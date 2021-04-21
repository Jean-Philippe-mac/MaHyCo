﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"

#include <arcane/MathUtils.h>
#include <arcane/IParallelMng.h>
#include <arcane/ITimeLoopMng.h>

#include <arcane/geometry/IGeometry.h>
#include <arcane/mesh/GhostLayerMng.h>

using namespace Arcane;
using namespace Arcane::Materials;


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
hydroStartInit()
{
   IParallelMng* m_parallel_mng = subDomain()->parallelMng();
   my_rank = m_parallel_mng->commRank();
   
  
  pinfo() <<  "Mon rang " << my_rank << " et mon nombre de mailles " << allCells().size();
  pinfo() <<  " Mes mailles pures : " << ownCells().size();
  pinfo() <<  " Mes mailles frantomes : " << allCells().size() - ownCells().size();
  
  info() << " Check donnees ";
  if ((options()->ordreProjection == 3) && (mesh()->ghostLayerMng()->nbGhostLayer() != 3) && (m_parallel_mng->isParallel() == true)) {
      pinfo() << " mode parallele : " << m_parallel_mng->isParallel();
      pinfo() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
      pinfo() << " incompatible avec la projection d'ordre " << options()->ordreProjection;
      pinfo() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
      subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  if ((options()->withProjection == true) && (mesh()->ghostLayerMng()->nbGhostLayer() < 2) && (m_parallel_mng->isParallel() == true)) {
      pinfo() << " mode parallele : " << m_parallel_mng->isParallel();
      pinfo() << " nombre de couches de mailles fantomes : " << mesh()->ghostLayerMng()->nbGhostLayer();
      pinfo() << " incompatible avec la projection ";
      pinfo() << " ----------------------------- fin du calcul à la fin de l'init ---------------------------------------------";
      subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  
  m_cartesian_mesh = ICartesianMesh::getReference(mesh());
  // Dimensionne les variables tableaux
  m_cell_cqs.resize(8);
  m_cell_cqs_n.resize(8);
  
    // Initialise le delta-t
  Real deltat_init = options()->deltatInit();
  m_global_deltat = deltat_init;

  info() << " Initialisation des environnements";
  hydroStartInitEnvAndMat();
 
  // Initialise les données géométriques: volume, cqs, longueurs caractéristiques
  computeGeometricValues(); 
  
  info() << " Initialisation des groupes de faces";
  PrepareFaceGroup();
  
  info() << " Initialisation des variables";
  // Initialises les variables (surcharge l'init d'arcane)
  hydroStartInitVar();
  
  if (!options()->sansLagrange) {
    for( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
        IMeshEnvironment* ienv = mm->environments()[i];
        // Initialise l'énergie et la vitesse du son
        options()->environment[i].eosModel()->initEOS(ienv);
    }

    CellToAllEnvCellConverter all_env_cell_converter(mm);

    // valeur moyenne
    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if (all_env_cell.nbEnvironment() !=1) {
            m_internal_energy[cell] = 0.;
            ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
                EnvCell ev = *ienvcell;        
                m_internal_energy[cell] += m_internal_energy[ev] * m_mass_fraction[ev];
                m_sound_speed[cell] = std::max(m_sound_speed[ev], m_sound_speed[cell]);
            }
        }
    }
  }
  // initialisation du volume Euler
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    // initialisation du volume Euler
    m_euler_volume[cell] = m_cell_volume[cell];
  }
  
  Arcane::Numerics::IGeometryMng* gm = options()->geometry();
  gm->init();
  auto g = gm->geometry();
  bool is_verbose = false;
  if (is_verbose){
    ENUMERATE_CELL(icell,allCells()){
      Cell c = *icell;
      // initialisation du volume Euler
      info() << "Volume = " << g->computeMeasure(c);  
      for (NodeEnumerator inode(c.nodes()); inode.hasNext(); ++inode)
        info() << c.localId() << " avec " << inode.localId();
    }
  }
  
  InitGeometricValues();
  
  // deplacé de hydrocontinueInit
  // calcul des volumes via les cqs, differents deu calcul dans computeGeometricValues ?
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = *icell;
    // Calcule le volume de la maille et des noeuds
    Real volume = 0.0;
    for (Integer inode = 0; inode < 8; ++inode) {
      volume += math::dot(m_node_coord[cell.node(inode)], m_cell_cqs[icell] [inode]);
      m_node_volume[cell.node(inode)] += volume; 
    }
    volume /= 3.;
    m_cell_volume[icell] = volume;
  }
}

/**
 *******************************************************************************
 * \file computeCellMass()
 * \brief Calcul de la masse des mailles
 *
 * \param  m_cell_volume, m_density, m_mass_fraction_env
 * \return m_cell_mass, m_cell_mass_env
 *******************************************************************************
 */
void MahycoModule::
computeCellMass()
{
  debug() << " Entree dans computeCellMass()";
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    m_cell_mass[cell] = m_density[cell] * m_cell_volume[cell];
  }  
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
      m_cell_mass[ev] = m_mass_fraction[ev] * m_cell_mass[cell];
    }
  }
}
/**
 *******************************************************************************
 * \file computeNodeMass()
 * \brief Calcul de la masse nodale
 *
 * \param  m_cell_mass
 * \return m_node_mass
 *******************************************************************************
 */
void MahycoModule::
computeNodeMass() 
{
  debug() << " Entree dans computeNodeMass()";
   // Initialisation ou reinitialisation de la masse nodale
  m_node_mass.fill(0.);
  ENUMERATE_CELL(icell, allCells()){    
    Cell cell = * icell;
    Real contrib_node_mass = 0.125 * m_cell_mass[cell];
    for( NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode){
      m_node_mass[inode] += contrib_node_mass; 
    }
  }
  m_node_mass.synchronize();
}
/**
 *******************************************************************************
 * \file hydroContinueInit()
 * \brief Initialisation suite à une reprise
 *
 * \param  
 * \return m_nb_env, m_nb_vars_to_project, m_sens_projection
 *******************************************************************************
 */

void MahycoModule::
hydroContinueInit()
{
  if (subDomain()->isContinue()) {
    
    debug() << " Entree dans hydroContinueInit()";
    // en reprise 
    m_cartesian_mesh = ICartesianMesh::getReference(mesh());
    mm = IMeshMaterialMng::getReference(defaultMesh());
    info() << mm->name();
  
    mm->recreateFromDump();
    m_nb_env = mm->environments().size();
    m_nb_vars_to_project = 3 * m_nb_env + 3 + 1 + 1;
    
    m_global_old_deltat = m_old_deltat;
    // mise a jour nombre iteration 
    m_global_iteration = m_global_iteration() +1;
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void MahycoModule::
saveValuesAtN()
{
  debug() << " Entree dans saveValuesAtN()";
  
  // le pas de temps a été mis a jour a la fin dunpas de temps precedent et arcanne met dans m_global_old_deltat ce pas de temps ?
  // donc nous ont remet le bon old pas de temps
  m_global_old_deltat = m_old_deltat;
  
  // synchronisation debut de pas de temps (avec projection nécéssaire ?)
  m_pseudo_viscosity.synchronize();
  m_density.synchronize();
  m_internal_energy.synchronize();
  m_cell_volume.synchronize();
  m_pressure.synchronize();
  m_cell_cqs.synchronize();
  m_velocity.synchronize();
  
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = *icell;
    m_pseudo_viscosity_n[cell] = m_pseudo_viscosity[cell];
    m_pressure_n[cell] = m_pressure[cell];
    m_cell_volume_n[cell] = m_cell_volume[cell];
    m_density_n[cell] = m_density[cell];
    m_internal_energy_n[cell] = m_internal_energy[cell];
    // sauvegarde des cqs
    for (NodeEnumerator inode(cell.nodes()); inode.index() < 8; ++inode) {
      m_cell_cqs_n[icell] [inode.index()] = m_cell_cqs[icell] [inode.index()];
    }
  }
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      m_pseudo_viscosity_n[ev] = m_pseudo_viscosity[ev];
      m_pressure_n[ev] = m_pressure[ev];
      m_cell_volume_n[ev] = m_cell_volume[ev];
      m_density_n[ev] = m_density[ev];
      m_internal_energy_n[ev] = m_internal_energy[ev];
    }
  }
  if (!options()->sansLagrange) {
    ENUMERATE_NODE(inode, allNodes()){
        m_velocity_n[inode] = m_velocity[inode];
    }
  }
  if (options()->withProjection) {
    ENUMERATE_NODE(inode, allNodes()){
      m_node_coord[inode] = m_node_coord_0[inode];
    }
    // pas besoin de m_cell_coord car recalculé dans hydroContinueInit
  }
}    
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeArtificialViscosity()
{
  if (options()->sansLagrange) return;
  debug() << " Entree dans computeArtificialViscosity()";
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
      m_pseudo_viscosity[ev] = 0.;     
      if (m_div_u[cell] < 0.0) {
        m_pseudo_viscosity[ev] = 1. / m_tau_density[ev]
          * (-0.5 * m_caracteristic_length[cell] * m_sound_speed[cell] * m_div_u[cell]
             + (adiabatic_cst + 1) / 2.0 * m_caracteristic_length[cell] * m_caracteristic_length[cell]
             * m_div_u[cell] * m_div_u[cell]);
      }
    }
  }
  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_pseudo_viscosity[cell] = 0.;     
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        m_pseudo_viscosity[cell] += m_pseudo_viscosity[ev] * m_fracvol[ev];
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
updateVelocity()
{
  if (options()->sansLagrange) {
    updateVelocityWithoutLagrange();
    return;
  }
  debug() << " Entree dans updateVelocity()";
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::zero());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode) {
      m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
     }
  }

  const Real dt(0.5 * (m_global_old_deltat() + m_global_deltat()));
  // Calcule l'impulsion aux noeuds
  ENUMERATE_NODE(inode, allNodes()){
    Real node_mass = m_node_mass[inode];
    Real3 old_velocity = m_velocity_n[inode];
    Real3 new_velocity = old_velocity + ( dt / node_mass) * m_force[inode];
    m_velocity[inode] = new_velocity;
  }

  m_velocity.synchronize();
}
/**
 *******************************************************************************
 * \file updateVelocityBackward()
 * \brief Calcul de la vitesse de n a n-1/2
 *
 * \param  gt->deltat_n, m_node_velocity_n
 *         m_pressure_n, m_pseudo_viscosity_n, m_cqs_n
 * \return m_node_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateVelocityBackward()
{
  if (options()->sansLagrange) return;
  debug() << " Entree dans updateVelocityBackward()";
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::null());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure_n[icell] + m_pseudo_viscosity_n[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      m_force[inode] += pressure * m_cell_cqs_n[icell] [inode.index()];
  }
  const Real dt(-0.5 * m_global_old_deltat());
  // Calcule l'impulsion aux noeuds
  ENUMERATE_NODE(inode, allNodes()){
    Real node_mass = m_node_mass[inode];
    Real3 old_velocity = m_velocity_n[inode];
    Real3 new_velocity = old_velocity + ( dt / node_mass) * m_force[inode];
    m_velocity_n[inode] = new_velocity;
  }

  m_velocity_n.synchronize();
  
}
/*******************************************************************************
 * \file updateVelocityForward()
 * \brief Calcul de la vitesse de n+1/2 a n+1
 *
 * \param  gt->deltat_nplus, m_node_velocity_n
 *         m_pressure_nplus, m_pseudo_viscosity_nplus, m_cqs (calcule juste avant)
 * \return m_node_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateVelocityForward()
{
  if (options()->sansLagrange) return;
  debug() << " Entree dans updateVelocityForward()";
  // Remise à zéro du vecteur des forces.
  m_force.fill(Real3::null());

  // Calcul pour chaque noeud de chaque maille la contribution
  // des forces de pression
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real pressure = m_pressure[icell] + m_pseudo_viscosity[icell];
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      m_force[inode] += pressure * m_cell_cqs[icell] [inode.index()];
  }
  const Real dt(0.5 * m_global_deltat());
  // Calcule l'impulsion aux noeuds
  ENUMERATE_NODE(inode, allNodes()){
    Real node_mass = m_node_mass[inode];
    Real3 old_velocity = m_velocity[inode];
    Real3 new_velocity = old_velocity + ( dt / node_mass) * m_force[inode];
    m_velocity[inode] = new_velocity;
  }
  m_velocity.synchronize();
}
/**
 *******************************************************************************
 * \file updateVelocityWithoutLagrange()
 * \brief Calcul de la vitesse pour les cas d'advection pure
 *
 * \param  gt->t_nplus1, m_node_velocity_n0, m_node_velocity_n, m_node_coord_n
 * \return m_node_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateVelocityWithoutLagrange()
{
  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    m_velocity[node] = m_velocity_n[node];

    if (options()->casTest == RiderVortexTimeReverse ||
        options()->casTest == MonoRiderVortexTimeReverse ) {
      m_velocity[node].x =
          m_velocity_n[node].x * cos(Pi * m_global_time() / 4.);
      m_velocity[node].y =
          m_velocity_n[node].y * cos(Pi * m_global_time() / 4.);
    } else if ( options()->casTest == RiderDeformationTimeReverse ||
		options()->casTest == MonoRiderDeformationTimeReverse) {
      m_velocity[node].x =
          m_velocity_n[node].x * cos(Pi * m_global_time());
      m_velocity[node].y =
          m_velocity_n[node].y * cos(Pi * m_global_time());
    }
  }
  m_velocity.synchronize();
}
/**
 *******************************************************************************
 * \file updatePosition()
 * \brief Calcul de la position
 *
 * \param  gt->t_nplus1, m_node_velocity_nplus1, m_node_coord_n
 * \return m_node_coord_nplus1
 *******************************************************************************
 */
void MahycoModule::
updatePosition()
{
  debug() << " Entree dans updatePosition()";
  Real deltat = m_global_deltat();
  ENUMERATE_NODE(inode, allNodes()){
    Node node = *inode;
    if (((options()->sansLagrange) && (node.nbCell() == 4)) || (!options()->sansLagrange))
        m_node_coord[inode] += deltat * m_velocity[inode];
  }
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = 0.125 * somme;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
applyBoundaryCondition()
{
  debug() << " Entree dans applyBoundaryCondition()";
  for (Integer i = 0, nb = options()->boundaryCondition.size(); i < nb; ++i){
    String NomBC = options()->boundaryCondition[i]->surface;
    FaceGroup face_group = mesh()->faceFamily()->findGroup(NomBC);
    Real value = options()->boundaryCondition[i]->value();
    TypesMahyco::eBoundaryCondition type = options()->boundaryCondition[i]->type();

    // boucle sur les faces de la surface
    ENUMERATE_FACE(j, face_group){
      Face face = * j;
      Integer nb_node = face.nbNode();

      // boucle sur les noeuds de la face
      for (Integer k = 0; k < nb_node; ++k){
        Node node = face.node(k);
        Real3& velocity = m_velocity[node];

        switch (type){
        case TypesMahyco::VelocityX:
          velocity.x = value;
          break;
        case TypesMahyco::VelocityY:
          velocity.y = value;
          break;
        case TypesMahyco::VelocityZ:
          velocity.z = value;
          break;
        case TypesMahyco::Unknown:
          break;
        }
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
InitGeometricValues()
{
  debug() << " Entree dans InitGeometricValues() ";
  ENUMERATE_NODE(inode, allNodes()){
      m_node_coord_0[inode] = m_node_coord[inode];
  }
  ENUMERATE_FACE (iFace, allFaces()) {
    Face face = *iFace;
    Real3 face_vec1 = m_node_coord[face.node(2)] - m_node_coord[face.node(0)]; 
    Real3 face_vec2 = m_node_coord[face.node(3)] - m_node_coord[face.node(1)];
    m_face_normal[iFace].x = produit(face_vec1.y, face_vec2.z, face_vec1.z, face_vec2.y);
    m_face_normal[iFace].y = - produit(face_vec2.x, face_vec1.z, face_vec2.z, face_vec1.x);
    m_face_normal[iFace].z = produit(face_vec1.x, face_vec2.y, face_vec1.y, face_vec2.x);
    m_face_normal[iFace] /= m_face_normal[iFace].abs();
    m_face_coord[face] = 0.;
    for (Integer inode = 0; inode < 4; ++inode) 
      m_face_coord[face] +=  0.25 * m_node_coord[face.node(inode)];
  }
  ENUMERATE_CELL(icell,allCells()) {
    ENUMERATE_FACE(iface, (*icell).faces()){
      const Face& face = *iface;
      Integer index = iface.index(); 
        m_outer_face_normal[icell][index] = (m_face_coord[face]-m_cell_coord[icell]) 
            / (m_face_coord[face]-m_cell_coord[icell]).abs();
    }
  }
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeGeometricValues()
{
  debug() << my_rank << " : " << " Entree dans computeGeometricValues() ";
  // Copie locale des coordonnées des sommets d'une maille
  Real3 coord[8];
  // Coordonnées des centres des faces
  Real3 face_coord[6];
  
  m_node_coord.synchronize();
  
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    // Recopie les coordonnées locales (pour le cache)
    for (NodeEnumerator inode(cell.nodes()); inode.index() < 8; ++inode) {
      coord[inode.index()] = m_node_coord[inode];
    }
    // Calcul les coordonnées des centres des faces
    face_coord[0] = 0.25 * (coord[0] + coord[3] + coord[2] + coord[1]);
    face_coord[1] = 0.25 * (coord[0] + coord[4] + coord[7] + coord[3]);
    face_coord[2] = 0.25 * (coord[0] + coord[1] + coord[5] + coord[4]);
    face_coord[3] = 0.25 * (coord[4] + coord[5] + coord[6] + coord[7]);
    face_coord[4] = 0.25 * (coord[1] + coord[2] + coord[6] + coord[5]);
    face_coord[5] = 0.25 * (coord[2] + coord[3] + coord[7] + coord[6]);

    // Calcule les résultantes aux sommets
    computeCQs(coord, face_coord, cell);
  }
  m_cell_cqs.synchronize();
  
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;    
    // Calcule le volume de la maille
    {
      Real volume = 0.;
      
      for (Integer inode = 0; inode < 8; ++inode) {
        volume += math::dot(m_node_coord[cell.node(inode)], m_cell_cqs[icell] [inode]);
        m_node_volume[cell.node(inode)] += volume;
      }
      volume /= 3.;
      
      m_cell_volume[cell] = volume;
    }
    
    // Calcule la longueur caractéristique de la maille.
    {
        if (options()->longueurCaracteristique() == "faces-opposees")
        {
            // Recopie les coordonnées locales (pour le cache)
            for (NodeEnumerator inode(cell.nodes()); inode.index() < 8; ++inode) {
                coord[inode.index()] = m_node_coord[inode];
            }
            // Calcul les coordonnées des centres des faces
            face_coord[0] = 0.25 * (coord[0] + coord[3] + coord[2] + coord[1]);
            face_coord[1] = 0.25 * (coord[0] + coord[4] + coord[7] + coord[3]);
            face_coord[2] = 0.25 * (coord[0] + coord[1] + coord[5] + coord[4]);
            face_coord[3] = 0.25 * (coord[4] + coord[5] + coord[6] + coord[7]);
            face_coord[4] = 0.25 * (coord[1] + coord[2] + coord[6] + coord[5]);
            face_coord[5] = 0.25 * (coord[2] + coord[3] + coord[7] + coord[6]);
            
            Real3 median1 = face_coord[0] - face_coord[3];
            Real3 median2 = face_coord[2] - face_coord[5];
            Real3 median3 = face_coord[1] - face_coord[4];
            Real d1 = median1.abs();
            Real d2 = median2.abs();
            Real d3 = median3.abs();
            
            Real dx_numerator = d1 * d2 * d3;
            Real dx_denominator = d1 * d2 + d1 * d3 + d2 * d3;
            m_caracteristic_length[icell] = dx_numerator / dx_denominator;
        }
        else if (options()->longueurCaracteristique() == "racine-cubique-volume")
        {
            m_caracteristic_length[icell] = std::pow(m_cell_volume[icell], 1./3.);
        }
        else
        {
            info() << " pas de longeur caractéritique definie dans le .arc " << options()->longueurCaracteristique(); 
            subDomain()->timeLoopMng()->stopComputeLoop(true);
        }
    }
 
  }

  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        Cell cell = ev.globalCell();
        m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
      }
    }
  }
}
/**
 *******************************************************************************
 * \file computeDivU()
 * \brief Calcul de la densité
 * \brief Calcul de la divergence de la vitesse
 * \brief Calcul de la variation du volume specifique (variation de 1/densite)
 *
 * pour les grandeurs moyennes et pour chaque materiau
 *
 * \param  m_density_nplus1, m_density_n, gt->deltat_nplus1,
 *        m_tau_density_nplus1 
 * \return m_divu_nplus1
 *******************************************************************************
 */
void MahycoModule::
updateDensity()
{
  if (options()->sansLagrange) return;
  debug() << my_rank << " : " << " Entree dans updateDensity() ";
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
       // pinfo() << my_rank << " : " << cell.uniqueId() << " " << m_cell_volume[ev];
       Real new_density = m_cell_mass[ev] / m_cell_volume[ev];
       // pinfo() << my_rank << " : " << cell.uniqueId() << " " << m_density_n[ev];
       // nouvelle density
       m_density[ev] = new_density;
       // volume specifique de l'environnement au temps n+1/2
       m_tau_density[ev] = 
        0.5 * (1.0 / m_density_n[ev] + 1.0 / m_density[ev]);
    }
  }
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = * icell;
    Real new_density = m_cell_mass[icell] / m_cell_volume[icell];
    // nouvelle density
    m_density[icell] = new_density;
    // volume specifique moyen au temps n+1/2
    m_tau_density[icell] = 
        0.5 * (1.0 / m_density_n[icell] + 1.0 / m_density[icell]);
    // divergence de la vitesse mode A1
    m_div_u[icell] =
      1.0 / m_global_deltat()  * ( 1.0 / m_density[icell] - 1.0 / m_density_n[icell] )
      / m_tau_density[icell];
  }
  
  m_density.synchronize();
  m_tau_density.synchronize();
  m_div_u.synchronize();
}
/**
 *******************************************************************************
 * \file updateEnergy()
 * \brief Calcul de l'energie interne (seul le cas du gaz parfait est codé)
 *
 * \param  m_density_env_nplus1, m_density_env_n,
 *         m_pseudo_viscosity_env_nplus1, m_pseudo_viscosity_env_n
 *         m_pressure_env_n, m_mass_fraction_env
 *
 * \return m_internal_energy_env_nplus1, m_internal_energy_nplus1
 *******************************************************************************
 */
/** */
/* la fonction dont on cherche un zero */
/** */
inline double fvnr(double e, double p, double dpde,
		double en, double qnn1, double pn, double rn1,
		double rn) 
{return e-en+0.5*(p+pn+2.*qnn1)*(1./rn1-1./rn);};

/** */
/* la derivee de f */
/** */
inline double fvnrderiv(double e, double dpde, double rn1, double rn)
{return 1.+0.5*dpde*(1./rn1-1./rn);};
/*
 *******************************************************************************
 */
void MahycoModule::
updateEnergyAndPressure()
{
  if (options()->sansLagrange) return;
  debug() << " Entree dans updateEnergyAndPressure()";
  bool csts = options()->schemaCsts();
  bool pseudo_centree = options()->pseudoCentree();
  // Calcul de l'énergie interne
  if (!csts) {
    ENUMERATE_ENV(ienv,mm){
      IMeshEnvironment* env = *ienv;
      Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Real pseudo(0.);
        if (pseudo_centree &&
            ((m_pseudo_viscosity_n[ev] + m_pseudo_viscosity[ev]) * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
          pseudo = 0.5 * (m_pseudo_viscosity[ev] + m_pseudo_viscosity_n[ev]);
        if (!pseudo_centree &&
            (m_pseudo_viscosity[ev] * (1.0 / m_density[ev] - 1.0 / m_density_n[ev]) < 0.))
          pseudo = m_pseudo_viscosity[ev];
          
        Real denom_accrois_nrj(1 + 0.5 * (adiabatic_cst - 1.0) *
                               m_density[ev] *
                               (1.0 / m_density[ev] -
                                1.0 / m_density_n[ev]));
        Real numer_accrois_nrj(m_internal_energy_n[ev] -
                               (0.5 * m_pressure[ev] + pseudo) *
                               (1.0 / m_density[ev] - 1.0 / m_density_n[ev]));
        m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
      }
    }
  } else {
    ENUMERATE_ENV(ienv,mm){
      IMeshEnvironment* env = *ienv;
      Real adiabatic_cst = options()->environment[env->id()].eosModel()->getAdiabaticCst(env);
      ENUMERATE_ENVCELL(ienvcell,env){
        EnvCell ev = *ienvcell;
        Cell cell = ev.globalCell();
        Real cqs_v_nplus1(0.);
        Real cqs_v_n(0.);
        for (Integer inode = 0; inode < 8; ++inode) {
          cqs_v_nplus1 += math::dot(m_velocity[cell.node(inode)], m_cell_cqs[cell] [inode])
            * m_global_deltat();
          cqs_v_n += math::dot(m_velocity[cell.node(inode)], m_cell_cqs_n[cell] [inode])
            * m_global_deltat();
        }
        Real denom_accrois_nrj(1 + 0.5 * (adiabatic_cst - 1.0)
                               * m_density[ev] * cqs_v_n / m_cell_mass[ev]) ;
        Real numer_accrois_nrj(m_internal_energy[ev]
                               - (0.5 * (m_pressure[ev] + m_pseudo_viscosity[ev])
                                  * cqs_v_n / m_cell_mass[ev])
                               - (0.5 * m_pseudo_viscosity[ev]
                                  * cqs_v_nplus1 / m_cell_mass[ev]));	
        m_internal_energy[ev] = numer_accrois_nrj / denom_accrois_nrj;
      }
    }
  }
  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_internal_energy[cell] = 0.;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;
        m_internal_energy[cell] += m_mass_fraction[ev] * m_internal_energy[ev];
      }
    }
  }
  if (! options()->withProjection) {
    // Calcul de la Pression si on ne fait pas de projection 
    for( Integer i=0,n=options()->environment().size(); i<n; ++i ) {
        IMeshEnvironment* ienv = mm->environments()[i];
        // Calcul de la pression et de la vitesse du son
        options()->environment[i].eosModel()->applyEOS(ienv);
    }
    computePressionMoyenne();
  }
}   
/**
 *******************************************************************************
 * \file computePressionMoyenne()
 * \brief Calcul de la pression moyenne
 *
 * \param  m_fracvol_env, m_pressure_env_nplus1, m_speed_velocity_env_nplus1
 * \return m_pressure_nplus1, m_speed_velocity_nplus1
 *******************************************************************************
 */
void MahycoModule::
computePressionMoyenne()
{
  if (options()->sansLagrange) return;
  debug() << " Entree dans computePressionMoyenne() ";
  // maille mixte
  // moyenne sur la maille
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() !=1) {
      m_pressure[cell] = 0.;
      m_sound_speed[icell] = 1.e-20;
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;        
        m_pressure[cell] += m_fracvol[ev] * m_pressure[ev];
        m_sound_speed[cell] = std::max(m_sound_speed[ev], m_sound_speed[cell]);
      }
    }
  }
}     
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::
computeDeltaT()
{
  debug() << " Entree dans compute DT avec " << m_global_old_deltat()
         << " et " << options()->deltatInit()
         << " et " << m_global_deltat()
         << " et " << options()->deltatMax();
         
  m_global_old_deltat = m_global_deltat;
  m_old_deltat = m_global_old_deltat();
         
  Real new_dt = FloatInfo < Real >::maxValue();
  if (options()->sansLagrange) {
    // on garde le meme pas de temps
    new_dt = options()->deltatInit();
    
  } else {
    CellToAllEnvCellConverter all_env_cell_converter(mm);

    // Calcul du pas de temps pour le respect du critère de CFL
    Real minimum_aux = FloatInfo < Real >::maxValue();
    Integer cell_id(-1), nbenvcell(-1);
    Real cc(0.), ll(0.);

    ENUMERATE_CELL(icell, allCells()){
        Cell cell = * icell;
        Real cell_dx = m_caracteristic_length[icell];
        Real sound_speed = m_sound_speed[icell];
        Real vmax(0.);
        if (options()->withProjection)
          for (NodeEnumerator inode(cell.nodes()); inode.index() < 8; ++inode) {
            vmax = math::max(m_velocity[inode].abs(), vmax);
          }
        Real dx_sound = cell_dx / (sound_speed + vmax);
        minimum_aux = math::min(minimum_aux, dx_sound);
        if (minimum_aux == dx_sound) {
            cell_id = icell.localId();
            cc = m_sound_speed[icell];
            ll = m_caracteristic_length[icell];
            AllEnvCell all_env_cell = all_env_cell_converter[cell];
            nbenvcell = all_env_cell.nbEnvironment();
        }
    }

    new_dt = options()->cfl() * minimum_aux;
    new_dt = parallelMng()->reduce(Parallel::ReduceMin, new_dt);
    
    // respect de taux de croissance max
    new_dt = math::min(new_dt, 1.05*m_global_old_deltat());
    // respect de la valeur max imposée par le fichier de données .plt
    debug() << " nouveau pas de temps " << new_dt << " par " << cell_id << 
        " (avec " << nbenvcell << " envs) avec " << cc << " " << ll << " et min " << minimum_aux;
    new_dt = math::min(new_dt, options()->deltatMax());
    // respect du pas de temps minimum
    if (new_dt < options()->deltatMin()) {
        info() << " pas de temps minimum ";
        info() << " nouveau pas de temps " << new_dt << " par " << cell_id 
            << " (avec " << nbenvcell << " envs) avec vitson = " << cc << " et longeur  =  " << ll << " et min " << minimum_aux;
        exit(1);
    }
    
    debug() << " nouveau pas de temps2 " << new_dt;
        
  }
  
  // Mise à jour du pas de temps
  m_global_deltat = new_dt;
  
  // Le dernier calcul se fait exactement au temps stopTime()
  Real stop_time  = options()->finalTime();
  bool not_yet_finish = (m_global_time() < stop_time);
  bool finish = (m_global_time() > stop_time);
  bool too_much = ((m_global_time()+new_dt) > stop_time);

  
  debug() << " nouveau pas de temps " << new_dt;
  if ( not_yet_finish && too_much){
    new_dt = stop_time - m_global_time();
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
  if (finish) {
    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }
    
  debug() << " time " << m_global_time() << " et fin à " << stop_time;
  debug() << " not_yet_finish " << not_yet_finish;
  debug() << " finish " << finish;
  debug() << " too_much " << too_much;
  
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline void MahycoModule::
computeCQs(Real3 node_coord[8], Real3 face_coord[6], const Cell & cell)
{
  const Real3 c0 = face_coord[0];
  const Real3 c1 = face_coord[1];
  const Real3 c2 = face_coord[2];
  const Real3 c3 = face_coord[3];
  const Real3 c4 = face_coord[4];
  const Real3 c5 = face_coord[5];

  // Calcul des normales face 1 :
  const Real3 n1a04 = 0.5 * math::vecMul(node_coord[0] - c0, node_coord[3] - c0);
  const Real3 n1a03 = 0.5 * math::vecMul(node_coord[3] - c0, node_coord[2] - c0);
  const Real3 n1a02 = 0.5 * math::vecMul(node_coord[2] - c0, node_coord[1] - c0);
  const Real3 n1a01 = 0.5 * math::vecMul(node_coord[1] - c0, node_coord[0] - c0);

  // Calcul des normales face 2 :
  const Real3 n2a05 = 0.5 * math::vecMul(node_coord[0] - c1, node_coord[4] - c1);
  const Real3 n2a12 = 0.5 * math::vecMul(node_coord[4] - c1, node_coord[7] - c1);
  const Real3 n2a08 = 0.5 * math::vecMul(node_coord[7] - c1, node_coord[3] - c1);
  const Real3 n2a04 = 0.5 * math::vecMul(node_coord[3] - c1, node_coord[0] - c1);

  // Calcul des normales face 3 :
  const Real3 n3a01 = 0.5 * math::vecMul(node_coord[0] - c2, node_coord[1] - c2);
  const Real3 n3a06 = 0.5 * math::vecMul(node_coord[1] - c2, node_coord[5] - c2);
  const Real3 n3a09 = 0.5 * math::vecMul(node_coord[5] - c2, node_coord[4] - c2);
  const Real3 n3a05 = 0.5 * math::vecMul(node_coord[4] - c2, node_coord[0] - c2);

  // Calcul des normales face 4 :
  const Real3 n4a09 = 0.5 * math::vecMul(node_coord[4] - c3, node_coord[5] - c3);
  const Real3 n4a10 = 0.5 * math::vecMul(node_coord[5] - c3, node_coord[6] - c3);
  const Real3 n4a11 = 0.5 * math::vecMul(node_coord[6] - c3, node_coord[7] - c3);
  const Real3 n4a12 = 0.5 * math::vecMul(node_coord[7] - c3, node_coord[4] - c3);

  // Calcul des normales face 5 :
  const Real3 n5a02 = 0.5 * math::vecMul(node_coord[1] - c4, node_coord[2] - c4);
  const Real3 n5a07 = 0.5 * math::vecMul(node_coord[2] - c4, node_coord[6] - c4);
  const Real3 n5a10 = 0.5 * math::vecMul(node_coord[6] - c4, node_coord[5] - c4);
  const Real3 n5a06 = 0.5 * math::vecMul(node_coord[5] - c4, node_coord[1] - c4);

  // Calcul des normales face 6 :
  const Real3 n6a03 = 0.5 * math::vecMul(node_coord[2] - c5, node_coord[3] - c5);
  const Real3 n6a08 = 0.5 * math::vecMul(node_coord[3] - c5, node_coord[7] - c5);
  const Real3 n6a11 = 0.5 * math::vecMul(node_coord[7] - c5, node_coord[6] - c5);
  const Real3 n6a07 = 0.5 * math::vecMul(node_coord[6] - c5, node_coord[2] - c5);

  // Calcul des résultantes aux sommets :
  m_cell_cqs[cell] [0] = (5. * (n1a01 + n1a04 + n2a04 + n2a05 + n3a05 + n3a01) +
                          (n1a02 + n1a03 + n2a08 + n2a12 + n3a06 + n3a09)) * (1. / 12.);
  m_cell_cqs[cell] [1] = (5. * (n1a01 + n1a02 + n3a01 + n3a06 + n5a06 + n5a02) +
                          (n1a04 + n1a03 + n3a09 + n3a05 + n5a10 + n5a07)) * (1. / 12.);
  m_cell_cqs[cell] [2] = (5. * (n1a02 + n1a03 + n5a07 + n5a02 + n6a07 + n6a03) +
                          (n1a01 + n1a04 + n5a06 + n5a10 + n6a11 + n6a08)) * (1. / 12.);
  m_cell_cqs[cell] [3] = (5. * (n1a03 + n1a04 + n2a08 + n2a04 + n6a08 + n6a03) +
                          (n1a01 + n1a02 + n2a05 + n2a12 + n6a07 + n6a11)) * (1. / 12.);
  m_cell_cqs[cell] [4] = (5. * (n2a05 + n2a12 + n3a05 + n3a09 + n4a09 + n4a12) +
                          (n2a08 + n2a04 + n3a01 + n3a06 + n4a10 + n4a11)) * (1. / 12.);
  m_cell_cqs[cell] [5] = (5. * (n3a06 + n3a09 + n4a09 + n4a10 + n5a10 + n5a06) +
                          (n3a01 + n3a05 + n4a12 + n4a11 + n5a07 + n5a02)) * (1. / 12.);
  m_cell_cqs[cell] [6] = (5. * (n4a11 + n4a10 + n5a10 + n5a07 + n6a07 + n6a11) +
                          (n4a12 + n4a09 + n5a06 + n5a02 + n6a03 + n6a08)) * (1. / 12.);
  m_cell_cqs[cell] [7] = (5. * (n2a08 + n2a12 + n4a12 + n4a11 + n6a11 + n6a08) +
                          (n2a04 + n2a05 + n4a09 + n4a10 + n6a07 + n6a03)) * (1. / 12.);
}

/*---------------------------------------------------------------------------*/
inline Real MahycoModule::produit(Real A, Real B, Real C, Real D)
  {
  return (A*B-C*D);
  }
/*---------------------------------------------------------------------------*/
inline Real MahycoModule::norme(Real E, Real F, Real G)
  {
  return (sqrt((E*E)+(F*F)+(G*G)));
  }

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_MAHYCO(MahycoModule);
