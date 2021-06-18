// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "MahycoModule.h"


/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void MahycoModule::hydroStartInitEnvAndMat()
{
  info() << " Preparation des env  ";
  IMeshMaterialMng* m_material_mng;
  mm = IMeshMaterialMng::getReference(defaultMesh());
  // lecture des environnements
  MeshBlockBuildInfo mbbi("BLOCK1",allCells());
  UniqueArray<IMeshEnvironment*> saved_envs;
  IMeshMaterial* m_mat1;
  IMeshMaterial* m_mat2;
  IMeshEnvironment* m_env1;
  IMeshEnvironment* m_env2;

  info() << "Lit les infos des matériaux du JDD ";
  // Lit les infos des matériaux du JDD et les enregistre dans le gestionnaire
  for( Integer i=0,n=options()->material().size(); i<n; ++i ){
    String mat_name = options()->material[i].name;
    info() << "Found material name=" << mat_name;
    mm->registerMaterialInfo(mat_name);
  }

  for( Integer i=0,n=options()->environment().size(); i<n; ++i ){
    String env_name = options()->environment[i].name;
    info() << "Found environment name=" << env_name;
    Materials::MeshEnvironmentBuildInfo env_build(env_name);
    info() << " nbmat " << options()->environment[i].material.size();
    for( Integer k=0,kn=options()->environment[i].material.size(); k<kn; ++k ){
      String mat_name = options()->environment[i].material[k];
      info() << " k = " << k << " Add material " << mat_name << " for environment " << env_name;
      env_build.addMaterial(mat_name);
    }
    info() << "Materiau cree";
    IMeshEnvironment* env = mm->createEnvironment(env_build);
    info() << "Environment cree";
    saved_envs.add(env);
    // Le bloc ne contient que 2 milieux
    if (i<2){
      info() << "Add environment " << env_name << " to block1";
      mbbi.addEnvironment(env);
    }
  }
  
  info() << " Rangement des mailles  ";
  
  IMeshBlock* m_block1 = mm->createBlock(mbbi);
  
  mm->endCreate(subDomain()->isContinue());
  
  Integer nb_cell = allCells().size();
  Integer nb_env = m_block1->nbEnvironment();
  m_nb_env = nb_env;
  m_nb_vars_to_project = 3 * nb_env + 3 + 1 + 1;
  m_sens_projection = 0;
  // (volumes, ma;sses, energies internes) * nbmatmax
  // + vitesses + energie cinétique + pseudo*
  //
  // on redimensionne les tableaux de la projection en fonction
  // du nombre total d'environnements
  m_is_dir_face.resize(3); // dimension 3
  m_outer_face_normal.resize(6); // dimension 6 faces par mailles
  
  options()->remap()->resizeRemapVariables(m_nb_vars_to_project, m_nb_env);
  
  Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
   ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = one_over_nbnode * somme ;
  }
  
  info() << " Initialisation du cas test";
  
  options()->casModel()->initMat();
      
  m_mat1 = mm->environments()[0]->materials()[0];
  m_mat2 = mm->environments()[1]->materials()[0];
  m_env1 = mm->environments()[0];
  m_env2 = mm->environments()[1];
  
  
  Int32UniqueArray mat1_indexes;
  Int32UniqueArray mat2_indexes;
  info() << " Nb environments " << nb_env ;
  
  info() << " Trie par environnements  ";
  ENUMERATE_CELL(icell,allCells()){
    Cell cell = *icell;
    if (m_materiau[icell] == 0.) {
      mat1_indexes.add(icell.itemLocalId());
    } else if (m_materiau[icell] == 1.) {
      mat2_indexes.add(icell.itemLocalId());
    } else {
      mat1_indexes.add(icell.itemLocalId());
      mat2_indexes.add(icell.itemLocalId());
    }
        
  }

  info() << " Ajout de mailles mixtes " ;
 
  
  info() << " fin de la boucle sur les mailles  ";
  info() << " mat 1 " << mat1_indexes.size();
  info() << " mat 2 " << mat2_indexes.size();
  if (m_mat1){
    Materials::MeshMaterialModifier modifier(mm);
    modifier.addCells(m_mat1, mat1_indexes);
    info() << " mat1 modifie ";
  }    
  if (m_mat2){
    Materials::MeshMaterialModifier modifier(mm);
    modifier.addCells(m_mat2, mat2_indexes);
    info() << " mat2 modifie ";
  }
  
}
/**
 *******************************************************************************
 * \file PrepareFaceGroup()
 * \brief Creation des groupes de faces suivant X, Y et Z
 *******************************************************************************
 */
void MahycoModule::PrepareFaceGroup() {
  Int32UniqueArray face_x0_lid;
  Int32UniqueArray face_y0_lid;
  Int32UniqueArray face_z0_lid;
  Int32UniqueArray face_xmax_lid;
  Int32UniqueArray face_ymax_lid;
  Int32UniqueArray face_zmax_lid;
  Real3 ex = {1. , 0. , 0.};
  Real3 ey = {0. , 1. , 0.};
  Real3 ez = {0. , 0. , 1.};
  Real3 maxCoor= {-1. , -1. , -1.};
  ENUMERATE_NODE(inode, allNodes()){
      maxCoor.x = std::max(maxCoor.x, m_node_coord[inode].x);
      maxCoor.y = std::max(maxCoor.y, m_node_coord[inode].y);
      maxCoor.z = std::max(maxCoor.z, m_node_coord[inode].z);
  }
  ENUMERATE_FACE(iface, allFaces()){
     const Face& face = *iface;
     Integer face_local_id = face.localId();
     for (Integer idir = 0 ; idir <  mesh()->dimension() ; ++idir) {
         m_is_dir_face[face][idir] = false;
     }
     bool flag_x0(true);
     bool flag_y0(true);
     bool flag_z0(true);
     bool flag_xmax(true);
     bool flag_ymax(true);
     bool flag_zmax(true);
     for (NodeEnumerator inode(face.nodes()); inode.index() < face.nbNode(); ++inode) {
         if (m_node_coord[inode].x > options()->threshold)  flag_x0 = false;
         if (m_node_coord[inode].y > options()->threshold)  flag_y0 = false;
         if (m_node_coord[inode].z > options()->threshold)  flag_z0 = false;
         if (math::abs(m_node_coord[inode].x - maxCoor.x) > options()->threshold) flag_xmax = false;
         if (math::abs(m_node_coord[inode].y - maxCoor.y) > options()->threshold) flag_ymax = false;
         if (math::abs(m_node_coord[inode].z - maxCoor.z) > options()->threshold) flag_zmax = false;
     }
     if (flag_x0 == true) face_x0_lid.add(face_local_id);
     if (flag_y0 == true) face_y0_lid.add(face_local_id);
     if (flag_z0 == true) face_z0_lid.add(face_local_id);
     if (flag_xmax == true) face_xmax_lid.add(face_local_id);
     if (flag_ymax == true) face_ymax_lid.add(face_local_id);
     if (flag_zmax == true) face_zmax_lid.add(face_local_id);
   }
   
   mesh()->faceFamily()->createGroup("XMIN", face_x0_lid,true);
   mesh()->faceFamily()->createGroup("YMIN", face_y0_lid,true);
   mesh()->faceFamily()->createGroup("ZMIN", face_z0_lid,true);
   FaceGroup facexmin = mesh()->faceFamily()->findGroup("XMIN");
   FaceGroup faceymin = mesh()->faceFamily()->findGroup("YMIN");
   FaceGroup facezmin = mesh()->faceFamily()->findGroup("ZMIN");
   info() << " taille x 0 " << facexmin.size();
   info() << " taille y 0 " << faceymin.size();
   info() << " taille z 0 " << facezmin.size();
   
   mesh()->faceFamily()->createGroup("XMAX", face_xmax_lid,true);
   mesh()->faceFamily()->createGroup("YMAX", face_ymax_lid,true);
   mesh()->faceFamily()->createGroup("ZMAX", face_zmax_lid,true);
   FaceGroup facexmax = mesh()->faceFamily()->findGroup("XMAX");
   FaceGroup faceymax = mesh()->faceFamily()->findGroup("YMAX");
   FaceGroup facezmax = mesh()->faceFamily()->findGroup("ZMAX");
   info() << " taille x max " << facexmax.size();
   info() << " taille y max " << faceymax.size();
   info() << " taille z max " << facezmax.size();
   
 
   info() << " nombre total de face " << allFaces().size();
   
   info() << " creation des groupes de dimension " << m_dimension;
} 
