#include <iostream>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <cmath>
#include <string>
#include "MahycoModule.h"

using namespace :: std;

 

void MahycoModule::calcCSV() {
  //calcul du centre de la maille : m_c_2  
  Real one_over_nbnode = m_dimension == 2 ? .25  : .125 ;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = one_over_nbnode * somme;
  }
  // calcul des volumes : m_V2
  computeGeometricValues();
  // calcul des normales : m_s_2
  one_over_nbnode = m_dimension == 2 ? .5  : .25 ;
  ENUMERATE_FACE (iFace, allFaces()) {
      Face face = *iFace;
      m_face_coord[face] = 0.;
      for (Integer inode = 0; inode < face.nbNode(); ++inode) 
        m_face_coord[face] +=  one_over_nbnode * m_node_coord[face.node(inode)];
  }
  ENUMERATE_CELL(icell, allCells()) {
    ENUMERATE_FACE(iface, (*icell).faces()){
      const Face& face = *iface;
      Integer index = iface.index(); 
      m_outer_face_normal[icell][index] = (m_face_coord[face]-m_cell_coord[icell]) 
            / (m_face_coord[face]-m_cell_coord[icell]).abs();
    }    
  }
  // plus petit longeur de maille est dans m_caracteristic_length[icell] : m_min_l
}
void MahycoModule::calcSigVrl_2() {
  //calcul des decentrement  en n+1/2 : m_sig_2
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    Real Vdc(0.);
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
       EnvCell ev = *ienvcell;   
       Integer index_env = ev.environmentId();     
       for (Integer index = 0; index < cell.nbFace(); ++index) {
         // Calcul des volumes des materiaux (environnements)
         m_cell_volume[ev] = max_0(m_outer_face_normal[icell][index].x * m_matcell_relativevelocity[ev].x +
          m_outer_face_normal[icell][index].y * m_matcell_relativevelocity[ev].y);
         // Calcul des coefficient de decentrement standart : m_mat_sigma (m_sig)
         // Face de cet index 
         Face face = cell.face(index);
         // Cell voisine vis-à-vis de cette face 
         Cell cellvois = getCellVoisine(cell, face);
         // index local dans cellvois de la face face
         Integer indexvois = getindiceface(cellvois, face);
         // obtenir l'envcell de cette maille voisine 
         EnvCell evvois = getenvcell(cellvois, index_env);
         Vdc = max_0(m_outer_face_normal[cellvois][indexvois].x * m_matcell_relativevelocity[evvois].x -
                   m_outer_face_normal[cellvois][indexvois].y * m_matcell_relativevelocity[evvois].y);
         if (m_cell_volume[ev] == 0. && Vdc == 0.) {
           m_mat_sigma[ev][index] = 0.5;
         } else {
          m_mat_sigma[ev][index] = m_cell_volume[ev] / ( m_cell_volume[ev] + Vdc);
         }
       }
      }
    }
  }
}
void MahycoModule::calcSigT() {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  //calcul des decentrement globaux en n+1/2 : m_sigT_2
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        m_total_mat_sigma[ev][index] = m_mat_sigma[ev][index];
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        m_total_mat_sigma[ev][index] += m_fracvol[evvois] * m_mat_sigma[evvois][indexvois];
      }
    }
  }
}
void MahycoModule::calcQ() {
  Real a1(0.), a2(0.);
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  // Calcul de la divergence du volume en n-1/2
  // utilisation des décentrement en n-1/2 : m_sigT soit m_total_mat_sigma_n
  m_div_u_n.fill(0.);
  Real r = 0.0;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        m_total_mat_sigma_n[ev][index] = m_mat_sigma_n[ev][index];
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        m_div_u_n[cell] += 0.5 * math::dot( m_outer_face_normal_n[cell][index], 
           (m_fracvol[ev] * m_total_mat_sigma_n[ev][index] * m_matcell_velocity_n[ev] + 
            m_fracvol[evvois] * m_total_mat_sigma_n[evvois][indexvois] * m_matcell_velocity_n[evvois]));
      }
      r +=  m_alpha_rho[ev];
    }
    //Calcul de la viscosite artificielle        
    m_pseudo_viscosity[cell] = r * (a2 * math::pow(m_div_u_n[cell], 2.) - a1 * m_sound_speed[cell] * m_div_u_n[cell]);
  }
}
//--------------------------------------------------------------------
//----------Evolution-------------------------------------------------
//--------------------------------------------------------------------	
void MahycoModule::evolueU() {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  Real deltat = m_global_deltat(); // dtnp12 de n a n+1
  const Real dt(0.5 * (m_global_old_deltat() + m_global_deltat())); // dtn de n-1/2 a n+1/2
  RealArray* coeff = new RealUniqueArray(m_nb_env);
  Real3Array* X = new Real3UniqueArray(m_nb_env);
  // m_total_mat_sigma_n.copy(m_mat_sigma_n);
  //Stokes coefficient
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    Real density0;
    Real3 matcellvelocity0(Real3(0., 0., 0.));
    Real3 cumul_stockes(Real3(0., 0., 0.));
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      //debut de calcul des decentrement globaux en n-1/2 
      m_total_mat_sigma_n[ev] = m_mat_sigma_n[ev];
      Integer index_env = ev.environmentId();   
      if (index_env == 0) {
        matcellvelocity0 = m_matcell_velocity[ev];
        density0 = m_density[ev];
      } else {
        (*coeff)[index_env] = dt * m_coef_s * 9 * m_fracvol[ev] * density0 * m_vu_air / (2 * m_r_eau * m_r_eau);
        m_stockes[ev] = - m_cell_volume[ev] * (*coeff)[index_env] * (m_matcell_velocity[ev] - matcellvelocity0);
        cumul_stockes += m_stockes[ev];
      }      
    }
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      if (index_env == 0) m_stockes[ev] = cumul_stockes;
    }
  }
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    Real3 transport(Real3(0., 0., 0.));
    Real3 grad_P(Real3(0., 0., 0.));
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        m_total_mat_sigma_n[ev][index] = m_mat_sigma_n[ev][index];
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        m_total_mat_sigma_n[ev][index] += m_fracvol[evvois] * m_mat_sigma_n[evvois][indexvois];
      }
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        transport += m_volumic_transfert[ev][index] * m_alpha_rho_old[ev] * m_matcell_velocity_n[ev] 
                   - m_volumic_transfert[evvois][indexvois] * m_alpha_rho_old[evvois] * m_matcell_velocity_n[evvois];
        grad_P += 0.5 * m_fracvol[ev] * m_total_mat_sigma_n[ev][index] * m_outer_face_normal[cell][index] 
              * ( m_pressure[evvois] - m_pressure[ev] + m_pseudo_viscosity[evvois] - m_pseudo_viscosity[ev]);
      }
      (*X)[index_env] = (m_alpha_rho_old[ev] * m_matcell_velocity_n[ev] * m_cell_volume_old[ev] - deltat * transport - dt * grad_P) / (m_cell_volume_n[ev] * m_alpha_rho[ev]);
    }
    //Avec Stokes version
	Real temp_1(0.);
	Real3 temp_2(Real3(0., 0., 0.));
    
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      temp_1 += m_alpha_rho[ev] * (*coeff)[index_env] / (m_alpha_rho[ev] * (*coeff)[index_env]);
      temp_2 += m_alpha_rho[ev] * (*coeff)[index_env] * (*X)[index_env] / (m_alpha_rho[ev] * (*coeff)[index_env]);
    }
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      if (index_env == 0) temp_1 += m_alpha_rho[ev];
    }
    Real3 matcellvelocity0(Real3(0., 0., 0.));
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      if (index_env == 0) {
	    // Calcul de la vitesse de la phase porteuse
        m_matcell_velocity[ev] = (temp_2 + m_alpha_rho[ev] * (*X)[index_env]) / temp_1;
        matcellvelocity0 = m_matcell_velocity[ev];
      } else {
	    // Calcul de la vitesse des fluides disperses
        m_matcell_velocity[ev] = (m_alpha_rho[ev] * (*X)[index_env] * m_matcell_velocity[ev]) 
        / (m_alpha_rho[ev] + (*coeff)[index_env]);
      }
      m_matcell_relativevelocity[ev] = m_matcell_velocity[ev] - (m_cell_coord[cell] - m_cell_coord_n[cell]) / dt;
    }
  }
}
void MahycoModule::evolueAR() {
  Real deltat = m_global_deltat(); // dtnp12 de n a n+1
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      m_alpha_rho[ev] = m_cell_volume_n[ev] * m_alpha_rho_n[ev] / m_cell_volume[ev];
      Integer index_env = ev.environmentId();   
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        m_total_mat_sigma_n[ev][index] = m_mat_sigma_n[ev][index];
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        m_alpha_rho[ev] +=  deltat * (m_volumic_transfert[ev][index] * m_alpha_rho_n[ev] -
            m_volumic_transfert[evvois][indexvois] * m_alpha_rho_n[evvois]) / m_cell_volume[ev];
      }      
    }
  }
}

void MahycoModule::evolueE() {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  Real deltat = m_global_deltat(); // dtnp12 de n a n+1
  const Real dt(0.5 * (m_global_old_deltat() + m_global_deltat())); // dtn de n-1/2 a n+1/2
  m_div_u.fill(0.);
  Real ech_P, ech_W = 0.0;
  Real transport;
  Real gamPhi, gamF;
  Real gruneisen, gruneisen_phi, p_ref, pphi_ref;
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
	Real temp(0.);
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();  
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      temp += m_fracvol[ev] / gruneisen;
    }
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      m_dissipation[ev] = - dt * m_fracvol[ev] / (temp * gruneisen) * m_pseudo_viscosity[ev] * m_div_u[cell];
      // Calcul de la divergence du volume
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        
        // Calcul de la divergence du volume
        m_div_u[cell] += 0.25 * math::dot(m_outer_face_normal_n[cell][index],
        ( m_fracvol[ev] * m_total_mat_sigma[ev][index] * 
        (m_matcell_velocity_n[ev] + m_matcell_velocity[ev]) +
        m_fracvol[evvois] * m_total_mat_sigma[evvois][indexvois] * 
        (m_matcell_velocity_n[evvois] + m_matcell_velocity[evvois])));
        
        // Calcul de l'advection du gradient de pression
        m_pressure_advection[ev] += 0.25 * m_total_mat_sigma[ev][index] * math::dot((m_matcell_velocity[ev] 
                 + m_matcell_velocity_n[ev]), m_outer_face_normal_n[cell][index])
                 * (m_pressure[evvois] - m_pressure[ev]);
        // Calcul de la dissipation
        m_dissipation[ev] += deltat * m_volumic_transfert[evvois][indexvois] *  0.5 * m_alpha_rho_old[evvois] 
        * math::dot((m_matcell_velocity[ev] - m_matcell_velocity_n[evvois]),  
                    (m_matcell_velocity_n[ev] - m_matcell_velocity_n[evvois]));
      }
      // On rajoute la dissipation due aux forces dissipatives
      m_dissipation[ev] += -0.5 * math::dot(m_stockes[ev], (m_matcell_velocity[ev] + m_matcell_velocity_n[ev]));
    }
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();   
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      p_ref = options()->environment[index_env].eosModel()->getPressionRef(ev.environment());
      gamF = gruneisen + 1. + p_ref/m_pressure[ev];
      for (Integer index = 0; index < cell.nbFace(); ++index) {
        Face face = cell.face(index);
        // Cell voisine vis-à-vis de cette face 
        Cell cellvois = getCellVoisine(cell, face);
        // index local dans cellvois de la face face
        Integer indexvois = getindiceface(cellvois, face);
        // obtenir l'envcell de cette maille voisine 
        EnvCell evvois = getenvcell(cellvois, index_env);
        
        transport += m_volumic_transfert[ev][index] * m_alpha_rho_n[ev] * m_internal_energy_n[ev] 
        - m_volumic_transfert[evvois][indexvois] * m_alpha_rho_n[evvois] * m_internal_energy_n[evvois];
      }
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev2 = *ienvcell;  
        Integer i_env = ev2.environmentId();   
        gruneisen_phi = options()->environment[i_env].eosModel()->getGruneisen(ev2.environment());
        pphi_ref = options()->environment[i_env].eosModel()->getPressionRef(ev2.environment());
        gamPhi = gruneisen_phi + 1. + p_ref/m_pressure[ev2];
        ech_P += (m_compress[ev2] * m_fracvol[ev] / gamF) * 
            (m_pressure_advection[ev] - m_pressure_advection[ev2]);
        ech_W += (m_compress[ev2]/gamF) * gruneisen * m_dissipation[ev] 
               - (m_compress[ev]/gamPhi) * gruneisen_phi * m_dissipation[ev2];
      }
      m_internal_energy[ev] = (m_internal_energy_n[ev] * m_alpha_rho_n[ev] * m_cell_volume_n[ev] 
                    - deltat * transport 
                    - dt * m_compress[ev] * m_pressure[ev] * m_cell_volume_n[ev] 
                    + dt * ech_P - ech_W + m_dissipation[ev])
                    / (m_alpha_rho[ev] * m_cell_volume[ev]);
    }
  }
}

void MahycoModule::calcVarState() {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  //Calcul de la pression
  double temp_1 = 0.0;
  double temp_2 = 0.0;
  double b, c, delta, gruneisen, p_ref;
  bool error(false);
  Real pref_1(0.);
  Real pref_2(2100000000.0);
  //pression anisentrope
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();  
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      p_ref = options()->environment[index_env].eosModel()->getPressionRef(ev.environment());
      if (ev.environmentId() == 0) 
          temp_1 += gruneisen * m_alpha_rho[ev] * m_internal_energy[ev];
      else 
          temp_2 += gruneisen * m_alpha_rho[ev] * m_internal_energy[ev];
    }
    b = pref_1 + pref_2 - temp_1 - temp_2;
    c = pref_1 * pref_2 - temp_1 * pref_2 - temp_2 * pref_1;
    delta = pow(b, 2) - 4 * c;
    if (delta < 0) {
      delta = 0;
      error = true;
    }
    m_pressure[cell] = (-b + sqrt(delta))*0.5;
    if (error) {
      cerr << " -- discriminant négatif dans la pression maille :" << cell.uniqueId() << endl;
      cerr << "coordonne de la maille= " << m_cell_coord[cell] << endl;
    }
    // Calcul de la masse volumique, de la fraction volumique
    // des coefficients de compression relatifs et de la vitesse du son
    double gamF; 
    temp_1 = temp_2 = 0.0;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      p_ref = options()->environment[index_env].eosModel()->getPressionRef(ev.environment());
      m_density[ev] = (m_pressure[cell] + p_ref) / (gruneisen * m_internal_energy[ev]);
      m_fracvol[ev] = m_alpha_rho[ev] / m_density[ev];
      gamF = gruneisen + 1 + p_ref / m_pressure[cell];
      temp_1 += m_fracvol[ev] / gamF;
    }
    double csF = 0.0;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      Integer index_env = ev.environmentId();
      gruneisen = options()->environment[index_env].eosModel()->getGruneisen(ev.environment());
      p_ref = options()->environment[index_env].eosModel()->getPressionRef(ev.environment());
      gamF = gruneisen + 1 + p_ref / m_pressure[cell];
      m_compress[ev] = m_fracvol[ev] / ( gamF * temp_1);
      csF = sqrt(gamF*m_pressure[cell] / m_density[ev]);
      m_sound_speed[cell] += m_compress[ev] * csF;
    }
  }
}
