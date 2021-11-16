#include "RIDERGEEKService.h"


void RIDERGEEKService::initMatMono(Integer dim)  {

}
void RIDERGEEKService::initMat(Integer dim)  {

  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
      // toutes les mailles ont les 2 materiaux
      m_materiau[cell] = 0.5;
  }
}
void RIDERGEEKService::initVarMono(Integer dim)  {
 
}
void RIDERGEEKService::initVar(Integer dim)  {
  CellToAllEnvCellConverter all_env_cell_converter(IMeshMaterialMng::getReference(mesh()));
  //caracteristiques du tests
  Real ampli = 0.1;
  Real matelas = 0.1;
  Real P_atm = 100000.0;
  Real r_0 = 1.0;
  Real r_1 = 1000;
  Real3 u_1(0.0, 0.0, 0.0);
  if ( options()->casTest == RiderTx) u_1.x = 100.;
  if ( options()->casTest == RiderTy) u_1.y = 100.;
  if ( options()->casTest == RiderT45) u_1.x = u_1.y = 100.;
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell]; 
    m_cell_gaussien1[cell] =   ampli * math::exp(-math::pow(m_cell_coord[cell].x - 0.25, 2.0) + math::pow(m_cell_coord[cell].y - 0.25, 2.0) / 0.005);
    m_cell_gaussien2[cell] =   ampli * math::exp(-math::pow(m_cell_coord[cell].x - 0.25, 2.0) + math::pow(m_cell_coord[cell].y - 0.25, 2.0) / 0.005);
    // pseudo-viscosité 
    m_pseudo_viscosity[cell] = 0.;
    m_density[cell] = 0.;
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
    Real fv0 = (1.0 - m_cell_gaussien1[cell] - matelas);
    Real fv1 = (m_cell_gaussien1[cell] - matelas);
    if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;  
        Integer index_env = ev.environmentId();  
        m_fracvol[ev] = fv0 * (1-index_env) + fv1 * index_env ;
        m_pressure[ev] = P_atm;
        m_density[ev] = r_0 * (1-index_env) + r_1 * index_env;
        // vitesse par materiaux
        m_matcell_velocity[ev] = - ( fv1 * u_1 / fv0) * (1-index_env) + u_1 * index_env;
        // densitee moyenne
        m_density[cell] +=  m_fracvol[ev] * m_density[ev];
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

bool RIDERGEEKService::hasReverseOption() { return options()->reverseOption;}
Real RIDERGEEKService::getReverseParameter() { return options()->parametre;}

ARCANE_REGISTER_SERVICE_RIDERGEEK(RIDERGEEK, RIDERGEEKService);
