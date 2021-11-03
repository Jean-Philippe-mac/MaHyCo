#include "POINTTRIPLEService.h"


void POINTTRIPLEService::initMatMono(Integer dim)  {
  pinfo() << " Initialisation du materiau du POINT-TRIPLE MonoMat";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    m_materiau[cell] = 0;      
  }
}
void POINTTRIPLEService::initMat(Integer dim)  {
  if (options()->casTest == MonoTriplePoint) {
    initMatMono(dim);
    return;
  }
  pinfo() << " Initialisation dzs materiaux du POINT-TRIPLE Multi-Mat";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;

    double x = m_cell_coord[cell].x;
    double y = m_cell_coord[cell].y;
    if (x <= 0.1) {
      m_materiau[cell] = 0;
    } else {
        if (y <= 0.15) {
            m_materiau[cell] = 1;
        } else {
            m_materiau[cell] = 2;
        }
    }
  }
}
void POINTTRIPLEService::initVarMono(Integer dim)  {
  pinfo() << " Initialisation des donnees  du POINT-TRIPLE MonoMat";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    // pseudo-viscosité 
    m_pseudo_viscosity[cell] = 0.;
    // que des mailles pures à l'init
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
    
    // vitesse eucclhy
    m_cell_velocity[cell] = {0.0, 0.0, 0.0};   
    
    // parametres maille
    double x = m_cell_coord[cell].x;
    double y = m_cell_coord[cell].y;
    if (x <= 0.1) {
        m_density[cell] = 1.0;
        m_pressure[cell] = 1.0;
    } else {
        if (y <= 0.15) {
            m_density[cell] = 1.0;
            m_pressure[cell] = 0.1;
        } else {
            m_density[cell] = 0.1;
            m_pressure[cell] = 0.1;
        }
    }
  }
  pinfo() << " fin de la boucle sur les mailles ";
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    //m_velocity[inode].y = 1.;
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
  pinfo() << " fin de la boucle sur les noeuds ";
}
void POINTTRIPLEService::initVar(Integer dim)  {
  if (options()->casTest == MonoTriplePoint) {
    initVarMono(dim);
    return;
  }
  pinfo() << " Initialisation des donnees  du POINT-TRIPLE Multi-Mat";
  ENUMERATE_CELL(icell,allCells()) {
    Cell cell = *icell;
    
    // vitesse eucclhy
    m_cell_velocity[cell] = {0.0, 0.0, 0.0};   
    // pseudo-viscosité 
    m_pseudo_viscosity[cell] = 0.;
    // que des mailles pures à l'init
    m_fracvol[cell] = 1.;
    m_mass_fraction[cell] = 1.;
    // parametres maille
    double x = m_cell_coord[cell].x;
    double y = m_cell_coord[cell].y;
    if (x <= 0.1) {
        m_density[cell] = 1.0;
        m_pressure[cell] = 1.0;
    } else {
        if (y <= 0.15) {
            m_density[cell] = 1.0;
            m_pressure[cell] = 0.1;
        } else {
            m_density[cell] = 0.1;
            m_pressure[cell] = 0.1;
        }
    }
  }
  pinfo() << " fin de la boucle sur les mailles ";
  ENUMERATE_NODE(inode, allNodes()){
    m_velocity[inode] = {0.0, 0.0, 0.0};    
    //m_velocity[inode].y = 1.;
    // sauvegarde des valeurs initiales mises dans m_velocity_n
    m_velocity_n[inode] = m_velocity[inode];
  }
  pinfo() << " fin de la boucle sur les noeuds ";
}

/*---------------------------------------------------------------------------*/

bool POINTTRIPLEService::hasReverseOption() { return options()->reverseOption;}
Real POINTTRIPLEService::getReverseParameter() { return options()->parametre;}

ARCANE_REGISTER_SERVICE_POINTTRIPLE(POINTTRIPLE, POINTTRIPLEService); 
