#include <iostream>
#include <vector>
#include <typeinfo>
#include <math.h>
#include <cmath>
#include <string>
#include "MahycoModule.h"

using namespace :: std;

void MahycoModule::calcCSV() {
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real3 somme = {0. , 0. , 0.};
    for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode)
      somme += m_node_coord[inode];
    m_cell_coord[cell] = one_over_nbnode * somme;
  }
  // calcul des volumes
  computeGeometricValues();
  // calcul des normales
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
  // plus petit longeur de maille est dans m_caracteristic_length[icell]
}
