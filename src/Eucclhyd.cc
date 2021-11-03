#include "MahycoModule.h"

/**
 * Job computeCornerNormal
 * In variables: m_node_coord
 * Out variables: m_lminus, m_lpc, m_lplus, m_nminus, m_nplus
 */
void MahycoModule::computeCornerNormal()  {
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      Node node = cell.node(inode);
      Integer kcell(-1);
      for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
        if (node->cell(jcell) == node) {
          kcell = jcell;
          break;
        }
      }
      Node nodemoins1;
      if (inode != 0) nodemoins1 = cell.node(inode-1);
      else nodemoins1 = cell.node(3);
      Node nodeplus1;
      if (inode != 3) nodeplus1 = cell.node(inode+1);
      else nodeplus1 = cell.node(0);
      Real3 xp = m_node_coord[node];
      Real3 xp_plus = m_node_coord[nodeplus1];
      Real3 xp_moins = m_node_coord[nodemoins1];
      Real3 npc_plus( 0.5 * (xp_plus.y - xp.y) , 0.5 * (xp.x - xp_plus.x) , 0.);
      Real lpc_plus = npc_plus.abs();
      npc_plus = npc_plus / lpc_plus;
      m_nplus[node][kcell] = npc_plus;
      m_lplus[node][kcell] = lpc_plus;
      Real3 npc_minus( 0.5 * (xp.y - xp_moins.y), 0.5 * (xp_moins.x - xp.x) , 0.);
      Real lpc_minus = npc_minus.abs();
      npc_minus = npc_minus / lpc_minus;
      m_nminus[node][kcell] = npc_minus;
      m_lminus[node][kcell] = lpc_minus;
      m_lpc[node][kcell] = (lpc_plus * npc_plus) + (lpc_minus * npc_minus);
    }
  }
}
/**
 * Job computeGradients 
 * In variables: m_node_force_n, m_node_velocity_n, m_lpc, spaceOrder, v
 * Out variables: m_velocity_gradient, m_pressure_gradient,
 * m_pressure_gradient_env
 */
void MahycoModule::computeGradients() {
  
  m_node_fracvol.fill(0.);
  ENUMERATE_ENV(ienv,mm){
    IMeshEnvironment* env = *ienv;
    ENUMERATE_ENVCELL(ienvcell,env){
      EnvCell ev = *ienvcell;
      Cell cell = ev.globalCell();
      for (NodeEnumerator inode(cell.nodes()); inode.hasNext(); ++inode) 
        m_node_fracvol[inode][env->id()] += 0.25 * m_fracvol[ev];
    }
  }
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  if (options()->getOrdreLagrange() == 2)
    // computeGradients
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      Real3x3 somme(Real3(0.0, 0.0, 0.0),
                 Real3(0.0, 0.0, 0.0),
                 Real3(0.0, 0.0, 0.0));
      for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
        somme += math::prodTens(m_velocity[cell.node(inode)], m_cell_cqs_n[cell][inode]);
      }
      m_velocity_gradient[cell] = somme / m_cell_volume_n[cell];
      
      Real3 somme_force(Real3(0.0,0.0,0.0));
      for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
        Node node = cell.node(inode);
        Integer kcell(-1);
        for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
          if (node->cell(jcell) == node) {
            kcell = jcell;
            break;
          }
        }
        somme_force += m_node_force[cell->node(inode)][kcell]; 
      }
      m_pressure_gradient[cell] = somme_force / m_cell_volume_n[cell];
      AllEnvCell all_env_cell = all_env_cell_converter[cell];
      if (all_env_cell.nbEnvironment() !=1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;  
        somme_force = {0.0,0.0,0.0};
        for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
         // somme_force += m_node_force[ev][inode];
        }
        m_pressure_gradient[ev] = somme_force / m_cell_volume_n[cell];
      }
    }
    }
}

/**
 * Job computeDissipationMatrix 
 * In variables: m_sound_speed, m_lminus, m_lplus, m_nminus, m_nplus, m_density_n
 * Out variables: M
 */
void MahycoModule::computeDissipationMatrix() {
  
  computeCellMass();

  ENUMERATE_NODE(node, allNodes()) {
       for (Integer icell = 0; icell < node->nbCell(); ++icell) {
           Cell cell=node->cell(icell);
           Real3x3 cornerMatrix = 
                math::prodTensScal(math::prodTens(m_nplus[node][icell],
                               m_nplus[node][icell]), m_lplus[node][icell])
                +
                math::prodTensScal(math::prodTens(m_nminus[node][icell],
                             m_nminus[node][icell]), m_lminus[node][icell]);

            m_dissipation_matrix[node][icell] =
                m_density_n[cell] * m_sound_speed[cell] * cornerMatrix;

//             m_dissipation_matrix_env(pNodes, cCellsOfNodeP, 0) =
//                 m_density_env_n(cCells)[0] * m_speed_velocity(cCells) *
//                 cornerMatrix;
//             m_dissipation_matrix_env(pNodes, cCellsOfNodeP, 1) =
//                 m_density_env_n(cCells)[1] * m_speed_velocity(cCells) *
//                 cornerMatrix;
//             m_dissipation_matrix_env(pNodes, cCellsOfNodeP, 2) =
//                 m_density_env_n(cCells)[2] * m_speed_velocity(cCells) *
//                 cornerMatrix;
          }
        }
}

/*
 * Job extrapolateValue 
 * In variables: m_cell_velocity_n, , m_cell_coord, m_velocity_gradient,
 * m_pressure_gradient, p, spaceOrder Out variables: m_cell_velocity_extrap,
 * m_pressure_extrap, m_pressure_env_extrap
 */
void MahycoModule::extrapolateValue()  {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  if (options()->getOrdreLagrange() == 1) {
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
        m_cell_velocity_extrap[cell][inode] = m_cell_velocity_n[cell];
        m_pressure_extrap[cell][inode] = m_pressure[cell];
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        if (all_env_cell.nbEnvironment() !=1) {
          ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;  
            m_pressure_extrap[ev][inode] = m_pressure[ev];
          }
        }
      }
    }
  } else {
   ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
     Real min_reduction(MAXFLOAT);
     Real min_reduction_ev(MAXFLOAT);
     for (Integer icell = 0; icell < node->nbCell(); ++icell) {
       Cell cell=node->cell(icell);   
       min_reduction = math::min(min_reduction, m_pressure[cell]);
       AllEnvCell all_env_cell = all_env_cell_converter[cell];
       //if (all_env_cell.nbEnvironment() !=1) {
         ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;  
            min_reduction_ev = math::min(min_reduction_ev, m_pressure[ev]);
         }
       //}
     }
     Real minp(min_reduction);
     Real minp_ev(min_reduction_ev);
      
     Real max_reduction(0.); // MINFLOAT
     Real max_reduction_ev(0.);
     for (Integer icell = 0; icell < node->nbCell(); ++icell) {
       Cell cell=node->cell(icell);   
       max_reduction = math::max(max_reduction, m_pressure[cell]);
       AllEnvCell all_env_cell = all_env_cell_converter[cell];
       //if (all_env_cell.nbEnvironment() !=1) {
         ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
            EnvCell ev = *ienvcell;  
            max_reduction_ev = math::max(max_reduction_ev, m_pressure[ev]);
         }
       //}
      }
      Real maxp(max_reduction);
      Real maxp_ev(max_reduction_ev);
      
      Real3 min_reduction_vec(Real3(0.0, 0.0, 0.0));
      for (Integer icell = 0; icell < node->nbCell(); ++icell) {
        Cell cell=node->cell(icell);   
        min_reduction_vec = math::min(min_reduction_vec, m_cell_velocity_n[cell]);
      }
      Real3 minv(min_reduction_vec);
      
      Real3 max_reduction_vec(Real3(0.0, 0.0, 0.0));
      for (Integer icell = 0; icell < node->nbCell(); ++icell) {
        Cell cell=node->cell(icell);   
        max_reduction_vec = math::min(max_reduction_vec, m_cell_velocity_n[cell]);
      }
      Real3 maxv(max_reduction_vec);
      
      for (Integer icell = 0; icell < node->nbCell(); ++icell) {
        Cell cell=node->cell(icell);  
        Integer knode(-1);
        for (Integer jnode = 0; jnode < cell->nbNode(); ++jnode) {
            if (cell->node(jnode) == node) {
                knode = jnode;
                break;
            }
        }
        AllEnvCell all_env_cell = all_env_cell_converter[cell];
        ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
          EnvCell ev = *ienvcell;  
          Real ptmp_ev = m_pressure[ev] + math::dot(m_pressure_gradient[ev], (m_node_coord[node] - m_cell_coord[cell]));
          m_pressure_extrap[ev][knode] = math::max(math::min(maxp_ev, ptmp_ev),                                             minp_ev);
          m_pressure_extrap[cell][knode] += m_fracvol[ev] * m_pressure_extrap[ev][knode];
         }
         Real3 vtmp = m_cell_velocity_n[cell] 
         + math::prodTensVec(m_velocity_gradient[cell], 
                             (m_node_coord[node] - m_cell_coord[cell]));
         m_cell_velocity_extrap[cell][knode] = math::max(math::min(maxv, vtmp), minv);
      }
    }
  }
}
/**
 * Job computeG.
 * In variables: M, m_cell_velocity_extrap, m_lpc, m_pressure_extrap
 * Out variables: G
 */
void MahycoModule::computeG() {    
   ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
    Real3 reduction(Real3(0.0, 0.0, 0.0));
    for (Integer icell = 0; icell < node->nbCell(); ++icell) {
      Cell cell=node->cell(icell);  
      Integer knode(-1);
      for (Integer jnode = 0; jnode < cell->nbNode(); ++jnode) {
        if (cell->node(jnode) == node) {
          knode = jnode;
          break;
        }
      }
      reduction += math::prodTensVec(m_dissipation_matrix[node][icell],
                                     m_cell_velocity_extrap[cell][knode])
      + m_pressure_extrap[cell][knode] * m_lpc[node][icell];
    
    }
    m_node_G[node] = reduction;
  }
}

/**
 * Job computeNodeDissipationMatrix
 * In variables: M
 * Out variables: m_node_dissipation
 */
void MahycoModule::computeNodeDissipationMatrix() {
  ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
    Real3x3 reduction(Real3(0.0, 0.0, 0.0),
                 Real3(0.0, 0.0, 0.0),
                 Real3(0.0, 0.0, 0.0));
    for (Integer icell = 0; icell < node->nbCell(); ++icell) {
      Cell cell=node->cell(icell);  
      reduction += m_dissipation_matrix[node][icell];
    }
    m_node_dissipation[node] = reduction;
  }
}
/**
 * Job computeNodeVelocity.
 * In variables: G, m_node_dissipation
 * Out variables: m_node_velocity
 */
void MahycoModule::computeNodeVelocity() {
  ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
    if (node.nbCell() != 4) { 
     m_velocity[node] = math::prodTensVec(math::inverseMatrix(m_node_dissipation[node], 1), m_node_G[node]);
    }
  }
}

/**
 * Job computeLagrangePosition.
 * In variables: m_node_velocity_nplus1, m_node_coord, deltat_n
 * Out variables: XLagrange
 */
void MahycoModule::computeLagrangePosition()  {  
  Real deltat = m_global_deltat();
  ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
    m_node_coord[node] += m_velocity[node] * deltat;
  }
//   auto faces(mesh->getFaces());
//   Kokkos::parallel_for(
//       "computeLagrangePosition", nbFaces, KOKKOS_LAMBDA(const int& fFaces) {
//         int fId(faces[fFaces]);
//         int n1FirstNodeOfFaceF(mesh->getFirstNodeOfFace(fId));
//         int n1Id(n1FirstNodeOfFaceF);
//         int n1Nodes(n1Id);
//         int n2SecondNodeOfFaceF(mesh->getSecondNodeOfFace(fId));
//         int n2Id(n2SecondNodeOfFaceF);
//         int n2Nodes(n2Id);
//         RealArray1D<dim> X_face =
//             0.5 * (varlp->XLagrange(n1Nodes) + varlp->XLagrange(n2Nodes));
//         RealArray1D<dim> face_vec =
//             varlp->XLagrange(n2Nodes) - varlp->XLagrange(n1Nodes);
//         varlp->XfLagrange(fFaces) = X_face;
//         varlp->faceLengthLagrange(fFaces) = MathFunctions::norm(face_vec);
//       });
  if (options()->getWithProjection() == 0) {
    ENUMERATE_CELL(icell, allCells()){
      Cell cell = * icell;
      Real reduction(0.); 
      for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
        Node node = cell.node(inode);
        Node nodeplus1;
        if (inode != 3) nodeplus1 = cell.node(inode+1);
        else nodeplus1 = cell.node(0);
        reduction = reduction + (m_node_coord[node] - m_node_coord[nodeplus1]).abs();
      }
      m_cell_perimeter[cell] = reduction;
    }
  }
}

/**
 * Job computeSubCellForce 
 * In variables: M, m_cell_velocity_extrap, m_node_velocity_nplus1, m_lpc, m_pressure_extrap 
 * Out variables: m_node_force_nplus1
 */
void MahycoModule::computeSubCellForce() {
  ENUMERATE_NODE(inode, allNodes()) {
    Node node = *inode;
    for (Integer icell = 0; icell < node->nbCell(); ++icell) {
      Cell cell=node->cell(icell);  
      Integer knode(-1);
      for (Integer jnode = 0; jnode < cell->nbNode(); ++jnode) {
        if (cell->node(jnode) == node) {
          knode = jnode;
          break;
        }
      }
      m_node_force[node][icell] = (- m_pressure_extrap[cell][knode] * m_lpc[node][icell]) + 
          math::prodTensVec(m_dissipation_matrix[node][icell],
                (m_velocity[node]-m_cell_velocity_extrap[cell][knode]));

           /*
            m_node_force_env_nplus1(pNodes, cCellsOfNodeP, 2) =
                (-m_pressure_env_extrap(cCells, pNodesOfCellC)[2] *
                 m_lpc(pNodes, cCellsOfNodeP)) +
                MathFunctions::matVectProduct(
                    m_dissipation_matrix_env(pNodes, cCellsOfNodeP, 2),
                    m_node_velocity_nplus1(pNodes) -
                        m_cell_velocity_extrap(cCells, pNodesOfCellC));                  
           */
      }
    }
}
/**
 * Job computeLagrangeVolumeAndCenterOfGravity called @6.0 in executeTimeLoopN
 * method. In variables: XLagrange Out variables: XcLagrange, vLagrange
 */
void MahycoModule::computeLagrangeVolumeAndCenterOfGravity()  {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    Real reduction(0.0);
    Real3 reduction_vec(Real3(0.0, 0.0, 0.0));
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      Node node = cell.node(inode);
      Integer kcell(-1);
      for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
        if (node->cell(jcell) == node) {
          kcell = jcell;
          break;
        }
      }
      Node nodeplus1;
      if (inode != 3) nodeplus1 = cell.node(inode+1);
      else nodeplus1 = cell.node(0);
      reduction += math::cross2D(m_node_coord[node], m_node_coord[nodeplus1]);
      reduction_vec += math::crossProduct3(m_node_coord[node], m_node_coord[nodeplus1]) * 
        (m_node_coord[node] + m_node_coord[nodeplus1]);
    }
    Real vol = 0.5 * reduction;
    m_cell_volume[cell] = vol;
    m_cell_coord[cell] = (1.0 / (6.0 * vol) * reduction_vec);
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    if (all_env_cell.nbEnvironment() > 1) {
      ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
        EnvCell ev = *ienvcell;  
        m_cell_volume[ev] = m_fracvol[ev] * m_cell_volume[cell];
      }
    }
  }
}
/**
 * Job updateCellCenteredLagrangeVariables
 * method. In variables: m_node_force_nplus1, m_cell_velocity_n,
 * m_node_velocity_nplus1, deltat_n, m_internal_energy_n, m_lpc, m, m_density_n,
 * vLagrange Out variables: ULagrange
 */
void MahycoModule::updateCellCenteredLagrangeVariables()  {
  CellToAllEnvCellConverter all_env_cell_converter(mm);
  Real deltat = m_global_deltat();
  ENUMERATE_CELL(icell, allCells()){
    Cell cell = * icell;
    AllEnvCell all_env_cell = all_env_cell_converter[cell];
    // 
    // densite
    // 
    Real reduction(0.0);
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      Node node = cell.node(inode);
      Integer kcell(-1);
      for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
        if (node->cell(jcell) == node) {
          kcell = jcell;
          break;
        }
      }
      reduction += math::dot(m_lpc[node][kcell], m_velocity[node]);
    }
    Real rhoLagrange = 1. / (1. / m_density_n[cell] +
                 deltat / m_cell_mass[cell] * reduction);
    m_density[cell] = rhoLagrange;

    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      if (m_fracvol[ev] > options()->getThreshold()) 
      m_density[ev] = m_mass_fraction[ev] *  rhoLagrange / m_fracvol[ev];
    }

    // 
    // Vitesse
    // 
    Real3 reduction_vec(Real3(0.0, 0.0, 0.0));
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      Node node = cell.node(inode);
      Integer kcell(-1);
      for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
        if (node->cell(jcell) == node) {
          kcell = jcell;
          break;
        }
      }
      reduction_vec += m_node_force[node][kcell];
    }
    m_cell_velocity[cell] = m_cell_velocity_n[cell] 
            + reduction_vec * deltat / m_cell_mass[cell];
    // 
    // Energie interne 
    // 
    Real reduction2(0.0);
    for (Integer inode = 0; inode < cell.nbNode(); ++inode) {
      Node node = cell.node(inode);
      Integer kcell(-1);
      for (Integer jcell = 0; jcell < node->nbCell(); ++jcell) {
        if (node->cell(jcell) == node) {
          kcell = jcell;
          break;
        }
      }
      reduction2 += math::dot(m_node_force[node][kcell], 
       m_velocity[node] - 0.5*(m_cell_velocity_n[cell] + m_cell_velocity[cell]));
    } 
    Real eLagrange = m_internal_energy_n[cell] +
            deltat / m_cell_mass[cell] * reduction2;
    ENUMERATE_CELL_ENVCELL(ienvcell,all_env_cell) {
      EnvCell ev = *ienvcell;  
      m_internal_energy[ev] = m_internal_energy[ev] + m_fracvol[ev] 
              * deltat / m_cell_mass[ev] * reduction2;
      m_internal_energy[cell] +=  m_mass_fraction[ev] * m_internal_energy[ev];
    }

  }
}

/* A TRAITER PLUS TARD 
        // m_total_energy_L(cCells) = (rhoLagrange * vLagrange(cCells)) *
        // (eLagrange + 0.5 * (cell_velocity_L[0] * cell_velocity_L[0] +
        // cell_velocity_L[1] * cell_velocity_L[1]));
        m_global_masse_L(cCells) = 0.;
        m_total_energy_L(cCells) =
            (rhoLagrange * varlp->vLagrange(cCells)) *
            (m_mass_fraction_env(cCells)[0] * peLagrange[0] +
             m_mass_fraction_env(cCells)[1] * peLagrange[1] +
             m_mass_fraction_env(cCells)[2] * peLagrange[2] +
             0.5 * (cell_velocity_L[0] * cell_velocity_L[0] +
                    cell_velocity_L[1] * cell_velocity_L[1]));
        for (int imat = 0; imat < nbmat; imat++) {
          m_global_masse_L(cCells) += m_mass_fraction_env(cCells)[imat] *
                                      (rhoLagrange * varlp->vLagrange(cCells));
        }
      });
  double reductionE(0.), reductionM(0.);
  {
    Kokkos::Sum<double> reducerE(reductionE);
    Kokkos::parallel_reduce("reductionE", nbCells,
                            KOKKOS_LAMBDA(const int& cCells, double& x) {
                              reducerE.join(x, m_total_energy_L(cCells));
                            },
                            reducerE);
    Kokkos::Sum<double> reducerM(reductionM);
    Kokkos::parallel_reduce("reductionM", nbCells,
                            KOKKOS_LAMBDA(const int& cCells, double& x) {
                              reducerM.join(x, m_global_masse_L(cCells));
                            },
                            reducerM);
  }
  m_global_total_energy_L = reductionE;
  m_total_masse_L = reductionM;
}
void Eucclhyd::switchalpharho_rho() noexcept {
  Kokkos::parallel_for(
      "updateParticleCoefficient", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        m_density_n(cCells) /=
            particules->m_cell_particle_volume_fraction(cCells);
        for (int imat = 0; imat < nbmatmax; imat++)
          m_density_env_n(cCells)[imat] /=
              particules->m_cell_particle_volume_fraction(cCells);
      });
}

void Eucclhyd::switchrho_alpharho() noexcept {
  Kokkos::parallel_for(
      "updateParticleCoefficient", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        m_density_n(cCells) *=
            particules->m_cell_particle_volume_fraction(cCells);
        for (int imat = 0; imat < nbmatmax; imat++)
          m_density_env_n(cCells)[imat] *=
              particules->m_cell_particle_volume_fraction(cCells);
      });
}
void Eucclhyd::PreparecellvariablesForParticles() noexcept {
  Kokkos::parallel_for(
      "copycellvariables", nbCells, KOKKOS_LAMBDA(const int& cCells) {
        for (int imat = 0; imat < nbmatmax; imat++) {
          particules->m_particlecell_fracvol_env(cCells)[imat] =
              m_fracvol_env(cCells)[imat];
          particules->m_particlecell_density_env_n(cCells)[imat] =
              m_density_env_n(cCells)[imat];
        }
        particules->m_particlecell_density_n(cCells) = m_density_n(cCells);
        particules->m_particlecell_euler_volume(cCells) =
            m_euler_volume(cCells);
        particules->m_particlecell_velocity_n(cCells) =
            m_cell_velocity_n(cCells);
        particules->m_particlecell_velocity_nplus1(cCells) =
            m_cell_velocity_nplus1(cCells);
        particules->m_particlecell_mass(cCells) = m_cell_mass(cCells);
      });
  std::cout << " fin prepare " << std::endl;
}
*/
