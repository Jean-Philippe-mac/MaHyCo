﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef IREMAPALESERVICE_H
#define IREMAPALESERVICE_H

#include "TypesMahyco.h"

#include "Remap/IRemap.h"
#include "Remap/RemapALE_axl.h"
#include <arcane/ItemTypes.h>
#include "arcane/IMesh.h"
#include "arcane/IItemFamily.h"
#include "arcane/materials/IMeshMaterialMng.h"
#include "arcane/materials/IMeshMaterial.h"
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/IMeshBlock.h"
#include "arcane/materials/MeshMaterialModifier.h"
#include "arcane/materials/MeshMaterialVariableRef.h"
#include "arcane/materials/MeshEnvironmentVariableRef.h"
#include "arcane/materials/MaterialVariableBuildInfo.h"
#include "arcane/materials/MeshBlockBuildInfo.h"
#include "arcane/materials/MeshEnvironmentBuildInfo.h"
#include "arcane/materials/MeshMaterialVariableDependInfo.h"
#include "arcane/materials/CellToAllEnvCellConverter.h"
#include "arcane/materials/MatCellVector.h"
#include "arcane/materials/EnvCellVector.h"
#include "arcane/materials/MatConcurrency.h"
#include "arcane/materials/MeshMaterialIndirectModifier.h"
#include "arcane/materials/MeshMaterialVariableSynchronizerList.h"
#include "arcane/materials/ComponentSimd.h"

#include "arcane/cea/ICartesianMesh.h"
#include "arcane/cea/CellDirectionMng.h"
#include "arcane/cea/FaceDirectionMng.h"
#include "arcane/cea/NodeDirectionMng.h"
#include "arcane/cea/CartesianConnectivity.h"

#include "AcceleratorUtils.h"


using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le service de Remap version ALE
 */
class RemapALEService 
: public ArcaneRemapALEObject
{
public:
  /** Constructeur de la classe */
  RemapALEService(const ServiceBuildInfo & sbi)
    : ArcaneRemapALEObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~RemapALEService() {};
  
  struct interval {
    double inf, sup;
  };
  ICartesianMesh* m_cartesian_mesh;

public:
  
  virtual void initGpu();
  
   /**
   * main du remap
   **/
  virtual void appliRemap(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env);
  /**
   * resize les variables du remap
   **/
  virtual void resizeRemapVariables(Integer nb_vars_to_project, Integer nb_env);
  /**
   * synchronisation des valeurs aux cellules 
   **/
  virtual void synchronizeUremap();

  
  virtual Integer getOrdreProjection();
  virtual bool hasProjectionPenteBorne();
  virtual bool hasConservationEnergieTotale(); 
  virtual bool isEuler();   
  
  virtual void remapVariables(Integer dimension, Integer withDualProjection, Integer nb_vars_to_project, Integer nb_env);
  
private:
    
  void ComputeNodeGroupToRelax();
  void computeLissage();
  void computeVolumes();
  void computeFlux();
  void computeApproPhi(MaterialVariableCellArrayReal, VariableCellArrayReal);
  void computeNewPhi(MaterialVariableCellReal, MaterialVariableCellReal, VariableCellArrayReal);
  
  
  Real m_arithmetic_thresold = 1.e-300;
  
  // Pour l'utilisation des accélérateurs
  ax::Runner m_runner;
  
  
};

#endif
