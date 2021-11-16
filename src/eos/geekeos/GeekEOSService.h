﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#ifndef GEEKEOSSERVICE_H
#define GEEKEOSSERVICE_H

#include "eos/IEquationOfState.h"
#include "eos/geekeos/GeekEOS_axl.h"
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
using namespace Arcane;
using namespace Arcane::Materials;

/**
 * Représente le modèle d'équation d'état <em>Geek</em>
 */
class GeekEOSService 
: public ArcaneGeekEOSObject
{
public:
  /** Constructeur de la classe */
  GeekEOSService(const ServiceBuildInfo & sbi)
    : ArcaneGeekEOSObject(sbi) {}
  
  /** Destructeur de la classe */
  virtual ~GeekEOSService() {};

public:
  /** 
   *  Initialise l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et l'énergie interne. 
   */
  virtual void initEOS(IMeshEnvironment* env);
  /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression. 
   */
  virtual void applyEOS(IMeshEnvironment* env);
   /** 
   *  Applique l'équation d'état au groupe de mailles passé en argument
   *  et calcule la vitesse du son et la pression pour une cellule
   */
  virtual void applyOneCellEOS(IMeshEnvironment* env, EnvCell ev);
  /** 
   *  Renvoie la constante adiabatic de l'environnement. 
   */
  virtual Real getAdiabaticCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la constante tension limit de l'environnement. 
   */
  virtual Real getTensionLimitCst(IMeshEnvironment* env);
  /** 
   *  Renvoie la pression de reference
   */
  virtual Real getPressionRef(IMeshEnvironment* env);
  /** 
   *  Renvoie le coefficeint de Gruneisen
   */
  virtual Real getGruneisen(IMeshEnvironment* env);
};

#endif