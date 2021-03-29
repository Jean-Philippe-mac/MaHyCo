#ifndef TYPESMAHYCO_H
#define TYPESMAHYCO_H

#include <arcane/ItemGroup.h>
#include "eos/IEquationOfState.h"

struct TypesMahyco
{
  enum eBoundaryCondition
  {
    VelocityX, //!< Vitesse X fix�e
    VelocityY, //!< Vitesse Y fix�e
    VelocityZ, //!< Vitesse Z fix�e
    Unknown //!< Type inconnu
  };
};

#endif

