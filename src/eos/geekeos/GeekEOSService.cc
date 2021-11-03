// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include "GeekEOSService.h"

using namespace Arcane;
using namespace Arcane::Materials;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void GeekEOSService::initEOS(IMeshEnvironment* env)
{
 
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void GeekEOSService::applyEOS(IMeshEnvironment* env)
{
 
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void GeekEOSService::applyOneCellEOS(IMeshEnvironment* env, EnvCell ev)
{
 
}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
Real GeekEOSService::getAdiabaticCst(IMeshEnvironment* env) { return options()->adiabaticCst();}
Real GeekEOSService::getTensionLimitCst(IMeshEnvironment* env) { return options()->limitTension();}
Real GeekEOSService::getPressionRef(IMeshEnvironment* env) { return options()->pref();}
Real GeekEOSService::getGruneisen(IMeshEnvironment* env) { return options()->gruneisen();}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_SERVICE_PERFECTGASEOS(PerfectGas, PerfectGasEOSService);
