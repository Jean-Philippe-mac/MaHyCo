#ifndef ACC_ENV_MULTI_ENV_UTILS_H
#define ACC_ENV_MULTI_ENV_UTILS_H

#include "accenv/AcceleratorUtils.h"

#include "arcane/MeshVariableScalarRef.h"
#include "arcane/MeshVariableArrayRef.h"
#include "arcane/accelerator/Views.h"
#include <arcane/IMesh.h>
#include <arcane/VariableBuildInfo.h>
#include <arcane/VariableView.h>

#include <arcane/materials/ComponentPartItemVectorView.h>
#include "arcane/materials/IMeshEnvironment.h"
#include "arcane/materials/MeshMaterialVariableRef.h"

/*---------------------------------------------------------------------------*/
/* Pour créer une vue sur les valeurs d'un environnement                     */
/*---------------------------------------------------------------------------*/
using namespace Arcane;
using namespace Arcane::Materials;

template<typename value_type>
ArrayView<value_type> envView(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshEnvironment* env) {
  // Pour rappel, [0] fait référence à .globalVariable() d'où le +1
  return var_menv._internalValue()[env->id()+1];
}

/*---------------------------------------------------------------------------*/
/* Pour repérer un EnvCell (fortement inspiré de MatVarIndex)                */
/*---------------------------------------------------------------------------*/
class EnvVarIndex {
 public:
  ARCCORE_HOST_DEVICE EnvVarIndex(Int32 array_index, Int32 value_index) :
    m_array_index (array_index),
    m_value_index (value_index)
  {
  }
  ARCCORE_HOST_DEVICE Int32 arrayIndex() const { return m_array_index; }
  ARCCORE_HOST_DEVICE Int32 valueIndex() const { return m_value_index; }
 protected:
  Int32 m_array_index;
  Int32 m_value_index;
};

/*---------------------------------------------------------------------------*/
/* Vues sur les variables multi-environnement                                */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvView {
 public:
  MultiEnvView(const Span< Span<value_type> >& spn) :
    m_var_menv_views (spn)
  {
  }

  ARCCORE_HOST_DEVICE MultiEnvView(const MultiEnvView<value_type>& rhs) :
    m_var_menv_views (rhs.m_var_menv_views)
  {
  }

  ARCCORE_HOST_DEVICE value_type operator[](const EnvVarIndex& evi) const {
    return m_var_menv_views[evi.arrayIndex()][evi.valueIndex()];
  }

  ARCCORE_HOST_DEVICE void setValue(const EnvVarIndex& evi, value_type val) const {
    return m_var_menv_views[evi.arrayIndex()].setItem(evi.valueIndex(), val);
  }

 public:
  Span< Span<value_type> > m_var_menv_views;
};

/*---------------------------------------------------------------------------*/
/* Pour créer des vues sur les variables multi-environnement                 */
/*---------------------------------------------------------------------------*/
template<typename value_type>
class MultiEnvVar {
 public:
  MultiEnvVar(CellMaterialVariableScalarRef<value_type>& var_menv, IMeshMaterialMng* mm) :
   m_var_menv_impl(platform::getAcceleratorHostMemoryAllocator(), mm->environments().size()+1)
  {
    m_var_menv_impl[0].setArray(Span<value_type>(var_menv._internalValue()[0]));
    ENUMERATE_ENV(ienv, mm) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();
      m_var_menv_impl[env_id+1].setArray(Span<value_type>(envView(var_menv,env)));
    }
  }

  //! Pour y accéder en lecture/écriture
  auto span() {
    return MultiEnvView<value_type>(m_var_menv_impl);
  }

 protected:
  UniqueArray< Span<value_type> > m_var_menv_impl;
};

/*---------------------------------------------------------------------------*/
/* Vue sur stockage du multi-env                                             */
/*---------------------------------------------------------------------------*/
class MultiEnvCellViewIn {
 public:
  MultiEnvCellViewIn(ax::RunCommand& command, Integer max_nb_env,
      const VariableCellInteger& v_nb_env,
      const UniqueArray<Int16>& v_l_env_arrays_idx,
      const VariableCellArrayInteger& v_l_env_values_idx) :
    m_max_nb_env (max_nb_env),
    m_nb_env_in (ax::viewIn(command,v_nb_env)),
    m_l_env_arrays_idx_in (v_l_env_arrays_idx.constSpan()),
    m_l_env_values_idx_in (ax::viewIn(command,v_l_env_values_idx))
  {}

  ARCCORE_HOST_DEVICE MultiEnvCellViewIn(const MultiEnvCellViewIn& rhs) : 
    m_max_nb_env (rhs.m_max_nb_env),
    m_nb_env_in (rhs.m_nb_env_in),
    m_l_env_arrays_idx_in (rhs.m_l_env_arrays_idx_in),
    m_l_env_values_idx_in (rhs.m_l_env_values_idx_in)
  {}

  ARCCORE_HOST_DEVICE Integer nbEnv(CellLocalId cid) const {
    return m_nb_env_in[cid];
  }

  ARCCORE_HOST_DEVICE EnvVarIndex envCell(CellLocalId cid, Integer ienv) const {
    return EnvVarIndex(
        m_l_env_arrays_idx_in[cid.localId()*m_max_nb_env+ienv],
        m_l_env_values_idx_in[cid][ienv]);
  }

 protected:
  Integer m_max_nb_env;
  ax::VariableCellInt32InView m_nb_env_in;
  Span<const Int16> m_l_env_arrays_idx_in;
  ax::ItemVariableArrayInViewT<Cell,Int32> m_l_env_values_idx_in;
};

/*---------------------------------------------------------------------------*/
/* Stockage du multi-env                                                     */
/*---------------------------------------------------------------------------*/
class MultiEnvCellStorage {
 public:
  MultiEnvCellStorage(IMeshMaterialMng* mm, AccMemAdviser* acc_mem_adv) :
    m_mesh_material_mng (mm),
    m_max_nb_env (mm->environments().size()),
    m_nb_env(VariableBuildInfo(mm->mesh(), "NbEnv" , IVariable::PNoDump| IVariable::PNoNeedSync)),
    m_l_env_arrays_idx(platform::getAcceleratorHostMemoryAllocator()),
    m_l_env_values_idx(VariableBuildInfo(mm->mesh(), "LEnvValuesIdx" , IVariable::PNoDump| IVariable::PNoNeedSync))
#ifdef ARCANE_DEBUG
    , m_env_id(VariableBuildInfo(mm->mesh(), "EnvId" , IVariable::PNoDump| IVariable::PNoNeedSync))
#endif
  {
    m_l_env_arrays_idx.resize(m_max_nb_env*mm->mesh()->allCells().size());
    acc_mem_adv->setReadMostly(m_l_env_arrays_idx.view());
    m_l_env_values_idx.resize(m_max_nb_env);
  }

  //! Remplissage
  void buildStorage(Materials::MaterialVariableCellInteger& v_global_cell) {
    ENUMERATE_CELL(icell, m_mesh_material_mng->mesh()->allCells()){
      // Init du nb d'env par maille qui va être calculé à la boucle suivante
      m_nb_env[icell] = 0;
    }

    auto in_nb_env  = Arcane::viewIn(m_nb_env);
    auto out_nb_env = Arcane::viewOut(m_nb_env);
    auto out_l_env_pos = Arcane::viewOut(m_l_env_values_idx);

    ENUMERATE_ENV(ienv, m_mesh_material_mng) {
      IMeshEnvironment* env = *ienv;
      Integer env_id = env->id();

      // Mailles mixtes
      Span<const Integer> in_global_cell(envView(v_global_cell, env));

      Integer nb_imp = env->impureEnvItems().nbItem();
      for(Integer imix = 0 ; imix < nb_imp ; ++imix) {
        CellLocalId cid(in_global_cell[imix]); // on récupère l'identifiant de la maille globale

        auto index_cell = in_nb_env[cid];

        // On relève le numéro de l'environnement 
        // et l'indice de la maille dans la liste de mailles mixtes env
        m_l_env_arrays_idx[cid*m_max_nb_env+index_cell] = env_id+1; // décalage +1 car 0 est pris pour global
        out_l_env_pos[cid][index_cell] = imix;

        out_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
      }

      // Mailles pures
      // Pour les mailles pures, valueIndexes() est la liste des ids locaux des mailles
      Span<const Int32> in_cell_id(env->pureEnvItems().valueIndexes());

      // Nombre de mailles pures de l'environnement
      Integer nb_pur = env->pureEnvItems().nbItem();

      for(Integer ipur = 0 ; ipur < nb_pur ; ++ipur) {
        CellLocalId cid(in_cell_id[ipur]); // accés indirect à la valeur de la maille

        auto index_cell = in_nb_env[cid];
        ARCANE_ASSERT(index_cell==0, ("Maille pure mais index_cell!=0"));

        m_l_env_arrays_idx[cid*m_max_nb_env+index_cell] = 0; // 0 référence le tableau global
        out_l_env_pos[cid][index_cell] = cid.localId();

        // Equivalent à affecter 1
        out_nb_env[cid] = index_cell+1; // ++ n'est pas supporté
      }
    }

    checkStorage(v_global_cell);
  }

  //! Verification
  void checkStorage(Materials::MaterialVariableCellInteger& v_global_cell) {
#ifdef ARCANE_DEBUG
    // m_env_id doit être calculé
    MultiEnvVar<Integer> menv_global_cell(v_global_cell, m_mesh_material_mng);
    auto in_menv_global_cell(menv_global_cell.span());

    ENUMERATE_CELL(icell, m_mesh_material_mng->mesh()->allCells()){
      Cell cell = * icell;
      Integer cid = cell.localId();

      Integer nb_env = m_nb_env[icell];
      Integer nb_env_bis = (m_env_id[cell]>=0 ? 1 : -m_env_id[cell]-1);
      ARCANE_ASSERT(nb_env_bis==nb_env, ("nb_env_bis!=nb_env"));

      for(Integer ienv=0 ; ienv<nb_env ; ++ienv) {
        Integer i=m_l_env_arrays_idx[cid*m_max_nb_env+ienv];
        Integer j=m_l_env_values_idx[icell][ienv];
        EnvVarIndex evi(i,j);

        Integer cid_bis=in_menv_global_cell[evi];
        ARCANE_ASSERT(cid_bis==cid, ("cid_bis!=cid"));
      }
    }
#endif
  }

  //! Vue sur le stockage pour utilisation en lecture sur GPU
  MultiEnvCellViewIn viewIn(ax::RunCommand& command) {
    return MultiEnvCellViewIn(command, m_max_nb_env, m_nb_env, m_l_env_arrays_idx, m_l_env_values_idx);
  }

 protected:
  IMeshMaterialMng* m_mesh_material_mng=nullptr;
  Integer m_max_nb_env;
  VariableCellInteger m_nb_env;  //! Nb d'env par maille
  UniqueArray<Int16> m_l_env_arrays_idx; //! liste des indexes des env par maille
  VariableCellArrayInteger m_l_env_values_idx;
#ifdef ARCANE_DEBUG
  VariableCellInteger m_env_id;
#endif
};

#endif

