/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(AMRREDUCTIONS_HPP)
#error "This file should only be included through AMRReductions.hpp"
#endif

#ifndef AMRREDUCTIONS_IMPL_HPP
#define AMRREDUCTIONS_IMPL_HPP

template <VariableType var_t>
AMRReductions<var_t>::AMRReductions(const GRAMR &a_gramr,
                                    const int a_base_level)
    : m_base_level(a_base_level),
      m_coarsest_dx(a_gramr.get_gramrlevels()[0]->get_dx())
{
    set_level_data_vect(a_gramr);
    set_ref_ratios_vect(a_gramr);
    set_domain_volume();
}

template <VariableType var_t>
void AMRReductions<var_t>::set_level_data_vect(const GRAMR &a_gramr)
{
    // this function already checks if the level pointers are null
    const auto gramrlevel_ptrs = a_gramr.get_gramrlevels();
    int num_levels = gramrlevel_ptrs.size();
    m_level_data_ptrs.resize(num_levels);

    for (int ilev = 0; ilev < num_levels; ++ilev)
    {
        m_level_data_ptrs[ilev] = const_cast<GRLevelData *>(
            &gramrlevel_ptrs[ilev]->getLevelData(var_t));
    }
}

template <VariableType var_t>
void AMRReductions<var_t>::set_ref_ratios_vect(const GRAMR &a_gramr)
{
    // this function already checks if the level pointers are null
    const auto gramrlevel_ptrs = a_gramr.get_gramrlevels();
    int num_levels = gramrlevel_ptrs.size();
    m_ref_ratios.resize(num_levels);

    for (int ilev = 0; ilev < num_levels; ++ilev)
    {
        m_ref_ratios[ilev] = gramrlevel_ptrs[ilev]->refRatio();
    }
}

template <VariableType var_t> void AMRReductions<var_t>::set_domain_volume()
{
    // first check if m_level_data_ptrs has been set
    CH_assert((m_level_data_ptrs.size() > 0) &&
              (m_level_data_ptrs[0] != nullptr));

    // first calculate the volume assuming each cell on the coarsest level has
    // unit length
    int cell_volume =
        m_level_data_ptrs[0]->disjointBoxLayout().physDomain().size().product();

    // multiply by dx_coarsest to get real volume
    m_domain_volume =
        pow(m_coarsest_dx, CH_SPACEDIM) * static_cast<double>(cell_volume);
}

template <VariableType var_t>
Real AMRReductions<var_t>::min(const Interval &a_vars) const
{
    CH_assert(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    CH_TIME("AMRReductions::min");
    return computeMin(m_level_data_ptrs, m_ref_ratios, a_vars, m_base_level);
}

template <VariableType var_t>
Real AMRReductions<var_t>::min(const int a_var) const
{
    return min(Interval(a_var, a_var));
}

template <VariableType var_t>
Real AMRReductions<var_t>::max(const Interval &a_vars) const
{
    pout() << a_vars.begin() << " : " << a_vars.end() << " : " << m_num_vars << std::endl;
    CH_assert(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    CH_TIME("AMRReductions::max");
    return computeMax(m_level_data_ptrs, m_ref_ratios, a_vars, m_base_level);
}

template <VariableType var_t>
Real AMRReductions<var_t>::max(const int a_var) const
{
    return max(Interval(a_var, a_var));
}

template <VariableType var_t>
Real AMRReductions<var_t>::norm(const Interval &a_vars,
                                const int a_norm_exponent,
                                const bool a_normalize_by_volume) const
{
    CH_assert(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    CH_TIME("AMRReductions::norm");
    Real norm = computeNorm(m_level_data_ptrs, m_ref_ratios, m_coarsest_dx,
                            a_vars, a_norm_exponent, m_base_level);
    if (a_normalize_by_volume)
    {
        norm /=
            pow(m_domain_volume, 1.0 / static_cast<double>(a_norm_exponent));
    }

    return norm;
}

template <VariableType var_t>
Real AMRReductions<var_t>::norm(const int a_var, const int a_norm_exponent,
                                const bool a_normalize_by_volume) const
{
    return norm(Interval(a_var, a_var), a_norm_exponent, a_normalize_by_volume);
}

template <VariableType var_t>
Real AMRReductions<var_t>::sum(const Interval &a_vars) const
{
    CH_assert(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    CH_TIME("AMRReductions::sum");
    return computeSum(m_level_data_ptrs, m_ref_ratios, m_coarsest_dx, a_vars,
                      m_base_level);
}

template <VariableType var_t>
Real AMRReductions<var_t>::sum(const int a_var) const
{
    return sum(Interval(a_var, a_var));
}

template <VariableType var_t>
Real AMRReductions<var_t>::get_domain_volume() const
{
    return m_domain_volume;
}

template <VariableType var_t>
RealVect AMRReductions<var_t>::maxIndex(const Interval &a_vars) const
{
    CH_assert(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    CH_TIME("AMRReductions::max");
    return computeMaxIndex(m_level_data_ptrs, m_ref_ratios, a_vars, m_base_level);
}

template <VariableType var_t>
RealVect AMRReductions<var_t>::maxIndex(const int a_var) const
{
    return maxIndex(Interval(a_var, a_var));
}

template <VariableType var_t>
RealVect AMRReductions<var_t>::computeMaxIndex(
    const Vector<LevelData<FArrayBox> *> &a_phi, const Vector<int> &a_nRefFine,
    const Interval a_comps, const int a_lBase) const {
  int numLevels = a_phi.size();

  // it is often the case that while a_phi has many possible
  // levels, only a subset of them are defined -- check that
  // just to be sure
  if (a_phi[numLevels - 1] == NULL) {
    int lev = numLevels - 1;

    while (a_phi[lev] == NULL) {
      lev--;
    }

    numLevels = lev + 1;
  }

  Real max;
  Real maxOnLevel;

  // Initializing these to -1 because indices are always > 0
  IntVect maxIndex = -1 * IntVect::Unit;
  IntVect maxIndexOnLevel = -1 * IntVect::Unit;

  RealVect max_L(-1 * IntVect::Unit);

  max = -HUGE_VAL;

  // loop over levels
  for (int lev = a_lBase; lev < numLevels; lev++) {
    // in case there are extra levels which are not defined
    if (a_phi[lev] != NULL) {
      CH_assert(a_phi[lev]->isDefined());

      LevelData<FArrayBox> &thisPhi = *(a_phi[lev]);
      const DisjointBoxLayout *finerGridsPtr = NULL;

      if (lev < numLevels - 1) {
        finerGridsPtr = &(a_phi[lev + 1]->getBoxes());

        maxOnLevel =
            computeMax(thisPhi, finerGridsPtr, a_nRefFine[lev], a_comps);
        maxIndexOnLevel =
            computeMaxIndex(thisPhi, finerGridsPtr, a_nRefFine[lev], a_comps);
      }
      // Not clear to me what this else is doing, but cba to find out
      else {
        int bogusRefRatio = 100000;
        maxOnLevel = computeMax(thisPhi, finerGridsPtr, bogusRefRatio, a_comps);
        maxIndexOnLevel =
            computeMaxIndex(thisPhi, finerGridsPtr, bogusRefRatio, a_comps);
      }

      if (maxOnLevel > max) {
        max = maxOnLevel;
        maxIndex = maxIndexOnLevel;
        // Below, the IntVec maxIndex turns into the RealVect max_L, which
        // gives length instead of gridpoint. The + .5 is because we like
        // to be cell-centered
        max_L = (RealVect(maxIndex) + .5) * m_coarsest_dx * pow(.5, lev);
      }
    }
  }

  // shouldn't need to do broadcast/gather thing
  return max_L;
}
template <VariableType var_t>
IntVect AMRReductions<var_t>::computeMaxIndex(
    const LevelData<FArrayBox> &a_phi, const DisjointBoxLayout *a_finerGridsPtr,
    const int a_nRefFine, const Interval a_comps) const {
  // we need this struct to accommodate the MPI_Reduce with MPI_MAXLOC later
  struct {
    Real levelMax;
    int rank;
  } levelMaxStruct, recv;

  Real thisMax;

  levelMaxStruct.levelMax = -99999999999.9;
  int my_rank = 0;
#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &my_rank);
#endif
  levelMaxStruct.rank = my_rank;

  IntVect LevelmaxIndex = -1 * IntVect::Unit;
  IntVect thismaxIndex = -1 * IntVect::Unit;

  const DisjointBoxLayout &levelGrids = a_phi.getBoxes();

  LevelData<FArrayBox> temp(levelGrids, a_comps.size());

  Interval tempComps(0, a_comps.size() - 1);

  DataIterator dit = a_phi.dataIterator();

  for (dit.reset(); dit.ok(); ++dit) {
    temp[dit()].copy(temp[dit()].box(), tempComps, temp[dit()].box(),
                     a_phi[dit()], a_comps);

    if (a_finerGridsPtr != NULL) {
      LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

      // now loop over fine boxes and set covered regions to 0
      // I don't understand why we do this but not in the mood
      // to find out rn
      for (litFine.reset(); litFine.ok(); ++litFine) {
        Box coveredBox(a_finerGridsPtr->get(litFine()));
        coveredBox.coarsen(a_nRefFine);
        coveredBox &= temp[dit()].box();

        if (!coveredBox.isEmpty()) {
          temp[dit()].setVal(levelMaxStruct.levelMax, coveredBox, 0,
                             tempComps.size());
        }
      } // end loop over fine-grid boxes
    }   // end if there is a finer level

    // while we're looping over the grids, get max as well
    // need to loop over comps
    for (int comp = tempComps.begin(); comp <= tempComps.end(); ++comp) {
      thisMax = temp[dit()].max(comp);
      thismaxIndex = temp[dit()].maxIndex(comp);
      if (thisMax > levelMaxStruct.levelMax) {
        levelMaxStruct.levelMax = thisMax;
        LevelmaxIndex = thismaxIndex;
      }
    }
  } // end loop over this level's grids

  // This is where we do an MPI_Reduce. The MPI_MAXLOC operation will tell us
  // where the maximum value of the variable is, so we can query its location
  // from that specific rank
#ifdef CH_MPI
  int root_rank = 0; // this is where we send the maximum's location
  int result = MPI_Reduce(&levelMaxStruct, &recv, 1, MPI_DOUBLE_INT, MPI_MAXLOC,
                          root_rank, Chombo_MPI::comm);

  if (result != MPI_SUCCESS) {
    MayDay::Error("Sorry, but I had a communication error in AMRReductions");
  }
  int rank_with_max = recv.rank;

  // The root rank must let the other ranks know on which rank the max value is
  result = MPI_Bcast(&rank_with_max, 1, MPI_INT, root_rank, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("Sorry, but I had a communication error in AMRReductions");
  }

  int LevelmaxIndex_x = LevelmaxIndex[0];
  int LevelmaxIndex_y = LevelmaxIndex[1];
#if CH_SPACEDIM == 3
  int LevelmaxIndex_z = LevelmaxIndex[2];
#endif

  // Now we can send the maximum's location from rank_with_max to the other
  // ranks. I'm doing one Bcast per coordinate for now, I'm sure you could
  // define a custom MPI type but no
  result =
      MPI_Bcast(&LevelmaxIndex_x, 1, MPI_INT, rank_with_max, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("Sorry, but I had a communication error in AMRReductions");
  }
  result =
      MPI_Bcast(&LevelmaxIndex_y, 1, MPI_INT, rank_with_max, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("Sorry, but I had a communication error in AMRReductions");
  }

  // Reconstruct the RealVect and done
  LevelmaxIndex[0] = LevelmaxIndex_x;
  LevelmaxIndex[1] = LevelmaxIndex_y;

  // Want to do the same for z, but only if SpaceDim == 3
#if CH_SPACEDIM == 3
  result =
      MPI_Bcast(&LevelmaxIndex_z, 1, MPI_INT, rank_with_max, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("Sorry, but I had a communication error in AMRReductions");
  }
  LevelmaxIndex[2] = LevelmaxIndex_z;
#endif
#endif // CH_MPI

  return LevelmaxIndex;
}

#endif
