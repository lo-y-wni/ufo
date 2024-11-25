/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_MPIRANK_H_
#define UFO_FILTERS_OBSFUNCTIONS_MPIRANK_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the MPIRank ObsFunction
class MPIRankParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MPIRankParameters, Parameters)
    // No options.
};

// -----------------------------------------------------------------------------

/// \brief Write MPI rank to the output variable.
///
/// For non-overlapping ObsSpace distributions such as the Round Robin,
/// there is an unambiguous rank that can be assigned to each location.
/// For overlapping ObsSpace distributions such as the InefficientDistribution,
/// the rank associated with the observation 'patch' is assigned to each location.
/// The patch rank is always unique for each location; for example,
/// the patch rank for the InefficientDistribution is always rank 0.
/// Further details can be found in the source code for each distribution.
class MPIRank : public ObsFunctionBase<int> {
 public:
  explicit MPIRank(const eckit::LocalConfiguration &);
  ~MPIRank();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  MPIRankParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_MPIRANK_H_
