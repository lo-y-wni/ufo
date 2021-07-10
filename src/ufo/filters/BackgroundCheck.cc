/* * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/BackgroundCheck.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/interface/ObsFilter.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "BackgroundCheck constructor" << std::endl;

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = parameters_.test_hofx.value();

  if (parameters_.functionAbsoluteThreshold.value()) {
    for (const Variable & var : *(parameters_.functionAbsoluteThreshold.value()))
      allvars_ += var;
  }
  allvars_ += Variables(filtervars_, test_hofx);
  ASSERT(parameters_.threshold.value() ||
         parameters_.absoluteThreshold.value() ||
         parameters_.functionAbsoluteThreshold.value());
  if (parameters_.functionAbsoluteThreshold.value()) {
    ASSERT(!parameters_.threshold.value() &&
           !parameters_.absoluteThreshold.value());
    ASSERT(!parameters_.functionAbsoluteThreshold.value()->empty());
  }
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void BackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;
  const oops::Variables observed = obsdb_.obsvariables();
  const float missing = util::missingValue(missing);
  oops::Log::debug() << "BackgroundCheck obserr: " << *obserr_;

  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsVariables(), "ObsValue");
  ioda::ObsDataVector<float> obsbias(obsdb_, filtervars.toOopsVariables(), "ObsBias", false);
  std::string test_hofx = parameters_.test_hofx.value();

// Get function absolute threshold

  if (parameters_.functionAbsoluteThreshold.value()) {
//  Get function absolute threshold info from configuration
    const Variable &rtvar = parameters_.functionAbsoluteThreshold.value()->front();
    ioda::ObsDataVector<float> function_abs_threshold(obsdb_, rtvar.toOopsVariables());
    data_.get(rtvar, function_abs_threshold);

    Variables varhofx(filtervars_, test_hofx);
    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      size_t iv = observed.find(filtervars.variable(jv).variable());
//    H(x)
      std::vector<float> hofx;
      data_.get(varhofx.variable(jv), hofx);
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass &&
            (*obserr_)[iv][jobs] != util::missingValue((*obserr_)[iv][jobs])) {
          ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
          ASSERT(hofx[jobs] != util::missingValue(hofx[jobs]));

//        Threshold for current observation
          float zz = function_abs_threshold[jv][jobs];
//        Check distance from background
          if (std::abs(static_cast<float>(hofx[jobs]) - obs[jv][jobs]) > zz) {
            flagged[jv][jobs] = true;
          }
        }
      }
    }
  } else {
    Variables varhofx(filtervars_, test_hofx);
    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      size_t iv = observed.find(filtervars.variable(jv).variable());
//    H(x)
      std::vector<float> hofx;
      data_.get(varhofx.variable(jv), hofx);

//    Threshold for current variable
      std::vector<float> abs_thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
      std::vector<float> thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
      std::vector<float> bc_factor(obsdb_.nlocs(), 0.0);

      if (parameters_.absoluteThreshold.value())
        abs_thr = getScalarOrFilterData(*parameters_.absoluteThreshold.value(), data_);
      if (parameters_.threshold.value())
        thr = getScalarOrFilterData(*parameters_.threshold.value(), data_);

//    Bias Correction parameter
      if (parameters_.BiasCorrectionFactor.value())
        bc_factor = getScalarOrFilterData(*parameters_.BiasCorrectionFactor.value(), data_);

      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass &&
            (*obserr_)[iv][jobs] != util::missingValue((*obserr_)[iv][jobs])) {
          ASSERT(obs[jv][jobs] != util::missingValue(obs[jv][jobs]));
          if (parameters_.BiasCorrectionFactor.value()) {
            ASSERT(obsbias[jv][jobs] != util::missingValue(obsbias[jv][jobs]));
            bc_factor[jobs] = bc_factor[jobs]*obsbias[jv][jobs];
          }
          ASSERT(hofx[jobs] != util::missingValue(hofx[jobs]));

//        Threshold for current observation
          float zz = (thr[jobs] == std::numeric_limits<float>::max()) ? abs_thr[jobs] :
            std::min(abs_thr[jobs], thr[jobs] * (*obserr_)[iv][jobs]);
          ASSERT(zz < std::numeric_limits<float>::max() && zz > 0.0);

//        Check distance from background
          if (std::abs(static_cast<float>(hofx[jobs]) - obs[jv][jobs] - bc_factor[jobs]) > zz) {
            flagged[jv][jobs] = true;
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
