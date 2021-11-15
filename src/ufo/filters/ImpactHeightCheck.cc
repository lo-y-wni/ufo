/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ImpactHeightCheck.h"

#include <Eigen/Core>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// ImpactHeightCheck: Calculate the impact height for model profiles, and reject
/// any observations which are outside the range of model impact heights.  Check
/// for any sharp refractivity gradients, and reject any observations below them.

ImpactHeightCheck::ImpactHeightCheck(
        ioda::ObsSpace & obsdb,
        const Parameters_ & parameters,
        std::shared_ptr<ioda::ObsDataVector<int> > flags,
        std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ImpactHeightCheck constructor" << std::endl;
  allvars_ += Variable("ObsDiag/refractivity");
  allvars_ += Variable("ObsDiag/model_heights");
  allvars_ += Variable("MetaData/impact_parameter");
  allvars_ += Variable("MetaData/earth_radius_of_curvature");

  // It is essential for observations to be grouped according to (e.g.) profile number
  // (unless there is only one profile in the sample, which would be very unusual).
  // Throw an exception if the "group variable" configuration option is missing.
  if (obsdb.obs_group_vars().empty())
    throw eckit::BadParameter("group variables configuration is empty.", Here());
}

// -----------------------------------------------------------------------------

ImpactHeightCheck::~ImpactHeightCheck() {
  oops::Log::trace() << "ImpactHeightCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Apply the filter: remove observations which are outside the range of
/// impact heights

void ImpactHeightCheck::applyFilter(const std::vector<bool> & apply,
                                    const Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ImpactHeightCheck post-filter" << std::endl;
  const oops::Variables observed = obsdb_.obsvariables();
  const float missingFloat = util::missingValue(missingFloat);

  // Get the refractivity from the obs diagnostics, including the number of
  // vertical levels on which the refractivity has been calculated (nRefLevels)
  Variable refractivityVariable = Variable("ObsDiag/refractivity");
  oops::Log::debug() << data_.nlevs(refractivityVariable) << std::endl;
  const size_t nRefLevels = data_.nlevs(refractivityVariable);
  std::vector<std::vector<float>> refractivity;

  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    std::vector<float> inputData;
    data_.get(refractivityVariable, static_cast<int>(iLevel), inputData);
    refractivity.push_back(inputData);
  }

  // For the benefits of debugging, output the refractivity for the first
  // observation
  oops::Log::debug() << "Refractivity(first ob) ";
  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    oops::Log::debug() << refractivity[iLevel][0] << " ";
  }
  oops::Log::debug() << std::endl;

  // Get the height of the levels on which the refractivity has been calculated.
  // Must be the same length as the array defining the refractivity.
  Variable modelHeightsVariable = Variable("ObsDiag/model_heights");
  oops::Log::debug() << data_.nlevs(modelHeightsVariable) << std::endl;
  if (data_.nlevs(modelHeightsVariable) != nRefLevels) {
    throw eckit::BadValue("Model heights and refractivity must have the same number of levels",
                          Here());
  }
  std::vector<std::vector<float>> modelHeights;

  // Read the heights of the refractivity levels
  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    std::vector<float> inputData;
    data_.get(modelHeightsVariable, static_cast<int>(iLevel), inputData);
    modelHeights.push_back(inputData);
  }

  // For debugging, output the heights of the refractivity levels for the first
  // observation.
  oops::Log::debug() << "Model heights (first ob) ";
  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    oops::Log::debug() << modelHeights[iLevel][0] << " ";
  }
  oops::Log::debug() << std::endl;

  // Read in the observation impact parameter for each observation
  Variable impactVariable = Variable("MetaData/impact_parameter");
  std::vector<float> impactParameter;
  data_.get(impactVariable, impactParameter);

  // Read in the earth's radius of curvature for each observation
  Variable radiusCurvatureParameter = Variable("MetaData/earth_radius_of_curvature");
  std::vector<float> radiusCurvature;
  data_.get(radiusCurvatureParameter, radiusCurvature);

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // For each variable, perform the filter
  for (size_t iFilterVar = 0; iFilterVar < filtervars.nvars(); ++iFilterVar) {
    const size_t iVar = observed.find(filtervars.variable(iFilterVar).variable());

    // Loop over the unique profiles
    for (size_t iProfile : record_numbers) {
      const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);

      // Find the observation with the smallest impact parameter (i.e. the
      // bottom of the profile).  If all impactParameters are missing, then this
      // should choose the first observation.
      auto it_nonmissing = std::min_element(obs_numbers.begin(), obs_numbers.end(),
                                            [missingFloat, &impactParameter]
                                            (const size_t & a, const size_t & b)
                                            {return impactParameter[a] == missingFloat ?
                                                    false :
                                                    (impactParameter[b] == missingFloat ?
                                                     true :
                                                     impactParameter[a] < impactParameter[b]);});
      const size_t bottomOb = *it_nonmissing;

      // Load the refractivity profile and the associated model heights for this observation,
      // cleansing any missing data
      std::vector<float> refracProfile;
      std::vector<float> heightProfile;
      for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel)
        if (refractivity[iLevel][bottomOb] != missingFloat &&
            modelHeights[iLevel][bottomOb] != missingFloat) {
          refracProfile.push_back(refractivity[iLevel][bottomOb]);
          heightProfile.push_back(modelHeights[iLevel][bottomOb]);
        }

      if (refracProfile.size() < 2) {
        oops::Log::error() << "Should have at least two valid points in every profile:" <<
                              std::endl << "size = " << refracProfile.size() << "  " <<
                              "bottomOb = " << bottomOb << std::endl;
        flagged[iFilterVar][bottomOb] = true;
        continue;
      }

      oops::Log::debug() << "Min and max refrac for profile " << iProfile <<
                            " is " << refracProfile[0] << " to " << refracProfile.back() <<
                            std::endl;

      const std::vector<float> & gradient = calcVerticalGradient(refracProfile, heightProfile);
      // Output the calculated refractivity gradient for the first profile
      if (iProfile == record_numbers[0]) {
          oops::Log::debug() << "Gradient found to be" << std::endl;
          for (float grad : gradient)
              oops::Log::debug() << grad << "  ";
          oops::Log::debug() << std::endl;
      }

      float sharpGradientImpact = std::numeric_limits<float>::lowest();
      // Search for sharp gradients (super-refraction) starting at the top of the profile
      for (int iLevel = static_cast<int>(nRefLevels)-2; iLevel >= 0; --iLevel) {
        size_t thisLevel = static_cast<size_t>(iLevel);
        if (gradient[thisLevel] != missingFloat &&
            gradient[thisLevel] < parameters_.gradientThreshold.value()) {
          sharpGradientImpact = calcImpactHeight(refracProfile[thisLevel],
                                                 heightProfile[thisLevel],
                                                 radiusCurvature[bottomOb]);
          oops::Log::info() << "Sharp refractivity gradient of " << gradient[thisLevel] <<
                               " found at " << thisLevel << "  " << sharpGradientImpact <<
                               std::endl;
          oops::Log::debug() << thisLevel << "   " << refracProfile[thisLevel] << "   " <<
                                heightProfile[thisLevel] << "   " <<
                                radiusCurvature[thisLevel] << std::endl;
          break;
        }
      }

      // Loop over all observations in the profile
      for (size_t jobs : obs_numbers) {
        // Check that this observation should be considered in this routine
        if (apply[jobs] && (*flags_)[iVar][jobs] == QCflags::pass) {
          // Reject observation if it is below the minimum (either surface or sharp gradient)
          const float obsImpactHeight = impactParameter[jobs] - radiusCurvature[jobs];
          if (parameters_.verboseOutput.value())
            oops::Log::debug() << "Checking minimum height " << obsImpactHeight << "   " <<
                                  sharpGradientImpact + parameters_.sharpGradientOffset.value() <<
                                  "   " << calcImpactHeight(refracProfile[0], heightProfile[0],
                                                            radiusCurvature[jobs]) +
                                  parameters_.surfaceOffset.value() << "   " <<
                                  (obsImpactHeight < sharpGradientImpact +
                                        parameters_.sharpGradientOffset.value()) << "   " <<
                                  (obsImpactHeight < calcImpactHeight(refracProfile[0],
                                                                      heightProfile[0],
                                                                      radiusCurvature[jobs]) +
                                      parameters_.surfaceOffset.value()) << std::endl;
          if (obsImpactHeight < sharpGradientImpact + parameters_.sharpGradientOffset.value() ||
              obsImpactHeight < calcImpactHeight(refracProfile[0], heightProfile[0],
                                                 radiusCurvature[jobs]) +
                                                 parameters_.surfaceOffset.value())
            flagged[iFilterVar][jobs] = true;

          if (parameters_.verboseOutput.value())
            oops::Log::debug() << "Checking maximum height " << obsImpactHeight << "   " <<
                                  calcImpactHeight(refracProfile.back(), heightProfile.back(),
                                                   radiusCurvature[jobs]) << std::endl;
          // Reject observation if it is above the maximum
          if (obsImpactHeight > calcImpactHeight(refracProfile.back(), heightProfile.back(),
                                                 radiusCurvature[jobs]))
            flagged[iFilterVar][jobs] = true;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Calculate the vertical gradient of refractivity, assuming that any missing
/// data has already been removed
std::vector<float> ImpactHeightCheck::calcVerticalGradient(
        const std::vector<float> & refrac,
        const std::vector<float> & height) const {
  std::vector<float> gradient;
  for (size_t iLevel = 0; iLevel < refrac.size()-1; ++iLevel) {
    gradient.push_back((refrac[iLevel+1] - refrac[iLevel]) / (height[iLevel+1] - height[iLevel]));
  }
  return gradient;
}

// -----------------------------------------------------------------------------
/// Calculate the impact parameter from the refractivity, impact height and
/// radius of curvature
float ImpactHeightCheck::calcImpactHeight(float refractivity,
                                          float modelHeight,
                                          float radiusCurv) const {
  return 1.0E-6f * refractivity * (radiusCurv + modelHeight) + modelHeight;
}


// -----------------------------------------------------------------------------
void ImpactHeightCheck::print(std::ostream & os) const {
  os << "ImpactHeightCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
