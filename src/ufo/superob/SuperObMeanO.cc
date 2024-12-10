/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObMeanO.h"

namespace ufo {

static SuperObMaker<SuperObMeanO> makerSuperObMeanO_("mean obs");

SuperObMeanO::SuperObMeanO(const GenericSuperObParameters & params,
                           const ObsFilterData & obsdb,
                           const std::vector<bool> & apply,
                           const Variables & filtervars,
                           const ioda::ObsDataVector<int> & flags,
                           std::vector<std::vector<bool>> & flagged)
  : SuperObBase(params, obsdb, apply, filtervars, flags, flagged)
{}

void SuperObMeanO::computeSuperOb(const std::vector<std::size_t> & locs,
                                  const std::vector<float> & obs,
                                  const std::vector<float> & /*hofx*/,
                                  const ioda::ObsDataRow<int> & flags,
                                  std::vector<float> & superobs,
                                  std::vector<bool> & flagged) const {
  const float missing = util::missingValue<float>();

  float sumO = 0.0f;
  int count = 0;
  // The superob location is set to the first valid entry in the record.
  int superobloc = -1;

  for (std::size_t jloc : locs) {
    const float obsValue = obs[jloc];
    // Only consider locations which have valid observation values
    // and are either passing QC or have been marked as passive (i.e. H(x) is computed
    // but the observation is not assimilated).
    if ((flags[jloc] != QCflags::pass && flags[jloc] != QCflags::passive) ||
        obsValue == missing) {
      continue;
    }
    sumO += obsValue;
    count++;
    if (superobloc < 0) {
      superobloc = jloc;
    }
  }
  if (count > 0) {
    superobs[superobloc] = sumO / count;
    // The superob location is not flagged.
    flagged[superobloc] = false;
  }
}

}  // namespace ufo
