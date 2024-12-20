/*
 * (C) Crown Copyright 2024, Met Office
 * (C) Copyright 2024 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
#define UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// enum type for surface correction type, and ParameterTraitsHelper for it
enum class SfcCorrectionType {
  UKMO, WRFDA, GSL
};
struct SfcCorrectionTypeParameterTraitsHelper {
  typedef SfcCorrectionType EnumType;
  static constexpr char enumTypeName[] = "SfcCorrectionType";
  static constexpr util::NamedEnumerator<SfcCorrectionType> namedValues[] = {
    { SfcCorrectionType::UKMO, "UKMO" },
    { SfcCorrectionType::WRFDA, "WRFDA" },
    { SfcCorrectionType::GSL, "GSL" }
  };
};

}  // namespace ufo

namespace oops {

/// Extraction of SfcCorrectionType parameters from config
template <>
struct ParameterTraits<ufo::SfcCorrectionType> :
    public EnumParameterTraits<ufo::SfcCorrectionTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the SfcCorrected operator.
class ObsSfcCorrectedParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSfcCorrectedParameters, ObsOperatorParametersBase)

 public:
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
      "variables",
      "List of variables to be simulated which must be a subset of the simulated variables "
      "in the ObsSace",
      this};

  oops::Parameter<std::string> geovarGeomZ{
      "geovar_geomz",
      "Model variable for height of vertical levels, geopotential heights will be converted",
      "height_above_mean_sea_level",
      this};

  oops::Parameter<std::string> geovarSfcGeomZ{
      "geovar_sfc_geomz",
      "Model variable for surface height, geopotential heights will be converted",
      "height_above_mean_sea_level_at_surface",
      this};

  oops::Parameter<std::string> obsHeightName{
      "station_altitude",
      "stationElevation",
      this};

  // WRFDA method uses model surface and level 1 data.
  // UKMO method uses model surface and 2000m data.
  // GSL method ??
  oops::Parameter<SfcCorrectionType> correctionType{
      "correction scheme to use",
      "Scheme used for correction ('WRFDA' or 'UKMO' or 'GSL')",
      SfcCorrectionType::WRFDA,
      this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_SFCCORRECTED_OBSSFCCORRECTEDPARAMETERS_H_
