/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFILTERS_H_
#define TEST_UFO_OBSFILTERS_H_

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/testing/Test.h"
#include "oops/base/ObsFilters.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/FinalCheck.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsBiasParameters.h"
#include "ufo/ObsTraits.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  // TODO(wsmigaj) Report this (especially for floats) as a bug in eckit?
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::Duration> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

/// \brief Options used to configure comparison of a variable generated by a filter
/// against a reference variable loaded from the input IODA file.
class CompareVariablesParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CompareVariablesParameters, Parameters)

 public:
  /// Variable that should be contained against reference values.
  oops::RequiredParameter<ufo::Variable> test{"test", this};
  /// Variable containing the reference values.
  oops::RequiredParameter<ufo::Variable> reference{"reference", this};

  /// If set, the comparison will fail if any corresponding elements of the test and reference
  /// variables differ by more than `absTol`. If neither `absTol` nor `relTol` is set, the
  /// comparison will fail if any corresponding elements do not match exactly.
  oops::OptionalParameter<float> absTol{"absTol", this};

  /// If set, the comparison will fail if the relative difference of any pair of corresponding
  /// elements of the test and reference variables exceeds `relTol`. If neither `absTol` nor
  /// `relTol` is set, the comparison will fail if any corresponding elements do not match exactly.
  oops::OptionalParameter<float> relTol{"relTol", this};
};

// -----------------------------------------------------------------------------

/// \brief Options used to configure a test running a sequence of filters on observations
/// from a single obs space.
///
/// Note: at least one of the options whose names end in `Benchmark` or the `compareVariables`
/// option needs to be set (otherwise the test won't test much).
class ObsTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsTypeParameters, Parameters)

 public:
  /// Options used to configure the observation space.
  oops::Parameter<ioda::ObsTopLevelParameters> obsSpace{"obs space", {}, this};

  /// Options used to configure observation filters.
  oops::ObsFiltersParameters<ObsTraits> filtersParams{this};

  /// Options passed to the observation operator that will be applied during the test. If not set,
  /// no observation operator will be applied. To speed up tests of filters that depend on the
  /// values produced by the observation operator (model equivalents), these values can be
  /// precalculated and stored in the IODA file used to initialize the ObsSpace. In that case the
  /// `obs operator` keyword should be omitted and instead the `HofX` option should be set to the
  /// name of the group of ObsSpace variables containing the precalculated model equivalents.
  oops::OptionalParameter<ObsOperatorParametersWrapper> obsOperator{"obs operator", this};

  /// Group of variables storing precalculated model equivalents of observations. See the
  /// description of the `obs operator` option for more information.
  oops::OptionalParameter<std::string> hofx{"HofX", this};

  /// Options used to load GeoVaLs from a file. Required if any observation filters depend on
  /// GeoVaLs or of the `obs operator` option is set.
  oops::Parameter<GeoVaLsParameters> geovals{"geovals", {}, this};

  /// Options used to load observation diagnostics from a file. Required if any observation filters
  /// depend on observation diagnostics.
  oops::Parameter<GeoVaLsParameters> obsDiagnostics{"obs diagnostics", {}, this};

  /// Options used to configure the observation bias.
  oops::Parameter<ObsBiasParameters> obsBias{"obs bias", {}, this};

  /// Indices of observations expected to pass quality control.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> passedObservationsBenchmark{
    "passedObservationsBenchmark", this};
  /// Number of observations expected to pass quality control.
  oops::OptionalParameter<size_t> passedBenchmark{"passedBenchmark", this};

  /// Indices of observations expected to fail quality control.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> failedObservationsBenchmark{
    "failedObservationsBenchmark", this};
  /// Number of observations expected to fail quality control.
  oops::OptionalParameter<size_t> failedBenchmark{"failedBenchmark", this};

  /// An integer corresponding to one of the constants in the QCflags namespace.
  oops::OptionalParameter<int> benchmarkFlag{"benchmarkFlag", this};

  /// Indices of observations expected to receive the `benchmarkFlag` flag.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> flaggedObservationsBenchmark{
    "flaggedObservationsBenchmark", this};
  /// Number of observations expected to receive the `benchmarkFlag` flag.
  oops::OptionalParameter<size_t> flaggedBenchmark{"flaggedBenchmark", this};

  /// A list of options indicating variables whose final values should be compared against reference
  /// values loaded from the input IODA file.
  oops::OptionalParameter<std::vector<CompareVariablesParameters>> compareVariables{
    "compareVariables", this};

  /// A list of names of variables expected not to exist after all filters finish operation.
  oops::OptionalParameter<std::vector<Variable>> expectVariablesNotToExist{
    "expectVariablesNotToExist", this};

  /// If set to a string, the test will pass only if the filters produce an exception whose message
  /// contains that string.
  oops::OptionalParameter<std::string> expectExceptionWithMessage{
    "expectExceptionWithMessage", this};
};

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ObsFilters test.
class ObsFiltersParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFiltersParameters, Parameters)

 public:
  /// Time window configuration options.
  oops::RequiredParameter<eckit::LocalConfiguration> timeWindow{"time window", this};
  /// A list whose elements are used to configure tests running sequences of filters on
  /// observations from individual observation spaces.
  oops::Parameter<std::vector<ObsTypeParameters>> observations{"observations", {}, this};
};

// -----------------------------------------------------------------------------

//!
//! \brief Run the FinalCheck filter.
//!
//! This needs to be done manually if post-filters aren't run because the HofX vector
//! is not available.
//!
void runFinalCheck(oops::ObsSpace<ufo::ObsTraits> &obsspace,
                   oops::ObsDataVector<ufo::ObsTraits, int> &qcflags,
                   oops::ObsDataVector<ufo::ObsTraits, float> &obserr) {
  FinalCheck finalCheck(obsspace.obsspace(), FinalCheckParameters(),
                        qcflags.obsdatavectorptr(),
                        std::make_shared<ioda::ObsDataVector<float>>(obserr.obsdatavector()));
  finalCheck.doFilter();
}

// -----------------------------------------------------------------------------

//!
//! \brief Convert indices of observations held by this process to global observation indices.
//!
//! \param[inout] indices
//!   On input: local indices of observations held by this process. On output: corresponding
//!   global observation indices.
//! \param globalIdxFromLocalIdx
//!   A vector whose ith element is the global index of ith observation held by this process.
//!
void convertLocalObsIndicesToGlobal(std::vector<size_t> &indices,
                                    const std::vector<size_t> &globalIdxFromLocalIdx) {
  for (size_t &index : indices)
    index = globalIdxFromLocalIdx[index];
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flags satisfy the
//! predicate in at least one variable.
//!
//! \param qcFlags
//!   Vector of quality control flags for all observations.
//! \param obsDisttribution
//!   The MPI distribution used by the ObsSpace.
//! \param globalIdxFromLocalIdx
//!   A vector whose ith element is the global index of ith observation held by this process.
//! \param predicate
//!   A function object taking an argument of type int and returning bool.
//!
template <typename Predicate>
std::vector<size_t> getObservationIndicesWhere(
    const ObsTraits::ObsDataVector<int> &qcFlags,
    const eckit::mpi::Comm &comm,
    const std::vector<size_t> &globalIdxFromLocalIdx,
    const Predicate &predicate) {
  // Among the locations held on the calling process, identify those that satisfy the predicate.
  std::vector<size_t> indices;
  for (size_t locIndex = 0; locIndex < qcFlags.nlocs(); ++locIndex) {
    bool satisfied = false;
    for (size_t varIndex = 0; varIndex < qcFlags.nvars(); ++varIndex) {
      if (predicate(qcFlags[varIndex][locIndex])) {
        satisfied = true;
        break;
      }
    }
    if (satisfied) {
      indices.push_back(locIndex);
    }
  }

  // Convert the local location indices to global ones.
  convertLocalObsIndicesToGlobal(indices, globalIdxFromLocalIdx);
  // Concatenate the lists of global indices produced on all processes.
  oops::mpi::allGatherv(comm, indices);
  // Remove any duplicates (which can occur if some locations are held on more than one process).
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  return indices;
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have passed quality control in
//! at least one variable.
//!
std::vector<size_t> getPassedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                const eckit::mpi::Comm &comm,
                                                const std::vector<size_t> &globalIdxFromLocalIdx) {
  return getObservationIndicesWhere(qcFlags, comm, globalIdxFromLocalIdx,
                                    [](int qcFlag) { return qcFlag == 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have failed quality control in
//! at least one variable.
//!
std::vector<size_t> getFailedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                const eckit::mpi::Comm &comm,
                                                const std::vector<size_t> &globalIdxFromLocalIdx) {
  return getObservationIndicesWhere(qcFlags, comm, globalIdxFromLocalIdx,
                                    [](int qcFlag) { return qcFlag != 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flag is set to \p flag in
//! at least one variable.
//!
std::vector<size_t> getFlaggedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                 const eckit::mpi::Comm &comm,
                                                 const std::vector<size_t> &globalIdxFromLocalIdx,
                                                 int flag) {
  return getObservationIndicesWhere(qcFlags, comm, globalIdxFromLocalIdx,
                                    [flag](int qcFlag) { return qcFlag == flag; });
}

// -----------------------------------------------------------------------------

//!
//! Return the number of nonzero elements of \p data (on all MPI ranks, but counting each location
//! only once even if it is held on multiple ranks).
//!
size_t numNonzero(const ObsTraits::ObsDataVector<int> & data,
                  const ioda::Distribution &dist) {
  auto accumulator = dist.createAccumulator<size_t>();
  // Local reduction
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    size_t numNonzerosAtLocation = 0;
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] != 0)
        ++numNonzerosAtLocation;
    }
    accumulator->addTerm(locIndex, numNonzerosAtLocation);
  }
  // Global reduction
  return accumulator->computeResult();
}

// -----------------------------------------------------------------------------
//!
//! Return the number of elements of \p data equal to \p value (on all MPI ranks, but counting each
//! location only once even if it is held on multiple ranks).
//!
size_t numEqualTo(const ObsTraits::ObsDataVector<int> & data, int value,
                  const ioda::Distribution &dist) {
  auto accumulator = dist.createAccumulator<size_t>();
  // Local reduction
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    size_t numHitsAtLocation = 0;
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] == value)
        ++numHitsAtLocation;
    }
    accumulator->addTerm(locIndex, numHitsAtLocation);
  }
  // Global reduction
  return accumulator->computeResult();
}

// -----------------------------------------------------------------------------

template <typename T>
void expectVariablesEqual(const ObsTraits::ObsSpace &obsspace,
                          const ufo::Variable &referenceVariable,
                          const ufo::Variable &testVariable)
{
  std::vector<T> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<T> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT_EQUAL(reference, test);
}

// -----------------------------------------------------------------------------

void expectVariablesApproximatelyEqual(const ObsTraits::ObsSpace &obsspace,
                                       const ufo::Variable &referenceVariable,
                                       const ufo::Variable &testVariable,
                                       float absTol)
{
  std::vector<float> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<float> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT(oops::are_all_close_absolute(reference, test, absTol));
}

// -----------------------------------------------------------------------------

void expectVariablesRelativelyEqual(const ObsTraits::ObsSpace &obsspace,
                                    const ufo::Variable &referenceVariable,
                                    const ufo::Variable &testVariable,
                                    float relTol)
{
  std::vector<float> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<float> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT(oops::are_all_close_relative(reference, test, relTol));
}

// -----------------------------------------------------------------------------

void testFilters(size_t obsSpaceIndex, oops::ObsSpace<ufo::ObsTraits> &obspace,
                 const ObsTypeParameters &params) {
  typedef oops::GeoVaLs<ufo::ObsTraits>           GeoVaLs_;
  typedef oops::ObsDiagnostics<ufo::ObsTraits>    ObsDiags_;
  typedef oops::ObsAuxControl<ufo::ObsTraits>     ObsAuxCtrl_;
  typedef oops::ObsFilters<ufo::ObsTraits>        ObsFilters_;
  typedef oops::ObsOperator<ufo::ObsTraits>       ObsOperator_;
  typedef oops::ObsVector<ufo::ObsTraits>         ObsVector_;
  typedef oops::ObsSpace<ufo::ObsTraits>          ObsSpace_;
  typedef oops::ObsDataVector<ufo::ObsTraits, float> ObsDataVector_;

/// init QC and error
  ObsDataVector_ obserrfilter(obspace, obspace.obsvariables(), "ObsError");
  std::shared_ptr<oops::ObsDataVector<ufo::ObsTraits, int> >
    qcflags(new oops::ObsDataVector<ufo::ObsTraits, int>  (obspace, obspace.obsvariables()));

//  Create filters and run preProcess
  ObsFilters_ filters(obspace,
                      params.filtersParams,
                      qcflags, obserrfilter);
  filters.preProcess();
/// call priorFilter and postFilter if hofx is available
  oops::Variables geovars = filters.requiredVars();
  oops::Variables diagvars = filters.requiredHdiagnostics();

/// initialize zero bias
  ObsVector_ bias(obspace);
  bias.zero();

  if (params.hofx.value() != boost::none) {
///   read GeoVaLs from file if required
    const GeoVaLs_ gval(params.geovals.value(), obspace, geovars);
    if (geovars.size() > 0) {
      filters.priorFilter(gval);
      } else {
      oops::Log::info() << "Filters don't require geovals, priorFilter not called" << std::endl;
    }
///   read H(x) and ObsDiags from file
    oops::Log::info() << "HofX section specified, reading HofX from file" << std::endl;
    const std::string &hofxgroup = *params.hofx.value();
    ObsVector_ hofx(obspace, hofxgroup);
    if (diagvars.size() > 0) {
      if (params.obsDiagnostics.value().filename.value() == boost::none)
        throw eckit::UserError("Element #" + std::to_string(obsSpaceIndex) +
                     " of the 'observations' list requires an 'obs diagnostics.filename' section",
                     Here());
      oops::Log::info() << "Obs diagnostics section specified, reading obs diagnostics from file"
                        << std::endl;
    }
    const ObsDiags_ diags(params.obsDiagnostics.value(), obspace, diagvars);
    filters.postFilter(gval, hofx, bias, diags);
  } else if (params.obsOperator.value() != boost::none) {
///   read GeoVaLs, compute H(x) and ObsDiags
    oops::Log::info() << "ObsOperator section specified, computing HofX" << std::endl;
    ObsOperator_ hop(obspace, *params.obsOperator.value());
    const ObsAuxCtrl_ ybias(obspace, params.obsBias.value());
    ObsVector_ hofx(obspace);
    oops::Variables vars = hop.requiredVars();
    oops::Variables reducedVars = filters.requiredVars();
    reducedVars += ybias.requiredVars();
    vars += reducedVars;  // the reduced format is derived from the sampled format
    GeoVaLs_ gval(params.geovals.value(), obspace, vars);
    hop.computeReducedVars(reducedVars, gval);
    oops::Variables diagvars;
    diagvars += filters.requiredHdiagnostics();
    diagvars += ybias.requiredHdiagnostics();
    ObsDiags_ diags(obspace, hop.locations(), diagvars);
    filters.priorFilter(gval);
    hop.simulateObs(gval, hofx, ybias, bias, diags);
    filters.postFilter(gval, hofx, bias, diags);
    hofx.save("hofx");
  } else if (geovars.size() > 0) {
///   Only call priorFilter
    const GeoVaLs_ gval(params.geovals.value(), obspace, geovars);
    filters.priorFilter(gval);
    oops::Log::info() << "HofX or ObsOperator sections not provided for filters, " <<
                         "postFilter not called" << std::endl;
///   apply the FinalCheck filter (which should always be run after all other filters).
    runFinalCheck(obspace, *qcflags, obserrfilter);
    obserrfilter.mask(*qcflags);
  } else {
///   no need to run priorFilter or postFilter
    oops::Log::info() << "GeoVaLs not required, HofX or ObsOperator sections not " <<
                         "provided for filters, only preProcess was called" << std::endl;
///   apply the FinalCheck filter (which should always be run after all other filters).
    runFinalCheck(obspace, *qcflags, obserrfilter);
    obserrfilter.mask(*qcflags);
  }

  qcflags->save("EffectiveQC");
  const std::string errname = "EffectiveError";
  obserrfilter.save(errname);

//  Compare with known results
  bool atLeastOneBenchmarkFound = false;
  const ObsTraits::ObsSpace &ufoObsSpace = obspace.obsspace();

  if (params.passedObservationsBenchmark.value() != boost::none) {
    atLeastOneBenchmarkFound = true;
    const std::vector<size_t> &passedObsBenchmark =
        *params.passedObservationsBenchmark.value();
    const std::vector<size_t> passedObs = getPassedObservationIndices(
          qcflags->obsdatavector(), ufoObsSpace.comm(), ufoObsSpace.index());
    EXPECT_EQUAL(passedObs, passedObsBenchmark);
  }

  if (params.passedBenchmark.value() != boost::none) {
    atLeastOneBenchmarkFound = true;
    const size_t passedBenchmark = *params.passedBenchmark.value();
    size_t passed = numEqualTo(qcflags->obsdatavector(), ufo::QCflags::pass,
                               *ufoObsSpace.distribution());
    EXPECT_EQUAL(passed, passedBenchmark);
  }

  if (params.failedObservationsBenchmark.value() != boost::none) {
    atLeastOneBenchmarkFound = true;
    const std::vector<size_t> &failedObsBenchmark =
        *params.failedObservationsBenchmark.value();
    const std::vector<size_t> failedObs = getFailedObservationIndices(
          qcflags->obsdatavector(), ufoObsSpace.comm(), ufoObsSpace.index());
    EXPECT_EQUAL(failedObs, failedObsBenchmark);
  }

  if (params.failedBenchmark.value() != boost::none) {
    atLeastOneBenchmarkFound = true;
    const size_t failedBenchmark = *params.failedBenchmark.value();
    size_t failed = numNonzero(qcflags->obsdatavector(), *ufoObsSpace.distribution());
    EXPECT_EQUAL(failed, failedBenchmark);
  }

  if (params.benchmarkFlag.value() != boost::none) {
    const int flag = *params.benchmarkFlag.value();

    if (params.flaggedObservationsBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> &flaggedObsBenchmark =
          *params.flaggedObservationsBenchmark.value();
      const std::vector<size_t> flaggedObs =
          getFlaggedObservationIndices(qcflags->obsdatavector(), ufoObsSpace.comm(),
                                       ufoObsSpace.index(), flag);
      EXPECT_EQUAL(flaggedObs, flaggedObsBenchmark);
    }

    if (params.flaggedBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const size_t flaggedBenchmark = *params.flaggedBenchmark.value();
      size_t flagged = numEqualTo(qcflags->obsdatavector(), flag, *ufoObsSpace.distribution());
      EXPECT_EQUAL(flagged, flaggedBenchmark);
    }
  }

  if (params.compareVariables.value() != boost::none) {
    for (const CompareVariablesParameters &compareVariablesParams :
         *params.compareVariables.value()) {
      atLeastOneBenchmarkFound = true;

      const ufo::Variable &referenceVariable = compareVariablesParams.reference;
      const ufo::Variable &testVariable = compareVariablesParams.test;

      switch (ufoObsSpace.dtype(referenceVariable.group(), referenceVariable.variable())) {
      case ioda::ObsDtype::Integer:
        expectVariablesEqual<int>(ufoObsSpace, referenceVariable, testVariable);
        break;
      case ioda::ObsDtype::Integer_64:
        // If this variable contains epoch metadata, that means it is a DateTime represented as
        // seconds from epoch, so we reinterpret and compare the int64 as a DateTime. Otherwise,
        // fall back to comparing as int64s.
        try {
          const std::string varname = referenceVariable.fullName();
          const ioda::Variable &iodavar = ufoObsSpace.getObsGroup().vars.open(varname);
          const util::DateTime epoch = ioda::getEpochAsDtime(iodavar);
          expectVariablesEqual<util::DateTime>(ufoObsSpace, referenceVariable, testVariable);
        } catch (ioda::Exception&) {
          expectVariablesEqual<int64_t>(ufoObsSpace, referenceVariable, testVariable);
        }
        break;
      case ioda::ObsDtype::String:
        expectVariablesEqual<std::string>(ufoObsSpace, referenceVariable, testVariable);
        break;
      case ioda::ObsDtype::DateTime:
        expectVariablesEqual<util::DateTime>(ufoObsSpace, referenceVariable, testVariable);
        break;
      case ioda::ObsDtype::Bool:
        expectVariablesEqual<bool>(ufoObsSpace, referenceVariable, testVariable);
        break;
      case ioda::ObsDtype::Float:
        if (compareVariablesParams.absTol.value() == boost::none &&
            compareVariablesParams.relTol.value() == boost::none) {
          expectVariablesEqual<float>(ufoObsSpace, referenceVariable, testVariable);
        } else {
          if (compareVariablesParams.absTol.value() != boost::none) {
            const float tol = *compareVariablesParams.absTol.value();
            expectVariablesApproximatelyEqual(ufoObsSpace, referenceVariable, testVariable, tol);
          } else if (compareVariablesParams.relTol.value() != boost::none) {
            const float tol = *compareVariablesParams.relTol.value();
            expectVariablesRelativelyEqual(ufoObsSpace, referenceVariable, testVariable, tol);
          }
        }
        break;
      case ioda::ObsDtype::None:
        ASSERT_MSG(false, "Reference variable not found in observation space");
        break;
      case ioda::ObsDtype::Empty:
        ASSERT_MSG(false, "Test function is not set up to handle empty observation space");
      }
    }
  }

  if (params.expectVariablesNotToExist.value() != boost::none) {
    for (const Variable &var : *params.expectVariablesNotToExist.value()) {
      atLeastOneBenchmarkFound = true;
      EXPECT_NOT(ufoObsSpace.has(var.group(), var.variable()));
    }
  }

  if (params.expectExceptionWithMessage.value() != boost::none) {
    // We shouldn't get here, but if we do, we want the test to fail with a message that an
    // expected exception hasn't been thrown rather than that no benchmarks have been found.
    atLeastOneBenchmarkFound = true;
  }

  EXPECT(atLeastOneBenchmarkFound);
}

// -----------------------------------------------------------------------------

void runTest() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::ObsSpace<ufo::ObsTraits>     ObsSpace_;

  ObsFiltersParameters params;
  params.validateAndDeserialize(::test::TestEnvironment::config());

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
/// identify parameters used for this group of observations
    const ObsTypeParameters &typeParams = params.observations.value()[jj];

    if (typeParams.expectExceptionWithMessage.value() != boost::none) {
      const char *expectedMessage = typeParams.expectExceptionWithMessage.value()->c_str();
      EXPECT_THROWS_MSG(testFilters(jj, Test_::obspace()[jj], typeParams),
                        expectedMessage);
    } else {
      testFilters(jj, Test_::obspace()[jj], typeParams);
    }
  }
}

// -----------------------------------------------------------------------------

class ObsFilters : public oops::Test {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
 public:
  ObsFilters() {}
  virtual ~ObsFilters() {}
 private:
  std::string testid() const override {return "test::ObsFilters";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFilters/testFilters")
      { runTest(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFILTERS_H_
