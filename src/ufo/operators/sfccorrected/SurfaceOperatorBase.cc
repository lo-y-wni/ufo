/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/sfccorrected/SurfaceOperatorBase.h"

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

namespace ufo {

SurfaceOperatorBase::SurfaceOperatorBase(const std::string &, const Parameters_ & params)
  : params_(params)
{}

SurfaceOperatorFactory::SurfaceOperatorFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name + " already registered in ufo::SurfaceOperatorFactory.", Here());
  getMakers()[name] = this;
}

std::unique_ptr<SurfaceOperatorBase>
SurfaceOperatorFactory::create(const std::string & name,
                               const Parameters_ & params) {
  oops::Log::trace() << "SurfaceOperatorBase::create starting" << std::endl;
  typename std::map<std::string, SurfaceOperatorFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::SurfaceOperatorFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<SurfaceOperatorBase> ptr = jloc->second->make(name, params);
  oops::Log::trace() << "SurfaceOperatorBase::create done" << std::endl;
  return ptr;
}

}  // namespace ufo
