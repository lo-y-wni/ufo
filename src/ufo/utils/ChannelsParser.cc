/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/ChannelsParser.h"

#include <algorithm>

#include <sstream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace ufo {

/// Function to split string on delimiter
std::vector<std::string> splitString(const std::string& str, char delim)
{
  std::vector<std::string> result;
  std::stringstream s(str);
  std::string substr;
  while (std::getline(s, substr, delim)) {
    result.push_back(substr);
  }
  return result;
}

/// Function to parse channels (supports commas for separating channels
//  and channel ranges and dashes for channel ranges).
//  For example: 1-5, 9, 13-45
std::vector<int> parseChannels(const std::string& str)
{
  // split string by commas to get individual channels or ranges
  std::vector<std::string> ranges = splitString(str, ',');

  std::vector<int> channels;
  for (int irange = 0; irange < ranges.size(); irange++) {
    // split the element by dashes (in case it is a range)
    std::vector<std::string> range = splitString(ranges[irange], '-');
    ASSERT((range.size() == 1) || (range.size() == 2));
    // add a single channel
    if (range.size() == 1) {
      // add a single channel
      channels.push_back(std::stoi(range[0]));
    } else if (range.size() == 2) {
      // add a range
      int start = std::stoi(range[0]);
      int stop  = std::stoi(range[1]);
      for (int ch = start; ch <= stop; ch++) {
        channels.push_back(ch);
      }
    }
  }

  // sort and remove duplicates
  std::sort(channels.begin(), channels.end());
  channels.erase(std::unique(channels.begin(), channels.end()), channels.end());

  return channels;
}

}  // namespace ufo
