/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "sep_injectors_common.h"

using namespace std;

namespace sep {
   
   particleinjector::frame::Frame getInjectionFrame(const std::string& frame) {
      if (frame == "simulation") {
	 return particleinjector::frame::SIMULATION;
      } else if (frame == "antiparallel") {
	 return particleinjector::frame::ANTIPARALLEL;
      } else if (frame == "parallel") {
	 return particleinjector::frame::PARALLEL;
      } else if (frame == "plasma") {
	   return particleinjector::frame::PLASMA;
      } else {
	 return particleinjector::frame::UNKNOWN;
      }
      return particleinjector::frame::UNKNOWN;
   }
   
   particleinjector::region::Region getInjectionRegion(const std::string& region) {
      if (region == "upstream") {
	 return particleinjector::region::UPSTREAM;
      } else if (region == "downstream") {
	 return particleinjector::region::DOWNSTREAM;
      } else if (region == "shock_upstream") {
	 return particleinjector::region::SHOCK_UPSTREAM;
      } else if (region == "shock_downstream") {
	 return particleinjector::region::SHOCK_DOWNSTREAM;
      } else {
	 return particleinjector::region::UNKNOWN;
      }
      return particleinjector::region::UNKNOWN;
   }

} // namespace sep
