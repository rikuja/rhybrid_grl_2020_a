/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2013 Finnish Meteorological Institute
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

#include <cstdlib>
#include <iostream>

#include "sep_distrib_pitch_container.h"

using namespace std;

namespace sep {

   PitchDistribution::PitchDistribution() {
      finalize = NULL;
      getPitch = NULL;
      initialize = NULL;
   }
   
   PitchDistribution::PitchDistribution(finalizePitchDistrib finalize,getPitchFunction getPitch,initializePitchDistrib initialize):
     finalize(finalize), getPitch(getPitch), initialize(initialize) { 
     
   }

   DistribPitchContainer::DistribPitchContainer() { }
   
   DistribPitchContainer::~DistribPitchContainer() { 
      for (map<string,PitchDistribution>::iterator it=pitchDistribs.begin(); it!=pitchDistribs.end(); ++it) {
	 (*it->second.finalize)();
      }
   }
   
   bool DistribPitchContainer::getDistribution(const std::string& name,finalizePitchDistrib& finalize,
					       getPitchFunction& getPitch,initializePitchDistrib& initialize) {
      map<string,PitchDistribution>::iterator it = pitchDistribs.find(name);
      if (it == pitchDistribs.end()) return false;
      
      finalize   = it->second.finalize;
      getPitch   = it->second.getPitch;
      initialize = it->second.initialize;
      return true;
   }
   
   bool DistribPitchContainer::registerDistribution(const std::string& name,finalizePitchDistrib finalize,
						    getPitchFunction getPitch,initializePitchDistrib initialize) {
      if (finalize == NULL) return false;
      if (getPitch == NULL) return false;
      if (initialize == NULL) return false;
      if (pitchDistribs.find(name) != pitchDistribs.end()) return false;
      pitchDistribs[name] = PitchDistribution(finalize,getPitch,initialize);
      return true;
   }

} // namespace sep
