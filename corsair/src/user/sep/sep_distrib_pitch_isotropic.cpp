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
#include <cmath>

#include "sep_distrib_pitch_isotropic.h"

using namespace std;

namespace sep {

   static PitchDistribIsotropic pitchIsotropic;
   
   PitchDistribIsotropic::PitchDistribIsotropic() {
      initialized = false;
      simClasses = NULL;
   }
   
   bool pitchDistribIsotropicFinalize() {return true;}
   
   Real pitchDistribIsotropicGetPitch() {
      return 1.0 - 2.0*pitchIsotropic.simClasses->random.uniform();
   }
   
   bool pitchDistribIsotropicInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName) {
      // Prevent multiple initializations:
      if (pitchIsotropic.initialized == true) return true;
      
      pitchIsotropic.simClasses = &simClasses;
      pitchIsotropic.initialized = true;
      return true;
   }

}
