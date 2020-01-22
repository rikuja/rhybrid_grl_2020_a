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

#include "sep_distrib_pitch_mono.h"

using namespace std;

namespace sep {

   static PitchDistribMono pitchMono;
   
   PitchDistribMono::PitchDistribMono() {
      initialized = false;
      injectionPitch = NAN;
      simClasses = NULL;
   }
   
   bool pitchDistribMonoFinalize() {return true;}
   
   Real pitchDistribMonoGetPitch() {
      return pitchMono.injectionPitch;
   }
   
   bool pitchDistribMonoInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName) {
      // Prevent multiple initializations:
      if (pitchMono.initialized == true) return true;
      pitchMono.simClasses = &simClasses;
      
      // Read injection pitch from config file:
      const Real DEF_VALUE = numeric_limits<Real>::infinity();
      cr.add(regionName+".injection_pitch","Injection pitch (float).",DEF_VALUE);
      cr.parse();
      cr.get(regionName+".injection_pitch",pitchMono.injectionPitch);
      
      if (pitchMono.injectionPitch == DEF_VALUE) {
	 simClasses.logger << "(PITCH DISTRIB MONO) ERROR: Parameter '" << regionName+".injection_pitch' was not found" << endl << write;
	 return false;
      }
      
      pitchMono.initialized = true;
      return true;
   }

}
