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

#include "sep_distrib_wave_energy_base_class.h"

using namespace std;

namespace sep {

   /** Default constructor. Sets member variable pointers to NULL
    * and WaveEnergySpectrumBaseClass::initialized to value false.*/
   WaveEnergySpectrumBaseClass::WaveEnergySpectrumBaseClass() {
      sim = NULL;
      simClasses = NULL;
      initialized = false;
   }
   
   /** Virtual destructor, sets member variable pointers to NULL.*/
   WaveEnergySpectrumBaseClass::~WaveEnergySpectrumBaseClass() {
      sim = NULL;
      simClasses = NULL;
   }

   /** Initializer function, sets member variable pointers to 
    * correct addresses.
    * @param sim Struct containing variables controlling execution of simulation.
    * @param simClasses Struct containing generic use simulation classes.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains parameters 
    * for this instance of class.
    * @return If true, base class initialized successfully. This does not mean that 
    * the derived class initialized correctly.*/
   bool WaveEnergySpectrumBaseClass::initialize(Simulation& sim,SimulationClasses& simClasses,
						ConfigReader& cr,const std::string& regionName) {
      this->sim = &sim;
      this->simClasses = &simClasses;
      return true;
   }
   
} // namespace sep
