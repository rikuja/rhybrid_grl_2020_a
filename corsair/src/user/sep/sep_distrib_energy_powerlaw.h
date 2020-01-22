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

#ifndef SEP_DISTRIP_ENERGY_POWERLAW_H
#define SEP_DISTRIB_ENERGY_POWERLAW_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_distrib_energy_container.h"

namespace sep {

   struct EnergyDistribPowerlaw {
      EnergyDistribPowerlaw();
      
      bool initialized;
      Real injectionEnergyMin;
      Real injectionEnergyMax;
      Real spectralIndex;
      SimulationClasses* simClasses;
   };
   
   bool energyDistribPowerlawFinalize(void*& parameters);
   bool energyDistribPowerlawGetDistrib(distrib::InjectionEnergy& injectionEnergy,void* parameters);
   bool energyDistribPowerlawGetEnergy(distrib::InjectionEnergy& injectionEnergy,void* parameters);
   bool energyDistribPowerlawInitialize(Simulation& sim,SimulationClasses& simClasses,
					ConfigReader& cr,const std::string& regionName,
					void*& parameters);
   
} // namespace sep

#endif
