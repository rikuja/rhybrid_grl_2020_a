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

#include "sep_distrib_energy_container.h"

using namespace std;

namespace sep {

   namespace distrib {
      InjectionEnergy::InjectionEnergy() {
	 superResolution = 1;
	 N_super = 0;
      }
   }
   
   EnergyDistribution::EnergyDistribution() {
      finalize = NULL;
      getDistrib = NULL;
      getEnergy = NULL;
      initialize = NULL;
   }
   
   EnergyDistribution::EnergyDistribution(finalizeEnergyDistrib finalize,getEnergyDistribFunction getDistrib,
					  getEnergyFunction getEnergy,initializeEnergyDistrib initialize):
     finalize(finalize), getDistrib(getDistrib), getEnergy(getEnergy), initialize(initialize) { 
     
   }

   DistribEnergyContainer::DistribEnergyContainer() { }
   
   DistribEnergyContainer::~DistribEnergyContainer() { }
   
   bool DistribEnergyContainer::getDistribution(const std::string& name,finalizeEnergyDistrib& finalize,
						getEnergyDistribFunction& getDistrib,
						getEnergyFunction& getEnergy,initializeEnergyDistrib& initialize) {
      map<string,EnergyDistribution>::iterator it = energyDistribs.find(name);
      if (it == energyDistribs.end()) return false;
      
      finalize   = it->second.finalize;
      getDistrib = it->second.getDistrib;
      getEnergy  = it->second.getEnergy;
      initialize = it->second.initialize;
      return true;
   }
   
   bool DistribEnergyContainer::registerDistribution(const std::string& name,finalizeEnergyDistrib finalize,
						     getEnergyDistribFunction getDistrib,
						     getEnergyFunction getEnergy,initializeEnergyDistrib initialize) {
      if (finalize == NULL) return false;
      if (getDistrib == NULL) return false;
      if (getEnergy == NULL) return false;
      if (initialize == NULL) return false;
      if (energyDistribs.find(name) != energyDistribs.end()) return false;
      energyDistribs[name] = EnergyDistribution(finalize,getDistrib,getEnergy,initialize);
      return true;
   }

} // namespace sep
