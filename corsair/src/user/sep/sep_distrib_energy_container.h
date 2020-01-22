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

#ifndef SEP_DISTRIB_ENERGY_CONTAINER_H
#define SEP_DISTRIB_ENERGY_CONTAINER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include <map>

namespace sep {

   typedef bool (*finalizeEnergyDistrib)(void*& parameters);
   
   namespace distrib {
      struct InjectionEnergy {
	 Real thermalEnergy;
	 Real energy;
	 Real weight;
	 int32_t N;
	 int32_t N_particles;
	 int32_t superResolution;
	 int32_t N_super;

	 InjectionEnergy();
      };
   }

   typedef bool (*getEnergyDistribFunction)(distrib::InjectionEnergy& energy,void* parameters);
   
   /** Get injection energy and the value of distribution function.
    * Note that energy is given per amu.
    * @param thermalEnergy Maxwell-Boltzmann thermal energy per amu.
    * @param energy Injection energy is written here.
    * @param weight Value of distribution function at injection energy is written here.
    * @param N
    * @param N_particles
    * @return If true, injection energy was calculated successfully.*/
   typedef bool (*getEnergyFunction)(distrib::InjectionEnergy& injectionEnergy,void* parameters);
   
   typedef bool (*initializeEnergyDistrib)(Simulation& sim,SimulationClasses& simClasses,
					   ConfigReader& cr,const std::string& regionName,
					   void*& parameters);
   
   struct EnergyDistribution {
      finalizeEnergyDistrib finalize;
      getEnergyDistribFunction getDistrib;
      getEnergyFunction getEnergy;
      initializeEnergyDistrib initialize;
      
      EnergyDistribution();
      EnergyDistribution(finalizeEnergyDistrib finalize,getEnergyDistribFunction getDistrib,
			 getEnergyFunction getEnergy,initializeEnergyDistrib initialize);
   };
   
   class DistribEnergyContainer {
    public:
      DistribEnergyContainer();
      ~DistribEnergyContainer();
      
      bool getDistribution(const std::string& name,finalizeEnergyDistrib& finalize,
			   getEnergyDistribFunction& getEnergyDistrib,
			   getEnergyFunction& getEnergy,initializeEnergyDistrib& initialize);
      
      bool registerDistribution(const std::string& name,finalizeEnergyDistrib finalize,
				getEnergyDistribFunction getEnergyDistrib,
				getEnergyFunction getEnergy,initializeEnergyDistrib initialize);
    private:
      std::map<std::string,EnergyDistribution> energyDistribs;
	
   };

} // namespace sep

#endif
