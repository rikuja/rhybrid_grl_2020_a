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

#include <constants.h>
#include "sep_distrib_energy_container.h"
#include "sep_distrib_energy_powerlaw.h"

using namespace std;

namespace sep {

   EnergyDistribPowerlaw::EnergyDistribPowerlaw() {
      initialized = false;
      injectionEnergyMin = NAN;
      injectionEnergyMax = NAN;
      spectralIndex = NAN;
      simClasses = NULL;
   }
   
   bool energyDistribPowerlawFinalize(void*& parameters) {
      EnergyDistribPowerlaw* e = reinterpret_cast<EnergyDistribPowerlaw*>(parameters);
      delete e; e = NULL;
      return true;
   }
   
   bool energyDistribPowerlawGetDistrib(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      injectionEnergy.weight = NAN;
      return false;
   }

   bool energyDistribPowerlawGetEnergy(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      EnergyDistribPowerlaw& energyPowerlaw = *reinterpret_cast<EnergyDistribPowerlaw*>(parameters);
      
      injectionEnergy.energy = energyPowerlaw.injectionEnergyMin 
	+ (energyPowerlaw.injectionEnergyMax-energyPowerlaw.injectionEnergyMin)
	* energyPowerlaw.simClasses->random.uniform();
      
      injectionEnergy.weight = pow(injectionEnergy.energy/energyPowerlaw.injectionEnergyMin,
				   energyPowerlaw.spectralIndex);
      return true;
   }
   
   bool energyDistribPowerlawInitialize(Simulation& sim,SimulationClasses& simClasses,
					ConfigReader& cr,const std::string& regionName,
					void*& parameters) {
      simClasses.logger << "(ENERGY DISTRIB POWERLAW) Starting initialization" << endl;
      
      // Prevent multiple initializations:
      if (parameters != NULL) {
	 simClasses.logger << write;
	 return true;
      }
      
      EnergyDistribPowerlaw* energyPowerlaw = new EnergyDistribPowerlaw();
      energyPowerlaw->initialized = true;
      parameters = energyPowerlaw;
      
      energyPowerlaw->simClasses = &simClasses;
      
      // Read injection energy from config file:
      const Real DEF_VALUE = numeric_limits<Real>::infinity();
      cr.add(regionName+".injection_energy_min","Minimum injection energy (float).",DEF_VALUE);
      cr.add(regionName+".injection_energy_max","Maximum injection energy (float).",DEF_VALUE);
      cr.add(regionName+".spectral_index","Power law spectral index (float)",DEF_VALUE);
      cr.add(regionName+".units","Units in which injection energy is given (eV/keV/MeV/GeV).",string(""));
      
      string energyUnits;
      cr.parse();
      cr.get(regionName+".injection_energy_min",energyPowerlaw->injectionEnergyMin);
      cr.get(regionName+".injection_energy_max",energyPowerlaw->injectionEnergyMax);
      cr.get(regionName+".spectral_index",energyPowerlaw->spectralIndex);
      cr.get(regionName+".units",energyUnits);
      
      if (energyPowerlaw->injectionEnergyMin == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".injection_energy_min' was not found" << endl;
	 energyPowerlaw->initialized = false;
      }
      if (energyPowerlaw->injectionEnergyMax == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".injection_energy_max' was not found" << endl;
	 energyPowerlaw->initialized = false;
      }
      if (energyPowerlaw->spectralIndex == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".spectral_index' was not found" << endl;
	 energyPowerlaw->initialized = false;
      }
      
      const Real energyFactor = simClasses.constants.getEnergyInSI(energyUnits);
      if (energyFactor == numeric_limits<Real>::infinity()) {
	 simClasses.logger << "\t ERROR: Unknown energy units in parameter '";
	 simClasses.logger << regionName << ".units'. Must be (eV/keV/MeV/GeV) !" << endl;
	 energyPowerlaw->initialized = false;
      }      
      energyPowerlaw->injectionEnergyMin *= energyFactor;
      energyPowerlaw->injectionEnergyMax *= energyFactor;

      // Write initialization status and exit:
      simClasses.logger << "\t Initialization complete, status is ";
      if (energyPowerlaw->initialized == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return energyPowerlaw->initialized;
   }

}
