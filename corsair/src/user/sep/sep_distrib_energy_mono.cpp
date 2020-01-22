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
#include "sep_distrib_energy_mono.h"

using namespace std;

namespace sep {

   EnergyDistribMono::EnergyDistribMono() {
      initialized = false;
      injectionEnergy = NAN;
      simClasses = NULL;
   }
   
   bool energyDistribMonoFinalize(void*& parameters) {
      EnergyDistribMono* e = reinterpret_cast<EnergyDistribMono*>(parameters);
      delete e; e = NULL;
      return true;
   }
   
   bool energyDistribMonoGetDistrib(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      injectionEnergy.energy = 0.0;
      return true;
   }

   bool energyDistribMonoGetEnergy(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      EnergyDistribMono& energyMono = *reinterpret_cast<EnergyDistribMono*>(parameters);
      
      injectionEnergy.energy = energyMono.injectionEnergy;
      injectionEnergy.weight = 1.0;
      return true;
   }

   bool energyDistribMonoInitialize(Simulation& sim,SimulationClasses& simClasses,
				    ConfigReader& cr,const std::string& regionName,
				    void*& parameters) {
      
      // Prevent multiple initializations:
      if (parameters != NULL) return true;

      EnergyDistribMono* energyMono = new EnergyDistribMono();
      parameters = energyMono;
      energyMono->simClasses = &simClasses;
      energyMono->initialized = true;

      // Read injection energy from config file:
      const Real DEF_VALUE = numeric_limits<Real>::infinity();
      cr.add(regionName+".injection_energy","Injection energy (float).",DEF_VALUE);
      cr.add(regionName+".units","Units in which injection energy is given (eV/keV/MeV/GeV).",string(""));
      
      string energyUnits;
      cr.parse();
      cr.get(regionName+".injection_energy",energyMono->injectionEnergy);
      cr.get(regionName+".units",energyUnits);
      
      if (energyMono->injectionEnergy == DEF_VALUE) {
	 simClasses.logger << "(ENERGY DISTRIB MONO) ERROR: Parameter '" << regionName+".injection_energy' was not found" << endl << write;
	 energyMono->initialized = false;
      }
      
      const Real energyFactor = simClasses.constants.getEnergyInSI(energyUnits);
      if (energyFactor == numeric_limits<Real>::infinity()) {
	 simClasses.logger << "(SEP ENERGY DISTRIB MONO) ERROR: Unknown energy units in parameter '";
	 simClasses.logger << regionName << ".units'. Must be (eV/keV/MeV/GeV) !" << endl;
	 energyMono->initialized = false;
      }      
      energyMono->injectionEnergy *= energyFactor;
      
      simClasses.logger << "(SEP ENERGY DISTRIB MONO) Initialization complete, status is ";
      if (energyMono->initialized == true) simClasses.logger << "SUCCESS" << endl;
      else simClasses.logger << "FAILURE" << endl;
      
      simClasses.logger << "\t Injection energy: " << energyMono->injectionEnergy << endl;
      simClasses.logger << "\t Units           : '" << energyUnits << "'" << endl << write;
      return true;
   }

}
