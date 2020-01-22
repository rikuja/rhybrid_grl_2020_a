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

#include <climits>

#include "sep_particle_species.h"

using namespace std;

namespace sep {
   
   bool Species::finalize() {
      return true;   
   }

   const std::string& Species::getName() const {return name;}
   
   const std::string& Species::getSpeciesType() const {return speciesType;}

   bool Species::readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& name) {
      bool success = true;
      this->name = name;
      this->speciesType = "particle";

      const Real defValue = numeric_limits<Real>::infinity();
      cr.add(name+".mass","Mass of species '"+name+"' in proton masses (float).",defValue);
      cr.add(name+".charge","Charge of species '"+name+"' in elementary charges (float)",defValue);
      cr.parse();
      cr.get(name+".mass",mass);
      cr.get(name+".charge",charge);
      
      if (mass == defValue) {
	 simClasses.logger << "(SPECIES) ERROR: Could not read mass of species '" << name << "' !" << endl << write;
	 success = false;
      }
      if (charge == defValue) {
	 simClasses.logger << "(SPECIES) ERROR: Could not read charge of species '" << name << "' !" << endl << write;
	 success = false;
      }
   
      mass *= constants::MASS_PROTON;
      charge *= constants::CHARGE_ELEMENTARY;
      q_per_m = charge/mass;
   
      return success;
   }

} // namespace sep
