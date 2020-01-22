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

#include "sep_lagr_species.h"

using namespace std;

namespace sep {
   
   bool LagrangianSpecies::finalize() {
      return true;   
   }

   const std::string& LagrangianSpecies::getName() const {return name;}
   
   const std::string& LagrangianSpecies::getSpeciesType() const {return speciesType;}

   bool LagrangianSpecies::readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& name) {
      bool success = true;
      this->name = name;
      this->speciesType = "lagrangian";

      cr.add(name+".direction","Direction of Alfven wave propagation, negative for antiparallel, positive for parallel (int).",(int)0);
      cr.parse();
      cr.get(name+".direction",propagationDirection);
      
      if (propagationDirection == 0) {
	 simClasses.logger << "(SEP LAGR SPECIES) ERROR: Could not read parameter '" << name + ".direction' value from config file" << endl << write;
	 success = false;
      } else if (propagationDirection < 0) {
	 propagationDirection = -1;
      } else {
	 propagationDirection = +1;
      }
      
      return success;
   }

} // namespace sep
