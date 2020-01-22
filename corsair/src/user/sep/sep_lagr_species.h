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

#ifndef SEP_LAGR_SPECIES_H
#define SEP_LAGR_SPECIES_H

#include <cstdlib>

#include <definitions.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace sep {
   
   struct LagrangianSpecies {
      std::string name;                 /**< Name of species.*/
      std::string speciesType;          /**< What type of particle this species represents.*/
      int32_t propagationDirection;     /**< Direction of Alfven wave propagation, -1 for antiparallel and +1 for parallel.*/
      
      bool finalize();
      const std::string& getName() const;
      const std::string& getSpeciesType() const;
      bool readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& name);   
   };
   
} // namespace sep
   
#endif
