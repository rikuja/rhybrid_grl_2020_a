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

#ifndef SEP_INIT_H
#define SEP_INIT_H

#include <user.h>

namespace sep {
   
   bool earlyInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,std::vector<ParticleListBase*>& particleLists);
   bool finalize(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
   bool lateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const ObjectFactories& objectFactories,std::vector<ParticleListBase*>& particleLists);
   bool registerObjectMakers(ObjectFactories& objectFactories);
   bool runTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
   
} // namespace sep

#endif