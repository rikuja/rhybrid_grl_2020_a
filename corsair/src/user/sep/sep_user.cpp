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

/** This file includes definitions for callback functions defined in user.h header.
 * These callbacks are called by Corsair main loop. Functions here simply call
 * functions defined in xxx .
 */

#include <cstdlib>
#include <iostream>

#include <user.h>

#include <sep_init.h>
#include <sep_propagate.h>
#include "sep_register_dataoperators.h"

using namespace std;

bool propagate(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   return sep::propagate(sim,simClasses,particleLists);
}

bool userEarlyInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,vector<ParticleListBase*>& particleLists) {
   return sep::earlyInitialization(sim,simClasses,cr,particleLists);
}
   
bool userLateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			    const ObjectFactories& objectFactories,vector<ParticleListBase*>& particleLists) {
   return sep::lateInitialization(sim,simClasses,cr,objectFactories,particleLists);
}

bool userFinalization(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   return sep::finalize(sim,simClasses,particleLists);
}

bool userRunTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists) {
   return sep::runTests(sim,simClasses,particleLists);
}

bool registerObjectMakers(ObjectFactories& objectFactories) {
   return sep::registerObjectMakers(objectFactories);
}

bool registerDataOperators() {
   return sep::registerDataOperators();
}
