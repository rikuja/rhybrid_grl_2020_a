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

#ifndef SEP_PROPAGATE_H
#define SEP_PROPAGATE_H

#include <vector>

#include <simulation.h>
#include <simulationclasses.h>
#include <particle_list_base.h>

namespace sep {
   
   bool propagate(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
   void recalculateCellVolumes(Simulation& sim,SimulationClasses& simClasses);
   void recalculateTimestep(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists);
   
} // namespace sep

#endif
