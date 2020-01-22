/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#include <main.h>

#include "sep_typedefs.h"
#include "sep_register_particle_propagators.h"

#include "sep_particle_propagator_rk2.h"
#include "sep_particle_propagator_shock_rk2.h"                  
#include "sep_particle_propagator_coronal_rk2.h"

// Temporarily here
#include "sep_object_wrapper.h"
#include "sep_particle_splitter.h"
// End temp

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * particle propagators.
    * @return If true, propagators were registered successfully.*/
   bool registerParticlePropagators() {
      bool success = true;

      ObjectFactories& objectFactories = corsair::getObjectWrapper().objectFactories;

      if (objectFactories.particlePropagators.registerMaker("ParticlePropagator",
							    ParticlePropagMaker<PARTICLE_SPECIES,PARTICLE,ACCELERATOR>) == false) 
	success = false;
      if (objectFactories.particlePropagators.registerMaker("ParticlePropagatorShock",
							    ParticlePropagShockMaker<PARTICLE_SPECIES,PARTICLE,ACCELERATOR>) == false) 
	success = false;
      if (objectFactories.particlePropagators.registerMaker("ParticlePropagatorCoronal",
							    ParticlePropagCoronalMaker<PARTICLE_SPECIES,PARTICLE,ACCELERATOR>) == false) 
	success = false;
      
      #warning Create registering file for splitters, good for now
      if (sep::getObjectWrapper().splitterFactory.registerMaker("ParticleSplitter",
								ParticleSplitterMaker<PARTICLE_SPECIES,PARTICLE>) == false) 
	success = false;

      return success;
   }
   
} // namespace sep
