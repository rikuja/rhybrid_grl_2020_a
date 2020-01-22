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
#include "sep_register_lagr_propagators.h"

#include "sep_lagr_propagator_rk2.h"
#include "sep_lagr_propagator_shock_rk2.h"

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * particle propagators.
    * @return If true, propagators were registered successfully.*/
   bool registerWavePacketPropagators() {
      bool success = true;

      ObjectFactories& objectFactories = corsair::getObjectWrapper().objectFactories;

      if (objectFactories.particlePropagators.registerMaker("LagrPropagator",
							    LagrPropagMaker<LAGR_SPECIES,LAGR_PARTICLE,LAGR_ACCELERATOR>) == false) 
	success = false;
      if (objectFactories.particlePropagators.registerMaker("LagrPropagatorShock",
							    LagrPropagShockMaker<LAGR_SPECIES,LAGR_PARTICLE,LAGR_ACCELERATOR>) == false) 
	success = false;
      
      return success;
   }
   
} // namespace sep
