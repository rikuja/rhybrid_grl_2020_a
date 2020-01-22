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
#include "sep_register_particle_injectors.h"

#include "sep_particle_injector_cell.h"
#include "sep_particle_injector_homogeneous.h"
#include "sep_particle_injector_ng_reames_1994.h"
#include "sep_particle_injector_shock.h"

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * particle injectors.
    * @return If true, injectors were registered successfully.*/
   bool registerParticleInjectors() {
      bool success = true;

      ObjectFactories& objectFactories = corsair::getObjectWrapper().objectFactories;
      
      if (objectFactories.particleInjectors.registerMaker("IonInjectorShock",PIShockMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleInjectors.registerMaker("IonInjectorCell",PICellMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleInjectors.registerMaker("IonInjectorHomog",PIHomogMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleInjectors.registerMaker("IonInjectorNgReames1994",PINgReames1994Maker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      
      return success;
   }
   
} // namespace sep
