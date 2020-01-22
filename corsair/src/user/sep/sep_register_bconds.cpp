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
#include "sep_register_bconds.h"

#include "sep_particle_boundary_cond.h"
#include "sep_particle_boundary_cond_double_shock.h"
#include "sep_particle_boundary_cond_shock.h"
#include "sep_lagr_boundary_cond.h"
#include "sep_particle_boundary_cond_none.h"

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * particle and wave packet boundary condition functions.
    * @return If true, boundary condition functions were registered successfully.*/
   bool registerBoundaryConditions() {
      bool success = true;

      ObjectFactories& objectFactories = corsair::getObjectWrapper().objectFactories;
      if (objectFactories.particleBoundaryConditions.registerMaker("IonBoundaryCondRemove",PBCondMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleBoundaryConditions.registerMaker("IonBoundaryCondDoubleShock",
								   PBCondDoubleShockMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleBoundaryConditions.registerMaker("IonBoundaryCondShock",PBCondShockMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleBoundaryConditions.registerMaker("IonBoundaryCondNone",ParticleBCondNoneMaker<PARTICLE_SPECIES,PARTICLE>) == false) success = false;
      if (objectFactories.particleBoundaryConditions.registerMaker("LagrBoundaryCondRemove",LagrBCondMaker<LAGR_SPECIES,LAGR_PARTICLE>) == false) success = false;
      if (objectFactories.particleBoundaryConditions.registerMaker("LagrBoundaryCondNone",ParticleBCondNoneMaker<LAGR_SPECIES,LAGR_PARTICLE>) == false) success = false;
      return success;
   }
   
} // namespace sep
