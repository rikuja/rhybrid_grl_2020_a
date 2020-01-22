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

#ifndef SEP_TYPEDEFS_H
#define SEP_TYPEDEFS_H

#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_particle_accelerator.h"
#include "sep_lagr_definition.h"
#include "sep_lagr_species.h"
#include "sep_lagr_accelerator.h"

namespace sep {

   typedef sep::Particle<Real> PARTICLE;
   typedef sep::Species PARTICLE_SPECIES;
   typedef sep::ParticleAccelerator<PARTICLE_SPECIES,PARTICLE> ACCELERATOR;
   
   typedef sep::LagrangianParticle<Real> LAGR_PARTICLE;
   typedef sep::LagrangianSpecies LAGR_SPECIES;
   typedef sep::LagrangianAccelerator<LAGR_SPECIES,LAGR_PARTICLE> LAGR_ACCELERATOR;
   //typedef sep::ParticleList<LAGR_SPECIES,LAGR_PARTICLE> LAGR_LIST;
   
} // namespace sep

#endif
