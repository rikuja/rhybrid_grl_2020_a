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

#ifndef SEP_PARTICLE_BOUNDARY_COND_NONE_H
#define SEP_PARTICLE_BOUNDARY_COND_NONE_H

#include <cstdlib>
#include <climits>
#include <stdint.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_boundary_condition.h>
#include <particle_list_util.h>

#include "sep_simcontrol.h"
#include "sep_particle_definition.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class ParticleBoundaryCondNone: public ParticleBoundaryCondBase {
    public: 
      ParticleBoundaryCondNone();
      
      bool apply(pargrid::DataID particleDataID,unsigned int* N_particles,const std::vector<pargrid::CellID>& exteriorBlocks);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      
    private:
      Simulation* sim;
      SimulationClasses* simClasses;
      SPECIES species;
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondBase* ParticleBCondNoneMaker() {return new ParticleBoundaryCondNone<SPECIES,PARTICLE>();}

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondNone<SPECIES,PARTICLE>::ParticleBoundaryCondNone(): ParticleBoundaryCondBase() { }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondNone<SPECIES,PARTICLE>::apply(pargrid::DataID particleDataID,unsigned int* N_particles,
							  const std::vector<pargrid::CellID>& exteriorBlocks) {
      return true;
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondNone<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
							       const std::string& regionName,const ParticleListBase* plist) {
      bool success = ParticleBoundaryCondBase::initialize(sim,simClasses,cr,regionName,plist);
      this->sim = &sim;
      this->simClasses = &simClasses;
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      return success;
   }
   
} // namespace sep
   
#endif
