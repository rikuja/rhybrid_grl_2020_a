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

#ifndef SEP_PARTICLE_BOUNDARY_COND_DOUBLE_SHOCK_H
#define SEP_PARTICLE_BOUNDARY_COND_DOUBLE_SHOCK_H

#include <cstdlib>
#include <climits>
#include <stdint.h>

#include <main.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_boundary_condition.h>
#include <particle_list_util.h>

#include "sep_simcontrol.h"
#include "sep_particle_definition.h"
#include "sep_lagr_species.h"
#include "sep_lagr_definition.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class ParticleBoundaryCondDoubleShock: public ParticleBoundaryCondBase {
    public: 
      ParticleBoundaryCondDoubleShock();
      
      bool apply(pargrid::DataID particleDataID,unsigned int* N_particles,const std::vector<pargrid::CellID>& exteriorBlocks);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      
    private:
      SPECIES species;
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondBase* PBCondDoubleShockMaker() {return new ParticleBoundaryCondDoubleShock<SPECIES,PARTICLE>();}

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondDoubleShock<SPECIES,PARTICLE>::ParticleBoundaryCondDoubleShock(): ParticleBoundaryCondBase() { }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondDoubleShock<SPECIES,PARTICLE>::apply(pargrid::DataID particleDataID,unsigned int* N_particles,
								 const std::vector<pargrid::CellID>& exteriorBlocks) {
      // Get particles:
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      
      // Get parallel and antiparallel wave packets:
      pargrid::DataID parDataID = pargrid::INVALID_DATAID;
      pargrid::DataID antiparDataID = pargrid::INVALID_DATAID;
      corsair::ObjectWrapper& objectWrapper = corsair::getObjectWrapper();

      for (size_t plist=0; plist<objectWrapper.particleLists.size(); ++plist) {
	 if (objectWrapper.particleLists[plist]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
	 const sep::LagrangianSpecies* species = reinterpret_cast<const sep::LagrangianSpecies*>(objectWrapper.particleLists[plist]->getSpecies());
	 if (species->propagationDirection < 0) objectWrapper.particleLists[plist]->getParticles(antiparDataID);
	 else objectWrapper.particleLists[plist]->getParticles(parDataID);
      }

      typedef sep::LagrangianParticle<Real> LAGR_PARTICLE;
      const pargrid::DataWrapper<LAGR_PARTICLE>& antiparWrapper = objectWrapper.simClasses.pargrid.getUserDataDynamic<LAGR_PARTICLE>(antiparDataID);
      const pargrid::DataWrapper<LAGR_PARTICLE>& parWrapper = objectWrapper.simClasses.pargrid.getUserDataDynamic<LAGR_PARTICLE>(parDataID);

      // Remove particles from all cells where there are no wave packets
      // (spectral wave intensity should be zero):
      for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 bool keepParticles = false;
	 if (parDataID != pargrid::INVALID_DATAID) {
	    if (parWrapper.size(blockLID) > 0) keepParticles = true;
	 }
	 if (antiparDataID != pargrid::INVALID_DATAID) {
	    if (antiparWrapper.size(blockLID) > 0) keepParticles = true;
	 }

	 if (keepParticles == false) {
	    wrapper.resize(blockLID,0);
	    N_particles[blockLID] = 0;
	 }
      }

      return true;
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondDoubleShock<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								      const std::string& regionName,const ParticleListBase* plist) {
      bool success = ParticleBoundaryCondBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      if (simControl.shock == NULL) {
	 simClasses.logger << "(SEP ION BCOND SHOCK) ERROR: There's no shock in simulation" << std::endl << write;
	 success = false;
      }
      
      return success;
   }
   
} // namespace sep
   
#endif
