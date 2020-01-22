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

#ifndef SEP_PARTICLE_BOUNDARY_COND_SHOCK_H
#define SEP_PARTICLE_BOUNDARY_COND_SHOCK_H

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
   class ParticleBoundaryCondShock: public ParticleBoundaryCondBase {
    public: 
      ParticleBoundaryCondShock();
      
      bool apply(pargrid::DataID particleDataID,unsigned int* N_particles,const std::vector<pargrid::CellID>& exteriorBlocks);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      
    private:
      SPECIES species;
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondBase* PBCondShockMaker() {return new ParticleBoundaryCondShock<SPECIES,PARTICLE>();}

   template<class SPECIES,class PARTICLE> inline
   ParticleBoundaryCondShock<SPECIES,PARTICLE>::ParticleBoundaryCondShock(): ParticleBoundaryCondBase() { }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondShock<SPECIES,PARTICLE>::apply(pargrid::DataID particleDataID,unsigned int* N_particles,
						      const std::vector<pargrid::CellID>& exteriorBlocks) {
      // Get particles:
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);

      // Temporary fix to periodic boundaries -- iterate over all blocks 
      // and check that particle coordinates are ok:
      uint32_t i_block,j_block,k_block;
      for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 
	 // Skip blocks that are not on mesh boundaries:
	 bool checkCoordinates = false;
	 if (i_block == 0 || i_block == sim->x_blocks-1) checkCoordinates = true;
	 if (j_block == 0 || j_block == sim->y_blocks-1) checkCoordinates = true;
	 if (k_block == 0 || k_block == sim->z_blocks-1) checkCoordinates = true;
	 if (checkCoordinates == false) continue;
	  
	 // Correct particles' coordinates if they crossed periodic boundaries:
	 PARTICLE* particles = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    // Clamp particle coordinates to be inside the block. Note that 
	    // after addition or substraction x >= i_block*block::WIDTH_X and 
	    // x < i_block*block::WIDTH_X must still apply to precision of Real
	    // (yes, there was a bug related to finite floating point precision here):
	    if (particles[p].state[0] < 0.0) {
	       particles[p].state[0] +=
		 (sim->x_blocks*block::WIDTH_X-simControl.BLOCK_EPSILON);
	    }
	    if (particles[p].state[0] >= sim->x_blocks*block::WIDTH_X) {
	       particles[p].state[0] = simControl.BLOCK_EPSILON
		 + (particles[p].state[0]-sim->x_blocks*block::WIDTH_X);
	    }
	    if (particles[p].state[1] < 0.0) {
	       particles[p].state[1] +=
		 (sim->y_blocks*block::WIDTH_Y-simControl.BLOCK_EPSILON);
	    }
	    if (particles[p].state[1] >= sim->y_blocks*block::WIDTH_Y) {
	       particles[p].state[1] = simControl.BLOCK_EPSILON
		 + (particles[p].state[1]-sim->y_blocks*block::WIDTH_Y);
	    }
	    if (particles[p].state[2] < 0.0) {
	       particles[p].state[2] +=
		 (sim->z_blocks*block::WIDTH_Z-simControl.BLOCK_EPSILON);
	    }
	    if (particles[p].state[2] >= sim->z_blocks*block::WIDTH_Z) {
	       particles[p].state[2] = simControl.BLOCK_EPSILON
		 + (particles[p].state[2]-sim->z_blocks*block::WIDTH_Z);
	    }
	    
	    #ifndef NDEBUG
	    bool ok = true;
	    if (particles[p].state[0] < 0)                            ok = false;
	    if (particles[p].state[0] > sim->x_blocks*block::WIDTH_X) ok = false;
	    if (particles[p].state[1] < 0)                            ok = false;
	    if (particles[p].state[1] > sim->y_blocks*block::WIDTH_Y) ok = false;
	    if (particles[p].state[2] < 0)                            ok = false;
	    if (particles[p].state[2] > sim->z_blocks*block::WIDTH_Z) ok = false;
	    if (ok == false) {
	       std::cerr << "PLIST error in coordinates" << std::endl;
	       std::cerr << "\t" << particles[p].state[0] << '\t' << particles[p].state[1] << '\t' << particles[p].state[2] << std::endl;
	       std::cerr << "\t min: " << i_block*block::WIDTH_X << ' ' << j_block*block::WIDTH_Y << ' ' << k_block*block::WIDTH_Z << ' ' << std::endl;
	       std::cerr << "\t max: " << (i_block+1)*block::WIDTH_X << ' ' << (j_block+1)*block::WIDTH_Y << ' ' << (k_block+1)*block::WIDTH_Z << ' ' << std::endl;
	       exit(1);
	    }
            #endif
	 }
      }

      // Remove particles that are in shock downstream region. 
      // Only do this after t >= t_shock:
      if (sim->t >= simControl.t_shock) {
      //if (simControl.includeShock == true && species.getSpeciesType() != simControl.lagrangianSpeciesTypename) {
	 std::vector<pargrid::ArraySizetype> removedParticles;
	 for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	    // Skip empty blocks:
	    if (wrapper.size(blockLID) == 0) continue;
	    
	    // Mark particles in downstream for removal, note that
	    // we are doing the checking after propagation at time t+dt:
	    removedParticles.clear();
	    PARTICLE* particles = wrapper.data()[blockLID];
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       if (simControl.shock->inDownstream(sim->t+sim->dt,particles[p].state) == true) {
		  removedParticles.push_back(p);
	       }
	    }
	    
	    particlelist::removeParticles(blockLID,removedParticles,wrapper,N_particles);

	    #ifndef NDEBUG
	    for (pargrid::ArraySizetype current=0; current<wrapper.size(blockLID); ++current) {
	       bool ok = true;
	       if (simControl.shock->inDownstream(sim->t+sim->dt,particles[current].state) == true) {
		  std::cerr << "ERROR: LID#" << blockLID << " particle in position " << current << " is in downstream, should be in upstream" << std::endl;
		  ok = false;
	       }
	       if (ok == false) exit(1);
	    }
	    #endif
	 }
      }

      // Remove particles on exterior cells:
      #warning FIXME sep_particle_list.h Lagrangian boundary conditions
      if (species.getSpeciesType() != "lagrangian") {
	 const std::vector<pargrid::CellID>& exteriorBlocks = simClasses->pargrid.getExteriorCells();
	 for (size_t c=0; c<exteriorBlocks.size(); ++c) {
	    const pargrid::CellID blockID = exteriorBlocks[c];
	    N_particles[blockID] = 0;
	    wrapper.resize(blockID,0);
	 }
      }

      return true;
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleBoundaryCondShock<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
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
