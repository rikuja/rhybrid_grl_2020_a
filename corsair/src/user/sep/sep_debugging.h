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

#ifndef SEP_DEBUGGING_H
#define SEP_DEBUGGING_H

#include <stdint.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <particle_list_base.h>

namespace sep {

   template<class SPECIES,class PARTICLE> inline
   bool checkCoordinates(Simulation& sim,SimulationClasses& simClasses,const SPECIES& species,
			 ParticleListBase* particleList,const std::string& step) {
      bool success = true;
      pargrid::DataID particleID = pargrid::INVALID_DATAID;
      if (particleList->getParticles(particleID) == false) return false;

      uint32_t* N_particles = particleList->getParticleNumberArray();
      if (N_particles == NULL) {
	 std::cerr << "ERROR: N_particles is NULL in checkCoordinates" << std::endl;
	 exit(1);
      }

      const pargrid::DataWrapper<PARTICLE>& wrapper = simClasses.pargrid.getUserDataDynamic<PARTICLE>(particleID);

      for (pargrid::CellID blockLID=0; blockLID<simClasses.pargrid.getNumberOfLocalCells(); ++blockLID) {
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);
	 
	 // Check that N_particles and wrapper.size() match:
	 if (N_particles[blockLID] != wrapper.size()[blockLID]) {
	    std::stringstream ss;
	    ss << "STEP " << sim.timestep << " t " << sim.t << " P#" << sim.mpiRank << " LID " << blockLID << " GID " << blockGID;
	    ss << " ERROR in particle lists: ";
	    ss << " N_particles " << N_particles[blockLID] << " wrapper.size() " << wrapper.size()[blockLID];
	    ss << ' ' << step << std::endl;
	    std::cerr << ss.str();
	    exit(1);
	 }

	 // Min,max acceptable x,y,z coordinate limits for block:
	 const uint32_t x_min = i_block*block::WIDTH_X;
	 const uint32_t x_max = (i_block+1)*block::WIDTH_X;
	 const uint32_t y_min = j_block*block::WIDTH_Y;
	 const uint32_t y_max = (j_block+1)*block::WIDTH_Y;
	 const uint32_t z_min = k_block*block::WIDTH_Z;
	 const uint32_t z_max = (k_block+1)*block::WIDTH_Z;
	 
	 PARTICLE* particles = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    // Check that logical coordinates are inside the block:
	    success = true;
	    if (particles[p].state[0]  < x_min) success = false;
	    if (particles[p].state[0] >= x_max) success = false;
	    if (particles[p].state[1]  < y_min) success = false;
	    if (particles[p].state[1] >= y_max) success = false;
	    if (particles[p].state[2]  < z_min) success = false;
	    if (particles[p].state[2] >= z_max) success = false;
	    
	    if (success == false) {
	       std::cerr << "STEP #" << sim.timestep << " P#" << sim.mpiRank;
	       std::cerr << " Block LID#" << blockLID << " GID#" << blockGID << " error in particle '";
	       std::cerr << species.getName() << "' coordinates" << std::endl;
               std::cerr << "\t indices: " << i_block << ' ' << j_block << ' ' << k_block << std::endl;
	       std::cerr << "\t particle coords: ";
               for (int i=0; i<3; ++i) std::cerr << particles[p].state[i] << '\t';
               std::cerr << std::endl;
	       std::cerr << "OCCURRED IN STEP '" << step << "'" << std::endl;
	       break;
            }

	    // Check offset indices, should ensure that accumulators work:
	    Real pos[3];
	    pos[0] = particles[p].state[0] - i_block*block::WIDTH_X + 1;
	    pos[1] = particles[p].state[1] - j_block*block::WIDTH_Y + 1;
	    pos[2] = particles[p].state[2] - k_block*block::WIDTH_Z + 1;
	    
	    int32_t indices[3];
	    for (int i=0; i<3; ++i) indices[i] = static_cast<int32_t>(pos[i]);
	    if (indices[0] < 1 || indices[0] > block::WIDTH_X) success = false;
	    if (indices[1] < 1 || indices[1] > block::WIDTH_Y) success = false;
	    if (indices[2] < 1 || indices[2] > block::WIDTH_Z) success = false;
	    if (success == false) {
	       std::cerr << "STEP #" << sim.timestep;
	       std::cerr << " Block LID#" << blockLID << " GID#" << blockGID << " accum indices invalid" << std::endl;
	       std::cerr << "\t indices: " << i_block << ' ' << j_block << ' ' << k_block << std::endl;
	       std::cerr << "\t coords : ";
	       for (int i=0; i<3; ++i) std::cerr << particles[p].state[i] << '\t';
	       std::cerr << std::endl;
	       std::cerr << "\t indices: " << indices[0] << ' ' << indices[1] << ' ' << indices[2] << std::endl;
	       std::cerr << "OCCURRED IN STEP '" << step << "'" << std::endl;
	       break;
	    }
	 }
      }
      
      return success;
   }
   
} // namespace sep

#endif
