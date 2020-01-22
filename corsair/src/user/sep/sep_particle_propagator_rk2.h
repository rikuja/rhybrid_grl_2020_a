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

#ifndef SEP_PARTICLE_PROPAGATOR_RK2_H
#define SEP_PARTICLE_PROPAGATOR_RK2_H

#include <cstdlib>
#include <iostream>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <linear_algebra.h>
#include <accumulators.h>
#include <base_class_particle_propagator.h>

#include "sep_particle_definition.h"
#include "sep_simcontrol.h"

namespace sep {

   template<class SPECIES,class PARTICLE,class ACCELERATOR>
   class ParticlePropagRK2: public ParticlePropagatorBase {
    public:
      ParticlePropagRK2();
      
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
			 const double* const coordinates,unsigned int N_particles);
      
    private:
      ACCELERATOR accelerator;
      bool initialized;
      SPECIES species;
   };

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   ParticlePropagatorBase* ParticlePropagMaker() {return new ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>();}
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>::ParticlePropagRK2(): ParticlePropagatorBase() {
      initialized = false;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
      return true;
   }
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>::finalize() {
      bool success = true;
      initialized = false;
      return success;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
								       const double* const coordinates,unsigned int N_particles) {
      bool success = initialized;
      const int32_t VARIABLES = 5;
      
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE* particleList = wrapper.data()[blockLID];
      
      uint32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

      // Calculate cell indices that correctly wrap around periodic boundaries:
      int32_t i_index[block::WIDTH_X+2];
      int32_t j_index[block::WIDTH_Y+2];
      int32_t k_index[block::WIDTH_Z+2];
      
      for (int32_t i=0; i<block::WIDTH_X+2; ++i) {
	 int32_t i_fetch = i_block*block::WIDTH_X + i - 1;
	 if (i_fetch < 0)
	   i_fetch += sim->x_blocks*block::WIDTH_X;
	 else if (i_fetch >= (int32_t)sim->x_blocks*block::WIDTH_X)
	   i_fetch -= sim->x_blocks*block::WIDTH_X;
	 i_index[i] = i_fetch;
      }
      
      for (int32_t j=0; j<block::WIDTH_Y+2; ++j) {
	 int32_t j_fetch = j_block*block::WIDTH_Y + j - 1;
	 if (j_fetch < 0)
	   j_fetch += sim->y_blocks*block::WIDTH_Y;
	 else if (j_fetch >= (int32_t)sim->y_blocks*block::WIDTH_Y)
	   j_fetch -= sim->y_blocks*block::WIDTH_Y;
	 j_index[j] = j_fetch;
      }
      
      for (int32_t k=0; k<block::WIDTH_Z+2; ++k) {
	 int32_t k_fetch = k_block*block::WIDTH_Z + k - 1;
	 if (k_fetch < 0)
	   k_fetch += sim->z_blocks*block::WIDTH_Z;
	 else if (k_fetch >= (int32_t)sim->z_blocks*block::WIDTH_Z)
	   k_fetch -= sim->z_blocks*block::WIDTH_Z;
	 k_index[k] = k_fetch;
      }
      
      Real x[sep::particle::DATASIZE];
      for (int i=0; i<sep::particle::DATASIZE; ++i) x[i] = 0.0;
      Real y[VARIABLES];
      for (uint32_t p=0; p<N_particles; ++p) {
	 #ifndef NDEBUG
	 bool ok = true;
	 if (particleList[p].state[particle::XCRD] < (i_block  )*block::WIDTH_X) ok = false;
	 if (particleList[p].state[particle::XCRD] > (i_block+1)*block::WIDTH_X) ok = false;
	 if (particleList[p].state[particle::YCRD] < (j_block  )*block::WIDTH_Y) ok = false;
	 if (particleList[p].state[particle::YCRD] > (j_block+1)*block::WIDTH_Y) ok = false;
	 if (particleList[p].state[particle::ZCRD] < (k_block  )*block::WIDTH_Z) ok = false;
	 if (particleList[p].state[particle::ZCRD] > (k_block+1)*block::WIDTH_Z) ok = false;
	 if (ok == false) {
	    std::cerr << "(SEP RK2) ERROR: Invalid particle coordinates prior to propagation" << std::endl;
	    std::cerr << "\t STEP: " << sim->timestep << std::endl;
	    std::cerr << "\t BLOCK: " << i_block << ' ' << j_block << ' ' << k_block << " GID: " << blockGID << std::endl;
	    std::cerr << "\t CRDS:  " << particleList[p].state[particle::XCRD] << '\t';
	    std::cerr << particleList[p].state[particle::YCRD] << '\t' << particleList[p].state[particle::ZCRD] << std::endl;
	    exit(1);
	 }
	 #endif

	 // Calculate i,j,k indices:
	 const uint32_t I_old = static_cast<uint32_t>(particleList[p].state[particle::XCRD]);
	 const uint32_t J_old = static_cast<uint32_t>(particleList[p].state[particle::YCRD]);
	 const uint32_t K_old = static_cast<uint32_t>(particleList[p].state[particle::ZCRD]);
	 
	 #ifndef NDEBUG
	 ok = true;
	 if (I_old < i_block*block::WIDTH_X || I_old >= (i_block+1)*block::WIDTH_X) ok = false;
	 if (J_old < j_block*block::WIDTH_Y || J_old >= (j_block+1)*block::WIDTH_Y) ok = false;
	 if (K_old < k_block*block::WIDTH_Z || K_old >= (k_block+1)*block::WIDTH_Z) ok = false;
	 if (ok == false) {
	    std::cerr << "(SEP RK2) ERROR: Invalid cell indices" << std::endl;
	    std::cerr << "\t BLOCK: " << i_block << ' ' << j_block << ' ' << k_block << " GID: " << blockGID << std::endl;
	    std::cerr << "\t I,J,K: " << I_old << ' ' << J_old << ' ' << K_old << std::endl;
	    exit(1);
	 }
	 #endif

	 // Evaluate acceleration at time t:
	 //std::cerr << "INI acc" << std::endl;
	 accelerator.get(blockLID,sim->t,species,particleList[p].state,y);	 
	 for (int i=0; i<VARIABLES; ++i) y[i] *= (0.5*sim->dt); // FIXME

	 #ifndef NDEBUG
	 ok = true;
	 if (y[0]/sim->dx_cell[I_old] >= 1.0) ok = false;
	 if (y[1]/sim->dy_cell[J_old] >= 1.0) ok = false;
	 if (y[2]/sim->dz_cell[K_old] >= 1.0) ok = false;
	 if (ok == false) {
	    std::cerr << "(SEP RK2) CFL VIOLATION TRIAL STEP = " << sim->timestep << std::endl;
	    std::cerr << "\t BLOCK GID#" << blockGID << '\t' << i_block << ' ' << j_block << ' ' << k_block << std::endl;
	    std::cerr << "\t CRDS: " << particleList[p].state[particle::XCRD] << '\t' << particleList[p].state[particle::YCRD] << '\t';
	    std::cerr << particleList[p].state[particle::ZCRD] << std::endl;
	    std::cerr << "\t ACC = " << y[0]/sim->dx_cell[I_old] << '\t' << y[1]/sim->dy_cell[J_old] << '\t' << y[2]/sim->dz_cell[K_old] << std::endl;
	    exit(1);
	 }
	 #endif
	 
	 // Advance to intermediate state at time t + 0.5*dt:
	 // NOTE: coordinates need to be scaled into logical units:
	 int32_t I_star1 = I_old + static_cast<int32_t>(particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_old);
	 int32_t J_star1 = J_old + static_cast<int32_t>(particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_old);
	 int32_t K_star1 = K_old + static_cast<int32_t>(particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_old);

	 Real dx_L = particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_star1;
	 Real dy_L = particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_star1;
	 Real dz_L = particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_star1;

	 const int32_t I_dot = static_cast<int32_t>(1 + y[0]/sim->dx_cell[I_old]);
	 const int32_t J_dot = static_cast<int32_t>(1 + y[1]/sim->dy_cell[J_old]);
	 const int32_t K_dot = static_cast<int32_t>(1 + y[2]/sim->dz_cell[K_old]);
	 
	 int32_t I_off = static_cast<int32_t>(1 + particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_old);
	 int32_t J_off = static_cast<int32_t>(1 + particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_old);
	 int32_t K_off = static_cast<int32_t>(1 + particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_old);
	 x[0] = i_index[I_off] + (1-I_dot) + dx_L * sim->dx_cell[I_old]/sim->dx_cell[i_index[I_off]];
	 x[1] = j_index[J_off] + (1-J_dot) + dy_L * sim->dy_cell[J_old]/sim->dy_cell[j_index[J_off]];
	 x[2] = k_index[K_off] + (1-K_dot) + dz_L * sim->dz_cell[K_old]/sim->dz_cell[k_index[K_off]];
	 x[particle::V_PAR] = particleList[p].state[particle::V_PAR] += y[3];
	 x[particle::MU]    = particleList[p].state[particle::MU];
	 
	 // Evaluate acceleration at intermediate time step:
	 accelerator.get(blockLID,sim->t+0.5*sim->dt,species,x,y);
	 for (int i=0; i<VARIABLES; ++i) y[i] *= sim->dt;
	 
	 #ifndef NDEBUG
	 if (y[0]/sim->dx_cell[I_old] >= 1.0) ok = false;
	 if (y[1]/sim->dy_cell[J_old] >= 1.0) ok = false;
	 if (y[2]/sim->dz_cell[K_old] >= 1.0) ok = false;
	 if (ok == false) {
	    std::cerr << "(SEP RK2) CFL VIOLATION FINAL STEP = " << sim->timestep << std::endl;
	    std::cerr << "\t BLOCK GID#" << blockGID << '\t' << i_block << ' ' << j_block << ' ' << k_block << std::endl;
	    std::cerr << "\t CRDS: " << particleList[p].state[particle::XCRD] << '\t' << particleList[p].state[particle::YCRD] << '\t';
	    std::cerr << particleList[p].state[particle::ZCRD] << std::endl;
	    std::cerr << "\t ACC = " << y[0]/sim->dx_cell[I_old] << '\t' << y[1]/sim->dy_cell[J_old] << '\t' << y[2]/sim->dz_cell[K_old] << std::endl;
	    exit(1);
	 }
	 #endif

	 I_star1 = I_old + static_cast<int32_t>(particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_old);
	 J_star1 = J_old + static_cast<int32_t>(particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_old);
	 K_star1 = K_old + static_cast<int32_t>(particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_old);
	 
	 dx_L = particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_star1;
	 dy_L = particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_star1;
	 dz_L = particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_star1;
	 
	 I_off = static_cast<int32_t>(1 + particleList[p].state[particle::XCRD] + y[0]/sim->dx_cell[I_old] - I_old);
	 J_off = static_cast<int32_t>(1 + particleList[p].state[particle::YCRD] + y[1]/sim->dy_cell[J_old] - J_old);
	 K_off = static_cast<int32_t>(1 + particleList[p].state[particle::ZCRD] + y[2]/sim->dz_cell[K_old] - K_old);
	 
	 particleList[p].state[particle::XCRD]   = I_star1 + dx_L * sim->dx_cell[I_old]/sim->dx_cell[i_index[I_off]];
	 particleList[p].state[particle::YCRD]   = J_star1 + dy_L * sim->dy_cell[J_old]/sim->dy_cell[j_index[J_off]];
	 particleList[p].state[particle::ZCRD]   = K_star1 + dz_L * sim->dz_cell[K_old]/sim->dz_cell[k_index[K_off]];
	 particleList[p].state[particle::V_PAR] += y[particle::V_PAR];
	 particleList[p].state[particle::MU]    += y[particle::MU];
      }

      return success;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagRK2<SPECIES,PARTICLE,ACCELERATOR>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								    const std::string& regionName,const ParticleListBase* plist) {
      initialized = ParticlePropagatorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());

      if (accelerator.initialize(sim,simClasses,cr) == false) {
	 simClasses.logger << "(SEP RK2) ERROR: Accelerator class failed to initialize" << std::endl << write;
	 initialized = false;
      }
      
      return initialized;
   }

} // namespace sep
   
#endif
