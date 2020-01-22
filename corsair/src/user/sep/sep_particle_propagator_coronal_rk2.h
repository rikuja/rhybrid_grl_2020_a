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

#ifndef SEP_PARTICLE_PROPAGATOR_CORONAL_RK2_H
#define SEP_PARTICLE_PROPAGATOR_CORONAL_RK2_H

#include <cstdlib>
#include <iostream>
#include <climits>
#include <map>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <linear_algebra.h>
#include <accumulators.h>
#include <base_class_particle_propagator.h>

#include "sep_particle_definition.h"
#include "sep_simcontrol.h"
#include "sep_shock_drift_acceleration.h"
#include "sep_shock_accelerator.h"
#include "sep_injection_buffer.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR>
   class ParticlePropagCoronalRK2: public ParticlePropagatorBase {
    public:
      ParticlePropagCoronalRK2();
      
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
			 const double* const coordinates,unsigned int N_particles);
      bool propagateParticle(pargrid::CellID blockLID,Real t,PARTICLE& particle);
    private:
      ACCELERATOR accelerator;
      bool initialized;
      SPECIES species;

      uint32_t i_block;
      uint32_t j_block;
      uint32_t k_block;
      int32_t i_index[block::WIDTH_X+2];
      int32_t j_index[block::WIDTH_Y+2];
      int32_t k_index[block::WIDTH_Z+2];

      pargrid::CellID blockLID;
      Real t;
      Real t_upstream;
      Real dt;
      Real maxShockDistance;
      int32_t splitAmount;
      Real splitThreshold;

      bool calculateShockDriftAcceleration(pargrid::CellID blockLID,Real t_crossing,PARTICLE& particle,PARTICLE& upstream,const Real d_shock);
      void binarySearchShockCrossing(PARTICLE& particle,PARTICLE& upstream,Real* x,Real* y);
      void propagateParticle(PARTICLE& particle,Real* x,Real* y);
   };

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   ParticlePropagatorBase* ParticlePropagCoronalMaker() {return new ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>();}
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::ParticlePropagCoronalRK2(): ParticlePropagatorBase() {
      initialized = false;
      splitAmount = 1;
      splitThreshold = 1.0;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
      return true;
   }
   
   /** Search approximate time of shock crossing. After this function returns, current 
    * simulation time for the particle, given by ParticlePropagCoronalRK2::t, is set to 
    * the approximate time of shock crossing, and ParticlePropagCoronalRK2::dt is set to 
    * sim->dt - dt, where dt is the timestep taken by this function.
    */
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   void ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::binarySearchShockCrossing(PARTICLE& initialState,PARTICLE& upstream,Real* x,Real* y) {
      // Shock crossing is in the interval [t,t+dt], where dt=sim->dt. 
      // Current simulation time is t=sim->t:
      Real t_min = t;
      Real t_max = t+dt;
      t_upstream = t;

      // Make a copy of particle:
      PARTICLE particle = initialState;
      upstream = initialState;

      #ifndef NDEBUG
         int iters=0;
      #endif
      
      // Binary search time of shock crossing. After do-while loop finishes,
      // approximate time of shock crossing is 0.5*(t_min+t_max);
      const Real d_shock_initial = simControl.shock->getSquaredDistanceToShock(t,particle.state);
      Real d_shock2 = d_shock_initial;
      do {
	 // Calculate time and timestep:
	 dt = 0.5*(t_min + t_max) - t;

	 // Propagate:
	 propagateParticle(particle,x,y);

	 // Calculate squared distance to shock:
	 d_shock2 = simControl.shock->getSquaredDistanceToShock(t+dt,particle.state);

	 if (d_shock_initial >= 0.0) {
	    if (d_shock2 < -maxShockDistance) {
	       // Current timestep is too large, reduce t_max:
	       t_max = 0.5*(t_min+t_max);
	    } else if (d_shock2 > 0.0) {
	       // Current timestep is too small, distance to shock did not 
	       // change sign (particle did not cross shock), increase t_min:
	       upstream = particle;
	       t_upstream = t+dt;
	       t_min = 0.5*(t_min+t_max);	       
	    } else {
	       // Current timestep is within accepted limit, exit:
	       break;
	    }
	 } else {
	    if (d_shock2 > maxShockDistance) {
	       // Current timestep is too large, reduce t_max:
	       t_max = 0.5*(t_min+t_max);
	    } else if (d_shock2 < 0.0) {
	       // Current timestep is too small, distance to shock did not
	       // change sign (particle did not cross shock), increase t_min:
	       t_min = 0.5*(t_min+t_max);
	    } else {
	       // Current timestep is within accepted limit, exit:
	       break;
	    }
	 }
	 
	 #ifndef NDEBUG
	 ++iters;
	 if (iters > 1000) {
	    std::cerr << "(SEP PARTICLE PROPAG SHOCK) Too many binary search iterations" << std::endl;
	    std::cerr << "\t t=" << t << std::endl;
	    std::cerr << "\t ini coords: " << initialState.state[0] << '\t' << initialState.state[1];
	    std::cerr << '\t' << initialState.state[2] << std::endl;
	    std::cerr << "\t cur coords: " << particle.state[0] << '\t' << particle.state[1] << '\t' << particle.state[2] << std::endl;
	    std::cerr << "\t V_par     : " << particle.state[sep::particle::V_PAR] << std::endl;
	    std::cerr << "\t ini dist is: " << d_shock_initial << std::endl;
	    std::cerr << "\t cur dist is: " << d_shock2;
	    std::cerr << " max shock distance is: " << maxShockDistance << std::endl;
	    exit(1);
	 }
	 #endif
	 
	 // Revert back to initial state:
	 particle = initialState;
      } while (true);

      this->t  += dt;
      this->dt = sim->dt - this->dt;
      initialState = particle;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::calculateShockDriftAcceleration(pargrid::CellID blockLID,Real t_crossing,
												PARTICLE& particle,PARTICLE& upstream,
												const Real d_shock) {
      // Get magnetic field at particle position:
      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(blockLID,t_crossing,particle.state,E,B,gradB);
      
      // Get plasma state at particle position:
      #ifndef NDEBUG
         if (simControl.fieldsGetPlasmaState == NULL) {
	    std::cerr << "getPlasmaState is NULL" << std::endl;
	    exit(1);
	 }
      #endif
      PlasmaState plasmaState;
      (*simControl.fieldsGetPlasmaState)(blockLID,t_crossing,particle.state,plasmaState);

      Real V_wave[3];
      Real dV_wave,V1_alfven;
      (*simControl.fieldsGetState)(blockLID,t_crossing,particle.state,B,V_wave,dV_wave,V1_alfven,+1);
      
      // Get shock normal:
      Real shockNormal[3];
      simControl.shock->getShockNormal(t,particle.state,shockNormal);

      // Get local shock velocity in simulation frame:
      Real V_shock_SIM[3];
      simControl.shock->getLocalShockVelocity(t,particle.state,V_shock_SIM);

      // Get transformation velocity to local SNIF:
      Real V_SIM_2_SNIF[3];
      getTransformSimToLocalSNIF(plasmaState.V_plasma_SIM,V_shock_SIM,shockNormal,V_SIM_2_SNIF);

      // Calculate parallel component of transform velocity:
      Real V_trans_par = dotProduct<3>(V_SIM_2_SNIF,B) / vectorMagnitude<3>(B);
      
      // Calculate plasma velocity in SNIF:
      Real V_plasma_SNIF[3];
      for (int i=0; i<3; ++i) V_plasma_SNIF[i] = plasmaState.V_plasma_SIM[i] - V_SIM_2_SNIF[i];
      
      // Calculate tangential component of B:
      Real B_tang[3];
      const Real B_norm = dotProduct<3>(B,shockNormal);
      for (int i=0; i<3; ++i) B_tang[i] = B[i] - B_norm*shockNormal[i];

      PlasmaParameters plasma;
      plasma.B1_norm        = -B_norm;
      plasma.B1_tang        = vectorMagnitude<3>(B_tang);
      plasma.V1_plasma_norm = -dotProduct<3>(V_plasma_SNIF,shockNormal);
      plasma.V1_plasma_tang = 0.0;
      plasma.R_gas          = simControl.shock->getGasCompressionRatio(blockLID,t,particle.state);
      plasma.R_magn         = simControl.shock->getMagneticCompressionRatio(blockLID,t,particle.state);
      plasma.L_shock        = 1000.0;

      // Shock crossing algorithm (binarySearchShockCrossing) propagates 
      // particle to the opposite side of shock, i.e., if particle started
      // out in upstream region (d_shock>0) it is propagated to downstream 
      // position (d_shock<0). So, we have to switch dowstream position to 
      // upstream position and vice versa before calling shock accelerator:
      ParticleParameters p;
      if (d_shock >= 0.0) p.state[shockaccelerator::XPOS] = +0.5*plasma.L_shock;
      else p.state[shockaccelerator::XPOS] = -0.5*plasma.L_shock;
      p.state[shockaccelerator::YPOS]  = 0.0;
      p.state[shockaccelerator::ZPOS]  = 0.0;
      p.state[shockaccelerator::V_PAR] = particle.state[sep::particle::V_PAR] - V_trans_par;
      p.state[shockaccelerator::MU]    = particle.state[sep::particle::MU];

      // Solve shock encounter:
      simControl.shockAccelerator.setSpecies(species);
      Real dummy;
      Real value = simControl.shockAccelerator.getReturnedParticle(plasma,p,V1_alfven,dummy);

      // If particle reflected or was returned back to upstream, 
      // set particle position to pre-encounter position, and 
      // reset time and time step to pre-encounter values:
      if (p.state[shockaccelerator::XPOS] < 0.0) {
	 particle.state[sep::particle::XCRD] = upstream.state[sep::particle::XCRD];
	 particle.state[sep::particle::YCRD] = upstream.state[sep::particle::YCRD];
	 particle.state[sep::particle::ZCRD] = upstream.state[sep::particle::ZCRD];
	 t = t_upstream;
	 dt = sim->dt - (t_upstream-sim->t);
      }

      // Set parallel speed and magnetic moment to post-encounter values:
      particle.state[sep::particle::V_PAR] = p.state[shockaccelerator::V_PAR] + V_trans_par;
      particle.state[sep::particle::MU]    = p.state[shockaccelerator::MU];
      particle.state[sep::particle::WEIGHT] *= value;

      #ifndef NDEBUG
         if (dt < 0 || sim->dt - dt < -1e-4) {
	    std::cerr << "(SEP PARTICLE PROPAG SHOCK) ERROR: current dt exceeds sim.dt!" << std::endl;
	    std::cerr << "\t sim.dt: " << sim->dt << "\t current dt: " << dt << std::endl;
	    exit(1);
	 }
      #endif

      if (p.state[shockaccelerator::XPOS] < 0.0) {
	 return true;
      }
      return false;
   }
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::finalize() {
      bool success = true;
      initialized = false;
      return success;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
									   const std::string& regionName,const ParticleListBase* plist) {
      // Initialize propagator base class:
      initialized = ParticlePropagatorBase::initialize(sim,simClasses,cr,regionName,plist);
      
      // Make a copy of propagated particle species:
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());

      // Init accelerator:
      if (accelerator.initialize(sim,simClasses,cr) == false) {
	 simClasses.logger << "(SEP PARTICLE SHOCK RK2) ERROR: Accelerator class failed to initialize" << std::endl << write;
	 initialized = false;
      }

      // Check that shock exists in simulation:
      if (simControl.shock == NULL) {
	 simClasses.logger << "(SEP PARTICLE SHOCK RK2) ERROR: Shock class is NULL" << std::endl << write;
	 initialized = false;
      } else {
	 maxShockDistance = simControl.shock->getAllowedDistance();
      }

      return initialized;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
									    const double* const coordinates,unsigned int N_particles) {
      bool success = initialized;
      
      // Create injection buffer:
      extern std::map<std::string,InjectionBuffer<PARTICLE> > particleInjectionBuffers;
      particleInjectionBuffers[species.getName()];
      typename std::map<std::string,InjectionBuffer<PARTICLE> >::iterator it
	= particleInjectionBuffers.find(species.getName());
      InjectionBuffer<PARTICLE>* injBuffer = &(it->second);
      
      #ifndef NDEBUG
         if (it == particleInjectionBuffers.end()) {
	    std::cerr << "(SEP) Particle propag coronal rk2 received NULL injection buffer" << std::endl;
	    exit(1);
	 }
         if (injBuffer == NULL) {
	    std::cerr << "(SEP) Particle propag coronal rk2 received NULL injection buffer" << std::endl;
	    exit(1);
	 }
      #endif

      const pargrid::CellID* const RESTRICT nbrIDs = simClasses->pargrid.getCellNeighbourIDs(blockLID);

      // Get array containing particles:
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE* particleList = wrapper.data()[blockLID];

      // Calculate block indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      // Calculate cell indices that correctly wrap around periodic boundaries:
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
      
      // Propagate all particles in block:
      this->blockLID = blockLID;
      Real x[sep::particle::DATASIZE];
      Real y[5];
      PARTICLE upstream;
      for (uint32_t p=0; p<N_particles; ++p) {
	 t  = sim->t;
	 dt = sim->dt;

	 Real E[3];
	 Real B[3];
	 Real gradB[9];

	 // Take a trial step of length dt. We need to calculate particle's 
	 // distance to shock here so that we can detect shock crossing:
	 const int32_t initialRegion = simControl.shock->getShockRegion(t,particleList[p].state);	 
	 const PARTICLE initialState = particleList[p];

	 propagateParticle(particleList[p],x,y);
	 if (sim->t+sim->dt < simControl.t_shock) continue;

	 // Calculate distance to shock at trial position. If d_shock changes 
	 // sign (d_shock2_initial*d_shock < 0.0) particle crossed the shock and 
	 // we need to apply shock drift acceleration:
	 const int32_t trialRegion = simControl.shock->getShockRegion(t+dt,particleList[p].state);
	 if (trialRegion != initialRegion) {
	    // Revert back to initial state:
	    particleList[p] = initialState;

	    // Find time of shock crossing. This will propagate particle 
	    // to t=t_crossing and set dt = dt - dt_taken, where dt_taken 
	    // is the time step taken by the binary search. The particle will 
	    // end up on the opposite side of the shock:
	    binarySearchShockCrossing(particleList[p],upstream,x,y);

	    // Calculate particle energy prior to shock encounter.
	    (*simControl.fieldsGetFields)(blockLID,t,particleList[p].state,E,B,gradB);
	    Real B_mag = vectorMagnitude<3>(B);
	    const Real U_old = B_mag*particleList[p].state[particle::MU] 
	      + 0.5*species.mass*particleList[p].state[particle::V_PAR]*particleList[p].state[particle::V_PAR];

	    // Apply shock-drift acceleration. Only modifies particle's parallel speed:
	    const Real d_shock2_initial = simControl.shock->getSquaredDistanceToShock(t,particleList[p].state);
	    const bool reflected = calculateShockDriftAcceleration(blockLID,t,particleList[p],upstream,d_shock2_initial);
	    
	    // Propagate to time t+dt:
	    propagateParticle(particleList[p],x,y);

	    // Calculate particle energy after shock encounter.
	    (*simControl.fieldsGetFields)(blockLID,t,particleList[p].state,E,B,gradB);
	    B_mag = vectorMagnitude<3>(B);
	    const Real U_new = B_mag*particleList[p].state[particle::MU]
	      + 0.5*species.mass*particleList[p].state[particle::V_PAR]*particleList[p].state[particle::V_PAR];	    

	    // Calculate energy gain. If it exceeds threshold value, split 
	    // particle into N new particles at same position:
	    const Real energyGain = (U_new-U_old)/U_old;

	    if (energyGain >= splitThreshold && reflected == true) {	       
	       // Calculate particle's i,j,k offset indices in block-based mesh:
	       const int32_t I_off = static_cast<uint32_t>(1 + (particleList[p].state[0]-i_block*block::WIDTH_X) / block::WIDTH_X) - 1;
	       const int32_t J_off = static_cast<uint32_t>(1 + (particleList[p].state[1]-j_block*block::WIDTH_Y) / block::WIDTH_Y) - 1;
	       const int32_t K_off = static_cast<uint32_t>(1 + (particleList[p].state[2]-k_block*block::WIDTH_Z) / block::WIDTH_Z) - 1;

	       // Insert splitted particles directly to target block:
	       const pargrid::NeighbourID nbrOffset = simClasses->pargrid.calcNeighbourTypeID(I_off,J_off,K_off);
	       const pargrid::CellID newCell = nbrIDs[nbrOffset];
	       
	       if (newCell != pargrid::INVALID_CELLID) {
		  particleList[p].state[particle::WEIGHT] /= splitAmount;
		  for (int cntr=1; cntr<splitAmount; ++cntr) injBuffer->insert(newCell,particleList[p]);
	       }
	    }
	 }
      }

      return success;
   }
   
   /**
    *
    * @return If true, particle encountered shock.*/
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateParticle(pargrid::CellID blockLID,Real t,PARTICLE& particle) {
      std::cerr << "ERROR" << std::endl;
      exit(1);
      
      // Calculate block indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      // Calculate cell indices that correctly wrap around periodic boundaries:
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
      
      this->blockLID = blockLID;
      Real x[sep::particle::DATASIZE];
      Real y[sep::particle::DATASIZE];
      
      this->t = t;
      dt = sim->dt;
      
      // Take a trial step of length dt. We need to calculate particle's 
      // distance to shock here so that we can detect shock crossing:
      const int32_t initialRegion = simControl.shock->getShockRegion(t,particle.state);
      const PARTICLE initialState = particle;
      
      propagateParticle(particle,x,y);

      // Calculate distance to shock at trial position. If d_shock changes 
      // sign (d_shock2_initial*d_shock < 0.0) particle crossed the shock and 
      // we need to apply shock drift acceleration:
      const int32_t trialRegion = simControl.shock->getShockRegion(t+dt,particle.state);
      if (trialRegion != initialRegion) {
	 // Revert back to initial state:
	 this->t = t;
	 dt      = sim->dt;
	 particle = initialState;
	 
	 // Find time of shock crossing. This will propagate particle 
	 // to t=t_crossing and set dt = dt - dt_taken, where dt_taken 
	 // is the time step taken by the binary search. The particle will 
	 // end up on the opposite side of the shock:
	 PARTICLE upstream;
	 binarySearchShockCrossing(particle,upstream,x,y);

	 // Apply shock-drift acceleration. Only modifies particle's parallel speed:
	 const Real d_shock2_initial = simControl.shock->getSquaredDistanceToShock(this->t,particle.state);
	 calculateShockDriftAcceleration(blockLID,this->t,particle,upstream,d_shock2_initial);

	 // Propagate to time t+dt and return value 'true' if 
	 // particle reflected off the shock:
	 propagateParticle(particle,x,y);

	 if (simControl.shock->getSquaredDistanceToShock(this->t+this->dt,particle.state) >= 0.0) return true;
	 return false;
      }
      
      // No shock encounter:
      return false;
   }
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   void ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateParticle(PARTICLE& particle,Real* x,Real* y) {
      // This clear is unnecessary but speeds up code by a tiny bit:
      for (int i=0; i<sep::particle::DATASIZE; ++i) x[i] = 0.0;
      for (int i=0; i<sep::particle::DATASIZE; ++i) y[i] = 0.0;

      // Calculate correct particle position to temporary array defined below.
      // This is needed because shock drift acceleration causes propagator to 
      // be called twice for some particles, and cell indices will be incorrect 
      // for some of these particles due to periodic boundaries:
      Real correctState[sep::particle::DATASIZE];
      int32_t I_fake = static_cast<int32_t>(particle.state[particle::XCRD]);
      int32_t J_fake = static_cast<int32_t>(particle.state[particle::YCRD]);
      int32_t K_fake = static_cast<int32_t>(particle.state[particle::ZCRD]);
      int32_t I_offset = I_fake - i_block*block::WIDTH_X + 1;
      int32_t J_offset = J_fake - j_block*block::WIDTH_Y + 1;
      int32_t K_offset = K_fake - k_block*block::WIDTH_Z + 1;
      correctState[0] = i_index[I_offset] + (particle.state[0]-I_fake);
      correctState[1] = j_index[J_offset] + (particle.state[1]-J_fake);
      correctState[2] = k_index[K_offset] + (particle.state[2]-K_fake);
      for (int i=3; i<sep::particle::DATASIZE; ++i) correctState[i] = particle.state[i];
      
      // Evaluate acceleration at time t:
      //accelerator.get(blockLID,t,species,particle.state,y);
      accelerator.get(blockLID,t,species,correctState,y);
      for (int i=0; i<sep::particle::DATASIZE; ++i) y[i] *= (0.5*dt);

      #ifndef NDEBUG
      bool violation = false;
      if (y[0]/sim->dx_cell[0] > 0.5) violation = true;
      if (y[1]/sim->dy_cell[0] > 0.5) violation = true;
      if (y[2]/sim->dz_cell[0] > 0.5) violation = true;
      if (violation == true) {
	 std::cerr << "(SEP RK2 SHOCK) CFL condition violated at trial step!" << std::endl;
	 std::cerr << "                " << y[0]/sim->dx_cell[0] << '\t' << y[1]/sim->dy_cell[0] << '\t' << y[2]/sim->dz_cell[0] << std::endl;
	 std::cerr << "                V_PAR: " << correctState[sep::particle::V_PAR] << std::endl;
	 std::cerr << "                I,J,K: " << i_index[I_offset] << ' ' << j_index[J_offset] << ' ' << k_index[K_offset] << std::endl;
	 std::cerr << "                x,y,z: " << particle.state[particle::XCRD] << '\t' << particle.state[particle::YCRD] << '\t' << particle.state[particle::ZCRD] << std::endl;
	 std::cerr << "                dx,dy,dz: " << sim->dx_cell[0] << '\t' << sim->dy_cell[0] << '\t' << sim->dz_cell[0] << std::endl;
	 std::cerr << "                dt,dt_sim: " << dt << '\t' << sim->dt << std::endl;
	 exit(1);
      }
      #endif

      x[0]               = correctState[0] + y[0] / sim->dx_cell[0];
      x[1]               = correctState[1] + y[1] / sim->dy_cell[0];
      x[2]               = correctState[2] + y[2] / sim->dz_cell[0];
      x[particle::V_PAR] = correctState[particle::V_PAR] + y[3];
      for (int i=4; i<sep::particle::DATASIZE; ++i) x[i] = correctState[i];
      
      // Evaluate acceleration at intermediate time step:
      accelerator.get(blockLID,t+0.5*dt,species,x,y);
      for (int i=0; i<sep::particle::DATASIZE; ++i) y[i] *= dt;

      #ifndef NDEBUG
      violation = false;
      if (y[0]/sim->dx_cell[0] > 1.0) violation = true;
      if (y[1]/sim->dy_cell[0] > 1.0) violation = true;
      if (y[2]/sim->dz_cell[0] > 1.0) violation = true;
      if (violation == true) {
	 std::cerr << "(SEP RK2 SHOCK) CFL condition violated at propagation step!" << std::endl;
	 std::cerr << "                " << y[0]/sim->dx_cell[0] << '\t' << y[1]/sim->dy_cell[0] << '\t' << y[2]/sim->dz_cell[0] << std::endl;
	 exit(1);
      }
      #endif
      
      particle.state[0] += y[0] / sim->dx_cell[0];
      particle.state[1] += y[1] / sim->dy_cell[0];
      particle.state[2] += y[2] / sim->dz_cell[0];
      particle.state[particle::V_PAR] += y[particle::V_PAR];
      particle.state[particle::MU]    += y[particle::MU];
      for (int i=5; i<sep::particle::DATASIZE; ++i) particle.state[i] = x[i];
   }

} // namespace sep
   
#endif
