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

#ifndef PARTICLE_PROPAGATOR_RK2_GC_H
#define PARTICLE_PROPAGATOR_RK2_GC_H

#include <cstdlib>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_propagator.h>

#include "particle_definition.h"

template<class FIELD,class SPECIES,class PARTICLE>
class ParticlePropagatorRk2GC: public ParticlePropagatorBase {
 public:
   
   ParticlePropagatorRk2GC();

   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist);
   bool propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
		      const double* const coordinates,unsigned int N_particles);
     
      
 private:
   Real acc[4];
   const SPECIES* species;
   PARTICLE state1;
   
   void stepForward(pargrid::CellID block,Real t,Real dt,PARTICLE& p);
};

// Maker function for constructing a new propagator:
template<class FIELD,class SPECIES,class PARTICLE> inline
ParticlePropagatorBase* RK2GCMaker() {return new ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>();}

template<class FIELD,class SPECIES,class PARTICLE> inline
ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::ParticlePropagatorRk2GC(): ParticlePropagatorBase() { }

template<class FIELD,class SPECIES,class PARTICLE> inline
bool ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::finalize() {
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
								    const double* const coordinates,unsigned int N_particles) {
   const Real t = sim->t;
   const Real dt = sim->dt;
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   PARTICLE* particleList = wrapper.data()[blockLID];
   for (unsigned int p=0; p<N_particles; ++p) {
      Real t_current = t;
      Real dt_current = dt;
      
      #ifdef USE_INDEX_OPERATOR
         state1[particle::MU] = particleList[p][particle::MU];
      #else
         state1.state[particle::MU] = particleList[p].state[particle::MU];
      #endif
      
      // Change to absolute coordinates:
      #ifdef USE_INDEX_OPERATOR
         particleList[p][particle::XPOS] += coordinates[particle::XPOS];
         particleList[p][particle::YPOS] += coordinates[particle::YPOS];
         particleList[p][particle::ZPOS] += coordinates[particle::ZPOS];
      #else
         particleList[p].state[particle::XPOS] += coordinates[particle::XPOS];
         particleList[p].state[particle::YPOS] += coordinates[particle::YPOS];
         particleList[p].state[particle::ZPOS] += coordinates[particle::ZPOS];
      #endif
      while (t_current < t+dt) {
	 // Get maximum timestep:
         #ifdef USE_INDEX_OPERATOR
	    dt_current = FIELD::getMaximumTimestep(blockLID,t_current,dt_current,*species,&(particleList[p][0]));
	 #else
	    dt_current = FIELD::getMaximumTimestep(blockLID,t_current,dt_current,*species,particleList[p].state);
	 #endif
	 
	 // Get acceleration at (t0,r0,v0):
	 #ifdef USE_INDEX_OPERATOR
	    FIELD::getAcceleration(blockLID,t_current,0.5*dt_current,*species,&(particleList[p][0]),acc);
	 #else
	    FIELD::getAcceleration(blockLID,t_current,0.5*dt_current,*species,particleList[p].state,acc);
	 #endif

	 // Step forward in time by dt_current:
	 stepForward(blockLID,t_current,dt_current,particleList[p]);
	 
	 t_current += dt_current;
	 dt_current = t+dt-t_current;
      }
      
      // Change back to relative coordinates:
      #ifdef USE_INDEX_OPERATOR
         particleList[p][particle::XPOS] -= coordinates[particle::XPOS];
         particleList[p][particle::YPOS] -= coordinates[particle::YPOS];
         particleList[p][particle::ZPOS] -= coordinates[particle::ZPOS];
      #else
         particleList[p].state[particle::XPOS] -= coordinates[particle::XPOS];
         particleList[p].state[particle::YPOS] -= coordinates[particle::YPOS];
         particleList[p].state[particle::ZPOS] -= coordinates[particle::ZPOS];
      #endif
   }
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								 const std::string& regionName,const ParticleListBase* plist) {
   // Call base class initializer:
   bool success = ParticlePropagatorBase::initialize(sim,simClasses,cr,regionName,plist);
   
   // Get propagated particle species:
   species = reinterpret_cast<const SPECIES*>(plist->getSpecies());
   
   // Init EM field:
   if (FIELD::initialize(sim,simClasses,cr) == false) success = false;
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
void ParticlePropagatorRk2GC<FIELD,SPECIES,PARTICLE>::stepForward(pargrid::CellID block,Real t,Real dt,PARTICLE& p) {
   const Real HALF = 0.5;
   bool inirefl = false;
   Real vpar = 0.0;
   
   // Propagate to half-step t = t0 + dt/2:
   #ifdef USE_INDEX_OPERATOR
      state1[particle::XPOS] = p[particle::XPOS] + acc[particle::XPOS]*HALF*dt;
      state1[particle::YPOS] = p[particle::YPOS] + acc[particle::YPOS]*HALF*dt;
      state1[particle::ZPOS] = p[particle::ZPOS] + acc[particle::ZPOS]*HALF*dt;
      state1[particle::VPAR] = p[particle::VPAR] + acc[particle::VPAR]*HALF*dt;
      state1[particle::REFLECTED] = p[particle::REFLECTED];
   #else
      state1.state[particle::XPOS] = p.state[particle::XPOS] + acc[particle::XPOS]*HALF*dt;
      state1.state[particle::YPOS] = p.state[particle::YPOS] + acc[particle::YPOS]*HALF*dt;
      state1.state[particle::ZPOS] = p.state[particle::ZPOS] + acc[particle::ZPOS]*HALF*dt;
      state1.state[particle::VPAR] = p.state[particle::VPAR] + acc[particle::VPAR]*HALF*dt;
      state1.state[particle::REFLECTED] = p.state[particle::REFLECTED];
   #endif
   
   #ifdef USE_INDEX_OPERATOR
      if (p[particle::REFLECTED] > 0.0) {
	 inirefl = true;
	 vpar = state1[particle::VPAR];
      }
   #else
      if (p.state[particle::REFLECTED] > 0.0) {
	 inirefl = true;
	 vpar = state1.state[particle::VPAR];
      }
   #endif
   
   // Get acceleration at half-step:
   #ifdef USE_INDEX_OPERATOR
      FIELD::getAcceleration(block,t+HALF*dt,dt,*species,&(state1[0]),acc);
   #else
      FIELD::getAcceleration(block,t+HALF*dt,dt,*species,state1.state,acc);
   #endif
   
   // Propagate to t = t0 + dt:
   #ifdef USE_INDEX_OPERATOR
      p[particle::XPOS] += acc[particle::XPOS]*dt;
      p[particle::YPOS] += acc[particle::YPOS]*dt;
      p[particle::ZPOS] += acc[particle::ZPOS]*dt;
      p[particle::VPAR] += acc[particle::VPAR]*dt;
      p[particle::REFLECTED] = state1[particle::REFLECTED];
      if (inirefl == true) p[particle::VPAR] = vpar;
   #else
      p.state[particle::XPOS] += acc[particle::XPOS]*dt;
      p.state[particle::YPOS] += acc[particle::YPOS]*dt;
      p.state[particle::ZPOS] += acc[particle::ZPOS]*dt;
      p.state[particle::VPAR] += acc[particle::VPAR]*dt;
      p.state[particle::REFLECTED] = state1.state[particle::REFLECTED];
      if (inirefl == true) p.state[particle::VPAR] = vpar;
   #endif
}

#endif
