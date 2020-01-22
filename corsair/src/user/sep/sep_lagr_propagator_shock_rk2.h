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

#ifndef SEP_LAGR_PROPAGATOR_SHOCK_RK2_H
#define SEP_LAGR_PROPAGATOR_SHOCK_RK2_H

#include <cstdlib>
#include <iostream>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <linear_algebra.h>
#include <accumulators.h>
#include <base_class_particle_propagator.h>

#include "sep_lagr_definition.h"
#include "sep_simcontrol.h"
#include "sep_injection_buffer.h"

namespace sep {
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR>
   class LagrangianPropagShockRK2: public ParticlePropagatorBase {
    public:
      LagrangianPropagShockRK2();

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
      
      uint32_t i_block;
      uint32_t j_block;
      uint32_t k_block;
      uint32_t i_index[block::WIDTH_X+2];
      uint32_t j_index[block::WIDTH_Y+2];
      uint32_t k_index[block::WIDTH_Z+2];

      Real t;
      Real dt;
      pargrid::CellID blockLID;
      Real maxShockDistance;

      void applyShockBoundaryConditions(Real t,const PARTICLE& incident,PARTICLE& reflected,PARTICLE& transmitted);
      void binarySearchShockCrossing(int propagationDirection,PARTICLE& initialState,PARTICLE& opposite);
      void propagateParticle(PARTICLE& particle,int propagationDirection,const Real& t,const Real& dt);
   };

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   ParticlePropagatorBase* LagrPropagShockMaker() {return new LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>();}
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::LagrangianPropagShockRK2(): ParticlePropagatorBase() {
      initialized = false;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
      return true;
   }

   /**
    * 
    * NOTE: incident state is defined at time t=sim->t, reflected and transmitted
    * states are defined at t=t_crossing.
    */
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   void LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::applyShockBoundaryConditions(Real t_crossing,
											     const PARTICLE& incident,
											     PARTICLE& reflected,
											     PARTICLE& transmitted) {
      // Calculate correct transmitted wave packet position:
      const int32_t VARIABLES=4;
      Real correctTransmitted[VARIABLES];
      const int32_t I_fake1 = static_cast<int32_t>(transmitted.state[lagr::XCRD]);
      const int32_t J_fake1 = static_cast<int32_t>(transmitted.state[lagr::YCRD]);
      const int32_t K_fake1 = static_cast<int32_t>(transmitted.state[lagr::ZCRD]);
      const int32_t I_offset1 = I_fake1 - i_block*block::WIDTH_X + 1;
      const int32_t J_offset1 = J_fake1 - j_block*block::WIDTH_Y + 1;
      const int32_t K_offset1 = K_fake1 - k_block*block::WIDTH_Z + 1;
      correctTransmitted[0] = i_index[I_offset1] + (transmitted.state[0]-I_fake1);
      correctTransmitted[1] = j_index[J_offset1] + (transmitted.state[1]-J_fake1);
      correctTransmitted[2] = k_index[K_offset1] + (transmitted.state[2]-K_fake1);

      // Get wave speeds of reflected and transmitted waves at reflected position:
      Real B2[3];
      Real V_wave_refl[3];
      Real V_wave_trans[3];
      Real dV_wave,V_alfven;
      (*simControl.fieldsGetState)(blockLID,t_crossing,correctTransmitted,B2,V_wave_refl ,dV_wave,V_alfven,-1*species.propagationDirection);
      (*simControl.fieldsGetState)(blockLID,t_crossing,correctTransmitted,B2,V_wave_trans,dV_wave,V_alfven,   species.propagationDirection);
      
      // Calculate correct indicent wave packet position:
      const int32_t I_fake2 = static_cast<int32_t>(incident.state[lagr::XCRD]);
      const int32_t J_fake2 = static_cast<int32_t>(incident.state[lagr::YCRD]);
      const int32_t K_fake2 = static_cast<int32_t>(incident.state[lagr::ZCRD]);
      const int32_t I_offset2 = I_fake2 - i_block*block::WIDTH_X + 1;
      const int32_t J_offset2 = J_fake2 - j_block*block::WIDTH_Y + 1;
      const int32_t K_offset2 = K_fake2 - k_block*block::WIDTH_Z + 1;
      Real incidentPos[3];
      incidentPos[0] = i_index[I_offset2] + (incident.state[0]-I_fake2);
      incidentPos[1] = j_index[J_offset2] + (incident.state[1]-J_fake2);
      incidentPos[2] = k_index[K_offset2] + (incident.state[2]-K_fake2);

      // Get gas compression ratio (this needs to be done on upstream side):
      #warning Intentionally passing invalid block LID, should work for now
      const pargrid::CellID blockLID = pargrid::INVALID_CELLID;
      const Real sqrt_R_gas = sqrt(simControl.shock->getGasCompressionRatio(blockLID,sim->t,incidentPos));

      // Get plasma state on upstream side (needed for reflection and transmission coefficients):
      PlasmaState plasmaState;
      (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,incidentPos,plasmaState);
      
      // Get wave speed at incident position:
      Real B1[3];
      Real V_wave1[3];
      (*simControl.fieldsGetState)(blockLID,sim->t,incidentPos,B1,V_wave1,dV_wave,V_alfven,species.propagationDirection);

      // Convert incident wavelength to physical units:
      Real lambdaIncident = (*simControl.getPhysicalWavelength)(incident.state[lagr::LAMBDA]);

      // Frequency (in local shock frame) is same for incident, reflected, and transmitted waves.
      // However, wave speeds differ, so there is a jump in wavelengths:
      Real V_shock[3];
      simControl.shock->getLocalShockVelocity(t_crossing,transmitted.state,V_shock);
      const Real lambdaReflected   = (V_wave_refl[0] - V_shock[0])  / (V_wave1[0] - V_shock[0]) * lambdaIncident;
      const Real lambdaTransmitted = (V_wave_trans[0] - V_shock[0]) / (V_wave1[0] - V_shock[0]) * lambdaIncident;
      reflected.state[lagr::LAMBDA]   = lambdaReflected;
      transmitted.state[lagr::LAMBDA] = lambdaTransmitted;

      //std::cerr << V_wave1[0] - V_shock[0] << '\t' << V_wave_trans[0] - V_shock[0] << '\t' << V_wave_refl[0] - V_shock[0] << std::endl;
      //std::cerr << V_wave1[0]/V_wave_refl[0] << '\t' << V_wave1[0]/V_wave_trans[0] << '\t' << V_shock[0] << std::endl;
      
      // Convert transmitted wavelength to logical units:
      transmitted.state[lagr::LAMBDA] = (*simControl.getLogicalWavelength)(transmitted.state[lagr::LAMBDA]);

      // Convert reflected wavelength to logical units:
      reflected.state[lagr::XCRD] = transmitted.state[lagr::XCRD];
      reflected.state[lagr::YCRD] = transmitted.state[lagr::YCRD];
      reflected.state[lagr::ZCRD] = transmitted.state[lagr::ZCRD];
      reflected.state[lagr::LAMBDA] = (*simControl.getLogicalWavelength)(reflected.state[lagr::LAMBDA]);
      
      // Calculate reflection and transmission coefficients:
      const Real V_alfvenn = sqrt(plasmaState.alfvenSpeed2);
      const Real Mach = fabs(plasmaState.V_plasma_SIM[0]) / V_alfvenn;
      Real coeffTransmission,coeffReflection;
      if (lambdaIncident*lambdaTransmitted > 0) {
	 //coeffTransmission = 0.5*sqrt_R_gas*(sqrt_R_gas+1)*(Mach + 1)/(Mach + sqrt_R_gas);
	 //coeffReflection   = 0.5*sqrt_R_gas*(sqrt_R_gas-1)*(Mach + 1)/(Mach - sqrt_R_gas);
	 coeffTransmission = (sqrt_R_gas+1)/2/sqrt_R_gas*lambdaIncident/lambdaTransmitted;
	 coeffReflection   = (sqrt_R_gas-1)/2/sqrt_R_gas*lambdaIncident/lambdaReflected;
      } else {
	 coeffTransmission = 0.0;
	 coeffReflection = 0.0;
      }

      // Calculate reflected and transmitted wave energies:
      reflected.state[lagr::ENERGY] = transmitted.state[lagr::ENERGY] 
	                            * coeffReflection*coeffReflection 
	                            * lambdaReflected/(lambdaIncident + 1e-10);
      transmitted.state[lagr::ENERGY] = transmitted.state[lagr::ENERGY]
	                              * coeffTransmission*coeffTransmission 
	                              * lambdaTransmitted/(lambdaIncident + 1e-10);

      #ifndef NDEBUG
      if (reflected.state[lagr::ENERGY] != reflected.state[lagr::ENERGY]) {
	 std::cerr << "NAN R=" << coeffReflection << '\t' << lambdaReflected << '\t' << lambdaIncident << '\t' << incident.state[lagr::LAMBDA] << std::endl;
	 exit(1);
      }
      
      if (transmitted.state[lagr::ENERGY] != transmitted.state[lagr::ENERGY]) {
	 std::cerr << "NAN T=" << coeffTransmission << '\t' << lambdaTransmitted << '\t' << lambdaIncident << std::endl;
	 exit(1);
      }
      #endif
   }

   /** Search approximate time of shock crossing. After function returns incident contains 
    * a particle that is on the initial side of the shock, and transmitted is on the opposite 
    * side of the shock.
    * @param incident Incident particle.
    * @param transmitted Transmitted particle.
    */
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   void LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::binarySearchShockCrossing(int propagationDirection,PARTICLE& incident,PARTICLE& transmitted) {
      // Shock crossing is in the interval [t,t+dt], where dt=sim->dt. 
      // Current simulation time is t=sim->t:
      Real t_min = sim->t;
      Real t_max = sim->t+sim->dt;
      
      // Make a copy of particle:
      const PARTICLE initialState = incident;
      const Real d_shock_initial = simControl.shock->getSquaredDistanceToShock(t,initialState.state);
      
      // Binary search time of shock crossing. After do-while loop finishes,
      // approximate time of shock crossing is 0.5*(t_min+t_max);
      Real d_shock2 = d_shock_initial;
      int iters = 0;
      do {
	 // Calculate time and timestep:
	 dt = 0.5*(t_min + t_max) - sim->t;
	 
	 // Propagate:
	 propagateParticle(transmitted,propagationDirection,sim->t,dt);

	 // Calculate squared distance to shock:
	 d_shock2 = simControl.shock->getSquaredDistanceToShock(sim->t+dt,transmitted.state);
	 
	 if (d_shock_initial >= 0.0) {
	    if (d_shock2 < -maxShockDistance) {
	       // Current timestep is too large, reduce t_max:
	       t_max = 0.5*(t_min+t_max);
	    } else if (d_shock2 > 0.0) {
	       // Current timestep is too small, distance to shock did not 
	       // change sign (particle did not cross shock), increase t_min:
	       t_min = 0.5*(t_min+t_max);
	       incident = transmitted;
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
	       incident = transmitted;
	    } else {
	       // Current timestep is within accepted limit, exit:
	       break;
	    }
	 }
	 
	 if (iters == 60) {
	    std::cerr << "(SEP LAGR PROPAGATOR SHOCK RK2) ERROR: Too many binary search iterations" << std::endl;
	    std::cerr << "\t x_incident    " << initialState.state[0] << "\t d_shock(initial) " << d_shock_initial << std::endl;
	    std::cerr << "\t x_transmitted " << transmitted.state[0] << "\t d_shock(current) " << d_shock2 << std::endl;
	    std::cerr << "\t sim.t " << sim->t << "\t t_search " << sim->t+dt << "\t dt_sim " << sim->dt << std::endl;
	    std::cerr << "\t max allowed distance to shock is " << maxShockDistance << std::endl;
	    exit(1);
	 }
	 
	 // Revert back to initial state:
	 transmitted = initialState;	 
	 ++iters;
      } while (true);

      t  = sim->t + dt;
      dt = std::max(0.0,sim->dt - dt);
   }
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::finalize() {
      bool success = true;
      initialized = false;
      return success;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateCell(pargrid::CellID blockLID,pargrid::DataID particleDataID,
									 const double* const coordinates,unsigned int N_particles) {
      bool success = initialized;
      
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE* particleList = wrapper.data()[blockLID];
      
      // Get Alfven waves travelling to opposite direction, needed for wave reflection below:
      sep::InjectionBuffer<PARTICLE>* injBuffer = NULL;
      extern sep::InjectionBuffer<PARTICLE> antiparInjectionBuffer;
      extern sep::InjectionBuffer<PARTICLE> parInjectionBuffer;
      if (species.propagationDirection < 0) {
	 injBuffer = &parInjectionBuffer;
      } else {
	 injBuffer = &antiparInjectionBuffer;
      }

      const pargrid::CellID* const RESTRICT nbrIDs = simClasses->pargrid.getCellNeighbourIDs(blockLID);

      this->blockLID = blockLID;
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
      
      for (uint32_t p=0; p<N_particles; ++p) {
	 // Take a trial step of length dt. We need to calculate particle's 
	 // distance to shock here so that we can detect shock crossing:
	 t  = sim->t;
	 dt = sim->dt;
	 const int32_t initialRegion = simControl.shock->getShockRegion(sim->t,particleList[p].state);
	 PARTICLE initialState = particleList[p];
	 propagateParticle(particleList[p],species.propagationDirection,sim->t,sim->dt);

	 // Calculate distance to shock at trial position. If d_shock changes 
	 // sign (d_shock2_initial*d_shock < 0.0) particle crossed the shock and 
	 // we need to apply boundary conditions:
	 int32_t trialRegion = simControl.shock->getShockRegion(sim->t+sim->dt,particleList[p].state);
	 if (trialRegion != initialRegion) {
	    // Revert back to initial state:
	    particleList[p] = initialState;

	    // Find time of shock crossing. This will propagate particle 
	    // to t=t_crossing and set dt = dt - dt_taken, where dt_taken 
	    // is the time step taken by the binary search.
	    // Particle initialState is overwritten to contain incident state, particleList[p] 
	    // will end up on the opposite side of the shock (it becomes transmitted particle):
	    binarySearchShockCrossing(species.propagationDirection,initialState,particleList[p]);
	    
	    // Apply boundary conditions at shock. Particle "reflected" is a new 
	    // wave produced during shock crossing, it needs to be added to correct 
	    // particleList and propagated:
	    PARTICLE reflected;
	    applyShockBoundaryConditions(t,initialState,reflected,particleList[p]);

	    // Propagate reflected and transmitted wave packets to time t+dt:
	    propagateParticle(particleList[p],   species.propagationDirection,t,dt);

	    // If reflected particle cannot be created, don't bother propagating it:
	    propagateParticle(reflected      ,-1*species.propagationDirection,t,dt);
	    
	    // Calculate particle's i,j,k offset indices in block-based mesh:
	    const int32_t I_off = static_cast<uint32_t>(1 + (reflected.state[0]-i_block*block::WIDTH_X) / block::WIDTH_X) - 1;
	    const int32_t J_off = static_cast<uint32_t>(1 + (reflected.state[1]-j_block*block::WIDTH_Y) / block::WIDTH_Y) - 1;
	    const int32_t K_off = static_cast<uint32_t>(1 + (reflected.state[2]-k_block*block::WIDTH_Z) / block::WIDTH_Z) - 1;

	    // Insert reflected particle directly to target block. This ensures
	    // that ParticleList class will correctly move reflected particles
	    // between processes:
	    const pargrid::NeighbourID nbrOffset = simClasses->pargrid.calcNeighbourTypeID(I_off,J_off,K_off);
	    const pargrid::CellID newCell = nbrIDs[nbrOffset];

	    if (reflected.state[0] > 1.0) {
	       if (newCell != pargrid::INVALID_CELLID) {
		  injBuffer->insert(newCell,reflected);
	       }
	    }
	 }
      }

      return success;
   }

   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   bool LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								      const std::string& regionName,const ParticleListBase* plist) {
      initialized = ParticlePropagatorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());

      if (accelerator.initialize(sim,simClasses,cr) == false) {
	 simClasses.logger << "(SEP LAGR RK2) ERROR: Accelerator class failed to initialize" << std::endl << write;
	 initialized = false;
      }

      if (simControl.shock == NULL) {
	 simClasses.logger << "(SEP LAGR RK2) ERROR: Shock class in NULL" << std::endl << write;
	 initialized = false;
      } else {
	 maxShockDistance = simControl.shock->getAllowedDistance();
      }

      return initialized;
   }
   
   template<class SPECIES,class PARTICLE,class ACCELERATOR> inline
   void LagrangianPropagShockRK2<SPECIES,PARTICLE,ACCELERATOR>::propagateParticle(PARTICLE& particle,int propagationDirection,const Real& t,const Real& dt) {
      // Create data structures for containing particle state at intermediate time step:
      const int32_t VARIABLES = sep::lagr::DATASIZE;
      const int32_t NWL = simControl.N_wavelengthMeshCells;

      Real sign = 1.0;
      if (particle.state[lagr::LAMBDA] < simControl.N_wavelengthMeshCells/2) sign = -1;
      
      // Calculate correct particle position to temporary array defined below.
      // This is needed because shock drift acceleration causes propagator to 
      // be called twice for some particles, and cell indices will be incorrect 
      // for some of these particles due to periodic boundaries:
      Real correctState[VARIABLES];
      const int32_t I_fake1 = static_cast<int32_t>(particle.state[lagr::XCRD]);
      const int32_t J_fake1 = static_cast<int32_t>(particle.state[lagr::YCRD]);
      const int32_t K_fake1 = static_cast<int32_t>(particle.state[lagr::ZCRD]);
      const int32_t I_offset1 = I_fake1 - i_block*block::WIDTH_X + 1;
      const int32_t J_offset1 = J_fake1 - j_block*block::WIDTH_Y + 1;
      const int32_t K_offset1 = K_fake1 - k_block*block::WIDTH_Z + 1;
      correctState[0] = i_index[I_offset1] + (particle.state[0]-I_fake1);
      correctState[1] = j_index[J_offset1] + (particle.state[1]-J_fake1);
      correctState[2] = k_index[K_offset1] + (particle.state[2]-K_fake1);

      // Evaluate acceleration at time t:      
      Real y[VARIABLES];
      accelerator.get(blockLID,t,propagationDirection,correctState,y);
      const int region = simControl.shock->getShockRegion(t,correctState);
      
      // TEST
      if (y[0]*propagationDirection < 0) y[lagr::LAMBDA] *= -1;
      // END TEST

      // Propagate to intermediate time step:
      Real x[VARIABLES];
      x[0]            = correctState[0] + y[0]*0.5*dt / sim->dx_cell[0];
      x[1]            = correctState[1] + y[1]*0.5*dt / sim->dy_cell[0];
      x[2]            = correctState[2] + y[2]*0.5*dt / sim->dz_cell[0];
      if (fabs(particle.state[lagr::LAMBDA] - NWL/2) < 1.0) {
	 x[lagr::LAMBDA] = particle.state[lagr::LAMBDA] + (particle.state[lagr::LAMBDA] - NWL/2)*y[lagr::LAMBDA]*0.5*dt;
      } else {
	 x[lagr::LAMBDA] = particle.state[lagr::LAMBDA] + sign*log(1+y[lagr::LAMBDA]*0.5*dt)/simControl.logicalWavelengthCellSize;
      }

      if (simControl.shock->getShockRegion(t,x) != region) {
	 particle.state[0] += y[0]*dt / sim->dx_cell[0];
	 particle.state[1] += y[1]*dt / sim->dy_cell[0];
	 particle.state[2] += y[2]*dt / sim->dz_cell[0];
	 particle.state[lagr::LAMBDA] += sign*log(1+y[lagr::LAMBDA]*dt)/simControl.logicalWavelengthCellSize;
	 return;
      }

      // Evaluate acceleration at intermediate time step:
      accelerator.get(blockLID,t+0.5*dt,propagationDirection,x,y);

      // TEST
      if (y[0]*propagationDirection < 0) y[lagr::LAMBDA] *= -1;
      // END TEST

      // Propagate to final time step:
      particle.state[0] += y[0]*dt / sim->dx_cell[0];
      particle.state[1] += y[1]*dt / sim->dy_cell[0];
      particle.state[2] += y[2]*dt / sim->dz_cell[0];
      if (fabs(x[lagr::LAMBDA] - NWL/2) < 1.0) {
	 particle.state[lagr::LAMBDA] += (x[lagr::LAMBDA] - NWL/2)*y[lagr::LAMBDA]*dt;
      } else {
	 particle.state[lagr::LAMBDA] += sign*log(1+y[lagr::LAMBDA]*dt)/simControl.logicalWavelengthCellSize;
      }
   }

} // namespace sep
   
#endif
