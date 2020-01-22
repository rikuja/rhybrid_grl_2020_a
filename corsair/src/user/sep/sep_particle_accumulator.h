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

#ifndef SEP_PARTICLE_ACCUMULATOR_H
#define SEP_PARTICLE_ACCUMULATOR_H

#include <cstdlib>
#include <simulation.h>
#include <simulationclasses.h>
#include <base_class_particle_accumulator.h>

#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_simcontrol.h"
#include "sep_accumulation_stretched.h"
#include "sep_fields_container.h"

namespace sep {

   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE,int ORDER>
   class IonAccumulator: public ParticleAccumulatorBase {
    public:
      IonAccumulator();
      ~IonAccumulator();
   
      bool accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles);
      bool accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles);
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool addRemoteUpdates();
      bool clearAccumulationArrays();
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool sendUpdates();
      bool wait();
   
    private:
      #if PROFILE_LEVEL > 0
         int accumulationID;
         int arrayClearID;
         int remoteUpdatesID;
      #endif
   
      #if PROFILE_LEVEL > 1
         int accumID;
         int copyID;
         int tempArrayClear;
      #endif

      SPECIES species;                        /**< Particle species accumulated by this IonAccumulator.*/
      static int speciesCounter;              /**< All particle species accumulate their weights to same array, 
					       * having pargrid::DataID SimControl.particleWeightDataID. This 
					       * counter is used to make sure that IonAccumulator only clears 
					       * accumulation array once per timestep etc.*/
      
      void accumulate(pargrid::CellID block,Real* array,unsigned int N_particles,const PARTICLE* particleList);
      Real* getAccumulationArray();
   };

   template<class SPECIES,class PARTICLE,int ORDER> inline
   ParticleAccumulatorBase* IonAccumMaker() {return new IonAccumulator<SPECIES,PARTICLE,ORDER>();}
   
   // ***** START CLASS DEFINITION ***** //
   
   template<class SPECIES,class PARTICLE,int ORDER>
   int IonAccumulator<SPECIES,PARTICLE,ORDER>::speciesCounter = 0;
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   IonAccumulator<SPECIES,PARTICLE,ORDER>::IonAccumulator(): ParticleAccumulatorBase() {
      sim = NULL;
      simClasses = NULL;

      speciesCounter = 0;
      
      #if PROFILE_LEVEL > 0
         accumulationID  = -1;
	 arrayClearID    = -1;
	 remoteUpdatesID = -1;
      #endif
   
      #if PROFILE_LEVEL > 1
         accumID = -1;
         copyID = -1;
         tempArrayClear = -1;
      #endif
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   IonAccumulator<SPECIES,PARTICLE,ORDER>::~IonAccumulator() { }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   void IonAccumulator<SPECIES,PARTICLE,ORDER>::accumulate(pargrid::CellID blockLID,Real* RESTRICT array,
							   unsigned int N_particles,const PARTICLE* RESTRICT particleList) {
      uint32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      Real B[3];
      Real V_wave[3];
      Real V_wave_par,dV_wave;
      for (unsigned int p=0; p<N_particles; ++p) {
	 // Calculate particle's resonant wavelength in logical coordinates:
	 #warning sep_particle_accumulator.h Accumulator assumes that simulation time is sim->t
	 Real V_alfven;
	 (*simControl.fieldsGetState)(blockLID,sim->t,particleList[p].state,B,V_wave,dV_wave,V_alfven,simControl.alfvenSign);

	 // Calculate parallel plasma speed:
	 const Real B_mag = vectorMagnitude<3>(B);
	 
	 // Calculate resonant wavelength in selected (antiparallel, parallel) Alfven wave rest frame:
	 V_wave_par = dotProduct<3>(V_wave,B)/B_mag;
	 Real lambda_res = 2.0*M_PI/species.q_per_m/B_mag * (particleList[p].state[sep::particle::V_PAR] - V_wave_par);
	 lambda_res = (*simControl.getLogicalWavelength)(lambda_res);

	 // Do not accumulate particles that are outside wavelength mesh:
	 if ((lambda_res < 0) || (lambda_res >= simControl.N_wavelengthMeshCells)) continue;
	 
	 Real pos[4];
	 pos[0] = particleList[p].state[sep::particle::XCRD] - i_block*block::WIDTH_X + 1;
	 pos[1] = particleList[p].state[sep::particle::YCRD] - j_block*block::WIDTH_Y + 1;
	 pos[2] = particleList[p].state[sep::particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	 pos[3] = lambda_res + 1;

	 // TEST
	 Real shapeFactors[12];
	 int32_t indices[4];
	 // END TEST
	 const Real weight = particleList[p].state[sep::particle::WEIGHT];
	 switch (ORDER) {
	  case 0:
	    accumScalarCentroidLogisticNGP_4D(array,pos,weight,simControl.N_wavelengthMeshCells+2);
	    break;
	  case 1:
	    accumScalarCentroidLogisticCIC_4D(array,pos,weight,simControl.N_wavelengthMeshCells+2);
	    break;
	  case 2:
	    // TEST
	    getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	    accumScalarCentroidLogicalTSC_4D(array,indices,shapeFactors,weight,simControl.N_wavelengthMeshCells+2);
	    // END TEST
	    //accumScalarCentroidLogisticTSC_4D(array,pos,weight,simControl.N_wavelengthMeshCells+2);
	    break;
	  default:
	    simClasses->logger << "(PARTICLE ACCUM) ERROR: Unknown accumulation order of accuracy!" << std::endl << write;
	    exit(1);
	    break;
	 }
	 
      }
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** particleLists = wrapper.data();
      
      // Get array where particle weights are accumulated to:
      Real* phaseSpaceWeight = getAccumulationArray();
      if (phaseSpaceWeight == NULL) return false;

      #if PROFILE_LEVEL > 0
         profile::start("Accumulation (total)",accumulationID);
      #endif

      bool success = true;
      Real t_propag = 0.0;
      const uint32_t SIZE = block::SIZE_PLUS_ONE_LAYER*(simControl.N_wavelengthMeshCells+2);
      Real* array = new Real[SIZE];
      const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(sim->inverseStencilID);
      for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	 // Measure block accumulation time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	 
	 const pargrid::CellID blockLID = boundaryBlocks[b];
	 
	 // Clear temp array:
	 for (unsigned int i=0; i<SIZE; ++i) array[i] = 0.0;
	 
	 // Accumulate all particles on this block:
	 accumulate(blockLID,array,N_particles[blockLID],particleLists[blockLID]);
	 
	 // Copy values from temp array to phaseSpaceWeight:
	 addValues4D(simClasses,blockLID,array,phaseSpaceWeight,simControl.N_wavelengthMeshCells);
	 
	 // Store block accumulation time:
	 if (sim->countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
      delete [] array; array = NULL;
      
      #if PROFILE_LEVEL > 0
         profile::stop(); // Accumulation (total)
      #endif
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
	         
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** particleLists = wrapper.data();
	
      // Get array where particle weights are accumulated to:
      Real* phaseSpaceWeight = getAccumulationArray();
      if (phaseSpaceWeight == NULL) return false;
      #if PROFILE_LEVEL > 0
         profile::start("Accumulation (total)",accumulationID);
      #endif

      const uint32_t SIZE = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2)*(simControl.N_wavelengthMeshCells+2);

      Real t_propag = 0.0;
      Real* array = new Real[SIZE];
      const std::vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(sim->inverseStencilID);
      for (size_t b=0; b<innerBlocks.size(); ++b) {
	 // Measure block accumulation time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	 
	 const pargrid::CellID blockLID = innerBlocks[b];
      
	 // Clear temp array:
	 for (unsigned int i=0; i<SIZE; ++i) array[i] = 0.0;

	 // Accumulate all particles on this block:
	 accumulate(blockLID,array,N_particles[blockLID],particleLists[blockLID]);      

	 // Copy values from temp array to wavePowers:
	 addValues4D(simClasses,blockLID,array,phaseSpaceWeight,simControl.N_wavelengthMeshCells);
	 
	 // Store block accumulation time:
	 if (sim->countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
      delete [] array; array = NULL;
      
      #if PROFILE_LEVEL > 0
         profile::stop(); // accumulation
      #endif
      return true;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
      return true;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::addRemoteUpdates() {
      // Only the last particle species adds remote updates:
      if (speciesCounter != simControl.N_particleSpecies) return true;
      if (speciesCounter == simControl.N_particleSpecies) speciesCounter = 0;
      
      bool success = true;

      // Get accumulation array from ParGrid:
      Real* phaseSpaceWeight = getAccumulationArray();
      if (phaseSpaceWeight == NULL) return false;
      
      // Get remote updates:
      unsigned int* offsets = NULL;
      Real* remoteUpdateArray = NULL;
      if (simClasses->pargrid.getRemoteUpdates<Real>(sim->inverseStencilID,simControl.particleWeightDataID,
						     offsets,remoteUpdateArray) == false) {
	 simClasses->logger << "(SEP PARTICLE ACCUM) ERROR: Failed to get remote updates from ParGrid!" << std::endl << write;
	 success = false;
	 return success;
      }

      #if PROFILE_LEVEL > 0
         profile::start("add remote updates",remoteUpdatesID);
      #endif

      // Add updates received from remote neighbours. If weight factors were
      // accumulated their values are not divided by cell volumes:
      const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(sim->inverseStencilID);
      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      const uint32_t SIZE_BLOCK = block::SIZE*NWL;
      for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	 const pargrid::CellID blockLID = boundaryBlocks[b];
	 for (unsigned int offset=offsets[b]; offset<offsets[b+1]; ++offset) {
	    for (uint32_t l=0; l<SIZE_BLOCK; ++l) phaseSpaceWeight[blockLID*SIZE_BLOCK+l] += remoteUpdateArray[offset*SIZE_BLOCK+l];						       
	 }
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::clearAccumulationArrays() {
      // Only the first particle species clears accumulation array:
      ++speciesCounter;
      if (speciesCounter > 1) {
	 if (speciesCounter == simControl.N_particleSpecies) speciesCounter = 0;
	 return true;
      }
      
      // Clear species counter if there is only one particle species:
      if (speciesCounter == simControl.N_particleSpecies) {
	 speciesCounter = 0;
      }
      
      Real* phaseSpaceWeight = getAccumulationArray();
      if (phaseSpaceWeight == NULL) return false;
   
      // Clear wave powers array:
      const size_t SIZE = block::SIZE*simControl.N_wavelengthMeshCells;
      for (size_t i=0; i<simClasses->pargrid.getNumberOfAllCells()*SIZE; ++i) phaseSpaceWeight[i] = 0.0;
      return true;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::finalize() {
      return true;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   Real* IonAccumulator<SPECIES,PARTICLE,ORDER>::getAccumulationArray() {
      Real* phaseSpaceWeight = simClasses->pargrid.getUserDataStatic<Real>(simControl.particleWeightDataID);
      if (phaseSpaceWeight == NULL) {
	 simClasses->logger << "(SEP PARTICLE ACCUM) ERROR: Failed to get accumulation array!" << std::endl << write;
      }
      return phaseSpaceWeight;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
							   const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP PARTICLE ACCUM) Starting init" << std::endl << write;
      
      // Init base class:
      bool success = ParticleAccumulatorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::sendUpdates() {
      // Only last particle species sends updates to neighbor processes:
      ++speciesCounter;
      if (speciesCounter != simControl.N_particleSpecies) return true;
      speciesCounter = 0;
      
      bool success = true;
      if (simClasses->pargrid.startNeighbourExchange(sim->inverseStencilID,simControl.particleWeightDataID) == false) {
	 simClasses->logger << "(SEP PARTICLE ACCUM) ERROR: Failed to start data exchange with neighbours!" << std::endl << write;
	 success = false;
      }
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool IonAccumulator<SPECIES,PARTICLE,ORDER>::wait() {
      // Only last particle species waits for MPI:
      ++speciesCounter;
      if (speciesCounter != simControl.N_particleSpecies) return true;
      
      bool success = true;
            
      if (simClasses->pargrid.wait(sim->inverseStencilID,simControl.particleWeightDataID,"Particle Accumulator wait") == false) {
         simClasses->logger << "(SEP PARTICLE ACCUM) ERROR: ParGrid::wait returned false!" << std::endl << write;
         success = false;
      }
      
      return success;
   }

} // namespace sep
   
#endif
