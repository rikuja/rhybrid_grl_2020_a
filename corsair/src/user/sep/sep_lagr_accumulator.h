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

#ifndef SEP_LAGR_ACCUMULATOR_H
#define SEP_LAGR_ACCUMULATOR_H

#include <cstdlib>
#include <simulation.h>
#include <simulationclasses.h>
#include <base_class_particle_accumulator.h>

#include "sep_lagr_definition.h"
#include "sep_lagr_species.h"
#include "sep_simcontrol.h"
#include "sep_accumulation_stretched.h"
#include "sep_wavelength_mesh_builder.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE,int ORDER>
   class LagrangianAccumulator: public ParticleAccumulatorBase {
    public:
      LagrangianAccumulator();
      ~LagrangianAccumulator();
   
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
         int mpiWaitID;
         int remoteUpdatesID;
      #endif

      SPECIES species;

      bool accumulate(const std::vector<pargrid::CellID>& blockLIDs,const unsigned int* N_particles,PARTICLE** particleList);
      Real* getAccumulationArray();
      Real* getAccumulationArray(pargrid::DataID& dataID);
   };

   template<class SPECIES,class PARTICLE,int ORDER> inline
   ParticleAccumulatorBase* LagrAccumMaker() {return new LagrangianAccumulator<SPECIES,PARTICLE,ORDER>();}
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::LagrangianAccumulator(): ParticleAccumulatorBase() {
      #if PROFILE_LEVEL > 0
         accumulationID  = -1;
	 arrayClearID    = -1;
	 mpiWaitID       = -1;
	 remoteUpdatesID = -1;
      #endif
   
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::~LagrangianAccumulator() { }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::accumulate(const std::vector<pargrid::CellID>& blockLIDs,
								  const unsigned int* N_particles,PARTICLE** particles) {
      // Get correct wave energy array:
      Real* RESTRICT waveEnergyGlobal = getAccumulationArray();
      if (waveEnergyGlobal == NULL) return false;

      // Get array containing min/max valid wavelength bins:
      wmesh::validDatatype* validBins = NULL;
      pargrid::DataID dataID = pargrid::INVALID_DATAID;
      if (species.propagationDirection < 0) {
	 dataID = simControl.maxAntiparValidWavelengthBinsDataID;
	 validBins = simClasses->pargrid.getUserDataStatic<wmesh::validDatatype>(dataID);
      } else {
	 dataID = simControl.maxParValidWavelengthBinsDataID;
	 validBins = simClasses->pargrid.getUserDataStatic<wmesh::validDatatype>(dataID);
      }
      wmesh::validDatatype tmpValidBins[block::SIZE*2];

      bool success = true;
      const uint32_t SIZE = block::SIZE_PLUS_ONE_LAYER*(simControl.N_wavelengthMeshCells+2);
      
      Real t_propag = 0.0;
      Real* array = new Real[SIZE];
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 // Measure block accumulation time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	 
	 // Clear temp arrays:
	 for (unsigned int i=0; i<SIZE; ++i) array[i] = 0.0;
	 for (int32_t i=0; i<block::SIZE; ++i) {
	    tmpValidBins[2*i+0] = simControl.N_wavelengthMeshCells;
	    tmpValidBins[2*i+1] = 0;
	 }
	  
	 // Calculate block indices:
	 const pargrid::CellID blockLID = blockLIDs[b];
	 uint32_t blockIndices[3];
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(*sim,blockGID,blockIndices[0],blockIndices[1],blockIndices[2]);

	 // If accumulation is upwinded/downwinded near shock, classify cells:
	 bool shockedBlock = false;
	 Real pos[4];
	 int32_t shockRegions[block::WIDTH_X+2];
	 if (simControl.useShockUpwinding == true) {
	    shockedBlock = classifyShockedCells(sim->t,shockRegions,blockIndices);
	 }

	 // Accumulate all particles on this block:
	 for (unsigned int p=0; p<N_particles[blockLID]; ++p) {
	    // Skip particles that are outside wavelength mesh in lambda direction:
	    if (particles[blockLID][p].state[lagr::LAMBDA] < 0.0) continue;
	    if (particles[blockLID][p].state[lagr::LAMBDA] > simControl.N_wavelengthMeshCells) continue;

	    // Calculate wave packet's offset in block:
	    pos[0] = particles[blockLID][p].state[lagr::XCRD]   - blockIndices[0]*block::WIDTH_X + 1;
	    pos[1] = particles[blockLID][p].state[lagr::YCRD]   - blockIndices[1]*block::WIDTH_Y + 1;
	    pos[2] = particles[blockLID][p].state[lagr::ZCRD]   - blockIndices[2]*block::WIDTH_Z + 1;
	    pos[3] = particles[blockLID][p].state[lagr::LAMBDA] + 1;

	    const Real weight = particles[blockLID][p].state[lagr::ENERGY];
	    
	    int32_t index;
	    int32_t indices[4];
	    Real shapeFactors[12];
	    
	    #ifndef NDEBUG
	       bool accumSuccess = true;
	    #endif
	    switch (ORDER) {
	     case 0:
	       getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	       
	       // Calculate upwinded/downwinded shape factors (if necessary):
	       if (shockedBlock == true) {
		  calculateUpwindedWaveShapeFactors(shockRegions,indices,shapeFactors,
						    simControl.shock->getShockRegion(sim->t,particles[blockLID][p].state));
	       }
	       
	       accumScalarCentroidLogicalNGP_4D(array,indices,shapeFactors,weight,simControl.N_wavelengthMeshCells+2);
	       if (weight > 1e-20) {
		  index = block::index(indices[0]-1,indices[1]-1,indices[2]-1);
		  if (indices[3]-1 < tmpValidBins[index*2+0]) tmpValidBins[index*2+0] = indices[3]-1;
		  if (indices[3]-1 > tmpValidBins[index*2+1]) tmpValidBins[index*2+1] = indices[3]-1;
	       }
	       break;
	     case 1:
	       getShapeFactorsCIC_4D(pos,indices,shapeFactors);

	       // Calculate upwinded/downwinded shape factors (if necessary):
	       if (shockedBlock == true) {
		  calculateUpwindedWaveShapeFactors(shockRegions,indices,shapeFactors,
						    simControl.shock->getShockRegion(sim->t,particles[blockLID][p].state));
	       }

	       accumScalarCentroidLogicalCIC_4D(array,indices,shapeFactors,weight,simControl.N_wavelengthMeshCells+2);

	       if (weight > 1e-20) {
		  index = block::index(static_cast<int32_t>(pos[0])-1,static_cast<int32_t>(pos[1])-1,static_cast<int32_t>(pos[2])-1);
		  const int32_t L_index = static_cast<int32_t>(pos[3])-1;
		  if (L_index < tmpValidBins[index*2+0]) tmpValidBins[index*2+0] = L_index;
		  if (L_index > tmpValidBins[index*2+1]) tmpValidBins[index*2+1] = L_index;
	       }
	       break;
	     case 2:
	       getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	       
	       // Calculate upwinded/downwinded shape factors (if necessary):
	       if (shockedBlock == true) {
		  calculateUpwindedWaveShapeFactors(shockRegions,indices,shapeFactors,
						    simControl.shock->getShockRegion(sim->t,particles[blockLID][p].state));
	       }
	       
	       accumScalarCentroidLogicalTSC_4D(array,indices,shapeFactors,weight,simControl.N_wavelengthMeshCells+2);

	       if (weight > 1e-20) {
		  index = block::index(indices[0]-1,indices[1]-1,indices[2]-1);
		  if (indices[3]-1 < tmpValidBins[index*2+0]) tmpValidBins[index*2+0] = indices[3]-1;
		  if (indices[3]-1 > tmpValidBins[index*2+1]) tmpValidBins[index*2+1] = indices[3]-1;
	       }
	       break;
	     default:
	       simClasses->logger << "(LAGR ACCUM) ERROR: Unknown accumulation order of accuracy!" << std::endl << write;
	       exit(1);
	       break;
	    }
	    
	    #ifndef NDEBUG
	    // Debugging only -- check that accumulation succeeded:
	    if (accumSuccess != true) {
	       const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	       uint32_t i_block,j_block,k_block;
	       block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	       
	       std::cerr << "(SEP LAGR ACCUM) ERROR: Wave energy accumulation failed" << std::endl;
	       std::cerr << "\t STEP: " << sim->timestep << std::endl;
	       std::cerr << "\t block LID#" << blockLID << " GID#" << blockGID << std::endl;
	       std::cerr << "\t block indices: " << i_block << ' ' << j_block << ' ' << k_block << std::endl;
	       std::cerr << "\t wave pos: ";
	       for (int i=0; i<4; ++i) std::cerr << particles[blockLID][p].state[i] << '\t';
	       std::cerr << std::endl;
	       exit(1);
	    }
	    #endif
	 }
	  
	 // Copy values from temp array to wavePowers:
	 addValues4D(simClasses,blockLID,array,waveEnergyGlobal,simControl.N_wavelengthMeshCells);	 
	 for (int i=0; i<block::SIZE*2; ++i) {
	    validBins[blockLID*block::SIZE*2+i] = tmpValidBins[i];
	 }
	 // Store block accumulation time:
	 if (sim->countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
      delete [] array; array = NULL;
      
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
      #if PROFILE_LEVEL > 0
         profile::start("Accumulation (total)",accumulationID);
      #endif

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** particleLists = wrapper.data();
      const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(sim->inverseStencilID);
      bool success = accumulate(boundaryBlocks,N_particles,particleLists);

      #if PROFILE_LEVEL > 0
         profile::stop(); // Accumulation (total)
      #endif
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
      #if PROFILE_LEVEL > 0
         profile::start("Accumulation (total)",accumulationID);
      #endif

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      PARTICLE** particleLists = wrapper.data();
      const std::vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(sim->inverseStencilID);
      bool success = accumulate(innerBlocks,N_particles,particleLists);

      #if PROFILE_LEVEL > 0
         profile::stop(); // accumulation (total)
      #endif
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
      return true;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::addRemoteUpdates() {
      bool success = true;
      
      // Get correct wave energy array:
      pargrid::DataID waveEnergyDataID;
      Real* waveEnergy = getAccumulationArray(waveEnergyDataID);
      if (waveEnergy == NULL) return false;

      // Get remote updates:
       unsigned int* offsets = NULL;
      Real* remoteUpdateArray = NULL;
      if (simClasses->pargrid.getRemoteUpdates<Real>(sim->inverseStencilID,waveEnergyDataID,offsets,remoteUpdateArray) == false) {
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: Failed to get remote updates from ParGrid!" << std::endl << write;
	 success = false;
	 return success;
      }

      // Add updates received from remote neighbours:
      const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(sim->inverseStencilID);
      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      const uint32_t SIZE_BLOCK = block::SIZE*NWL;
      for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	 const pargrid::CellID blockLID = boundaryBlocks[b];
	 for (unsigned int offset=offsets[b]; offset<offsets[b+1]; ++offset) {
	    for (uint32_t l=0; l<SIZE_BLOCK; ++l) waveEnergy[blockLID*SIZE_BLOCK+l] += remoteUpdateArray[offset*SIZE_BLOCK+l];
	 }
      }

      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::clearAccumulationArrays() {
      // Get correct wave energy array:
      Real* waveEnergy = getAccumulationArray();
      if (waveEnergy == NULL) return false;
   
      // Clear wave energy array:
      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      const uint32_t SIZE_BLOCK = block::SIZE*NWL;
      for (size_t b=0; b<simClasses->pargrid.getNumberOfAllCells()*SIZE_BLOCK; ++b) waveEnergy[b] = 0.0;
      return true;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::finalize() {
      return true;
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   Real* LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::getAccumulationArray() {
      pargrid::DataID dummy;
      return getAccumulationArray(dummy);
   }
   
   template<class SPECIES,class PARTICLE,int ORDER> inline
   Real* LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::getAccumulationArray(pargrid::DataID& dataID) {
      // Check that accumulator is called for Lagrangian species:
      if (species.getSpeciesType() != simControl.lagrangianSpeciesTypename) {
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: Accumulator called for non-Lagrangian species" << std::endl << write;
	 dataID = pargrid::INVALID_DATAID;
	 return NULL;
      }
      
      // Get correct wave energy array:
      dataID = pargrid::INVALID_DATAID;
      if (species.propagationDirection < 0) {
	 dataID = simControl.antiparAlfvenWaveEnergyDataID;
      } else if (species.propagationDirection > 0) {
	 dataID = simControl.parAlfvenWaveEnergyDataID;
      } else {
	 dataID = pargrid::INVALID_DATAID;
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: Could not get wave energy array" << std::endl << write;
	 return NULL;
      }
      Real* waveEnergyGlobal = simClasses->pargrid.getUserDataStatic<Real>(dataID);
      if (waveEnergyGlobal == NULL) {
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: Wave energy array is NULL" << std::endl << write;
	 dataID = pargrid::INVALID_DATAID;
	 return NULL;
      }
      
      return waveEnergyGlobal;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								  const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP LAGR ACCUM) Starting init" << std::endl << write;
      bool success = ParticleAccumulatorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::sendUpdates() {
      bool success = true;

      // Send updates from remote cells to neighbor processes:
      pargrid::DataID waveEnergyDataID = pargrid::INVALID_DATAID;
      getAccumulationArray(waveEnergyDataID);
      
      if (simClasses->pargrid.startNeighbourExchange(sim->inverseStencilID,waveEnergyDataID) == false) {
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: Failed to start data exchange with neighbours!" << std::endl << write;
	 success = false;
      }
      
      return success;
   }

   template<class SPECIES,class PARTICLE,int ORDER> inline
   bool LagrangianAccumulator<SPECIES,PARTICLE,ORDER>::wait() {
      bool success = true;

      pargrid::DataID waveEnergyDataID;
      getAccumulationArray(waveEnergyDataID);
      
      if (simClasses->pargrid.wait(sim->inverseStencilID,waveEnergyDataID,"Lagrangian accumulator wait") == false) {
	 simClasses->logger << "(SEP LAGR ACCUM) ERROR: ParGrid::wait returned false!" << std::endl << write;
	 success = false;
      }

      return success;
   }

} // namespace sep
   
#endif
