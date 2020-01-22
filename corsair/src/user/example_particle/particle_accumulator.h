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

#ifndef PARTICLE_ACCUMULATOR_H
#define PARTICLE_ACCUMULATOR_H

#include <cstdlib>
#include <simulation.h>
#include <simulationclasses.h>
#include <accumulators.h>
#include <base_class_particle_accumulator.h>

template<class SPECIES,class PARTICLE,int ORDER>
class Accumulator: public ParticleAccumulatorBase {
 public:
   Accumulator();
   ~Accumulator();
   
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
      int addRemoteUpdatesID;
      int particleAccumulation;
      int dataCopying;
      int mpiWaits;
      int arrayClearing;
      int tempArrayClearing;
      int overhead;
   #endif
   
   const SPECIES* species;
};

// Maker function that returns new Accumulator:
template<class SPECIES,class PARTICLE,int ORDER> inline
ParticleAccumulatorBase* AccumulatorMaker() {return new Accumulator<SPECIES,PARTICLE,ORDER>();}

template<class SPECIES,class PARTICLE,int ORDER> inline
Accumulator<SPECIES,PARTICLE,ORDER>::Accumulator(): ParticleAccumulatorBase() {
   #if PROFILE_LEVEL > 0
      arrayClearing = -1;
      tempArrayClearing = -1;
      addRemoteUpdatesID = -1;
      particleAccumulation = -1;
      dataCopying = -1;
      mpiWaits = -1;
      overhead = -1;
   #endif
}

template<class SPECIES,class PARTICLE,int ORDER> inline
Accumulator<SPECIES,PARTICLE,ORDER>::~Accumulator() {
   finalize();
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
   bool success = true;
   
   #if PROFILE_LEVEL > 0
      profile::start("overhead",overhead);
   #endif
   
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   PARTICLE** particleLists = wrapper.data();
   Real* densityArray = simClasses->pargrid.getUserDataStatic<Real>(SimControl::densityDataID);
   const int SIZE = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[SIZE];
   Real cellSizes[3];
   const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(SimControl::densityStencilID);
   
   #if PROFILE_LEVEL > 0
      profile::stop(); // overhead
   #endif
   
   Real t_propag = 0.0;
   for (size_t b=0; b<boundaryBlocks.size(); ++b) {
      // Measure accumulation time if we are preparing for mesh repartitioning:
      if (sim->countPropagTime == true) {
	 t_propag = MPI_Wtime();
      }
      
      #if PROFILE_LEVEL > 1
         profile::start("temp array clearing",tempArrayClearing);
      #endif
      const pargrid::CellID block = boundaryBlocks[b];
      for (int i=0; i<SIZE; ++i) array[i] = 0.0;
      getBlockCellSize(*simClasses,*sim,block,cellSizes);
      
      #if PROFILE_LEVEL > 1
         profile::stop(); // temp array clearing
         profile::start("accumulation",particleAccumulation);
      #endif
      
      switch (ORDER) {
       case 0:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarNGP_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       case 1:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarCIC_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       case 2:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarTSC_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       default:
	 std::cerr << "(PARTICLE ACCUMULATOR) ERROR: Unknown order of accuracy!" << std::endl;
	 exit(1);
	 break;
      }
	  
      #if PROFILE_LEVEL > 1
         profile::stop(); // accumulation
         profile::start("data copying",dataCopying);
      #endif
      
      // Add values to mesh:
      block::addValues3D(*simClasses,block,array,densityArray);
      
      #if PROFILE_LEVEL > 1
         profile::stop(); // data copying
      #endif
      
      if (sim->countPropagTime == true) {
	 simClasses->pargrid.getCellWeights()[block] += (MPI_Wtime() - t_propag);
      }
   }
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
   bool success = true;
   #if PROFILE_LEVEL > 0
      profile::start("overhead",overhead);
   #endif

   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   PARTICLE** particleLists = wrapper.data();
   Real* densityArray = simClasses->pargrid.getUserDataStatic<Real>(SimControl::densityDataID);
   Real cellSizes[3];
   const int SIZE = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[SIZE];
   const std::vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(SimControl::densityStencilID);

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   
   Real t_propag = 0.0;
   for (size_t b=0; b<innerBlocks.size(); ++b) {
      // Measure accumulation time if we are preparing for mesh repartitioning:
      if (sim->countPropagTime == true) {
	 t_propag = MPI_Wtime();
      }
      
      #if PROFILE_LEVEL > 1
         profile::start("temp array clearing",tempArrayClearing);
      #endif
      const pargrid::CellID block = innerBlocks[b];
      for (int i=0; i<SIZE; ++i) array[i] = 0.0;
      getBlockCellSize(*simClasses,*sim,block,cellSizes);

      #if PROFILE_LEVEL > 1
         profile::stop(); // temp array clearing
         profile::start("accumulation",particleAccumulation);
      #endif

      switch (ORDER) {
       case 0:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarNGP_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       case 1:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarCIC_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       case 2:
	 for (unsigned int p=0; p<N_particles[block]; ++p) {
	    const Real WEIGHT = particleLists[block][p].state[particle::WEIGHT];
	    accumulateScalarTSC_XYZ(cellSizes,array,WEIGHT,particleLists[block][p]);
	 }
	 break;
       default:
	 std::cerr << "(PARTICLE ACCUMULATOR) ERROR: Unknown order of accuracy!" << std::endl;
	 exit(1);
	 break;
      }

      #if PROFILE_LEVEL > 1
         profile::stop(); // accumulation
         profile::start("data copying",dataCopying);
      #endif
      
      // Add values to mesh:
      block::addValues3D(*simClasses,block,array,densityArray);

      #if PROFILE_LEVEL > 1
         profile::stop(); // data copying
      #endif
      
      if (sim->countPropagTime == true) {
	 simClasses->pargrid.getCellWeights()[block] += (MPI_Wtime() - t_propag);
      }
   }
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
   return true;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::addRemoteUpdates() {
   bool success = true;

   #if PROFILE_LEVEL > 0
      profile::start("overhead",overhead);
   #endif

   // Get density array from ParGrid:
   Real* densityArray = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(SimControl::densityDataID));
   if (densityArray == NULL) {
      simClasses->logger << "(ACCUMULATOR) ERROR: Failed to get densityArray from ParGrid" << std::endl << write;
      success = false;
   }
   
   // Get remote updates from ParGrid. "bufferPtr" is a pointer to an array containing 
   // all updates received from remote neighbours for all blocks. "offsets" is an array that 
   // tells which elements in "bufferPtr" belong to each boundary block:
   unsigned int* offsets = NULL;
   char* bufferPtr = NULL;
   if (simClasses->pargrid.getRemoteUpdates(SimControl::densityStencilID,SimControl::densityDataID,offsets,bufferPtr) == false) {
      simClasses->logger << "(ACCUMULATOR) ERROR: Failed to obtain remote updates from ParGrid!" << std::endl << write;
      success = false;
   }
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   
   if (success == true) {
      #if PROFILE_LEVEL > 0
         profile::start("add remote updates",addRemoteUpdatesID);
      #endif
      
      // Add updates received from remote process(es) to local cells:
      Real* buffer = reinterpret_cast<Real*>(bufferPtr);
      const std::vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(SimControl::densityStencilID);
      for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	 const pargrid::CellID blockLID = boundaryBlocks[b];
	 for (unsigned int i=offsets[b]; i<offsets[b+1]; ++i) {
	    for (int j=0; j<block::SIZE; ++j) {
	       densityArray[blockLID*block::SIZE+j] += buffer[i*block::SIZE+j];
	    }
	 }
      }
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
   }
   
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::clearAccumulationArrays() {
   bool success = true;

   #if PROFILE_LEVEL > 0
      profile::start("array clearing",arrayClearing);
   #endif
      
   Real* densityArray = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(SimControl::densityDataID));
   if (densityArray == NULL) return false;
   for (size_t i=0; i<simClasses->pargrid.getNumberOfAllCells()*block::SIZE; ++i) densityArray[i] = 0.0;
   
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::finalize() {
      bool success = true;
      sim = NULL;
      simClasses = NULL;
      return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
						     const std::string& regionName,const ParticleListBase* plist) {
   bool success = ParticleAccumulatorBase::initialize(sim,simClasses,cr,regionName,plist);
   species = reinterpret_cast<const SPECIES*>(plist->getSpecies());
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::sendUpdates() {
   bool success = true;
   
   // Send locally calculated updates on remote cells to neighbour process(es):
   if (simClasses->pargrid.startNeighbourExchange(SimControl::densityStencilID,SimControl::densityDataID) == false) {
      simClasses->logger << "(ACCUMULATOR) ERROR: Failed to send updates to neighbour process(es)!" << std::endl << write;
      success = false;
   }
   
   return success;
}

template<class SPECIES,class PARTICLE,int ORDER> inline
bool Accumulator<SPECIES,PARTICLE,ORDER>::wait() {
   bool success = true;
   #if PROFILE_LEVEL > 0
      profile::start("MPI waits",mpiWaits);
   #endif
   
   if (simClasses->pargrid.wait(SimControl::densityStencilID,SimControl::densityDataID) == false) {
      simClasses->logger << "(ACCUMULATOR) ERROR: MPI wait failed!" << std::endl << write;
      success = false;
   }
   
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   
   return success;
}

#endif
