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

#include <cstdlib>
#include <iostream>
#include <sstream>

#include <map>
#include <set>
#include <pargrid_buffers.h>
#include <pargrid_copy_protocol_mpi.h>

#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_accumulation_base.h"
#include "sep_particle_definition.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   
   OperatorAccumulationBase::OperatorAccumulationBase(uint32_t vectorSize,uint32_t order): DataOperator(),order(order),vectorSize(vectorSize) {
      initialized = false;
      accumArray = NULL;
      
      #if PROFILE_LEVEL > 0
         accumulationID = -1;
         bufferCopyID = -1;
         densityOperatorID = -1;
         mpiOverheadID = -1;
         mpiWaitID = -1;
         mpiWaitSendsID = -1;
      #endif
   }
   
   OperatorAccumulationBase::~OperatorAccumulationBase() {
      finalize();
   }

   bool OperatorAccumulationBase::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }
   
   uint32_t OperatorAccumulationBase::getOrder() const {return order;}

   uint32_t OperatorAccumulationBase::getVectorSize() const {return vectorSize;}
   
   bool OperatorAccumulationBase::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP EM FIELDS) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }

      return initialized;
   }
   
   bool OperatorAccumulationBase::postProcessData(Real* array) {
      return true;
   }
   
   bool OperatorAccumulationBase::preProcessData(Real* array) {
      return true;
   }

   #if PROFILE_LEVEL > 0
   void OperatorAccumulationBase::setProfileName(const std::string& name) {
      profileName = name;
   }
   #endif

   bool OperatorAccumulationBase::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         profile::start(profileName,densityOperatorID);
      #endif

      // Get number of arrays written by this Operator:
      N_arrays = getNumberOfArrays();
      
      // Write each array to output file:
      for (uint32_t i=0; i<N_arrays; ++i) {
	 currentArray = i;
	 if (writeData(i,spatMeshName,particleLists) == false) success = false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

   bool OperatorAccumulationBase::writeData(uint32_t arrayIndex,const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      bool success = setAccumulatedArray(arrayIndex);
      if (success == false) return success;

      // Get local cell IDs:
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();

      // Allocate global accumulation array in ParGrid:
      const pargrid::DataID accumDataID = simClasses->pargrid.addUserData<Real>("densityOperator",block::SIZE*vectorSize,false);
      if (accumDataID == pargrid::INVALID_DATAID) {
	 simClasses->logger << "(SEP OP DENSITY) ERROR: Failed to create accumulation array" << endl << write;
	 return false;
      }
      accumArray = simClasses->pargrid.getUserDataStatic<Real>(accumDataID);
      for (size_t i=0; i<simClasses->pargrid.getNumberOfAllCells()*block::SIZE*vectorSize; ++i) accumArray[i] = 0.0;

      if (preProcessData(NULL) == false) return false;
      
      // Add data transfer TBD
      #if PROFILE_LEVEL > 0
         profile::start("MPI overhead",mpiOverheadID);
      #endif
      typedef pargrid::Buffer<Real> BUFFER;
      vector<BUFFER> inBuffers;
      vector<BUFFER> outBuffers;
      std::vector<pargrid::CopyProtocolMPI<BUFFER> > inBufferCopy;
      std::vector<pargrid::CopyProtocolMPI<BUFFER> > outBufferCopy;
      std::vector<std::vector<pargrid::CellID> > recvLIDs;
      std::vector<std::vector<pargrid::CellID> > sendLIDs;
      
      const set<pargrid::MPI_processID>& neighbourProcesses = simClasses->pargrid.getNeighbourProcesses();
      inBuffers.resize(neighbourProcesses.size());
      outBuffers.resize(neighbourProcesses.size());
      inBufferCopy.resize(neighbourProcesses.size());
      outBufferCopy.resize(neighbourProcesses.size());
      recvLIDs.resize(neighbourProcesses.size());
      sendLIDs.resize(neighbourProcesses.size());
      
      size_t counter = 0;
      for (set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
	 const std::size_t N_recvs = simClasses->pargrid.getNumberOfReceives(sim->inverseStencilID,*it);
	 inBuffers[counter].resize(N_recvs);
	 if (inBufferCopy[counter].set(&(inBuffers[counter]),false,sim->comm,*it) == false) success = false;
	 
	 const std::size_t N_sends = simClasses->pargrid.getNumberOfSends(sim->inverseStencilID,*it);
	 outBuffers[counter].resize(N_sends);
	 if (outBufferCopy[counter].set(&(outBuffers[counter]),true,sim->comm,*it) == false) success = false;
      
	 const map<pargrid::MPI_processID,set<pargrid::CellID> >::const_iterator recvs
	   = simClasses->pargrid.getReceives(sim->inverseStencilID).find(*it);
	 recvLIDs[counter].resize(recvs->second.size());
	 size_t c=0;
	 for (set<pargrid::CellID>::const_iterator cellGID=recvs->second.begin(); cellGID!=recvs->second.end(); ++cellGID) {
	    const pargrid::CellID cellLID = simClasses->pargrid.getLocalID(*cellGID);
	    recvLIDs[counter][c] = cellLID;
	    ++c;
	 }
	 
	 const map<pargrid::MPI_processID,set<pargrid::CellID> >::const_iterator sends
	   = simClasses->pargrid.getSends(sim->inverseStencilID).find(*it);
	 sendLIDs[counter].resize(sends->second.size());
	 c=0;
	 for (set<pargrid::CellID>::const_iterator cellGID=sends->second.begin(); cellGID!=sends->second.end(); ++cellGID) {
	    const pargrid::CellID cellLID = simClasses->pargrid.getLocalID(*cellGID);
	    sendLIDs[counter][c] = cellLID;
	    ++c;
	 }
	 
	 ++counter;
      }
      
      // Post receives for data:
      MPITypes::rank bufferCounter = 0;
      for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
	 if (inBufferCopy[bufferCounter].start() == false) success = false;
	 ++bufferCounter;
      }
      #if PROFILE_LEVEL > 0
         profile::stop();
         profile::start("Accumulation",accumulationID);
      #endif
      
      // Accumulate in boundary blocks
      const vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(sim->inverseStencilID);
      for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	 accumulateBlock(boundaryBlocks[b],accumArray,particleLists);
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
         profile::start("Buffer Copy",bufferCopyID);
      #endif
      
      // Iterate over all neighbor processes who to send data:
      MPITypes::rank processCounter = 0;
      for (map<pargrid::MPI_processID,set<pargrid::CellID> >::const_iterator
	   it  = simClasses->pargrid.getSends(sim->inverseStencilID).begin();
	   it != simClasses->pargrid.getSends(sim->inverseStencilID).end();
	   ++it) {
	 
	 // Iterate over all blocks sent to neighbour process nbrHostID:
	 size_t blockCounter = 0;
	 vector<Real>& buffer = outBuffers[processCounter].getBuffer();
	 buffer.clear();
	 for (set<pargrid::CellID>::const_iterator block=it->second.begin(); block!=it->second.end(); ++block) {
	    const pargrid::CellID blockLID = sendLIDs[processCounter][blockCounter];

	    outBuffers[processCounter].getBlockSizes()[blockCounter] = block::SIZE*vectorSize;
	    ++blockCounter;

	    // Copy values from remote cell to output buffer:
	    for (uint32_t i=0; i<block::SIZE*vectorSize; ++i) {
	       buffer.push_back( accumArray[blockLID*block::SIZE*vectorSize+i] );
	    }
	 }

	 // Post send for particles sent to process nbrHostID:
	 if (outBufferCopy[processCounter].start() == false) success = false;
	 
	 ++processCounter;
      }
      #if PROFILE_LEVEL > 0
         profile::stop();                               // Buffer Copy
         profile::start("Accumulation",accumulationID);
      #endif
      
      // Accumulate in inner blocks
      const vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(sim->inverseStencilID);
      for (size_t b=0; b<innerBlocks.size(); ++b) {
	 accumulateBlock(innerBlocks[b],accumArray,particleLists);
      }
      
      // Wait for remote updates:
      #if PROFILE_LEVEL > 0
         profile::stop();                              // Accumulation
         profile::start("MPI Wait (recvs)",mpiWaitID);
      #endif
      bufferCounter = 0;
      for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
	 inBufferCopy[bufferCounter].wait();
	 ++bufferCounter;
      }
      #if PROFILE_LEVEL > 0
         profile::stop();                               // MPI Wait (recsv)
         profile::start("Buffer Copy",bufferCopyID);
      #endif

      // Add remote updates to local values:
      bufferCounter = 0;
      for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
	 size_t offset = 0;
	 const std::vector<Real>& buffer = inBuffers[bufferCounter].getBuffer();
	 for (size_t b=0; b<inBuffers[bufferCounter].getNumberOfBlocks(); ++b) {
	    const pargrid::CellID blockLID = recvLIDs[bufferCounter][b];
	    
	    for (uint32_t i=0; i<block::SIZE*vectorSize; ++i) {
	       accumArray[blockLID*block::SIZE*vectorSize+i] += buffer[offset];
	       ++offset;
	    }
	 }
	 ++bufferCounter;
      }

      #if PROFILE_LEVEL > 0
         profile::stop();                               // Buffer Copy
      #endif

      // Post process accumulated data:
      if (postProcessData(accumArray) == false) return false;
      
      // Write arrays to output file:
      std::map<std::string,std::string> attribs;
      attribs["name"] = spatMeshName+'/'+getOutputName(arrayIndex);
      attribs["mesh"] = spatMeshName;
      attribs["unit"] = getOutputUnits(arrayIndex);
      attribs["centering"] = "zone";
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
       case sep::CARTESIAN:
	 attribs["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
       case sep::CYLINDRICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 break;
       case sep::SPHERICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	 break;
       default:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
      }      
      
      const uint64_t arraySize = N_localBlocks*block::SIZE;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,accumArray) == false) success = false;

      // Deallocate memory and exit:
      if (simClasses->pargrid.removeUserData(accumDataID) == false) success = false;

      // Wait for sends to complete:
      #if PROFILE_LEVEL > 0
         profile::start("MPI Wait (sends)",mpiWaitSendsID);
      #endif
      
      bufferCounter = 0;
      for (std::set<pargrid::MPI_processID>::const_iterator it=neighbourProcesses.begin(); it!=neighbourProcesses.end(); ++it) {
	 outBufferCopy[bufferCounter].wait();
	 ++bufferCounter;
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif

      return success;
   }
      
   
} // namespace sep
