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

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_operator_wavelength_mesh_density.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const string PREFIX = "WavelengthMeshDensity";
   
   WavelengthMeshDensityOP::WavelengthMeshDensityOP(): SpatialSliceOP() { 

   }
   
   WavelengthMeshDensityOP::~WavelengthMeshDensityOP() { }

   bool WavelengthMeshDensityOP::finalize() { 
      return true;
   }
   
   std::string WavelengthMeshDensityOP::getName() const { 
      return "WLM_DensityOP";
   }
   
   bool WavelengthMeshDensityOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) { 
      // Init base class:
      if (SpatialSliceOP::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP WLM VOLUME) Base class failed to init" << endl << write;
	 return false;
      }
      
      bool success = true;
      return success;
   }
   
   bool WavelengthMeshDensityOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) { 
      bool success = true;

      // Get array where phase-space densities are accumulated to:
      Real* phaseSpaceWeight = simClasses->pargrid.getUserDataStatic<Real>(simControl.particleWeightDataID);
      if (phaseSpaceWeight == NULL) {
	 simClasses->logger << "(SEP OP WLM DENSITY) ERROR: Could not find accumulation array" << endl << write;
	 return false;
      }
      
      // Phase-space weights are calculates on parallel Alfven rest frame:
      simControl.alfvenSign = +1;

      #if PROFILE_LEVEL > 0
         static int profTotalTime = -1;
         profile::start("WLM Density",profTotalTime);
      #endif

      recalculateCellVolumes(*sim,*simClasses);

      for (size_t p=0; p<particleLists.size(); ++p) {
	 // Skip wave packets:
	 if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 
	 // Clear accumulation array:
	 if (particleLists[p]->clearAccumulationArrays() == false) success = false;
	 
	 // Accumulate particle phase-space weights:
	 if (particleLists[p]->accumulateBoundaryCells() == false) success = false;
	 if (particleLists[p]->accumulateInnerCells() == false) success = false;	 

	 // Write phase-space density to all slices:
	 for (size_t slice=0; slice<getNumberOfSlices(); ++slice) {
	    stringstream ss;
	    ss << "WaveMesh" << slice;
	    const string meshName = ss.str();

	    // Get IDs of blocks in slice:
	    prepareSlice(slice,simControl.N_wavelengthMeshCells);
	    size_t counter = 0;
	    vector<pargrid::CellID> blockLIDs;
	    vector<pargrid::CellID> newBlockGIDs;
	    getAcceptedBlocks(slice,newBlockGIDs,blockLIDs);

	    int32_t sliceIndex = numeric_limits<int32_t>::max();
	    for (size_t b=0; b<blockLIDs.size(); ++b) {
	       int32_t i_block,j_block,k_block;
	       const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLIDs[b]];
	       block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	       switch (getSlicedCoordinate(slice)) {
		case 0:
		  for (int32_t i=0; i<block::WIDTH_X; ++i) {
		     if (sim->x_crds_node[i_block*block::WIDTH_X+i] > getSliceOrigin(slice)) continue;
		     if (getSliceOrigin(slice) > sim->x_crds_node[i_block*block::WIDTH_X+i+1]) continue;
		     sliceIndex = i;
		     break;
		  }
		  break;
		case 1:
		  for (int32_t j=0; j<block::WIDTH_Y; ++j) {
		     if (sim->y_crds_node[j_block*block::WIDTH_Y+j] > getSliceOrigin(slice)) continue;
		     if (getSliceOrigin(slice) > sim->y_crds_node[j_block*block::WIDTH_Y+j+1]) continue;
		     sliceIndex = j;
		     break;
		  }
		  break;
		case 2:
		  for (int32_t k=0; k<block::WIDTH_Z; ++k) {
		  if (sim->z_crds_node[k_block*block::WIDTH_Z+k] > getSliceOrigin(slice)) continue;
		     if (getSliceOrigin(slice) > sim->z_crds_node[k_block*block::WIDTH_Z+k+1]) continue;
		     sliceIndex = k;
		     break;
		  }
		  break;
		default:
		  continue;
		  break;
	       }
	       break;
	    }
	 
	    if (blockLIDs.size() > 0 && sliceIndex == numeric_limits<int32_t>::max()) {
	       simClasses->logger << "(SEP OP WLM VOLUME) ERROR: Could not find sliceIndex" << endl << write;
	       exit(1);
	    }

	    // Create a buffer for phase-space volume:
	    size_t arraySize;
	    switch (getSlicedCoordinate(slice)) {
	     case 0:
	       arraySize = block::WIDTH_Y*block::WIDTH_Z*simControl.N_wavelengthMeshCells;
	       break;
	     case 1:
	       arraySize = block::WIDTH_X*block::WIDTH_Z*simControl.N_wavelengthMeshCells;
	       break;
	     case 2:
	       arraySize = block::WIDTH_X*block::WIDTH_Y*simControl.N_wavelengthMeshCells;
	       break;
	     default:
	       continue;
	       break;
	    }
	    arraySize *= blockLIDs.size();
	    vector<Real> buffer;
	    buffer.resize(arraySize);

	    const Real* RESTRICT dl = simControl.wavelengthMeshCellSizes;

	    counter = 0;
	    for (size_t b=0; b<blockLIDs.size(); ++b) {
	       // Get block i,j,k indices:
	       const pargrid::CellID blockLID = blockLIDs[b];
	       const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	       uint32_t i_block,j_block,k_block;
	       block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

	       // Divide particle phase-space weights by phase-space volumes:
	       const uint32_t NWL = simControl.N_wavelengthMeshCells;
	       switch (getSlicedCoordinate(slice)) {
		case 0:
		  for (uint32_t l=0; l<NWL; ++l) for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) {
		     const uint64_t index1 = blockLID*block::SIZE+block::index(sliceIndex,j,k);
		     const uint64_t index2 = blockLID*block::SIZE*NWL + block::index(sliceIndex,j,k)*NWL + l;
		     buffer[counter] = phaseSpaceWeight[index2] / (simControl.cellVolumes[index1] * dl[l]);
		     ++counter;
		  }
		  break;
		case 1:
		  for (uint32_t l=0; l<NWL; ++l) for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) {
		     const uint64_t index1 = blockLID*block::SIZE+block::index(i,sliceIndex,k);
		     const uint64_t index2 = blockLID*block::SIZE*NWL + block::index(i,sliceIndex,k)*NWL + l;
		     buffer[counter] = phaseSpaceWeight[index2] / (simControl.cellVolumes[index1] * dl[l]);
		     ++counter;
		  }
		  break;
		case 2:
		  for (uint32_t l=0; l<NWL; ++l) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
		     const uint64_t index1 = blockLID*block::SIZE+block::index(i,j,sliceIndex);
		     const uint64_t index2 = blockLID*block::SIZE*NWL + block::index(i,j,sliceIndex)*NWL + l;
		     buffer[counter] = phaseSpaceWeight[index2] / (simControl.cellVolumes[index1] * dl[l]);
		     ++counter;
		  }
		  break;
		default:
		  exit(1);
		  break;
	       }
	    }
	 
	    // Write phase-space volume to output file:
	    map<string,string> xmlAttributes;
	    xmlAttributes["name"] = meshName + "/Density_" + particleLists[p]->getSpeciesName();
	    xmlAttributes["mesh"] = meshName;
	    xmlAttributes["centering"] = "zone";
	    switch (simControl.coordinateSystem) {
	     case sep::UNKNOWN:
	       xmlAttributes["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	       break;
	     case sep::CARTESIAN:
	       xmlAttributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	       break;
	     case sep::CYLINDRICAL:
	       xmlAttributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	       break;
	     case sep::SPHERICAL:
	       xmlAttributes["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	       break;
	     default:
	       xmlAttributes["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	       break;
	    }
	 
	    if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,arraySize,1,&(buffer[0])) == false) {
	       simClasses->logger << "(SEP OP WLM VOLUME) ERROR: Could not write phase-space volume on slice #" << slice << endl << write;
	       success = false;
	    }
	 }
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
