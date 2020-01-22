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
#include "sep_operator_wavelength_mesh.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   
   WavelengthMeshOP::WavelengthMeshOP(): SpatialSliceOP() { 
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   WavelengthMeshOP::~WavelengthMeshOP() { }
   
   bool WavelengthMeshOP::finalize() {
      bool success = true;
      return success;
   }
   
   std::string WavelengthMeshOP::getName() const {
      return "WavelengthMeshWriter";
   }
   
   bool WavelengthMeshOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (SpatialSliceOP::initialize(cr,sim,simClasses) == false) return false;
      return success;
   }
   
   bool WavelengthMeshOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      if (baseClassInitialized == false) {
	 simClasses->logger << "(SEP OP WLM) ERROR: Base class is not initialized" << endl << write;
	 return false;
      }

      #if PROFILE_LEVEL > 0
         profile::start("WLM Slices",profTotalTime);
      #endif
      for (size_t slice=0; slice<getNumberOfSlices(); ++slice) {
	 // Tell spatial slice operator that we're processing new slice:
	 prepareSlice(slice,simControl.N_wavelengthMeshCells);

	 // Create a name for the wavelength mesh slice:
	 stringstream ss;
	 ss << "WaveMesh" << slice;
	 const string meshName = ss.str();

	 // Get sliced coordinate (0-2):
	 const uint8_t slicedCoordinate = getSlicedCoordinate(slice);
	 
	 // Create mesh bounding box:
	 uint32_t bbox[6];
	 bbox[0] = sim->x_blocks;  // Number of mesh blocks in x-direction in mesh bounding box.
	 bbox[1] = sim->y_blocks;  // Number of mesh blocks in y-direction in mesh bounding box.
	 bbox[2] = sim->z_blocks;  // Number of mesh blocks in z-direction in mesh bounding box.
	 bbox[3] = block::WIDTH_X; // Number of cells in each block in x-direction.
	 bbox[4] = block::WIDTH_Y; // Number of cells in each block in y-direction.
	 bbox[5] = block::WIDTH_Z; // Number of cells in each block in z-direction.

	 if (isCylindrical(slice) == true) {
	    bbox[slicedCoordinate+0] = bbox[2];
	    bbox[slicedCoordinate+3] = bbox[5];
	    bbox[2] = simControl.N_wavelengthMeshCells;
	    bbox[5] = 1;
	 } else {
	    bbox[slicedCoordinate+0] = simControl.N_wavelengthMeshCells;
	    bbox[slicedCoordinate+3] = 1;
	 }

	 // Get number of nodes:
	 uint32_t N_nodes[3];
	 N_nodes[0] = sim->x_blocks*block::WIDTH_X + 1;
	 N_nodes[1] = sim->y_blocks*block::WIDTH_Y + 1;
	 N_nodes[2] = sim->z_blocks*block::WIDTH_Z + 1;
	 
	 if (isCylindrical(slice) == true) {
	    N_nodes[slicedCoordinate] = N_nodes[2];
	    N_nodes[2] = simControl.N_wavelengthMeshCells+1;
	 } else {
	    N_nodes[slicedCoordinate] = simControl.N_wavelengthMeshCells+1;
	 }
	 
	 // Create array containing wavelength mesh coordinates:
	 Real* wmesh = new Real[simControl.N_wavelengthMeshCells+1];
	 for (size_t i=0; i<simControl.N_wavelengthMeshCells+1; ++i) {
	    Real lambda = simControl.wavelengthMeshNodeCoordinates[i];
	    if (lambda < 0.0) {
	       wmesh[i] = min(0.0,-log10(-lambda));
	    } else {
	       wmesh[i] = max(0.0,+log10(lambda));
	    }
	 }

	 // Get pointers to mesh node coordinates:
	 Real* crds[3];
	 crds[0] = sim->x_crds_node;
	 crds[1] = sim->y_crds_node;
	 crds[2] = sim->z_crds_node;
	 
	 if (isCylindrical(slice) == true) {
	    crds[slicedCoordinate] = crds[2];
	    crds[2] = wmesh;
	 } else {
	    crds[slicedCoordinate] = wmesh;
	 }

	 uint64_t arraySize = 6;
	 uint32_t* ptr = bbox;
	 if (sim->mpiRank != sim->MASTER_RANK) {
	    arraySize = 0;
	    ptr = NULL;
	    for (int i=0; i<3; ++i) crds[i]    = NULL;
	    for (int i=0; i<3; ++i) N_nodes[i] = 0;
	 }

	 vector<pargrid::CellID> validLocalIDs;
	 vector<pargrid::CellID> newGlobalIDs;
	 getAcceptedBlocks(slice,newGlobalIDs,validLocalIDs);

	 if (newGlobalIDs.size() == 0) {
	    simClasses->logger << "(SEP OP WAVELENGTH MESH) WARNING: Slice has no valid cells, slice: " << slice;
	    simClasses->logger << " sliced coordinate: " << (int)slicedCoordinate << endl;
	    simClasses->logger << "\t slice origin: " << getSliceOrigin(slice) << endl;
	 }

	 // Write mesh bounding box and node coordinates (master process only):
	 map<string,string> xmlAttributes;
	 xmlAttributes["mesh"] = meshName;
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,arraySize,1,ptr) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttributes,N_nodes[0],1,crds[0]) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttributes,N_nodes[1],1,crds[1]) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttributes,N_nodes[2],1,crds[2]) == false) success = false;

	 #warning IMPROVE ME Add ghost cells to wavelength mesh
	 const uint64_t N_blocks = newGlobalIDs.size(); // This is locals + ghosts
	 const uint64_t N_ghosts = 0;

	 // Write mesh global IDs:
	 xmlAttributes.clear();
	 xmlAttributes["name"] = meshName;
	 xmlAttributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
	 stringstream sss;
	 sss << (pargrid::MPI_processID)sim->mpiProcesses;
	 xmlAttributes["domains"] = sss.str();

	 if (isCylindrical(slice) == true) {
	    xmlAttributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 } else {
	    xmlAttributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 }

	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,N_blocks,1,&(newGlobalIDs[0])) == false) {
	    simClasses->logger << "(OP WAVELENGTH MESH) ERROR: Failed to write block global IDs" << endl << write;
	    success = false;
	 }

	 // Write domain size:
	 xmlAttributes.clear();
	 xmlAttributes["mesh"] = meshName;
	 int domainSize[2];
	 domainSize[0] = N_blocks;
	 domainSize[1] = N_ghosts;
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,2,domainSize) == false) {
	    simClasses->logger << "(OP WAVELENGTH MESH) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }

	 vector<uint32_t> validNeighbors;
	 vector<pargrid::MPI_processID> validHosts;

	 // Write ghost block local IDs:
	 if (simClasses->vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttributes,N_ghosts,1,&(validNeighbors[0])) == false) {
	    simClasses->logger << "(OP WAVELENGTH MESH) ERROR: Failed to write ghost block local IDs!" << endl << write;
	    success = false;
	 }

	 // Write ghost block domain IDs:
	 if (simClasses->vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttributes,N_ghosts,1,&(validHosts[0])) == false) {
	    simClasses->logger << "(OP WAVELENGTH MESH) ERROR: Failed to write block ghost domain IDs!" << endl << write;
	    success = false;
	 }

	 delete [] wmesh; wmesh = NULL;
      }
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
