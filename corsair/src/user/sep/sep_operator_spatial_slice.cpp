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
#include <sstream>

#include "sep_simcontrol.h"
#include "sep_operator_spatial_slice.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const std::string prefix = "SpatialSliceOperator";

   SpatialSliceOP::SpatialSliceOP(): DataOperator() { 
      baseClassInitialized = false;
   }
   
   SpatialSliceOP::~SpatialSliceOP() { 
      finalize();
   }

   bool SpatialSliceOP::addConfigFileItems(ConfigReader& cr) {
      cr.addComposed(prefix+".sliced_coordinate","Sliced coordinate: x/y/z (string).");
      cr.addComposed(prefix+".slice_origin","Origin of 2D spatial mesh slice (float).");
      cr.addComposed(prefix+".slice_geometry","Geometry of output mesh cartesian/cylindrical (string)");
      return true;
   }

   uint32_t SpatialSliceOP::calculateNewGlobalID(uint32_t* indices) const {
      return indices[2]*y_blocks*x_blocks + indices[1]*x_blocks + indices[0];
   }
   
   bool SpatialSliceOP::finalize() {
      bool success = true;
      return success;
   }
   
   void SpatialSliceOP::getAcceptedBlocks(uint32_t slice,std::vector<pargrid::CellID>& blockGIDs,std::vector<pargrid::CellID>& blockLIDs) {
      blockLIDs.clear();
      blockGIDs.clear();
      if (slice >= getNumberOfSlices()) return;
      
      const vector<pargrid::CellID>& globalIDs = simClasses->pargrid.getGlobalIDs();
      const double* const blockCoordinateArray = getBlockCoordinateArray(*sim,*simClasses);
      
      Real blockCoords[3];
      uint32_t blockIndices[3];
      Real blockSizes[3];
      for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 // Calculate spatial block indices:
	 block::calculateBlockIndices(*sim,globalIDs[blockLID],blockIndices[0],blockIndices[1],blockIndices[2]);
	 
	 // Get block coordinates:
	 blockCoords[0] = blockCoordinateArray[blockLID*3+0];
	 blockCoords[1] = blockCoordinateArray[blockLID*3+1];
	 blockCoords[2] = blockCoordinateArray[blockLID*3+2];
	 
	 // Get block size:
	 blockSizes[0] = sim->dx_block[blockIndices[0]];
	 blockSizes[1] = sim->dy_block[blockIndices[1]];
	 blockSizes[2] = sim->dz_block[blockIndices[2]];
	 
	 // Check if slice intersects block:
	 const Real sliceOrigin = sliceOrigins[slice];
	 const Real CRD = blockCoords[sliceCoordinates[slice]];
	 const Real DX  = blockSizes[sliceCoordinates[slice]];

	 if (CRD > sliceOrigin || CRD+DX < sliceOrigin) continue;
	 sliceIndices[slice] = blockIndices[sliceCoordinates[slice]];

	 if (isCylindrical(slice) == true) {
	    blockIndices[sliceCoordinates[slice]] = blockIndices[2];
	 }

	 // Calculate new global ID for cells in the accepted block:
	 for (uint32_t cell=0; cell<N_cells; ++cell) {
	    if (isCylindrical(slice) == true) {
	       blockIndices[2] = cell;
	    } else {
	       blockIndices[sliceCoordinates[slice]] = cell;
	    }
	    blockGIDs.push_back( calculateNewGlobalID(blockIndices) );
	 }
	 blockLIDs.push_back( blockLID );
      }
   }

   size_t SpatialSliceOP::getNumberOfSlices() const {return sliceCoordinates.size();}
   
   uint8_t SpatialSliceOP::getSlicedCoordinate(uint32_t slice) const {
      return sliceCoordinates[slice];
   }
   
   /*std::string SpatialSliceOP::getSliceName(uint32_t slice) const {
      // Create unique name for the slice:
      string meshName;
      stringstream ss;
      ss << "WaveMesh" << slice;
      ss >> meshName;
      return meshName;
   }*/

   uint32_t SpatialSliceOP::getSliceIndex(uint32_t slice) const {
      return sliceIndices[slice];
   }
   
   Real SpatialSliceOP::getSliceOrigin(uint32_t slice) const {
      return sliceOrigins[slice];
   }
   
   bool SpatialSliceOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      if (baseClassInitialized == true) return baseClassInitialized;
      baseClassInitialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      // Add config file items:
      if (addConfigFileItems(cr) == false) {
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      vector<string> slicedCoordinateStrings;
      vector<string> sliceOriginStrings;
      vector<string> sliceGeometryString;
      cr.parse();
      cr.get(prefix+".sliced_coordinate",slicedCoordinateStrings);
      cr.get(prefix+".slice_origin",sliceOriginStrings);
      cr.get(prefix+".slice_geometry",sliceGeometryString);
      
      // Check that config file contains enough parameters for all slices:
      if (slicedCoordinateStrings.size() != sliceOriginStrings.size()) {
	 simClasses.logger << "(OP SPATIAL SLICE) ERROR: Number of sliced coordinates and number of slice origins differ." << endl << write;
	 baseClassInitialized = false; return baseClassInitialized;
      }
      if (slicedCoordinateStrings.size() != sliceGeometryString.size()) {
	 simClasses.logger << "(OP SPATIAL SLICE) ERROR: Number of sliced coordinates and number of slice geometries." << endl << write;
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      for (size_t i=0; i<slicedCoordinateStrings.size(); ++i) {
	 if (slicedCoordinateStrings[i] == "x") sliceCoordinates.push_back(0);
	 else if (slicedCoordinateStrings[i] == "y") sliceCoordinates.push_back(1);
	 else if (slicedCoordinateStrings[i] == "z") sliceCoordinates.push_back(2);
	 else {
	    simClasses.logger << "(OP SPATIAL SLICE) ERROR: Unknown slice coordinate '" << slicedCoordinateStrings[i];
	    simClasses.logger << "' given in config file, should be one of 'x', 'y' or 'z'." << endl << write;
	    baseClassInitialized = false; return baseClassInitialized;
	 }
      }
      
      for (size_t i=0; i<sliceGeometryString.size(); ++i) {
	 bool ok = false;
	 if (sliceGeometryString[i] == "cartesian") {ok = true; sliceGeometries.push_back(false);}
	 if (sliceGeometryString[i] == "cylindrical") {ok = true; sliceGeometries.push_back(true);}

	 if (ok == false) {
	    simClasses.logger << "(OP SPATIAL SLICE) ERROR: Unknown slice geometry '" << sliceGeometryString[i] << "'";
	    simClasses.logger << "' given in config file, should be 'cartesian' or 'cylindrical'." << endl << write;
	    baseClassInitialized = false; return baseClassInitialized;
	 }
      }
      
      for (size_t i=0; i<sliceOriginStrings.size(); ++i)
	sliceOrigins.push_back(atof(sliceOriginStrings[i].c_str()));

      sliceIndices.resize(sliceOrigins.size());
      for (size_t i=0; i<sliceIndices.size(); ++i) {
	 sliceIndices[i] = numeric_limits<uint32_t>::max();
      }

      return baseClassInitialized;
   }
   
   bool SpatialSliceOP::isCylindrical(uint32_t slice) const {
      return sliceGeometries[slice];
   }

   void SpatialSliceOP::prepareSlice(uint32_t slice,uint32_t N_cells) {
      this->N_cells = N_cells;
      
      if (slice >= getNumberOfSlices()) {
	 x_blocks = numeric_limits<uint32_t>::max();
	 y_blocks = numeric_limits<uint32_t>::max();
	 z_blocks = numeric_limits<uint32_t>::max();
	 return;
      }
      
      uint32_t bbox[3];
      bbox[0] = sim->x_blocks;  // Number of mesh blocks in x-direction in mesh bounding box.
      bbox[1] = sim->y_blocks;  // Number of mesh blocks in y-direction in mesh bounding box.
      bbox[2] = sim->z_blocks;  // Number of mesh blocks in z-direction in mesh bounding box.
      if (isCylindrical(slice) == true) {
	 bbox[getSlicedCoordinate(slice)+0] = bbox[2];
	 bbox[2] = N_cells;
      } else {
	 bbox[getSlicedCoordinate(slice)+0] = N_cells;
      }

      x_blocks = bbox[0];
      y_blocks = bbox[1];
      z_blocks = bbox[2];
   }

} // namespace sep

