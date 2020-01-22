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
#include "sep_operator_spatial_lineout.h"
#include "sep_mesh_logical.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const std::string prefix = "SpatialLineoutOperator";

   OperatorSpatialLineout::OperatorSpatialLineout(): DataOperator() { 
      baseClassInitialized = false;
   }
   
   OperatorSpatialLineout::~OperatorSpatialLineout() { 
      finalize();
   }

   bool OperatorSpatialLineout::addConfigFileItems(ConfigReader& cr) {
      cr.addComposed(prefix+".lineout_coordinate","Coordinate direction of lineout: x/y/z (string).");
      cr.addComposed(prefix+".lineout_origin_x","Line origin, x-coordinate (float).");
      cr.addComposed(prefix+".lineout_origin_y","Line origin, y-coordinate (float).");
      cr.addComposed(prefix+".lineout_origin_z","Line origin, z-coordinate (float).");
      return true;
   }
   
   bool OperatorSpatialLineout::finalize() {
      bool success = true;
      return success;
   }
   
   bool OperatorSpatialLineout::getAcceptedBlocks(size_t line,std::vector<pargrid::CellID>& inner,std::vector<pargrid::CellID>& boundary) {
      bool success = true;
      if (line >= lineoutCoordinates.size()) return false;
      if (lineoutCoordinates[line] >= 3) return false;

      // Get logical coordinates of line:
      Real lineCoords[3];
      lineCoords[0] = lineoutOriginsX[line];
      lineCoords[1] = lineoutOriginsY[line];
      lineCoords[2] = lineoutOriginsZ[line];
      
      uint32_t lineIndices[3];
      for (int i=0; i<3; ++i) lineIndices[i] = numeric_limits<uint32_t>::max();
      if (lineoutCoordinates[line] != 0) lineIndices[0] = static_cast<uint32_t>(lineCoords[0]) / block::WIDTH_X;
      if (lineoutCoordinates[line] != 1) lineIndices[1] = static_cast<uint32_t>(lineCoords[1]) / block::WIDTH_Y;
      if (lineoutCoordinates[line] != 2) lineIndices[2] = static_cast<uint32_t>(lineCoords[2]) / block::WIDTH_Z;
      
      const vector<pargrid::MPI_processID>& hosts = simClasses->pargrid.getHosts();
      for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 // Calculate block indices:
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 
	 // Check if lineout intersects the block:
	 if (lineoutCoordinates[line] != 0) if (i_block != lineIndices[0]) continue;
	 if (lineoutCoordinates[line] != 1) if (j_block != lineIndices[1]) continue;
	 if (lineoutCoordinates[line] != 2) if (k_block != lineIndices[2]) continue;
	 
	 // Block is accepted, figure out whether it is inner or boundary block:
	 pargrid::CellID* nbrs = simClasses->pargrid.getCellNeighbourIDs(blockLID);
	 	 
	 pargrid::CellID posNeighborID,negNeighborID;
	 switch (lineoutCoordinates[line]) {
	  case 0:
	    negNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,0,0)];
	    posNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,0,0)];
	    break;
	  case 1:
	    negNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(0,-1,0)];
	    posNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(0,+1,0)];
	    break;
	  case 2:
	    negNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(0,0,-1)];
	    posNeighborID = nbrs[simClasses->pargrid.calcNeighbourTypeID(0,0,+1)];
	    break;
	  default:
	    negNeighborID = pargrid::INVALID_CELLID;
	    posNeighborID = pargrid::INVALID_CELLID;
	    break;
	 }
	 
	 // Block is boundary block if either of its neighbor to the direction 
	 // of lineout is hosted on another MPI process:
	 bool isBoundary = false;
	 if (negNeighborID != pargrid::INVALID_CELLID) if (hosts[negNeighborID] != sim->mpiRank) isBoundary = true;
	 if (posNeighborID != pargrid::INVALID_CELLID) if (hosts[posNeighborID] != sim->mpiRank) isBoundary = true;
	 if (isBoundary == true) boundary.push_back(blockLID);
	 else inner.push_back(blockLID);
      }

      return success;
   }
   
   uint8_t OperatorSpatialLineout::getLineoutCoordinate(size_t line) const {
      if (line >= lineoutCoordinates.size()) return numeric_limits<uint8_t>::max();
      return lineoutCoordinates[line];
   }
   
   std::string OperatorSpatialLineout::getLineoutName(size_t line) const {
      if (line >= lineoutCoordinates.size()) return "";
      stringstream ss;
      ss << "lineout" << line;
      return ss.str();
   }

   size_t OperatorSpatialLineout::getNumberOfLineouts() const {return lineoutCoordinates.size();}
   
   bool OperatorSpatialLineout::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      if (baseClassInitialized == true) return baseClassInitialized;
      baseClassInitialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      // Add config file items:
      if (addConfigFileItems(cr) == false) {
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      vector<string> lineoutCoordinateStrings;
      vector<string> lineoutOriginStringsX;
      vector<string> lineoutOriginStringsY;
      vector<string> lineoutOriginStringsZ;
      cr.parse();
      cr.get(prefix+".lineout_coordinate",lineoutCoordinateStrings);
      cr.get(prefix+".lineout_origin_x",lineoutOriginStringsX);
      cr.get(prefix+".lineout_origin_y",lineoutOriginStringsY);
      cr.get(prefix+".lineout_origin_z",lineoutOriginStringsZ);
      
      if (lineoutCoordinateStrings.size() != lineoutOriginStringsX.size()) baseClassInitialized = false;
      if (lineoutCoordinateStrings.size() != lineoutOriginStringsY.size()) baseClassInitialized = false;
      if (lineoutCoordinateStrings.size() != lineoutOriginStringsZ.size()) baseClassInitialized = false;
      if (baseClassInitialized == false) {
	 simClasses.logger << "(SEP OP SPATIAL LINEOUT) ERROR: Number of lineout coordinates and number of lineout origins differ." << endl << write;
	 baseClassInitialized = false; return baseClassInitialized;
      }
      
      for (size_t i=0; i<lineoutCoordinateStrings.size(); ++i) {
	 if (lineoutCoordinateStrings[i] == "x") lineoutCoordinates.push_back(0);
	 else if (lineoutCoordinateStrings[i] == "y") lineoutCoordinates.push_back(1);
	 else if (lineoutCoordinateStrings[i] == "z") lineoutCoordinates.push_back(2);
	 else {
	    simClasses.logger << "(SEP OP SPATIAL LINEOUT) ERROR: Unknown slice coordinate '" << lineoutCoordinateStrings[i];
	    simClasses.logger << "' given in config file, should be one of 'x', 'y' or 'z'." << endl << write;
	    baseClassInitialized = false; return baseClassInitialized;
	 }
      }

      // Convert origin coordinates to logical units:
      for (size_t i=0; i<lineoutOriginStringsX.size(); ++i) {
	 Real phys[3];
	 Real logical[3];
	 phys[0] = atof(lineoutOriginStringsX[i].c_str());
	 phys[1] = atof(lineoutOriginStringsY[i].c_str());
	 phys[2] = atof(lineoutOriginStringsZ[i].c_str());
	 
	 // Check that coordinates are valid:
	 if (lineoutCoordinates[i] != 0) {
	    if (phys[0] < sim.x_crds_node[0] || phys[0] >= sim.x_crds_node[sim.x_blocks*block::WIDTH_X]) baseClassInitialized = false;
	    else phys[0] = sim.x_crds_node[0] + 0.5*sim.dx_cell[0];
	 }
	 if (lineoutCoordinates[i] != 1) {
	    if (phys[1] < sim.y_crds_node[0] || phys[1] >= sim.y_crds_node[sim.y_blocks*block::WIDTH_Y]) baseClassInitialized = false;
	    else {
	       //phys[1] = sim.y_crds_node[0] + 0.5*sim.dy_cell[0];
	    }
	 }
	 if (lineoutCoordinates[i] != 2) {
	    if (phys[2] < sim.z_crds_node[0] || phys[2] >= sim.z_crds_node[sim.z_blocks*block::WIDTH_Z]) baseClassInitialized = false;
	    else {
	       //phys[2] = sim.z_crds_node[0] + 0.5*sim.dz_cell[0];
	    }
	 }
	 if (baseClassInitialized == false) {
	    simClasses.logger << "(SEP OP SPATIAL LINEOUT) ERROR: Origin coordinates for lineout #" << i << " are invalid" << endl;
	    simClasses.logger << "x_min/max: " << sim.x_crds_node[0] << '\t' << sim.x_crds_node[sim.x_blocks*block::WIDTH_X] << endl;
	    simClasses.logger << "y_min/max: " << sim.y_crds_node[0] << '\t' << sim.y_crds_node[sim.y_blocks*block::WIDTH_Y] << endl;
	    simClasses.logger << "z_min/max: " << sim.z_crds_node[0] << '\t' << sim.z_crds_node[sim.z_blocks*block::WIDTH_Z] << endl;
	    simClasses.logger << "Given origin: " << phys[0] << '\t' << phys[1] << '\t' << phys[2] << endl;
	    simClasses.logger << write;
	    break;
	 }

	 getLogicalCoordinates(&sim,phys,logical);
	 lineoutOriginsX.push_back(logical[0]);
	 lineoutOriginsY.push_back(logical[1]);
	 lineoutOriginsZ.push_back(logical[2]);
	 
	 if (logical[0] == numeric_limits<Real>::infinity()) {
	    simClasses.logger << "(SEP OP SPATIAL LINEOUT) ERROR: Invalid logical origin coordinates" << endl;
	    simClasses.logger << "x,y,z coords:  " << phys[0] << '\t' << phys[1] << '\t' << phys[2] << endl;
	    simClasses.logger << "x,y,z logical: " << logical[0] << '\t' << logical[1] << '\t' << logical[2] << endl << write;
	    baseClassInitialized = false;
	    break;
	 }
      }
      
      return baseClassInitialized;
   }
      
} // namespace sep

