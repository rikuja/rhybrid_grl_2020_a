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

#include <constants.h>
#include "sep_coordinate_transform.h"
#include "sep_simcontrol.h"
#include "sep_operator_shock_mesh.h"
#include "sep_base_class_shock.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   const string shockMeshName = "shock";

   namespace shockmesh {
      void getBlockIDs(Simulation* sim,SimulationClasses* simClasses,const Real* nodeCoords,
		       pargrid::CellID& blockLID,pargrid::CellID& blockGID) {
	 const uint32_t i_block = nodeCoords[0] / block::WIDTH_X;
	 const uint32_t j_block = nodeCoords[1] / block::WIDTH_Y;
	 const uint32_t k_block = nodeCoords[2] / block::WIDTH_Z;
	 blockGID = block::calculateGlobalIndex(*sim,i_block,j_block,k_block);
	 
	 // Hack: for now it is safe to ignore this, but in future if plasma state 
	 // or field data is read from mesh, the local IDs need to be correct and 
	 // this DataOperator properly parallelized.
         #warning HACK: Intentionally returning incorrect cell local IDs, for now safe to ignore
	 blockLID = simClasses->pargrid.getLocalID(blockGID);
	 
	 #ifndef NDEBUG
	 /*if (blockLID == pargrid::INVALID_CELLID) {
	    cerr << "\t in sep_operator shock mesh: " << nodeCoords[0] << '\t' << nodeCoords[1] << '\t' << nodeCoords[2] << endl;
	    blockLID = 0;
	 }*/
	 #endif
      }

      void transformToCartesian(Simulation* sim,const Real* RESTRICT position,const Real* RESTRICT vectorIn,Real* vectorOut) {
	 const uint32_t I = static_cast<uint32_t>(position[0]);
	 const uint32_t J = static_cast<uint32_t>(position[1]);
	 const uint32_t K = static_cast<uint32_t>(position[2]);
	 Real pos[3];
	 pos[0] = sim->x_crds_node[I] + (position[0]-I)*sim->dx_cell[I];
	 pos[1] = sim->y_crds_node[J] + (position[1]-J)*sim->dy_cell[J];
	 pos[2] = sim->z_crds_node[K] + (position[2]-K)*sim->dz_cell[K];
	 transformVector(pos,vectorIn,vectorOut,simControl.coordinateSystem,sep::CARTESIAN);
      }
      
      bool writeVariable(Simulation* sim,SimulationClasses* simClasses,bool nodeCentered,
			 const std::vector<Real>& values,size_t vectorSize,const std::string& name) {
	 bool success = true;
	 const size_t N_points = values.size() / vectorSize;
	 
	 map<string,string> xmlAttributes;
	 xmlAttributes["mesh"] = shockMeshName;
	 xmlAttributes["units"] = "";
	 xmlAttributes["name"] = shockMeshName+"/"+name;
	 if (nodeCentered == true) xmlAttributes["centering"] = "node";
	 else xmlAttributes["centering"] = "zone";

         #warning Only writes data on master process
	 if (sim->mpiRank == sim->MASTER_RANK) {
	    if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,N_points,vectorSize,&(values[0])) == false) {
	       simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write variable '" << name << "'" << endl << write;
	       success = false;
	    }
	 } else {
	    if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,0,vectorSize,&(values[0])) == false) {
	       simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write variable '" << name << "'" << endl << write;
	       success = false;
	    }
	 }
	 
	 return success;
      }
   }
   
   OperatorShockMesh::OperatorShockMesh(): DataOperator() {
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorShockMesh::~OperatorShockMesh() { }
   
   bool OperatorShockMesh::finalize() {
      bool success = true;
      return success;
   }
   
   std::string OperatorShockMesh::getName() const {
      return "ShockMeshWriter";
   }
   
   bool OperatorShockMesh::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (DataOperator::initialize(cr,sim,simClasses) == false) success = false;
      return success;
   }
   
   bool OperatorShockMesh::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      if (simControl.includeShock == false) {
	 simClasses->logger << "(SEP SHOCK MESH OP) No shock in simulation, exiting." << endl << write;
	 return success;
      }
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Shock pointer is NULL while includeShock==true" << endl << write;
	 return false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::start("Shock Mesh",profTotalTime);
      #endif

      // Allocate shock mesh:
      if (simControl.shock->initializeMesh(sim->t,simControl.N_shockSurfaceRefinements) == false) {
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	 #endif
	 simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Shock mesh initialization failed" << endl << write;
	 return false;
      }
      
      // Write cell connectivity array:
      vector<uint32_t>& cellConnectivity = simControl.shock->getCellConnectivity();
      const uint32_t N_cells = simControl.shock->getNumberOfCells();
      const uint32_t N_nodes = simControl.shock->getNumberOfNodes();
      
      string N_cells_string,N_nodes_string;
      stringstream ss1,ss2;
      ss1 << N_cells;
      ss1 >> N_cells_string;
      ss2 << N_nodes;
      ss2 >> N_nodes_string;
      
      map<string,string> xmlAttributes;
      xmlAttributes["name"] = shockMeshName;
      xmlAttributes["type"] = vlsv::mesh::STRING_UCD_GENERIC_MULTI;
      xmlAttributes["domains"] = "1";
      xmlAttributes["cells"] = N_cells_string;
      xmlAttributes["nodes"] = N_nodes_string;

      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,cellConnectivity.size(),1,&(cellConnectivity[0])) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write cell connectivity" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,0,1,&(cellConnectivity[0])) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write cell connectivity" << endl << write;
	    success = false;
	 }
      }
      
      // Write node coordinates:
      xmlAttributes.clear();
      xmlAttributes["mesh"] = shockMeshName;
      
      if (sim->mpiRank == sim->MASTER_RANK) {
	 vector<Real>& nodeCoordinates = simControl.shock->getNodeCoordinates();
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS",xmlAttributes,N_nodes,3,&(nodeCoordinates[0])) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write node coordinates" << endl << write;
	    success = false;
	 }
      } else {
	 vector<Real> nodeCoordinates(1);
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS",xmlAttributes,0,1,&(nodeCoordinates[0])) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write node coordinates" << endl << write;
	    success = false;
	 }
      }
      
      // Write cell offsets:
      //vector<uint32_t>& cellOffsets = simControl.shock->getCellOffsets();
      //const uint32_t N_connectivityEntries = cellOffsets[cellOffsets.size()-1];
      uint32_t offsetEntries[2];
      //offsetEntries[0] = cellOffsets[cellOffsets.size()-1];
      offsetEntries[0] = cellConnectivity.size();
      offsetEntries[1] = N_nodes;
      
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_OFFSETS",xmlAttributes,1,2,offsetEntries) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write offsets" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH_OFFSETS",xmlAttributes,0,2,offsetEntries) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write offsets" << endl << write;
	    success = false;
	 }
      }
      
      // Write mesh bounding box (master process only):
      if (sim->mpiRank == sim->MASTER_RANK) {
	 unsigned int bbox[6];
	 bbox[0] = 1; // Number of blocks in x-direction in mesh bounding box.
	 bbox[1] = 1; // Number of blocks in y-direction in mesh bounding box.
	 bbox[2] = 1; // Number of blocks in z-direction in mesh bounding box.
	 bbox[3] = 1; // Number of cells in each block in x-direction.
	 bbox[4] = 1; // Number of cells in each block in y-direction.
	 bbox[5] = 1; // Number of cells in each block in z-direction.
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,6,1,bbox) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write mesh bounding box" << endl << write;
	    success = false;
	 }
      } else {
	 unsigned int* ptr = NULL;
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,0,1,ptr) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write mesh bounding box" << endl << write;
	    success = false;
	 }
      }
      
      // Write domain sizes:
      int domainSize[4];
      const int N_ghosts = 0;
      domainSize[0] = N_cells;
      domainSize[1] = N_ghosts;
      domainSize[2] = N_nodes;
      domainSize[3] = 0;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,4,domainSize) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,0,4,domainSize) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }
      }

      // Deallocate shock mesh:
      if (simControl.shock->finalizeMesh() == false) {
	 simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to finalize shock mesh" << endl << write;
	 success = false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
