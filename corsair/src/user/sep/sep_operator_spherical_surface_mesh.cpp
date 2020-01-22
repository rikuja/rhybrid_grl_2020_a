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
#include <map>
#include <vector>

#include <triangulated_sphere.h>

#include "sep_simcontrol.h"
#include "sep_operator_spherical_surface_mesh.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   OperatorSphericalSurfaceMesh::OperatorSphericalSurfaceMesh(): DataOperator() {
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorSphericalSurfaceMesh::~OperatorSphericalSurfaceMesh() { }
   
   bool OperatorSphericalSurfaceMesh::finalize() {
      bool success = true;
      return success;
   }

   std::string OperatorSphericalSurfaceMesh::getName() const {
      return "SphericalSurfaceMeshWriter";
   }
   
   bool OperatorSphericalSurfaceMesh::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (DataOperator::initialize(cr,sim,simClasses) == false) success = false;
      radius0 = constants::DIST_SOLAR_RADIUS;
      return success;
   }
   
   bool OperatorSphericalSurfaceMesh::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         profile::start("Spherical Surface Mesh",profTotalTime);
      #endif

      vector<ucdmesh::Node> nodes;
      vector<ucdmesh::TriangularFace> faces;
      Real centroid[3] = {0.0,0.0,0.0};
      const int refinementLevel = 2;
      size_t N_connections;
      ucdmesh::createSphere(ucdmesh::spheregenerator::icosahedron,refinementLevel,centroid,radius0,nodes,faces,N_connections);
      
      stringstream ss1;
      ss1 << nodes.size();
      stringstream ss2;
      ss2 << faces.size();
      
      string N_cells_string,N_nodes_string;
      ss1 >> N_nodes_string;
      ss2 >> N_cells_string;
      
      map<string,string> xmlAttributes;
      xmlAttributes["name"] = "Sun";
      xmlAttributes["type"] = vlsv::mesh::STRING_UCD_GENERIC_MULTI;
      xmlAttributes["domains"] = "1";
      xmlAttributes["cells"] = N_cells_string;
      xmlAttributes["nodes"] = N_nodes_string;

      const uint32_t connectivitySize = N_connections;
      
      uint32_t* facePointer = &(faces[0].cellType);
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,connectivitySize,1,facePointer) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write cell connectivity" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,0,1,facePointer) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write cell connectivity" << endl << write;
	    success = false;
	 }
      }
      
      // Write node coordinates:
      xmlAttributes.clear();
      xmlAttributes["mesh"] = "Sun";
      
      Real* nodePointer = &(nodes[0].x);
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS",xmlAttributes,nodes.size(),3,nodePointer) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write node coordinates" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS",xmlAttributes,0,1,nodePointer) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write node coordinates" << endl << write;
	    success = false;
	 }
      }
      
      // Write cell offsets:
      uint32_t offsetEntries[vlsv::ucdgenericmulti::offsets::SIZE];
      offsetEntries[vlsv::ucdgenericmulti::offsets::ZONE_ENTRIES] = connectivitySize;
      offsetEntries[vlsv::ucdgenericmulti::offsets::NODE_ENTRIES] = nodes.size();
      
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_OFFSETS",xmlAttributes,1,vlsv::ucdgenericmulti::offsets::SIZE,offsetEntries) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write offsets" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH_OFFSETS",xmlAttributes,0,vlsv::ucdgenericmulti::offsets::SIZE,offsetEntries) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write offsets" << endl << write;
	    success = false;
	 }
      }
      
      // Write mesh bounding box (master process only):
      if (sim->mpiRank == sim->MASTER_RANK) {
	 unsigned int bbox[vlsv::ucdgenericmulti::bbox::SIZE];
	 bbox[vlsv::ucdgenericmulti::bbox::X_BLOCKS] = 1; // Number of blocks in x-direction in mesh bounding box.
	 bbox[vlsv::ucdgenericmulti::bbox::Y_BLOCKS] = 1; // Number of blocks in y-direction in mesh bounding box.
	 bbox[vlsv::ucdgenericmulti::bbox::Z_BLOCKS] = 1; // Number of blocks in z-direction in mesh bounding box.
	 bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X] = 1; // Number of cells in each block in x-direction.
	 bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y] = 1; // Number of cells in each block in y-direction.
	 bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z] = 1; // Number of cells in each block in z-direction.
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,vlsv::ucdgenericmulti::bbox::SIZE,1,bbox) == false) {
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
      int domainSize[vlsv::ucdgenericmulti::domainsizes::SIZE];
      const int N_ghosts = 0;
      domainSize[vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS] = faces.size();
      domainSize[vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS] = N_ghosts;
      domainSize[vlsv::ucdgenericmulti::domainsizes::TOTAL_NODES] = nodes.size();
      domainSize[vlsv::ucdgenericmulti::domainsizes::GHOST_NODES] = 0;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,vlsv::ucdgenericmulti::domainsizes::SIZE,domainSize) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }
      } else {
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,0,vlsv::ucdgenericmulti::domainsizes::SIZE,domainSize) == false) {
	    simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
