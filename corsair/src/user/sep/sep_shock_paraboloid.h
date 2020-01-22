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

#ifndef SEP_SHOCK_PARABOLOID_H
#define SEP_SHOCK_PARABOLOID_H

#include <vector>

#include <vlsv_common.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace sep {
   /*
   struct UCDMeshNode {
      vlsv::celltype::type cellType;
      uint32_t* cell;
      uint32_t* face0;
      uint32_t* face1;
      uint32_t* face2;
      Real coords[3];
      uint16_t sizes[4];
      
      UCDMeshNode();
      UCDMeshNode(const UCDMeshNode& mn);
      ~UCDMeshNode();      
   };
   */
   class ShockParaboloid {
    public:
      ShockParaboloid();
      ~ShockParaboloid();

      bool finalize();
      std::vector<uint32_t>& getCellConnectivity();
      std::vector<uint32_t>& getCellOffsets();
      std::vector<Real>& getNodeCoordinates();
      uint32_t getNumberOfCells() const;
      uint32_t getNumberOfNodes() const;
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      //bool writeMesh(Simulation& sim,SimulationClasses& simClasses,const std::string& meshName,bool writeMesh);
      
    private:
      Real maxHeight;
      Real height0;
      Real radius0_x;
      Real radius0_y;
      Real x0;
      uint32_t N_nodes_u;
      uint32_t N_nodes_v;
      
      //std::vector<UCDMeshNode> nodes;
      std::vector<uint32_t> cellConnectivity;
      std::vector<uint32_t> cellOffsets;
      std::vector<Real> nodeCoordinates;
      
      bool createNodes();
      uint32_t globalIndex(uint32_t i,uint32_t j);
   };
   
} // namespace sep

#endif
