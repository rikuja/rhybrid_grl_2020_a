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
#include <cmath>

#include "sep_mesh_logical.h"
#include "sep_simcontrol.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   uint32_t binarySearchCellIndex(const Real* nodeCoords,uint32_t N_cells,Real x) {
      // Define search range (whole node coordinates array):
      uint32_t i_min = 0;
      uint32_t i_max = N_cells;

      // Check that given coordinate is actually inside the search range:
      if (x < nodeCoords[i_min]) return numeric_limits<uint32_t>::max();
      if (x > nodeCoords[i_max]) return numeric_limits<uint32_t>::max();

      // Iterate until given physical coordinate is between node i and node i+1:
      while (i_max > i_min+1) {
	 const uint32_t i_middle = (i_min+i_max)/2;
	 if (x < nodeCoords[i_middle]) i_max = i_middle;
	 else i_min = i_middle;
      }

      #ifndef NDEBUG
         if (i_max != i_min+1) {
	    cerr << "(SEP MESH LOGICAL) ERROR: i_max != i_min+1 in binarySearchCellIndex" << endl;
	    exit(1);
	 }
         if (nodeCoords[i_min] > x || nodeCoords[i_max] < x) {
	    cerr << "(SEP MESH LOGICAL) ERROR: x " << x << " is not in range [" << nodeCoords[i_min];
	    cerr << '\t' << nodeCoords[i_max] << "]" << endl;
	    exit(1);
	 }
      #endif
      
      return i_min;
   }
   
   /** Get logical coordinates corresponsing to given physical coordinates. If 
    * physical coordinates are outside simulation domain, all logical coordinates 
    * will contain value 'inf'.
    * @param sim Generic simulation control variables.
    * @param physical Physical coordinates in simulation coordinate system.
    * @param logical Array where logical coordinates are written to.*/
   void getLogicalCoordinates(Simulation* sim,const Real* RESTRICT physical,Real* logical) {
      uint32_t index,N_cells;

      /* //#ifndef NDEBUG
      bool success=true;
      if (physical[0] < sim->x_crds_node[0] || physical[0] > sim->x_crds_node[sim->x_blocks*block::WIDTH_X]) success=false;
      if (physical[1] < sim->y_crds_node[0] || physical[1] > sim->y_crds_node[sim->y_blocks*block::WIDTH_Y]) success=false;
      if (physical[2] < sim->z_crds_node[0] || physical[2] > sim->z_crds_node[sim->z_blocks*block::WIDTH_Z]) success=false;
      if (success == false) {
	 cerr << "(SEP MESH LOGICAL) ERROR: Invalid physical position." << endl;
	 cerr << "pos: ";
	 for (int i=0; i<3; ++i) cerr << physical[i] << '\t';
	 cerr << endl;
	 cerr << "x limits: " << sim->x_crds_node[0] << '\t' << sim->x_crds_node[sim->x_blocks*block::WIDTH_X] << std::endl;
	 cerr << "y limits: " << sim->y_crds_node[0] << '\t' << sim->y_crds_node[sim->y_blocks*block::WIDTH_Y] << std::endl;
	 cerr << "z limits: " << sim->z_crds_node[0] << '\t' << sim->z_crds_node[sim->z_blocks*block::WIDTH_Z] << std::endl;
	 exit(1);
      }
      //#endif */
      
      // Find logical x-coordinate:
      N_cells = sim->x_blocks * block::WIDTH_X;
      index = binarySearchCellIndex(sim->x_crds_node,N_cells,physical[0]);
      if (index == numeric_limits<uint32_t>::max()) {
	 for (int i=0; i<3; ++i) logical[i] = numeric_limits<Real>::infinity();
	 //cerr << "x inf " << sim->x_crds_node[0] << ' ' << sim->x_crds_node[1] << ' ' << physical[0] << endl;
	 //exit(1);
	 return;
      } else {
	 logical[0] = index + (physical[0] - sim->x_crds_node[index]) / sim->dx_cell[index];
      }
      
      // Find logical y-coordinate:
      N_cells = sim->y_blocks * block::WIDTH_Y;
      index = binarySearchCellIndex(sim->y_crds_node,N_cells,physical[1]);
      if (index == numeric_limits<uint32_t>::max()) {
	 for (int i=0; i<3; ++i) logical[i] = numeric_limits<Real>::infinity();
	 //cerr << "y inf " << sim->y_crds_node[0] << ' ' << sim->y_crds_node[1] << ' ' << physical[1] << endl;
	 //exit(1);
	 return;
      } else {
	 logical[1] = index + (physical[1] - sim->y_crds_node[index]) / sim->dy_cell[index];
      }
      
      // Find logical z-coordinate:
      N_cells = sim->z_blocks * block::WIDTH_Z;
      index = binarySearchCellIndex(sim->z_crds_node,N_cells,physical[2]);
      if (index == numeric_limits<uint32_t>::max()) {
	 for (int i=0; i<3; ++i) logical[i] = numeric_limits<Real>::infinity();
	 //cerr << "z inf " << sim->z_crds_node[0] << ' ' << sim->z_crds_node[1] << ' ' << physical[2] << endl;
	 //cerr << physical[0] << '\t' << physical[1] << '\t' << physical[2] << endl;
	 //exit(1);
	 return;
      } else {
	 logical[2] = index + (physical[2] - sim->z_crds_node[index]) / sim->dz_cell[index];
      }
   }
   
   Real lambdaScalingInvRadius(const Real* pos) {
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 return 1.0;
	 break;
       case sep::CARTESIAN:
	 return 1.0;
	 break;
       case sep::CYLINDRICAL:
	 return simControl.R_reference/pos[0];
	 break;
       case sep::SPHERICAL:
	 return simControl.R_reference/pos[0];
	 break;
       default:
	 return 1.0;
	 break;
      }
   }
   
   Real lambdaScalingNone(const Real* pos) {
      return 1.0;
   }

   Real lambdaScalingRadius(const Real* pos) {
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 return 1.0;
	 break;
       case sep::CARTESIAN:
	 return 1.0;
	 break;
       case sep::CYLINDRICAL:
	 return pos[0]/simControl.R_reference;
	 break;
       case sep::SPHERICAL:
	 return pos[0]/simControl.R_reference;
	 break;
       default:
	 return 1.0;
	 break;
      }
   }
   
   Real lambdaScalingRadius2(const Real* pos) {
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 return 1.0;
	 break;
       case sep::CARTESIAN:
	 return 1.0;
	 break;
       case sep::CYLINDRICAL:
	 return (pos[0]*pos[0])/(simControl.R_reference*simControl.R_reference);
	 break;
       case sep::SPHERICAL:
	 return (pos[0]*pos[0])/(simControl.R_reference*simControl.R_reference);
	 break;
       default:
	 return 1.0;
	 break;
      }
   }

} // namespace sep
