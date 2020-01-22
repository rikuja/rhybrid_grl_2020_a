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

#include "sep_shock_paraboloid.h"
#include "sep_simcontrol.h"

using namespace std;

namespace sep {

   static const string PREFIX = "ShockParaboloid";
   extern sep::SimControl simControl;
   
   ShockParaboloid::ShockParaboloid() { }
   
   ShockParaboloid::~ShockParaboloid() { }

   bool ShockParaboloid::createNodes() {
      bool success = true;
      const Real du = maxHeight / (N_nodes_u - 1);
      const Real dv = 2.0*M_PI / N_nodes_v;
      
      for (uint32_t i=0; i<N_nodes_u; ++i) for (uint32_t j=0; j<N_nodes_v; ++j) {
	 // Write node coordinates:
	 const Real u = i*du + 0.01*maxHeight;
	 const Real v = j*dv;
	 nodeCoordinates.push_back( x0 - u ) ;
	 nodeCoordinates.push_back( radius0_y * sqrt(u/height0) * sin(v) );
	 nodeCoordinates.push_back( radius0_x * sqrt(u/height0) * cos(v) );
      }
      
      uint32_t counter = 0;
      for (uint32_t i=0; i<N_nodes_u-1; ++i) for (uint32_t j=0; j<N_nodes_v; ++j) {
	 cellConnectivity.push_back( vlsv::celltype::QUAD );
	 cellConnectivity.push_back( 4 );
	 cellConnectivity.push_back( globalIndex(i  ,j  ) );
	 
	 if (j < N_nodes_v-1) {
	    cellConnectivity.push_back( globalIndex(i  ,j+1) );
	    cellConnectivity.push_back( globalIndex(i+1,j+1) );
	 } else {
	    cellConnectivity.push_back( globalIndex(i  ,0  ) );
	    cellConnectivity.push_back( globalIndex(i+1,0  ) );
	 }
	 
	 cellConnectivity.push_back( globalIndex(i+1,j  ) );
	 
	 cellOffsets.push_back(counter*6);
	 ++counter;
      }
      cellOffsets.push_back(counter*6);
      
      return success;
   }
   
   bool ShockParaboloid::finalize() { 
      bool success = true;
      return success;
   }
   
   std::vector<uint32_t>& ShockParaboloid::getCellConnectivity() {
      return cellConnectivity;
   }
   
   std::vector<uint32_t>& ShockParaboloid::getCellOffsets() {
      return cellOffsets;
   }
   
   std::vector<Real>& ShockParaboloid::getNodeCoordinates() {
      return nodeCoordinates;
   }
   
   uint32_t ShockParaboloid::getNumberOfCells() const {
      return cellOffsets.size()-1;
   }
   
   uint32_t ShockParaboloid::getNumberOfNodes() const {
      return N_nodes_u*N_nodes_v;
   }
   
   uint32_t ShockParaboloid::globalIndex(uint32_t i,uint32_t j) {
      return i*N_nodes_v + j;
   }
   
   bool ShockParaboloid::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) { 
      bool success = true;
      
      const Real defValue = numeric_limits<Real>::infinity();
      cr.add(PREFIX+".length_units","Units in which lengths are given, defaults to 'RS' (string)",string("RS"));
      cr.add(PREFIX+".reference_height","Height at which reference radius is given (float)",defValue);
      cr.add(PREFIX+".reference_radius1","Shock radius1 at reference height (float)",defValue);
      cr.add(PREFIX+".reference_radius2","Shock radius2 at reference height (float)",defValue);
      cr.add(PREFIX+".size_radius","Number of nodes in circular direction (int)",(uint32_t)0);
      cr.add(PREFIX+".size_height","Number of nodes in height direction (int)",(uint32_t)0);
      cr.add(PREFIX+".maximum_height","Maximum height (float)",defValue);
      cr.add(PREFIX+".initial_x_position","Position of shock at t=0 (float)",defValue);
      cr.parse();
      
      string lengthUnitsString;
      cr.get(PREFIX+".length_units",lengthUnitsString);
      cr.get(PREFIX+".reference_height",height0);
      cr.get(PREFIX+".reference_radius1",radius0_x);
      cr.get(PREFIX+".reference_radius2",radius0_y);
      cr.get(PREFIX+".size_radius",N_nodes_v);
      cr.get(PREFIX+".size_height",N_nodes_u);
      cr.get(PREFIX+".maximum_height",maxHeight);
      cr.get(PREFIX+".initial_x_position",x0);
      
      // Check input values for sanity:
      const double lengthUnits = simClasses.constants.getDistanceInSI(lengthUnitsString);
      if (lengthUnits == numeric_limits<double>::infinity()) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Unsupported length units '" << lengthUnitsString << "'" << endl << write;
	 success = false;
      }
      
      if (height0 == defValue) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Reference height was not given in config file" << endl << write;
	 success = false;
      }
      if (radius0_x == defValue) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Reference radius was not given in config file" << endl << write;
	 success = false;
      }
      if (radius0_y == defValue) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Reference radius was not given in config file" << endl << write;
	 success = false;
      }
      if (maxHeight == defValue) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Maximum height was not given in config file" << endl << write;
	 success = false;
      }
      if (N_nodes_u == 0) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Number of radial nodes was not given in config file" << endl << write;
	 success = false;
      }
      if (N_nodes_v == 0) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Number of height nodes was not given in config file" << endl << write;
	 success = false;
      }
      if (x0 == defValue) {
	 simClasses.logger << "(SEP SHOCK) ERROR: Shock initial position was not given in config file" << endl << write;
	 success = false;
      }

      // Exit if error(s) have occurred:
      if (success == false) return success;
      
      // Scale distances to SI units:
      height0 *= lengthUnits;
      radius0_x *= lengthUnits;
      radius0_y *= lengthUnits;
      maxHeight *= lengthUnits;
      x0 *= lengthUnits;
      ++N_nodes_u;

      success = createNodes();
      
      return success;
   }

} // namespace sep
