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

#include <linear_algebra.h>
#include "sep_simcontrol.h"
#include "sep_fields_cylindrical_linecurrent.h"

using namespace std;

namespace sep {

   static CylindricalLineCurrent cylindricalLineCurrent;
   static const string PREFIX = "CylindricalLineCurrent";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   extern sep::SimControl simControl;
   
   CylindricalLineCurrent::CylindricalLineCurrent() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      blockCoordinates = NULL;
      B_axial = NAN;
      B_reference = NAN;
      R_reference = NAN;
      rho_mass = NAN;
   }
   
   static bool addConfigFileOptions(ConfigReader& cr) {
      //cr.add(PREFIX+".Ex","Electric field, x-component (float).",DEF_VALUE);
      //cr.add(PREFIX+".Ey","Electric field, y-component (float).",DEF_VALUE);
      //cr.add(PREFIX+".Ez","Electric field, z-component (float).",DEF_VALUE);
      cr.add(PREFIX+".B_axial","Magnitude of axial/z magnetic field (float)",DEF_VALUE);
      cr.add(PREFIX+".B_reference","Magnitude of magnetic field (float)",DEF_VALUE);
      cr.add(PREFIX+".R_reference","Reference distance at which B is B_magnitude.",DEF_VALUE);
      cr.add(PREFIX+".mass_density","Mass density (float).",DEF_VALUE);
      return true;
   }
   
   bool cylindricalLineCurrentFieldFinalize() {return true;}

   void cylindricalLineCurrentFieldGetFields(pargrid::CellID blockLID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]) {
      #ifndef NDEBUG
      /*int32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = cylindricalLineCurrent.simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*cylindricalLineCurrent.sim,blockGID,i_block,j_block,k_block);
      bool ok = true;
      if (position[0] < i_block*block::WIDTH_X || position[0] > (i_block+1)*block::WIDTH_X) ok = false;
      if (position[1] < j_block*block::WIDTH_Y || position[1] > (j_block+1)*block::WIDTH_Y) ok = false;
      if (position[2] < k_block*block::WIDTH_Z || position[2] > (k_block+1)*block::WIDTH_Z) ok = false;
      if (ok == false) {
	 cerr << "(SEP LINECURRENT) ERROR: Invalid coordinates" << endl;
	 cerr << "\t TIME: " << cylindricalLineCurrent.sim->t << " STEP: " << cylindricalLineCurrent.sim->timestep << endl;
	 cerr << "\t BLOCK: " << i_block << ' ' << j_block << ' ' << k_block << " GID: " << blockGID << endl;
	 cerr << "\t CRDS:  " << position[0] << '\t' << position[1] << '\t' << position[2] << endl;
	 exit(1);
      }*/
      #endif
      
      // Calculate radial coordinate:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const double R = cylindricalLineCurrent.sim->x_crds_node[I] + (position[0]-I)*cylindricalLineCurrent.sim->dx_cell[I];

      for (int i=0; i<3; ++i) E[i] = 0.0;
      B[0] = 0.0;
      B[1] = cylindricalLineCurrent.B_reference * cylindricalLineCurrent.R_reference / R;
      B[2] = cylindricalLineCurrent.B_axial;
      
      gradB[matrixIndex(0,0)] = 0.0;
      gradB[matrixIndex(0,1)] = -B[1]/R;
      gradB[matrixIndex(0,2)] = 0.0;
      gradB[matrixIndex(1,0)] = -B[1]/R;
      gradB[matrixIndex(1,1)] = 0.0;
      gradB[matrixIndex(1,2)] = 0.0;
      gradB[matrixIndex(2,0)] = 0.0;
      gradB[matrixIndex(2,1)] = 0.0;
      gradB[matrixIndex(2,2)] = 0.0;
   }

   void cylindricalLineCurrentFieldGetState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,Real* B,
					    Real* RESTRICT V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign) {
      // Calculate radial coordinate:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const double R = cylindricalLineCurrent.sim->x_crds_node[I] + (position[0]-I)*cylindricalLineCurrent.sim->dx_cell[I];
      
      B[0] = 0.0;
      B[1] = cylindricalLineCurrent.B_reference * cylindricalLineCurrent.R_reference / R;
      B[2] = cylindricalLineCurrent.B_axial;
      
      //const Real B_mag = vectorMagnitude<3>(B);
      //const Real parallelPlasmaSpeed = dotProduct<3>(B,cylindricalLineCurrent.V_plasma)/B_mag;
      const Real alfvenConst = 1.0 / sqrt(constants::PERMEABILITY * cylindricalLineCurrent.rho_mass);
      V_alfven = alfvenSign*alfvenConst*vectorMagnitude<3>(B);

      V_wave[0] = cylindricalLineCurrent.V_plasma[0] + alfvenSign*alfvenConst*B[0];
      V_wave[1] = cylindricalLineCurrent.V_plasma[1] + alfvenSign*alfvenConst*B[1];
      V_wave[2] = cylindricalLineCurrent.V_plasma[2] + alfvenSign*alfvenConst*B[2];
      dV_wave = 0.0;
   }

   void cylindricalLineCurrentFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,
						  PlasmaState& plasmaState) {
      // Calculate radial coordinate:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const double R = cylindricalLineCurrent.sim->x_crds_node[I] + (position[0]-I)*cylindricalLineCurrent.sim->dx_cell[I];
      
      plasmaState.B[0] = 0.0;
      plasmaState.B[1] = cylindricalLineCurrent.B_reference * cylindricalLineCurrent.R_reference / R;
      plasmaState.B[2] = cylindricalLineCurrent.B_axial;
      
      for (int i=0; i<3; ++i) plasmaState.E[i] = 0.0;

      for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cylindricalLineCurrent.V_plasma[i];
      plasmaState.ionMassDensity = cylindricalLineCurrent.rho_mass;
   }

   bool cylindricalLineCurrentFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Prevent duplicate initialization:
      if (cylindricalLineCurrent.initialized == true) return true;
      
      cylindricalLineCurrent.initialized = true;
      cylindricalLineCurrent.sim = &sim;
      cylindricalLineCurrent.simClasses = &simClasses;

      // Parse config file:
      if (addConfigFileOptions(cr) == false) cylindricalLineCurrent.initialized = false;
      cr.parse();
      //cr.get(PREFIX+".Ex",E[0]);
      //cr.get(PREFIX+".Ey",E[1]);
      //cr.get(PREFIX+".Ez",E[2]);
      cr.get(PREFIX+".B_axial",cylindricalLineCurrent.B_axial);
      cr.get(PREFIX+".B_reference",cylindricalLineCurrent.B_reference);
      cr.get(PREFIX+".R_reference",cylindricalLineCurrent.R_reference);
      cr.get(PREFIX+".mass_density",cylindricalLineCurrent.rho_mass);

      // Get array containing block coordinates:
      cylindricalLineCurrent.blockCoordinates = getBlockCoordinateArray(sim,simClasses);
      if (cylindricalLineCurrent.blockCoordinates == NULL) {
	 simClasses.logger << "(SEP FIELDS CYLINDRICAL LC) ERROR: Block coordinate array is NULL" << endl << write;
	 cylindricalLineCurrent.initialized = false;
      }
      
      // Check input values:
      if (cylindricalLineCurrent.B_axial == DEF_VALUE) cylindricalLineCurrent.initialized = false;
      if (cylindricalLineCurrent.rho_mass == DEF_VALUE) cylindricalLineCurrent.initialized = false;
      if (cylindricalLineCurrent.B_reference == DEF_VALUE) cylindricalLineCurrent.initialized = false;
      if (cylindricalLineCurrent.R_reference == DEF_VALUE) cylindricalLineCurrent.initialized = false;
      
      if (cylindricalLineCurrent.initialized == false) {
	 cylindricalLineCurrent.simClasses->logger << "(SEP FIELDS CYLINDRICAL LC) ERROR: Missing config file options" << endl << write;
      }
      
      for (int i=0; i<3; ++i) cylindricalLineCurrent.V_plasma[i] = 0.0;

      // Set reference B:
      simControl.R_reference = cylindricalLineCurrent.R_reference;
      
      simClasses.logger << "(SEP FIELDS CYLINDRICAL LC) Initialization complete, status is ";
      if (cylindricalLineCurrent.initialized == true) simClasses.logger << "success" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return cylindricalLineCurrent.initialized;
   }
      
} // namespace sep
