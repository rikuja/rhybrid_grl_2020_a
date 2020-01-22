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
#include "sep_fields_spherical_parker.h"

using namespace std;

namespace sep {

   static SphericalParker sphericalParker;
   static const string PREFIX = "SphericalParker";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   extern sep::SimControl simControl;
   
   SphericalParker::SphericalParker() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      blockCoordinates = NULL;
      B_earth_AU2 = NAN;
      V_sw_r = NAN;
      rho_mass_AU2 = NAN;
   }
   
   static bool addConfigFileOptions(ConfigReader& cr) {
      cr.add(PREFIX+".B_reference","Magnitude of radial component of magnetic field at 1 AU (float)",DEF_VALUE);
      cr.add(PREFIX+".V_radial_sw","Radial solar wind speed at 1 AU (float).",DEF_VALUE);
      cr.add(PREFIX+".mass_density","Mass density at 1 AU (float).",DEF_VALUE);
      return true;
   }
   
   bool sphericalParkerFieldFinalize() {return true;}

   void sphericalParkerFieldGetFields(pargrid::CellID blockLID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]) {
      // Note: position is given in logical coordinates and needs to be 
      // transformed into physical spherical (r,phi,z) coordinates:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const uint32_t J = static_cast<uint32_t>(position[1]);
      const Real R = sphericalParker.sim->x_crds_node[I] + (position[0]-I)*sphericalParker.sim->dx_cell[I];
      //const double R = sphericalParker.blockCoordinates[3*blockLID+0] + position[0];
      const Real R2 = R*R;
      const Real THETA = sphericalParker.sim->y_crds_node[J] + (position[1]-J)*sphericalParker.sim->dy_cell[J];
      const Real SIN_THETA = sin(THETA);
      
      // Aliases for magnetic field components:
      const Real B_radial    = sphericalParker.B_earth_AU2 / R2;
      const Real B_azimuthal = -B_radial * R * constants::OMEGA_SOLAR_EQUATOR / sphericalParker.V_sw_r;
      
      B[0] = B_radial;
      B[1] = 0.0;
      B[2] = B_azimuthal*SIN_THETA;

      Real V_sw[3];
      V_sw[0] = sphericalParker.V_sw_r;
      V_sw[1] = 0.0;
      V_sw[2] = 0.0;
      crossProduct(B,V_sw,E);
      //for (int i=0; i<3; ++i) E[i] = 0.0;
      
      gradB[matrixIndex(0,0)] = -2*B_radial/R;
      gradB[matrixIndex(0,1)] = 0.0;
      gradB[matrixIndex(0,2)] = -B_azimuthal/R;
      gradB[matrixIndex(1,0)] = 0.0;
      gradB[matrixIndex(1,1)] = B_radial/R;
      gradB[matrixIndex(1,2)] = B_azimuthal/R/tan(THETA);
      gradB[matrixIndex(2,0)] = -B_azimuthal/R;
      gradB[matrixIndex(2,1)] = -B_azimuthal/R/tan(THETA);
      gradB[matrixIndex(2,2)] = B_radial/R;
      /*
      if (R > 1.0e9) {
	 std::cerr << sphericalParker.sim->timestep << '\t' << R << '\t' << B_radial << '\t' << B_azimuthal << std::endl;
      }*/
   }
   
   void sphericalParkerFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,
					   PlasmaState& plasmaState) {
      cerr << "sphericalParkerFieldGetPlasmaState NOT IMPLEMENTED" << endl;
      exit(1);
   }

   void sphericalParkerFieldGetState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,Real* RESTRICT B,
				     Real* RESTRICT V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign) {
      // Note: position is given in logical coordinates and needs to be 
      // transformed into physical spherical (r,phi,z) coordinates:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const Real R = sphericalParker.sim->x_crds_node[I] + (position[0]-I)*sphericalParker.sim->dx_cell[I];
      const Real R2 = R*R;
      
      // Aliases for magnetic field components:
      const Real B_radial    = sphericalParker.B_earth_AU2 / R2;
      const Real B_azimuthal = -B_radial * R * constants::OMEGA_SOLAR_EQUATOR / sphericalParker.V_sw_r;
      
      B[0] = B_radial;
      B[1] = 0.0;
      B[2] = B_azimuthal;
      
      //const Real B_mag = sqrt(B_radial*B_radial + B_azimuthal*B_azimuthal);
      
      // Calculate Alfven speed along B:
      //const Real parallelPlasmaSpeed = sphericalParker.V_sw_r * B_radial / B_mag;
      const Real rho_mass = sphericalParker.rho_mass_AU2 / R2;
      const Real alfvenConst = 1.0 / sqrt(constants::PERMEABILITY * rho_mass);
      V_alfven = alfvenSign*alfvenConst*vectorMagnitude<3>(B);
      V_wave[0] = sphericalParker.V_sw_r + alfvenSign*alfvenConst*B[0];
      V_wave[1] = alfvenSign*alfvenConst*B[1];
      V_wave[2] = alfvenSign*alfvenConst*B[2];

      #warning Derivative of alfven speed not calculated correctly
      dV_wave = 0.0;
   }
   
   /*
   void sphericalParkerFieldGetState(pargrid::CellID cellID,const Real* position,Real B[3],Real V[3],Real& rho_mass) {
      // Calculate radial coordinate:
      const uint32_t I = static_cast<uint32_t>(position[0]);
      const Real R = sphericalParker.sim->x_crds_node[I] + (position[0]-I)*sphericalParker.sim->dx_cell[I];
      const Real R2 = R*R;
      
      // Calculate B:
      const Real B_radial    = sphericalParker.B_earth_AU2 / R2;
      const Real B_azimuthal = -B_radial * R * constants::OMEGA_SOLAR_EQUATOR / sphericalParker.V_sw_r;
      
      B[0] = B_radial;
      B[1] = 0.0;
      B[2] = B_azimuthal;
      
      // Calculate V_plasma:
      V[0] = sphericalParker.V_sw_r;
      V[1] = 0.0;
      V[2] = 0.0;
      
      // Calculate mass density:
      rho_mass = sphericalParker.rho_mass_AU2 / R2;
   }
   */
   bool sphericalParkerFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Prevent duplicate initialization:
      if (sphericalParker.initialized == true) return true;
      
      sphericalParker.initialized = true;
      sphericalParker.sim = &sim;
      sphericalParker.simClasses = &simClasses;

      // Parse config file:
      if (addConfigFileOptions(cr) == false) sphericalParker.initialized = false;
      cr.parse();
      cr.get(PREFIX+".B_reference",sphericalParker.B_earth_AU2);
      cr.get(PREFIX+".V_radial_sw",sphericalParker.V_sw_r);
      cr.get(PREFIX+".mass_density",sphericalParker.rho_mass_AU2);

      // Get array containing block coordinates:
      sphericalParker.blockCoordinates = getBlockCoordinateArray(sim,simClasses);
      if (sphericalParker.blockCoordinates == NULL) {
	 simClasses.logger << "(SEP FIELDS SPHERICAL PARKER) ERROR: Block coordinate array is NULL" << endl << write;
	 sphericalParker.initialized = false;
      }
      
      // Check input values:
      if (sphericalParker.rho_mass_AU2 == DEF_VALUE) sphericalParker.initialized = false;
      if (sphericalParker.rho_mass_AU2 <= 0.0) sphericalParker.initialized = false;
      if (sphericalParker.B_earth_AU2 == DEF_VALUE) sphericalParker.initialized = false;
      if (sphericalParker.V_sw_r == DEF_VALUE) sphericalParker.initialized = false;
      
      if (sphericalParker.initialized == false) {
	 sphericalParker.simClasses->logger << "(SEP FIELDS SPHERICAL PARKER) ERROR: Missing config file options" << endl << write;
      }
      
      // Set reference value of B:
      simControl.R_reference = constants::DIST_ASTRONOMICAL_UNIT;
      
      sphericalParker.B_earth_AU2 *= constants::DIST_ASTRONOMICAL_UNIT*constants::DIST_ASTRONOMICAL_UNIT;
      sphericalParker.rho_mass_AU2 *= constants::DIST_ASTRONOMICAL_UNIT*constants::DIST_ASTRONOMICAL_UNIT;
      
      simClasses.logger << "(SEP FIELDS SPHERICAL PARKER) Initialization complete, status is ";
      if (sphericalParker.initialized == true) simClasses.logger << "success" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return sphericalParker.initialized;
   }
      
} // namespace sep
