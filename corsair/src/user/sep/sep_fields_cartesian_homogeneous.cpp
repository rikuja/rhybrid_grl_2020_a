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
#include "sep_fields_cartesian_homogeneous.h"

using namespace std;

namespace sep {

   static CartesianHomogeneous cartesianHomogeneous;
   static const string PREFIX = "CartesianHomogeneous";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   extern sep::SimControl simControl;
   
   CartesianHomogeneous::CartesianHomogeneous() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      rho_mass = NAN;
      for (int i=0; i<3; ++i) B[i] = NAN;
      for (int i=0; i<3; ++i) E[i] = NAN;
      for (int i=0; i<3; ++i) V_plasma[i] = NAN;
   }
   /*
   CartesianHomogeneous::~CartesianHomogeneous() {
      finalize();
   }*/
   
   static bool addConfigFileOptions(ConfigReader& cr) {
      cr.add(PREFIX+".B_magnitude","Magnitude of magnetic field (float)",DEF_VALUE);
      cr.add(PREFIX+".B_x","Magnetic field direction, x-component (float).",DEF_VALUE);
      cr.add(PREFIX+".B_y","Magnetic field direction, y-component (float).",DEF_VALUE);
      cr.add(PREFIX+".B_z","Magnetic field direction, z-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_magnitude","Plasma speed (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_x","Plasma velocity direction, x-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_y","Plasma velocity direction, y-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_z","Plasma velocity direction, z-component (float).",DEF_VALUE);
      cr.add(PREFIX+".ion_number_density","Mass density (float).",DEF_VALUE);
      return true;
   }
   
   bool cartesianHomogeneousFieldFinalize() {return true;}

   void cartesianHomogeneousFieldGetFields(pargrid::CellID blockLID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]) {
      for (int i=0; i<3; ++i) E[i] = cartesianHomogeneous.E[i];
      for (int i=0; i<3; ++i) B[i] = cartesianHomogeneous.B[i];
      for (int i=0; i<9; ++i) gradB[i] = 0.0;
   }

   void cartesianHomogeneousFieldGetPlasmaState(pargrid::CellID blockLID,Real t,const Real* position,PlasmaState& plasmaState) {
      for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianHomogeneous.V_plasma[i];
      plasmaState.ionMassDensity = cartesianHomogeneous.rho_mass;
   }

   void cartesianHomogeneousFieldGetState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,
					  Real* RESTRICT B,Real* RESTRICT V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign) {
      for (int i=0; i<3; ++i) B[i] = cartesianHomogeneous.B[i];
      
      const Real B_mag = vectorMagnitude<3>(cartesianHomogeneous.B);
      V_alfven = alfvenSign*B_mag / sqrt(constants::PERMEABILITY * cartesianHomogeneous.rho_mass);
      //const Real parallelPlasmaSpeed = dotProduct<3>(cartesianHomogeneous.V_plasma,cartesianHomogeneous.B) / B_mag;
      
      //V_alfven  = parallelPlasmaSpeed + alfvenSign * alfvenSpeed;
      V_wave[0] = cartesianHomogeneous.V_plasma[0] + V_alfven*B[0]/B_mag;
      V_wave[1] = cartesianHomogeneous.V_plasma[1] + V_alfven*B[1]/B_mag;
      V_wave[2] = cartesianHomogeneous.V_plasma[2] + V_alfven*B[2]/B_mag;
      dV_wave = 0.0;
   }

   bool cartesianHomogeneousFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Prevent duplicate initialization:
      if (cartesianHomogeneous.initialized == true) return true;
      
      cartesianHomogeneous.initialized = true;
      cartesianHomogeneous.sim = &sim;
      cartesianHomogeneous.simClasses = &simClasses;

      // Parse config file:
      if (addConfigFileOptions(cr) == false) cartesianHomogeneous.initialized = false;
      cr.parse();
      //cr.get(PREFIX+".Ex",E[0]);
      //cr.get(PREFIX+".Ey",E[1]);
      //cr.get(PREFIX+".Ez",E[2]);
      Real B_magnitude,V_magnitude;
      cr.get(PREFIX+".B_magnitude",B_magnitude);
      cr.get(PREFIX+".B_x",cartesianHomogeneous.B[0]);
      cr.get(PREFIX+".B_y",cartesianHomogeneous.B[1]);
      cr.get(PREFIX+".B_z",cartesianHomogeneous.B[2]);
      cr.get(PREFIX+".V_plasma_magnitude",V_magnitude);
      cr.get(PREFIX+".V_plasma_x",cartesianHomogeneous.V_plasma[0]);
      cr.get(PREFIX+".V_plasma_y",cartesianHomogeneous.V_plasma[1]);
      cr.get(PREFIX+".V_plasma_z",cartesianHomogeneous.V_plasma[2]);
      cr.get(PREFIX+".ion_number_density",cartesianHomogeneous.rho_mass);
      
      // Check input values:
      if (cartesianHomogeneous.rho_mass == DEF_VALUE) cartesianHomogeneous.initialized = false;
      if (B_magnitude == DEF_VALUE) cartesianHomogeneous.initialized = false;
      if (V_magnitude == DEF_VALUE) cartesianHomogeneous.initialized = false;
      for (int i=0; i<3; ++i) if (cartesianHomogeneous.B[i] == DEF_VALUE) cartesianHomogeneous.initialized = false;
      for (int i=0; i<3; ++i) if (cartesianHomogeneous.V_plasma[i] == DEF_VALUE) cartesianHomogeneous.initialized = false;
      
      if (cartesianHomogeneous.initialized == false) {
	 cartesianHomogeneous.simClasses->logger << "(SEP FIELDS CARTEASIAN) ERROR: Missing config file options" << endl << write;
      }
      
      // Calculate B and V_plasma fields:
      cartesianHomogeneous.rho_mass *= constants::MASS_PROTON;
      unitVector<3,Real>(cartesianHomogeneous.B);
      for (int i=0; i<3; ++i) cartesianHomogeneous.B[i] *= B_magnitude;
      unitVector<3,Real>(cartesianHomogeneous.V_plasma);
      for (int i=0; i<3; ++i) cartesianHomogeneous.V_plasma[i] *= V_magnitude;
      
      // Set reference value of B:
      simControl.R_reference = NAN;
      
      // Calculate E = B x V_plasma field:
      crossProduct<Real>(cartesianHomogeneous.V_plasma,cartesianHomogeneous.B,cartesianHomogeneous.E);

      simClasses.logger << "(SEP FIELDS CARTESIAN) Initialization complete, status is ";
      if (cartesianHomogeneous.initialized == true) simClasses.logger << "success" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return cartesianHomogeneous.initialized;
   }
      
} // namespace sep