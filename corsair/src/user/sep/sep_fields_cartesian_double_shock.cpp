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
#include "sep_fields_cartesian_double_shock.h"
#include "sep_mesh_logical.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

namespace sep {

   static CartesianDoubleShock cartesianDoubleShock;
   static const string PREFIX = "CartesianDoubleShock";
   static const string SHOCKPREFIX = "Shock.planar.double";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   extern sep::SimControl simControl;

   CartesianDoubleShock::CartesianDoubleShock() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      ionNumberDensity = NAN;
      ionPolytropicIndex = NAN;
      ionTemperature1 = NAN;
      massIon = NAN;
      for (int i=0; i<3; ++i) B1[i] = NAN;
      for (int i=0; i<3; ++i) E1_SIM[i] = NAN;
      for (int i=0; i<3; ++i) V1_plasma_SIM[i] = NAN;
   }
   
   bool initialize(Real B_magnitude,Real V_magnitude) {
      bool success = true;

      // ************************************************************************** //
      // ***** CALCULATE PLASMA VARIABLES IN UPSTREAM REGION OF LEADING SHOCK ***** //
      // ************************************************************************** //

      // Calculate E, B, and V_plasma in upstream region of leading shock:
      unitVector<3,Real>(cartesianDoubleShock.B1);
      for (int i=0; i<3; ++i) cartesianDoubleShock.B1[i] *= B_magnitude;

      unitVector<3,Real>(cartesianDoubleShock.V1_plasma_SIM);
      for (int i=0; i<3; ++i) cartesianDoubleShock.V1_plasma_SIM[i] *= V_magnitude;

      crossProduct(cartesianDoubleShock.B1,cartesianDoubleShock.V1_plasma_SIM,cartesianDoubleShock.E1_SIM);

      // Calculate Alfven velocity:
      Real ionMassDensity = cartesianDoubleShock.ionNumberDensity * cartesianDoubleShock.massIon;
      for (int i=0; i<3; ++i) cartesianDoubleShock.V1_alfven[i] = cartesianDoubleShock.B1[i] / sqrt(constants::PERMEABILITY * ionMassDensity);
       
      // Get leading shock gas and magnetic compression ratios.
      Real shockCentroid[3];
      for (int i=0; i<3; ++i) shockCentroid[i] = cartesianDoubleShock.leadingShockInitialCentroid[i];
      getLogicalCoordinates(cartesianDoubleShock.sim,shockCentroid,cartesianDoubleShock.leadingShockCentroid);
      for (int i=0; i<3; ++i) shockCentroid[i] = cartesianDoubleShock.leadingShockCentroid[i];
      
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 shockCentroid[0] -= 0.001;
      } else {
	 shockCentroid[0] += 0.001;
      }

      uint32_t i_block = static_cast<uint32_t>(shockCentroid[0] / block::WIDTH_X);
      uint32_t j_block = static_cast<uint32_t>(shockCentroid[1] / block::WIDTH_Y);
      uint32_t k_block = static_cast<uint32_t>(shockCentroid[2] / block::WIDTH_Z);
      
      pargrid::CellID blockGID = block::calculateGlobalIndex(*cartesianDoubleShock.sim,i_block,j_block,k_block);
      pargrid::CellID blockLID = cartesianDoubleShock.simClasses->pargrid.getLocalID(blockGID);

      // If block local ID is valid, this process can calculate compression ratios:
      Real R_gas = 1.0, R_magn = 1.0;
      if (blockLID != pargrid::INVALID_CELLID) {
	 R_gas  = simControl.shock->getGasCompressionRatio(blockLID,cartesianDoubleShock.sim->t,shockCentroid);
	 R_magn = simControl.shock->getMagneticCompressionRatio(blockLID,cartesianDoubleShock.sim->t,shockCentroid);
      }

      // Reduce compression ratios to master process which will broadcast them to everyone:
      MPI_Reduce(&R_gas,&cartesianDoubleShock.leadingShockGasCompressionRatio,1,MPI_Type<Real>(),MPI_MAX,
		 cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Reduce(&R_magn,&cartesianDoubleShock.leadingShockMagnCompressionRatio,1,MPI_Type<Real>(),MPI_MAX,
		 cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Bcast(&cartesianDoubleShock.leadingShockGasCompressionRatio,1,MPI_Type<Real>(),
		cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Bcast(&cartesianDoubleShock.leadingShockMagnCompressionRatio,1,MPI_Type<Real>(),
		cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);

      R_gas  = cartesianDoubleShock.leadingShockGasCompressionRatio;
      R_magn = cartesianDoubleShock.leadingShockMagnCompressionRatio;
      
      cartesianDoubleShock.ionPressure1 = cartesianDoubleShock.ionNumberDensity 
	                                * constants::BOLTZMANN * cartesianDoubleShock.ionTemperature1;

      // *************************************************************** //
      // ***** CALCULATE PLASMA VARIABLES IN REGION BETWEEN SHOCKS ***** //
      // *************************************************************** //
       
      // Transformation velocity from simulation frame to leading shock's normal incidence frame:
      Real V_T_SNIF[3];
      getTransformSimToLocalSNIF(cartesianDoubleShock.V1_plasma_SIM,cartesianDoubleShock.leadingShockVelocity,
				 cartesianDoubleShock.leadingShockNormal,V_T_SNIF);

      // Plasma velocity in leading shock's normal incidence frame:
      Real V_plasma_SNIF[3];
      for (int i=0; i<3; ++i) V_plasma_SNIF[i] = cartesianDoubleShock.V1_plasma_SIM[i] - V_T_SNIF[i];

      // Transformation velocity from SNIF to de Hoffmann-Teller (HT) frame:
      Real V_T_HT[3];
      getTransformLocalSNIFToLocalHT(V_plasma_SNIF,cartesianDoubleShock.B1,cartesianDoubleShock.leadingShockNormal,V_T_HT);
      
      // Plasma velocity in local HT frame:
      //Real V_plasma_HT[3];
      //for (int i=0; i<3; ++i) V_plasma_HT[i] = V_plasma_SNIF[i] - V_T_HT[i];
       
      // Calculate B in region between shocks:
      Real B_norm = dotProduct<3>(cartesianDoubleShock.B1,cartesianDoubleShock.leadingShockNormal);
      Real B_tang[3];
      for (int i=0; i<3; ++i) B_tang[i] = cartesianDoubleShock.B1[i] - B_norm*cartesianDoubleShock.leadingShockNormal[i];
      for (int i=0; i<3; ++i) cartesianDoubleShock.B2[i] = B_norm*cartesianDoubleShock.leadingShockNormal[i]
	+ cartesianDoubleShock.leadingShockMagnCompressionRatio*B_tang[i];
       
      // Calculate E and Plasma velocity in region between shocks:      
      for (int i=0; i<3; ++i) {
	 cartesianDoubleShock.V2_plasma_SIM[i] = V_plasma_SNIF[i]/cartesianDoubleShock.leadingShockGasCompressionRatio + V_T_SNIF[i]
	   + V_T_HT[i] * (1 - cartesianDoubleShock.leadingShockMagnCompressionRatio/
                                 cartesianDoubleShock.leadingShockGasCompressionRatio);
      }
      
      crossProduct(cartesianDoubleShock.B2,cartesianDoubleShock.V2_plasma_SIM,cartesianDoubleShock.E2_SIM);
      
      // Calculate Alfven velocity:
      ionMassDensity = cartesianDoubleShock.ionNumberDensity * cartesianDoubleShock.massIon
	             * cartesianDoubleShock.leadingShockGasCompressionRatio;
      for (int i=0; i<3; ++i) cartesianDoubleShock.V2_alfven[i] = cartesianDoubleShock.B2[i] / sqrt(constants::PERMEABILITY * ionMassDensity);

      // Calculate pressure jump and temperature:
      const Real gamma = cartesianDoubleShock.ionPolytropicIndex;
      Real pressureJump = ((gamma+1)*R_gas - (gamma-1)) / ((gamma+1) - (gamma-1)*R_gas);
      Real temperatureJump = pressureJump / R_gas;
      
      cartesianDoubleShock.ionTemperature2 = cartesianDoubleShock.ionTemperature1 * temperatureJump;
      cartesianDoubleShock.ionPressure2 = cartesianDoubleShock.ionNumberDensity*R_gas
	* constants::BOLTZMANN * cartesianDoubleShock.ionTemperature2;

      // ******************************************************************* //
      // ***** CALCULATE PLASMA VARIABLES DOWNSTREAM OF TRAILING SHOCK ***** //
      // ******************************************************************* //
      
      // Get trailing shock gas and magnetic compression ratios. We can calculate these 
      // using the plasma state just downstream of leading shock:
      for (int i=0; i<3; ++i) shockCentroid[i] = cartesianDoubleShock.leadingShockCentroid[i];
      
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 shockCentroid[0] += 0.002;
      } else {
	 shockCentroid[0] -= 0.002;
      }
      i_block = static_cast<uint32_t>(shockCentroid[0] / block::WIDTH_X);
      j_block = static_cast<uint32_t>(shockCentroid[1] / block::WIDTH_Y);
      k_block = static_cast<uint32_t>(shockCentroid[2] / block::WIDTH_Z);
      blockGID = block::calculateGlobalIndex(*cartesianDoubleShock.sim,i_block,j_block,k_block);
      blockLID = cartesianDoubleShock.simClasses->pargrid.getLocalID(blockGID);

      // If block local ID is valid, this process can calculate compression ratios:
      R_gas = 1.0, R_magn = 1.0;
      if (blockLID != pargrid::INVALID_CELLID) {
	 R_gas  = simControl.shock->getGasCompressionRatio(blockLID,cartesianDoubleShock.sim->t,shockCentroid);
	 R_magn = simControl.shock->getMagneticCompressionRatio(blockLID,cartesianDoubleShock.sim->t,shockCentroid);
      }

      // Reduce compression ratios to master process which will broadcast them to everyone:
      MPI_Reduce(&R_gas,&cartesianDoubleShock.trailingShockGasCompressionRatio,1,MPI_Type<Real>(),MPI_MAX,
		 cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Reduce(&R_magn,&cartesianDoubleShock.trailingShockMagnCompressionRatio,1,MPI_Type<Real>(),MPI_MAX,
		 cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Bcast(&cartesianDoubleShock.trailingShockGasCompressionRatio,1,MPI_Type<Real>(),
		cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);
      MPI_Bcast(&cartesianDoubleShock.trailingShockMagnCompressionRatio,1,MPI_Type<Real>(),
		cartesianDoubleShock.sim->MASTER_RANK,cartesianDoubleShock.sim->comm);

      R_gas  = cartesianDoubleShock.trailingShockGasCompressionRatio;
      R_magn = cartesianDoubleShock.trailingShockMagnCompressionRatio;
      cartesianDoubleShock.doubleGasCompressionRatio = cartesianDoubleShock.leadingShockGasCompressionRatio 
	                                             * cartesianDoubleShock.trailingShockGasCompressionRatio;
      cartesianDoubleShock.doubleMagnCompressionRatio = cartesianDoubleShock.leadingShockMagnCompressionRatio
	                                              * cartesianDoubleShock.trailingShockMagnCompressionRatio;
      
      // Transformation velocity from simulation frame to leading shock's normal incidence frame:
      getTransformSimToLocalSNIF(cartesianDoubleShock.V2_plasma_SIM,cartesianDoubleShock.trailingShockVelocity,
				 cartesianDoubleShock.trailingShockNormal,V_T_SNIF);
      
      // Plasma velocity in trailing shock's normal incidence frame:
      for (int i=0; i<3; ++i) V_plasma_SNIF[i] = cartesianDoubleShock.V2_plasma_SIM[i] - V_T_SNIF[i];
      
      // Transformation velocity from SNIF to de Hoffmann-Teller (HT) frame:
      getTransformLocalSNIFToLocalHT(V_plasma_SNIF,cartesianDoubleShock.B2,cartesianDoubleShock.trailingShockNormal,V_T_HT);
      
      // Calculate B in downstream region of trailing shock:
      B_norm = dotProduct<3>(cartesianDoubleShock.B2,cartesianDoubleShock.trailingShockNormal);
      for (int i=0; i<3; ++i) B_tang[i] = cartesianDoubleShock.B2[i] - B_norm*cartesianDoubleShock.trailingShockNormal[i];
      for (int i=0; i<3; ++i) cartesianDoubleShock.B3[i] = B_norm*cartesianDoubleShock.trailingShockNormal[i]
	                                                 + cartesianDoubleShock.trailingShockMagnCompressionRatio*B_tang[i];
      
      // Calculate E and Plasma velocity in region between shocks:      
      for (int i=0; i<3; ++i) {
	 cartesianDoubleShock.V3_plasma_SIM[i] = V_plasma_SNIF[i]/cartesianDoubleShock.trailingShockGasCompressionRatio + V_T_SNIF[i]
	                                       + V_T_HT[i] * (1 - cartesianDoubleShock.trailingShockMagnCompressionRatio/
							      cartesianDoubleShock.trailingShockGasCompressionRatio);
      }

      crossProduct(cartesianDoubleShock.B3,cartesianDoubleShock.V3_plasma_SIM,cartesianDoubleShock.E3_SIM);
      
      // Calculate Alfven velocity:
      ionMassDensity *= cartesianDoubleShock.trailingShockGasCompressionRatio;
      for (int i=0; i<3; ++i) cartesianDoubleShock.V3_alfven[i] = cartesianDoubleShock.B3[i] / sqrt(constants::PERMEABILITY * ionMassDensity);

      // Calculate pressure jump and temperature:
      pressureJump = ((gamma+1)*R_gas - (gamma-1)) / ((gamma+1) - (gamma-1)*R_gas);
      temperatureJump = pressureJump / R_gas;
      cartesianDoubleShock.ionTemperature3 = cartesianDoubleShock.ionTemperature2 * temperatureJump;
      cartesianDoubleShock.ionPressure3 = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.doubleGasCompressionRatio
	* constants::BOLTZMANN * cartesianDoubleShock.ionTemperature3;

      SimulationClasses* simClasses = cartesianDoubleShock.simClasses;
      simClasses->logger << "(SEP FIELDS DOUBLE SHOCK) Plasma variables after initialization:" << endl;
      simClasses->logger << "Leading shock gas compression ratio      : " << cartesianDoubleShock.leadingShockGasCompressionRatio << endl;
      simClasses->logger << "Leading shock magnetic compression ratio : " << cartesianDoubleShock.leadingShockMagnCompressionRatio << endl;
      simClasses->logger << "Trailing shock gas compression ratio     : " << cartesianDoubleShock.trailingShockGasCompressionRatio << endl;
      simClasses->logger << "Trailing shock magnetic compression ratio: " << cartesianDoubleShock.trailingShockMagnCompressionRatio << endl;
      
      simClasses->logger << "Plasma velocity V1_SIM   : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V1_plasma_SIM[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "                V2_SIM   : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V2_plasma_SIM[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "                V3_SIM   : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V3_plasma_SIM[i] << '\t';
      simClasses->logger << endl;
      
      simClasses->logger << "Alfven velocity V1_alfven: ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V1_alfven[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "                V2_alfven: ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V2_alfven[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "                V3_alfven: ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.V3_alfven[i] << '\t';
      simClasses->logger << endl;
      
      simClasses->logger << "Electric field E1_SIM    : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.E1_SIM[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "               E2_SIM    : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.E2_SIM[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "               E3_SIM    : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.E3_SIM[i] << '\t';
      simClasses->logger << endl;
      
      simClasses->logger << "Magnetic field B1        : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.B1[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "               B2        : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.B2[i] << '\t';
      simClasses->logger << endl;
      simClasses->logger << "               B3        : ";
      for (int i=0; i<3; ++i) simClasses->logger << cartesianDoubleShock.B3[i] << '\t';
      simClasses->logger << endl;

      simClasses->logger << "Pressure (ion) 1         : " << cartesianDoubleShock.ionPressure1 << endl;
      simClasses->logger << "Pressure (ion) 2         : " << cartesianDoubleShock.ionPressure2 << endl;
      simClasses->logger << "Pressure (ion) 3         : " << cartesianDoubleShock.ionPressure3 << endl;
      simClasses->logger << "Temperature (ion) 1      : " << cartesianDoubleShock.ionTemperature1 << endl;
      simClasses->logger << "Temperature (ion) 2      : " << cartesianDoubleShock.ionTemperature2 << endl;
      simClasses->logger << "Temperature (ion) 3      : " << cartesianDoubleShock.ionTemperature3 << endl;
      
      simClasses->logger << write;
      return success;
   }
   
   static bool addConfigFileOptions(ConfigReader& cr) {
      cr.add(PREFIX+".B_magnitude","Magnitude of magnetic field in T (float)",DEF_VALUE);
      cr.add(PREFIX+".B_x","Magnetic field direction, x-component (float).",DEF_VALUE);
      cr.add(PREFIX+".B_y","Magnetic field direction, y-component (float).",DEF_VALUE);
      cr.add(PREFIX+".B_z","Magnetic field direction, z-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_magnitude","Plasma speed in m/s in simulation frame (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_x","Plasma velocity direction, x-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_y","Plasma velocity direction, y-component (float).",DEF_VALUE);
      cr.add(PREFIX+".V_plasma_z","Plasma velocity direction, z-component (float).",DEF_VALUE);
      cr.add(PREFIX+".ion_number_density","Ion number density in 1/m3 (float).",DEF_VALUE);
      cr.add(PREFIX+".ion_mass","Ion mass in proton masses (float).",DEF_VALUE);
      cr.add(PREFIX+".ion_polytropic_index","Ion polytropic index (float).",DEF_VALUE);
      cr.add(PREFIX+".ion_temperature","Ion temperature in Kelvins (float).",DEF_VALUE);
      cr.add(PREFIX+".reference_distance","Distance in length units where plasma parameters are given (float).",1.0);

      cr.add(SHOCKPREFIX+".leading_shock.speed","Speed of leading shock in simulation frame (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".leading_shock.centroid_x","Leading shock centroid position at t=0, x-component (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".leading_shock.centroid_y","Leading shock centroid position at t=0, y-component (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".leading_shock.centroid_z","Leading shock centroid position at t=0, z-component (float).",DEF_VALUE);
      
      cr.add(SHOCKPREFIX+".trailing_shock.speed","Speed of trailing shock in simulation frame (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".trailing_shock.centroid_x","Trailing shock centroid position at t=0, x-component (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".trailing_shock.centroid_y","Trailing shock centroid position at t=0, y-component (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".trailing_shock.centroid_z","Trailing shock centroid position at t=0, z-component (float).",DEF_VALUE);
      cr.add(SHOCKPREFIX+".length_units","Units in which positions are given, defaults to 'm' (string).",string("m"));
      return true;
   }
   
   bool cartesianDoubleShockFieldFinalize() {return true;}

   void cartesianDoubleShockFieldGetFields(pargrid::CellID blockLID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]) {
      #ifndef NDEBUG
         if (cartesianDoubleShock.initialized == false) {
	    cerr << "(SEP FIELDS DOUBLE SHOCK) ERROR: Fields have not been initialized" << endl;
	    exit(1);
	 }
      #endif

      // Calculate shock centroids in SIM (comoving) frame:
      Real leadingCentroid[3];
      Real trailingCentroid[3];
      const Real dt = max((Real)0.0,cartesianDoubleShock.sim->t - simControl.t_setup);
      for (int i=0; i<3; ++i) {
	 leadingCentroid[i] = cartesianDoubleShock.leadingShockInitialCentroid[i]
		            + (cartesianDoubleShock.leadingShockVelocity[i] - simControl.V_frame[i])*dt;
      }
      for (int i=0; i<3; ++i) {
	 trailingCentroid[i] = leadingCentroid[i] + cartesianDoubleShock.relativeShockCentroid[i]
	                     + cartesianDoubleShock.relativeShockVelocity[i]*t;
      }

      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 // Convert centroid coordinates to logical coordinates
	 // so that they can be efficiently compared to particle coordinates:
	 getLogicalCoordinates(cartesianDoubleShock.sim,leadingCentroid,cartesianDoubleShock.leadingShockCentroid);
	 getLogicalCoordinates(cartesianDoubleShock.sim,trailingCentroid,cartesianDoubleShock.trailingShockCentroid);
	 
	 if (position[0] > cartesianDoubleShock.trailingShockCentroid[0]) {
	    for (int i=0; i<3; ++i) E[i] = cartesianDoubleShock.E3_SIM[i];
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B3[i];
	 } else if (position[0] > cartesianDoubleShock.leadingShockCentroid[0]) {
	    for (int i=0; i<3; ++i) E[i] = cartesianDoubleShock.E2_SIM[i];
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B2[i];
	 } else {
	    for (int i=0; i<3; ++i) E[i] = cartesianDoubleShock.E1_SIM[i];
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B1[i];
	 }
	 for (int i=0; i<9; ++i) gradB[i] = 0.0;
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 // Convert logical coordinates to physical ones, we need to 
	 // guard against logical coordinates being out of simulation box here:
	 uint32_t I;
	 Real x_logical = position[0];
	 if (x_logical < 0.0) {
	    I         = 0;
	    x_logical = 0.0;
	 } else if (x_logical > cartesianDoubleShock.sim->x_blocks*block::WIDTH_X) {
	    I         = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	    x_logical = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	 } else {
	    I         = static_cast<uint32_t>(position[0]);
	 }

	 // Calculate radial coordinate in SIM (comoving) frame:
	 Real R_SIM = cartesianDoubleShock.sim->x_crds_node[I]
	            + (x_logical-I)*cartesianDoubleShock.sim->dx_cell[I];

	 // Calculate radial coordinate in fixed frame:
	 Real R_FIXED = R_SIM + simControl.V_frame[0]*dt;
	 const Real R_limited = std::max(cartesianDoubleShock.sim->x_crds_node[1],R_FIXED);

	 const Real R0 = cartesianDoubleShock.sim->x_crds_node[0];
	 const Real B_radial = cartesianDoubleShock.B1[0] * (R0/R_FIXED) * (R0/R_FIXED);

	 B[0] = B_radial;
	 B[1] = 0.0;
	 B[2] = 0.0;
	 for (int i=0; i<3; ++i) E[i] = 0.0;

	 if (position[0] > cartesianDoubleShock.trailingShockCentroid[0]) {
	    
	 } else if (position[0] > cartesianDoubleShock.leadingShockCentroid[0]) {
	    
	 } else {
	    
	 }

	 gradB[matrixIndex(0,0)] = -2*B_radial/R_limited;
	 gradB[matrixIndex(0,1)] = 0.0;
	 gradB[matrixIndex(0,2)] = 0.0;
	 gradB[matrixIndex(1,0)] = 0.0;
	 gradB[matrixIndex(1,1)] = 0.0;
	 gradB[matrixIndex(1,2)] = 0.0;
	 gradB[matrixIndex(2,0)] = 0.0;
	 gradB[matrixIndex(2,1)] = 0.0;
	 gradB[matrixIndex(2,2)] = 0.0;
      }
   }

   void cartesianDoubleShockFieldGetPlasmaState(pargrid::CellID blockLID,Real t,const Real* position,PlasmaState& plasmaState) {
      #ifndef NDEBUG
         if (cartesianDoubleShock.initialized == false) {
	    cerr << "(SEP FIELDS DOUBLE SHOCK) ERROR: Fields have not been initialized" << endl;
	    exit(1);
	 }
      #endif

      // Calculate shock centroids in SIM (comoving) frame:
      Real leadingCentroid[3];
      Real trailingCentroid[3];
      const Real dt = max((Real)0.0,cartesianDoubleShock.sim->t - simControl.t_setup);
      for (int i=0; i<3; ++i) {
	 leadingCentroid[i] = cartesianDoubleShock.leadingShockInitialCentroid[i]
	   + (cartesianDoubleShock.leadingShockVelocity[i] - simControl.V_frame[i])*dt;
      }
      for (int i=0; i<3; ++i) {
	 trailingCentroid[i] = leadingCentroid[i] + cartesianDoubleShock.relativeShockCentroid[i]
	                     + cartesianDoubleShock.relativeShockVelocity[i]*t;
      }

      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 // Convert centroid coordinates to logical coordinates
	 // so that they can be efficiently compared to particle coordinates:
	 getLogicalCoordinates(cartesianDoubleShock.sim,leadingCentroid,cartesianDoubleShock.leadingShockCentroid);
	 getLogicalCoordinates(cartesianDoubleShock.sim,trailingCentroid,cartesianDoubleShock.trailingShockCentroid);
	 
	 if (position[0] > cartesianDoubleShock.trailingShockCentroid[0]) {
	    const Real gasCompressionRatio = cartesianDoubleShock.doubleGasCompressionRatio;
	    for (int i=0; i<3; ++i) plasmaState.B[i] = cartesianDoubleShock.B3[i];
	    for (int i=0; i<3; ++i) plasmaState.E[i] = cartesianDoubleShock.E3_SIM[i];
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V3_plasma_SIM[i];
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * gasCompressionRatio;
	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature3;
	 } else if (position[0] > cartesianDoubleShock.leadingShockCentroid[0]) {
	    const Real gasCompressionRatio = cartesianDoubleShock.leadingShockGasCompressionRatio;
	    for (int i=0; i<3; ++i) plasmaState.B[i] = cartesianDoubleShock.B2[i];
	    for (int i=0; i<3; ++i) plasmaState.E[i] = cartesianDoubleShock.E2_SIM[i];
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V2_plasma_SIM[i];
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * gasCompressionRatio;
	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature2;
	 } else {
	    for (int i=0; i<3; ++i) plasmaState.B[i] = cartesianDoubleShock.B1[i];
	    for (int i=0; i<3; ++i) plasmaState.E[i] = cartesianDoubleShock.E1_SIM[i];
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V1_plasma_SIM[i];
	    
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon;
	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature1;
	 }
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 // Convert logical coordinates to physical ones, we need to 
	 // guard against logical coordinates being out of simulation box here:
	 uint32_t I;
	 Real x_logical = position[0];
	 if (x_logical < 0.0) {
	    I         = 0;
	    x_logical = 0.0;
	 } else if (x_logical > cartesianDoubleShock.sim->x_blocks*block::WIDTH_X) {
	    I         = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	    x_logical = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	 } else {
	    I         = static_cast<uint32_t>(position[0]);
	 }
	 
	 // Calculate radial coordinate in SIM (comoving) frame:
	 Real R_SIM = cartesianDoubleShock.sim->x_crds_node[I] 
	            + (x_logical-I)*cartesianDoubleShock.sim->dx_cell[I];
	 
	 // Calculate radial coordinate in fixed frame:
	 Real R_FIXED = R_SIM + simControl.V_frame[0]*dt;
	 //const Real R_limited = std::max(cartesianDoubleShock.sim->x_crds_node[1],R_FIXED);

	 const Real R0 = cartesianDoubleShock.sim->x_crds_node[0];
	 const Real B_radial = cartesianDoubleShock.B1[0] * (R0/R_FIXED) * (R0/R_FIXED);

	 // Set E,B field values:
	 plasmaState.B[0] = B_radial;// * (R0/R_FIXED) * (R0/R_FIXED);
	 plasmaState.B[1] = 0.0;
	 plasmaState.B[2] = 0.0;
	 for (int i=0; i<3; ++i) plasmaState.E[i] = 0.0;

	 if (R_SIM < trailingCentroid[0]) {
	 //if (position[0] < cartesianDoubleShock.trailingShockCentroid[0]) {
	    /*
	    if (cartesianDoubleShock.sim->timestep > 1900) {
	      cerr << position[0] << '\t' << cartesianDoubleShock.trailingShockCentroid[0] << " in double ds" << endl;
	    }*/
	    
	    /*
	    // In SIM frame, relative speed between trailing shock front and downstream plasma:
	    const Real V_relative_trail = cartesianDoubleShock.leadingShockVelocity[0] - cartesianDoubleShock.V3_plasma_SIM[0];
	    
	    // Calculate time t_dot_trail when plasma at position R was at trailing shock front:
	    const Real t = max(0.0,cartesianDoubleShock.sim->t - simControl.t_setup);
	    const Real t_dot_trail = t - (trailingCentroid[0] - R) / V_relative_trail;
	    
	    // Calculate position R_dot_trail of trailing shock front at time t_dot_trail:
	    Real R_dot_trail = trailingCentroid[0] - cartesianDoubleShock.trailingShockVelocity[0]*t_dot_trail;
	    R_dot_trail = max(R0,R_dot_trail);
	    
	    // ***** //
	    const Real V_relative = cartesianDoubleShock.leadingShockVelocity[0] - cartesianDoubleShock.V2_plasma_SIM[0];
	    const Real t_dot = t - (leadingCentroid[0] - R) / V_relative;
	    Real R_dot = leadingCentroid[0] - cartesianDoubleShock.leadingShockVelocity[0]*t_dot;
	    R_dot = max(R0,R_dot);
	    Real rho0 = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * (R0/R_dot) * (R0/R_dot);
	    //rho0 *= cartesianDoubleShock.leadingShockGasCompressionRatio * (R_dot/R) * (R_dot/R);
	    // ***** //
	    
	    // Plasma density drops as (R_dot_trail/R)^2 downstream of trailing shock:
	    //plasmaState.ionMassDensity = rho0 * cartesianDoubleShock.trailingShockGasCompressionRatio * (R_dot_trail/R) * (R_dot_trail/R);
	    plasmaState.ionMassDensity = rho0;
	     */
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon
	      * cartesianDoubleShock.doubleGasCompressionRatio * (R0/R_FIXED) * (R0/R_FIXED);
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V3_plasma_SIM[i]
	                                                        - simControl.V_frame[i];
	    
	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature3;
	 } else if (R_SIM < leadingCentroid[0]) {
	 //} else if (position[0] < cartesianDoubleShock.leadingShockCentroid[0]) {
	    /*
	    // Detailed calculation gives result that plasma density is simply 
	     * unshocked plasma density multiplied by gas compression ratio.
	    // In SIM frame, relative speed between leading shock front and downstream plasma:
	    const Real V_relative = cartesianDoubleShock.leadingShockVelocity[0] - cartesianDoubleShock.V2_plasma_SIM[0];
	    
	    // Calculate how long time ago (t_dot) plasma at position R was at leading shock front:
	    const Real t = cartesianDoubleShock.sim->t - simControl.t_setup;
	    const Real t_dot = -(leadingCentroid[0] - R) / V_relative;
	    
	    // Calculate position R_dot of leading shock front at time t_dot:
	    Real R_dot = leadingCentroid[0] + cartesianDoubleShock.leadingShockVelocity[0]*t_dot;
	    
	    // Calculate plasma density at position R_dot:
	    Real rho0 = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * (R0/R_dot) * (R0/R_dot);
	    if (R_dot < R0) rho0 = 0;

	    // Plasma density drops as (R_dot/R)^2 downstream of leading shock:
	    const Real gasCompressionRatio = cartesianDoubleShock.leadingShockGasCompressionRatio;
	    plasmaState.ionMassDensity = rho0 * gasCompressionRatio * (R_dot/R) * (R_dot/R);
	    */
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon
	      * cartesianDoubleShock.leadingShockGasCompressionRatio * (R0/R_FIXED) * (R0/R_FIXED);
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V2_plasma_SIM[i]
	                                                        - simControl.V_frame[i];
	    
	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature2;
	 } else {
	    for (int i=0; i<3; ++i) plasmaState.V_plasma_SIM[i] = cartesianDoubleShock.V1_plasma_SIM[i]
	                                                        - simControl.V_frame[i];
	    plasmaState.ionMassDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon
	      * (R0/R_FIXED) * (R0/R_FIXED);

	    plasmaState.ionTemperature = cartesianDoubleShock.ionTemperature1;	    
	 }
      }

      plasmaState.ionPolytropicIndex = cartesianDoubleShock.ionPolytropicIndex;
      plasmaState.alfvenSpeed2 = vectorMagnitude2<3>(plasmaState.B) / constants::PERMEABILITY / plasmaState.ionMassDensity;
      plasmaState.soundSpeed2 = cartesianDoubleShock.ionPolytropicIndex*constants::BOLTZMANN*plasmaState.ionTemperature / cartesianDoubleShock.massIon;
   }

   void cartesianDoubleShockFieldGetState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,
					  Real* RESTRICT B,Real* RESTRICT V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign) {
      #ifndef NDEBUG
         if (cartesianDoubleShock.initialized == false) {
	    cerr << "(SEP FIELDS DOUBLE SHOCK) ERROR: Fields have not been initialized" << endl;
	    exit(1);
	 }
      #endif

      const Real dt = max((Real)0.0,t - simControl.t_setup);
      Real leadingCentroid[3];
      Real trailingCentroid[3];
      for (int i=0; i<3; ++i) {
	 leadingCentroid[i] = cartesianDoubleShock.leadingShockInitialCentroid[i]
	   + (cartesianDoubleShock.leadingShockVelocity[i] - simControl.V_frame[i])*dt;
      }
      for (int i=0; i<3; ++i) {
	 trailingCentroid[i] = leadingCentroid[i] + cartesianDoubleShock.relativeShockCentroid[i] 
	                     + cartesianDoubleShock.relativeShockVelocity[i]*t;
      }
      
      getLogicalCoordinates(cartesianDoubleShock.sim,leadingCentroid,cartesianDoubleShock.leadingShockCentroid);
      getLogicalCoordinates(cartesianDoubleShock.sim,trailingCentroid,cartesianDoubleShock.trailingShockCentroid);

      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 if (position[0] > cartesianDoubleShock.trailingShockCentroid[0]) {
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B3[i];
	    for (int i=0; i<3; ++i) V_wave[i] = cartesianDoubleShock.V3_plasma_SIM[i] + alfvenSign*cartesianDoubleShock.V3_alfven[i];
	    V_alfven = alfvenSign*fabs(cartesianDoubleShock.V3_alfven[0]);
	 } else if (position[0] > cartesianDoubleShock.leadingShockCentroid[0]) {
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B2[i];
	    for (int i=0; i<3; ++i) V_wave[i] = cartesianDoubleShock.V2_plasma_SIM[i] + alfvenSign*cartesianDoubleShock.V2_alfven[i];
	    V_alfven = alfvenSign*fabs(cartesianDoubleShock.V2_alfven[0]);
	 } else {	 
	    for (int i=0; i<3; ++i) B[i] = cartesianDoubleShock.B1[i];
	    for (int i=0; i<3; ++i) V_wave[i] = cartesianDoubleShock.V1_plasma_SIM[i] + alfvenSign*cartesianDoubleShock.V1_alfven[i];
	    V_alfven = alfvenSign*fabs(cartesianDoubleShock.V1_alfven[0]);
	 }
	 dV_wave = 0.0;
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 // Convert logical coordinates to physical ones, we need to 
	 // guard against logical coordinates being out of simulation box here:
	 uint32_t I;
	 Real x_logical = position[0];
	 if (x_logical < 0.0) {
	    I         = 0;
	    x_logical = 0.0;
	 } else if (x_logical > cartesianDoubleShock.sim->x_blocks*block::WIDTH_X) {
	    I         = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	    x_logical = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	 } else {
	    I         = static_cast<uint32_t>(position[0]);
	 }
	 
	 // Calculate SIM frame radial coordinate:
	 Real R_SIM = cartesianDoubleShock.sim->x_crds_node[I] + (x_logical-I)*cartesianDoubleShock.sim->dx_cell[I];

	 // Calculate fixed frame radial coordinate:
	 uint32_t I_max = cartesianDoubleShock.sim->x_blocks*block::WIDTH_X;
	 Real R_FIXED = R_SIM + simControl.V_frame[0]*dt;
	 Real R_limited = std::max(R_FIXED,cartesianDoubleShock.sim->x_crds_node[1]);
	 R_limited = std::min(R_limited,cartesianDoubleShock.sim->x_crds_node[I_max]);

	 const Real R0 = cartesianDoubleShock.sim->x_crds_node[0];
	 const Real B_radial = cartesianDoubleShock.B1[0] * (R0/R_FIXED) * (R0/R_FIXED);
	 B[0] = B_radial;
	 B[1] = 0;
	 B[2] = 0;

	 if (R_SIM < trailingCentroid[0]) {
	    V_alfven = alfvenSign * cartesianDoubleShock.V_alfven_surface 
	             * (R0/R_limited) / sqrt(cartesianDoubleShock.doubleGasCompressionRatio);
	    V_wave[0] = cartesianDoubleShock.V3_plasma_SIM[0] + V_alfven;
	    dV_wave = -V_alfven / R_FIXED;
	 } else if (R_SIM < leadingCentroid[0]) {
	    V_alfven = alfvenSign * cartesianDoubleShock.V_alfven_surface 
	             * (R0/R_limited) / sqrt(cartesianDoubleShock.leadingShockGasCompressionRatio);
	    V_wave[0] = cartesianDoubleShock.V2_plasma_SIM[0] + V_alfven;
	    dV_wave = -V_alfven / R_FIXED;
	 } else {
	    V_alfven = alfvenSign * cartesianDoubleShock.V_alfven_surface * (R0/R_limited);
	    V_wave[0] = cartesianDoubleShock.V1_plasma_SIM[0] + V_alfven;
	    dV_wave = -V_alfven / R_FIXED;
	 }

	 V_wave[0] -= simControl.V_frame[0];
	 V_wave[1] = 0.0;
	 V_wave[2] = 0.0;
      }      
   }

   bool cartesianDoubleShockFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Prevent duplicate initialization:
      if (cartesianDoubleShock.initialized == true) return true;
      simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) Starting initialization" << endl << write;
      
      cartesianDoubleShock.initialized = true;
      cartesianDoubleShock.sim = &sim;
      cartesianDoubleShock.simClasses = &simClasses;

      // Parse config file:
      if (addConfigFileOptions(cr) == false) {
	 cartesianDoubleShock.initialized = false;
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Failed to add config file items" << endl << write;
      }
      cr.parse();
      Real B_magnitude,V_magnitude;
      cr.get(PREFIX+".B_magnitude",B_magnitude);
      cr.get(PREFIX+".B_x",cartesianDoubleShock.B1[0]);
      cr.get(PREFIX+".B_y",cartesianDoubleShock.B1[1]);
      cr.get(PREFIX+".B_z",cartesianDoubleShock.B1[2]);
      cr.get(PREFIX+".V_plasma_magnitude",V_magnitude);
      cr.get(PREFIX+".V_plasma_x",cartesianDoubleShock.V1_plasma_SIM[0]);
      cr.get(PREFIX+".V_plasma_y",cartesianDoubleShock.V1_plasma_SIM[1]);
      cr.get(PREFIX+".V_plasma_z",cartesianDoubleShock.V1_plasma_SIM[2]);
      cr.get(PREFIX+".ion_number_density",cartesianDoubleShock.ionNumberDensity);
      cr.get(PREFIX+".ion_mass",cartesianDoubleShock.massIon);
      cr.get(PREFIX+".ion_polytropic_index",cartesianDoubleShock.ionPolytropicIndex);
      cr.get(PREFIX+".ion_temperature",cartesianDoubleShock.ionTemperature1);
      cr.get(PREFIX+".reference_distance",cartesianDoubleShock.referenceDistance);
      
      // Check input values:
      if (cartesianDoubleShock.ionNumberDensity == DEF_VALUE) cartesianDoubleShock.initialized = false;
      if (cartesianDoubleShock.massIon == DEF_VALUE) cartesianDoubleShock.initialized = false;
      if (B_magnitude == DEF_VALUE) cartesianDoubleShock.initialized = false;
      if (V_magnitude == DEF_VALUE) cartesianDoubleShock.initialized = false;
      if (cartesianDoubleShock.ionTemperature1 == DEF_VALUE || cartesianDoubleShock.ionTemperature1 < 0.0) cartesianDoubleShock.initialized = false;
      if (cartesianDoubleShock.ionPolytropicIndex == DEF_VALUE || cartesianDoubleShock.ionPolytropicIndex < 1.0) cartesianDoubleShock.initialized = false;
      for (int i=0; i<3; ++i) if (cartesianDoubleShock.B1[i] == DEF_VALUE) cartesianDoubleShock.initialized = false;
      for (int i=0; i<3; ++i) if (cartesianDoubleShock.V1_plasma_SIM[i] == DEF_VALUE) cartesianDoubleShock.initialized = false;

      // Scale ion mass to kg:
      cartesianDoubleShock.massIon *= constants::MASS_PROTON;
      
      // Upstream magnetic field must be in xz-plane:
      cartesianDoubleShock.B1[1] = 0.0;

      if (cartesianDoubleShock.initialized == false) {
	 cartesianDoubleShock.simClasses->logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Missing config file options" << endl << write;
      }
      
      // Read shock parameters:
      cr.get(SHOCKPREFIX+".leading_shock.speed",cartesianDoubleShock.leadingShockVelocity[0]);
      cr.get(SHOCKPREFIX+".leading_shock.centroid_x",cartesianDoubleShock.leadingShockInitialCentroid[0]);
      cr.get(SHOCKPREFIX+".leading_shock.centroid_y",cartesianDoubleShock.leadingShockInitialCentroid[1]);
      cr.get(SHOCKPREFIX+".leading_shock.centroid_z",cartesianDoubleShock.leadingShockInitialCentroid[2]);
      
      cr.get(SHOCKPREFIX+".trailing_shock.speed",cartesianDoubleShock.trailingShockVelocity[0]);
      cr.get(SHOCKPREFIX+".trailing_shock.centroid_x",cartesianDoubleShock.trailingShockInitialCentroid[0]);
      cr.get(SHOCKPREFIX+".trailing_shock.centroid_y",cartesianDoubleShock.trailingShockInitialCentroid[1]);
      cr.get(SHOCKPREFIX+".trailing_shock.centroid_z",cartesianDoubleShock.trailingShockInitialCentroid[2]);

      string lengthUnitsString;
      cr.get(SHOCKPREFIX+".length_units",lengthUnitsString);
      const double lengthScale = simClasses.constants.getDistanceInSI(lengthUnitsString);
      if (lengthScale == numeric_limits<double>::infinity()) {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Unknown shock length units" << endl << write;
	 cartesianDoubleShock.initialized = false;
      }
      cartesianDoubleShock.referenceDistance *= lengthScale;
      
      // Scale magnetic field and number density to 1 solar radii (spherical geometry only):
      if (simControl.coordinateSystem == sep::SPHERICAL) {
	 //Real V_alfven = sqrt(cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * constants::PERMEABILITY);
	 //V_alfven = B_magnitude / V_alfven;
	 //cerr << "V_a " << V_alfven << endl;
	 
	 Real ratio2 = cartesianDoubleShock.referenceDistance / constants::DIST_SOLAR_RADIUS;
	 ratio2 *= ratio2;
	 B_magnitude *= ratio2;
	 cartesianDoubleShock.ionNumberDensity *= ratio2;

	 //V_alfven = sqrt(cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon * constants::PERMEABILITY);
	 //V_alfven = B_magnitude / V_alfven;
	 //cerr << "V_a " << V_alfven << endl;
	 
	 // Scale values to mesh lower -x boundary:
	 ratio2 = constants::DIST_SOLAR_RADIUS / sim.x_crds_node[0];
	 ratio2 *= ratio2;
	 B_magnitude *= ratio2;
	 cartesianDoubleShock.ionNumberDensity *= ratio2;
      }

      // Shocks only propagate to +x direction:
      for (int i=1; i<3; ++i) cartesianDoubleShock.leadingShockVelocity[i] = 0.0;
      for (int i=1; i<3; ++i) cartesianDoubleShock.trailingShockVelocity[i] = 0.0;
      for (int i=0; i<3; ++i) cartesianDoubleShock.leadingShockNormal[i] = 0.0;
      for (int i=0; i<3; ++i) cartesianDoubleShock.trailingShockNormal[i] = 0.0;

      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Unknown mesh geometry" << endl << write;
	 cartesianDoubleShock.initialized = false;
	 break;
       case sep::CARTESIAN:
	 for (int i=0; i<3; ++i) cartesianDoubleShock.leadingShockInitialCentroid[i] *= lengthScale;
	 for (int i=0; i<3; ++i) cartesianDoubleShock.trailingShockInitialCentroid[i] *= lengthScale;
	 cartesianDoubleShock.leadingShockNormal[0] = -1.0;
	 cartesianDoubleShock.trailingShockNormal[0] = -1.0;
	 cartesianDoubleShock.trailingShockVelocity[0] *= -1.0;
	 break;
       case sep::CYLINDRICAL:
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: cylindrical geometry not implemented" << endl << write;
	 cartesianDoubleShock.initialized = false;
	 break;
       case sep::SPHERICAL:
	 cartesianDoubleShock.leadingShockInitialCentroid[0]  *= lengthScale;
	 cartesianDoubleShock.trailingShockInitialCentroid[0] *= lengthScale;
	 cartesianDoubleShock.leadingShockNormal[0]  = +1.0;
	 cartesianDoubleShock.trailingShockNormal[0] = +1.0;
	 break;
       default:
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Unknown mesh geometry" << endl << write;
	 cartesianDoubleShock.initialized = false;
	 break;
      }
      // Check shock parameters:
      if (cartesianDoubleShock.trailingShockVelocity[0] < 0.0) {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Shock speeds must be positive" << endl << write;
	 cartesianDoubleShock.initialized = false;
      }
      
      for (int i=0; i<3; ++i) cartesianDoubleShock.relativeShockVelocity[i]
	= cartesianDoubleShock.trailingShockVelocity[i] 
	- cartesianDoubleShock.leadingShockVelocity[i];
      
      for (int i=0; i<3; ++i) cartesianDoubleShock.relativeShockCentroid[i] 
	= cartesianDoubleShock.trailingShockInitialCentroid[i]
	- cartesianDoubleShock.leadingShockInitialCentroid[i];
      
      // Initialize gas and magnetic compression ratios to NANs. We need to get 
      // these values from shock wave class, but it has not been initialized yet. 
      // We will calculate these later.
      cartesianDoubleShock.leadingShockGasCompressionRatio = NAN;
      cartesianDoubleShock.trailingShockGasCompressionRatio = NAN;
      cartesianDoubleShock.leadingShockMagnCompressionRatio = NAN;
      cartesianDoubleShock.trailingShockMagnCompressionRatio = NAN;
      cartesianDoubleShock.doubleGasCompressionRatio = NAN;
      cartesianDoubleShock.doubleMagnCompressionRatio = NAN;

      // Set reference value of B (dummy value here):
      simControl.R_reference = NAN;

      // Attempt to initialize plasma variables including compression ratios, which 
      // need to be calculated using shock class:
      if (simControl.shock == NULL) {
	 cartesianDoubleShock.initialized = false;
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) WARNING: Could not initialize plasma variables because shock class is NULL." << endl << write;
	 cartesianDoubleShock.initialized = false;
      } else if (simControl.shock->isInitialized() == false) {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) WARNING: Could not initialize plasma variables because shock class is uninitialized." << endl << write;
	 cartesianDoubleShock.initialized = false;
      }

      // Calculate surface Alfven speed (spherical geometry only):
      if (simControl.coordinateSystem == sep::SPHERICAL) {
	 Real massDensity = cartesianDoubleShock.ionNumberDensity*cartesianDoubleShock.massIon;	 
	 cartesianDoubleShock.V_alfven_surface = fabs(B_magnitude) / sqrt(constants::PERMEABILITY * massDensity);
      }
      
      // Check that all processes have succeeded in initialization so far:
      if (simClasses.pargrid.checkSuccess(cartesianDoubleShock.initialized) == false) {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: One or more processes failed in initialization" << endl << write;
	 return false;
      }
      
      if (initialize(B_magnitude,V_magnitude) == false) {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) ERROR: Failed to initialize plasma variables." << endl << write;
	 cartesianDoubleShock.initialized = false;
	 return false;
      } else {
	 simClasses.logger << "(SEP FIELDS DOUBLE SHOCK) Initialization successful." << endl << write;
	 return true;
      }

   }
      
} // namespace sep