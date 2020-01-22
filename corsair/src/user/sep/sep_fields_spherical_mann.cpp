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
#include <root_finding.h>
#include "sep_simcontrol.h"
#include "sep_fields_spherical_mann.h"

using namespace std;

namespace sep {

   SphericalMann sphericalMann;
   static const string PREFIX = "SphericalMann";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   extern sep::SimControl simControl;

   static Real rhs = 0.0;
   
   Real mannModel(Real x) {
      return x*x - log(x*x) - rhs;
   }
   
   void getPlasmaState(Real* pos,Real* velocity,Real& rho_mass) {
      const Real tol_diff = 1e-10;
      const Real tol_min = 1e-3;
      const Real v_min_global = 1e-8;
      const Real v_max_global = 4.0;
      
      // For now, assume spherical geometry:
      const Real radius = pos[0];

      // Below 1.8 solar radii Newkirk 1961 model is used for electron density:
      if (radius / constants::DIST_SOLAR_RADIUS <= 1.8) {
	 Real R = radius / constants::DIST_SOLAR_RADIUS;
	 Real numberDensity = (1e6 * 4.2e4)*pow(10.0,4.32/R);
	 velocity[0] = 0.0;
	 velocity[1] = 0.0;
	 velocity[2] = 0.0;
	 rho_mass = 1.92 * sphericalMann.meanMolecularWeight * constants::MASS_PROTON * numberDensity;
	 return;
      } 
      
      Real v_min = v_min_global;
      Real v_max = v_max_global;
      if (radius/sphericalMann.criticalRadius < 1.0) v_max = 1.0;
      else v_min = 1.0;
	 
      rhs = 4.0*log(radius/sphericalMann.criticalRadius) + 4*sphericalMann.criticalRadius/radius - 3.0;
      velocity[0] = rootf::bisection(mannModel,v_min,v_max,tol_diff,tol_min,1000) * sphericalMann.criticalSpeed;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
      const Real numberDensity = sphericalMann.mannModelConstant / velocity[0] / (radius*radius);
      rho_mass = 1.92 * sphericalMann.meanMolecularWeight * constants::MASS_PROTON * numberDensity;
   }

   void calculatePlasmaState() {
      delete [] sphericalMann.V_alfven;
      delete [] sphericalMann.V_radial;
      sphericalMann.N_points = sphericalMann.resolution * sphericalMann.sim->x_blocks*block::WIDTH_X + 1;
      sphericalMann.V_alfven = new Real[sphericalMann.N_points];
      sphericalMann.V_radial = new Real[sphericalMann.N_points];
      sphericalMann.r_min = sphericalMann.sim->x_crds_node[0];
      sphericalMann.r_max = sphericalMann.sim->x_crds_node[sphericalMann.sim->x_blocks*block::WIDTH_X];
      sphericalMann.dr = (sphericalMann.r_max-sphericalMann.r_min) / (sphericalMann.N_points-1);
      
      Real pos[3] = {0.0,0.0,0.0};
      Real vel[3];
      Real rho_mass;
      for (int i=0; i<sphericalMann.N_points-1; ++i) {
	 // Solve radial solar wind speed and mass density:
	 pos[0] = sphericalMann.r_min + i * sphericalMann.dr;
	 getPlasmaState(pos,vel,rho_mass);
	 sphericalMann.V_radial[i] = vel[0];
	 //sphericalMann.V_alfven[i] = rho_mass;
	 
	 // Calculate Alfven speed:
	 const Real R2 = pos[0]*pos[0];
	 const Real B_radial = sphericalMann.B_radial_R2 / R2;
	 sphericalMann.V_alfven[i] = B_radial / sqrt(constants::PERMEABILITY * rho_mass);
      }

      sphericalMann.V_radial[0] = sphericalMann.V_radial[1];
      sphericalMann.V_alfven[0] = sphericalMann.V_alfven[1];
      sphericalMann.V_radial[sphericalMann.N_points-1] = sphericalMann.V_radial[sphericalMann.N_points-2];
      sphericalMann.V_alfven[sphericalMann.N_points-1] = sphericalMann.V_alfven[sphericalMann.N_points-2];
   }

   void writeOutput() {
      if (sphericalMann.sim->mpiRank != sphericalMann.sim->MASTER_RANK) return;
      
      // Open output file on master process:
      fstream out("solar_wind_model.txt", iostream::out);
      if (out.good() == false) {
	 sphericalMann.simClasses->logger << "(SEP FIELDS SPHERICAL PARKER) ERROR: Failed to write output file" << endl << write;
	 return;
      }
      
      
      out << "#R/R_Sun\tR_AU\tV(solar wind)\tV(Alfven)\tN_e\tf_plasma(Hz)\tB(nt)" << endl;
      out << "#Speeds are in km/s" << endl;

      Real pos[3];
      Real rho_mass;
      Real V_plasma[3];
      for (int i=0; i<sphericalMann.N_points; ++i) {
	 Real R = sphericalMann.r_min + i*sphericalMann.dr;
	 pos[0] = R;
	 getPlasmaState(pos,V_plasma,rho_mass);

	 Real V_plasma = sphericalMann.V_radial[i]/1000.0;
	 Real V_alfven = sphericalMann.V_alfven[i]/1000.0;
	 const Real n_e = rho_mass/(1.92*sphericalMann.meanMolecularWeight*constants::MASS_PROTON);
	 const Real omega_pe = sqrt(n_e/constants::PERMITTIVITY/constants::MASS_ELECTRON)*constants::CHARGE_ELEMENTARY;
	 
	 out << R/constants::DIST_SOLAR_RADIUS << '\t';
	 out << R/constants::DIST_ASTRONOMICAL_UNIT << '\t';
	 out << V_plasma << '\t';
	 out << V_alfven << '\t';
	 out << n_e << '\t';
	 out << omega_pe*0.5/M_PI << '\t';
	 out << sphericalMann.B_radial_R2 / (R*R) << '\t';
	 
	 out << rho_mass << '\t';

	 out << endl;	 
      }
      
      out.close();
   }
   
   SphericalMann::SphericalMann() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      blockCoordinates = NULL;
      resolution = 1;
      V_alfven = NULL;
      V_radial = NULL;
   }
   
   SphericalMann::~SphericalMann() {
      delete [] V_alfven;
      delete [] V_radial;
   }
   
   void SphericalMann::interpolate(Real radius,Real& V_sw,Real& V_alfven,Real& dV_alfven,int alfvenSign) {
      const int i = static_cast<int>((radius - sphericalMann.r_min) / sphericalMann.dr);
      const Real R_i = sphericalMann.r_min + i*sphericalMann.dr;
      const Real coeff = (radius - R_i) / sphericalMann.dr;

      #ifndef NDEBUG
      if (i < 0 || i > sphericalMann.sim->x_crds_node[sphericalMann.sim->x_blocks*block::WIDTH_X]) {
	 cerr << "(SEP FIELDS SPHERICAL MANN) ERROR: Invalid coordinates in interpolate" << endl;
	 cerr << "\t radius: " << radius << endl;
	 cerr << "\t i     : " << i << endl;
	 cerr << "\t R_i   : " << R_i << endl;
	 exit(1);
      }
      #endif

      if (i == 0) {
	 V_sw = sphericalMann.V_radial[0];
	 V_alfven = sphericalMann.V_alfven[0];
	 dV_alfven = 0.0;
	 return;
      }
      
      V_sw     = sphericalMann.V_radial[i] + coeff * (sphericalMann.V_radial[i+1] - sphericalMann.V_radial[i]);
      V_alfven = sphericalMann.V_alfven[i] + coeff * (sphericalMann.V_alfven[i+1] - sphericalMann.V_alfven[i]);

      const Real V_n = sphericalMann.V_radial[i  ] + alfvenSign * sphericalMann.V_alfven[i  ];
      const Real V_p = sphericalMann.V_radial[i+1] + alfvenSign * sphericalMann.V_alfven[i+1];
      dV_alfven = (V_p - V_n) / sphericalMann.dr;
   }
   
   static bool addConfigFileOptions(ConfigReader& cr) {
      cr.add(PREFIX+".radial_B","Magnitude of radial component of magnetic field at reference radius (float)",DEF_VALUE);
      cr.add(PREFIX+".radial_velocity","Solar wind radial velocity at reference radius (float).",DEF_VALUE);
      cr.add(PREFIX+".number_density","Solar wind number density at reference radius (float).",DEF_VALUE);
      cr.add(PREFIX+".temperature","Plasma temperature in K (float).",DEF_VALUE);
      cr.add(PREFIX+".mean_molecular_weight","Mean molecular weight (float).",DEF_VALUE);
      cr.add(PREFIX+".reference_radius","Reference radius in given distance units, defaults to 1 (float).",(Real)1.0);
      cr.add(PREFIX+".distance_units","Units in which reference radius is given, defaults to AU.",string("AU"));
      cr.add(PREFIX+".polytropic_index","Ion polytropic index, defaults to isothermal value 1.0 (float).",(Real)1.0); 
      return true;
   }
   
   bool sphericalMannFieldFinalize() {return true;}

   void sphericalMannFieldGetFields(pargrid::CellID blockLID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]) {
      #ifndef NDEBUG
      if (sphericalMann.sim == NULL) {
	 cerr << "getFields SIM is NULL "<< endl;
	 exit(1);
      }
      #endif
      
      // Note: position is given in logical coordinates and needs to be 
      // transformed into physical spherical (r,phi,z) coordinates:
      uint32_t I;
      Real x_logical = position[0];
      if (x_logical < 0.0) {
	 I         = 0;
	 x_logical = 0.0;
      } else if (x_logical > sphericalMann.sim->x_blocks*block::WIDTH_X) {
	 I         = sphericalMann.sim->x_blocks*block::WIDTH_X;
	 x_logical = sphericalMann.sim->x_blocks*block::WIDTH_X;
      } else {
	 I         = static_cast<uint32_t>(position[0]);
      }

      #ifndef NDEBUG
      /*if (position[0] < 0.0 || position[0] > sphericalMann.sim->x_blocks*block::WIDTH_X+1 || I > sphericalMann.sim->x_blocks*block::WIDTH_X) {
	    cerr << "(SEP FIELDS SPHERICAL MANN) STEP " << sphericalMann.sim->timestep;
	    cerr << " P#" << sphericalMann.sim->mpiRank << " LID " << blockLID << " GID " << sphericalMann.simClasses->pargrid.getGlobalIDs()[blockLID];
	    cerr << " ERROR: Invalid coordinates in GetFields" << endl;
	    cerr << "\t pos: " << position[0] << '\t' << position[1] << '\t' << position[2] << endl;
	    exit(1);
	 }*/
      #endif

      const Real R = sphericalMann.sim->x_crds_node[I] + (x_logical-I)*sphericalMann.sim->dx_cell[I];
      const Real R2 = R*R;
      
      // Aliases for magnetic field components:
      const Real B_radial = sphericalMann.B_radial_R2 / R2;
      
      B[0] = B_radial;
      B[1] = 0.0;
      B[2] = 0.0;

      E[0] = 0.0;
      E[1] = 0.0;
      E[2] = 0.0;
      
      gradB[matrixIndex(0,0)] = -2*B_radial/R;
      gradB[matrixIndex(0,1)] = 0.0;
      gradB[matrixIndex(0,2)] = 0.0;
      gradB[matrixIndex(1,0)] = 0.0;
      gradB[matrixIndex(1,1)] = 0.0;
      gradB[matrixIndex(1,2)] = 0.0;
      gradB[matrixIndex(2,0)] = 0.0;
      gradB[matrixIndex(2,1)] = 0.0;
      gradB[matrixIndex(2,2)] = 0.0;
   }
   
   void sphericalMannFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* position,PlasmaState& plasmaState) {
      // Calculate radial coordinate:
      //const uint32_t I = static_cast<uint32_t>(position[0]);
      uint32_t I;
      Real x_logical = position[0];
      if (x_logical < 0.0) {
	 I         = 0;
	 x_logical = 0.0;
      } else if (x_logical > sphericalMann.sim->x_blocks*block::WIDTH_X) {
	 I         = sphericalMann.sim->x_blocks*block::WIDTH_X;
	 x_logical = sphericalMann.sim->x_blocks*block::WIDTH_X;
      } else {
	 I         = static_cast<uint32_t>(position[0]);
      }
      
      #ifndef NDEBUG
         bool ok = true;
         if (I < 0 || I > sphericalMann.sim->x_blocks*block::WIDTH_X) ok = false;
         if (position[0] != position[0] || (position[0] == -numeric_limits<Real>::infinity() || position[0] == numeric_limits<Real>::infinity())) ok = false;
         if (ok == false) {
	    stringstream ss;
	    ss << "(SEP FIELDS SPHERICAL MANN) ERROR: Invalid coordinates in GetPlasmaState" << endl;
	    ss << "\t P#" << sphericalMann.sim->mpiRank << " LID " << cellID << " GID " << sphericalMann.simClasses->pargrid.getGlobalIDs()[cellID] << endl;
	    ss << "\t pos: " << position[0] << '\t' << position[1] << '\t' << position[2] << endl;
	    
	    cerr << ss.str();
	    exit(1);
	 }
      #endif
      
      const Real R = sphericalMann.sim->x_crds_node[I] + (x_logical-I)*sphericalMann.sim->dx_cell[I];
      Real pos[3];
      pos[0] = R;
      
      getPlasmaState(pos,plasmaState.V_plasma_SIM,plasmaState.ionMassDensity);
      plasmaState.electronTemperature = sphericalMann.temperature;
      plasmaState.electronNumberDensity = plasmaState.ionMassDensity / (1.92 * sphericalMann.meanMolecularWeight * constants::MASS_PROTON);
      plasmaState.soundSpeed2 = constants::BOLTZMANN*sphericalMann.temperature/(sphericalMann.meanMolecularWeight*constants::MASS_PROTON);
      plasmaState.ionPolytropicIndex = sphericalMann.ionPolytropicIndex;
      plasmaState.ionTemperature = sphericalMann.temperature;

      const Real B_radial = sphericalMann.B_radial_R2 / (R*R);
      plasmaState.B[0] = B_radial;
      plasmaState.B[1] = 0.0;
      plasmaState.B[2] = 0.0;
      plasmaState.E[0] = 0.0;
      plasmaState.E[1] = 0.0;
      plasmaState.E[2] = 0.0;
      plasmaState.alfvenSpeed2 = B_radial*B_radial / constants::PERMEABILITY / plasmaState.ionMassDensity;
   }

   void sphericalMannFieldGetState(pargrid::CellID cellID,Real t,const Real* RESTRICT position,Real* RESTRICT B,
				     Real* RESTRICT V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign) {
      // Calculate radial coordinate:
      uint32_t I;
      Real x_logical = position[0];
      if (x_logical < 0.0) {
	 I         = 0;
	 x_logical = 0.0;
      } else if (x_logical > sphericalMann.sim->x_blocks*block::WIDTH_X) {
	 I         = sphericalMann.sim->x_blocks*block::WIDTH_X;
	 x_logical = sphericalMann.sim->x_blocks*block::WIDTH_X;
      } else {
	 I         = static_cast<uint32_t>(position[0]);
      }
      
      /* #ifndef NDEBUG
         if (position[0] < 0.0 || position[0] > sphericalMann.sim->x_blocks*block::WIDTH_X || I > sphericalMann.sim->x_blocks) {
	    cerr << "(SEP FIELDS SPHERICAL MANN) ERROR: Invalid coordinates in GetState" << endl;
	    cerr << "\t pos: " << position[0] << '\t' << position[1] << '\t' << position[2] << endl;
	    exit(1);
	 }
      #endif */

      const Real R = sphericalMann.sim->x_crds_node[I] + (x_logical-I)*sphericalMann.sim->dx_cell[I];
      const Real R2 = R*R;
      
      // Calculate B:
      const Real B_radial = sphericalMann.B_radial_R2 / R2;
      B[0] = B_radial;
      B[1] = 0.0;
      B[2] = 0.0;

      sphericalMann.interpolate(R,V_wave[0],V_alfven,dV_wave,alfvenSign);
      V_alfven *= alfvenSign;
      V_wave[0] += V_alfven;
      V_wave[1] = 0.0;
      V_wave[2] = 0.0;
   }

   bool sphericalMannFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Prevent duplicate initialization:
      if (sphericalMann.initialized == true) return true;
      
      sphericalMann.initialized = true;
      sphericalMann.sim = &sim;
      sphericalMann.simClasses = &simClasses;

      // Parse config file:
      if (addConfigFileOptions(cr) == false) sphericalMann.initialized = false;
      cr.parse();
      string referenceRadiusUnitsString;
      cr.get(PREFIX+".radial_B",sphericalMann.radialMagneticField);
      cr.get(PREFIX+".radial_velocity",sphericalMann.radialSpeed);
      cr.get(PREFIX+".number_density",sphericalMann.numberDensity);
      cr.get(PREFIX+".temperature",sphericalMann.temperature);
      cr.get(PREFIX+".mean_molecular_weight",sphericalMann.meanMolecularWeight);
      cr.get(PREFIX+".reference_radius",sphericalMann.referenceRadius);
      cr.get(PREFIX+".distance_units",referenceRadiusUnitsString);
      cr.get(PREFIX+".polytropic_index",sphericalMann.ionPolytropicIndex);
      
      // Get array containing block coordinates:
      sphericalMann.blockCoordinates = getBlockCoordinateArray(sim,simClasses);
      if (sphericalMann.blockCoordinates == NULL) {
	 simClasses.logger << "(SEP FIELDS SPHERICAL PARKER) ERROR: Block coordinate array is NULL" << endl << write;
	 sphericalMann.initialized = false;
      }
      
      // Check input values:
      if (sphericalMann.numberDensity == DEF_VALUE) sphericalMann.initialized = false;
      if (sphericalMann.numberDensity <= 0.0) sphericalMann.initialized = false;
      if (sphericalMann.radialMagneticField == DEF_VALUE) sphericalMann.initialized = false;
      if (sphericalMann.radialSpeed == DEF_VALUE) sphericalMann.initialized = false;
      if (sphericalMann.temperature == DEF_VALUE) sphericalMann.initialized = false;
      if (sphericalMann.ionPolytropicIndex < 0.0) sphericalMann.initialized = false;
      
      if (sphericalMann.initialized == false) {
	 sphericalMann.simClasses->logger << "(SEP FIELDS SPHERICAL PARKER) ERROR: Missing config file options" << endl << write;
      }

      // Set reference value of B:
      const double distanceUnits = sphericalMann.simClasses->constants.getDistanceInSI(referenceRadiusUnitsString);
      if (distanceUnits == numeric_limits<double>::max()) sphericalMann.initialized = false;
      sphericalMann.referenceRadius *= distanceUnits;
      simControl.R_reference = sphericalMann.referenceRadius;

      // Calculate derived parameters:
      sphericalMann.mannModelConstant = sphericalMann.referenceRadius*sphericalMann.referenceRadius
	                                * sphericalMann.numberDensity * sphericalMann.radialSpeed;
      sphericalMann.criticalSpeed  = sqrt(constants::BOLTZMANN*sphericalMann.temperature/sphericalMann.meanMolecularWeight/constants::MASS_PROTON);
      sphericalMann.criticalRadius = 0.5*constants::GRAVITY*constants::MASS_SOLAR/(sphericalMann.criticalSpeed*sphericalMann.criticalSpeed);

      sphericalMann.B_radial_R2 = sphericalMann.radialMagneticField * sphericalMann.referenceRadius*sphericalMann.referenceRadius;

      if (sphericalMann.initialized == true) {
	 calculatePlasmaState();
	 writeOutput();
      }
      
      simClasses.logger << "(SEP FIELDS SPHERICAL PARKER) Initialization complete, status is ";
      if (sphericalMann.initialized == true) simClasses.logger << "success" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return sphericalMann.initialized;
   }
      
} // namespace sep
