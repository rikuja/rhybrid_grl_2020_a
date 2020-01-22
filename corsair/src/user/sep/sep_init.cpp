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

#include <logically_cartesian_builder.h>
#include <object_factory_generic.h>

#include "sep_typedefs.h"
#include "sep_init.h"
#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_wavelength_mesh_builder.h"
#include "sep_mesh_logical.h"

// Field function includes
#include "sep_fields_cartesian_double_shock.h"
#include "sep_fields_cartesian_homogeneous.h"
#include "sep_fields_cylindrical_linecurrent.h"
#include "sep_fields_spherical_mann.h"
#include "sep_fields_spherical_parker.h"

// Particle propagators:
#include "sep_particle_scatterer.h"

// Particle includes:
#include "sep_particle_list.h"
#include "sep_particle_accumulator.h"

// Lagrangian mesh includes:
#include "sep_lagr_accumulator.h"
#include "sep_injection_buffer.h"

// Shock includes:
#include "sep_shock_planar.h"
#include "sep_shock_spherical.h"
#include "sep_shock_drift_acceleration.h"

// Includes for wrapper functions that register injectors, distributions 
// et cetera to ObjectFactories and such. These were added to speed up 
// compilation of sep_init.cpp file (lot of template functions).
#include "sep_object_wrapper.h"
#include "sep_register_bconds.h"
#include "sep_register_distributions.h"
#include "sep_register_lagr_injectors.h"
#include "sep_register_lagr_propagators.h"
#include "sep_register_particle_injectors.h"
#include "sep_register_particle_propagators.h"

using namespace std;

namespace sep {

   static sep::ObjectWrapper sepObjectWrapper;

   const int ORDER = 1;

   typedef sep::ParticleList<LAGR_SPECIES,LAGR_PARTICLE> LAGR_LIST;
   
   static Simulation* sim = NULL;
   static SimulationClasses* simClasses = NULL;
   SimControl simControl;

   InjectionBuffer<LAGR_PARTICLE> antiparInjectionBuffer;
   InjectionBuffer<LAGR_PARTICLE> parInjectionBuffer;
   map<string,InjectionBuffer<PARTICLE> > particleInjectionBuffers;

   const string writeBoolean(const bool& value) {
      if (value == true) return "YES";
      else return "NO";
   }

   bool createParticleSpecies(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const ObjectFactories& objectFactories,
			      vector<ParticleListBase*>& particleLists,const vector<string>& speciesNames) {
      bool success = true;
      
      typedef sep::ParticleList<PARTICLE_SPECIES,PARTICLE> PARTICLE_LIST;

      for (vector<string>::const_iterator it=speciesNames.begin(); it!=speciesNames.end(); ++it) {
	 particleLists.push_back(new PARTICLE_LIST);
	 if (particleLists[particleLists.size()-1]->initialize(sim,simClasses,cr,objectFactories,*it) == false) success = false;
	 ++simControl.N_particleSpecies;
      }

      return success;
   }
   
   bool earlyInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      sep::sim = &sim;
      sep::simClasses = &simClasses;
      simClasses.logger << "(SEP EARLY INIT) Starting init" << endl;

      // Add known field functions to container:
      if (sep::getObjectWrapper().fieldsContainer.registerField("CartesianHomogeneous",
								sep::cartesianHomogeneousFieldFinalize,
								sep::cartesianHomogeneousFieldGetFields,
								sep::cartesianHomogeneousFieldGetPlasmaState,
								sep::cartesianHomogeneousFieldGetState,
								sep::cartesianHomogeneousFieldInitialize) == false) success = false;
      if (sep::getObjectWrapper().fieldsContainer.registerField("CartesianDoubleShock",
								sep::cartesianDoubleShockFieldFinalize,
								sep::cartesianDoubleShockFieldGetFields,
								sep::cartesianDoubleShockFieldGetPlasmaState,
								sep::cartesianDoubleShockFieldGetState,
								sep::cartesianDoubleShockFieldInitialize) == false) success = false;
      if (sep::getObjectWrapper().fieldsContainer.registerField("CylindricalLineCurrent",
								sep::cylindricalLineCurrentFieldFinalize,
								sep::cylindricalLineCurrentFieldGetFields,
								sep::cylindricalLineCurrentFieldGetPlasmaState,
								sep::cylindricalLineCurrentFieldGetState,
								sep::cylindricalLineCurrentFieldInitialize) == false) success = false;
      if (sep::getObjectWrapper().fieldsContainer.registerField("SphericalMann",
								sep::sphericalMannFieldFinalize,
								sep::sphericalMannFieldGetFields,
								sep::sphericalMannFieldGetPlasmaState,
								sep::sphericalMannFieldGetState,
								sep::sphericalMannFieldInitialize) == false) success = false;
      if (sep::getObjectWrapper().fieldsContainer.registerField("SphericalParker",
								sep::sphericalParkerFieldFinalize,
								sep::sphericalParkerFieldGetFields,
								sep::sphericalParkerFieldGetPlasmaState,
								sep::sphericalParkerFieldGetState,
								sep::sphericalParkerFieldInitialize) == false) success = false;
      if (success == false) {
	 simClasses.logger << "\t ERROR: Failed to register one or more field functions" << endl;
      }

      // Get GridBuilder name:
      string gridBuilderName;
      cr.add("gridbuilder","Name of GridBuilder (string).",string(""));
      cr.parse();
      cr.get("gridbuilder",gridBuilderName);
      if (gridBuilderName == "") {
	 simClasses.logger << "\t ERROR: Config file does not contain GridBuilder name" << endl;
	 success = false;
      }
      
      // Get used coordinate system:
      string geometryString;
      cr.add(gridBuilderName+".geometry","Used coordinate system (string)",string(""));
      cr.parse();
      cr.get(gridBuilderName+".geometry",geometryString);
      if (geometryString == "") {
	 simClasses.logger << "\t ERROR: Config file does not contain grid geometry" << endl;
	 success = false;
      }
      
      // Get name of EM field function:
      cr.add(geometryString+".field","Name of electromagnetic field function",string(""));
      cr.parse();
      cr.get(geometryString+".field",simControl.fieldFunctionName);
      if (simControl.fieldFunctionName == "") {
	 simClasses.logger << "\t ERROR: Config file does not contain electromagnetic field function name" << endl;
	 success = false;
      }
      
      simControl.coordinateSystem = sep::UNKNOWN;
      if (geometryString == "cartesian") simControl.coordinateSystem = sep::CARTESIAN;
      else if (geometryString == "cylindrical") simControl.coordinateSystem = sep::CYLINDRICAL;
      else if (geometryString == "spherical") simControl.coordinateSystem = sep::SPHERICAL;
      
      simClasses.logger << "\t Initialization status is ";
      if (success == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      
      return success;
   }
   
   bool finalize(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;

      // Finalize electromagnetic fields:
      if (simControl.fieldsFinalize != NULL) if ((*simControl.fieldsFinalize)() == false) success = false;
      
      // Remove ParGrid arrays:
      if (simClasses.pargrid.removeUserData(simControl.particleWeightDataID) == false) success = false;
      if (simClasses.pargrid.removeUserData(simControl.antiparAlfvenWaveEnergyDataID) == false) success = false;
      if (simClasses.pargrid.removeUserData(simControl.parAlfvenWaveEnergyDataID) == false) success = false;
      if (simClasses.pargrid.removeUserData(simControl.waveGrowthDataID) == false) success = false;
      
      return success;
   }
   
   sep::ObjectWrapper& getObjectWrapper() {
      return sepObjectWrapper;
   }

   bool lateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			   const ObjectFactories& objectFactories,vector<ParticleListBase*>& particleLists) {
      bool success = true;

      // Read config file parameters:
      const string PREFIX="SEP";
      string propagateWaveString,applyWaveGrowthString,scatterParticlesString;
      string inclAntiparWavesString,inclParWavesString,shockNameString,splitterNameString;
      string limitWaveEnergyString,antiparWavelengthScalerName,parWavelengthScalerName;
      string shockUpwindingString;
      cr.add(PREFIX+".propagate_alfven_waves","If 'yes', Alfven waves are propagated (string)",string("yes"));
      cr.add(PREFIX+".alfven_wave_growth","If 'yes', Alfven wave growth is applied to make simulation self-consistent (string)",string("yes"));
      cr.add(PREFIX+".alfven_wave_growth_antiparallel","If 'no', wave growth is not applied to antiparallel waves (string)",string("yes"));
      cr.add(PREFIX+".alfven_wave_growth_parallel","If 'no', wave growth is not applied to parallel waves (string)",string("yes"));
      cr.add(PREFIX+".scatter_particles","If 'yes', particles are scattered off waves (string)",string("yes"));
      cr.add(PREFIX+".scatter_particles_antiparallel","If 'no' particles do not scatter off antiparallel waves (string)",string("yes"));
      cr.add(PREFIX+".scatter_particles_parallel","If 'no' particles do not scatter off parallel waves (string)",string("yes"));
      cr.add(PREFIX+".setup_time","Code is run until t >= t_setup until particles are injected (float)",(Real)0.0);
      cr.add(PREFIX+".shock_time","Shock is not turned on until t >= t_shock (float)",(Real)0.0);
      cr.add(PREFIX+".include_antiparallel_alfven_waves","If 'yes' antiparallel propagating Alfven waves are included (string)",string("yes"));
      cr.add(PREFIX+".include_parallel_alfven_waves","If 'yes' parallel propagating Alfven waves are included (string)",string("yes"));
      cr.add(PREFIX+".timestep_recalculate_interval","Interval (in timesteps) between maximum timestep re-evaluations, defaults to value 0 (int)",(uint32_t)0);
      cr.add(PREFIX+".timestep_maximum","Maximum allowed timestep, defaults to -1.0 (no maximum) (float).",(Real)-1.0);
      cr.add(PREFIX+".shock_name","Name of shock wave class, if left empty there is no shock (string)",string(""));
      cr.add(PREFIX+".timestep_particle_factor","Maximum particle timestep is multiplied by this factor, defaults to value 1 (float)",(Real)1.0);
      cr.add(PREFIX+".limit_wave_energy","If 'yes' a limiter is applied to available wave energy for scattering, defaults to 'no' (string)",string("no"));
      cr.add(PREFIX+".antiparallel_wavelength_scaler","Name of antiparallel wavelength mesh scaler class, 'none' (default) if no mesh is not scaled (string)",string("none"));
      cr.add(PREFIX+".parallel_wavelength_scaler","Name of antiparallel wavelength mesh scaler class, 'none' (default) if no mesh is not scaled (string)",string("none"));
      cr.add(PREFIX+".max_scattering_delta_mu","Maximum allowed change in pitch in scattering (Real)",(Real)0.1); 
      cr.add(PREFIX+".frame_speed.x","Simulation frame speed, x-component (Real)",(Real)0.0);
      cr.add(PREFIX+".frame_speed.y","Simulation frame speed, y-component (Real)",(Real)0.0);
      cr.add(PREFIX+".frame_speed.z","Simulation frame speed, z-component (Real)",(Real)0.0);
      cr.add(PREFIX+".courant_number","CFL number used to calculate simulation time step (Real)",(Real)0.5);
      cr.add(PREFIX+".lambda0","Minimum possible wavelength (Real)",(Real)1000.0);
      cr.add(PREFIX+".particle_splitter","Name of particle splitter class (string)",string(""));
      cr.add(PREFIX+".use_shock_upwinding","If 'yes' waves and particles upwind/downwind interpolations near shock (string)",string("no"));
      cr.add(PREFIX+".use_upwinding","If 'yes' accumulations and interpolations are up/downwinded everywhere (string)",string("no"));
      cr.add(PREFIX+".maximum_wave_loss_percentage","Wave packet can only lose this many percentages of its energy during one time step (float)",(Real)(0.9));
      cr.add(PREFIX+".bohm_coefficient","Proton mean free paths are limited to this coefficient times gyro radius (float)",(Real)0.05);
      cr.add(PREFIX+".shock_surface_refinements","How many times shock mesh surface elements are refined (int)",(uint32_t)3);
      cr.parse();
      cr.get(PREFIX+".propagate_alfven_waves",propagateWaveString);
      cr.get(PREFIX+".alfven_wave_growth",applyWaveGrowthString);
      cr.get(PREFIX+".scatter_particles",scatterParticlesString);
      cr.get(PREFIX+".include_antiparallel_alfven_waves",inclAntiparWavesString);
      cr.get(PREFIX+".include_parallel_alfven_waves",inclParWavesString);
      cr.get(PREFIX+".timestep_recalculate_interval",simControl.timestepRecalculateInterval);
      cr.get(PREFIX+".timestep_maximum",simControl.dt_maximum);
      cr.get(PREFIX+".shock_name",shockNameString);
      cr.get(PREFIX+".timestep_particle_factor",simControl.dt_particle_coeff);
      cr.get(PREFIX+".limit_wave_energy",limitWaveEnergyString);
      cr.get(PREFIX+".setup_time",simControl.t_setup);
      cr.get(PREFIX+".shock_time",simControl.t_shock);
      cr.get(PREFIX+".antiparallel_wavelength_scaler",antiparWavelengthScalerName);
      cr.get(PREFIX+".parallel_wavelength_scaler",parWavelengthScalerName);
      cr.get(PREFIX+".max_scattering_delta_mu",simControl.maxScatteringDeltaMu);
      cr.get(PREFIX+".frame_speed.x",simControl.V_frame[0]);
      cr.get(PREFIX+".frame_speed.y",simControl.V_frame[1]);
      cr.get(PREFIX+".frame_speed.z",simControl.V_frame[2]);
      cr.get(PREFIX+".courant_number",simControl.courantNumber);
      cr.get(PREFIX+".lambda0",simControl.lambda0);
      cr.get(PREFIX+".particle_splitter",splitterNameString);
      cr.get(PREFIX+".use_shock_upwinding",shockUpwindingString);
      cr.get(PREFIX+".maximum_wave_loss_percentage",simControl.maximumWaveEnergyLossPercentage);
      cr.get(PREFIX+".bohm_coefficient",simControl.bohmLimitCoefficient);
      cr.get(PREFIX+".shock_surface_refinements",simControl.N_shockSurfaceRefinements);

      if (propagateWaveString == "no") simControl.propagateAlfvenWaves = false;
      
      // Wave growth is turned on and by default growth factors are applied to both 
      // wave modes. Check if either (or both) were separately turned off in config file:
      if (applyWaveGrowthString == "no") simControl.applyWaveGrowth = false;
      if (simControl.applyWaveGrowth == true) {
	 simControl.applyAntiparallelWaveGrowth = true;
	 simControl.applyParallelWaveGrowth = true;
	 cr.get(PREFIX+".alfven_wave_growth_antiparallel",applyWaveGrowthString);
	 if (applyWaveGrowthString == "no") simControl.applyAntiparallelWaveGrowth = false;
	 cr.get(PREFIX+".alfven_wave_growth_parallel",applyWaveGrowthString);
	 if (applyWaveGrowthString == "no") simControl.applyParallelWaveGrowth = false;
      } else {
	 simControl.applyAntiparallelWaveGrowth = false;
	 simControl.applyParallelWaveGrowth = false;
      }

      // Scattering is turned on and by default both wave modes scatter particles.
      // Check if either (or both) were separately turned off in config file:
      if (scatterParticlesString == "no") simControl.scatterParticles = false;
      if (simControl.scatterParticles == true) {
	 simControl.scatterAntiparallel = true;
	 simControl.scatterParallel = true;
	 cr.get(PREFIX+".scatter_particles_antiparallel",scatterParticlesString);
	 if (scatterParticlesString == "no") simControl.scatterAntiparallel = false;
	 cr.get(PREFIX+".scatter_particles_parallel",scatterParticlesString);
	 if (scatterParticlesString == "no") simControl.scatterParallel = false;
      } else {
	 simControl.scatterAntiparallel = false;
	 simControl.scatterParallel = false;
      }
      
      if (inclAntiparWavesString == "no") simControl.includeAntiparWaves = false;
      if (inclParWavesString == "no") simControl.includeParWaves = false;
      if (shockUpwindingString == "yes") simControl.useShockUpwinding = true;

      // Check that maximum wave energy loss percentage has a reasonable value:
      if (simControl.maximumWaveEnergyLossPercentage < 0 || simControl.maximumWaveEnergyLossPercentage > 1) {
	 simClasses.logger << "(SEP INIT) ERROR: Value of parameter '" << PREFIX+".maximum_wave_loss_percentage' must be";
	 simClasses.logger << " between 0 and 1." << endl << write;
	 success = false;
      }
      simControl.maximumWaveEnergyLossPercentage = 1 - simControl.maximumWaveEnergyLossPercentage;

      // Check if accumulations and interpolations should be up/downwinded everywhere:
      cr.get(PREFIX+".use_upwinding",shockUpwindingString);
      if (shockUpwindingString == "yes") simControl.useUpwinding = true;
      if (simControl.useShockUpwinding == true && simControl.useUpwinding == true) {
	 simClasses.logger << "(SEP INIT) ERROR: upwinding and shock upwinding cannot currently be applied simultaneously" << endl << write;
	 success = false;
      }
      
      if (simControl.courantNumber <= 0.0 || simControl.courantNumber >= 1.0) {
	 simClasses.logger << "(SEP INIT) ERROR: Config file parameter '" << PREFIX+".courant_number' must ";
	 simClasses.logger << "have value between 0 and 1." << endl << write;
	 success = false;
      }
      
      if (simControl.maxScatteringDeltaMu >= 1.0) {
	 simClasses.logger << "(SEP INIT) ERROR: Config file parameter '" << PREFIX+".max_scattering_delta_mu' must have value ";
	 simClasses.logger << "below unity" << endl << write;
	 success = false;
      }
      
      if (simControl.dt_particle_coeff > 1.0 || simControl.dt_particle_coeff <= 0.0) {
	 simClasses.logger << "(SEP INIT) ERROR: Config file parameter '" << PREFIX+".timestep_particle_factor' has illegal value" << endl;
	 simClasses.logger << "\t It must have value between 0.0 and 1.0" << endl << write;
	 success = false;
      }
      
      simControl.order = ORDER;
      
      if (buildWaveMesh(sim,simClasses,cr) == false) success = false;

      // Init shock class if one is included in simulation:
      simControl.t_setup = max((Real)0.0,simControl.t_setup);
      simControl.shock = sep::getObjectWrapper().shockFactory.create(shockNameString);
      if (simControl.shock == NULL) {
	 simControl.includeShock = false;
      } else {
	 if (simControl.shock->initialize(sim,simClasses,cr) == false) success = false;
	 if (success == true) simControl.includeShock = true;
	 else {
	    delete simControl.shock;
	    simControl.shock = NULL;
	    simControl.includeShock = false;
	 }
      }

      // Get electromagnetic field functions and initialize them:
      if (sep::getObjectWrapper().fieldsContainer.getField(simControl.fieldFunctionName,simControl.fieldsFinalize,simControl.fieldsGetFields,
							   simControl.fieldsGetPlasmaState,simControl.fieldsGetState,simControl.fieldsInitialize) == false) {
         simClasses.logger << "(SEP LATE INIT) ERROR: Failed to get EM field functions." << endl << write;
         success = false;
      } else {
         if ((*simControl.fieldsInitialize)(sim,simClasses,cr) == false) {
            simClasses.logger << "(SEP LATE INIT) ERROR: EM field failed to initialize." << endl << write;
            success = false;
         }
      }
      
      // ************************************** //
      // ***** CREATE PARGRID DATA ARRAYS ***** //
      // ************************************** //
      
      const size_t SIZE_WAVEBLOCK = block::SIZE*simControl.N_wavelengthMeshCells;

      // Array for wave growth factors:
      simControl.waveGrowthName = "waveGrowth";
      simControl.waveGrowthDataID = 
	simClasses.pargrid.addUserData<Real>(simControl.waveGrowthName,SIZE_WAVEBLOCK);
      if (simControl.waveGrowthDataID == pargrid::INVALID_DATAID) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create array '" << simControl.waveGrowthName << "'" << endl << write;
	 success = false;
      }
      if (simClasses.pargrid.addDataTransfer(simControl.waveGrowthDataID,sim.defaultStencilID) == false) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to add transfer to array '" << simControl.waveGrowthName << "'" << endl << write;
	 success = false;
      }
      if (simClasses.pargrid.addDataTransfer(simControl.waveGrowthDataID,sim.inverseStencilID) == false) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to add transfer to array '" << simControl.waveGrowthName << "'" << endl << write;
	 success = false;
      }
      
      // Array for 4D (r,lambda) phase-space particle weight:
      simControl.particleWeightName = "particleWeight4D";
      simControl.particleWeightDataID =
	simClasses.pargrid.addUserData<Real>(simControl.particleWeightName,SIZE_WAVEBLOCK);
      if (simControl.particleWeightDataID == pargrid::INVALID_DATAID) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create array '" << simControl.particleWeightName << "'" << endl << write;
	 success = false;
      }
      if (simClasses.pargrid.addDataTransfer(simControl.particleWeightDataID,sim.inverseStencilID) == false) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to add transfer to array '" << simControl.particleWeightName << "'" << endl << write;
	 success = false;
      }
      
      simControl.temporaryWaveEnergyName = "temporaryWaveEnergy";
      simControl.temporaryWaveEnergyDataID = simClasses.pargrid.addUserData<Real>(simControl.temporaryWaveEnergyName,SIZE_WAVEBLOCK);
      if (simControl.temporaryWaveEnergyDataID == pargrid::INVALID_DATAID) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create temporary wave energy array" << endl << write;
	 success = false;
      }
      if (simClasses.pargrid.addDataTransfer(simControl.temporaryWaveEnergyDataID,sim.defaultStencilID) == false) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to transfer to temporary wave energy array" << endl << write;
	 success = false;
      }
      
      // Arrays for L,R (anti)parallel propagating Alfven waves:
      simControl.antiparAlfvenName     = "antiparAlfven";
      simControl.antiparWaveGrowthName = "antiparWaveGrowth";
      simControl.parAlfvenName         = "parAlfven";
      simControl.parWaveGrowthName     = "parWaveGrowth";
      
      // Create arrays and check for success:
      simControl.antiparAlfvenWaveEnergyDataID = 
	simClasses.pargrid.addUserData<Real>(simControl.antiparAlfvenName    ,SIZE_WAVEBLOCK);
      simControl.antiparWaveGrowthDataID =
	simClasses.pargrid.addUserData<Real>(simControl.antiparWaveGrowthName,SIZE_WAVEBLOCK);
      simControl.parAlfvenWaveEnergyDataID =
	simClasses.pargrid.addUserData<Real>(simControl.parAlfvenName        ,SIZE_WAVEBLOCK);
      simControl.parWaveGrowthDataID = 
	simClasses.pargrid.addUserData<Real>(simControl.parWaveGrowthName    ,SIZE_WAVEBLOCK);

      if (success == true) {
	 if (simControl.antiparAlfvenWaveEnergyDataID == pargrid::INVALID_DATAID) success = false;
	 if (simControl.antiparWaveGrowthDataID       == pargrid::INVALID_DATAID) success = false;
	 if (simControl.parAlfvenWaveEnergyDataID     == pargrid::INVALID_DATAID) success = false;
	 if (simControl.parWaveGrowthDataID           == pargrid::INVALID_DATAID) success = false;
	 if (success == false) {
	    simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create one or more wave mesh array(s)" << endl << write;
	 }
      }

      // Array containing min/max valid wavelength bins for scattering:
      simControl.maxAntiparValidWavelengthBinsName = "maxAntiparValidWavelengthBins";
      simControl.maxParValidWavelengthBinsName     = "maxParValidWavelengthBins";
      simControl.maxAntiparValidWavelengthBinsDataID =
	simClasses.pargrid.addUserData<sep::wmesh::validDatatype>(simControl.maxAntiparValidWavelengthBinsName,block::SIZE*2);
      simControl.maxParValidWavelengthBinsDataID =
	simClasses.pargrid.addUserData<sep::wmesh::validDatatype>(simControl.maxParValidWavelengthBinsName,block::SIZE*2);
      if (simControl.maxAntiparValidWavelengthBinsDataID == pargrid::INVALID_DATAID) {
	 success = false;
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create maxAntiparValidwavelengthBins array" << endl << write;
      }
      if (simControl.maxParValidWavelengthBinsDataID == pargrid::INVALID_DATAID) {
	 success = false;
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create maxParValidwavelengthBins array" << endl << write;
      }

      // Add transfers to arrays:
      if (success == true) {
	 if (simClasses.pargrid.addDataTransfer(simControl.antiparAlfvenWaveEnergyDataID,
						sim.defaultStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.antiparAlfvenWaveEnergyDataID,
						sim.inverseStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.antiparWaveGrowthDataID,
						sim.defaultStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.antiparWaveGrowthDataID,
						sim.inverseStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.parAlfvenWaveEnergyDataID,
						sim.defaultStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.parAlfvenWaveEnergyDataID,
						sim.inverseStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.parWaveGrowthDataID,
						sim.defaultStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.parWaveGrowthDataID,
					     sim.inverseStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.maxAntiparValidWavelengthBinsDataID,
						sim.defaultStencilID) == false) success = false;
	 if (simClasses.pargrid.addDataTransfer(simControl.maxParValidWavelengthBinsDataID,
						sim.defaultStencilID) == false) success = false;
	 if (success == false) {
	    simClasses.logger << "(SEP LATE INIT) ERROR: Failed to add transfer(s) to one or more wave array(s)."  << endl << write;
	 }
      }
      simClasses.logger << "(SEP LATE INIT) ParGrid array(s) created" << endl << write;

      // Calculate spatial cell volumes:
      sep::recalculateCellVolumes(sim,simClasses);
      
      // ********************************************* //
      // ***** CREATE PARTICLE SPECIES AND LISTS ***** //
      // ********************************************* //
      
#warning FIXME Hard-coded species type names
      simControl.lagrangianSpeciesTypename = "lagrangian";
      simControl.particleSpeciesTypename = "particle";
      
      // Get names of propagated particle species:
      vector<string> speciesNames;
      cr.addComposed("Simulation.particle.species","Names of simulated particle species (string)");
      cr.parse();
      cr.get("Simulation.particle.species",speciesNames);
      
      // Erase empty entries in vector speciesNames:
      bool erased = false;
      do {
	 erased = false;
	 for (vector<string>::iterator it=speciesNames.begin(); it!=speciesNames.end(); ++it) {
	    if ((*it).size() == 0) {
	       speciesNames.erase(it);
	       erased = true;
	       break;
	    }
	 }
      } while (erased == true);
      
      // Create particle lists:
      if (createParticleSpecies(sim,simClasses,cr,objectFactories,particleLists,speciesNames) == false) {
	 simClasses.logger << "(SEP LATE INIT) ERROR: Failed to create one or more particle species" << endl << write;
	 success = false;
      }

      // Create Lagrangian particle lists:
      if (simControl.includeAntiparWaves == true) {
	 particleLists.push_back(new LAGR_LIST);
	 if (particleLists[particleLists.size()-1]->initialize(sim,simClasses,cr,objectFactories,"AlfvenMeshAP") == false) success = false;
	 particleLists[particleLists.size()-1]->getParticles(simControl.antiparAlfvenDataID);
	 ++simControl.N_lagrangianSpecies;
      }
      if (simControl.includeParWaves == true) {
	 particleLists.push_back(new LAGR_LIST);
	 if (particleLists[particleLists.size()-1]->initialize(sim,simClasses,cr,objectFactories,"AlfvenMesh") == false) success = false;
	 particleLists[particleLists.size()-1]->getParticles(simControl.parAlfvenDataID);
	 ++simControl.N_lagrangianSpecies;
      }

      simClasses.logger << "(SEP LATE INIT) Macroparticle lists created, status is ";
      if (success == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      
      // Init shock drift acceleration class:
      if (simControl.shockAccelerator.initialize(sim,simClasses,cr) == false) {
	 simClasses.logger << "(SEP LATE INIT) Shock accelerator failed to init" << endl << write;
	 success = false;
      }

      // Check that all processes have succeeded in initialization:
      if (simClasses.pargrid.checkSuccess(success) == false) {
	 simClasses.logger << "(SEP LATE INIT) One or more processes failed to initialize" << endl << write;
	 success = false;
      }

      // Calculate value of simulation dt:
      if (success == true) recalculateTimestep(sim,simClasses,particleLists);

      // Get particle splitter:
      simControl.particleSplitter = sep::getObjectWrapper().splitterFactory.create(splitterNameString);
      if (simControl.particleSplitter == NULL && splitterNameString.size() > 0) {
	 simClasses.logger << "(SEP LATE INIT) Failed to create Particle Splitter" << endl << write;
	 success = false;
      }
      if (simControl.particleSplitter != NULL) {
	 if (simControl.particleSplitter->initialize(sim,simClasses,cr) == false) {
	    simClasses.logger << "(SEP LATE INIT) Particle splitter failed to initialize" << endl << write;
	    success = false;
	 }
      }
      
      // Write log file message:
      simClasses.logger << "(SEP LATE INIT) Initialization completed" << endl;
      simClasses.logger << "Scatter particles               : " << writeBoolean(simControl.scatterParticles) << endl;
      simClasses.logger << "Scatter particles (antiparallel): " << writeBoolean(simControl.scatterAntiparallel) << endl;
      simClasses.logger << "Scatter particles   (parallel)  : " << writeBoolean(simControl.scatterParallel) << endl;
      simClasses.logger << "Apply wave growth               : " << writeBoolean(simControl.applyWaveGrowth) << endl;
      simClasses.logger << "Apply wave growth (antiparallel): " << writeBoolean(simControl.applyAntiparallelWaveGrowth) << endl;
      simClasses.logger << "Apply wave growth   (parallel)  : " << writeBoolean(simControl.applyParallelWaveGrowth) << endl;
      simClasses.logger << write;

      // If simulation was restarted, skip creation of initial state:
      if (success == false) return success;
      if (sim.restarted == true) return success;

      // Accumulate Alfven waves and particles to get correct initial state:
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->clearAccumulationArrays() == false) success = false;
      }
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->accumulateBoundaryCells() == false) success = false;
      }
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->accumulateInnerCells() == false) success = false;
      }
      
      return success;
   }

   bool registerObjectMakers(ObjectFactories& objectFactories) {
      bool success = true;
      if (objectFactories.gridBuilders.exists("LogicallyCartesian") == false) {
	 if (objectFactories.gridBuilders.registerMaker("LogicallyCartesian",LCCreator) == false) success = false;
      }

      if (objectFactories.particleAccumulators.registerMaker("IonAccumulator",IonAccumMaker<PARTICLE_SPECIES,PARTICLE,ORDER>) == false) success = false;
      if (objectFactories.particleAccumulators.registerMaker("LagrAccumulator",LagrAccumMaker<LAGR_SPECIES,LAGR_PARTICLE,ORDER>) == false) success = false;

      if (registerBoundaryConditions() == false) success = false;
      if (registerDistributions() == false) success = false;
      if (registerParticleInjectors() == false) success = false;
      if (registerParticlePropagators() == false) success = false;
      if (registerWavePacketInjectors() == false) success = false;
      if (registerWavePacketPropagators() == false) success = false;

      if (sep::getObjectWrapper().shockFactory.registerMaker("Spherical",sep::SphericalShockMaker) == false) success = false;
      if (sep::getObjectWrapper().shockFactory.registerMaker("Planar",sep::PlanarShockMaker) == false) success = false;

      return success;
   }
   
   bool runTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists) {
      bool success = true;
      //if (sep::oblique::runTests() == false) success = false;
      if (simControl.shockAccelerator.runTests() == false) success = false;
      sep::runScattererTests(sim,simClasses);

      return success;
   }
   
} // namespace sep
