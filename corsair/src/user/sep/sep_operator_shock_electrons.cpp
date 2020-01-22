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
#include <map>
#include <vector>
#include <typeinfo>
#include <omp.h>

#include <main.h>
#include <constants.h>
#include <triangulated_sphere.h>

#include "sep_object_wrapper.h"
#include "sep_mesh_logical.h"
#include "sep_coordinate_transform.h"
#include "sep_simcontrol.h"
#include "sep_operator_shock_mesh.h"
#include "sep_operator_shock_electrons.h"
#include "sep_base_class_shock.h"
#include "sep_shock_spherical.h"
#include "sep_shock_drift_acceleration.h"
#include "sep_shock_accelerator.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   const string shockMeshName = "shock";
   const string regionName = "DataOperators.ShockElectrons";
   
   OperatorShockElectrons::OperatorShockElectrons(): DataOperator() {
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
      
      xDistribCells = 80;

      getEnergy = NULL;
      finalizeEnergy = NULL;
      energyParams = NULL;
      
      testsRun = false;
   }

   OperatorShockElectrons::~OperatorShockElectrons() { }

   bool OperatorShockElectrons::addConfigFileItems(ConfigReader& cr,const std::string& region) {
      cr.add(region+".energy_distribution","Name of electron energy distribution function (string)","");
      cr.add(region+".energy_distribution_parameters","Name of config file region containing parameters for electron energy distribution (string)","");
      cr.add(region+".injection_speed_max","Maximum value for injection speed (float)",+1.0e7);
      cr.add(region+".display_speed_min","Minimum value for speed in output data (float)",-2e7);
      cr.add(region+".display_speed_max","Maximum value for speed in output data (float)",+2e7);
      cr.add(region+".speed_samples","How many points are used to sample injection speed range (int)",(uint32_t)100);
      cr.add(region+".pitch_samples","How many points are used to sample injection pitch range (int)",(uint32_t)100);
      cr.addComposed(region+".instrument_names","Names of energy instruments (string)");
      return true;
   }
   
   bool OperatorShockElectrons::finalize() {
      bool success = true;
      
      // Finalize energy distribution:
      if (finalizeEnergy != NULL) {
	 (*finalizeEnergy)(energyParams);
      }
      finalizeEnergy = NULL;
      getEnergy = NULL;
      energyParams = NULL;

      return success;
   }

   std::string OperatorShockElectrons::getName() const {
      return "ShockMeshElectrons";
   }
   
   bool OperatorShockElectrons::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (DataOperator::initialize(cr,sim,simClasses) == false) success = false;
      
      // Initialize electron energy distribution:
      //std::string energyDistribString,energyDistribParamsString;
      if (addConfigFileItems(cr,regionName) == false) success = false;
      cr.parse();
      cr.get(regionName+".energy_distribution",energyDistribString);
      cr.get(regionName+".energy_distribution_parameters",energyDistribParamsString);
      cr.get(regionName+".injection_speed_max",V_injection_max);
      cr.get(regionName+".display_speed_min",V_par_min);
      cr.get(regionName+".display_speed_max",V_par_max);
      cr.get(regionName+".speed_samples",N_speedSamples);
      cr.get(regionName+".pitch_samples",N_pitchSamples);

      // Cell size in output distribution:
      d_V_par = (V_par_max-V_par_min) / xDistribCells;

      vector<string> instrumentNames;
      cr.get(regionName+".instrument_names",instrumentNames);
      size_t N_channels = 0;
      for (size_t i=0; i<instrumentNames.size(); ++i) {
	 // Skip empty lines
      	 if (instrumentNames[i].size() == 0) continue;
	 
	 // Read config file items for instrument and create it:
	 spacecraft::Instrument instrument;
	 if (spacecraft::createInstrument(simClasses,cr,regionName,instrumentNames[i],instrument) == false) {
	    simClasses.logger << "(SEP SHOCK ELECTRONS OP) ERROR: Failed to read required config file items for ";
	    simClasses.logger << "instrument '" << instrumentNames[i] << "'" << endl << write;
	    initialized = false;
	 } else {
	    instruments.push_back(instrument);
	    N_channels += instrument.minValues.size();
	 }
      }
      instrumentData.resize(N_channels);
      longitudeInstrData.resize(N_channels);
      filesCreated.resize(N_channels);
      for (size_t i=0; i<filesCreated.size(); ++i) filesCreated[i] = false;
      
      return success;
   }

   int32_t OperatorShockElectrons::parDistribIndex(const int32_t& i,const int32_t& j) {
      return j*xDistribCells + i;
   }

   bool OperatorShockElectrons::runTests() {
      if (testsRun == true) return testsRun;
      if (sim->mpiRank != sim->MASTER_RANK) return true;

      testsRun = true;
      return true;
   }
   
   bool OperatorShockElectrons::solveReflection(std::vector<Real>& distributionThermal,
						std::vector<Real>& distributionReflected,
						std::vector<Real>& distributionTransmitted,
						std::vector<Real>& reflectionEfficiency) {
      bool success = true;

      // Get shock surface elements on this process and solve distribution on those blocks:
      std::vector<LocalSurface> localSurfaces;
      if (simControl.shock->getLocalSurfaces(sim->t,localSurfaces,N_shockMeshCells,simControl.N_shockSurfaceRefinements) == false) {
	 success = false;
      } else {
	 for (size_t i=0; i<instrumentData.size(); ++i) {
	    instrumentData[i].resize(N_shockMeshCells);
	    for (size_t j=0; j<instrumentData[i].size(); ++j) instrumentData[i][j] = 0;
	 }
	 for (size_t i=0; i<longitudeInstrData.size(); ++i) {
	    longitudeInstrData[i].resize(yDistribCells);
	    for (size_t j=0; j<longitudeInstrData[i].size(); ++j) longitudeInstrData[i][j] = 0;
	 }
	 surfaceAreaSum.resize(yDistribCells);
	 for (size_t i=0; i<surfaceAreaSum.size(); ++i) surfaceAreaSum[i] = 0;

	 // TEST
	 const size_t SIZE = xDistribCells*yDistribCells;
	 vector<vector<Real> > cachedEfficiency(omp_get_max_threads());
	 vector<vector<Real> > cachedReflected(omp_get_max_threads());
	 vector<vector<Real> > cachedThermal(omp_get_max_threads());
	 for (size_t i=0; i<cachedThermal.size(); ++i) {
	    cachedEfficiency[i].resize(SIZE);
	    cachedReflected[i].resize(SIZE);
	    cachedThermal[i].resize(SIZE);
	 }
	 
         #pragma omp parallel
	   {	      
	      const int tid = omp_get_thread_num();
	      #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) cachedEfficiency[tid][i] = 0.0;
	      #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) cachedReflected[tid][i] = 0.0;
	      #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) cachedThermal[tid][i] = 0.0;
	   }
	 // END TEST
	 
	 // Iterate over all surface elements on this process 
	 // and solve electron reflection:
	 for (size_t s=0; s<localSurfaces.size(); ++s) {
	    // Index of the phi-cell where surface element resides (based on centroid position):
	    int32_t phiIndex = static_cast<int32_t>(localSurfaces[s].position[2]);

	    // Get plasma state:
	    PlasmaState plasmaState;
	    (*simControl.fieldsGetPlasmaState)(localSurfaces[s].localID,sim->t,localSurfaces[s].position,plasmaState);

	    // Get local shock velocity in simulation frame:
	    Real V_shock_SIM[3];
	    simControl.shock->getLocalShockVelocity(sim->t,localSurfaces[s].position,V_shock_SIM);

	    // Get shock normal:
	    Real shockNormal[3];
	    simControl.shock->getShockNormal(sim->t,localSurfaces[s].position,shockNormal);
	    if (shockNormal[0] < 0) continue;

	    // Get gas and magnetic compression ratios:
	    Real R_gas  = simControl.shock->getGasCompressionRatio(localSurfaces[s].localID,sim->t,localSurfaces[s].position);
	    Real R_magn = simControl.shock->getMagneticCompressionRatio(localSurfaces[s].localID,sim->t,localSurfaces[s].position);
	    if (fabs(R_gas) < 1.001) continue;

	    // Add surface area to area sum:
	    surfaceAreaSum[phiIndex] += localSurfaces[s].area;

	    // Measure block computation time if we're testing for repartitioning:
	    Real t_propag = 0.0;
	    if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	    
	    cerr << "surface " << s << " / " << localSurfaces.size() << endl;
	    if (solveReflection(localSurfaces[s].position,phiIndex,localSurfaces[s].zoneIndex,plasmaState.electronNumberDensity,
				plasmaState.electronTemperature,R_gas,R_magn,V_shock_SIM,shockNormal,
				plasmaState.V_plasma_SIM,plasmaState.B,
				cachedThermal,cachedReflected,distributionTransmitted,cachedEfficiency) == false) success = false;
				//distributionThermal,distributionReflected,distributionTransmitted,reflectionEfficiency) == false) success = false;

	    // Measure block computation time if we're testing for repartitioning:
	    if (sim->countPropagTime == true) {
	       t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	       simClasses->pargrid.getCellWeights()[localSurfaces[s].localID] += t_propag;
	    }
	 }
	 
	 #pragma omp parallel 
	   {
	      #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) {
		 for (int thread=0; thread<omp_get_max_threads(); ++thread) {
		    reflectionEfficiency[i] += cachedEfficiency[thread][i];
		 }
	      }
              #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) {
		 for (int thread=0; thread<omp_get_max_threads(); ++thread) {
		    distributionReflected[i] += cachedReflected[thread][i];
		 }
	      }
              #pragma omp for nowait
	      for (size_t i=0; i<SIZE; ++i) {
		 for (int thread=0; thread<omp_get_max_threads(); ++thread) {
		    distributionThermal[i] += cachedThermal[thread][i];
		 }
	      }
	   }

	 size_t counter = 0;
	 for (size_t i=0; i<instruments.size(); ++i) {
	    for (size_t c=0; c<instruments[i].minValues.size(); ++c) {
	       for (size_t s=0; s<localSurfaces.size(); ++s) {
		  int32_t phiIndex = static_cast<int32_t>(localSurfaces[s].position[2]);
		  longitudeInstrData[counter][phiIndex] += instrumentData[counter][localSurfaces[s].zoneIndex]
		                                         * localSurfaces[s].area;
	       }
	       ++counter;
	    }
	 }

	 for (size_t i=0; i<longitudeInstrData.size(); ++i) {
	    for (size_t j=0; j<longitudeInstrData[i].size(); ++j) {
	       longitudeInstrData[i][j] /= (1e-10 + surfaceAreaSum[j]);
	    }
	 }
      }

      return success;
   }
   
   bool OperatorShockElectrons::solveReflection(Real* centroid,int32_t phiIndex,uint32_t zoneIndex,Real numberDensity,Real temperature,Real R_gas,Real R_magn,
						Real* V_shock_SIM,Real* shockNormal,Real* V_plasma_SIM,Real* B,
						std::vector<std::vector<Real> >& cachedThermal,
						std::vector<std::vector<Real> >& cachedReflected,
						std::vector<Real>& distributionTransmitted,
						std::vector<std::vector<Real> >& cachedEfficiency) {
      
      Real B_mag = vectorMagnitude<3>(B);
      Real B_unit[3];
      for (int i=0; i<3; ++i) B_unit[i] = B[i];
      unitVector<3>(B_unit);
      
      // Get transform velocities:
      Real V_SIM_2_SNIF[3];
      Real V_SIM_2_HT[3];
      getTransformSimToLocalSNIF(V_plasma_SIM,V_shock_SIM,shockNormal,V_SIM_2_SNIF);
      getTransformSimToLocalHT(V_plasma_SIM,B,V_shock_SIM,shockNormal,V_SIM_2_HT);

      // Calculate electric drift velocity in SNIF:
      Real E[3];
      Real V_plasma_SNIF[3];
      for (int i=0; i<3; ++i) V_plasma_SNIF[i] = V_plasma_SIM[i] - V_SIM_2_SNIF[i];
      crossProduct(B,V_plasma_SNIF,E);
      Real V_electric_SNIF[3];
      crossProduct(E,B,V_electric_SNIF);
      for (int i=0; i<3; ++i) V_electric_SNIF[i] /= (B_mag*B_mag);

      // Setup plasma parameters for shock accelerator:
      PlasmaParameters plasmaParameters;
      plasmaParameters.B1_norm = -dotProduct<3>(B,shockNormal);
      plasmaParameters.B1_tang = sqrt(vectorMagnitude2<3>(B) - plasmaParameters.B1_norm*plasmaParameters.B1_norm);
      plasmaParameters.V1_plasma_norm = fabs(dotProduct<3>(V_plasma_SNIF,shockNormal));
      plasmaParameters.V1_plasma_tang = 0.0;
      plasmaParameters.R_gas  = R_gas;
      plasmaParameters.R_magn = R_magn;
      plasmaParameters.L_shock = 10000.0;

      Species species;
      species.mass = constants::MASS_ELECTRON;
      species.charge = -constants::CHARGE_ELEMENTARY;
      species.q_per_m = species.charge / species.mass;
      simControl.shockAccelerator.setSpecies(species);
      
      const Real mass = constants::MASS_ELECTRON;
      
      Real pitch_min = -1.0;
      Real pitch_max = +1.0;
      Real d_pitch = (pitch_max-pitch_min) / N_pitchSamples;
      
      Real U_inj_min = 0.0;
      Real U_inj_max = 0.5*constants::MASS_ELECTRON*V_injection_max*V_injection_max;
      Real dU = (U_inj_max-U_inj_min) / N_speedSamples;

      Real shapeFactors[3];
      shapeFactors[0] = 0.5*(phiIndex+1-centroid[2])*(phiIndex+1-centroid[2]);
      shapeFactors[2] = 0.5*(centroid[2]-phiIndex)*(centroid[2]-phiIndex);
      shapeFactors[1] = 1 - shapeFactors[0] - shapeFactors[2];

      int tid;
      #pragma omp parallel private(tid)
	{
	   tid = omp_get_thread_num();
	   distrib::InjectionEnergy injEnergy;
	   injEnergy.thermalEnergy = constants::BOLTZMANN*temperature;

           #pragma omp for
	   for (size_t sample=0; sample<N_speedSamples*N_pitchSamples; ++sample) {
	      const uint32_t i = sample / N_pitchSamples;
	      const uint32_t j = sample - i*N_pitchSamples;

	      // Value of energy distribution at injection energy 
	      // (getDistrib is thread-safe):
	      injEnergy.energy = U_inj_min + (i+0.5)*dU;
	      (*getDistrib)(injEnergy,energyParams);

	      // Value of (energy,pitch) distribution:
	      const Real f = 0.5*numberDensity*injEnergy.weight;

	      // Electron injection speed, pitch, gyro speed:
	      const Real speed = sqrt(2*injEnergy.energy/mass);
	      const Real pitch = pitch_min + (j+0.5)*d_pitch;
	      const Real V_gyro            = speed * sqrt(1 - pitch*pitch);

	      // Calculate incident GC velocity vector in SNIF and HT frames. 
	      // Note that these include ExB drifts:
	      Real V_incident_SNIF[3];
	      Real V_incident_HT[3];
	      for (int k=0; k<3; ++k) V_incident_SNIF[k] = speed*pitch*B_unit[k] + V_plasma_SIM[k] - V_SIM_2_SNIF[k];
	      for (int k=0; k<3; ++k) V_incident_HT[k]   = speed*pitch*B_unit[k] + V_plasma_SIM[k] - V_SIM_2_HT[k];
	 
	      ParticleParameters particle;
	      particle.state[shockaccelerator::XPOS]  = -0.501*plasmaParameters.L_shock;
	      particle.state[shockaccelerator::YPOS]  = 0.0;
	      particle.state[shockaccelerator::ZPOS]  = 0.0;
	      particle.state[shockaccelerator::V_PAR] = dotProduct<3>(V_incident_SNIF,B_unit);
	      particle.state[shockaccelerator::MU]    = 0.5*mass*V_gyro*V_gyro / B_mag;
	 
	      // Calculate pitch in HT frame:
	      const Real V_parallel_HT = dotProduct<3>(V_incident_HT,B_unit);
	      const Real V_HT = sqrt(V_parallel_HT*V_parallel_HT + V_gyro*V_gyro);
	      const Real cos_theta_HT = V_parallel_HT / V_HT;

	      // Add electron to thermal distribution:
	      Real I = (speed*pitch - V_par_min) / d_V_par;
	      int32_t I_index = static_cast<int32_t>(I);
	      
	      for (int32_t J=-1; J<2; ++J) {
		 if (phiIndex+J < 0) continue;
		 if (phiIndex+J >= yDistribCells) continue;
		 
		 Real W_L = 0.5*(I_index+1-I)*(I_index+1-I);
		 Real W_U = 0.5*(I-I_index)*(I-I_index);
		 Real W_C = 1 - W_L - W_U;
		 if ((I_index-1) >= 0 && (I_index-1) < (int64_t)xDistribCells) {
		    cachedThermal[tid][parDistribIndex(I_index-1,phiIndex+J)] += f*d_pitch*dU * W_L * shapeFactors[J+1];
		 }
		 if (I_index     >= 0 && I_index < (int64_t)xDistribCells) {
		    cachedThermal[tid][parDistribIndex(I_index  ,phiIndex+J)] += f*d_pitch*dU * W_C * shapeFactors[J+1];
		 }
		 if ((I_index+1) >= 0 && (I_index+1) < (int64_t)xDistribCells) {
		    cachedThermal[tid][parDistribIndex(I_index+1,phiIndex+J)] += f*d_pitch*dU * W_U * shapeFactors[J+1];
		 }
	      }

	      const int32_t I_inj = I_index;
	      
	      if ((plasmaParameters.B1_norm < 0 && cos_theta_HT < 0) || (plasmaParameters.B1_norm > 0 && cos_theta_HT > 0)) {
		 // Solve shock encounter:
		 simControl.shockAccelerator.solveShockEncounterGC(plasmaParameters,particle);
	    
		 if (particle.state[shockaccelerator::XPOS] < 0.0) { // Reflection
		    // Calculate reflected velocity in SNIF:
		    Real V_reflected_SNIF[3];
		    for (int k=0; k<3; ++k) V_reflected_SNIF[k] = particle.state[shockaccelerator::V_PAR]*B_unit[k] + V_electric_SNIF[k];
		    
		    // Calculate reflected velocity in plasma frame:
		    Real V_reflected_PLASMA[3];
		    for (int k=0; k<3; ++k) V_reflected_PLASMA[k] = V_reflected_SNIF[k] + V_SIM_2_SNIF[k] - V_plasma_SIM[k];
		    Real factor1 = fabs(dotProduct<3>(V_incident_SNIF,shockNormal));
		    Real factor2 = fabs(dotProduct<3>(V_reflected_PLASMA,shockNormal));
		    Real V_parallel_new_plasma = dotProduct<3>(V_reflected_PLASMA,B_unit);
	       
		    // Flux-weighted reflected number density:
		    Real amount = f*d_pitch*dU * factor1/factor2;
		    
		    // Velocity bin index as float and integer:
		    I = (V_parallel_new_plasma - V_par_min) / d_V_par;
		    I_index = static_cast<int32_t>(I);

		    for (int32_t J=-1; J<2; ++J) {
		       if (phiIndex+J < 0) continue;
		       if (phiIndex+J >= yDistribCells) continue;
		       
		       // Calculate TSC shape factors and accumulate to distribution:
		       Real W_L = 0.5*(I_index+1-I)*(I_index+1-I);
		       Real W_U = 0.5*(I-I_index)*(I-I_index);
		       Real W_C = 1 - W_L - W_U;
		    
		       if ((I_index-1) >= 0 && (I_index-1) < (int64_t)xDistribCells) {
			  cachedReflected[tid][parDistribIndex(I_index-1,phiIndex+J)] += amount * W_L * shapeFactors[J+1];
		       }
		       if (I_index     >= 0 && I_index < (int64_t)xDistribCells) {
			  cachedReflected[tid][parDistribIndex(I_index  ,phiIndex+J)] += amount * W_C * shapeFactors[J+1];
		       }
		       if ((I_index+1) >= 0 && (I_index+1) < (int64_t)xDistribCells) {
			  cachedReflected[tid][parDistribIndex(I_index+1,phiIndex+J)] += amount * W_U * shapeFactors[J+1];
		       }
		       if (I_inj >= 0 && I_inj < (int64_t)xDistribCells) {
			  cachedEfficiency[tid][parDistribIndex(I_inj,phiIndex)] += f*d_pitch*dU;
		       }
		    }

		    // Calculate energy of reflected particle in SIM frame:
		    Real V_reflected_SIM[3];
		    for (int k=0; k<3; ++k) V_reflected_SIM[k] = V_reflected_SNIF[k] + V_SIM_2_SNIF[k];
		    Real energy = 0.5 * species.mass * vectorMagnitude2<3>(V_reflected_SIM)
		      + particle.state[shockaccelerator::MU]*B_mag;

		    // Add electron to energy channel(s) which are stored on shock mesh:
		    uint32_t index=0;
		    for (size_t inst=0; inst<instruments.size(); ++inst) {
		       for (size_t chann=0; chann<instruments[inst].minValues.size(); ++chann) {
			  if (energy < instruments[inst].minValues[chann]) {++index; continue;}
			  if (energy >= instruments[inst].maxValues[chann]) {++index; continue;}
		     
                          #pragma omp critical
			    {
			       instrumentData[index][zoneIndex] += f*d_pitch*dU*factor1;
			    }

			  ++index;
		       }
		       ++index;
		    }

		 } else {
		    // Transmitted
		    // 
		 }
	      } else {
		 // No shock encounter
		 // 
	      }
	   } // end loop over pitch and energy samples
	} // end parallel region
      
      return true;
   }
   
   bool OperatorShockElectrons::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      if (simControl.includeShock == false) return success;
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP SHOCK MESH OP) ERROR: Shock pointer is NULL while includeShock==true" << endl << write;
	 return false;
      }
      
      // Do not plot electrons until shock is turned on:
      if (sim->t < simControl.t_shock) return true;

      // Init energy distribution unless it's already been done:
      if (getEnergy == NULL) {
	 corsair::ObjectWrapper& objectWrapper = corsair::getObjectWrapper();

	 sep::initializeEnergyDistrib initializeEnergy;
	 if (sep::getObjectWrapper().energyDistribContainer.getDistribution(energyDistribString,finalizeEnergy,
									    getDistrib,getEnergy,initializeEnergy) == false) {
	    objectWrapper.simClasses.logger << "(SEP SHOCK MESH OP) ERROR: Could not find energy distribution ";
	    objectWrapper.simClasses.logger << "function called '";
	    objectWrapper.simClasses.logger << energyDistribString << "'," << endl;
	    objectWrapper.simClasses.logger << "\t given with parameter '" << regionName+".energy_distribution'";
	    objectWrapper.simClasses.logger << endl << write;
	    success = false;
	    return success;
	 }
	 if (initializeEnergy(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader,
			      energyDistribParamsString,energyParams) == false) {
	    objectWrapper.simClasses.logger << "(SEP SHOCK MESH OP) ERROR: Energy distribution function failed ";
	    objectWrapper.simClasses.logger << "to initialize" << endl << write;
	    success = false;
	    finalizeEnergy(energyParams);
	    energyParams = NULL;
	 }
      }

      // Exit if error(s) occurred:
      if (success == false) return success;
      
      // Run tests if needed:
      //if (runTests() == false) success = false;
      
      #if PROFILE_LEVEL > 0
         profile::start("Shock Electrons",profTotalTime);
      #endif
      
      vector<Real> nodeCoordsLogical;
      simControl.shock->getLogicalNodeCoordinates(nodeCoordsLogical);

      if (writeDistributionMesh("ElectronDistribMesh") == false) success = false;
      //if (writeEmissionMesh("TypeIIMesh") == false) success = false;
      if (writeElectronDistribution(nodeCoordsLogical,"ElectronDistribMesh") == false) success = false;
      //if (writeElectronEmission(nodeCoordsLogical,"TypeIIMesh") == false) success = false;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
   bool OperatorShockElectrons::writeDistributionMesh(const std::string& meshName) {
      bool success = true;

      Real V_par_offset = 2 * fabs(V_par_min);
      yDistribCells = sim->z_blocks*block::WIDTH_Z;

      // Write node coordinates:
      map<string,string> xmlAttributes;
      xmlAttributes["mesh"] = meshName;
      vector<Real> coords(xDistribCells+1);
      for (size_t i=0; i<coords.size(); ++i) coords[i] = V_par_offset + V_par_min + i*d_V_par;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttributes,coords.size(),1,&(coords[0])) == false) success = false;
      } else {
	 Real* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      coords.resize(yDistribCells+1);
      for (size_t phi=0; phi<sim->z_blocks*block::WIDTH_Z; ++phi) {
	 coords[phi] = sim->z_crds_node[phi];
      }
      coords[coords.size()-1] = sim->z_crds_node[sim->z_blocks*block::WIDTH_Z];
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttributes,coords.size(),1,&(coords[0])) == false) success = false;
      } else {
	 Real* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      coords.resize(1);
      coords[0] = 0.0;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttributes,coords.size(),1,&(coords[0])) == false) success = false;
      } else {
	 Real* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      // Create mesh bounding box and write it:
      uint32_t bbox[6];
      bbox[0] = xDistribCells;
      bbox[1] = yDistribCells;
      bbox[2] = 0;
      bbox[3] = 1;
      bbox[4] = 1;
      bbox[5] = 1;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,6,1,bbox) == false) success = false;
      } else {
	 uint32_t* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      // Write valid global IDs:
      xmlAttributes.clear();
      xmlAttributes["name"] = meshName;
      xmlAttributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
      xmlAttributes["domains"] = "1";
      xmlAttributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
      xmlAttributes["spatial_dimension"] = "2";
      vector<uint32_t> globalIDs(xDistribCells*yDistribCells);
      for (int32_t j=0; j<yDistribCells; ++j) for (int32_t i=0; i<xDistribCells; ++i) {
	 globalIDs[j*xDistribCells+i] = j*xDistribCells+i;
      }
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,globalIDs.size(),1,&(globalIDs[0])) == false) success = false;
      } else {
	 uint32_t* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      // Write domain sizes:
      xmlAttributes.clear();
      xmlAttributes["mesh"] = meshName;
      int domainSize[2];
      domainSize[0] = xDistribCells*yDistribCells;
      domainSize[1] = 0;
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,2,domainSize) == false) success = false;
      } else {
	 int* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,0,2,NULLPTR) == false) success = false;
      }

      // Write ghosts:
      globalIDs.clear();
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (simClasses->vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttributes,0,1,&(globalIDs[0])) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttributes,0,1,&(globalIDs[0])) == false) success = false;
      } else {
	 uint32_t* NULLPTR = NULL;
	 if (simClasses->vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttributes,0,1,NULLPTR) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttributes,0,1,NULLPTR) == false) success = false;
      }

      return success;
   }

   bool OperatorShockElectrons::writeElectronDistribution(const std::vector<Real>& nodeCoords,const std::string& meshName) {
      bool success = true;

      // Check that the shock is a spherical shock:
      ShockSpherical* shockSpherical = dynamic_cast<ShockSpherical*>(simControl.shock);
      if (typeid(*simControl.shock) != typeid(*shockSpherical)) return true;
      
      #if PROFILE_LEVEL > 0
         static int profTotal = -1;
         profile::start("Electron distributions",profTotal);
      #endif
      
      // Each process creates local copies of distribution array(s):
      size_t SIZE = xDistribCells*yDistribCells;
      vector<Real> distributionThermal(SIZE);
      vector<Real> distributionReflected(SIZE);
      vector<Real> distributionTransmitted(SIZE);
      vector<Real> reflectionEfficiency(SIZE);
      for (size_t i=0; i<distributionThermal.size(); ++i) distributionThermal[i] = 0.0;
      for (size_t i=0; i<distributionReflected.size(); ++i) distributionReflected[i] = 0.0;
      for (size_t i=0; i<distributionTransmitted.size(); ++i) distributionTransmitted[i] = 0.0;
      for (size_t i=0; i<reflectionEfficiency.size(); ++i) reflectionEfficiency[i] = 0.0;

      // Each process solves electron reflection on local blocks:
      if (solveReflection(distributionThermal,distributionReflected,distributionTransmitted,reflectionEfficiency) == false) success = false;

      // Reduce values to master process:
      vector<Real> globalThermal(SIZE);
      vector<Real> globalReflected(SIZE);
      vector<Real> globalTransmitted(SIZE);
      vector<Real> globalEfficiency(SIZE);
      MPI_Reduce(&(distributionThermal[0])    ,&(globalThermal[0])    ,SIZE,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);
      MPI_Reduce(&(distributionReflected[0])  ,&(globalReflected[0])  ,SIZE,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);
      MPI_Reduce(&(distributionTransmitted[0]),&(globalTransmitted[0]),SIZE,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);
      MPI_Reduce(&(reflectionEfficiency[0])   ,&(globalEfficiency[0]) ,SIZE,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);

      // Post-process values:
      if (sim->mpiRank == sim->MASTER_RANK) {
	 for (size_t i=0; i<SIZE; ++i) globalEfficiency[i] /= (1e-30 + globalThermal[i]);
	 for (size_t i=0; i<SIZE; ++i) globalThermal[i] /= d_V_par;
	 for (size_t i=0; i<SIZE; ++i) globalReflected[i] /= d_V_par;
	 for (size_t i=0; i<SIZE; ++i) globalTransmitted[i] /= d_V_par;
      }

      // Master process writes values to file:
      map<string,string> xmlAttributes;
      xmlAttributes["mesh"] = meshName;
      xmlAttributes["units"] = "";
      xmlAttributes["name"] = meshName+"/distrib_thermal";
      xmlAttributes["centering"] = "zone";      
      const uint64_t vectorSize = 1;

      Real* ptr = NULL;
      if (sim->mpiRank != sim->MASTER_RANK) {
	 SIZE = 0;
      }

      if (sim->mpiRank == sim->MASTER_RANK) ptr = &(globalThermal[0]);
      if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,SIZE,vectorSize,ptr) == false) {
	 simClasses->logger << "(SEP SHOCK ELECTRONS OP) ERROR: Failed to write electron distribution" << endl << write;
	 success = false;
      }

      if (sim->mpiRank == sim->MASTER_RANK) ptr = &(globalReflected[0]);
      xmlAttributes["name"] = meshName+"/distrib_refl";
      if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,SIZE,vectorSize,ptr) == false) {
	 simClasses->logger << "(SEP SHOCK ELECTRONS OP) ERROR: Failed to write electron distribution" << endl << write;
	 success = false;
      }
      
      if (sim->mpiRank == sim->MASTER_RANK) ptr = &(globalTransmitted[0]);
      xmlAttributes["name"] = meshName+"/distrib_trans";
      if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,SIZE,vectorSize,ptr) == false) {
	 simClasses->logger << "(SEP SHOCK ELECTRONS OP) ERROR: Failed to write electron distribution" << endl << write;
	 success = false;
      }
      
      if (sim->mpiRank == sim->MASTER_RANK) ptr = &(globalEfficiency[0]);
      xmlAttributes["name"] = meshName+"/refl_efficiency";
      if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,SIZE,vectorSize,ptr) == false) {
	 simClasses->logger << "(SEP SHOCK ELECTRONS OP) ERROR: Failed to write electron distribution" << endl << write;
	 success = false;
      }

      vector<Real> globalInstrumentData(N_shockMeshCells);
      size_t counter=0; 
      for (size_t i=0; i<instruments.size(); ++i) {
	 for (size_t c=0; c<instruments[i].channelNames.size(); ++c) {
	    MPI_Reduce(&(instrumentData[counter][0]),&(globalInstrumentData[0]),N_shockMeshCells,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);
	    
	    if (instruments[i].divideByBinWidth == true) {
	       for (size_t k=0; k<globalInstrumentData.size(); ++k) {
		  globalInstrumentData[k] *= constants::CHARGE_ELEMENTARY/instruments[i].binWidths[c];
	       }
	    }

	    string name = instruments[i].name + '/' + instruments[i].channelNames[c];
	    shockmesh::writeVariable(sim,simClasses,false,globalInstrumentData,1,name);
	    ++counter;
	 }
      }

      counter = 0;
      vector<vector<Real> > globalLongiData = longitudeInstrData;
      for (size_t i=0; i<instruments.size(); ++i) {
	 for (size_t c=0; c<instruments[i].channelNames.size(); ++c) {
	    MPI_Reduce(&(longitudeInstrData[counter][0]),&(globalLongiData[counter][0]),yDistribCells,MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);

	    if (instruments[i].divideByBinWidth == true) {
	       for (size_t k=0; k<globalLongiData[counter].size(); ++k) {
		  globalLongiData[counter][k] *= constants::CHARGE_ELEMENTARY/instruments[i].binWidths[c];
	       }
	    }
	    ++counter;
	 }	 
      }

      // Write differential injection fluxes vs. longitude to output files:
      if (sim->mpiRank == sim->MASTER_RANK) {

	 for (size_t i=0; i<longitudeInstrData.size(); ++i) {
	    stringstream fname;
	    fname << "e_diff_flux_channel_" << i+1 << ".txt";
	    fstream out;
	    if (filesCreated[i] == false) {
	       // If output file has not been opened yet, do so and write header:
	       out.open(fname.str().c_str(),fstream::out);
	       filesCreated[i] = true;

	       out << "# ";
	       for (int j=0; j<yDistribCells; ++j) {
		  Real phi = sim->z_crds_node[j] + 0.5*sim->dz_cell[j];
		  out << phi << '\t';
	       }
	       out << endl;
	    } else {
	       // Append to existing file:
	       out.open(fname.str().c_str(),fstream::out | fstream::app);
	    }

	    out << sim->t << '\t';
	    for (int j=0; j<yDistribCells; ++j) {
	       out << globalLongiData[i][j] << '\t';
	    }
	    out << endl;
	    out.close();
	 }
      }
      
      // Write parallel speed distributions vs. longitude to output files:
      if (sim->mpiRank == sim->MASTER_RANK) {
	 stringstream fname;
	 fname << "e_par_distrib_" << sim->timestep << ".txt";
	 fstream out(fname.str().c_str(), fstream::out);

	 for (int i=0; i<xDistribCells; ++i) {
	    Real V_par = V_par_min + (i+0.5)*d_V_par;
	    out << V_par << '\t';

	    for (int j=0; j<yDistribCells; ++j) {
	       out << globalThermal[parDistribIndex(i,j)] << '\t';
	       out << globalReflected[parDistribIndex(i,j)] << '\t';
	       out << globalTransmitted[parDistribIndex(i,j)] << '\t';
	       out << globalEfficiency[parDistribIndex(i,j)] << '\t';
	    }
	    
	    out << endl;
	 }
	 out.close();
      }
      

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
