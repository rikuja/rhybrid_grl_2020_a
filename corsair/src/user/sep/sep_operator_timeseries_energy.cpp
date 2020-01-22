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

#include "sep_lagr_definition.h"
#include "sep_lagr_species.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_simcontrol.h"
#include "sep_fields_container.h"
#include "sep_guiding_center_theory.h"
#include "sep_operator_timeseries_energy.h"

using namespace std;

namespace sep {
   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   OperatorTimeseriesEnergy::OperatorTimeseriesEnergy(): DataOperator() {
      initialized = false;
      fileOut = NULL;
      namesWritten = false;
      
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorTimeseriesEnergy::~OperatorTimeseriesEnergy() {
      finalize();
   }
   
   bool OperatorTimeseriesEnergy::finalize() {
      if (initialized == false) return true;
      initialized = false;
      delete fileOut; fileOut = NULL;
      return true;
   }

   std::string OperatorTimeseriesEnergy::getName() const {return "sepOperatorTimeseriesEnergy";}
   
   bool OperatorTimeseriesEnergy::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP TIMESERIES ENERGY) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }

      const string PREFIX = "EnergyTimeseries";
      cr.add(PREFIX+".output_file_name","Name of energy timeseries output file",string("timeseries_energy.txt"));
      cr.parse();
      cr.get(PREFIX+".output_file_name",outputFileName);
      
      if (sim.mpiRank == sim.MASTER_RANK) {
	 fileOut = new fstream;
	 fileOut->open(outputFileName.c_str(),fstream::out);
	 if (fileOut->good() == false) initialized = false;
      }
      
      return initialized;
   }
   
   bool OperatorTimeseriesEnergy::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start("Energy Timeseries",profTotalTime);
      #endif
   
      // Get local cell IDs:
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
      
      // Allocate arrays for E,B vectors:
      Real B[3];
      Real E[3];
      Real gradB[9];

      // Get arrays containing (anti)parallel wave energies:
      //const uint32_t NWL = simControl.N_wavelengthMeshCells;
      Real* parWaveEnergy = simClasses->pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      Real* antiparWaveEnergy = simClasses->pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      if (parWaveEnergy == NULL) return false;
      if (antiparWaveEnergy == NULL) return false;

      vector<string> energyNames;
      vector<Real> energySums;
      
      /*
      // Calculate total energy in different wave modes:
      energyNames.push_back("AntiparallelAlfven");
      energySums.push_back(0.0);
      for (size_t i=0; i<N_localBlocks*block::SIZE*NWL; ++i) energySums[0] += antiparWaveEnergy[i];
      energyNames.push_back("ParallelAlfven");
      energySums.push_back(0.0);
      for (size_t i=0; i<N_localBlocks*block::SIZE*NWL; ++i) energySums[1] += parWaveEnergy[i];
      */
      
      // Calculate total energy in Alfven waves:
      for (size_t species=0; species<particleLists.size(); ++species) {
	 // Skip non-Lagrangian species:
	 if (particleLists[species]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
      
	 // Get Lagrangian species:
	 const sep::LagrangianSpecies* speciesData = reinterpret_cast<const sep::LagrangianSpecies*>(particleLists[species]->getSpecies());
	 if (speciesData == NULL) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Failed to get species '" << particleLists[species]->getName();
	    simClasses->logger << "' data, skipping it." << endl << write;
	    continue;
	 }
	 
	 // Get particles:
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if (particleLists[species]->getParticles(speciesDataID) == false) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Failed to get Lagrangian species '" << particleLists[species]->getName();
	    simClasses->logger << "' particles, skipping it." << endl << write;
	    continue;
	 }
	 
	 // Check that particle data is valid:
	 typedef sep::LagrangianParticle<Real> PARTICLE;
	 pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(speciesDataID);
	 if (wrapper.valid() == false) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Data for Lagrangian species '" << particleLists[species]->getName();
	    simClasses->logger << "' is invalid, skipping it." << endl << write;
	    continue;
	 }
	 
	 // Calculate total energy in this wave mode:
	 if (speciesData->propagationDirection < 0) {
	    energyNames.push_back("AntiparallelAlfven");
	    energySums.push_back(0.0);
	 } else {
	    energyNames.push_back("ParallelAlfven");
	    energySums.push_back(0.0);
	 }

	 const size_t INDEX = energySums.size()-1;
	 for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       energySums[INDEX] += wrapper.data()[blockLID][p].state[lagr::ENERGY];
	    }
	 }
      }
	 
      // Calculate total particle energy:
      for (size_t species=0; species<particleLists.size(); ++species) {
	 // Skip non-particle species:
	 if (particleLists[species]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
      
	 // Get particle species:
	 const sep::Species* speciesData = reinterpret_cast<const sep::Species*>(particleLists[species]->getSpecies());
	 if (speciesData == NULL) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Failed to get species '" << particleLists[species]->getName();
	    simClasses->logger << "' data, skipping it." << endl << write;
	    continue;
	 }
	 
	 // Get particles:
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if (particleLists[species]->getParticles(speciesDataID) == false) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Failed to get species '" << particleLists[species]->getName();
	    simClasses->logger << "' particles, skipping it." << endl << write;
	    continue;
	 }
	 
	 // Check that particle data is valid:
	 typedef sep::Particle<Real> PARTICLE;
	 pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(speciesDataID);
	 if (wrapper.valid() == false) {
	    simClasses->logger << "(SEP OP TIMESERIES ENERGY) ERROR: Data for species '" << particleLists[species]->getName();
	    simClasses->logger << "' is invalid, skipping it." << endl << write;
	    continue;
	 }

	 // TEST
	 Real U_drift_sum = 0.0;
	 // END TEST
	 
	 // Calculate total energy of particles of this species:
	 energyNames.push_back(particleLists[species]->getName());
	 energySums.push_back(0.0);
	 const size_t INDEX = energySums.size()-1;
	 for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       // Get fields at particle position:
	       (*simControl.fieldsGetFields)(blockLID,sim->t,wrapper.data()[blockLID][p].state,E,B,gradB);
	       
	       // Calculate GC drift velocity:
	       Real V_drift[3];
	       sep::calculateDriftVelocity(*speciesData,wrapper.data()[blockLID][p],E,B,gradB,V_drift);
	       
	       // Calculate total particle energy:
	       const Real U_gyro = wrapper.data()[blockLID][p].state[sep::particle::MU]*vectorMagnitude<3>(B);
	       const Real V_par2 = wrapper.data()[blockLID][p].state[sep::particle::V_PAR]*wrapper.data()[blockLID][p].state[sep::particle::V_PAR];
	       //const Real U_total = 0.5*speciesData->mass*(V_par2 + vectorMagnitude2<3>(V_drift)) + U_gyro;
	       const Real U_total = 0.5*speciesData->mass*(V_par2) + U_gyro;
	       energySums[INDEX] += U_total * wrapper.data()[blockLID][p].state[sep::particle::WEIGHT];
	       
	       // TEST
	       U_drift_sum += 0.5*speciesData->mass*vectorMagnitude2<3>(V_drift)*wrapper.data()[blockLID][p].state[sep::particle::WEIGHT];
	       // END TEST
	    }
	 }
      }
      
      // Collect total energy changes of particles and waves to master process:
      energyNames.push_back("dU_particles");
      energySums.push_back(simControl.dU_particles);
      energyNames.push_back("dU_waves");
      energySums.push_back(simControl.dU_waves);
      energyNames.push_back("dU_difference");
      energySums.push_back(simControl.dU_particles+simControl.dU_waves);

      // Collect total energies to master process:
      vector<Real> globalEnergySums;
      globalEnergySums.resize(energySums.size());
      if (energySums.size() > 0) {
	 MPI_Reduce(&(energySums[0]),&(globalEnergySums[0]),energySums.size(),MPI_Type<Real>(),MPI_SUM,sim->MASTER_RANK,sim->comm);
      }
      
      // Write out data (master process only):
      if (sim->mpiRank == sim->MASTER_RANK) {
	 if (namesWritten == false) {
	    (*fileOut) << "Time(s)" << '\t';
	    for (size_t i=0; i<energyNames.size(); ++i) (*fileOut) << energyNames[i] << '\t';
	    (*fileOut) << "Total Energy" << endl;
	    namesWritten = true;
	 }
	 
	 fileOut->precision(14);
	 Real totalEnergy = 0.0;
	 (*fileOut) << sim->t << '\t';
	 for (size_t i=0; i<globalEnergySums.size(); ++i) {
	    (*fileOut) << globalEnergySums[i] << '\t';
	    totalEnergy += globalEnergySums[i];
	 }
	 (*fileOut) << totalEnergy;
	 (*fileOut) << endl;
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
      
   
} // namespace sep
