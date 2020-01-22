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

#ifndef SEP_PARTICLE_INJECTOR_HOMOGENEOUS_H
#define SEP_PARTICLE_INJECTOR_HOMOGENEOUS_H

#include <cstdlib>
#include <climits>
#include <stdint.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include <linear_algebra.h>

#include "sep_object_wrapper.h"
#include "sep_namespaces.h"
#include "sep_simcontrol.h"
#include "sep_particle_definition.h"
#include "sep_fields_container.h"
#include "sep_injectors_common.h"

namespace sep {

   extern sep::SimControl simControl;

   template<class SPECIES,class PARTICLE>
   class ParticleInjectorHomog: public ParticleInjectorBase {
    public: 
      ParticleInjectorHomog();
      ~ParticleInjectorHomog();
   
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		       const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);
   
    private:
      bool initialInjection;
      uint32_t inflowBoundaryMask;
      bool initialized;
      bool initialInjectionDone;
      uint32_t N_particlesPerCell;        /**< Number of macroparticles injected to each cell.*/
      Real numberDensity;                 /**< Number density injected to each cell.*/
      SPECIES species;
      particleinjector::frame::Frame injectionFrame;
      
      void* energyParams;
      sep::finalizeEnergyDistrib finalizeEnergy;
      sep::getEnergyFunction getEnergy;
      sep::getPitchFunction getPitch;
      
      const Real DEF_VALUE;
      
      void injectParticles(pargrid::CellID block,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* PIHomogMaker() {return new ParticleInjectorHomog<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorHomog<SPECIES,PARTICLE>::ParticleInjectorHomog(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity())  {
      getEnergy = NULL;
      finalizeEnergy = NULL;
      getPitch = NULL;
      energyParams = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorHomog<SPECIES,PARTICLE>::~ParticleInjectorHomog() {
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorHomog<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".energy_distribution","Name of energy distribution function (string)",std::string(""));
      cr.add(PREFIX+".energy_distribution_parameters","Name of region containing energy distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution","Name of pitch distribution function (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution_parameters","Name of region containing pitch distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".macroparticles_per_cell","Number of macroparticles injected per cell (integer)",(uint32_t)0);
      cr.add(PREFIX+".number_density","Number density injected to each cell (float)",DEF_VALUE);
      cr.add(PREFIX+".inflow_boundary","Inflow boundary where particles are injected at every time step (string)",std::string(""));
      cr.add(PREFIX+".initial_injection","If 'yes' particles are injected to entire simulation domain before first time step (string)",std::string("yes"));
      cr.add(PREFIX+".injection_frame","Injection frame, antiparallel/parallel/simulation (string)",std::string("simulation"));
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorHomog<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      initialInjectionDone = false;
      sim = NULL;
      simClasses = NULL;
      if (finalizeEnergy != NULL) (*finalizeEnergy)(energyParams);
      finalizeEnergy = NULL;
      getEnergy = NULL;
      getPitch = NULL;
      return true;
   }

   /** Called by ParticleListSkeleton.*/
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorHomog<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      bool success = true;
      
      // Do not inject particles until setup time has passed:
      if (sim->t < simControl.t_setup) return success;

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      /*
      // TEST
      if (sim->t >= simControl.t_setup && initialInjectionDone == false) {
	 const std::vector<pargrid::CellID>& interiorBlocks = simClasses->pargrid.getInteriorCells();
	 for (size_t b=0; b<interiorBlocks.size(); ++b) {
	    const pargrid::CellID blockLID = interiorBlocks[b];
	    injectParticles(blockLID,N_particles,wrapper);
	 }
	 initialInjectionDone = true;
      }
      // END TEST
       */
      // Initial injection:
      Real t_propag = 0.0;
      if (initialInjection == true && initialInjectionDone == false) {
	 const std::vector<pargrid::CellID>& interiorBlocks = simClasses->pargrid.getInteriorCells();
	 for (size_t b=0; b<interiorBlocks.size(); ++b) {
	    // Measure block injection time if we are testing for repartitioning:
	    if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	    
	    const pargrid::CellID blockLID = interiorBlocks[b];
	    injectParticles(blockLID,N_particles,wrapper);
	    
	    // Store block injection time:
	    if (sim->countPropagTime == true) {
	       t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	       simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	    }
	 }
	 initialInjectionDone = true;
      }

      // Inject new particles to outflow boundaries:
      if (sim->initializing == true) return success;
      const std::vector<pargrid::CellID> exteriorBlocks = simClasses->pargrid.getExteriorCells();
      for (size_t block=0; block<exteriorBlocks.size(); ++block) {
	 // Measure block injection time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();

	 const pargrid::CellID blockLID = exteriorBlocks[block];
	 const uint32_t existingNbrFlags = simClasses->pargrid.getNeighbourFlags()[blockLID];
	 const uint32_t missingNbrFlags = (existingNbrFlags ^ pargrid::ALL_NEIGHBOURS_EXIST);
	 if ((missingNbrFlags & inflowBoundaryMask) > 0) {
	    injectParticles(blockLID,N_particles,wrapper);
	 }

	 // Store block injection time:
	 if (sim->countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
      
      return success;
   }

   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorHomog<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,unsigned int* N_particles,
								 pargrid::DataWrapper<PARTICLE>& wrapper) {
      // Resize block's particle list to have space for new particles:
      const pargrid::ArraySizetype oldSize = wrapper.size()[blockLID];
      N_particles[blockLID] += N_particlesPerCell*block::SIZE;
      wrapper.resize(blockLID,oldSize+N_particlesPerCell*block::SIZE);

      // Get block global ID and calculate bounding box (i,j,k) indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      uint32_t i_block,j_block,k_block;
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

      //Real E[3];
      Real B[3];
      //Real gradB[9];
      Real V_wave[3];
      Real dV_wave;
      distrib::InjectionEnergy injEnergy;
      injEnergy.thermalEnergy = constants::BOLTZMANN*1e6;
      
      int32_t alfvenSign = +1;
      if (injectionFrame == particleinjector::frame::ANTIPARALLEL) alfvenSign = -1;

      Real weightSum = 0.0;
      unsigned int counter = oldSize;
      PARTICLE* particles = wrapper.data()[blockLID];
      for (int32_t k=0; k<block::WIDTH_Z; ++k) for (int32_t j=0; j<block::WIDTH_Y; ++j) for (int32_t i=0; i<block::WIDTH_X; ++i) {
	 // Calculate cell (i,j,k) indices:
	 const uint32_t i_cell = i_block*block::WIDTH_X + i;
	 const uint32_t j_cell = j_block*block::WIDTH_Y + j;
	 const uint32_t k_cell = k_block*block::WIDTH_Z + k;
	 
	 // Inject particles to cell:
	 injEnergy.N_particles = N_particlesPerCell;
	 for (uint32_t n=0; n<N_particlesPerCell; ++n) {
	    injEnergy.N = n;
	    injEnergy.superResolution = 1;

	    // Calculate injection position in logical coordinates:
	    particles[counter].state[particle::XCRD] = i_cell + simClasses->random.uniform();
	    particles[counter].state[particle::YCRD] = j_cell + simClasses->random.uniform();
	    particles[counter].state[particle::ZCRD] = k_cell + simClasses->random.uniform();

	    // Get magnetic field and wave speed (significant only if injecting in wave rest frame) at injection position:
	    #warning Injector assumes simulation time sim->t
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,particles[counter].state,B,V_wave,dV_wave,V_alfven,alfvenSign);
	    const Real B_mag = vectorMagnitude<3>(B);
	    
	    // Determine injection frame parallel speed:
	    Real frameSpeed = 0.0;
	    if (injectionFrame == particleinjector::frame::PARALLEL || injectionFrame == particleinjector::frame::ANTIPARALLEL) {
	       frameSpeed = dotProduct<3>(V_wave,B)/B_mag;
	    }

	    // Calculate pitch, energy:
	    // NOTE: energy and speed are calculated in physical units
	    (*getEnergy)(injEnergy,energyParams);
	    #ifndef NDEBUG
	    if (injEnergy.weight != injEnergy.weight) {
	       std::cerr << "(SEP PARTICLE INJ HOMOG) ERROR: Received NAN particle weight" << std::endl;
	       exit(1);
	    }
	    #endif
	    particles[counter].state[particle::WEIGHT] = injEnergy.weight;
	    const Real pitch  = (*getPitch)();
	    const Real driftFreeEnergy = injEnergy.energy;
	    const Real driftFreeSpeed = sqrt(2*driftFreeEnergy/species.mass);
	    const Real parallelSpeed  = driftFreeSpeed * pitch;
	    const Real gyroSpeed      = driftFreeSpeed * sqrt(1 - pitch*pitch);
	    const Real magneticMoment = 0.5*species.mass*gyroSpeed*gyroSpeed / B_mag;

	    particles[counter].state[particle::V_PAR] = parallelSpeed + frameSpeed;
	    particles[counter].state[particle::MU]    = magneticMoment;
	    weightSum += particles[counter].state[particle::WEIGHT];

	    ++counter;
	 }

	 // Get cell volume:
	 const Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	 
	 // Normalize weights:
	 counter -= N_particlesPerCell;
	 Real sum = 0.0;
	 for (uint32_t n=0; n<N_particlesPerCell; ++n) {
	    particles[counter].state[particle::WEIGHT] *= (numberDensity/weightSum)*cellVolume;
	    sum += particles[counter].state[particle::WEIGHT];
	    ++counter;
	 }
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorHomog<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
							    const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP PARTICLE INJ HOMOG) Starting to init." << std::endl << write;
      
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      this->species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      // Parse config file options:
      std::string energyDistribString,energyDistribParamsString;
      std::string pitchDistribString,pitchDistribParamsString;
      std::string inflowBoundaryString,initialInjectionString;
      std::string injectionFrameString;
      cr.parse();
      cr.get(regionName+".energy_distribution",energyDistribString);
      cr.get(regionName+".energy_distribution_parameters",energyDistribParamsString);
      cr.get(regionName+".pitch_distribution",pitchDistribString);
      cr.get(regionName+".pitch_distribution_parameters",pitchDistribParamsString);
      cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell);
      cr.get(regionName+".number_density",numberDensity);
      cr.get(regionName+".inflow_boundary",inflowBoundaryString);
      cr.get(regionName+".initial_injection",initialInjectionString);
      cr.get(regionName+".injection_frame",injectionFrameString);

      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Parameter '" << regionName+".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
      }
      if (numberDensity == DEF_VALUE) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Parameter '" << regionName+".number_density' was not found" << std::endl << write;
	 initialized = false;
      }

      injectionFrame = getInjectionFrame(injectionFrameString);
      if (injectionFrame == sep::particleinjector::frame::UNKNOWN) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Unknown injection frame" << std::endl << write;
	 initialized = false;
      }

      // Determine inflow boundary (if any):
      inflowBoundaryMask = 0;
      if (inflowBoundaryString == "-x") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0));
      if (inflowBoundaryString == "+x") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0));
      if (inflowBoundaryString == "-y") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0));
      if (inflowBoundaryString == "+y") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0));
      if (inflowBoundaryString == "-z") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1));
      if (inflowBoundaryString == "+z") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1));

      // Check if initial state injection should be done:
      if (initialInjectionString == "yes") initialInjection = true;
      else initialInjection = false;

      // Attempt to get pitch distribution function and initialize it:
      sep::finalizePitchDistrib finalizePitch;
      sep::initializePitchDistrib initializePitch;
      if (sep::getObjectWrapper().pitchDistribContainer.getDistribution(pitchDistribString,finalizePitch,getPitch,initializePitch) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Could not find pitch distribution function called '";
	 simClasses.logger << pitchDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".pitch_distribution'" << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      if (initializePitch(sim,simClasses,cr,pitchDistribParamsString) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Pitch distribution function failed to initialize" << std::endl << write;
	 initialized = false;
      }
      
      // Attempt to get energy distribution function:
      sep::initializeEnergyDistrib initializeEnergy;
      sep::getEnergyDistribFunction getDistrib;
      if (sep::getObjectWrapper().energyDistribContainer.getDistribution(energyDistribString,finalizeEnergy,getDistrib,getEnergy,initializeEnergy) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Could not find energy distribution function called '";
	 simClasses.logger << energyDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".energy_distribution'" << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      if (initializeEnergy(sim,simClasses,cr,energyDistribParamsString,energyParams) == false) {
	 simClasses.logger << "(PARTICLE INJ HOMOG) ERROR: Energy distribution function failed to initialize" << std::endl << write;
	 (*finalizeEnergy)(energyParams);
	 initialized = false;
      }

      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }

      simClasses.logger << "(PARTICLE INJ HOMOG) Initialization complete, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << std::endl << write;
      else simClasses.logger << "FAILURE" << std::endl << write;

      return initialized;
   }

} // namespace sep
   
#endif
