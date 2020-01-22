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

#ifndef SEP_PARTICLE_INJECTOR_CELL_H
#define SEP_PARTICLE_INJECTOR_CELL_H

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
#include "sep_simcontrol.h"
#include "sep_particle_definition.h"
#include "sep_propagate.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class ParticleInjectorCell: public ParticleInjectorBase {
    public: 
      ParticleInjectorCell();
      ~ParticleInjectorCell();
      
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);
   
    private:
      bool initialized;
      bool initialInjectionDone;
      Real injectionPosition[3];
      uint32_t N_particlesPerCell;        /**< Number of macroparticles injected to each cell.*/
      Real numberDensity;                 /**< Number density injected to each cell.*/
      SPECIES species;

      void* energyParams;                 /**< Parameters to getEnergy function.*/
      sep::finalizeEnergyDistrib finalizeEnergy;
      sep::getEnergyFunction getEnergy;
      sep::getPitchFunction getPitch;
      
      const Real DEF_VALUE;
      
      void injectParticles(pargrid::CellID block,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper,
			  int i_cell,int j_cell,int k_cell);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* PICellMaker() {return new ParticleInjectorCell<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorCell<SPECIES,PARTICLE>::ParticleInjectorCell(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity())  {
      finalizeEnergy = NULL;
      getEnergy = NULL;
      getPitch = NULL;
      energyParams = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorCell<SPECIES,PARTICLE>::~ParticleInjectorCell() {
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorCell<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".energy_distribution","Name of energy distribution function (string)",std::string(""));
      cr.add(PREFIX+".energy_distribution_parameters","Name of region containing energy distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution","Name of pitch distribution function (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution_parameters","Name of region containing pitch distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".macroparticles_per_cell","Number of macroparticles injected per cell (integer)",(uint32_t)0);
      cr.add(PREFIX+".number_density","Number density injected to each cell (float)",DEF_VALUE);
      cr.add(PREFIX+".position_x","x-coordinate of injection position (float)",DEF_VALUE);
      cr.add(PREFIX+".position_y","y-coordinate of injection position (float)",DEF_VALUE);
      cr.add(PREFIX+".position_z","z-coordinate of injection position (float)",DEF_VALUE);
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorCell<SPECIES,PARTICLE>::finalize() {
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
   bool ParticleInjectorCell<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      bool success = true;

      recalculateCellVolumes(*sim,*simClasses);
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
            
      Real t_propag = 0.0;
      if (initialInjectionDone == false) {
	 const std::vector<pargrid::CellID>& interiorBlocks = simClasses->pargrid.getInteriorCells();
	 for (size_t b=0; b<interiorBlocks.size(); ++b) {
	    // Measure block injection time if we are testing for repartitioning:
	    if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	    
	    const pargrid::CellID blockLID = interiorBlocks[b];
	    const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	    uint32_t i_block,j_block,k_block;
	    block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	    
	    for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	       const Real x = sim->x_crds_node[i_block*block::WIDTH_X+i];
	       const Real y = sim->y_crds_node[j_block*block::WIDTH_Y+j];
	       const Real z = sim->z_crds_node[k_block*block::WIDTH_Z+k];
	       const Real dx = sim->dx_cell[i_block*block::WIDTH_X+i];
	       const Real dy = sim->dy_cell[j_block*block::WIDTH_Y+j];
	       const Real dz = sim->dz_cell[k_block*block::WIDTH_Z+k];
	       
	       bool inject = true;
	       if (x > injectionPosition[0] || x+dx <= injectionPosition[0]) inject = false; 
	       if (y > injectionPosition[1] || y+dy <= injectionPosition[1]) inject = false;
	       if (z > injectionPosition[2] || z+dz <= injectionPosition[2]) inject = false;
	       if (inject == true) {
		  injectParticles(blockLID,N_particles,wrapper,i,j,k);
	       }
	    }
	    // Store block injection time:
	    if (sim->countPropagTime == true) {
	       t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	       simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	    }
	 }
	 initialInjectionDone = true;
      }
	 
      return success;
   }

   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorCell<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,unsigned int* N_particles,
								 pargrid::DataWrapper<PARTICLE>& wrapper,int i_cell,int j_cell,int k_cell) {
      // Resize block's particle list to have space for new particles:
      const pargrid::ArraySizetype oldSize = wrapper.size()[blockLID];
      N_particles[blockLID] += N_particlesPerCell*block::SIZE;
      wrapper.resize(blockLID,oldSize+N_particlesPerCell*block::SIZE);

      // Get block global ID and calculate bounding box (i,j,k) indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      uint32_t i_block,j_block,k_block;
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      // Get spatial cell volume:
      const Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i_cell,j_cell,k_cell)];
      const Real weight = numberDensity*cellVolume/N_particlesPerCell;

      Real E[3];
      Real B[3];
      Real gradB[9];
      
      unsigned int counter = oldSize;
      PARTICLE* particles = wrapper.data()[blockLID];
      
      distrib::InjectionEnergy injEnergy;
      
      // Inject particles to cell:
      for (uint32_t n=0; n<N_particlesPerCell; ++n) {
	 // Calculate injection position in logical coordinates:
	 particles[counter].state[particle::XCRD] = i_block*block::WIDTH_X + i_cell + simClasses->random.uniform();
	 particles[counter].state[particle::YCRD] = j_block*block::WIDTH_Y + j_cell + simClasses->random.uniform();
	 particles[counter].state[particle::ZCRD] = k_block*block::WIDTH_Z + k_cell + simClasses->random.uniform();
	 //particles[counter].state[particle::XCRD] = i_cell + 0.5;
	 //particles[counter].state[particle::YCRD] = j_cell + 0.5;
	 //particles[counter].state[particle::ZCRD] = k_cell + 0.5;
	    
	 // Get fields at injection position:
	 #warning Injector assumes simulation time sim->t
	 (*simControl.fieldsGetFields)(blockLID,sim->t,particles[counter].state,E,B,gradB);

	 // Calculate pitch, energy:
	 // NOTE: energy and speed are calculated in physical units
	 (*getEnergy)(injEnergy,energyParams);
	 particles[counter].state[sep::particle::WEIGHT] = injEnergy.weight;
	 const Real pitch  = (*getPitch)();

	 // Transform to solar wind frame:
	 /*
	  Real V_GC[3];
	  crossProduct(E,B,V_GC);
	  const Real B_mag2 = vectorMagnitude2<3>(B);
	  V_GC[0] /= B_mag2;
	  V_GC[1] /= B_mag2;
	  V_GC[2] /= B_mag2;
	  const Real driftEnergy = 0.5*species.m*vectorMagnitude2<3>(V_GC);
	  */
	 const Real driftFreeEnergy = injEnergy.energy;
	 
	 const Real driftFreeSpeed = sqrt(2*driftFreeEnergy/species.mass);
	 const Real parallelSpeed  = driftFreeSpeed * pitch;
	 const Real gyroSpeed      = driftFreeSpeed * sqrt(1 - pitch*pitch);
	 
	 const Real B_mag = vectorMagnitude<3>(B);
	 const Real magneticMoment = 0.5*species.mass*gyroSpeed*gyroSpeed / B_mag;
	 
	 particles[counter].state[particle::V_PAR] = parallelSpeed;
	 particles[counter].state[particle::MU]    = magneticMoment;
	 particles[counter].state[particle::WEIGHT] = weight;
	 ++counter;
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorCell<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
							   const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP PARTICLE INJ CELL) Starting to init." << std::endl << write;
      
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      this->species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      // Parse config file options:
      std::string energyDistribString,energyDistribParamsString;
      std::string pitchDistribString,pitchDistribParamsString;
      cr.parse();
      cr.get(regionName+".energy_distribution",energyDistribString);
      cr.get(regionName+".energy_distribution_parameters",energyDistribParamsString);
      cr.get(regionName+".pitch_distribution",pitchDistribString);
      cr.get(regionName+".pitch_distribution_parameters",pitchDistribParamsString);
      cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell);
      cr.get(regionName+".number_density",numberDensity);
      cr.get(regionName+".position_x",injectionPosition[0]);
      cr.get(regionName+".position_y",injectionPosition[1]);
      cr.get(regionName+".position_z",injectionPosition[2]);

      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Parameter '" << regionName+".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
      }
      if (numberDensity == DEF_VALUE) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Parameter '" << regionName+".number_density' was not found" << std::endl << write;
	 initialized = false;
      }
      for (int i=0; i<3; ++i) if (injectionPosition[i] == DEF_VALUE) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Parameter(s) '" << regionName+".position_x/y/z' contained invalid value(s)" << std::endl << write;
	 initialized = false;
      }
      
      // Attempt to get pitch distribution function and initialize it:
      sep::finalizePitchDistrib finalizePitch;
      sep::initializePitchDistrib initializePitch;
      if (sep::getObjectWrapper().pitchDistribContainer.getDistribution(pitchDistribString,finalizePitch,getPitch,initializePitch) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Could not find pitch distribution function called '";
	 simClasses.logger << pitchDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".pitch_distribution'" << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      if (initializePitch(sim,simClasses,cr,pitchDistribParamsString) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Pitch distribution function failed to initialize" << std::endl << write;
	 initialized = false;
      }
      
      // Attempt to get energy distribution function:
      sep::initializeEnergyDistrib initializeEnergy;
      sep::getEnergyDistribFunction getDistrib;
      if (sep::getObjectWrapper().energyDistribContainer.getDistribution(energyDistribString,finalizeEnergy,
									 getDistrib,getEnergy,initializeEnergy) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Could not find energy distribution ";
	 simClasses.logger << "function called '";
	 simClasses.logger << energyDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".energy_distribution'"
	   << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      if (initializeEnergy(sim,simClasses,cr,energyDistribParamsString,energyParams) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ CELL) ERROR: Energy distribution function failed ";
	 simClasses.logger << "to initialize" << std::endl << write;
	 initialized = false;
	 finalizeEnergy(energyParams);
      }

      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }

      simClasses.logger << "(SEP PARTICLE INJ CELL) Initialization complete, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << std::endl << write;
      else simClasses.logger << "FAILURE" << std::endl << write;
      
      return initialized;
   }

} // namespace sep
   
#endif
