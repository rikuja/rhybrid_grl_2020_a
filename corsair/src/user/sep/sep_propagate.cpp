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

#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_particle_scatterer.h"
#include "sep_fields_container.h"
#include "sep_lagr_definition.h"
#include "sep_lagr_species.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_accumulation_stretched.h"
#include "sep_debugging.h"
#include "sep_injection_buffer.h"
#include "sep_particle_list.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   /** Helper function that stops profiling and exits a function with given status.
    * @param status Exit status.
    * @return Always equal the value of status.*/
   bool exitWithStatus(const bool& status) {
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return status;
   }
   
   bool accumulateParticlePhaseSpaceWeight(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      
      // Clear accumulation array(s):
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 if (particleLists[p]->clearAccumulationArrays() == false) success = false;
      }        
      
      // Accumulate particle quantities to simulation mesh:
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 if (particleLists[p]->accumulateBoundaryCells() == false) success = false;
      }
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 if (particleLists[p]->accumulateInnerCells() == false) success = false;
      }
      
      // Divide values by phase-space volume:
      Real* weight4D = simClasses.pargrid.getUserDataStatic<Real>(simControl.particleWeightDataID);
      const int NWL = simControl.N_wavelengthMeshCells;
      for (pargrid::CellID blockLID=0; blockLID<simClasses.pargrid.getNumberOfAllCells(); ++blockLID) {
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);
	 
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    const Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	    
	    for (int l=0; l<NWL; ++l) {
	       const Real factor = 1.0 / simControl.wavelengthMeshCellSizes[l] / cellVolume;
	       weight4D[blockLID*block::SIZE*NWL+block::index(i,j,k)*NWL+l] *= factor;
	    }
	 }
      }
      
      if (success == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: accumulateParticlePhaseSpaceWeight failed" << endl << write;
      }
      return success;
   }
   
   bool accumulateWaveEnergies(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      
      // Clear accumulation array(s):
      for (size_t p=0; p<particleLists.size(); ++p) {
	 // Skip non-Lagrangian species:
	 if (particleLists[p]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
	 
	 // Only accumulate current Alfven wave mode energy (parallel or antiparallel):
	 const sep::LagrangianSpecies* species = reinterpret_cast<const sep::LagrangianSpecies*>(particleLists[p]->getSpecies());
	 if (species->propagationDirection * simControl.alfvenSign < 0) continue;
	 
	 if (particleLists[p]->clearAccumulationArrays() == false) success = false;
      }

      // Accumulate particle quantities to simulation mesh:
      for (size_t p=0; p<particleLists.size(); ++p) {
	 // Skip non-Lagrangian species:
	 if (particleLists[p]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
	 
	 // Only accumulate current Alfven wave mode energy (parallel or antiparallel):
	 const sep::LagrangianSpecies* species = reinterpret_cast<const sep::LagrangianSpecies*>(particleLists[p]->getSpecies());
	 if (species->propagationDirection * simControl.alfvenSign < 0) continue;

	 if (particleLists[p]->accumulateBoundaryCells() == false) success = false;
      }
      for (size_t p=0; p<particleLists.size(); ++p) {
	 // Skip non-Lagrangian species:
	 if (particleLists[p]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;

	 // Only accumulate current Alfven wave mode energy (parallel or antiparallel):
	 const sep::LagrangianSpecies* species = reinterpret_cast<const sep::LagrangianSpecies*>(particleLists[p]->getSpecies());
	 if (species->propagationDirection * simControl.alfvenSign < 0) continue;

	 if (particleLists[p]->accumulateInnerCells() == false) success = false;
      }

      if (success == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: accumulateWaveEnergies failed" << endl << write;
      }

      return success;
   }
   
   void checkCoordinates(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists,
			 const string& step) {
      for (size_t p=0; p<particleLists.size(); ++p) {
	 bool success = true;
	 if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) {
	    typedef LagrangianParticle<Real> LAGR;
	    typedef LagrangianSpecies SPECIES;
	    const SPECIES* species = reinterpret_cast<const SPECIES*>(particleLists[p]->getSpecies());
	    if (sep::checkCoordinates<SPECIES,LAGR>(sim,simClasses,*species,particleLists[p],step) == false)
	      success = false;
	 } else {
	    typedef Particle<Real> PARTICLE;
	    typedef Species SPECIES;
	    const SPECIES* species = reinterpret_cast<const SPECIES*>(particleLists[p]->getSpecies());
	    if (sep::checkCoordinates<SPECIES,PARTICLE>(sim,simClasses,*species,particleLists[p],step) == false)
	      success = false;
	 }
	 if (success == false) exit(1);
      }
   }

   bool clearWaveGrowthArray(Simulation& sim,SimulationClasses& simClasses) {
      Real* waveGrowth = NULL;
      if (simControl.alfvenSign < 0) {
	 waveGrowth = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparWaveGrowthDataID);
      } else {
	 waveGrowth = simClasses.pargrid.getUserDataStatic<Real>(simControl.parWaveGrowthDataID);
      }
      if (waveGrowth == NULL) return false;
      
      #if PROFILE_LEVEL > 0
         static int profClearWaveGrowth = -1;
         profile::start("Clear Wave Growth Array",profClearWaveGrowth);
      #endif
      
      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      for (size_t i=0; i<simClasses.pargrid.getNumberOfAllCells()*block::SIZE*NWL; ++i) {
	 waveGrowth[i] = 0.0;
      }

      #if PROFILE_LEVEL > 0
         profile::stop(); // Clear wave growth array
      #endif
      
      return true;
   }

   /** During shock crossing incident waves will split into transmitted
    * and reflected waves in downstream region. This function injects 
    * reflected waves into simulation.*/
   bool injectExtraParticles() {
      bool success = true;
      /*
      #if PROFILE_LEVEL > 0
         static int profID = -1;
         profile::start("Inject Reflected Waves",profID);
      #endif
      */
      // Iterate over all particle lists and pick the ones that contain
      // Alfven wave packets:
      corsair::ObjectWrapper& objectWrapper = corsair::getObjectWrapper();
      for (size_t plist=0; plist<objectWrapper.particleLists.size(); ++plist) {
	 if (objectWrapper.particleLists[plist]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;

	 // Get Alfven wave packet species, needed to get 
	 // propagation direction below:
	 const sep::LagrangianSpecies* species = 
	   reinterpret_cast<const sep::LagrangianSpecies*>(objectWrapper.particleLists[plist]->getSpecies());

	 // Get ParGrid DataID of particles in particle list:
	 pargrid::DataID speciesID = pargrid::INVALID_DATAID;
	 if (objectWrapper.particleLists[plist]->getParticles(speciesID) == false) {
	    success = false;
	    #ifndef NDEBUG
	    cerr << "ERROR: Lagrangian species DataID is INVALID" << endl;
	    #endif
	    continue;
	 }
	 
	 // Get particle data:
	 pargrid::DataWrapper<LagrangianParticle<Real> > wrapper = 
	   objectWrapper.simClasses.pargrid.getUserDataDynamic<LagrangianParticle<Real> >(speciesID);
	 if (wrapper.valid() == false) {
	    success = false;
	    #ifndef NDEBUG
	    cerr << "ERROR: Lagrangian wrapper is invalid" << endl;
	    #endif
	    continue;
	 }
	 
	 // These two injection buffers contain reflected wave packets that need
	 // to be injected to simulation:
	 extern sep::InjectionBuffer<LagrangianParticle<Real> > antiparInjectionBuffer;
	 extern sep::InjectionBuffer<LagrangianParticle<Real> > parInjectionBuffer;
	 
	 // Copy all particles from injection buffer to parallel mesh:
	 uint32_t* N_particles = objectWrapper.particleLists[plist]->getParticleNumberArray();
	 if (species->propagationDirection < 0) {
	    for (size_t p=0; p<antiparInjectionBuffer.cellIDs.size(); ++p) {
	       wrapper.push_back(antiparInjectionBuffer.cellIDs[p],antiparInjectionBuffer.particles[p]);
	       if (antiparInjectionBuffer.cellIDs[p] < objectWrapper.simClasses.pargrid.getNumberOfLocalCells()) {
		  ++N_particles[antiparInjectionBuffer.cellIDs[p]];
	       }
	    }
	    antiparInjectionBuffer.clear();
	 } else {
	    for (size_t p=0; p<parInjectionBuffer.cellIDs.size(); ++p) {
	       wrapper.push_back(parInjectionBuffer.cellIDs[p],parInjectionBuffer.particles[p]);
	       if (parInjectionBuffer.cellIDs[p] < objectWrapper.simClasses.pargrid.getNumberOfLocalCells()) {
		  ++N_particles[parInjectionBuffer.cellIDs[p]];
	       }
	    }	       
	    parInjectionBuffer.clear();
	 }

	 // Send particles from buffered (remote) cells to neighbor processes:
	 sep::ParticleList<LagrangianSpecies,LagrangianParticle<Real> >* particleList =
	   reinterpret_cast<sep::ParticleList<LagrangianSpecies,LagrangianParticle<Real> >*>(objectWrapper.particleLists[plist]);
	 if (particleList->exchangeParticles() == false) success = false;
      }

      extern map<string,sep::InjectionBuffer<sep::Particle<Real> > > particleInjectionBuffers;
      for (size_t plist=0; plist<objectWrapper.particleLists.size(); ++plist) {
	 if (objectWrapper.particleLists[plist]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 
	 map<string,sep::InjectionBuffer<Particle<Real> > >::iterator it 
	   = particleInjectionBuffers.find(objectWrapper.particleLists[plist]->getName());
	 if (it == particleInjectionBuffers.end()) continue;
	 
	 // Get ParGrid DataID of particles in particle list:
	 pargrid::DataID speciesID = pargrid::INVALID_DATAID;
	 if (objectWrapper.particleLists[plist]->getParticles(speciesID) == false) {
	    success = false;
	    #ifndef NDEBUG
	    cerr << "ERROR: Particle DataID is INVALID" << endl;
	    #endif
	    continue;
	 }

	 // Get particle data:
	 pargrid::DataWrapper<Particle<Real> > wrapper =
	   objectWrapper.simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesID);
	 if (wrapper.valid() == false) {
	    success = false;
	    #ifndef NDEBUG
	    cerr << "ERROR: Particle wrapper is not valid" << endl;
            #endif
	    continue;
	 }
	 
	 // Copy all particles from injection buffer to parallel mesh:
	 uint32_t* N_particles = objectWrapper.particleLists[plist]->getParticleNumberArray();
	 for (size_t p=0; p<it->second.cellIDs.size(); ++p) {
	    wrapper.push_back(it->second.cellIDs[p],it->second.particles[p]);
	    
	    // If we are copying a particle to local cell, we need to increase N_particles 
	    // so that propagator knows to propagate all particles.
	    if (it->second.cellIDs[p] < objectWrapper.simClasses.pargrid.getNumberOfLocalCells()) {
	       ++N_particles[it->second.cellIDs[p]];
	    }
	 }
	 particleInjectionBuffers.erase(it);
	 
	 // Send particles from buffered (remote) cells to neighbor processes:
	 sep::ParticleList<Species,Particle<Real> >* particleList =
	   reinterpret_cast<sep::ParticleList<Species,Particle<Real> >*>(objectWrapper.particleLists[plist]);
	 if (particleList->exchangeParticles() == false) success = false;
      }

      /*
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif*/
      return success;
   }
   
   void recalculateCellVolumes(Simulation& sim,SimulationClasses& simClasses) {
      // Prevent multiple recalculations on same time step:
      if ((int64_t)sim.timestep == simControl.timestepCellVolumesCalculated) return;
      if (sim.meshRepartitioned == false) return;

      delete [] simControl.cellVolumes;
      simControl.cellVolumes = new Real[simClasses.pargrid.getNumberOfAllCells()*block::SIZE];
      
      size_t counter = 0;
      for (pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfAllCells(); ++b) {
	 // Calculate block's indices:
	 const pargrid::CellID blockLID = b;
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 uint32_t i_cell,j_cell,k_cell;
	 block::calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);

	 // Calculate spatial cell volumes for all cells in block:
	 Real dx = NAN;
	 Real dy = NAN;
	 Real dz = NAN;
	 Real X0,Y0,Z0,X1,Y1,Z1;
	 const Real ZERO = 0.0;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    i_cell = i_block*block::WIDTH_X + i;
	    j_cell = j_block*block::WIDTH_Y + j;
	    k_cell = k_block*block::WIDTH_Z + k;
	    switch (simControl.coordinateSystem) {
	     case sep::UNKNOWN:
	       break;
	     case sep::CARTESIAN:
	       dx = sim.x_crds_node[i_cell+1]-sim.x_crds_node[i_cell];
	       dy = sim.y_crds_node[j_cell+1]-sim.y_crds_node[j_cell];
	       dz = sim.z_crds_node[k_cell+1]-sim.z_crds_node[k_cell];
	       break;
	     case sep::CYLINDRICAL:
	       dx = 0.5*(sim.x_crds_node[i_cell+1]*sim.x_crds_node[i_cell+1] - sim.x_crds_node[i_cell]*sim.x_crds_node[i_cell]);
	       dy = sim.y_crds_node[j_cell+1]-sim.y_crds_node[j_cell];
	       dz = sim.z_crds_node[k_cell+1]-sim.z_crds_node[k_cell];
	       break;
	     case sep::SPHERICAL:
	       // Transform SIM frame node coordinates to fixed frame coordinates:
	       X0 = sim.x_crds_node[i_cell  ] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[0];
	       X1 = sim.x_crds_node[i_cell+1] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[0];
	       dx = (X1*X1*X1 - X0*X0*X0)/3.0;

	       Y0 = sim.y_crds_node[j_cell  ] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[1];
	       Y1 = sim.y_crds_node[j_cell+1] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[1];
	       dy = cos(Y0) - cos(Y1);
	       
	       Z0 = sim.z_crds_node[k_cell  ] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[2];
	       Z1 = sim.z_crds_node[k_cell+1] + max(ZERO,sim.t-simControl.t_setup)*simControl.V_frame[2];
	       dz = Z1-Z0;
	       break;
	     default:
	       dx = sim.x_crds_node[i_cell+1]-sim.x_crds_node[i_cell];
	       dy = sim.y_crds_node[j_cell+1]-sim.y_crds_node[j_cell];
	       dz = sim.z_crds_node[k_cell+1]-sim.z_crds_node[k_cell];
	       break;
	    }
	    simControl.cellVolumes[counter] = dx*dy*dz;
	    ++counter;
	 }
      }

      simControl.timestepCellVolumesCalculated = sim.timestep;
   }

   bool prepareWaveEnergy(Simulation& sim,SimulationClasses& simClasses,const vector<pargrid::CellID>& blockLIDs) {
      return true;
      // Get pointer to (anti)parallel wave energy array:
      Real* waveEnergy = NULL;
      if (simControl.alfvenSign > 0) {
	 waveEnergy = simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      } else {
	 waveEnergy = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      }
      if (waveEnergy == NULL) return false;
      
      // Get pointer to temporary wave energy array:
      Real* temporaryWaveEnergy = simClasses.pargrid.getUserDataStatic<Real>(simControl.temporaryWaveEnergyDataID);
      if (temporaryWaveEnergy == NULL) return false;

      const uint32_t NWL        = simControl.N_wavelengthMeshCells;
      const uint32_t SIZE_BLOCK = block::SIZE*NWL;

      // Wave energies are not limited, just copy values to temporary array
      // and divide by spatial cell volume:
	 #warning OPTIMIZE ME, this copy is unnecessary
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    const Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	    
	    for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	       temporaryWaveEnergy[blockLID*SIZE_BLOCK+block::index(i,j,k)*NWL+l] 
		 = waveEnergy[blockLID*SIZE_BLOCK+block::index(i,j,k)*NWL+l] / cellVolume;
	    }
	 }
      }

      return true;
   }
   
   bool scaleWaveEnergy(Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         static int profWaveScalingTotal = -1;
         profile::start("Scale Wave Energy",profWaveScalingTotal);
      #endif
      
      // Scale wave energy by total particle phase-space density on boundary blocks:
      const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(sim.defaultStencilID);
      prepareWaveEnergy(sim,simClasses,boundaryBlocks);
      
      pargrid::DataID waveEnergyDataID = pargrid::INVALID_DATAID;
      if (simControl.alfvenSign < 0) waveEnergyDataID = simControl.antiparAlfvenWaveEnergyDataID;
      else waveEnergyDataID = simControl.parAlfvenWaveEnergyDataID;

      // Send scaled wave energy to neighbor processes:
      if (simClasses.pargrid.startNeighbourExchange(sim.defaultStencilID,waveEnergyDataID) == false) exit(1);
      if (simClasses.pargrid.startNeighbourExchange(sim.defaultStencilID,simControl.temporaryWaveEnergyDataID) == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to send temporary wave energy" << endl << write;
	 return exitWithStatus(false);
      }
      
      // Scale wave energy by total particle phase-space density on inner blocks:
      const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(sim.defaultStencilID);
      prepareWaveEnergy(sim,simClasses,innerBlocks);

      // Wait for spectral wave energy density:
      #if PROFILE_LEVEL > 0
         static int profWaveEnergyWait = -1;
         profile::start("MPI Wait (Scaled Wave Energy)",profWaveEnergyWait);
      #endif
      if (simClasses.pargrid.wait(sim.defaultStencilID,waveEnergyDataID) == false) exit(1);
      if (simClasses.pargrid.wait(sim.defaultStencilID,simControl.temporaryWaveEnergyDataID,"scaleWaveEnergy") == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to receive temporary wave energy" << endl << write;
	 #if PROFILE_LEVEL > 0
	    profile::stop(); // MPI Wait (Scaled Wave Energy)
	 #endif
	 return exitWithStatus(false); // Scale wave energy
      }
      #if PROFILE_LEVEL > 0
         profile::stop();  // profWaveEnergyWait
      #endif
      return exitWithStatus(success);
   }

   bool applySmoothing(Simulation& sim,SimulationClasses& simClasses,const std::vector<pargrid::CellID>& blockLIDs) { 
      // Get wave growth factor array:
      Real* RESTRICT waveGrowth = NULL;
      if (simControl.alfvenSign < 0) {
	 waveGrowth = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparWaveGrowthDataID);
      } else {
	 waveGrowth = simClasses.pargrid.getUserDataStatic<Real>(simControl.parWaveGrowthDataID);
      }
      if (waveGrowth == NULL) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to get array(s) in applySmoothing" << endl << write;
	 return false;
      }
      
      #if PROFILE_LEVEL > 0
         static int profCounter = -1;
         profile::start("Smoothing",profCounter);
      #endif

      const int32_t NWL = simControl.N_wavelengthMeshCells;
      Real* diffFlux = new Real[NWL];
      const Real diffConst = 0.25;
      const uint32_t SIZE = block::SIZE*simControl.N_wavelengthMeshCells;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];
	 for (int32_t cell=0; cell<block::SIZE; ++cell) {
	    // Clear diffusive flux:
	    for (int32_t l=0; l<NWL; ++l) diffFlux[l] = 0.0;
	    
	    // Find boundaries in wavelength space where diffusive smoothing is applied:
	    int32_t l_min = 0;
	    int32_t l_max = NWL-1;
	    while (l_min < NWL && waveGrowth[blockLID*SIZE+cell*NWL+l_min] == 0) ++l_min;
	    while (l_max > -1 && waveGrowth[blockLID*SIZE+cell*NWL+l_max] == 0) --l_max;
	    l_min = min(NWL-1,l_min);
	    l_max = max(0,l_max);

	    for (int i=0; i<1; ++i) {
	       // Calculate diffusive fluxes:
	       for (int32_t l=0; l<NWL-1; ++l) {
		  Real U_left = waveGrowth[blockLID*SIZE+cell*NWL+l  ] / simControl.wavelengthMeshCellSizes[l  ];
		  Real U_rght = waveGrowth[blockLID*SIZE+cell*NWL+l+1] / simControl.wavelengthMeshCellSizes[l+1];
		  Real avgCellSize = 0.5*(simControl.wavelengthMeshCellSizes[l  ]+simControl.wavelengthMeshCellSizes[l+1]);
		  diffFlux[l] = -diffConst * (U_rght-U_left) * avgCellSize;
	       }

	       // Prevent energy from flowing out of wavelength mesh:
	       for (int32_t l=0; l<l_min; ++l) diffFlux[l] = 0.0;
	       for (int32_t l=l_max; l<NWL; ++l) diffFlux[l] = 0.0;

	       // Apply diffusive fluxes to growth factors:
	       waveGrowth[blockLID*SIZE+cell*NWL+0] -= diffFlux[0];
	       for (int32_t l=1; l<NWL; ++l) {
		  waveGrowth[blockLID*SIZE+cell*NWL+l  ] -= (diffFlux[l]-diffFlux[l-1]);
	       }
	    }
	 }
      }
      delete [] diffFlux; diffFlux = NULL;

      return exitWithStatus(true);
   }
   
   bool applyWaveGrowth(Simulation& sim,SimulationClasses& simClasses,
			const std::vector<pargrid::CellID>& blockLIDs,vector<ParticleListBase*>& particleLists) {
      // Get wave growth factor array:
      const Real* RESTRICT waveGrowthGlobal = NULL;
      const Real* RESTRICT waveEnergyGlobal = NULL;
      if (simControl.alfvenSign < 0) {
	 waveGrowthGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparWaveGrowthDataID);
	 waveEnergyGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      } else {
	 waveGrowthGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.parWaveGrowthDataID);
	 waveEnergyGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      }
      if (waveGrowthGlobal == NULL) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to get wave growth factor array in applyWaveGrowth" << endl << write;
	 return false;
      }
      
      // Find the correct Alfven wave packet population (parallel or antiparallel):
      pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
      for (size_t s=0; s<particleLists.size(); ++s) {
	 if (particleLists[s]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
	 const LagrangianSpecies* species = reinterpret_cast<const LagrangianSpecies*>(particleLists[s]->getSpecies());
	 if (simControl.alfvenSign == species->propagationDirection) {
	    if (particleLists[s]->getParticles(speciesDataID) == false) {
	       simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to get Lagrangian particles" << endl << write;
	       return false;
	    }
	    break;
	 }
      }
      if (speciesDataID == pargrid::INVALID_DATAID) return false;

      typedef LagrangianParticle<Real> PARTICLE;
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses.pargrid.getUserDataDynamic<PARTICLE>(speciesDataID);
      if (wrapper.valid() == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Got invalid wrapper from ParGrid in applyWaveGrowth" << endl << write;
	 return false;
      }

      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      const uint32_t SIZE_BLOCK = block::SIZE_PLUS_ONE_LAYER*(NWL+2);
      Real* waveGrowthFactors = new Real[SIZE_BLOCK];
      
      Real* waveEnergy = new Real[SIZE_BLOCK];
      Real t_propag = 0.0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 // Measure block accumulation time if we are testing for repartitioning:
	 if (sim.countPropagTime == true) t_propag = MPI_Wtime();
	 
	 const pargrid::CellID blockLID = blockLIDs[b];
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 
	 // Calculate block's indices in global mesh:
	 uint32_t blockIndices[3];
	 block::calculateBlockIndices(sim,blockGID,blockIndices[0],blockIndices[1],blockIndices[2]);
	
	 // If interpolations are upwinded/downwinded near the shock, 
	 // classify cells according to their shock regions:
	 bool shockedBlock = false;
	 int shockRegions[block::WIDTH_X+2];
	 if (simControl.useShockUpwinding == true) {
	    shockedBlock = classifyShockedCells(sim.t,shockRegions,blockIndices);
	 }

	 // Fetch wave growth factors:
	 for (uint32_t i=0; i<SIZE_BLOCK; ++i) waveGrowthFactors[i] = 0.0;
	 loadScalar4D(&simClasses,blockLID,waveGrowthGlobal,waveGrowthFactors,NWL);
	 for (int32_t i=0; i<block::SIZE_PLUS_ONE_LAYER; ++i) waveGrowthFactors[i*(NWL+2)    +0] = 0.0;
	 for (int32_t i=0; i<block::SIZE_PLUS_ONE_LAYER; ++i) waveGrowthFactors[i*(NWL+2)+NWL+1] = 0.0;

	 for (uint32_t i=0; i<SIZE_BLOCK; ++i) waveEnergy[i] = 0.0;
	 loadScalar4D(&simClasses,blockLID,waveEnergyGlobal,waveEnergy,NWL);
	 for (int32_t i=0; i<block::SIZE_PLUS_ONE_LAYER; ++i) waveEnergy[i*(NWL+2)    +0] = 0.0;
	 for (int32_t i=0; i<block::SIZE_PLUS_ONE_LAYER; ++i) waveEnergy[i*(NWL+2)+NWL+1] = 0.0;

	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    // Skip wave packets that are outside wavelength mesh in lambda direction:
	    if (wrapper.data()[blockLID][p].state[lagr::LAMBDA] < 0.0) continue;
	    if (wrapper.data()[blockLID][p].state[lagr::LAMBDA] > simControl.N_wavelengthMeshCells) continue;
	    
	    // Interpolate wave growth factor to wave packet position:
	    Real pos[4];
	    pos[0] = wrapper.data()[blockLID][p].state[lagr::XCRD] - blockIndices[0]*block::WIDTH_X + 1;
	    pos[1] = wrapper.data()[blockLID][p].state[lagr::YCRD] - blockIndices[1]*block::WIDTH_Y + 1;
	    pos[2] = wrapper.data()[blockLID][p].state[lagr::ZCRD] - blockIndices[2]*block::WIDTH_Z + 1;
	    pos[3] = wrapper.data()[blockLID][p].state[lagr::LAMBDA] + 1;
	    

	    int32_t indices[4];
	    Real sf[12];
	    const Real energy = wrapper.data()[blockLID][p].state[lagr::ENERGY];
	    const Real EPSILON = 1e-10;

	    Real waveGrowthFactor = 0.0;
	    switch (simControl.order) {
	     case 0:
	       #ifndef NDEBUG
	          cerr << "ERROR: NGP wave growth not implementedin applyWaveGrowth" << endl;
	          exit(1);
	       #endif
	       waveGrowthFactor = interpolateScalarLogicalNGP_4D(waveGrowthFactors,pos,NWL+2);
	       break;
	     case 1:
	       getShapeFactorsCIC_4D(pos,indices,sf);

	       // If upwinding/downwinding is used near the shock, modify shape 
	       // factors so that particles scattering in downstream region do not 
	       // grow waves in upstream region and vice versa:
	       if (shockedBlock == true) {
		  calculateUpwindedWaveShapeFactors(shockRegions,indices,sf,
						    simControl.shock->getShockRegion(sim.t,wrapper.data()[blockLID][p].state));
	       }

	       for (int k_off=0; k_off<2; ++k_off) for (int j_off=0; j_off<2; ++j_off) for (int i_off=0; i_off<2; ++i_off) for (int l_off=0; l_off<2; ++l_off) {
		  const int32_t index = block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*(NWL+2)+indices[3]+l_off;
		  const Real shapeFactor = sf[0+i_off] * sf[2+j_off] * sf[4+k_off] * sf[6+l_off] * (energy / (waveEnergy[index]+EPSILON));
		  waveGrowthFactor += shapeFactor * waveGrowthFactors[index];
	       }
	       break;
	     case 2:
	       getShapeFactorsTSC_4D(pos,indices,sf);
	       
	       // If upwinding/downwinding is used near the shock, modify shape 
	       // factors so that particles scattering in downstream region do not 
	       // grow waves in upstream region and vice versa:
	       if (shockedBlock == true) {
		  calculateUpwindedWaveShapeFactors(shockRegions,indices,sf,
						    simControl.shock->getShockRegion(sim.t,wrapper.data()[blockLID][p].state));
	       }
	       
	       for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) for (int l_off=-1; l_off<2; ++l_off) {
		  const int32_t index = block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*(NWL+2)+indices[3]+l_off;
		  const Real shapeFactor = sf[1+i_off] * sf[4+j_off] * sf[7+k_off] * sf[10+l_off] * (energy / (waveEnergy[index]+EPSILON));
		  waveGrowthFactor += shapeFactor * waveGrowthFactors[index];
	       }
	       break;
	     default:
	       simClasses.logger << "(SEP PROPAGATE) ERROR: Unknown order of accuracy in applyWaveGrowth" << endl << write;
	       exit(1);
	       break;
	    }

	    // If we are at data save step, collect total wave energy change:
	    if (sim.atDataSaveStep == true) {
	       simControl.dU_waves += waveGrowthFactor;
	    }

	    // Apply wave growth:
	    Real minEnergy = simControl.maximumWaveEnergyLossPercentage * wrapper.data()[blockLID][p].state[lagr::ENERGY];
	    wrapper.data()[blockLID][p].state[lagr::ENERGY] += waveGrowthFactor;
	    wrapper.data()[blockLID][p].state[lagr::ENERGY] = max(minEnergy,wrapper.data()[blockLID][p].state[lagr::ENERGY]);
	 }

	 // Store block propagation time:
	 if (sim.countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses.pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }

      delete [] waveGrowthFactors; waveGrowthFactors = NULL;
      delete [] waveEnergy; waveEnergy = NULL;
      return true;
   }
   
   /** Calculate wave growth on all local blocks. Wave growth factors (=gamma) must 
    * be correct on all local blocks prior to calling this function. After 
    * this function returns correctly, wave energy has been updated on all local blocks.
    * 
    * @return If true, wave growth was calculated successfully.*/
   bool calculateWaveGrowth(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      if (simControl.alfvenSign < 0 && simControl.applyAntiparallelWaveGrowth == false) return success;
      if (simControl.alfvenSign > 0 && simControl.applyParallelWaveGrowth == false) return success;

      #if PROFILE_LEVEL > 0
         static int profWaveGrowthTotal = -1;
         profile::start("Wave Growth Total",profWaveGrowthTotal);
      #endif

      const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(sim.defaultStencilID);
      if (applySmoothing(sim,simClasses,boundaryBlocks) == false) {
	 success = false;
      }

      // Get correct (antiparallel or parallel) wave growth array DataID:
      pargrid::DataID waveGrowthDataID = pargrid::INVALID_DATAID;
      if (simControl.alfvenSign < 0) waveGrowthDataID = simControl.antiparWaveGrowthDataID;
      else waveGrowthDataID = simControl.parWaveGrowthDataID;

      // Send wave growth factors to neighbor processes:
      if (simClasses.pargrid.startNeighbourExchange(sim.defaultStencilID,waveGrowthDataID) == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to start copy of wave growth factor" << endl;
	 simClasses.logger << "\t direction = " << simControl.alfvenSign << " array DataID = ";
	 simClasses.logger << waveGrowthDataID << endl << write;
	 success = false;
      }

      // TEST
      pargrid::DataID waveEnergyID = simControl.parAlfvenWaveEnergyDataID;
      if (simControl.alfvenSign < 0) waveEnergyID = simControl.antiparAlfvenWaveEnergyDataID;
      simClasses.pargrid.startNeighbourExchange(sim.defaultStencilID,waveEnergyID);
      // END TEST

      const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(sim.defaultStencilID);

      // Calculate wave growth on inner blocks:
      if (applySmoothing(sim,simClasses,innerBlocks) == false) success = false;
      if (applyWaveGrowth(sim,simClasses,innerBlocks,particleLists) == false) success = false;

      // Wait for wave growth factors:
      #if PROFILE_LEVEL > 0
         static int profWaitWaveGrowth = -1;
         profile::start("MPI Wait (wave growth factor)",profWaitWaveGrowth);
      #endif
      if (simClasses.pargrid.wait(sim.defaultStencilID,waveGrowthDataID,"calculateWaveGrowth") == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Wait for wave growth factors failed" << endl << write;
	 success = false;
      }
      simClasses.pargrid.wait(sim.defaultStencilID,waveEnergyID,"calculateWaveGrowth");
      #if PROFILE_LEVEL > 0
         profile::stop();  // profWaitWaveGrowth
      #endif

      // Calculate wave growth on boundary blocks:
      if (applyWaveGrowth(sim,simClasses,boundaryBlocks,particleLists) == false) success = false;
      return exitWithStatus(success);
   }

   /** Scatter particles on all local blocks and accumulate wave growth factors. 
    * After this function returns, wave growth factors are OK on all local blocks.
    * 
    * @return If true, particles were scattered successfully.*/
   bool scatterParticlesAccumulateGrowthFactor(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;

      #if PROFILE_LEVEL > 0
         static int profScatterAccumulateTotal = -1;
         profile::start("Scattering and Wave Growth",profScatterAccumulateTotal);
      #endif

      // Clear wave growth array:
      if (clearWaveGrowthArray(sim,simClasses) == false) success = false;
      
      // Scatter particles on boundary blocks:
      const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(sim.inverseStencilID);
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 if (scatterParticles(boundaryBlocks,sim,simClasses,particleLists[p],simControl.fieldsGetState) == false) success = false;
      }

      // Get correct (antiparallel or parallel) wave growth array DataID:
      pargrid::DataID waveGrowthDataID = pargrid::INVALID_DATAID;
      if (simControl.alfvenSign < 0) waveGrowthDataID = simControl.antiparWaveGrowthDataID;
      else waveGrowthDataID = simControl.parWaveGrowthDataID;
      
      if (simControl.applyWaveGrowth == true) {
	 // Send wave growth updates to neighbor processes:
	 if (simClasses.pargrid.startNeighbourExchange(sim.inverseStencilID,waveGrowthDataID) == false) {
	    simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to send wave growth factors in scatterParticlesAccumulateGrowthFactor" << endl << write;
	    success = false;
	 }
      }

      // Scatter particles on inner blocks:
      const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(sim.inverseStencilID);
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (particleLists[p]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 if (scatterParticles(innerBlocks,sim,simClasses,particleLists[p],simControl.fieldsGetState) == false) success = false;
      }

      if (simControl.applyWaveGrowth == true) {
	 // Wait for growth factors:
         #if PROFILE_LEVEL > 0
            static int profMPIWait = -1;
            profile::start("MPI Wait (wave growth updates)",profMPIWait);
         #endif
	 
	 if (simClasses.pargrid.wait(sim.inverseStencilID,waveGrowthDataID,"scatterParticlesAccumulateGrowthFactor") == false) {
	    simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to receive wave growth updates in ";
	    simClasses.logger << "scatterParticlesAccumulateGrowthFactor" << endl << write;
	    success = false;
	 }
         #if PROFILE_LEVEL > 0
            profile::stop();  // profMPIWait
         #endif
      
	 // Add remote updates to growth factors:
	 Real* waveGrowth = simClasses.pargrid.getUserDataStatic<Real>(waveGrowthDataID);
	 if (waveGrowth == NULL) {
	    simClasses.logger << "(SEP PROPAGATE) ERROR: Wave growth array is NULL in ";
	    simClasses.logger << "scatterParticlesAccumulateGrowthFactor" << endl << write;
	    success = false;
	 }
	 unsigned int* offsets = NULL;
	 Real* remoteUpdateArray = NULL;
	 if (simClasses.pargrid.getRemoteUpdates<Real>(sim.inverseStencilID,waveGrowthDataID,
						       offsets,remoteUpdateArray) == false) {
	    simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to get remote wave growth updates" << endl << write;
	    success = false;
	 }
	 if (success == true) {
	    const uint32_t NWL = simControl.N_wavelengthMeshCells;
	    const uint32_t SIZE_BLOCK = block::SIZE*NWL;
	    for (size_t b=0; b<boundaryBlocks.size(); ++b) {
	       const pargrid::CellID blockLID = boundaryBlocks[b];
	       for (unsigned int offset=offsets[b]; offset<offsets[b+1]; ++offset) {
		  for (uint32_t l=0; l<SIZE_BLOCK; ++l) waveGrowth[blockLID*SIZE_BLOCK+l] += remoteUpdateArray[offset*SIZE_BLOCK+l];
	       }
	    }
	 }
      }  // if (simControl.applyWaveGrowth == true)

      return exitWithStatus(success); // Scattering and Wave Growth
   }
   
   bool scatterStep(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         static int profScatterTotal = -1;
         profile::start("(SEP) Scatter Step",profScatterTotal);
      #endif
      
      // Accumulate Alfven wave energy:
      if (accumulateWaveEnergies(sim,simClasses,particleLists) == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to accumulate Alfven wave energies" << endl << write;
	 success = false;
      }

      if (sim.t < simControl.t_setup) return exitWithStatus(success);

      // Exit if particles do not scatter off the chosen wave mode:
      if (simControl.alfvenSign < 0 && simControl.scatterAntiparallel == false) return exitWithStatus(success);
      if (simControl.alfvenSign > 0 && simControl.scatterParallel == false) return exitWithStatus(success);

      // Accumulate particle 4D phase-space density on antiparallel or parallel Alfven
      // wave rest frame, as selected by the sign of simControl.alfvenSign variable. After 
      // this function returns, phase space weight is correct on all local blocks.
      /*if (sim.atDataSaveStep == true) {
	 if (accumulateParticlePhaseSpaceWeight(sim,simClasses,particleLists) == false) {
	    simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to accumulate particle phase-space weight" << endl << write;
	    return exitWithStatus(false);
	 }
      }*/
	 
      // Scale wave energy by particle phase-space weight:
      if (scaleWaveEnergy(sim,simClasses) == false) success = false;
	 
      // Scatter particles and accumulate wave growth factors:
      if (scatterParticlesAccumulateGrowthFactor(sim,simClasses,particleLists) == false) success = false;
	 
      // Calculate wave growth on all local blocks:
      if (simControl.applyWaveGrowth == true) {
	 if (calculateWaveGrowth(sim,simClasses,particleLists) == false) success = false;
	 
	 if (sim.atDataSaveStep == true) {
	    accumulateWaveEnergies(sim,simClasses,particleLists);
	 }
      }

      return exitWithStatus(success);
   }
   
   bool propagate(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
      bool success = true;
      // If partitioning has changed, recalculate spatial cell volumes:
      recalculateCellVolumes(sim,simClasses);

      // Set SIM frame speed to zero if setup time has passes:
      if (sim.t >= simControl.t_setup) {
	 for (int i=0; i<3; ++i) simControl.V_frame[i] = 0.0;
      }

      // Clear wave and particle energy change counters:
      simControl.dU_particles = 0.0;
      simControl.dU_waves = 0.0;

      // Reduce data save interval while simulation is being set up.
      static Real dataSaveInterval = NAN;
      static bool dataSaveIntervalStored = false;
      if (dataSaveIntervalStored == false && sim.dataIntervalIsTime == true) {
	 dataSaveInterval = sim.dataIntervalFloat;
	 dataSaveIntervalStored = true;
      }

      if (sim.dataIntervalIsTime) {
	 if (sim.t >= simControl.t_setup) {
	    sim.dataIntervalFloat = dataSaveInterval;
	 } else {
	    sim.dataIntervalFloat = 100.0 * sim.dt;
	 }
      }
      
      // Re-evaluate timestep dt:
      static bool timestepRecalculatedAfterSetup = false;
      if (simControl.timestepRecalculateInterval > 0) {
	 if (sim.timestep % simControl.timestepRecalculateInterval == 0) {
	    #ifndef NDEBUG
	       if (sim.mpiRank == sim.MASTER_RANK) {
		  simClasses.logger << "(SEP PROPAGATE) Recalculating dt" << endl << write;
	       }
	    #endif
	    recalculateTimestep(sim,simClasses,particleLists);
	 }
	 if (timestepRecalculatedAfterSetup == false) {
	    if (sim.t >= simControl.t_setup) {
	       #ifndef NDEBUG
	          if (sim.mpiRank == sim.MASTER_RANK) {
		     simClasses.logger << "(SEP PROPAGATE) Recalculating dt as t > t_setup" << endl << write;
		  }
               #endif
	       recalculateTimestep(sim,simClasses,particleLists);
	       timestepRecalculatedAfterSetup = true;
	    }
	 }
      }

      // Scatter particles using both wave modes, accumulations
      // are also done in scatterStep function:
      #ifndef NDEBUG
         if (sim.mpiRank == sim.MASTER_RANK) {
	    simClasses.logger << "(SEP PROPAGATE) Starting scatter step" << endl << write;
	 }
         checkCoordinates(sim,simClasses,particleLists,"scatter step");
      #endif

      simControl.currentMaxScatteringSubsteps = 1;

      // Master process randomizes which scattering (parallel or antiparallel waves) 
      // is solved first and distributes it to all processes:
      if (sim.mpiRank == sim.MASTER_RANK) {
	 if (simClasses.random.uniform() < 0.5) {
	    simControl.alfvenSign = -1;
	 } else {
	    simControl.alfvenSign = +1;
	 }
      }
      MPI_Bcast(&simControl.alfvenSign,1,MPI_Type<int32_t>(),sim.MASTER_RANK,sim.comm);
      for (int i=0; i<2; ++i) {
	 // Scatter particles using selected wave mode:
	 if (simControl.alfvenSign < 0 && simControl.includeAntiparWaves == true) 
	   if (scatterStep(sim,simClasses,particleLists) == false) success = false;
	 
	 if (simControl.alfvenSign > 0 && simControl.includeParWaves == true)
	   if (scatterStep(sim,simClasses,particleLists) == false) success = false;
	 
	 // Switch wave mode:
	 simControl.alfvenSign *= -1;
      }

      // Propagate all wave packets and particle species:
      #ifndef NDEBUG
         if (sim.mpiRank == sim.MASTER_RANK) {
	    simClasses.logger << "(SEP PROPAGATE) Starting propagation step" << endl << write;
	 }
         checkCoordinates(sim,simClasses,particleLists,"propagation");
      #endif
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (simControl.propagateAlfvenWaves == false) 
	   if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 if (particleLists[p]->propagateBoundaryCellParticles() == false) success = false;
      }      
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (simControl.propagateAlfvenWaves == false)
	   if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 if (particleLists[p]->propagateInnerCellParticles() == false) success = false;
      }
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (simControl.propagateAlfvenWaves == false)
	   if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 if (particleLists[p]->waitParticleSends() == false) success = false;
      }

      // Apply boundary conditions:
      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (simControl.propagateAlfvenWaves == false)
	   if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 if (particleLists[p]->applyBoundaryConditions() == false) success = false;
      }

      // Inject new particles and wave packets. Note that 
      // all existing particles have been propagated to t+dt 
      // at this point:
      #ifndef NDEBUG
      if (sim.mpiRank == sim.MASTER_RANK) {
	 simClasses.logger << "(SEP PROPAGATE) Injecting new particles" << endl << write;
      }
      checkCoordinates(sim,simClasses,particleLists,"before injection");
      #endif

      for (size_t p=0; p<particleLists.size(); ++p) {
	 if (simControl.propagateAlfvenWaves == false)
	   if (particleLists[p]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
	 if (particleLists[p]->injectParticles() == false) success = false;
      }

      // Inject "extra" wave packets (possible reflected waves) and 
      // particles (some injectors produce these) to simulation:
      if (injectExtraParticles() == false) {
	 simClasses.logger << "(SEP PROPAGATE) ERROR: Failed to inject extra particles" << endl << write;
	 success = false;
      }

      #ifndef NDEBUG
      checkCoordinates(sim,simClasses,particleLists,"after injection");
      #endif

      // Check for particle splitting:
      if (simControl.particleSplitter != NULL) {
	 if (simControl.particleSplitter->split() == false) success = false;
	 #ifndef NDEBUG
	    checkCoordinates(sim,simClasses,particleLists,"after splitting");
         #endif
      }
      
      return success;
   }
   
   void recalculateTimestep(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists) {
      #if PROFILE_LEVEL > 0
         static int profTimestep = -1;
         profile::start("Timestep Recalculation",profTimestep);
      #endif

      // Calculate maximum timestep for Alfven waves:
      Real dt_max_waves = numeric_limits<Real>::max();
      Real centroid[3];
      Real B[3];
      Real V_wave_minus[3], V_wave_plus[3];
      Real dV_wave;
      for (pargrid::CellID blockLID=0; blockLID<simClasses.pargrid.getNumberOfLocalCells(); ++blockLID) {
	 // Calculate block i,j,k indices:
	 uint32_t i_block,j_block,k_block;
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);
	 
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    const int64_t i_cell = i_block*block::WIDTH_X + i;
	    const int64_t j_cell = j_block*block::WIDTH_Y + j;
	    const int64_t k_cell = k_block*block::WIDTH_Z + k;
	    
	    // Calculate smallest cell sizes around cell (i,j,k). This ensures that 
	    // Courant condition is evaluated correctly for non-uniform meshes:
	    Real dx_min = sim.dx_cell[i_cell];
	    Real dy_min = sim.dy_cell[j_cell];
	    Real dz_min = sim.dz_cell[k_cell];
	    if (i_cell-1 >= 0) dx_min = min(dx_min,sim.dx_cell[i_cell-1]);
	    if (j_cell-1 >= 0) dy_min = min(dy_min,sim.dy_cell[j_cell-1]);
	    if (k_cell-1 >= 0) dz_min = min(dz_min,sim.dz_cell[k_cell-1]);
	    if (i_cell+1 < sim.x_blocks*block::WIDTH_X) dx_min = min(dx_min,sim.dx_cell[i_cell+1]);
	    if (j_cell+1 < sim.y_blocks*block::WIDTH_Y) dy_min = min(dy_min,sim.dy_cell[j_cell+1]);
	    if (k_cell+1 < sim.z_blocks*block::WIDTH_Z) dz_min = min(dz_min,sim.dz_cell[k_cell+1]);
	    
	    // Calculate cell centroid position in logical coordinates:
	    centroid[0] = i_cell + 0.5;
	    centroid[1] = j_cell + 0.5;
	    centroid[2] = k_cell + 0.5;

	    const Real R     = sim.x_crds_node[i_cell] + 0.5*sim.dx_cell[i_cell];
	    const Real THETA = sim.y_crds_node[j_cell] + 0.5*sim.dy_cell[j_cell];

	    // Evaluate Alfven speed in cell centroid:
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim.t,centroid,B,V_wave_minus,dV_wave,V_alfven,-1);
	    (*simControl.fieldsGetState)(blockLID,sim.t,centroid,B,V_wave_plus,dV_wave,V_alfven,+1);
	    Real vx_max = max(fabs(V_wave_minus[0]),fabs(V_wave_plus[0])) + simControl.EPSILON;
	    Real vy_max = max(fabs(V_wave_minus[1]),fabs(V_wave_plus[1])) + simControl.EPSILON;
	    Real vz_max = max(fabs(V_wave_minus[2]),fabs(V_wave_plus[2])) + simControl.EPSILON;
	    
	    switch (simControl.coordinateSystem) {
	     case sep::UNKNOWN:
	       exit(1);
	       break;
	     case sep::CARTESIAN:
	       dt_max_waves = min(dt_max_waves, dx_min / vx_max);
	       dt_max_waves = min(dt_max_waves, dy_min / vy_max);
	       dt_max_waves = min(dt_max_waves, dz_min / vz_max);
	       break;
	     case sep::CYLINDRICAL:
	       dt_max_waves = min(dt_max_waves, dx_min / vx_max);
	       dt_max_waves = min(dt_max_waves, dy_min * R / vy_max);
	       dt_max_waves = min(dt_max_waves, dz_min / vz_max);
	       break;
	     case sep::SPHERICAL:
	       dt_max_waves = min(dt_max_waves, dx_min / vx_max);
	       dt_max_waves = min(dt_max_waves, dy_min * R / vy_max);
	       dt_max_waves = min(dt_max_waves, dz_min * R*sin(THETA) / vz_max);
	       break;
	     default:
	       exit(1);
	       break;
	    }
	 }
      }

      // Calculate maximum timestep for particles:
      Real dt_max_particles = numeric_limits<Real>::max();
      for (size_t species=0; species<particleLists.size(); ++species) {
	 // Skip Lagrangian wave packets as their dt is evaluated differently:
	 if (particleLists[species]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;

	 // Get maximum injection speed:
	 Real v_injector_max = particleLists[species]->getInjector()->getMaximumSpeed();
	 
	 pargrid::DataID dataID = pargrid::INVALID_DATAID;
	 if (particleLists[species]->getParticles(dataID) == false) continue;
	 if (dataID == pargrid::INVALID_DATAID) continue;
	 
	 typedef sep::Particle<Real> PARTICLE;
	 const pargrid::DataWrapper<PARTICLE>& wrapper = simClasses.pargrid.getUserDataDynamic<PARTICLE>(dataID);
	 if (wrapper.valid() == false) continue;

	 for (pargrid::CellID blockLID=0; blockLID<simClasses.pargrid.getNumberOfLocalCells(); ++blockLID) {
	    //if (wrapper.size(blockLID) == 0) continue;
	    
	    // Calculate block i,j,k indices:
	    uint32_t i_block,j_block,k_block;
	    const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	    block::calculateBlockIndices(sim,blockGID,i_block,j_block,k_block);
	    
	    Real dx_min = numeric_limits<Real>::max();
	    Real dy_min = numeric_limits<Real>::max();
	    Real dz_min = numeric_limits<Real>::max();
	    for (int k=-1; k<block::WIDTH_Z+1; ++k) for (int j=-1; j<block::WIDTH_Y+1; ++j) for (int i=-1; i<block::WIDTH_X+1; ++i) {
	       const int64_t i_cell = i_block*block::WIDTH_X + i;
	       const int64_t j_cell = j_block*block::WIDTH_Y + j;
	       const int64_t k_cell = k_block*block::WIDTH_Z + k;
	       if (i_cell >= 0 && i_cell < sim.x_blocks*block::WIDTH_X) dx_min = min(dx_min,sim.dx_cell[i_cell]);
	       if (j_cell >= 0 && j_cell < sim.y_blocks*block::WIDTH_Y) dy_min = min(dy_min,sim.dy_cell[j_cell]);
	       if (k_cell >= 0 && k_cell < sim.z_blocks*block::WIDTH_Z) dz_min = min(dz_min,sim.dz_cell[k_cell]);
	    }
	    
	    Real pos[3] = {0,0,0};
	    Real v_par_max = 0.0;
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       if (fabs(wrapper.data()[blockLID][p].state[particle::V_PAR]) > v_par_max) {
		  v_par_max = fabs(wrapper.data()[blockLID][p].state[particle::V_PAR]);
		  pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD];
		  pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD];
		  pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD];
	       }
	    }
	    v_par_max = max(v_par_max,v_injector_max);
	    v_par_max += simControl.EPSILON;

	    Real R = 0.0;
	    Real THETA = 0.0;
	    const uint32_t i_max = static_cast<uint32_t>(pos[0]);
	    const uint32_t j_max = static_cast<uint32_t>(pos[1]);
	    
	    switch (simControl.coordinateSystem) {
	     case sep::UNKNOWN:
	       exit(1);
	       break;
	     case sep::CARTESIAN:
	       dt_max_particles = min(dt_max_particles,dx_min / v_par_max);
	       dt_max_particles = min(dt_max_particles,dy_min / v_par_max);
	       dt_max_particles = min(dt_max_particles,dz_min / v_par_max);
	       break;
	     case sep::CYLINDRICAL:
	       R = sim.x_crds_node[i_max] + (pos[0]-i_max)*sim.dx_cell[i_max];
	       dt_max_particles = min(dt_max_particles,dx_min / v_par_max);
	       dt_max_particles = min(dt_max_particles,dy_min*R / v_par_max);
	       dt_max_particles = min(dt_max_particles,dz_min / v_par_max);
	       break;
	     case sep::SPHERICAL:
	       R = sim.x_crds_node[i_max] + (pos[0]-i_max)*sim.dx_cell[i_max];
	       THETA = sim.y_crds_node[j_max] + (pos[1]-j_max)*sim.dy_cell[j_max];
	       
	       dt_max_particles = min(dt_max_particles,dx_min / v_par_max);
	       dt_max_particles = min(dt_max_particles,dy_min*R / v_par_max);
	       dt_max_particles = min(dt_max_particles,dz_min*R*sin(THETA) / v_par_max);
	       break;
	     default:
	       exit(1);
	       break;
	    }
	 }
      }
      Real dt_max = min(dt_max_waves,dt_max_particles*simControl.dt_particle_coeff);
      
      // Reduce timestep further if we're taking too many scattering substeps:
      if (simControl.currentMaxScatteringSubsteps >= 10) dt_max = min(dt_max, sim.dt / 2);

      // Do not increase dt too much so that scattering substeps does jump up:
      if (simControl.currentMaxScatteringSubsteps > 1) dt_max = min(dt_max, 2*sim.dt);

      // Reduce maximum timestep to master process:
      Real dt_max_global;
      MPI_Reduce(&dt_max,&dt_max_global,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
      
      // Reduce max scatter substeps to master process:
      uint32_t maxScatterSubsteps;
      MPI_Reduce(&simControl.currentMaxScatteringSubsteps,&maxScatterSubsteps,1,MPI_Type<uint32_t>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
      
      // Use Courant number 0.5:
      if (simControl.dt_maximum > 0.0) {
	 sim.dt = min(simControl.dt_maximum,simControl.courantNumber*dt_max_global);
      } else {
	 sim.dt = simControl.courantNumber*dt_max_global;
      }
      MPI_Bcast(&sim.dt,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);

      simClasses.logger << "(SEP PROPAGATE) Timestep adjusted at step " << sim.timestep << " time " << sim.t << " new dt = " << sim.dt << endl;
      simClasses.logger << "\t dt_max(waves): " << dt_max_waves << " dt_max(particles): " << dt_max_particles;
      simClasses.logger << " scatter substeps: " << maxScatterSubsteps << endl;
      simClasses.logger << write;
      
      #if PROFILE_LEVEL > 0
         profile::stop(); // Timestep Recalculation
      #endif
   }
   
} // namespace sep
