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

#ifndef SEP_LAGR_INJECTOR_DEFAULT_H
#define SEP_LAGR_INJECTOR_DEFAULT_H

#include <cstdlib>
#include <climits>
#include <stdint.h>

#include <main.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include <linear_algebra.h>
#include <object_factory_generic.h>

#include "sep_object_wrapper.h"
#include "sep_simcontrol.h"
#include "sep_lagr_definition.h"
#include "sep_fields_container.h"
#include "sep_distrib_wave_energy_base_class.h"

namespace sep {

   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class LagrangianInjectorDefault: public ParticleInjectorBase {
    public: 
      LagrangianInjectorDefault();

      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);
   
    private:
      bool inflowInjection;
      bool initialInjection;
      bool initialized;
      bool initialInjectionDone;
      uint32_t N_particlesPerCell;        /**< Number of macroparticles injected to each cell.*/
      uint32_t resolution;
      uint32_t resolutionLambda;          /**< Number of injected particles per lambda cell.*/
      SPECIES species;
      WaveEnergySpectrumBaseClass* waveEnergySpectrum;
      
      const Real DEF_VALUE;
      
      void injectParticles(pargrid::CellID block,Real t,Real B_mag,Real volume,Real* position,
			   unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* LIDefaultMaker() {return new LagrangianInjectorDefault<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   LagrangianInjectorDefault<SPECIES,PARTICLE>::LagrangianInjectorDefault(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity())  {
      waveEnergySpectrum = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorDefault<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".macroparticles_per_cell","Number of macroparticles injected per cell (integer)",(uint32_t)0);
      cr.add(PREFIX+".initial_injection","If 'yes' wave packets are injected to entire simulation domain before first time step (string)",std::string("yes"));
      cr.add(PREFIX+".inflow_boundary_injection","If 'yes' wave packets are injected to inflow boundaries (string)",std::string("yes"));
      cr.add(PREFIX+".energy_spectrum","Name of wave energy spectrum function (string)",std::string(""));
      cr.add(PREFIX+".energy_spectrum_parameters","Name of config file region containing parameters for wave energy spectrum function (string)",std::string(""));
      cr.add(PREFIX+".resolution_lambda","Number of injected particles per lambda cell, defaults to 1 (int)",(uint)1);
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorDefault<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      initialInjectionDone = false;
      if (waveEnergySpectrum != NULL) {
	 waveEnergySpectrum->finalize();
	 delete waveEnergySpectrum;
      }
      waveEnergySpectrum = NULL;
      return true;
   }

   /** Called by ParticleListSkeleton.*/
   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorDefault<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      if (waveEnergySpectrum == false) return false;
      bool success = true;
      
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      
      // Inject Lagrangian wave packets to whole simulation domain 
      // before simulation starts (if necessary):
      Real t_propag = 0.0;
      if (initialInjection == true && initialInjectionDone == false) {
	 for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	    // Measure block injection time if we are testing for repartitioning:
	    if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	    
	    // Calculate block indices:
	    const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	    uint32_t i_block,j_block,k_block;
	    block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	    
	    for (int32_t k=0; k<block::WIDTH_Z; ++k) for (int32_t j=0; j<block::WIDTH_Y; ++j) for (int32_t i=0; i<block::WIDTH_X; ++i) {
	       Real position[3];
	       position[0] = i_block*block::WIDTH_X + i + 0.5;
	       position[1] = j_block*block::WIDTH_Y + j + 0.5;
	       position[2] = k_block*block::WIDTH_Z + k + 0.5;
	       const Real volume = simControl.cellVolumes[blockLID*block::SIZE + block::index(i,j,k)];
	       
	       Real E[3];
	       Real B[3];
	       Real gradB[9];
	       (*simControl.fieldsGetFields)(blockLID,sim->t,position,E,B,gradB);
	       const Real B_mag = vectorMagnitude<3>(B);
	       
	       injectParticles(blockLID,sim->t,B_mag,volume,position,N_particles,wrapper);
	    }
	    
	    // Store block injection time:
	    if (sim->countPropagTime == true) {
	       t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	       simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	    }
	 }
      }
      initialInjectionDone = true;

      // Inject new waves to inflow boundaries:
      if (inflowInjection == false) return success;
      const std::vector<pargrid::CellID>& exteriorBlocks = simClasses->pargrid.getExteriorCells();
      for (size_t b=0; b<exteriorBlocks.size(); ++b) {
	 // Measure block injection time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	 
	 const pargrid::CellID blockLID = exteriorBlocks[b];
	 const uint32_t existingNbrFlags = simClasses->pargrid.getNeighbourFlags()[blockLID];
	 const uint32_t missingNbrFlags = (existingNbrFlags ^ pargrid::ALL_NEIGHBOURS_EXIST);

	 // Calculate block indices:
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 
	 // Count number of waves in each cell:
	 uint32_t N_waves[block::WIDTH_X*block::WIDTH_Y*block::WIDTH_Z];
	 for (uint32_t n=0; n<block::WIDTH_X*block::WIDTH_Y*block::WIDTH_Z; ++n) N_waves[n] = 0;
	 PARTICLE* waves = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    const int32_t i_index = static_cast<int32_t>(waves[p].state[0]) - i_block*block::WIDTH_X;
	    const int32_t j_index = static_cast<int32_t>(waves[p].state[1]) - j_block*block::WIDTH_Y;
	    const int32_t k_index = static_cast<int32_t>(waves[p].state[2]) - k_block*block::WIDTH_Z;
	    ++N_waves[k_index*block::WIDTH_Y*block::WIDTH_X+j_index*block::WIDTH_X+i_index];
	 }
	 
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    // Check if cell has enough wave particles in it:
	    if (N_waves[k*block::WIDTH_Y*block::WIDTH_X+j*block::WIDTH_X+i] >= resolution*resolution*resolution) continue;
	    
	    // Form a unit vector pointing to simulation domain:
	    Real normal[3];
	    for (int iii=0; iii<3; ++iii) normal[iii] = 0.0;
	    if (i == 0)                if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(-1,+0,+0))) > 0) normal[0] = +1.0;
	    if (i == block::WIDTH_X-1) if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+1,+0,+0))) > 0) normal[0] = -1.0;
	    if (j == 0)                if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,-1,+0))) > 0) normal[1] = +1.0;
	    if (j == block::WIDTH_Y-1) if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+1,+0))) > 0) normal[1] = -1.0;
	    if (k == 0)                if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+0,-1))) > 0) normal[2] = +1.0;
	    if (k == block::WIDTH_Z-1) if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+0,+1))) > 0) normal[2] = -1.0;
	    unitVector<3>(normal);
	    
	    // Check if particles should be injected on this cell:
	    if (vectorMagnitude<3>(normal) < 0.01) continue;
	    
	    const uint32_t i_cell = i_block*block::WIDTH_X + i;
	    const uint32_t j_cell = j_block*block::WIDTH_Y + j;
	    const uint32_t k_cell = k_block*block::WIDTH_Z + k;
	    
	    // Get wave velocity at cell centroid:
	    Real centroid[3];
	    centroid[0] = i_cell;
	    centroid[1] = j_cell;
	    centroid[2] = k_cell;
	    const Real volume = simControl.cellVolumes[blockLID*block::SIZE + block::index(i,j,k)];
	    
	    Real B[3];
	    Real V_wave[3];
	    Real dV_wave;
	    #warning Injector assumes simulation time sim->t
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,centroid,B,V_wave,dV_wave,V_alfven,species.propagationDirection);
	    const Real B_mag = vectorMagnitude<3>(B);
	    
	    // If this cell is an inflow cell, attempt to inject new wave packets:
	    if (dotProduct<3>(V_wave,normal) > 0.0) {
	       // Wave packets have been propagated to t+dt prior to calling injector,
	       // but sim->t has not yet been modified. That's why we need to use sim->t + sim->dt below:
	       Real offset[3];
	       offset[0] = 0.5 + V_wave[0]/sim->dx_cell[i_block*block::WIDTH_X+i] * (sim->t + sim->dt);
	       offset[1] = 0.5 + V_wave[1]/sim->dy_cell[i_block*block::WIDTH_Y+j] * (sim->t + sim->dt);
	       offset[2] = 0.5 + V_wave[2]/sim->dz_cell[i_block*block::WIDTH_Z+k] * (sim->t + sim->dt);
	       offset[0] -= static_cast<int32_t>(offset[0]);
	       offset[1] -= static_cast<int32_t>(offset[1]);
	       offset[2] -= static_cast<int32_t>(offset[2]);

	       // Integer cast above can return too small integer due to roundoff errors,
	       // i.e., if offset[0] = 9.999999999 we should get value 10 instead of 9.
	       // Here we correct these errors:
	       const Real EPS = 1.0e-8;
	       if (offset[0] > 1-EPS) offset[0] = 0.0;
	       if (offset[1] > 1-EPS) offset[1] = 0.0;
	       if (offset[2] > 1-EPS) offset[2] = 0.0;
	       
	       centroid[0] += offset[0];
	       centroid[1] += offset[1];
	       centroid[2] += offset[2];
	       
	       // If injection coordinate is within epsilon of cell node coordinates,
	       // i.e., x >= i_cell-EPS, clamp coordinates to lie within injection cell.
	       // Otherwise skip injection for this cell:
	       centroid[0] = std::max(1.0*i_cell,centroid[0]);
	       centroid[1] = std::max(1.0*j_cell,centroid[1]);
	       centroid[2] = std::max(1.0*k_cell,centroid[2]);
	       
	       #warning Antiparallel waves not injected to downstream region of trailing shock
	       if (i_block == 0) injectParticles(blockLID,sim->t,B_mag,volume,centroid,N_particles,wrapper);
	    }
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
   void LagrangianInjectorDefault<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,Real t_injection,Real B0_mag,
								     Real volume,Real* position, unsigned int* N_particles,
								     pargrid::DataWrapper<PARTICLE>& wrapper) {
      // Resize block's particle list to have space for new particles:
      const size_t N_injected = simControl.N_wavelengthMeshCells*block::SIZE*resolutionLambda;
      const pargrid::ArraySizetype oldSize = wrapper.size()[blockLID];
      N_particles[blockLID] += N_injected;
      wrapper.resize(blockLID,oldSize+N_injected);

      unsigned int counter = oldSize;
      PARTICLE* particles = wrapper.data()[blockLID];

      // Inject particles to cell:
      for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	 for (uint32_t p=0; p<resolutionLambda; ++p) {
	    const Real l_logical  = l + (1+2*p)*0.5/resolutionLambda;
	 
	    // Calculate injection position in logical coordinates:
	    particles[counter].state[lagr::XCRD]   = position[0];
	    particles[counter].state[lagr::YCRD]   = position[1];
	    particles[counter].state[lagr::ZCRD]   = position[2];
	    //particles[counter].state[lagr::LAMBDA] = l + 0.5;
	    particles[counter].state[lagr::LAMBDA] = l_logical;
	 
	    // This wave packet represents waves with logical lambda in interval [l,l+1[. 
	    // Energy is the total energy in waves in that interval, which we get from 
	    // wave energy spectrum function. We need to convert logical lambda coordinates 
	    // into physical ones, however:
	    //const Real minLambda = simControl.wavelengthMeshNodeCoordinates[l  ];
	    //const Real maxLambda = simControl.wavelengthMeshNodeCoordinates[l+1];
	    //Real l_avg = (*simControl.getLogicalWavelength)(0.5*(minLambda+maxLambda));
	    //if (fabs(l_avg - simControl.N_wavelengthMeshCells/2) < 1e-3) l_avg = l + 0.5;
	  
	    Real l_min = (*simControl.getPhysicalWavelength)(l_logical - 0.5/resolutionLambda);
	    Real l_max = (*simControl.getPhysicalWavelength)(l_logical + 0.5/resolutionLambda);
	    if (l == simControl.N_wavelengthMeshCells/2-1) l_max = 0.0;
	    if (l == simControl.N_wavelengthMeshCells/2  ) l_min = 0.0;
	    
	    Real l_logical_avg = (*simControl.getLogicalWavelength)(0.5*(l_min+l_max));
	    if (fabs(l_logical_avg - simControl.N_wavelengthMeshCells/2) < 1e-3) {
	       l_logical_avg = l_logical;
	    }
	    
	    const Real minLambda = l_min;
	    const Real maxLambda = l_max;	    

	    particles[counter].state[lagr::LAMBDA] = l_logical_avg;

	 #ifndef NDEBUG
	 if (fabs(l_logical_avg - simControl.N_wavelengthMeshCells/2) < 0.001) {
	    std::cerr << "(SEP LAGR INJ DEFAULT) ERROR: Invalid injection l " << l_logical_avg << std::endl;
	    std::cerr << "\t Trial value " << l + 0.5 << " min,max " << minLambda << ' ' << maxLambda << std::endl;
	    exit(1);
	 }
	 #endif
	 
	    particles[counter].state[lagr::ENERGY] 
	      = waveEnergySpectrum->getEnergyDensity(B0_mag,minLambda,maxLambda)*volume;
	 
	    ++counter;
	 }
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorDefault<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,
								ConfigReader& cr,const std::string& regionName,
								const ParticleListBase* plist) {
      simClasses.logger << "(SEP LAGR INJ DEFAULT) Starting to init." << std::endl << write;
      
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());

      if (addConfigFileItems(cr,regionName) == false) initialized = false;
      
      // Parse config file options:
      std::string initialInjectionString,inflowInjectionString;
      std::string energySpectrumName,energySpectrumParameters;
      cr.parse();
      if (cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell) == false) initialized = false;
      if (cr.get(regionName+".initial_injection",initialInjectionString) == false) initialized = false;
      if (cr.get(regionName+".inflow_boundary_injection",inflowInjectionString) == false) initialized = false;
      if (cr.get(regionName+".energy_spectrum",energySpectrumName) == false) initialized = false;
      if (cr.get(regionName+".energy_spectrum_parameters",energySpectrumParameters) == false) initialized = false;
      if (cr.get(regionName+".resolution_lambda",resolutionLambda) == false) initialized = false;
      
      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP LAGR INJ DEFAULT) ERROR: Parameter '" << regionName + 
	   ".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
      }
      
      // Check if initial state injection should be done:
      if (initialInjectionString == "yes") initialInjection = true;
      else initialInjection = false;

      // Check if wave packets are injected at inflow boundaries:
      if (inflowInjectionString == "yes") inflowInjection = true;
      else inflowInjection = false;
      
      // Get energy spectrum function and initialize it:
      waveEnergySpectrum = sep::getObjectWrapper().waveEnergySpectrumFactory.create(energySpectrumName);
      if (waveEnergySpectrum == NULL) {
	 simClasses.logger << "(SEP LAGR INJ DEFAULT) ERROR: Wave energy spectrum called '" << energySpectrumName 
	   << "' does not exist" << std::endl << write;
	 initialized = false;
      } else {
	 if (waveEnergySpectrum->initialize(sim,simClasses,cr,energySpectrumParameters) == false) {
	    simClasses.logger << "(SEP LAGR INJ DEFAULT) ERROR: Wave energy spectrum failed to initialize" 
	      << std::endl << write;
	    waveEnergySpectrum->finalize();
	    waveEnergySpectrum = NULL;
	    initialized = false;
	 }
	 
	 if (species.propagationDirection > 0) simControl.parWaveEnergy = waveEnergySpectrum;
	 else simControl.antiparWaveEnergy = waveEnergySpectrum;
      }
      
      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }

      resolution = 1;
      
      return initialized;
   }
   
} // namespace sep
   
#endif
