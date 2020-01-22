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

#ifndef SEP_LAGR_INJECTOR_NG_REAMES_1994_H
#define SEP_LAGR_INJECTOR_NG_REAMES_1994_H

#include <cstdlib>
#include <climits>
#include <stdint.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>

#include "sep_object_wrapper.h"
#include "sep_namespaces.h"
#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_lagr_definition.h"
#include "sep_fields_container.h"

namespace sep {

   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class LagrangianInjectorNgReames1994: public ParticleInjectorBase {
    public: 
      LagrangianInjectorNgReames1994();
   
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      //Real getMaximumSpeed();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);
   
    private:
      bool initialized;
      bool initialInjectionDone;
      bool initialInjection;
      bool inflowInjection;
      uint32_t N_particlesPerCell;        /**< Number of macroparticles injected to each cell.*/
      SPECIES species;

      lagrinjector::InflowBoundary inflowBoundary;
      pargrid::NeighbourID nbrMinusR;
      pargrid::NeighbourID nbrPlusR;

      Real maxInjectionTime;
      Real maxInjectionWavelength;
      Real minInjectionWavelength;
      
      uint32_t resolutionLambda;          /**< Number of injected particles per lambda cell.*/
      uint32_t resolutionRadial;          /**< Number of injected particles per radial cell.*/
      uint32_t superresolutionLambda;
      
      WaveEnergySpectrumBaseClass* waveEnergySpectrum;
      
      const Real DEF_VALUE;

      bool acceptCell(const pargrid::CellID& blockLID,const pargrid::CellID& blockGID,Real* normal);
      Real getWaveEnergy(Real l_min,Real l_max,Real lambda0,Real lambdaScaleFactor,Real B_mag2,Real U_relative,Real cellVolume);
      void injectParticles(pargrid::CellID block,Real* pos,Real B_mag,
			   unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* LINgReames1994Maker() {return new LagrangianInjectorNgReames1994<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::LagrangianInjectorNgReames1994(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity()) {
      initialInjection = true;
      waveEnergySpectrum = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::acceptCell(const pargrid::CellID& blockLID,const pargrid::CellID& blockGID,Real* normal) {
      const uint32_t nbrFlags = simClasses->pargrid.getNeighbourFlags(blockLID);

      Real B[3];
      Real injPos[3];
      Real V_wave[3];
      Real dV_wave;
      normal[0] = 0;
      normal[1] = 0;
      normal[2] = 0;
      uint32_t i_block,j_block,k_block;
      switch (inflowBoundary) {
       case lagrinjector::X_NEG:
	 // Inject to cell if it's missing -x neighbor:
	 if (((nbrFlags >> nbrMinusR) & 1) == 0) {normal[0] = 1; return true;}
	 break;
       case lagrinjector::X_POS:
	 // Inject to cell if it's missing +x neigbor:
	 if (((nbrFlags >> nbrPlusR) & 1) == 0) {normal[0] = -1; return true;}
	 break;
       case lagrinjector::VELOCITY_DIR:
	 // Inject to cell if wave velocity vector points to simulation box:
	 if (((nbrFlags >> nbrMinusR) & 1) == 0) normal[0] = +1.0;
	 if (((nbrFlags >> nbrPlusR) & 1) == 0) normal[0] = -1.0;
	 
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 injPos[0] = i_block*block::WIDTH_X + 0.5;
	 injPos[1] = j_block*block::WIDTH_Y + 0.5;
	 injPos[2] = k_block*block::WIDTH_Z + 0.5;

	 Real V_alfven;
	 (*simControl.fieldsGetState)(blockLID,sim->t,injPos,B,V_wave,dV_wave,V_alfven,species.propagationDirection);
	 if (dotProduct<3>(V_wave,normal) > 0) return true;
	 break;
       default:
	 break;
      }
      
      return false;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".macroparticles_per_cell","Number of macroparticles injected per cell (integer)",(uint32_t)0);
      cr.add(PREFIX+".initial_state_injection","If 'yes' initial state is created, defaults to 'yes' (string)",std::string("yes"));
      cr.add(PREFIX+".resolution_lambda","Number of injected particles per lambda cell, defaults to 1 (int)",(uint)1);
      cr.add(PREFIX+".superresolution_lambda","Number of injected wave packets near zero wavelength (int)",(uint)1);
      cr.add(PREFIX+".resolution_r","Number of injected particles per radial cell, defaults to 1 (int)",(uint)1);
      cr.add(PREFIX+".energy_spectrum","Name of wave energy spectrum function (string)",std::string(""));
      cr.add(PREFIX+".energy_spectrum_parameters","Name of config file region containing parameters for wave energy spectrum function (string)",std::string(""));
      cr.add(PREFIX+".min_injection_wavelength","Minimum wavelength where wave packets are injected at inflow boundary (float)",DEF_VALUE);
      cr.add(PREFIX+".max_injection_wavelength","Maximum wavelength where wave packets are injected at inflow boundary (float)",DEF_VALUE);
      cr.add(PREFIX+".max_injection_time","Maximum injection time (float)",std::numeric_limits<Real>::infinity());
      cr.add(PREFIX+".inflow_boundary_injection","If 'yes' wave packets are injected at inflow boundary (string)",std::string("no"));
      cr.add(PREFIX+".inflow_boundary","Either '-x', '+x', or 'velocity', determines which boundary is treated as inflow boundary (string)",std::string("velocity"));
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::finalize() {
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
   bool LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      bool success = true;
      if (sim->t >= maxInjectionTime) return initialized;

      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      recalculateCellVolumes(*sim,*simClasses);

      const std::vector<pargrid::CellID>& exteriorBlocks = simClasses->pargrid.getExteriorCells();

      Real injPos[3];
      Real B[3];
      Real V_wave[3];
      Real dV_wave;

      // Initial state injection -- inject wave packets to inflow boundary only:
      if (initialInjectionDone == false && inflowInjection == true) {
	 // Initial injection at time sim->t:
	 for (size_t b=0; b<exteriorBlocks.size(); ++b) {
	    // Get block global ID and calculate bounding box (i,j,k) indices:
	    const pargrid::CellID blockLID = exteriorBlocks[b];
	    const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	    Real normal[3];
	    if (acceptCell(blockLID,blockGID,normal) == false) continue;

	    uint32_t i_block,j_block,k_block;
	    block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	    
	    // Get plasma state at cell centroid:
	    injPos[0] = i_block*block::WIDTH_X + 0.5;
	    injPos[1] = j_block*block::WIDTH_Y + 0.5;
	    injPos[2] = k_block*block::WIDTH_Z + 0.5;
	    
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,injPos,B,V_wave,dV_wave,V_alfven,species.propagationDirection);
	    const Real B_mag = vectorMagnitude<3>(B);

	    if (normal[0] > 0) {
	       injPos[0] = i_block*block::WIDTH_X     + 1e-6;
	    } else {
	       injPos[0] = (i_block+1)*block::WIDTH_X - 1e-6;
	    }

	    injectParticles(blockLID,injPos,B_mag,N_particles,wrapper);
	 }

	 initialInjectionDone = true;
      }

      // Regular injection -- inject new wave packets to inflow boundary cells 
      // if the cells have too few packets:
      if (inflowInjection == false) {
	 return initialized;
      }
      Real t_propag = 0.0;
      for (size_t b=0; b<exteriorBlocks.size(); ++b) {
	 // Measure block injection time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();

	 const pargrid::CellID blockLID = exteriorBlocks[b];

	 // Get block global ID and calculate bounding box (i,j,k) indices:
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 Real normal[3];
	 if (acceptCell(blockLID,blockGID,normal) == false) continue;
	 
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

	 // HACK: Calculate the minimum x-coordinate of wave packet 
	 // that already exists in block:
	 Real x_min;
	 if (normal[0] > 0) {
	    x_min = i_block+2;
	    for (uint64_t i=0; i<wrapper.size(blockLID); ++i) {
	       if (wrapper.data()[blockLID][i].state[0] < x_min) x_min = wrapper.data()[blockLID][i].state[0];
	    }
	 } else {
	    x_min = i_block*block::WIDTH_X-1;
	    for (uint64_t i=0; i<wrapper.size(blockLID); ++i) {
	       if (wrapper.data()[blockLID][i].state[0] > x_min) x_min = wrapper.data()[blockLID][i].state[0];
	    }
	 }
	 // END HACK

	 // Get plasma state at cell centroid:
	 injPos[0] = i_block*block::WIDTH_X + 0.5;
	 injPos[1] = j_block*block::WIDTH_Y + 0.5;
	 injPos[2] = k_block*block::WIDTH_Z + 0.5;
	 Real V_alfven;
	 (*simControl.fieldsGetState)(blockLID,sim->t+sim->dt,injPos,B,V_wave,dV_wave,V_alfven,species.propagationDirection);
	 const Real B_mag = vectorMagnitude<3>(B);

	 // Calculate wavelength scale factor:
	 //const Real lambdaScaleFactor = wlengthScaler->getScaleFactor(sim->t,injPos);

	 // Calculate injection position:
	 Real dx = sim->dx_cell[i_block*block::WIDTH_X] / resolutionRadial;
	 Real injectionInterval = dx / V_wave[0];
	 int32_t intervalNumber = static_cast<int32_t>((sim->t+sim->dt)/injectionInterval);
	 Real t_offset = sim->t+sim->dt - intervalNumber*injectionInterval;
	 Real x_offset = V_wave[0]*t_offset / sim->dx_cell[i_block*block::WIDTH_X];

	 if (normal[0] > 0) {
	    injPos[0] = i_block*block::WIDTH_X + x_offset;
	    if (injPos[0] < i_block*block::WIDTH_X) continue;
	 } else {
	    injPos[0] = (i_block+1)*block::WIDTH_X + x_offset;
	    if (injPos[0] > (i_block+1)*block::WIDTH_X) continue;
	 }

	 // HACK: Skip injection if the injection position already contains 
	 // a wave packet. This is not pretty but works beautifully:
	 if (fabs(x_min-injPos[0]) > 1.0e-2/resolutionRadial) {
	    injectParticles(blockLID,injPos,B_mag,N_particles,wrapper);
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
   void LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,Real* pos,Real B_mag,
									  unsigned int* N_particles,
									  pargrid::DataWrapper<PARTICLE>& wrapper) {
      Real l_min = simControl.wavelengthMeshNodeCoordinates[0];
      Real l_max = simControl.wavelengthMeshNodeCoordinates[simControl.N_wavelengthMeshCells-1];
      if (minInjectionWavelength != DEF_VALUE) l_min = minInjectionWavelength;
      if (maxInjectionWavelength != DEF_VALUE) l_max = maxInjectionWavelength;

      Real l_min_logical = (*simControl.getLogicalWavelength)(l_min);
      Real l_max_logical = (*simControl.getLogicalWavelength)(l_max)+1;

      const int32_t l_min_index = std::max(0,static_cast<int32_t>(floor(l_min_logical)));
      const int32_t l_max_index = std::min((int32_t)simControl.N_wavelengthMeshCells,static_cast<int32_t>(floor(l_max_logical)));

      int count = 0;
      PARTICLE particle;
      const Real volume = simControl.cellVolumes[blockLID];
      particle.state[lagr::XCRD]   = pos[0];
      particle.state[lagr::YCRD]   = pos[1];
      particle.state[lagr::ZCRD]   = pos[2];
      
      const int32_t N = simControl.N_wavelengthMeshCells;
      for (int32_t l=l_min_index; l<N/2-1; ++l) {
	 for (uint32_t p=0; p<resolutionLambda; ++p) {
	    const Real l_logical  = l + (1+2*p)*0.5/resolutionLambda;
	    const Real l_min = (*simControl.getPhysicalWavelength)(l_logical - 0.5/resolutionLambda);
	    const Real l_max = (*simControl.getPhysicalWavelength)(l_logical + 0.5/resolutionLambda);
	    Real l_logical_avg = (*simControl.getLogicalWavelength)(0.5*(l_min+l_max));
	    if (fabs(l_logical_avg - simControl.N_wavelengthMeshCells/2) < 1e-3) {
	       l_logical_avg = l_logical;
	    }
	    particle.state[lagr::LAMBDA] = l_logical_avg;

	    const Real minLambda = l_min;
	    const Real maxLambda = l_max;
	    const Real injWaveEnergy = waveEnergySpectrum->getEnergyDensity(B_mag,minLambda,maxLambda)*volume;
	    particle.state[lagr::ENERGY] = injWaveEnergy/resolutionRadial;
	    
	    ++N_particles[blockLID];
	    wrapper.push_back(blockLID,particle);
	    ++count;
	 }
      }
      
      for (int32_t l=N/2-1; l<N/2+1; ++l) {
	 const uint32_t superres = superresolutionLambda;
	 for (uint32_t p=0; p<superres; ++p) {
	    
	    const Real dl = 1.0 / superres;
	    Real l_logical = l + 0.5/superres + p*dl;
	    Real l_phys_min = (*simControl.getPhysicalWavelength)(1e-6 + l_logical - 0.5*dl);
	    Real l_phys_max = (*simControl.getPhysicalWavelength)(1e-6 + l_logical + 0.5*dl);

	    if (l_logical < N/2) {
	       l_phys_max = std::min(l_phys_max,0.0);
	    } else {
	       l_phys_min = std::max(l_phys_min,0.0);
	    }
	    
	    particle.state[lagr::LAMBDA] = l_logical;
	    const Real injWaveEnergy = waveEnergySpectrum->getEnergyDensity(B_mag,l_phys_min,l_phys_max)*volume;
	    particle.state[lagr::ENERGY] = injWaveEnergy/resolutionRadial;
	    
	    ++N_particles[blockLID];
	    wrapper.push_back(blockLID,particle);
	    ++count;
	 }
      }

      for (int32_t l=N/2+1; l<l_max_index; ++l) {
	 for (uint32_t p=0; p<resolutionLambda; ++p) {
	    const Real l_logical  = l + (1+2*p)*0.5/resolutionLambda;
	    const Real l_min = (*simControl.getPhysicalWavelength)(l_logical - 0.5/resolutionLambda);
	    const Real l_max = (*simControl.getPhysicalWavelength)(l_logical + 0.5/resolutionLambda);
	    Real l_logical_avg = (*simControl.getLogicalWavelength)(0.5*(l_min+l_max));
	    if (fabs(l_logical_avg - simControl.N_wavelengthMeshCells/2) < 1e-3) {
	       l_logical_avg = l_logical;
	    }
	    particle.state[lagr::LAMBDA] = l_logical_avg;
	    
	    const Real minLambda = l_min;
	    const Real maxLambda = l_max;
	    const Real injWaveEnergy = waveEnergySpectrum->getEnergyDensity(B_mag,minLambda,maxLambda)*volume;
	    particle.state[lagr::ENERGY] = injWaveEnergy/resolutionRadial;
	    
	    ++N_particles[blockLID];
	    wrapper.push_back(blockLID,particle);
	    ++count;
	 }
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool LagrangianInjectorNgReames1994<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								     const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP LAGR INJ NG REAMES) Starting to init." << std::endl << write;
      
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      // Parse config file options:
      std::string injectInitialString,inflowInjectionString,inflowBoundaryString;
      std::string energySpectrumName,energySpectrumParameters;
      cr.parse();
      cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell);
      cr.get(regionName+".initial_state_injection",injectInitialString);
      cr.get(regionName+".resolution_lambda",resolutionLambda);
      cr.get(regionName+".resolution_r",resolutionRadial);
      if (cr.get(regionName+".energy_spectrum",energySpectrumName) == false) initialized = false;
      if (cr.get(regionName+".energy_spectrum_parameters",energySpectrumParameters) == false) initialized = false;
      if (cr.get(regionName+".min_injection_wavelength",minInjectionWavelength) == false) initialized = false;
      if (cr.get(regionName+".max_injection_wavelength",maxInjectionWavelength) == false) initialized = false;
      if (cr.get(regionName+".max_injection_time",maxInjectionTime) == false) initialized = false;
      if (cr.get(regionName+".inflow_boundary_injection",inflowInjectionString) == false) initialized = false;
      if (cr.get(regionName+".inflow_boundary",inflowBoundaryString) == false) initialized = false;
      if (cr.get(regionName+".superresolution_lambda",superresolutionLambda) == false) initialized = false;

      inflowInjection = false;
      if (inflowInjectionString == "yes") inflowInjection = true;

      inflowBoundary = lagrinjector::InflowBoundary::VELOCITY_DIR;
      if (inflowBoundaryString == "-x") inflowBoundary = lagrinjector::InflowBoundary::X_NEG;
      else if (inflowBoundaryString == "+x") inflowBoundary = lagrinjector::InflowBoundary::X_POS;
      else if (inflowBoundaryString == "velocity") inflowBoundary = lagrinjector::InflowBoundary::VELOCITY_DIR;
      
      // Get energy spectrum function and initialize it:
      waveEnergySpectrum = sep::getObjectWrapper().waveEnergySpectrumFactory.create(energySpectrumName);
      if (waveEnergySpectrum == NULL) {
	 simClasses.logger << "(SEP LAGR INJ NG REAMES) ERROR: Wave energy spectrum called '" << energySpectrumName
	   << "' does not exist" << std::endl << write;
	 initialized = false;
      } else {
	 if (waveEnergySpectrum->initialize(sim,simClasses,cr,energySpectrumParameters) == false) {
	    simClasses.logger << "(SEP LAGR INJ NG REAMES) ERROR: Wave energy spectrum failed to initialize"
	      << std::endl << write;
	    waveEnergySpectrum->finalize();
	    delete waveEnergySpectrum; waveEnergySpectrum = NULL;
	    initialized = false;
	 }
	 
	 if (species.propagationDirection > 0) simControl.parWaveEnergy = waveEnergySpectrum;
	 else simControl.antiparWaveEnergy = waveEnergySpectrum;
      }
      
      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP LAGR INJ NG REAMES) ERROR: Parameter '" << regionName+".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      
      if (injectInitialString == "yes") initialInjection = true;
      else initialInjection = false;

      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }
      
      nbrMinusR = simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0);
      nbrPlusR  = simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0);
      
      simClasses.logger << "(SEP LAGR INJ NG REAMES) Initialization done, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << std::endl << write;
      else simClasses.logger << "FAILURE" << std::endl << write;
      
      return initialized;
   }

} // namespace sep
   
#endif
