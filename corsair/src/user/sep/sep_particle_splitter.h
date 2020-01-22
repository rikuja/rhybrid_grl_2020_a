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

#ifndef SEP_PARTICLE_SPLITTER_H
#define SEP_PARTICLE_SPLITTER_H

#include <cstdlib>
#include <vector>
#include <limits>

#include <main.h>
#include <linear_algebra.h>

#include "sep_particle_definition.h"
#include "sep_simcontrol.h"
#include "sep_base_class_particle_splitter.h"

namespace sep {
   
   extern sep::SimControl simControl;
   static std::string prefix = "ParticleSplitter";

   template<class SPECIES,class PARTICLE>
   class ParticleSplitter: public ParticleSplitterBase {
    public:
      ParticleSplitter();
      virtual ~ParticleSplitter();

      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      bool split();
      
    private:
      bool checkWaveEnergies();
      bool highEnergySplitter();

      int32_t highEnergySplitterBinsPerInterval;
      Real highEnergySplitterMinEnergy;
      int32_t highEnergySplitterMinMacroparticlesPerCell;
      int32_t highEnergySplits;                           /**< If an energy bin has less than minMacroparticlesInBin 
							   * macroparticles in it, each macroparticle is split
							   * highEnergySplits times.*/
      uint32_t highEnergySplitterMinMacroparticlesPerBin; /**< Minimum number of macroparticles in each energy bin, 
							   * parameter for high energy splitter.*/
      
      int32_t waveEnergySplits;
      Real waveEnergySplitterMaxEnergyRatio;
      int32_t waveEnergySplitterMaxMacroparticlesPerCell; /**< If number of macroparticles in spatial cell is larger
							   * than this value, then wave energy splitter switches off.*/
      
      bool addConfigFileItems();
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleSplitterBase* ParticleSplitterMaker() {return new ParticleSplitter<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   ParticleSplitter<SPECIES,PARTICLE>::ParticleSplitter(): ParticleSplitterBase() {
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   ParticleSplitter<SPECIES,PARTICLE>::~ParticleSplitter() {
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::addConfigFileItems() {
      corsair::ObjectWrapper& wrapper = corsair::getObjectWrapper();

      wrapper.configReader.add(prefix+".high_energy_splitter.bins_per_interval","Number of bins per unit logarithmic interval (int)",(int32_t)3);
      wrapper.configReader.add(prefix+".high_energy_splitter.splits","How many particles are produced at each high energy split (int)?",(int32_t)5);
      wrapper.configReader.add(prefix+".high_energy_splitter.min_energy","Minimum energy for high energy splitter in MeV (float)",(Real)1.0);
      wrapper.configReader.add(prefix+".high_energy_splitter.min_macroparticles","Minimum number of macroparticles in cell until high energy splitter is turned on (int)",(int32_t)10);
      wrapper.configReader.add(prefix+".high_energy_splitter.min_macroparticles_per_bin","Macroparticles are not split if bin has at least this many particles",(int32_t)30);
      wrapper.configReader.add(prefix+".wave_energy_splitter.max_energy_ratio","Macroparticle is split if it's energy exceeds spectral wave energy times this value (float)",(Real)0.25);
      wrapper.configReader.add(prefix+".wave_energy_splitter.max_macroparticles_per_cell","Wave energy splitter is turned off is spatial cell contains at least this many macroparticles (int)",(int32_t)(1500));
      wrapper.configReader.add(prefix+".wave_energy_splitter.splits","How many particles are produces at each wave energy split (int)?",(int32_t)4);
      
      return true;
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::checkWaveEnergies() {
      // There are no particles in simulation until setup time has passed:
      corsair::ObjectWrapper& owrapper = corsair::getObjectWrapper();
      if (owrapper.sim.t < simControl.t_setup) return true;

      // Get wave energy arrays (may be NULL):
      const Real* RESTRICT waveEnergyAntipar = owrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      const Real* RESTRICT waveEnergyPar     = owrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      if (waveEnergyAntipar == NULL && waveEnergyPar == NULL) return true;
      
      // TEMP
      if (waveEnergyPar == NULL) return true;
      // END TEMP
      
      const int32_t NWL = simControl.N_wavelengthMeshCells;
      std::vector<Real> minWaveEnergyDens(NWL/2);
      
      Real t_propag = 0.0;
      for (pargrid::CellID blockLID=0; blockLID<owrapper.simClasses.pargrid.getNumberOfLocalCells(); ++blockLID) {
	 // Measure computation time if we are testing for repartitioning:
	 if (owrapper.sim.countPropagTime == true) t_propag = MPI_Wtime();

	 // Get minimum value of wave energy / wavelength in each wavelength bin (par waves only):
	 if (getMinimumWaveEnergies(blockLID,+1,minWaveEnergyDens) == false) continue;
	 
	 // Calculate wave speed and B in block centroid.
	 // We need these to approximate particle resonant wavelengths:
	 pargrid::CellID blockGID = owrapper.simClasses.pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(owrapper.sim,blockGID,i_block,j_block,k_block);
	 
	 Real pos[4];
	 pos[0] = i_block*block::WIDTH_X + 0.5*block::WIDTH_X;
	 pos[1] = j_block*block::WIDTH_Y + 0.5*block::WIDTH_Y;
	 pos[2] = k_block*block::WIDTH_Z + 0.5*block::WIDTH_Z;
	 pos[3] = 0.0;
	 
	 Real B[3];
	 Real V_wave[3];
	 Real dV_wave,V_alfven;
	 (*simControl.fieldsGetState)(blockLID,owrapper.sim.t,pos,B,V_wave,dV_wave,V_alfven,+1);
	 const Real waveSpeed = vectorMagnitude<3>(V_wave);
	 const Real parWaveSpeed = dotProduct<3>(V_wave,B) / vectorMagnitude<3>(B);
	 const Real B_mag = vectorMagnitude<3>(B);

	 for (size_t plist=0; plist<owrapper.particleLists.size(); ++plist) {
	    // Skip wave packets:
	    if (owrapper.particleLists[plist]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 
	    // Get particle data:
	    pargrid::DataID dataID;
	    if (owrapper.particleLists[plist]->getParticles(dataID) == false) continue;
	    pargrid::DataWrapper<Particle<Real> > wrapper = owrapper.simClasses.pargrid.getUserDataDynamic<Particle<Real> >(dataID);
	    if (wrapper.valid() == false) continue;

	    // Get species struct:
	    const SPECIES* species = reinterpret_cast<const SPECIES*>(owrapper.particleLists[plist]->getSpecies());

	    // Calculate gyro frequency:
	    const Real omega = species->q_per_m * vectorMagnitude<3>(B);
	    
	    // Buffer where splitted particles are temporarily inserted:
	    std::vector<Particle<Real> > splitBuffer;
	    std::vector<pargrid::ArraySizetype> indices;
	    
	    // Test all particles for splitting:
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       // Speed in wave rest frame:
	       Real parSpeed2 = wrapper.data()[blockLID][p].state[particle::V_PAR] - parWaveSpeed;
	       parSpeed2 *= parSpeed2;
	       const Real gyroSpeed2 = 2*wrapper.data()[blockLID][p].state[particle::MU]*B_mag/species->mass;
	       const Real speed = sqrt(parSpeed2 + gyroSpeed2);
	       
	       // Max. resonant wavelenght (physical and logical):
	       Real maxResonantLambdaPhys = -2*M_PI/omega*speed;
	       Real maxResonantLambdaLogi = (*simControl.getLogicalWavelength)(maxResonantLambdaPhys);
	       const int32_t l_index = static_cast<int32_t>(maxResonantLambdaLogi);
	       
	       // Change in wave energy / wavelength this macroparticle may cause:
	       Real energyChange = 0.5 * species->mass / M_PI * omega * waveSpeed;
	       energyChange *= wrapper.data()[blockLID][p].state[particle::WEIGHT];
	       
	       // If energy change exceeds threshold, insert a copy of particle to buffer:
	       if (energyChange/minWaveEnergyDens[l_index] > waveEnergySplitterMaxEnergyRatio) {
		  indices.push_back(p);
		  splitBuffer.push_back(wrapper.data()[blockLID][p]);
	       }
	    }

	    // Do not split particles if cell has too many macroparticles in it:
	    if (wrapper.size(blockLID) + splitBuffer.size() > (uint32_t)waveEnergySplitterMaxMacroparticlesPerCell) continue;

	    // Divide original macroparticle weights by waveEnergySplits:
	    for (size_t i=0; i<indices.size(); ++i) wrapper.data()[blockLID][indices[i]].state[particle::WEIGHT] /= waveEnergySplits;

	    // Inject waveEnergySplits-1 new particles to simulation:
	    for (size_t i=0; i<splitBuffer.size(); ++i) {
	       splitBuffer[i].state[particle::WEIGHT] /= waveEnergySplits;
	       for (int32_t s=1; s<waveEnergySplits; ++s) {
		  wrapper.push_back(blockLID,splitBuffer[i]);
	       }
	    }
	    owrapper.particleLists[plist]->getParticleNumberArray()[blockLID] += splitBuffer.size()*(waveEnergySplits-1);
	 }
	    
	 // Measure computation time if we are testing for repartitioning:
	 if (owrapper.sim.countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    owrapper.simClasses.pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
	 
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::finalize() {
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::highEnergySplitter() {
      // There are no particles in sim until setup time has passed:
      corsair::ObjectWrapper& owrapper = corsair::getObjectWrapper();
      if (owrapper.sim.t < simControl.t_setup) return true;
      
      int minLog10Energy = -2;
      int maxLog10Energy = +4;
      const Real d_log10_energy = 1.0 / highEnergySplitterBinsPerInterval;
      
      int N_bins = highEnergySplitterBinsPerInterval * (maxLog10Energy - minLog10Energy);
      std::vector<uint32_t> N_macroparticles(N_bins);
      std::vector<bool> splitParticlesInBin(N_bins);

      for (size_t plist=0; plist<owrapper.particleLists.size(); ++plist) {
	 // Skip wave packets:
	 if (owrapper.particleLists[plist]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
      
	 // Get particle data:
	 pargrid::DataID dataID;
	 if (owrapper.particleLists[plist]->getParticles(dataID) == false) continue;
	 pargrid::DataWrapper<Particle<Real> > wrapper = owrapper.simClasses.pargrid.getUserDataDynamic<Particle<Real> >(dataID);
	 if (wrapper.valid() == false) continue;
	 
	 // Get species struct:
	 const SPECIES* species = reinterpret_cast<const SPECIES*>(owrapper.particleLists[plist]->getSpecies());

	 std::vector<Particle<Real> > splitBuffer;
	 Real t_propag = 0.0;
	 for (pargrid::CellID blockLID=0; blockLID<owrapper.simClasses.pargrid.getNumberOfLocalCells(); ++blockLID) {
	    // Do not split particles in cells where there few macroparticles:
	    if (wrapper.size(blockLID) < (uint32_t)highEnergySplitterMinMacroparticlesPerCell) continue;

	    // Measure computation time if we are testing for repartitioning:
	    if (owrapper.sim.countPropagTime == true) t_propag = MPI_Wtime();
	    splitBuffer.clear();

	    // Clear N_macroparticles:
	    for (size_t i=0; i<N_macroparticles.size(); ++i) N_macroparticles[i] = 0;
	    for (size_t i=0; i<splitParticlesInBin.size(); ++i) splitParticlesInBin[i] = false;

	    // Calculate block centroid coordinates:
	    pargrid::CellID blockGID = owrapper.simClasses.pargrid.getGlobalIDs()[blockLID];
	    uint32_t i_block,j_block,k_block;
	    block::calculateBlockIndices(owrapper.sim,blockGID,i_block,j_block,k_block);
	                
	    Real pos[4];
	    pos[0] = i_block*block::WIDTH_X + 0.5*block::WIDTH_X;
	    pos[1] = j_block*block::WIDTH_Y + 0.5*block::WIDTH_Y;
	    pos[2] = k_block*block::WIDTH_Z + 0.5*block::WIDTH_Z;
	    pos[3] = 0.0;
	    
	    // Get fields in block centroid:
	    Real E[3];
	    Real B[3];
	    Real gradB[9];
	    (*simControl.fieldsGetFields)(blockLID,owrapper.sim.t,pos,E,B,gradB);
	    const Real B_mag = vectorMagnitude<3>(B);
	    
	    // Calculate number of macroparticles in each energy bin:
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       const Real parSpeed2 = wrapper.data()[blockLID][p].state[particle::V_PAR]
		                    * wrapper.data()[blockLID][p].state[particle::V_PAR];
	       const Real gyroSpeed2 = 2*wrapper.data()[blockLID][p].state[particle::MU]*B_mag/species->mass;
	       const Real energyPerAmu = 0.5*constants::MASS_PROTON*(parSpeed2 + gyroSpeed2)
		                       / constants::CHARGE_ELEMENTARY / 1e6;
	       
	       int index = static_cast<int>((log10(energyPerAmu)-minLog10Energy) / d_log10_energy);
	       if (index < 0) index = 0;
	       if (index > (int32_t)N_macroparticles.size()-1) index = N_macroparticles.size()-1;
	       ++N_macroparticles[index];
	    }

	    // Determine if particles in different energy intervals should be split:
	    for (size_t i=0; i<N_macroparticles.size(); ++i) {
	       if (N_macroparticles[i] == 0) continue;
	       if (minLog10Energy + i*d_log10_energy < highEnergySplitterMinEnergy) continue;
	       if (N_macroparticles[i] < highEnergySplitterMinMacroparticlesPerBin) splitParticlesInBin[i] = true;
	    }

	    // Split particles:
	    for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	       const Real parSpeed2 = wrapper.data()[blockLID][p].state[particle::V_PAR]
		 * wrapper.data()[blockLID][p].state[particle::V_PAR];
	       const Real gyroSpeed2 = 2*wrapper.data()[blockLID][p].state[particle::MU]*B_mag/species->mass;
	       const Real energyPerAmu = 0.5*constants::MASS_PROTON*(parSpeed2 + gyroSpeed2)
		 / constants::CHARGE_ELEMENTARY / 1e6;
	       
	       int index = static_cast<int>((log10(energyPerAmu)-minLog10Energy) / d_log10_energy);
	       if (index < 0) index = 0;
	       if (index > (int32_t)N_macroparticles.size()-1) index = N_macroparticles.size()-1;
	       if (splitParticlesInBin[index] == false) continue;

	       wrapper.data()[blockLID][p].state[particle::WEIGHT] /= highEnergySplits;
	       splitBuffer.push_back(wrapper.data()[blockLID][p]);
	    }
	    
	    for (size_t i=0; i<splitBuffer.size(); ++i) {
	       for (int32_t j=1; j<highEnergySplits; ++j) {
		  wrapper.push_back(blockLID,splitBuffer[i]);
	       }
	    }
	    owrapper.particleLists[plist]->getParticleNumberArray()[blockLID] += splitBuffer.size()*(highEnergySplits-1);

	    // Measure computation time if we are testing for repartitioning:
	    if (owrapper.sim.countPropagTime == true) {
	       t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	       owrapper.simClasses.pargrid.getCellWeights()[blockLID] += t_propag;
	    }
	 }
	 
      }
      
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      initialized = ParticleSplitterBase::initialize(sim,simClasses,cr);
      
      addConfigFileItems();
      cr.parse();
      cr.get(prefix+".high_energy_splitter.splits",highEnergySplits);
      cr.get(prefix+".high_energy_splitter.min_energy",highEnergySplitterMinEnergy);
      cr.get(prefix+".high_energy_splitter.min_macroparticles",highEnergySplitterMinMacroparticlesPerCell);
      cr.get(prefix+".high_energy_splitter.bins_per_interval",highEnergySplitterBinsPerInterval);
      cr.get(prefix+".high_energy_splitter.min_macroparticles_per_bin",highEnergySplitterMinMacroparticlesPerBin);
      cr.get(prefix+".wave_energy_splitter.splits",waveEnergySplits);
      cr.get(prefix+".wave_energy_splitter.max_macroparticles_per_cell",waveEnergySplitterMaxMacroparticlesPerCell);
      cr.get(prefix+".wave_energy_splitter.max_energy_ratio",waveEnergySplitterMaxEnergyRatio);

      if (highEnergySplitterMinMacroparticlesPerCell < 0) initialized = false;
      highEnergySplitterMinEnergy = log10(highEnergySplitterMinEnergy / 1e6);

      if (highEnergySplits < 1) {
	 simClasses.logger << "(SEP PARTICLE SPLITTER) ERROR: Parameter '" << prefix+".high_energy_splitter.splits' must have value >= 1" << std::endl << write;
	 initialized = false;
      }
      if (waveEnergySplitterMaxEnergyRatio <= 0.0 || waveEnergySplitterMaxEnergyRatio >= 1.0) {
	 simClasses.logger << "(SEP PARTICLE SPLITTER) ERROR: Parameter '" << prefix+".max_energy_ratio' must have a positive value below unity" << std::endl << write;
	 initialized = false;
      }
      if (waveEnergySplitterMaxMacroparticlesPerCell < 1) {
	 simClasses.logger << "(SEP PARTICLE SPLITTER) ERROR: Parameter '" << prefix+".wave_energy_splitter.max_macroparticles_per_cell' must have value >= 1" << std::endl << write;
	 initialized = false;
      }
      if (waveEnergySplits < 1) {
	 simClasses.logger << "(SEP PARTICLE SPLITTER) ERROR: Parameter '" << prefix+".wave_energy_splitter.splits' must have value >= 1" << std::endl << write;
	 initialized = false;
      }

      // Write log file message:
      simClasses.logger << "(SEP PARTICLE SPLITTER) Initialization complete, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << std::endl;
      else simClasses.logger << "FAILURE" << std::endl;
      
      simClasses.logger << "\t High energy splitter min energy log10(MeV): " << highEnergySplitterMinEnergy << std::endl;
      simClasses.logger << "\t                                     splits: " << highEnergySplits << std::endl;
      simClasses.logger << "\t                         min macroparticles: " << highEnergySplitterMinMacroparticlesPerCell << std::endl;
      simClasses.logger << "\t                      bins per log interval: " << highEnergySplitterBinsPerInterval << std::endl;
      simClasses.logger << "\t                 min macroparticles per bin: " << highEnergySplitterMinMacroparticlesPerBin << std::endl;
      simClasses.logger << "\t Wave energy splitter      max energy ratio: " << waveEnergySplitterMaxEnergyRatio << std::endl;
      simClasses.logger << "\t                max macroparticles per cell: " << waveEnergySplitterMaxMacroparticlesPerCell << std::endl;
      simClasses.logger << "\t                                     splits: " << waveEnergySplits << std::endl;
      simClasses.logger << write;

      return initialized;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleSplitter<SPECIES,PARTICLE>::split() {      
      bool success = true;
      if (corsair::getObjectWrapper().sim.timestep % 5 != 0) return success;
      
      #if PROFILE_LEVEL > 0
         static int profCounter = -1;
         profile::start("(SEP) Particle Splitter",profCounter);
      #endif
      
      if (checkWaveEnergies() == false) success = false;
      if (highEnergySplitter() == false) success = false;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

} // namespace sep
   
#endif
