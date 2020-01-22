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

#ifndef DEFAULT_INJECTOR_H
#define DEFAULT_INJECTOR_H

#include <cstdlib>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>

template<class FIELD,class SPECIES,class PARTICLE>
class DefaultInjector: public ParticleInjectorBase {
 public: 
   DefaultInjector();

   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);
   
 private:
   enum EnergyDistrib {Gaussian,PowerLaw};
   enum PitchDistrib {Isotropic};
   
   EnergyDistrib energyInjector;     /**< Injection energy distribution.*/
   PitchDistrib pitchInjector;       /**< Injection pitch distribution.*/
   int N_particlesPerCell;           /**< Number of macroparticles per cell.*/
   Real numberDensity;               /**< Number density of particles in each cell.*/
   const SPECIES* species;
   
   // Variables needed for power law injection:
   Real deltaEnergy;                 /**< Max. energy - Min. energy.*/
   Real maxEnergy;                   /**< Minimum injection energy.*/
   Real minEnergy;                   /**< Maximum injection energy.*/
   Real normalization;               /**< Normalization constant of the power law.*/
   Real sigma;                       /**< Index of the power law distribution.*/
   
   // Miscellaneous class variables:
   bool initialInjectionDone;
   bool initialized;

   Real distribPitchIsotropic(SimulationClasses& simClasses);
   bool injectParticles(pargrid::CellID block,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper);
};

// Maker function that returns a new DefaultInjector:
template<class FIELD,class SPECIES,class PARTICLE> inline
ParticleInjectorBase* DefInjectorMaker() {return new DefaultInjector<FIELD,SPECIES,PARTICLE>();}

template<class FIELD,class SPECIES,class PARTICLE> inline
DefaultInjector<FIELD,SPECIES,PARTICLE>::DefaultInjector(): ParticleInjectorBase() {
   initialInjectionDone = false;
   initialized = false;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool DefaultInjector<FIELD,SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
Real DefaultInjector<FIELD,SPECIES,PARTICLE>::distribPitchIsotropic(SimulationClasses& simClasses) {
   const Real pitch = 1.0 - 2.0*simClasses.random.uniform();
   return pitch;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool DefaultInjector<FIELD,SPECIES,PARTICLE>::finalize() {
   initialized = false;
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool DefaultInjector<FIELD,SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
   if (initialized == false) return initialized;
   bool success = true;
   
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   
   if (sim->timestep == 0 && initialInjectionDone == false) {
      for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
	 if (injectParticles(block,N_particles,wrapper) == false) success = false;
      }
      initialInjectionDone = true;
   } else {
      const std::vector<pargrid::CellID>& exteriorBlocks = simClasses->pargrid.getExteriorCells();
      for (size_t b=0; b<exteriorBlocks.size(); ++b) {
	 const pargrid::CellID block = exteriorBlocks[b];
	 if (injectParticles(block,N_particles,wrapper) == false) success = false;
      }
   }
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool DefaultInjector<FIELD,SPECIES,PARTICLE>::injectParticles(pargrid::CellID block,unsigned int* N_particles,
							      pargrid::DataWrapper<PARTICLE>& wrapper) {
   bool success = true;
   
   // Resize block's particle list to have space for injected particles:
   const pargrid::ArraySizetype oldSize = wrapper.size()[block];
   N_particles[block] += N_particlesPerCell;
   wrapper.resize(block,oldSize+N_particlesPerCell);
   
   Real energy = NAN;
   Real pitch = NAN;
   Real E[3];
   Real B[3];
   Real coords[3];
   Real V_drift[3];
   Real B_mag_inverse = 0.0;
   
   Real blockSize[3];
   getBlockSize(*simClasses,*sim,block,blockSize);
   Real cellSize[3];
   getBlockCellSize(*simClasses,*sim,block,cellSize);
   
   const double* coordinates = getBlockCoordinateArray(*sim,*simClasses);
   
   Real weightSum = 0.0;
   Real dE = (maxEnergy-minEnergy)/(N_particlesPerCell-1);
   
   unsigned int counter = oldSize;
   PARTICLE* particles = wrapper.data()[block];
   for (int p=0; p<N_particlesPerCell; ++p) {
      // First three elements returned by operator[] are reserved for particle coordinates:
      particles[counter][particle::XPOS] = simClasses->random.uniform()*blockSize[0];
      particles[counter][particle::YPOS] = simClasses->random.uniform()*blockSize[1];
      particles[counter][particle::ZPOS] = simClasses->random.uniform()*blockSize[2];
      
      // Get fields at particle position:
      coords[particle::XPOS] = particles[counter][particle::XPOS] + coordinates[3*block+particle::XPOS];
      coords[particle::YPOS] = particles[counter][particle::YPOS] + coordinates[3*block+particle::YPOS];
      coords[particle::ZPOS] = particles[counter][particle::ZPOS] + coordinates[3*block+particle::ZPOS];
      FIELD::getFields(block,sim->t,coords,E,B);
      B_mag_inverse = 1.0 / vectorMagnitude<3>(B);
      
      // Calculate injection pitch:
      switch (pitchInjector) {
       case (Isotropic):
	 pitch = distribPitchIsotropic(*simClasses);
	 break;
       default:
	 pitch = NAN;
	 break;
      }
      
      // Calculate injection energy:
      energy = minEnergy + p*dE;
      particles[counter][particle::WEIGHT] = pow(energy,sigma);
      weightSum += particles[counter][particle::WEIGHT];
      
      // Calculate energy in drift-free frame:
      #warning default_injector.h energy calculation only includes ExB drift!
      crossProduct(E,B,V_drift);
      const Real B_mag2 = vectorMagnitude2<3>(B);
      for (int i=0; i<3; ++i) V_drift[i] /= B_mag2;
      const Real energy_dot = energy - 0.5*species->m*vectorMagnitude2<3>(V_drift);
      
      // Calculate magnetic moment and parallel speed:
      particles[counter][particle::MU]   = energy_dot * B_mag_inverse * (1.0 - pitch*pitch);
      particles[counter][particle::VPAR] = pitch * sqrt(2.0*energy_dot/species->m);
      ++counter;
   }

   // Normalize particle weights so that the total number 
   // density equals the injection number density:
   counter = oldSize;
   const Real normalization = numberDensity*blockSize[0]*blockSize[1]*blockSize[2]/weightSum;
   for (int p=0; p<N_particlesPerCell; ++p) {
      particles[counter][particle::WEIGHT] *= normalization;
      ++counter;
   }   
   return success;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool DefaultInjector<FIELD,SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
							 const std::string& regionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
   species = reinterpret_cast<const SPECIES*>(plist->getSpecies());
   
   // Parse configuration file:
   std::string energyDistrib,pitchDistrib;
   const Real defValue = std::numeric_limits<Real>::infinity();
   cr.add(regionName+".particles_per_cell","Number of macroparticles to create per cell (int).",(int)-1);
   cr.add(regionName+".number_density","Number density of particles (float).",defValue);
   cr.add(regionName+".distribution_energy","Name of energy distribution function (gaussian/pitch).",std::string(""));
   cr.add(regionName+".distribution_pitch","Name of pitch distribution function (isotropic).",std::string(""));
   cr.parse();
   cr.get(regionName+".particles_per_cell",N_particlesPerCell);
   cr.get(regionName+".number_density",numberDensity);
   cr.get(regionName+".distribution_energy",energyDistrib);
   cr.get(regionName+".distribution_pitch",pitchDistrib);
   
   // Check input values for sanity:
   if (N_particlesPerCell < 0) {
      simClasses.logger << "(DEFAULTINJECTOR) ERROR: Parameter '" << regionName << ".particles_per_cell' must have value >= 0 !" << std::endl;
      simClasses.logger << "\t Obtainer value is '" << N_particlesPerCell << "'." << std::endl << write;
      initialized = false;
   }
   if (numberDensity == defValue || numberDensity < 0.0) {
      simClasses.logger << "(DEFAULTINJECTOR) ERROR: Parameter '" << regionName << ".number_density' must have a finite value > 0.0 !" << std::endl;
      simClasses.logger << "\t Obtained value is '" << numberDensity << "'." << std::endl << write;
      initialized = false;
   }
   if (pitchDistrib == "isotropic") {
      pitchInjector = Isotropic;
   } else {
      simClasses.logger << "(DEFAULTINJECTOR) ERROR: Could not obtain pitch distribution function given in parameter '";
      simClasses.logger << regionName << ".distribution_pitch' !" << std::endl;
      simClasses.logger << "\t Please check configuration file(s)." << std::endl << write;
      initialized = false;
   }
   
   if (energyDistrib == "powerlaw") {
      energyInjector = PowerLaw;

      // Define additional parameters reguired for power law injection and parse values:
      std::string energyUnitString;
      cr.add(regionName+".spectral_index","Spectral index of power law distribution (float).",defValue);
      cr.add(regionName+".energy_unit","Unit in which 'energy_min' and 'energy_max' are given (eV/keV/MeV/GeV).","");
      cr.add(regionName+".energy_max","Maximum injection energy in 'energy_unit' units (float).",defValue);
      cr.add(regionName+".energy_min","Minimum injection energy in 'energy_unit' units (float).",defValue);
      cr.parse();
      cr.get(regionName+".spectral_index",sigma);
      cr.get(regionName+".energy_unit",energyUnitString);
      cr.get(regionName+".energy_min",minEnergy);
      cr.get(regionName+".energy_max",maxEnergy);
      
      // Check power law parameters for sanity:
      const Real energyFactor = simClasses.constants.getEnergyInSI(energyUnitString);
      if (sigma != sigma) {
	 simClasses.logger << "(DEFAULTINJECTOR) ERROR: Power law spectral index has not been set with parameter '";
	 simClasses.logger << regionName << ".spectral_index' !" << std::endl << write;
	 initialized = false;
      }
      if (energyFactor != energyFactor) {
	 simClasses.logger << "(DEFAULTINJECTOR) ERROR: Unknown energy units in parameter '";
	 simClasses.logger << regionName << ".energy_unit'. Must be (eV/keV/MeV/GeV) !" << std::endl << write;
	 initialized = false;
      }
      if (minEnergy != minEnergy) {
	 simClasses.logger << "(DEFAULTINJECTOR) ERROR: Minimum energy not given with parameter '";
	 simClasses.logger << regionName << ".energy_min' !" << std::endl << write;
	 initialized = false;
      }
      if (maxEnergy != maxEnergy) {
	 simClasses.logger << "(DEFAULTINJECTOR) ERROR: Maximum energy not given with parameter '";
	 simClasses.logger << regionName << ".energy_max' !" << std::endl << write;
	 initialized = false;
      }
      if (minEnergy > maxEnergy) {
	 simClasses.logger << "(DEFAULTINJECTOR) ERROR: Min. energy '" << regionName;
	 simClasses.logger << ".energy_min' > max. energy '" << regionName << ".energy_max' !" << std::endl << write;
	 initialized = false;
      }
      
      // Scale min,max energies to SI units:
      maxEnergy *= energyFactor;
      minEnergy *= energyFactor;
      minEnergy += std::numeric_limits<Real>::min();
      deltaEnergy = maxEnergy - minEnergy;
      
      // Calculate normalization constant, case sigma=-1 is a special case
      // as the integral of power law is a logarithm:
      if (sigma == -1.0) {
	 normalization = 1.0 / (log(maxEnergy) - log(minEnergy));
      } else {
	 normalization = (sigma+1.0) / (pow(maxEnergy,sigma+1.0) - pow(minEnergy,sigma+1.0));
      }
   } else if (energyDistrib == "gaussian") {
      energyInjector = Gaussian;
   } else {
      simClasses.logger << "(DEFAULTINJECTOR) ERROR: Could not obtain energy distribution function given in '";
      simClasses.logger << regionName << ".distribution_energy' !" << std::endl;
      simClasses.logger << "\t Please check configuration file(s)." << std::endl << write;
      initialized = false;
   }
   
   return initialized;
}

#endif
