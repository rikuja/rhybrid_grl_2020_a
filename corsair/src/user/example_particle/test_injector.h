/** This file is part of Corsair simulation.
 *
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

#ifndef TEST_INJECTOR_H
#define TEST_INJECTOR_H

#include <cstdlib>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>

template<class FIELD,class SPECIES,class PARTICLE>
class TestInjector {
 public:
   TestInjector();
   
   bool finalize();
   bool inject(const SPECIES& species,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper);
   bool readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		       const std::string& regionName,const SPECIES& species);
   
 private:
   Real initialEnergy;
   Real initialPitch;
   bool initialized;
   bool injected;
   int N_injected;
   Real origin[3];
   Simulation* sim;
   SimulationClasses* simClasses;
};

template<class FIELD,class SPECIES,class PARTICLE> inline
TestInjector<FIELD,SPECIES,PARTICLE>::TestInjector() {
   injected = false;
   initialized = false;
   sim = NULL;
   simClasses = NULL;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool TestInjector<FIELD,SPECIES,PARTICLE>::finalize() {
   initialized = false;
   return true;
}

template<class FIELD,class SPECIES,class PARTICLE> inline
bool TestInjector<FIELD,SPECIES,PARTICLE>::inject(const SPECIES& species,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper) {
   // If particles have already been injected return immediately:
   if (initialized == false) return initialized;
   if (injected == true) return true;
   
   const double* blockCrds = getBlockCoordinateArray(*sim,*simClasses);
   Real blockSize[3];
   
   for (pargrid::CellID block=0; block<simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      // Get the size of the block:
      getBlockSize(*simClasses,*sim,block,blockSize);
      
      // Test if injection point is inside this block:
      if (origin[0] < blockCrds[3*block+0] || origin[0] > (blockCrds[3*block+0]+blockSize[0])) continue;
      if (origin[1] < blockCrds[3*block+1] || origin[1] > (blockCrds[3*block+1]+blockSize[1])) continue;
      if (origin[2] < blockCrds[3*block+2] || origin[2] > (blockCrds[3*block+2]+blockSize[2])) continue;
      
      // Resize the particle list of this block to have space
      // for N_injected extra particles:
      const pargrid::ArraySizetype currentSize = wrapper.size(block);
      wrapper.resize(block,currentSize+N_injected);
      N_particles[block] += N_injected;

      // Inject N_particles particles to this block. New particles
      // should be added to the end of the particle list:
      pargrid::ArraySizetype counter = currentSize;
      PARTICLE* const particles = wrapper.data()[block];
      Real E[3];
      Real B[3];
      Real coords[3];
      Real V_drift[3];
      for (int i=0; i<N_injected; ++i) {
	 // Particle coordinates are relative to block bottom lower left corner.
	 // First three data fields of class PARTICLE are reserved for coordinates:
	 particles[counter][XPOS] = simClasses->random.uniform()*blockSize[0];
	 particles[counter][YPOS] = simClasses->random.uniform()*blockSize[1];
	 particles[counter][ZPOS] = simClasses->random.uniform()*blockSize[2];
	 
	 // Fetch field values from correct position:
	 coords[XPOS] = particles[counter][XPOS] + blockCrds[3*block+0];
	 coords[YPOS] = particles[counter][YPOS] + blockCrds[3*block+0];
	 coords[ZPOS] = particles[counter][ZPOS] + blockCrds[3*block+0];	 
	 FIELD::getFields(block,sim->t,coords,E,B);
	 
	 // Calculate energy in drift-free frame:
	 crossProduct(E,B,V_drift);
	 const Real B_mag2 = vectorMagnitude2<3>(B);
	 for (int i=0; i<3; ++i) V_drift[i] /= B_mag2;
	 const Real energy_dot = initialEnergy - 0.5*species.m*vectorMagnitude2<3>(V_drift);
	 const Real B_mag_inverse = 1.0 / sqrt(B_mag2);
	 
	 particles[counter][particle::VPAR]   = initialPitch*sqrt(2.0*energy_dot/species.m);
	 particles[counter][particle::MU]     = energy_dot*B_mag_inverse*(1.0-initialPitch*initialPitch);
	 particles[counter][particle::WEIGHT] = 1.0;
	 
	 ++counter;
      }
   }
   injected = true;
   return true;
}

/** Read input parameters required by this particle injector from the given
 * config file region.
 * @param sim Struct containing generic simulation parameters.
 * @param simClasses Struct containing generic simulation classes.
 * @param cr Config file reader.
 * @param configName Name of config file region that contains injector's parameters.
 * @return If true, injector read all required values and initialized successfully.*/
template<class FIELD,class SPECIES,class PARTICLE> inline
bool TestInjector<FIELD,SPECIES,PARTICLE>::readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
					 const std::string& configName,const SPECIES& species) {
   bool success = true;
   this->sim = &sim;
   this->simClasses = &simClasses;
   
   const Real defValue = std::numeric_limits<Real>::infinity();

   // Define and read parameters from config file:
   cr.add(configName+".particles_per_cell","Number of particles to inject (int).",(int)0);
   cr.add(configName+".pitch","Value of pitch used in injection, -1 <= pitch <= 1 (float).",defValue);
   cr.add(configName+".energy","Value of energy used in injection, in given energy units (float).",defValue);
   cr.add(configName+".energy_unit","Energy unit (eV/keV/MeV/GeV) (string).",std::string(""));
   cr.add(configName+".origin_x","X-coordinate within the injection cell (float).",defValue);
   cr.add(configName+".origin_y","Y-coordinate within the injection cell (float).",defValue);
   cr.add(configName+".origin_z","Z-coordinate within the injection cell (float).",defValue);
   cr.parse();
   std::string energyUnit;
   cr.get(configName+".particles_per_cell",N_injected);
   cr.get(configName+".pitch",initialPitch);
   cr.get(configName+".energy",initialEnergy);
   cr.get(configName+".energy_unit",energyUnit);
   cr.get(configName+".origin_x",origin[0]);
   cr.get(configName+".origin_y",origin[1]);
   cr.get(configName+".origin_z",origin[2]);

   // Sanity check on input parameters:
   if (N_injected < 0) {
      simClasses.logger << "(TESTINJECTOR) ERROR: Number of injected particles was not given with parameter '";
      simClasses.logger << configName+".particles' " << std::endl;
      simClasses.logger << "\t or it has a negative value!" << std::endl << write;
      success = false;
   }
   if ((initialPitch < -1.0 || initialPitch > 1.0) || initialPitch == defValue) {
      simClasses.logger << "(TESTINJECTOR) ERROR: Injection pitch was not given with parameter '" << configName+".pitch' " << std::endl;
      simClasses.logger << "\t or it has a value outside range [-1,+1] !" << std::endl << write;
      success = false;
   }
   if (initialEnergy < 0.0 || initialEnergy == defValue) {
      simClasses.logger << "(TESTINJECTOR) ERROR: Injection energy was not given with parameter '" << configName+".energy' " << std::endl;
      simClasses.logger << "\t or it has a negative value !" << std::endl << write;
      success = false;
   }
   for (int i=0; i<3; ++i) if (origin[i] == defValue) {
      simClasses.logger << "(TESTINJECTOR) ERROR: Position within the injection cell was not given with parameters" << std::endl;
      simClasses.logger << "\t '" << configName+".origin_x', '" << configName+".origin_y', and '" << configName+".origin_z' !" << std::endl << write;
      success = false;
   }
   const Real energyFactor = simClasses.constants.getEnergyInSI(energyUnit);
   if (energyFactor == simClasses.constants.notFound()) {
      simClasses.logger << "(TESTINJECTOR) ERROR: Unknown energy unit given in parameter '" << configName+".unit' ";
      simClasses.logger << "or it was not defined at all!" << std::endl << write;
      success = false;
   }
   initialEnergy *= energyFactor;

   initialized = success;
   return initialized;
}

#endif
