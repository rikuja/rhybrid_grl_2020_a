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
#include <sstream>

#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_particle_energy.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_guiding_center_theory.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   extern sep::FieldsContainer fieldsContainer;
   
   OperatorParticleEnergy::OperatorParticleEnergy(int32_t order): OperatorAccumulationBase(1,order) {
      initialized = false;
   }
   
   OperatorParticleEnergy::~OperatorParticleEnergy() {
      finalize();
   }

   void OperatorParticleEnergy::accumulateBlock(pargrid::CellID blockLID,Real* RESTRICT accumArray,const std::vector<ParticleListBase*>& particleLists) {
      // Clear temp array:
      const int32_t SIZE = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real array[SIZE];
      for (int32_t i=0; i<SIZE; ++i) array[i] = 0.0;
      
      // Get block global ID:
      int32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      // Arrays for E, B, nabla B:
      Real E[3];
      Real B[3];
      Real gradB[9];
      
      for (size_t s=0; s<particleLists.size(); ++s) {
	 // Skip non-particle species:
	 if (particleLists[s]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 
	 pargrid::DataID speciesDataID;
	 if (particleLists[s]->getParticles(speciesDataID) == false) {
	    continue;
	 }

	 const sep::Species* species = reinterpret_cast<const sep::Species*>(particleLists[s]->getSpecies());
	 const pargrid::DataWrapper<sep::Particle<Real> >& wrapper = simClasses->pargrid.getUserDataDynamic<sep::Particle<Real> >(speciesDataID);
	 for (size_t p=0; p<wrapper.size(blockLID); ++p) {
	    (*simControl.fieldsGetFields)(blockLID,sim->t,wrapper.data()[blockLID][p].state,E,B,gradB);
	    
	    const Real B_mag = vectorMagnitude<3>(B);
	    
	    // Calculate gyro energy:
	    const Real U_gyro = wrapper.data()[blockLID][p].state[sep::particle::MU] * B_mag;
	    
	    // Calculate drift velocity:
	    Real V_drift[3];
	    sep::calculateDriftVelocity(*species,wrapper.data()[blockLID][p],E,B,gradB,V_drift);

	    Real energy = wrapper.data()[blockLID][p].state[sep::particle::V_PAR]*wrapper.data()[blockLID][p].state[sep::particle::V_PAR]
	                + vectorMagnitude2<3>(V_drift);
	    energy *= 0.5*species->mass;
	    energy += U_gyro;
	    energy *= wrapper.data()[blockLID][p].state[sep::particle::WEIGHT];
	    energy /= (constants::CHARGE_ELEMENTARY*1000.0);
	    
	    Real pos[3];
	    pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD] - i_block*block::WIDTH_X + 1;
	    pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD] - j_block*block::WIDTH_Y + 1;
	    pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	    
	    switch (OperatorAccumulationBase::getOrder()) {
	     case 0:
	       sep::accumScalarCentroidLogisticNGP_3D(array,pos,energy);
	       break;
	     case 1:
	       sep::accumScalarCentroidLogisticCIC_3D(array,pos,energy);
	       break;
	     case 2:
	       sep::accumScalarCentroidLogisticTSC_3D(array,pos,energy);
	       break;
	    }
	 }
      }
      
      // Add values from array to global accumulation array:
      block::addValues3D(*simClasses,blockLID,array,accumArray);
   }
   
   bool OperatorParticleEnergy::finalize() {
      if (initialized == false) return true;
      initialized = OperatorAccumulationBase::finalize();
      return true;
   }

   std::string OperatorParticleEnergy::getName() const {
      stringstream ss;
      ss << "sepOperatorParticleEnergy" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }
   
   uint32_t OperatorParticleEnergy::getNumberOfArrays() const {
      return 1;
   }
   
   std::string OperatorParticleEnergy::getOutputName(uint32_t arrayIndex) const {
      stringstream ss;
      ss << "TotalParticleEnergy" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }
   
   std::string OperatorParticleEnergy::getOutputUnits(uint32_t arrayIndex) const {
      return "keV";
   }
   
   bool OperatorParticleEnergy::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;

      // Init base class:
      if (OperatorAccumulationBase::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP PARTICLE ENERGY) ERROR: OperatorAccumulationBase failed to initialize!" << endl << write;
	 initialized = false;
      }

      #if PROFILE_LEVEL > 0
         stringstream ss;
         ss << "Energies" << OperatorAccumulationBase::getOrder();
         OperatorAccumulationBase::setProfileName(ss.str());
      #endif
      
      return initialized;
   }
   
   bool OperatorParticleEnergy::setAccumulatedArray(uint32_t arrayIndex) {
      return true;
   }
   
} // namespace sep
