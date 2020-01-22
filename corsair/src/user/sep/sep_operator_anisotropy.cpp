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

#include <linear_algebra.h>
#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_anisotropy.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   extern sep::FieldsContainer fieldsContainer;
   
   OperatorAnisotropy::OperatorAnisotropy(int32_t order): OperatorAccumulationBase(1,order) {
      initialized = false;
   }
   
   OperatorAnisotropy::~OperatorAnisotropy() {
      finalize();
   }

   void OperatorAnisotropy::accumulateBlock(pargrid::CellID blockLID,Real* RESTRICT accumArray,const std::vector<ParticleListBase*>& particleLists) {
      // Clear temp array:
      const int32_t SIZE = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real array[SIZE];
      for (int32_t i=0; i<SIZE; ++i) array[i] = 0.0;
      
      // Get block global ID:
      int32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      for (size_t s=0; s<particleLists.size(); ++s) {
	 // Skip non-particle species:
	 if (particleLists[s]->getSpeciesType() != simControl.particleSpeciesTypename) continue;
	 
	 pargrid::DataID speciesDataID;
	 if (particleLists[s]->getParticles(speciesDataID) == false) {
	    continue;
	 }

	 // Get Species struct:
	 const sep::Species* species = reinterpret_cast<const sep::Species*>(particleLists[s]->getSpecies());
	 
	 Real B[3];
	 Real V_wave[3];
	 Real dV_wave;
	 const pargrid::DataWrapper<sep::Particle<Real> >& wrapper = simClasses->pargrid.getUserDataDynamic<sep::Particle<Real> >(speciesDataID);
	 for (size_t p=0; p<wrapper.size(blockLID); ++p) {
	    // Get fields at particle position:
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,wrapper.data()[blockLID][p].state,B,V_wave,dV_wave,V_alfven,+1);
	    
	    // Calculate pitch (cosine of pitch angle):
	    const Real B_mag   = vectorMagnitude<3>(B);
	    const Real V_par   = wrapper.data()[blockLID][p].state[particle::V_PAR];

	    const Real V_gyro2 = 2.0*wrapper.data()[blockLID][p].state[particle::MU]*B_mag/species->mass;
	    const Real V_mag   = sqrt(V_par*V_par + V_gyro2);
	    const Real mu      = V_par / V_mag;

	    // Calculate offset relative to accumulation block's corner:
	    Real pos[3];
	    pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD] - i_block*block::WIDTH_X + 1;
	    pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD] - j_block*block::WIDTH_Y + 1;
	    pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	    
	    // Accumulate pitch x particle weight:
	    const Real value = wrapper.data()[blockLID][p].state[particle::WEIGHT] * mu;
	    switch (OperatorAccumulationBase::getOrder()) {
	     case 0:
	       sep::accumScalarCentroidLogisticNGP_3D(array,pos,value);
	       break;
	     case 1:
	       sep::accumScalarCentroidLogisticCIC_3D(array,pos,value);
	       break;
	     case 2:
	       sep::accumScalarCentroidLogisticTSC_3D(array,pos,value);
	       break;
	    }
	 } // for (size_t p=0; p<wrapper.size(blockLID); ++p)
      } // for (size_t s=0; s<particleLists.size(); ++s)
      
      // Add values from array to global accumulation array:
      block::addValues3D(*simClasses,blockLID,array,accumArray);
   }
   
   bool OperatorAnisotropy::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   std::string OperatorAnisotropy::getName() const {
      stringstream ss;
      ss << "sepOperatorAnisotropy" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }
   
   uint32_t OperatorAnisotropy::getNumberOfArrays() const {
      return 1;
   }
   
   std::string OperatorAnisotropy::getOutputName(uint32_t arrayIndex) const {
      stringstream ss;
      ss << "anisotropy" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }
   
   std::string OperatorAnisotropy::getOutputUnits(uint32_t arrayIndex) const {
      return "";
   }
   
   bool OperatorAnisotropy::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP ANISOTROPY) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }

      #if PROFILE_LEVEL > 0
         stringstream ss;
         ss << "Anisotropy" << OperatorAccumulationBase::getOrder();
         profileName = ss.str();
         OperatorAccumulationBase::setProfileName(ss.str());
      #endif
      
      return initialized;
   }
   
   bool OperatorAnisotropy::setAccumulatedArray(uint32_t arrayIndex) {
      return true;
   }
   
} // namespace sep
