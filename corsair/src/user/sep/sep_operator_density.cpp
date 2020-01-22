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

#include "sep_propagate.h"
#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_density.h"
#include "sep_particle_definition.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   
   OperatorDensity::OperatorDensity(int32_t order): OperatorAccumulationBase(1,order) {
      initialized = false;
   }
   
   OperatorDensity::~OperatorDensity() {
      finalize();
   }

   void OperatorDensity::accumulateBlock(pargrid::CellID blockLID,Real* RESTRICT accumArray,const std::vector<ParticleListBase*>& particleLists) {
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
	       
	 const pargrid::DataWrapper<sep::Particle<Real> >& wrapper = simClasses->pargrid.getUserDataDynamic<sep::Particle<Real> >(speciesDataID);
	 switch (OperatorAccumulationBase::getOrder()) {
	  case 0:
	    for (size_t p=0; p<wrapper.size(blockLID); ++p) {
	       Real pos[3];
	       pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD] - i_block*block::WIDTH_X + 1;
	       pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD] - j_block*block::WIDTH_Y + 1;
	       pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	       const Real weight = wrapper.data()[blockLID][p].state[particle::WEIGHT];
	       sep::accumScalarCentroidLogisticNGP_3D(array,pos,weight);
	    }
	    break;
	    
	  case 1:
	    for (size_t p=0; p<wrapper.size(blockLID); ++p) {
	       Real pos[3];
	       pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD] - i_block*block::WIDTH_X + 1;
	       pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD] - j_block*block::WIDTH_Y + 1;
	       pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	       const Real weight = wrapper.data()[blockLID][p].state[particle::WEIGHT];
	       sep::accumScalarCentroidLogisticCIC_3D(array,pos,weight);
	    }
	    break;
	    
	  case 2:
	    for (size_t p=0; p<wrapper.size(blockLID); ++p) {
	       Real pos[3];
	       pos[0] = wrapper.data()[blockLID][p].state[particle::XCRD] - i_block*block::WIDTH_X + 1;
	       pos[1] = wrapper.data()[blockLID][p].state[particle::YCRD] - j_block*block::WIDTH_Y + 1;
	       pos[2] = wrapper.data()[blockLID][p].state[particle::ZCRD] - k_block*block::WIDTH_Z + 1;
	       const Real weight = wrapper.data()[blockLID][p].state[particle::WEIGHT];
	       sep::accumScalarCentroidLogisticTSC_3D(array,pos,weight);
	    }
	    break;
	 }
      }
      
      // Add values from array to global accumulation array:
      block::addValues3D(*simClasses,blockLID,array,accumArray);
   }
   
   bool OperatorDensity::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   std::string OperatorDensity::getName() const {
      stringstream ss;
      ss << "sepOperatorDensity" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }

   uint32_t OperatorDensity::getNumberOfArrays() const {
      return 1;
   }
   
   std::string OperatorDensity::getOutputName(uint32_t arrayIndex) const {
      stringstream ss;
      ss << "density" << OperatorAccumulationBase::getOrder();
      return ss.str();
   }
   
   std::string OperatorDensity::getOutputUnits(uint32_t arrayIndex) const {
      return "1/m^3";
   }
   
   bool OperatorDensity::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP EM FIELDS) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }

      #if PROFILE_LEVEL > 0
         stringstream ss;
         ss << "density" << OperatorAccumulationBase::getOrder();
         profileName = ss.str();
         OperatorAccumulationBase::setProfileName(ss.str());
      #endif
      
      return initialized;
   }
   
   bool OperatorDensity::postProcessData(Real* array) {
      recalculateCellVolumes(*sim,*simClasses);

      // Divide values in each cell by spatial cell volumes:
      for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	    array[blockLID*block::SIZE+block::index(i,j,k)] /= cellVolume;
	 }
      }
      return true;
   }

   bool OperatorDensity::setAccumulatedArray(uint32_t arrayIndex) {
      return true;
   }
   
} // namespace sep
