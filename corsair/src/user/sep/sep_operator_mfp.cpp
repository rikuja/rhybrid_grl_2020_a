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
#include <sstream>

#include <linear_algebra.h>
#include <main.h>

#include "sep_propagate.h"
#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_mfp.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_guiding_center_theory.h"
#include "sep_particle_scatterer.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   const string PREFIX = "MeanFreePathOperator";
   
   OperatorMeanFreePath::OperatorMeanFreePath(int32_t order): OperatorDifferentialFlux(order) { }
   
   OperatorMeanFreePath::~OperatorMeanFreePath() { }

   void OperatorMeanFreePath::accumulateBlock(pargrid::CellID blockLID,Real* RESTRICT accumArray,
					      const std::vector<ParticleListBase*>& particleLists) {
      //return;
      // Create a temporary array for accumulating:
      const int32_t SIZE = block::SIZE_PLUS_ONE_LAYER;
      vector<Real> array(SIZE);
      for (int32_t i=0; i<SIZE; ++i) array[i] = 0.0;

      const int32_t NWL = simControl.N_wavelengthMeshCells+2;
      const int32_t WAVE_SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      vector<Real> parIntensity(WAVE_SIZE);

      // Load wave intensity:
      loadIntensity(blockLID,&(parIntensity[0]),+1);

      // Arrays for E,B,grad B values:
      Real E[3];
      Real B[3];
      Real gradB[9];
      
      // Get block global ID:
      int32_t i_block,j_block,k_block;
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

      // Iterate over all species accumulated to this channel:
      for (size_t s=0; s<instruments[currentInstrument].speciesNames.size(); ++s) {
	 // Get particle list for this species:
	 const std::string& speciesName = instruments[currentInstrument].speciesNames[s];
	 size_t speciesIndex = numeric_limits<size_t>::max();
	 for (size_t l=0; l<particleLists.size(); ++l) {
	    if (particleLists[l]->getSpeciesName() == speciesName) {
	       speciesIndex = l;
	       break;
	    }
	 }
	 if (speciesIndex == numeric_limits<size_t>::max()) continue;

	 // Get species struct:
	 const sep::Species& species = *reinterpret_cast<const sep::Species*>(particleLists[speciesIndex]->getSpecies());
	 
	 // Get fields at block centroid:
	 Real pos[4];
	 pos[0] = i_block*block::WIDTH_X + 0.5;
	 pos[1] = j_block*block::WIDTH_Y + 0.5;
	 pos[2] = k_block*block::WIDTH_Z + 0.5;
	 pos[3] = 0;
	 (*simControl.fieldsGetFields)(blockLID,corsair::getObjectWrapper().sim.t,pos,E,B,gradB);
	 const Real B_mag = vectorMagnitude<3>(B);
	 const Real omega = species.q_per_m * B_mag;
	 const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);

	 // Assume that energy is the average bin energy.
	 // Calculate maximum resonant wavelength:
	 Real energy = instruments[currentInstrument].minValues[currentChannel]
	             + instruments[currentInstrument].maxValues[currentChannel];
	 energy *= 0.5;
	 Real speed = sqrt(2*energy/species.mass);
	 Real maxResonantLambda = 2*M_PI/omega*speed;

	 pos[0] = 1.5;
	 pos[1] = 1.5;
	 pos[2] = 1.5;
	 
	 int32_t indices[4];
	 Real shapeFactors[12];
	 switch (simControl.order) {
	  case 0:
	    getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	    break;
	  case 1:
	    getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	    break;
	  case 2:
	    getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	    break;
	  default:
	    cerr << "ERROR: Unknown order of accuracy in sep_operator_mfp.cpp" << endl;
	    exit(1);
	    break;
	 }

	 // Integrate mean free path:
	 Real integral = 0.0;
	 Real pitchMin = -1.0;
	 Real pitchMax = +1.0;
	 int32_t N_pitch = 100;
	 Real d_pitch = (pitchMax-pitchMin) / N_pitch;

	 //lengthScale = speed/omega;

	 for (int32_t p=0; p<N_pitch; ++p) {
	    Real pitch = pitchMin + (p+0.5)*d_pitch;
	    Real D_mumu, d_D_mumu;
	    pos[3] = (*simControl.getLogicalWavelength)(maxResonantLambda * pitch);

	    evaluateDiffusionCoefficients(pitch,maxResonantLambda,diffConstant,pos,indices,shapeFactors,&(parIntensity[0]),D_mumu,d_D_mumu);
	    if (fabs(D_mumu) < 1e-15) continue;
	    integral += 3*speed/8 * (1-pitch*pitch)*(1-pitch*pitch)/D_mumu * d_pitch;
	 }
	 integral /= lengthScale;

	 switch (simControl.order) {
	  case 0:
	    accumScalarCentroidLogisticNGP_3D(&(array[0]),pos,integral);
	    break;
	  case 1:
	    accumScalarCentroidLogisticCIC_3D(&(array[0]),pos,integral);
	    break;
	  case 2:
	    accumScalarCentroidLogisticTSC_3D(&(array[0]),pos,integral);
	    break;
	  default:
	    cerr << "ERROR: Unknown order of accuracy in sep_operator_mfp.cpp" << endl;
	    exit(1);
	    break;
	 }
      }

      // Add values from array to global accumulation array:
      block::addValues3D(*simClasses,blockLID,&(array[0]),accumArray);
   }

   std::string OperatorMeanFreePath::getName() const {
      stringstream ss;
      ss << "MeanFreePath";
      return ss.str();
   }
   
   std::string OperatorMeanFreePath::getOutputName(uint32_t arrayIndex) const {
      const pair<size_t,size_t>& indices = instrumentIndices[arrayIndex];

      stringstream ss;
      ss << instruments[indices.first].name << "_mfp/";
      ss << instruments[indices.first].channelNames[indices.second];
      return ss.str();
   }
   
   std::string OperatorMeanFreePath::getOutputUnits(uint32_t arrayIndex) const {
      stringstream ss;
      ss << lengthScaleString;
      return ss.str();
   }
   
   bool OperatorMeanFreePath::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      simClasses.logger << "(SEP OP MFP) Starting initialization" << endl;

      // Init base class:
      if (OperatorDifferentialFlux::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP MFP) ERROR: Base class failed to initialize!" << endl << write;
	 initialized = false;
      }
      
      string gridBuilderName;
      cr.get("gridbuilder",gridBuilderName);
      cr.get(gridBuilderName+".input_units",lengthScaleString);
      lengthScale = simClasses.constants.getDistanceInSI(lengthScaleString);
      
      // Create a profiler section name for this DataOperator:
      #if PROFILE_LEVEL > 0
         stringstream ss;
         ss << "Mean Free Path" << OperatorAccumulationBase::getOrder();
         profileName = ss.str();
         OperatorAccumulationBase::setProfileName(ss.str());
      #endif

      // Write init status and exit:
      simClasses.logger << "(SEP OP MFP) Initialization complete, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      
      return initialized;
   }
   
   bool OperatorMeanFreePath::postProcessData(Real* array) {
      return true;
   }
   
   bool OperatorMeanFreePath::preProcessData(Real* array) {

      resetPitchLimits();
      //recalculateCellVolumes(*sim,*simClasses);
      
      #warning FIXME temporary solution, wave energy arrays are not synced
      bool success = true;
      if (simClasses->pargrid.startNeighbourExchange(sim->defaultStencilID,simControl.parAlfvenWaveEnergyDataID) == false) {
	 simClasses->logger << "(SEP OP MFP) Failed to start wave energy array sync" << std::endl << write;
	 success = false;
      }
      if (simClasses->pargrid.wait(sim->defaultStencilID,simControl.parAlfvenWaveEnergyDataID) == false) {
	 simClasses->logger << "(SEP OP MFP) Failed to sync wave energy arrays" << std::endl << write;
	 success = false;
      }

      return success;
   }
   
} // namespace sep
