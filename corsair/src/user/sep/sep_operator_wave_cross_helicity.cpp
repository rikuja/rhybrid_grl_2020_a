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

#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_operator_wave_cross_helicity.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   
   OperatorWaveCrossHelicity::OperatorWaveCrossHelicity(): DataOperator() {
      initialized = false;
      
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorWaveCrossHelicity::~OperatorWaveCrossHelicity() {
      finalize();
   }
   
   bool OperatorWaveCrossHelicity::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }
   
   bool OperatorWaveCrossHelicity::calculateCrossHelicity(pargrid::CellID blockLID,Real* crossHelicity) {
      // Get pointers to arrays containing wave energies:
      Real* antiparWaveEnergy = simClasses->pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      Real* parWaveEnergy     = simClasses->pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      if (antiparWaveEnergy == NULL) return false;
      if (parWaveEnergy == NULL) return false;
      
      for (int32_t cell=0; cell<block::SIZE; ++cell) {
	 // Calculate total antiparallel wave energy:
	 Real antiparEnergy = 0.0;
	 for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	    antiparEnergy += antiparWaveEnergy[blockLID*block::SIZE*simControl.N_wavelengthMeshCells + cell*block::SIZE + l];
	 }
      
	 // Calculate total parallel wave energy:
	 Real parEnergy = 0.0;
	 for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	    parEnergy += parWaveEnergy[blockLID*block::SIZE*simControl.N_wavelengthMeshCells + cell*block::SIZE + l];
	 }
	 
	 // Calculate normalized cross helicity:
	 crossHelicity[cell] = 2.0*parEnergy / (parEnergy + antiparEnergy) - 1.0;
      }

      return true;
   }

   std::string OperatorWaveCrossHelicity::getName() const {return "sepOperatorWaveCrossHelicity";}
   
   bool OperatorWaveCrossHelicity::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP WAVE CROSS HELICITY) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }
   
      return initialized;
   }

   bool OperatorWaveCrossHelicity::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start("Wave Cross Helicity",profTotalTime);
      #endif

      // Calculate cross helicity in each block:
      Real* crossHelicity = new Real[simClasses->pargrid.getNumberOfLocalCells()*block::SIZE];      
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
      for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	 Real* ptr = crossHelicity + blockLID*block::SIZE;
	 if (calculateCrossHelicity(blockLID,ptr) == false) {
	    simClasses->logger << "(SEP OP WAVE CROSS HELICITY) ERROR: Failed to calculate cross helicity in block LID#";
	    simClasses->logger << blockLID << " GID#" << simClasses->pargrid.getGlobalIDs()[blockLID] << std::endl << write;
	    success = false;
	    break;
	 }
      }
      
      // Write out data:
      map<std::string,std::string> attribs;
      attribs["name"] = spatMeshName + "/waves/crossHelicity";
      attribs["mesh"] = spatMeshName;
      attribs["unit"] = "";
      attribs["centering"] = "zone";

      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
       case sep::CARTESIAN:
	 attribs["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
       case sep::CYLINDRICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 break;
       case sep::SPHERICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	 break;
       default:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
      }
      
      const uint64_t arraySize = N_localBlocks*block::SIZE;
      const uint64_t vectorSize = 1;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,crossHelicity) == false) success = false;
      
      delete [] crossHelicity; crossHelicity = NULL;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      
      return success;
   }
      
   
} // namespace sep
