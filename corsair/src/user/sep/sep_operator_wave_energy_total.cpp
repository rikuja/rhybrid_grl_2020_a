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
#include "sep_operator_wave_energy_total.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   
   /** Default constructor, only calls base class constructor.*/
   OperatorWaveEnergyTotal::OperatorWaveEnergyTotal(): DataOperator() {
      initialized = false;
      
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   /** Destructor, calls OperatorWaveEnergyTotal::finalize member function.*/
   OperatorWaveEnergyTotal::~OperatorWaveEnergyTotal() {
      finalize();
   }
   
   /** Finalizer function, does not do anything.
    * @return If true, class was finalized successfully.*/
   bool OperatorWaveEnergyTotal::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   /** Return name of this DataOperator, must be unique, i.e., no other DataOperator 
    * can have the same name. This name is used when registering and requesting 
    * DataOperators from a container class.
    * @return Name of this DataOperator.*/
   std::string OperatorWaveEnergyTotal::getName() const {return "sepOperatorWaveEnergyTotal";}

   /** Class initializer, calls base class initializer DataOperator::initialize.
    * @param cr Configuration file reader.
    * @param sim Struct containing variables that control execution of simulation.
    * @param simClasses Struct containing generic use simulation classes.
    * @return If true, class initialized successfully and is ready for use.*/
   bool OperatorWaveEnergyTotal::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP WAVE ENERGY TOTAL) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }
   
      return initialized;
   }
   
   /** Write an array to spatial mesh to output file.
    * @param spatMeshName Name of spatial mesh.
    * @param name Name of the written array.
    * @param buffer Array containing values to be written.
    * @return If true, data was written successfully.*/
   bool OperatorWaveEnergyTotal::writeArray(const std::string& spatMeshName,const std::string& name,const Real* buffer) {
      bool success = true;
      
      map<std::string,std::string> attribs;
      attribs["name"] = name;
      attribs["mesh"] = spatMeshName;
      attribs["unit"] = "keV";
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

      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
      const uint64_t arraySize = N_localBlocks*block::SIZE;
      const uint64_t vectorSize = 1;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,buffer) == false) success = false;
      
      return success;
   }
   
   bool OperatorWaveEnergyTotal::writeWaveGrowth(const std::string& spatMeshName,const std::string& name,pargrid::DataID dataID) {
      bool success = true;
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
      Real* waveValues = simClasses->pargrid.getUserDataStatic<Real>(dataID);
      if (waveValues == NULL) return false;
      
      Real* buffer = new Real[N_localBlocks*block::SIZE];
      
      for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	 // Iterate over all spatial cells in block:
	 for (int32_t cell=0; cell<block::SIZE; ++cell) {
	    // Sum up total wave energy:
	    buffer[blockLID*block::SIZE + cell*block::SIZE] = 0.0;
	    for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	       buffer[blockLID*block::SIZE + cell]
		 += waveValues[blockLID*block::SIZE*simControl.N_wavelengthMeshCells + cell*block::SIZE+l];
	    }
	    
	    buffer[blockLID*block::SIZE + cell*block::SIZE] /= sim->dt;
	 }
      }
      
      if (writeArray(spatMeshName,spatMeshName+"/Waves/"+name,buffer) == false) success = false;
      delete [] buffer; buffer = NULL;
      return success;
   }

   /** Calculates total energy in given wave mode (parallel or antiparallel) and 
    * writes data to output file. Wave mode is defined by giving its ParGrid DataID.
    * @param spatMeshName Name of spatial mesh.
    * @param name Name of array written to output file, must be valid for Visit.
    * @param dataID ParGrid DataID of array containing wave energy in each wavelength bin.
    * @return If true, total wave energy densities were successfully calculated and written to output file.*/
   bool OperatorWaveEnergyTotal::writeWaveMode(const std::string& spatMeshName,const std::string& name,pargrid::DataID dataID) {
      bool success = true;
      
      // Get local cell IDs:
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();
      
      // Get pointer to array containing wave energies:
      Real* waveEnergy = simClasses->pargrid.getUserDataStatic<Real>(dataID);
      if (waveEnergy == NULL) return false;

      // Allocate buffer for total wave energy:
      Real* buffer = new Real[N_localBlocks*block::SIZE];
      
      // Iterate over all local blocks:
      for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	 // Iterate over all spatial cells in block:
	 for (int32_t cell=0; cell<block::SIZE; ++cell) {
	    // Sum up total wave energy:
	    buffer[blockLID*block::SIZE + cell*block::SIZE] = 0.0;
	    
	    for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	       buffer[blockLID*block::SIZE + cell] 
		 += waveEnergy[blockLID*block::SIZE*simControl.N_wavelengthMeshCells + cell*block::SIZE+l];
	    }
	    buffer[blockLID*block::SIZE + cell*block::SIZE] /= simControl.cellVolumes[blockLID*block::SIZE + cell];
	 }
      }

      if (writeArray(spatMeshName,spatMeshName+"/Waves/"+name,buffer) == false) success = false;
      delete [] buffer; buffer = NULL;
      return success;
   }
   
   bool OperatorWaveEnergyTotal::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start("Total Wave Energy",profTotalTime);
      #endif

      // Recalculate spatial cell volumes, ensures that the volumes 
      // are correct after mesh has been repartitioned:
      recalculateCellVolumes(*sim,*simClasses);

      if (writeWaveGrowth(spatMeshName,"WaveGrowth_Antiparallel",simControl.antiparWaveGrowthDataID) == false) {
	 simClasses->logger << "(SEP OPERATOR WAVE ENERGY TOTAL) ERROR: Failed to write parallel Alfven wave energy" << endl << write;
	 success = false;
      }
      if (writeWaveGrowth(spatMeshName,"WaveGrowth_Parallel",simControl.parWaveGrowthDataID) == false) {
	 simClasses->logger << "(SEP OPERATOR WAVE ENERGY TOTAL) ERROR: Failed to write parallel Alfven wave energy" << endl << write;
	 success = false;
      }

      // Calculate energy densities in antiparallel propagating Alfven waves (sum of L and R helicities):
      if (writeWaveMode(spatMeshName,"EnergyDens_Antiparallel",simControl.antiparAlfvenWaveEnergyDataID) == false) {
	 simClasses->logger << "(SEP OPERATOR WAVE ENERGY TOTAL) ERROR: Failed to write antiparallel Alfven wave energy" << endl << write;
	 success = false;
      }
      
      // Calculate energy densities is parallel propagating Alfven waves (sum of L and R helicities):
      if (writeWaveMode(spatMeshName,"EnergyDens_Parallel",simControl.parAlfvenWaveEnergyDataID) == false) {
	 simClasses->logger << "(SEP OPERATOR WAVE ENERGY TOTAL) ERROR: Failed to write parallel Alfven wave energy" << endl << write;
	 success = false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

} // namespace sep
