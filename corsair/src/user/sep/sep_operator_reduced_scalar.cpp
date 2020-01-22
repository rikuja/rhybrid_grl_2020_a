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
#include <map>

#include "sep_accumulation_stretched.h"
#include "sep_simcontrol.h"
#include "sep_operator_reduced_scalar.h"
#include "sep_particle_definition.h"

using namespace std;

namespace sep {
   extern sep::SimControl simControl;
   const string PREFIX = "ReducedScalar";
   
   OperatorReducedScalar::OperatorReducedScalar(): DataOperator() {
      initialized = false;
      sliceIndex = numeric_limits<uint32_t>::max();
      variablesSearched = false;
      
      #if PROFILE_LEVEL > 0
         totalTimeID = -1;
      #endif
   }
   
   OperatorReducedScalar::~OperatorReducedScalar() {
      finalize();
   }
   
   bool OperatorReducedScalar::addConfigFileItems(ConfigReader& cr) {
      cr.add(PREFIX+".slice_index","Which index of wavelength mesh is written (int).",numeric_limits<uint32_t>::max());
      cr.addComposed(PREFIX+".variable","Names of variables written to output file (string).");
      return true;
   }

   bool OperatorReducedScalar::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   std::string OperatorReducedScalar::getName() const {
      return "sepReducedScalar";
   }
   
   bool OperatorReducedScalar::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP REDUCED SCALAR) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }
      
      if (addConfigFileItems(cr) == false) {
	 simClasses.logger << "(SEP OP REDUCED SCALAR) ERROR: Failed to add config file items" << endl << write;
	 initialized = false;
	 return initialized;
      }
      
      cr.parse();
      cr.get(PREFIX+".slice_index",sliceIndex);
      cr.get(PREFIX+".variable",variableNames);
      
      if (sliceIndex == numeric_limits<uint32_t>::max()) {
	 simClasses.logger << "(SEP OP REDUCED SCALAR) ERROR: Could not find variable '" << PREFIX+".slice_index' from config file." << endl << write;
	 initialized = false;
      }      
      if (sliceIndex < 0 || sliceIndex >= simControl.N_wavelengthMeshCells) {
	 simClasses.logger << "(SEP OP REDUCED SCALAR) ERROR: Variable '" << PREFIX+".slice_index' has invalid value" << endl;
	 simClasses.logger << "\t Read value is '" << sliceIndex << "', it should be >= 0 and < " << simControl.N_wavelengthMeshCells << endl << write;
	 initialized = false;
      }
      
      #if PROFILE_LEVEL > 0
         profileName = "Reduced Scalar";
      #endif
      
      return initialized;
   }
   
   bool OperatorReducedScalar::searchVariableDataIDs() {
      bool success = true;
      
      for (size_t i=0; i<variableNames.size(); ++i) {
	 // Get ParGrid data ID of given variable:
	 bool arrayIsDynamic = true;
	 const pargrid::DataID dataID = simClasses->pargrid.getUserDataID(variableNames[i],arrayIsDynamic);
	 
	 // Check that array is static:
	 if (arrayIsDynamic == true) {
	    simClasses->logger << "(SEP OP REDUCED SCALAR) Skipping variable '" << variableNames[i] << "' because it is dynamic." << endl << write;
	    success = false;
	    continue;
	 }
	 
	 // Check that array contains a scalar in 4D wavelength mesh:
	 const size_t N_elements = simClasses->pargrid.getUserDataStaticElements(dataID);
	 if (N_elements != (size_t)block::SIZE*simControl.N_wavelengthMeshCells) {
	    simClasses->logger << "(SEP OP REDUCED SCALAR) Skipping variable '" << variableNames[i] << "' because ";
	    simClasses->logger << "it is not scalar." << endl << write;
	    success = false;
	    continue;
	 }
	 
	 variableOutputNames.push_back(variableNames[i]);
	 variableDataIDs.push_back(dataID);
      }
      return success;
   }
   
   bool OperatorReducedScalar::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start(profileName,totalTimeID);
      #endif

      if (variablesSearched == false) {
	 success = searchVariableDataIDs();
	 variablesSearched = true;
      }
      
      // Create map for XML tag attributes and write mesh geometry to it:
      map<string,string> xmlAttributes;
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 xmlAttributes["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
       case sep::CARTESIAN:
	 xmlAttributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
       case sep::CYLINDRICAL:
	 xmlAttributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 break;
       case sep::SPHERICAL:
	 xmlAttributes["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	 break;
       default:
	 xmlAttributes["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
      }
      
      // Create a buffer for data:
      const size_t SIZE_BUFFER = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
      const size_t NWL = simControl.N_wavelengthMeshCells;
      Real* buffer = new Real[SIZE_BUFFER];
      
      // Reduce all scalar quantities and write to file:
      for (size_t scalar=0; scalar<variableDataIDs.size(); ++scalar) {
	 // Get pointer to ParGrid data array:
	 const Real* const array = simClasses->pargrid.getUserDataStatic<Real>(variableDataIDs[scalar]);
	 if (array == NULL) {
	    simClasses->logger << "(SEP OP REDUCED SCALAR) WARNING: Skipping variable '" << variableOutputNames[scalar];
	    simClasses->logger << "' because I got NULL pointer from ParGrid" << endl << write;
	    continue;
	 }
	 
	 // Copy data to buffer:
	 for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	    const Real* const ptr = array + blockLID*block::SIZE*NWL;
	    for (int k=0; k<block::WIDTH_X; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	       buffer[blockLID*block::SIZE+block::index(i,j,k)] = ptr[block::index(i,j,k)*NWL + sliceIndex];
	    }
	 }
	 
	 // Write buffer to output file:
	 xmlAttributes["name"] = spatMeshName + '/' + variableOutputNames[scalar];
	 xmlAttributes["mesh"] = spatMeshName;
	 xmlAttributes["centering"] = "zone";
	 if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,SIZE_BUFFER,1,buffer) == false) success = false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      delete [] buffer; buffer = NULL;
      return success;
   }
      
   
} // namespace sep
