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
#include <map>
#include <vector>

#include "sep_simcontrol.h"
#include "sep_operator_wavelength_mesh_scalar.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const string PREFIX = "WavelengthMeshVariables";
   
   WavelengthMeshScalarOP::WavelengthMeshScalarOP(): SpatialSliceOP() { 
      variablesSearched = false;
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   WavelengthMeshScalarOP::~WavelengthMeshScalarOP() { }

   bool WavelengthMeshScalarOP::addConfigFileItems(ConfigReader& cr) { 
      cr.addComposed(PREFIX+".variable","Names of variables written to output file (string).");
      return true;
   }
   
   bool WavelengthMeshScalarOP::finalize() { 
      return true;
   }
   
   std::string WavelengthMeshScalarOP::getName() const { 
      return "WLM_ScalarOP";
   }
   
   bool WavelengthMeshScalarOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) { 
      // Init base class:
      if (SpatialSliceOP::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP WLM SCALAR) Base class failed to init" << endl << write;
	 return false;
      }
      
      bool success = true;

      addConfigFileItems(cr);
      cr.parse();
      cr.get(PREFIX+".variable",variableNames);
      
      return success;
   }
   
   bool WavelengthMeshScalarOP::searchVariableDataIDs() {
      for (size_t i=0; i<variableNames.size(); ++i) {
	 //cerr << "Checking " << variableNames[i] << endl;
	 
	 // Get ParGrid data ID of given variable:
	 bool arrayIsDynamic = true;
	 const pargrid::DataID dataID = simClasses->pargrid.getUserDataID(variableNames[i],arrayIsDynamic);
	 
	 // Check that array is static:
	 if (arrayIsDynamic == true) {
	    simClasses->logger << "(SEP OP WLM SCALAR) Skipping variable '" << variableNames[i] << "' because it is dynamic." << endl << write;
	    continue;
	 }
	 
	 // Check that array contains a scalar:
	 const size_t N_elements = simClasses->pargrid.getUserDataStaticElements(dataID);
	 if (N_elements != (size_t)block::SIZE*simControl.N_wavelengthMeshCells) {
	    simClasses->logger << "(SEP OP WLM SCALAR) Skipping variable '" << variableNames[i] << "' because ";
	    simClasses->logger << "it is not scalar." << endl << write;
	    continue;
	 }
	 
	 variableOutputNames.push_back(variableNames[i]);
	 variableDataIDs.push_back(dataID);
      }
      
      // Clear vector variableNames:
      vector<string> varDummy;
      variableNames.swap(varDummy);
      
      return true;
   }
   
   bool WavelengthMeshScalarOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) { 
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start("WLM Scalars",profTotalTime);
      #endif
      
      if (variablesSearched == false) {
	 success = searchVariableDataIDs();
	 variablesSearched = true;
      }

      // Write scalars on all slices:
      for (size_t slice=0; slice<getNumberOfSlices(); ++slice) {      
	 stringstream ss;
	 ss << "WaveMesh" << slice;
	 const string meshName = ss.str();

	 // Get IDs of blocks in slice:
	 prepareSlice(slice,simControl.N_wavelengthMeshCells);
	 size_t counter = 0;
	 vector<pargrid::CellID> blockLIDs;
	 vector<pargrid::CellID> newBlockGIDs;
	 getAcceptedBlocks(slice,newBlockGIDs,blockLIDs);

	 int32_t sliceIndex = numeric_limits<int32_t>::max();
	 for (size_t b=0; b<blockLIDs.size(); ++b) {
	    int32_t i_block,j_block,k_block;
	    const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLIDs[b]];
	    block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	    switch (getSlicedCoordinate(slice)) {
	     case 0:
	       for (int32_t i=0; i<block::WIDTH_X; ++i) {
		  if (sim->x_crds_node[i_block*block::WIDTH_X+i] > getSliceOrigin(slice)) continue;
		  if (getSliceOrigin(slice) > sim->x_crds_node[i_block*block::WIDTH_X+i+1]) continue;
		  sliceIndex = i;
		  break;
	       }
	       break;
	     case 1:
	       for (int32_t j=0; j<block::WIDTH_Y; ++j) {
		  if (sim->y_crds_node[j_block*block::WIDTH_Y+j] > getSliceOrigin(slice)) continue;
		  if (getSliceOrigin(slice) > sim->y_crds_node[j_block*block::WIDTH_Y+j+1]) continue;
		  sliceIndex = j;
		  break;
	       }
	       break;
	     case 2:
	       for (int32_t k=0; k<block::WIDTH_Z; ++k) {
		  if (sim->z_crds_node[k_block*block::WIDTH_Z+k] > getSliceOrigin(slice)) continue;
		  if (getSliceOrigin(slice) > sim->z_crds_node[k_block*block::WIDTH_Z+k+1]) continue;
		  sliceIndex = k;
		  break;
	       }
	       break;
	     default:
	       continue;
	       break;
	    }
	    break;
	 }

	 if (blockLIDs.size() > 0 && sliceIndex == numeric_limits<int32_t>::max()) {
	    simClasses->logger << "(SEP OP WLM SCALAR) ERROR: Could not find sliceIndex" << endl << write;
	    exit(1);
	 }
	 
	 size_t arraySize;
	 switch (getSlicedCoordinate(slice)) {
	  case 0:
	    arraySize = block::WIDTH_Y*block::WIDTH_Z*simControl.N_wavelengthMeshCells;
	    break;
	  case 1:
	    arraySize = block::WIDTH_X*block::WIDTH_Z*simControl.N_wavelengthMeshCells;
	    break;
	  case 2:
	    arraySize = block::WIDTH_X*block::WIDTH_Y*simControl.N_wavelengthMeshCells;
	    break;
	  default:
	    continue;
	    break;
	 }
	 arraySize *= blockLIDs.size();

	 // Create a buffer for WLM scalar:
	 vector<Real> buffer;
	 buffer.resize(arraySize);

	 // Write all scalars:
	 for (size_t v=0; v<variableDataIDs.size(); ++v) {
	    // Get pointer to user data:
	    Real* array = simClasses->pargrid.getUserDataStatic<Real>(variableDataIDs[v]);
	    const uint32_t NWL = simControl.N_wavelengthMeshCells;

	    counter = 0;
	    for (size_t b=0; b<blockLIDs.size(); ++b) {
	       const pargrid::CellID blockLID = blockLIDs[b];
	       const Real* const ptr = array + blockLID*block::SIZE*NWL;

	       switch (getSlicedCoordinate(slice)) {
		case 0:
		  for (uint32_t l=0; l<NWL; ++l) for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) {
		     buffer[counter] = ptr[block::index(sliceIndex,j,k)*simControl.N_wavelengthMeshCells+l];
		     ++counter;
		    }
		  break;
		case 1:
		  for (uint32_t l=0; l<NWL; ++l) for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) {
		     buffer[counter] = ptr[block::index(i,sliceIndex,k)*simControl.N_wavelengthMeshCells+l];
		     ++counter;
		    }
		  break;
		case 2:
		  for (uint32_t l=0; l<NWL; ++l) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
		     buffer[counter] = ptr[block::index(i,j,sliceIndex)*simControl.N_wavelengthMeshCells+l];
		     ++counter;
		  }
		  break;
		default:
		  exit(1);
		  break;
	       }	       
	    }
	    
	    // Write variable to file:
	    map<string,string> xmlAttributes;
	    xmlAttributes["name"] = meshName + '/' + variableOutputNames[v];
	    xmlAttributes["mesh"] = meshName;
	    xmlAttributes["centering"] = "zone";
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

	    if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,arraySize,1,&(buffer[0])) == false) success = false;
	 }
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
} // namespace sep
