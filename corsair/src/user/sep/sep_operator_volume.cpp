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
#include "sep_operator_volume.h"

using namespace std;

namespace sep {
   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   OperatorVolume::OperatorVolume(): DataOperator() {
      initialized = false;
      
      #if PROFILE_LEVEL > 0
         totalTime = -1;
      #endif
   }
   
   OperatorVolume::~OperatorVolume() {
      finalize();
   }
   
   bool OperatorVolume::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   std::string OperatorVolume::getName() const {return "sepOperatorVolume";}
   
   bool OperatorVolume::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP VOLUME) ERROR: DataOperator base class failed to initialize!" << endl << write;
	 initialized = false;
      }
   
      return initialized;
   }
   
   bool OperatorVolume::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;
      /*
      if (simControl.cellVolumes == NULL) {
	 simClasses->logger << "(SEP OP VOLUME) ERROR: Cell volumes array is NULL" << endl << write;
	 return false;
      }*/
      recalculateCellVolumes(*sim,*simClasses);
      
      #if PROFILE_LEVEL > 0
         profile::start("SpatialVolume",totalTime);
      #endif
      
      // Write cell volume array to output file:
      std::map<std::string,std::string> attribs;
      attribs["name"] = spatMeshName + "/CellVolume";
      attribs["mesh"] = spatMeshName;
      attribs["unit"] = "m^3";
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
      
      const uint64_t arraySize = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
      const uint64_t vectorSize = 1;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,simControl.cellVolumes) == false) success = false;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
      
   
} // namespace sep
