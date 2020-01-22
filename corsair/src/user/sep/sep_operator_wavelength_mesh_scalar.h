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

#ifndef SEP_OPERATOR_WAVELENGTH_MESH_SCALAR_H
#define SEP_OPERATOR_WAVELENGTH_MESH_SCALAR_H

#include <dataoperator.h>

#include "sep_operator_spatial_slice.h"

namespace sep {
   
   class WavelengthMeshScalarOP: public SpatialSliceOP {
    public:
      WavelengthMeshScalarOP();
      ~WavelengthMeshScalarOP();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);
      
    private:
      
      bool addConfigFileItems(ConfigReader& cr);      
      bool searchVariableDataIDs();
      
      std::vector<pargrid::DataID> variableDataIDs;
      std::vector<std::string> variableNames;
      std::vector<std::string> variableOutputNames;
      bool variablesSearched;
      
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif
   };
   
} // namespace sep

#endif
