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

#ifndef SEP_OPERATOR_SPHERICAL_SURFACE_MESH_H
#define SEP_OPERATOR_SPHERICAL_SURFACE_MESH_H

#include <vector>
#include <unordered_map>

#include <dataoperator.h>

namespace sep {
   class OperatorSphericalSurfaceMesh: public DataOperator {
    public:
      OperatorSphericalSurfaceMesh();
      ~OperatorSphericalSurfaceMesh();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

    private:

      Real radius0;

      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif
   };
   
} // namespace sep

#endif
