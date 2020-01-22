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

#ifndef SEP_OPERATOR_VOLUME_H
#define SEP_OPERATOR_VOLUME_H

#include <map>
#include <vector>

#include <dataoperator.h>

#include "sep_fields_container.h"

namespace sep {
   
   class OperatorVolume: public DataOperator {
    public:
      OperatorVolume();
      ~OperatorVolume();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);
   
    private:
      
      #if PROFILE_LEVEL > 0
         int totalTime;
      #endif
   };

} // namespace sep
   
#endif
