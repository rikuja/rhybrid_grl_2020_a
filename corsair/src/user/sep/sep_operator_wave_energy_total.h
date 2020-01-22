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

#ifndef SEP_OPERATOR_WAVE_ENERGY_TOTAL_H
#define SEP_OPERATOR_WAVE_ENERGY_TOTAL_H

#include <map>
#include <vector>

#include <dataoperator.h>

namespace sep {

   /** This DataOperator calculates total energy density in antiparallel and 
    * parallel propagating Alfven waves (sum of energies L and R helicity waves).*/
   class OperatorWaveEnergyTotal: public DataOperator {
    public:
      OperatorWaveEnergyTotal();
      ~OperatorWaveEnergyTotal();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);
   
    private:
      
      bool writeArray(const std::string& spatMeshName,const std::string& name,const Real* buffer);
      bool writeWaveGrowth(const std::string& spatMeshName,const std::string& name,pargrid::DataID dataID);
      bool writeWaveMode(const std::string& spatMeshName,const std::string& name,pargrid::DataID dataID);
      
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif
   };

} // namespace sep
   
#endif
