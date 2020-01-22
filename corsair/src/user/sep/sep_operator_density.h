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

#ifndef SEP_OPERATOR_DENSITY_H
#define SEP_OPERATOR_DENSITY_H

#include <stdint.h>
#include <map>
#include <vector>

#include <sep_operator_accumulation_base.h>

namespace sep {
   
   class OperatorDensity: public OperatorAccumulationBase {
    public:
      OperatorDensity(int32_t order=1);
      ~OperatorDensity();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   
    private:

      void accumulateBlock(pargrid::CellID blockLID,Real* accumArray,const std::vector<ParticleListBase*>& particleLists);
      uint32_t getNumberOfArrays() const;
      std::string getOutputName(uint32_t arrayIndex) const;
      std::string getOutputUnits(uint32_t arrayIndex) const;
      bool postProcessData(Real* array);
      bool setAccumulatedArray(uint32_t arrayIndex);
   };

} // namespace sep
   
#endif
