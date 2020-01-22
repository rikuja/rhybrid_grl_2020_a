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

#ifndef SEP_OPERATOR_MFP_H
#define SEP_OPERATOR_MFP_H

#include <stdint.h>
#include <map>
#include <vector>

#include "sep_operator_differential_flux.h"

namespace sep {
   
   class OperatorMeanFreePath: public OperatorDifferentialFlux {
    public:
      OperatorMeanFreePath(int32_t order=1);
      ~OperatorMeanFreePath();
      
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   
    private:
      void accumulateBlock(pargrid::CellID blockLID,Real* accumArray,const std::vector<ParticleListBase*>& particleLists);
      std::string getOutputName(uint32_t arrayIndex) const;
      std::string getOutputUnits(uint32_t arrayIndex) const;
      bool postProcessData(Real* array);
      bool preProcessData(Real* array);
      
      Real lengthScale;
      std::string lengthScaleString;
   };

} // namespace sep
   
#endif
