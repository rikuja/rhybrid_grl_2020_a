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

#ifndef SEP_OPERATOR_DMUMU_H
#define SEP_OPERATOR_DMUMU_H

#include <vector>

#include <dataoperator.h>

#include "sep_operator_spatial_slice.h"
#include "sep_operator_differential_flux.h"

namespace sep {
   
   class OperatorDmumu: public SpatialSliceOP {
    public:
      OperatorDmumu();
      ~OperatorDmumu();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

    private:
      bool writeDmumu(int32_t slice,size_t currentInstrument,size_t currentChannel,
		      const std::vector<pargrid::CellID>& validLocalIDs,
		      const std::vector<pargrid::CellID>& newGlobalIDs);
      
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif
      
      size_t currentInstrument;
      size_t currentChannel;
      std::string meshName;
      uint32_t N_pitchCells;
      uint32_t totalNumberOfArrays;
      std::vector<spacecraft::FluxInstrument> instruments;
      std::vector<std::pair<size_t,size_t> > instrumentIndices;
   };
   
} // namespace sep

#endif
