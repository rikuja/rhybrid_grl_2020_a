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

#ifndef SEP_OPERATOR_SPATIAL_SLICE_H
#define SEP_OPERATOR_SPATIAL_SLICE_H

#include <vector>

#include <dataoperator.h>

namespace sep {

   class SpatialSliceOP: public DataOperator {
    public:
      SpatialSliceOP();
      virtual ~SpatialSliceOP();
      
      virtual bool finalize();
      virtual bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      
    protected:
      bool baseClassInitialized;
      
      bool addConfigFileItems(ConfigReader& cr);
      uint32_t calculateNewGlobalID(uint32_t* indices) const;
      void getAcceptedBlocks(uint32_t slice,std::vector<pargrid::CellID>& newBlockGIDs,std::vector<pargrid::CellID>& blockLIDs);
      size_t getNumberOfSlices() const;
      uint8_t getSlicedCoordinate(uint32_t slice) const;
      uint32_t getSliceIndex(uint32_t slice) const;
      Real getSliceOrigin(uint32_t slice) const;
      bool isCylindrical(uint32_t slice) const;
      void prepareSlice(uint32_t slice,uint32_t N_cells);

    private:
      uint32_t x_blocks;
      uint32_t y_blocks;
      uint32_t z_blocks;
      std::vector<uint8_t> sliceCoordinates;
      std::vector<Real> sliceOrigins;
      std::vector<bool> sliceGeometries;
      std::vector<uint32_t> sliceIndices;
      uint32_t N_cells;                      /**< Number of cells in sliced dimension mesh.*/
   };

} // namespace sep
   
#endif
