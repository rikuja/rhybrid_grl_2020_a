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

#ifndef SEP_OPERATOR_SPATIAL_LINEOUT_H
#define SEP_OPERATOR_SPATIAL_LINEOUT_H

#include <vector>

#include <dataoperator.h>

namespace sep {

   class OperatorSpatialLineout: public DataOperator {
    public:
      OperatorSpatialLineout();
      virtual ~OperatorSpatialLineout();
      
      virtual bool finalize();
      virtual bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      
    protected:
      bool baseClassInitialized;
      
      bool addConfigFileItems(ConfigReader& cr);
      bool getAcceptedBlocks(size_t line,std::vector<pargrid::CellID>& inner,std::vector<pargrid::CellID>& boundary);
      size_t getNumberOfLineouts() const;
      uint8_t getLineoutCoordinate(size_t line) const;
      std::string getLineoutName(size_t line) const;
      bool getRemoteBlocks(std::vector<pargrid::CellID>& remote);

    private:
      std::vector<uint8_t> lineoutCoordinates;
      std::vector<Real> lineoutOriginsX;
      std::vector<Real> lineoutOriginsY;
      std::vector<Real> lineoutOriginsZ;
      std::vector<uint8_t> vectorSizes;
   };

} // namespace sep
   
#endif
