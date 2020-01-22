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

#ifndef SEP_BASE_CLASS_WAVELENGTH_SCALING_H
#define SEP_BASE_CLASS_WAVELENGTH_SCALING_H

#include <vector>

#include <main.h>

namespace sep {

   class WavelengthScalingBaseClass {
    public:
      WavelengthScalingBaseClass();
      virtual ~WavelengthScalingBaseClass();

      virtual bool finalize();
      virtual Real getScaleFactor(const Real& t,const Real* pos) =0;
      virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,
			      ConfigReader& cr,int propagationDirection);
      bool isInitialized() const;

    protected:
      bool initialized;
   };

} // namespace sep

#endif
