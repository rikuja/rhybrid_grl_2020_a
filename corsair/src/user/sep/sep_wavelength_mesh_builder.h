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

#ifndef SEP_WAVELENGTH_MESH_BUILDER_H
#define SEP_WAVELENGTH_MESH_BUILDER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace sep {

   typedef Real (*getLogicalWavelengthFunction)(const Real& wavelength);
   typedef Real (*getPhysicalWavelengthFunction)(const Real& l);
   typedef Real (*getWavelengthScaleFactorFunction)(const Real& wavelenght);
   
   namespace wmesh {
      typedef int16_t validDatatype;
      
      namespace linear {
	 Real getLogicalWavelength(const Real& wavelength);
	 Real getPhysicalWavelength(const Real& l);
	 Real getPhysicalWavelengthCellSize(const Real& l);
	 Real getWavelengthScaleFactor(const Real& wavelength);
      }
      
      namespace logarithmic {
	 Real getLogicalWavelength(const Real& wavelength);
	 Real getPhysicalWavelength(const Real& l);
	 Real getPhysicalWavelengthCellSize(const Real& l);
	 Real getWavelengthScaleFactor(const Real& wavelength);
      }
   }
   
   uint32_t waveIndex(uint32_t i,uint32_t j,uint32_t k,uint32_t l);
   
   bool buildWaveMesh(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   
} // namespace sep
   
#endif
