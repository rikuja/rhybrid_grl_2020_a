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

#include <cstdlib>
#include <iostream>
#include <limits>

#include "sep_base_class_particle_splitter.h"
#include "sep_simcontrol.h"

using namespace std;

namespace sep {
   
   extern sep::SimControl simControl;

   ParticleSplitterBase::ParticleSplitterBase() {
      finalize();
   }
   
   ParticleSplitterBase::~ParticleSplitterBase() {
      finalize();
   }
   
   bool ParticleSplitterBase::finalize() {
      initialized = false;
      return true;
   }

   bool ParticleSplitterBase::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      initialized = true;
      return initialized;
   }
   
   bool getMinimumWaveEnergies(const pargrid::CellID& blockLID,int waveMode,std::vector<Real>& minWaveEnergyDens) {
      corsair::ObjectWrapper& owrapper = corsair::getObjectWrapper();
      
      // Get wave energy array (may be NULL);
      const Real* RESTRICT waveEnergy = NULL;
      if (waveMode < 0) {
	 waveEnergy = owrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      } else {
	 waveEnergy = owrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      }
      if (waveEnergy == NULL) return false;
      
      const int32_t NWL = simControl.N_wavelengthMeshCells;
      const int32_t SIZE = block::SIZE*simControl.N_wavelengthMeshCells;
      minWaveEnergyDens.resize(NWL/2);

      // Invalidate minWaveEnergyDens:
      for (size_t i=0; i<minWaveEnergyDens.size(); ++i) minWaveEnergyDens[i] = std::numeric_limits<Real>::max();

      // Calculate minimum value of wave energy / wavelength in each wavelength bin (par waves):
      int32_t l_min = NWL/2-1;
      int32_t l_max = NWL/2;
      do {
	 for (int32_t cell=0; cell<block::SIZE; ++cell) {
	    const Real leftValue = waveEnergy[blockLID*SIZE+cell*NWL + l_min];
	    const Real rghtValue = waveEnergy[blockLID*SIZE+cell*NWL + l_max];
	    const Real smallerValue = std::min(leftValue,rghtValue) / simControl.wavelengthMeshCellSizes[l_min];
	    
	    if (smallerValue < minWaveEnergyDens[l_max-NWL/2]) minWaveEnergyDens[l_max-NWL/2] = smallerValue;
	 }

	 --l_min;
	 ++l_max;
	 if (l_min < 0) break;
      } while (true);
      
      for (size_t i=0; i<minWaveEnergyDens.size()-1; ++i) {
	 if (minWaveEnergyDens[i] < minWaveEnergyDens[i+1]) minWaveEnergyDens[i+1] = minWaveEnergyDens[i];
      }
      
      return true;
   }

} // namespace sep
