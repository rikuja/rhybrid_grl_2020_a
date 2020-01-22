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

#ifndef SEP_BASE_CLASS_PARTICLE_SPLITTER_H
#define SEP_BASE_CLASS_PARTICLE_SPLITTER_H

#include <vector>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <particle_list_base.h>

namespace sep {

   class ParticleSplitterBase {
    public:
      ParticleSplitterBase();
      virtual ~ParticleSplitterBase();

      virtual bool finalize();
      virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      virtual bool split() = 0;

    protected:
      bool initialized;
   };

   bool getMinimumWaveEnergies(const pargrid::CellID& blockLID,int waveMode,std::vector<Real>& minWaveEnergyDens);

} // namespace sep
   
#endif
