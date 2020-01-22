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

#include <cstdlib>
#include <iostream>

#include <pargrid_definitions.h>
#include <simulation.h>
#include "sep_simcontrol.h"

using namespace std;

sep::SimControl::SimControl(): EPSILON(numeric_limits<Real>::min()*100), BLOCK_EPSILON(1.0e-9) {
   alfvenSign = 0;
   cellVolumes = NULL;
   fieldsFinalize = NULL;
   fieldsGetFields = NULL;
   fieldsGetPlasmaState = NULL;
   fieldsGetState = NULL;
   fieldsInitialize = NULL;
   shock = NULL;
   getLogicalWavelength = NULL;
   getPhysicalWavelength = NULL;
   getWavelengthScaleFactor = NULL;
   particleSplitter = NULL;
   N_shockSurfaceRefinements = 3;

   parWaveEnergy = NULL;
   antiparWaveEnergy = NULL;

   propagateAlfvenWaves = true;
   scatterParticles     = true;
   applyWaveGrowth      = true;
   includeAntiparWaves  = true;
   includeParWaves      = true;
   includeShock         = false;
   useShockUpwinding    = false;
   useUpwinding         = false;

   scatterAntiparallel = scatterParticles;
   scatterParallel     = scatterParticles;
   applyAntiparallelWaveGrowth = applyWaveGrowth;
   applyParallelWaveGrowth     = applyWaveGrowth;

   courantNumber = NAN;
   N_lagrangianSpecies = 0;
   N_particleSpecies   = 0;
   timestepRecalculateInterval = 0;
   currentMaxScatteringSubsteps = 1;

   for (int i=0; i<3; ++i) V_frame[i] = 0.0;
   timestepCellVolumesCalculated = -1;
   
   // Invalidate wavelength mesh variables:
   N_wavelengthMeshCells = numeric_limits<unsigned int>::max();
   wavelengthMeshCellSizes = NULL;
   wavelengthMeshNodeCoordinates = NULL;
   
   // Invalidate pargrid DataIDs:
   waveGrowthDataID               = pargrid::getInvalidID<pargrid::DataID>();
   particleWeightDataID           = pargrid::getInvalidID<pargrid::DataID>();
   antiparAlfvenWaveEnergyDataID  = pargrid::getInvalidID<pargrid::DataID>();
   antiparWaveGrowthDataID        = pargrid::getInvalidID<pargrid::DataID>();
   parAlfvenWaveEnergyDataID      = pargrid::getInvalidID<pargrid::DataID>();
   parWaveGrowthDataID            = pargrid::getInvalidID<pargrid::DataID>();
   temporaryWaveEnergyDataID      = pargrid::getInvalidID<pargrid::DataID>();
   antiparAlfvenDataID            = pargrid::getInvalidID<pargrid::DataID>();
   parAlfvenDataID                = pargrid::getInvalidID<pargrid::DataID>();
   maxAntiparValidWavelengthBinsDataID = pargrid::getInvalidID<pargrid::DataID>();
   maxParValidWavelengthBinsDataID = pargrid::getInvalidID<pargrid::DataID>();
}

sep::SimControl::~SimControl() {
   delete [] cellVolumes; cellVolumes = NULL;
   delete [] wavelengthMeshCellSizes; wavelengthMeshCellSizes = NULL;
   delete [] wavelengthMeshNodeCoordinates; wavelengthMeshNodeCoordinates = NULL;
   delete shock; shock = NULL;
   delete particleSplitter; particleSplitter = NULL;
}

