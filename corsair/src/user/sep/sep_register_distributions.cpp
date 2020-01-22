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

#include "sep_object_wrapper.h"
#include "sep_register_distributions.h"

// Energy distributions
#include "sep_distrib_energy_kappa.h"
#include "sep_distrib_energy_mono.h"
#include "sep_distrib_energy_powerlaw.h"

// Pitch distributions
#include "sep_distrib_pitch_isotropic.h"
#include "sep_distrib_pitch_mono.h"

// Wave energy spectra:
#include "sep_distrib_wave_energy_base_class.h"
#include "sep_distrib_wave_energy_powerlaw.h"

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * particle and wave packet pitch, energy, etc. distribution functions.
    * @return If true, distributions were registered successfully.*/
   bool registerDistributions() {
      bool success = true;

      // Register energy distributions:
      if (sep::getObjectWrapper().energyDistribContainer.registerDistribution("kappa",sep::energyDistribKappaFinalize,
									      sep::energyDistribKappaGetDistrib,
									      sep::energyDistribKappaGetEnergy,
									      sep::energyDistribKappaInitialize) == false) success = false;
      if (sep::getObjectWrapper().energyDistribContainer.registerDistribution("mono",sep::energyDistribMonoFinalize,
									      sep::energyDistribMonoGetDistrib,
									      sep::energyDistribMonoGetEnergy,
									      sep::energyDistribMonoInitialize) == false) success = false;
      if (sep::getObjectWrapper().energyDistribContainer.registerDistribution("powerlaw",sep::energyDistribPowerlawFinalize,
									      sep::energyDistribPowerlawGetDistrib,
									      sep::energyDistribPowerlawGetEnergy,
									      sep::energyDistribPowerlawInitialize) == false) success = false;

      // Register pitch distributions:
      if (sep::getObjectWrapper().pitchDistribContainer.registerDistribution("isotropic",sep::pitchDistribIsotropicFinalize,
									     sep::pitchDistribIsotropicGetPitch,
									     sep::pitchDistribIsotropicInitialize) == false) success = false;
      if (sep::getObjectWrapper().pitchDistribContainer.registerDistribution("mono",sep::pitchDistribMonoFinalize,
									     sep::pitchDistribMonoGetPitch,
									     sep::pitchDistribMonoInitialize) == false) success = false;
      
      // Register wave energy distributions:
      if (sep::getObjectWrapper().waveEnergySpectrumFactory.registerMaker("PowerLaw",sep::waveEnergyPowerlawMaker) == false) success = false;
      
      return success;
   }
   
} // namespace sep
