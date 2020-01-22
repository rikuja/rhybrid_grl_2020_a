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

#ifndef SEP_DISTRIB_WAVE_ENERGY_POWERLAW_H
#define SEP_DISTRIB_WAVE_ENERGY_POWERLAW_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_distrib_wave_energy_base_class.h"

namespace sep {

   /** Derived class WaveEnergySpectrumPowerlaw defines the energy spectrum
    * of parallel or antiparallel propagating Alfven waves. Wave spectral 
    * intensity is assumed to be a simple powerlaw in wavelength interval 
    * [0, maxLambda] (interval [k_min,infinity] in wave number space),
    * I(lambda) = I_0 * (lambda/lambda_0)^spectralIndex.
    * 
    * Instead of assuming that lambda_0 and maxLambda have the same values 
    * everywhere, these parameters are instead defined by the corresponding 
    * resonant particle energies (or speeds), which in practice means that
    * lambda0, maxLambda are inversely proportional to magnitude of ambient 
    * magnetic field. Constant I_0 is defined using total energy in Alfven 
    * waves (of given helicity) with respect to energy in ambient magnetic field.
    * 
    * This class is used to calculate energy carried by Alfven wave packets 
    * when they are injected to simulation at inflow boundaries.*/
   class WaveEnergySpectrumPowerlaw: public WaveEnergySpectrumBaseClass {
    public:
      WaveEnergySpectrumPowerlaw();
      ~WaveEnergySpectrumPowerlaw();
      
      bool finalize();
      Real getEnergyDensity(Real B0_mag,Real minLambda,Real maxLambda);
      Real getIntensity(Real B0_mag,Real wavelength);
      Real getIntensityDerivated(Real B0_mag,Real wavelength);
      Real getMaximumWavelength(Real B0_mag,int helicity) const;
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName);
      
    private:
      Real I0_lambda_L;        /**< Wave intensity constant for L-helicity waves, roughly corresponds 
				* to I0 parameter in wave number power law spectrum.*/
      Real I0_lambda_R;        /**< Wave intensity constant for R-helicity waves, roughly corresponds 
				* to I0 parameter in wave number power law spectrum.*/
      Real refResonantSpeed_L; /**< Reference resonant speed with L-helicity waves, corresponds to k0 
				parameter in spectral intensity defined using wave number.*/
      Real refResonantSpeed_R; /**< Reference resonant speed with R-helicity waves, corresponds to k0
				parameter in spectral intensity defined using wave number.*/
      Real maxResonantSpeed_L; /**< Maximum resonant speed with L-helicity waves, corresponds to k_min
				* parameter in spectral intensity defined using wave number.*/
      Real maxResonantSpeed_R; /**< Maximum resonant speed with R-helicity waves, corresponds to k_min
				* parameter in spectral intensity defined using wave number.*/
      Real relativeEnergy_L;
      Real relativeEnergy_R;
      Real spectralIndex_L;    /**< Power law spectral index of L-helicity distribution, in wave number
				* spectrum the corresponding spectral index is two larger.*/
      Real spectralIndex_R;    /**< Power law spectral index of R-helicity distribution, in wave number
				* spectrum the corresponding spectral index is two larger.*/
   };
   
   WaveEnergySpectrumBaseClass* waveEnergyPowerlawMaker();
   
} // namespace sep

#endif
