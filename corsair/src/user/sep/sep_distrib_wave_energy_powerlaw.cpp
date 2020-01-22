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
#include <cmath>

#include "sep_simcontrol.h"
#include "sep_distrib_wave_energy_powerlaw.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const Real DEF_REAL = numeric_limits<Real>::infinity();

   /** Maker function for class WaveEnergySpectrumPowerLaw. This function 
    * is registered to an object factory that creates objects of type 
    * WaveEnergySpectrumBaseClass.
    * @return New instance of class WaveEnergySpectrumPowerLaw.
    * @see ObjectFactoryGeneric.*/
   WaveEnergySpectrumBaseClass* waveEnergyPowerlawMaker() {return new WaveEnergySpectrumPowerlaw;}

   /** Default constructor, initializes member variables to NAN values.*/
   WaveEnergySpectrumPowerlaw::WaveEnergySpectrumPowerlaw() {
      I0_lambda_L        = NAN;
      I0_lambda_R        = NAN;
      maxResonantSpeed_L = NAN;
      maxResonantSpeed_R = NAN;
      spectralIndex_L    = NAN;
      spectralIndex_R    = NAN;
   }
   
   /** Destructor, calls WaveEnergySpectrumPowerLaw::finalize() member function.*/
   WaveEnergySpectrumPowerlaw::~WaveEnergySpectrumPowerlaw() { 
      finalize();
   }

   /** Finalizer function, does not do anything.
    * @return If true, class finalized successfully.*/
   bool WaveEnergySpectrumPowerlaw::finalize() {return true;}

   /** Get wave energy density in wavelength interval [minLambda,maxLambda[. Positive 
    * (negative) wavelengths minLambda, maxLambda indicate L (R) helicity.
    * If NDEBUG preprocesses macro is not set, function checks if class initialized successfully 
    * and that minLambda and maxLambda have the same sign.
    * @param B0_mag Magnitude of ambient magnetic field.
    * @param minLambda Minimum wavelength.
    * @param maxLambda Maximum wavelength.
    * @return Energy density of waves with wavelengths in requested interval.*/
   Real WaveEnergySpectrumPowerlaw::getEnergyDensity(Real B0_mag,Real minLambda,Real maxLambda) {
      #ifndef NDEBUG
      if (initialized == false) {
	 simClasses->logger << "(SEP WAVE POWERLAW) ERROR: Wave energy function failed to init" << endl << write;
	 exit(1);
      }
      if (minLambda*maxLambda < 0.0) {
	 simClasses->logger << "(SEP WAVE POWERLAW) ERROR: minLambda and maxLambda must have same sign" << endl;
	 simClasses->logger << "\t Got values min: " << minLambda << "\t max: " << maxLambda << endl << write;
	 exit(1);
      }
      #endif
      
      // Note that both minLambda and maxLambda can be greater than maxWavelength,
      // in which case their values need to be set to maxWavelength.
      const Real omega = constants::CHARGE_ELEMENTARY*B0_mag/constants::MASS_PROTON;
      const Real EPSILON = numeric_limits<Real>::min();
      if (maxLambda > 0.0) {
	 // Waves with wavelength > 0 have L helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_L;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_L;
	 minLambda = max(+1.0,min(maxWavelength,minLambda));
	 maxLambda = max(+1.0,min(maxWavelength,maxLambda));

	 // Flatten spectral energy density at wavelengths shorten than LAMBDA0:
	 /*if (minLambda < simControl.lambda0 && maxLambda < simControl.lambda0) {
	    return 0.0;
	 } else if (minLambda < simControl.lambda0) {
	 if (minLambda < simControl.lambda0) {
	    const Real I_min_L = pow(maxWavelength/simControl.lambda0+EPSILON,spectralIndex_L-1);
	    const Real I_max_L = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_L-1);
	    return 0.5*0.5*relativeEnergy_L*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_L-I_min_L);
	 } else {
	    const Real I_min_L = pow(maxWavelength/minLambda+EPSILON,spectralIndex_L-1);
	    const Real I_max_L = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_L-1);
	    return 0.5*relativeEnergy_L*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_L-I_min_L);
	 }*/
	 
	 Real I_min_L,I_max_L;
	 Real rvalue = 0.0;
	 /*
	 if (minLambda < simControl.lambda0 && maxLambda < simControl.lambda0) {
	    I_min_L = getIntensity(B0_mag,simControl.lambda0);
	    rvalue += I_min_L * (maxLambda-minLambda)/constants::PERMEABILITY;
	    return rvalue;
	 } else if (minLambda < simControl.lambda0) {
	    I_min_L = getIntensity(B0_mag,simControl.lambda0);
	    rvalue += I_min_L * (simControl.lambda0 - minLambda)/constants::PERMEABILITY;
	    minLambda = simControl.lambda0;
	 }*/
	 minLambda = max(1.0,minLambda);
	 I_min_L = pow(maxWavelength/minLambda+EPSILON,spectralIndex_L-1);
	 I_max_L = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_L-1);
	 rvalue += 0.5*relativeEnergy_L*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_L-I_min_L);
	 return rvalue;
      } else {
	 // Waves with wavelength < 0 have R helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_R;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_R;
	 const Real tmp = fabs(minLambda);
	 minLambda = max(+1.0,min(maxWavelength,fabs(maxLambda)));
	 maxLambda = max(+1.0,min(maxWavelength,tmp));

	 // Flatten spectral energy density at wavelengths shorten than LAMBDA0:
	 /*if (minLambda < simControl.lambda0 && maxLambda < simControl.lambda0) {
	    return 0.0;
	 } else if (minLambda < simControl.lambda0) {
	    const Real I_min_R = pow(maxWavelength/simControl.lambda0+EPSILON,spectralIndex_R-1);
	    const Real I_max_R = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_R-1);
	    return 0.5*0.5*relativeEnergy_R*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_R-I_min_R);
	 } else {
	    const Real I_min_L = pow(maxWavelength/minLambda+EPSILON,spectralIndex_L-1);
	    const Real I_max_L = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_L-1);
	    return 0.5*relativeEnergy_L*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_L-I_min_L);
	 }*/

	 Real I_min_L,I_max_L;
	 Real rvalue = 0.0;
	 /*if (minLambda < simControl.lambda0 && maxLambda < simControl.lambda0) {
	    I_min_L = getIntensity(B0_mag,simControl.lambda0);
	    //rvalue += I_min_L * (simControl.lambda0 - minLambda)/constants::PERMEABILITY;
	    rvalue += I_min_L * (maxLambda-minLambda)/constants::PERMEABILITY;
	    return rvalue;
	 } else if (minLambda < simControl.lambda0) {
	    I_min_L = getIntensity(B0_mag,simControl.lambda0);
	    rvalue += I_min_L * (simControl.lambda0 - minLambda)/constants::PERMEABILITY;
	    minLambda = simControl.lambda0;
	 }*/
	 minLambda = max(1.0,minLambda);
	 I_min_L = pow(maxWavelength/minLambda+EPSILON,spectralIndex_R-1);
	 I_max_L = pow(maxWavelength/maxLambda+EPSILON,spectralIndex_R-1);
	 rvalue += 0.5*relativeEnergy_R*B0_mag*B0_mag/constants::PERMEABILITY*(I_max_L-I_min_L);
	 return rvalue;
      }
      
      cerr << "ERROR" << endl;
      exit(1);
   }
   
   
   Real WaveEnergySpectrumPowerlaw::getIntensity(Real B0_mag,Real wavelength) {
      #ifndef NDEBUG
      if (initialized == false) {
	 simClasses->logger << "(SEP WAVE POWERLAW) ERROR: class is not initialized" << endl << write;
	 exit(1);
      }
      #endif
      
      // Note that wavelength can be larger than maxWavelength, in which case zero value is returned:
      const Real omega = constants::CHARGE_ELEMENTARY*B0_mag/constants::MASS_PROTON;
      const Real EPSILON = numeric_limits<Real>::min();
      if (wavelength > 0.0) {
	 // Waves with wavelength > 0 have L helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_L;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_L;
	 if (wavelength > maxWavelength) return 0.0;
	 const Real I0_L = I0_lambda_L * B0_mag*B0_mag / maxWavelength;
	 return I0_L * pow(maxWavelength/wavelength+EPSILON,spectralIndex_L);
      } else {
	 // Waves with wavelength < 0 have R helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_R;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_R;
	 wavelength = -wavelength;
	 if (wavelength > maxWavelength) return 0.0;
	 const Real I0_R = I0_lambda_R * B0_mag*B0_mag / maxWavelength;
	 return I0_R * pow(maxWavelength/wavelength+EPSILON,spectralIndex_R);
      }
   }
   
   Real WaveEnergySpectrumPowerlaw::getIntensityDerivated(Real B0_mag,Real wavelength) { 
      #ifndef NDEBUG
      if (initialized == false) {
	 simClasses->logger << "(SEP WAVE POWERLAW) ERROR: class is not initialized" << endl << write;
	 exit(1);
      }
      #endif
      
      // Note that wavelength can be larger than maxWavelength, in which case zero value is returned:
      const Real omega = constants::CHARGE_ELEMENTARY*B0_mag/constants::MASS_PROTON;
      const Real EPSILON = numeric_limits<Real>::min();
      if (wavelength > 0.0) {
	 // Waves with wavelength > 0 have L helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_L;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_L;
	 if (wavelength > maxWavelength) return 0.0;
	 const Real I0_L = I0_lambda_L * B0_mag*B0_mag / maxWavelength;
	 return -I0_L * pow(maxWavelength/wavelength+EPSILON,spectralIndex_L+1) * spectralIndex_L/maxWavelength;
      } else {
	 // Waves with wavelength < 0 have R helicity. Convert reference and maximum
	 // resonant energies to wavelengths, and calculate intensity spectrum constant:
	 //const Real refWavelength = 2*M_PI/omega*refResonantSpeed_R;
	 const Real maxWavelength = 2*M_PI/omega*maxResonantSpeed_R;
	 wavelength = -wavelength;
	 if (wavelength > maxWavelength) return 0.0;
	 const Real I0_R = I0_lambda_R * B0_mag*B0_mag / maxWavelength;
	 return I0_R * pow(maxWavelength/wavelength+EPSILON,spectralIndex_R+1) * spectralIndex_R/maxWavelength;
      }
   }
   
   Real WaveEnergySpectrumPowerlaw::getMaximumWavelength(Real B0_mag,int helicity) const {
      const Real omega = constants::CHARGE_ELEMENTARY*B0_mag/constants::MASS_PROTON;
      if (helicity > 0) {
	 // L helicity:
	 return +2*M_PI/omega*maxResonantSpeed_L;
      } else {
	 // R helicity
	 return -2*M_PI/omega*maxResonantSpeed_R;
      }
   }

   /** Class initializer, reads values from configuration file and 
    * initializes member variables to correct values.
    * @param sim Struct containing variables controlling execution of simulation.
    * @param simClasses Struct containing generic use simulation classes.
    * @param cr Configuration file reader.
    * @param regionName Name of configuration file region that contains parameters 
    * for this instance of the class.
    * @return If true, class initialized successfully and is ready for use.*/
   bool WaveEnergySpectrumPowerlaw::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName) {
      // Prevent multiple initializations:
      if (WaveEnergySpectrumBaseClass::initialized == true) return true;
      bool success = true;

      // Init base class:
      if (WaveEnergySpectrumBaseClass::initialize(sim,simClasses,cr,regionName) == false) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Failed to initialize base class" << endl << write;
	 success = false;
      }

      cr.add(regionName+".relative_energy_L","Relative energy in L helicity waves (float).",DEF_REAL);
      cr.add(regionName+".relative_energy_R","Relative energy in R helicity waves (float).",DEF_REAL);
      cr.add(regionName+".spectral_index_L","Spectral index of I(k) spectrum of L helicity waves, must be positive (float).",DEF_REAL);
      cr.add(regionName+".spectral_index_R","Spectral index of I(k) spectrum of R helicity waves, must be positive (float).",DEF_REAL);
      cr.add(regionName+".max_resonant_energy_L","Maximum resonant energy of particles with L helicity waves in keV/amu (float).",DEF_REAL);
      cr.add(regionName+".max_resonant_energy_R","Maximum resonant energy of particles with R helicity waves in keV/amu (float).",DEF_REAL);
      cr.add(regionName+".ref_resonant_energy_L","Reference resonant energy of particles with L helicity waves in keV/amu (float).",DEF_REAL);
      cr.add(regionName+".ref_resonant_energy_R","Reference resonant energy of particles with R helicity waves in keV/amu (float).",DEF_REAL);
      cr.parse();
      
      Real maxResonantEnergy_L,maxResonantEnergy_R;
      Real refResonantEnergy_L,refResonantEnergy_R;
      cr.get(regionName+".relative_energy_L",relativeEnergy_L);
      cr.get(regionName+".relative_energy_R",relativeEnergy_R);
      cr.get(regionName+".spectral_index_L",spectralIndex_L);
      cr.get(regionName+".spectral_index_R",spectralIndex_R);
      cr.get(regionName+".max_resonant_energy_L",maxResonantEnergy_L);
      cr.get(regionName+".max_resonant_energy_R",maxResonantEnergy_R);
      cr.get(regionName+".ref_resonant_energy_L",refResonantEnergy_L);
      cr.get(regionName+".ref_resonant_energy_R",refResonantEnergy_R);
      
      // Sanity check on input parameters:
      if (spectralIndex_L <= 0.0 || spectralIndex_L == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Spectral index given with '" << regionName+".spectral_index_L' must be positive" << endl << write;
	 success = false;
      }
      if (spectralIndex_R <= 0.0 || spectralIndex_R == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Spectral index given with '" << regionName+".spectral_index_R' must be positive" << endl << write;
	 success = false;
      }
      if (relativeEnergy_L < 0.0 || relativeEnergy_L == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Relative energy given with '" << regionName+".relative_energy_L' >= 0 required" << endl << write;
	 success = false;
      }
      if (relativeEnergy_R < 0.0 || relativeEnergy_R == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Relative energy given with '" << regionName+".relative_energy_R' >= 0 required" << endl << write;
	 success = false;
      }
      if (maxResonantEnergy_L <= 0.0 || maxResonantEnergy_L == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Max. resonant energy/amu given with '" << regionName+".max_resonant_energy_L'" << endl;
	 simClasses.logger << "                           must be positive" << endl << write;
	 success = false;
      }
      if (maxResonantEnergy_R <= 0.0 || maxResonantEnergy_R == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Max. resonant energy/amu given with '" << regionName+".max_resonant_energy_R'" << endl;
	 simClasses.logger << "                           must be positive" << endl << write;
	 success = false;
      }
      if (refResonantEnergy_L <= 0.0 || refResonantEnergy_L == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Reference resonant energy/amu, given with parameter" << endl;
	 simClasses.logger << "                           '" << regionName+".ref_resonant_energy_L' must be greater than zero" << endl;
	 simClasses.logger << "                           Value read is '" << refResonantEnergy_L << "'" << endl << write;
	 success = false;
      }
      if (refResonantEnergy_R <= 0.0 || refResonantEnergy_R == DEF_REAL) {
	 simClasses.logger << "(SEP WAVE POWERLAW) ERROR: Reference resonant energy/amu, given with parameter" << endl;
	 simClasses.logger << "                           '" << regionName+".ref_resonant_energy_R' must be greater than zero" << endl;
	 simClasses.logger << "                           Value read is '" << refResonantEnergy_R << "'" << endl << write;
	 success = false;
      }
      
      // Spectral index of lambda power law distribution is 2 lower than in k distribution:
      spectralIndex_L = -spectralIndex_L + 2;
      spectralIndex_R = -spectralIndex_R + 2;

      // Convert resonant energies into speeds. Resonant speeds will be converted into 
      // wavelengths when injecting particles:
      const Real units = 1000.0*constants::CHARGE_ELEMENTARY;
      maxResonantSpeed_L = sqrt(2.0*maxResonantEnergy_L*units/constants::MASS_PROTON);
      maxResonantSpeed_R = sqrt(2.0*maxResonantEnergy_R*units/constants::MASS_PROTON);
      refResonantSpeed_L = sqrt(2.0*refResonantEnergy_L*units/constants::MASS_PROTON);
      refResonantSpeed_R = sqrt(2.0*refResonantEnergy_R*units/constants::MASS_PROTON);

      // Calculate spectral intensity constants. These still need to be multiplied by 
      // (B0^2 / maxResonantWavelength) but that is factored in later:
      I0_lambda_L = 0.5*(1-spectralIndex_L)*relativeEnergy_L;// * pow(refResonantSpeed_L/maxResonantSpeed_L,spectralIndex_L);
      I0_lambda_R = 0.5*(1-spectralIndex_R)*relativeEnergy_R;// * pow(refResonantSpeed_R/maxResonantSpeed_R,spectralIndex_R);

      if (success == true) {
	 simClasses.logger << "(SEP WAVE POWERLAW) Power law parameters for '" << regionName << "' are:" << endl;
	 simClasses.logger << "\t Spectral Index L: " << spectralIndex_L << endl;
	 simClasses.logger << "\t Spectral Index R: " << spectralIndex_R << endl;
	 simClasses.logger << "\t I_0L/B0^2*speed : " << I0_lambda_L << endl;
	 simClasses.logger << "\t I_0R/B0^2*speed : " << I0_lambda_R << endl;
	 simClasses.logger << "\t lambda_0_L*B0   : " << 2*M_PI*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*refResonantSpeed_L << endl;
	 simClasses.logger << "\t lambda_0_R*B0   : " << 2*M_PI*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*refResonantSpeed_R << endl;
	 simClasses.logger << "\t lambda_max_L*B0 : " << 2*M_PI*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*maxResonantSpeed_L << endl;
	 simClasses.logger << "\t lambda_max_R*B0 : " << 2*M_PI*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*maxResonantSpeed_R << endl;
	 simClasses.logger << "\t max. res speed L: " << maxResonantSpeed_L << endl;
	 simClasses.logger << "\t max. res speed R: " << maxResonantSpeed_R << endl;
	 simClasses.logger << write;
      }
      
      initialized = success;
      return initialized;
   }

}
