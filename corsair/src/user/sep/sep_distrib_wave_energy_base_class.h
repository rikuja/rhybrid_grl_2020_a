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

#ifndef SEP_DISTRIB_WAVE_ENERGY_BASE_H
#define SEP_DISTRIB_WAVE_ENERGY_BASE_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace sep {

   /** Definition of base class for Alfven wave energy spectra. 
    * Derived classes are used to calculate the energy carried by 
    * antiparallel or parallel propagating wave packets when they 
    * are injected to simulation at inflow boundaries.*/
   class WaveEnergySpectrumBaseClass {
    public:
      WaveEnergySpectrumBaseClass();
      virtual ~WaveEnergySpectrumBaseClass();
      
      /** Pure virtual class finalizer function. Must deallocate all member variables 
       * in derived class(es). Finalizer function must be safe to call multiple times.
       * @return It true, class finalized successfully.*/
      virtual bool finalize() =0;
      
      /** Pure virtual function that returns energy density in Alfven waves that have wavelenghts
       * in interval [minLambda,maxLambda].
       * @param B0_mag Magnitude of ambient magnetic field.
       * @param minLambda Lower limit for requested wavelength interval.
       * @param maxLambda Upper limit for requested wavelength interval.
       * @return Energy in Alfven waves in requested wavelength interval.
       * @see WaveEnergySpectrumBaseClass::getIntensity.*/
      virtual Real getEnergyDensity(Real B0_mag,Real minLambda,Real maxLambda) =0;
      
      /** Pure virtual function that returns the intensity of Alfven waves at given wavelength.
       * Integral of intensity is returned by member function getEnergy.
       * @param B0_mag Magnitude of ambient magnetic field.
       * @param wavelength Wavelength.
       * @return Intensity at given wavelength.
       * @see WaveEnergySpectrumBaseClass::getEnergy.*/
      virtual Real getIntensity(Real B0_mag,Real wavelength) =0;

      virtual Real getIntensityDerivated(Real B0_mag,Real wavelength) =0;
      
      virtual Real getMaximumWavelength(Real B0_mag,int helicity) const =0;
      
      virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName);
      
    protected:
      bool initialized;              /**< If true, (derived) class initialized successfully and is ready for use.*/
      Simulation* sim;               /**< Pointer to struct that contains variables controlling execution of simulation.*/
      SimulationClasses* simClasses; /**< Pointer to struct that contains generic use simulation classes.*/
   };
   
} // namespace sep

#endif
