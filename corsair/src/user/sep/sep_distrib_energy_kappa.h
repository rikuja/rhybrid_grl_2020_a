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

#ifndef SEP_DISTRIP_ENERGY_KAPPA_H
#define SEP_DISTRIB_ENERGY_KAPPA_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_distrib_energy_container.h"

namespace sep {
   
   struct EnergyDistribKappa {
      EnergyDistribKappa();

      Real energyIntegral;        /**< Integral of kappa distribution over the chosen 
				   * injection energy range.*/
      Real injectionEnergyMin;    /**< Minimum injection energy per amu.*/
      Real injectionEnergyMax;    /**< Maximum injection energy per amu.*/
      Real kappaIndex;            /**< Kappa-index of distribution.*/
      Real normalizationConstant; /**< Kappa-distribution normalization constant 
				   * without thermal speed. Equals to gamma(kappa+1) /
				   * gamma(kappa-1/2) / pi^(1.5).*/
      Real thermalConstant;       /**< Conversion factor from Maxwell-Boltzmann thermal
				   * energy to kappa-distribution thermal energy.*/
      
      bool initialized;
      bool binnedIntervals;       /**< If true, then injection energy is completely random 
				   * over the injection energy range. Otherwise it is random
				   * over binned intervals.*/
      SimulationClasses* simClasses;
   };

   bool energyDistribKappaFinalize(void*& parameters);
   bool energyDistribKappaGetDistrib(distrib::InjectionEnergy& energy,void* parameters);
   bool energyDistribKappaGetEnergy(distrib::InjectionEnergy& injectionEnergy,void* parameters);
   bool energyDistribKappaInitialize(Simulation& sim,SimulationClasses& simClasses,
				     ConfigReader& cr,const std::string& regionName,
				     void*& parameters);

} // namespace sep

#endif
