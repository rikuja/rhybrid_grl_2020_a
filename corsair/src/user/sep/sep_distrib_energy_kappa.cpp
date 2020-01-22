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

#include <main.h>
#include <constants.h>
#include "sep_distrib_energy_container.h"
#include "sep_distrib_energy_kappa.h"

using namespace std;

namespace sep {

   void runTests(Simulation& sim,SimulationClasses& simClasses);
   
   EnergyDistribKappa::EnergyDistribKappa() {
      initialized = false;
      injectionEnergyMin = NAN;
      injectionEnergyMax = NAN;
      kappaIndex = NAN;
      simClasses = NULL;
      binnedIntervals = false;
   }

   // Calculate normalization constant of kappa energy distribution:
   Real calculateNormalizationConstant(Real kappa) {
      return tgamma(kappa+1) / (tgamma(kappa-0.5)*tgamma(1.5));
   }

   bool energyDistribKappaFinalize(void*& parameters) {
      EnergyDistribKappa* e = reinterpret_cast<EnergyDistribKappa*>(parameters);
      delete e; e = NULL;
      return true;
   }
   
   bool energyDistribKappaGetDistrib(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      EnergyDistribKappa& energyKappa = *reinterpret_cast<EnergyDistribKappa*>(parameters);

      // Divide energy by kappa index and kappa distribution thermal energy:
      const Real kappaThermalEnergy = energyKappa.kappaIndex * energyKappa.thermalConstant * injectionEnergy.thermalEnergy;
      const Real x = injectionEnergy.energy / kappaThermalEnergy;

      // Return value of energy distribution:
      injectionEnergy.weight = energyKappa.normalizationConstant * sqrt(x) * pow(1.0 + x,-energyKappa.kappaIndex-1.0) / kappaThermalEnergy;
      return true;
   }

   Real getDistributionMB(const distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      Real x = injectionEnergy.energy / injectionEnergy.thermalEnergy;
      Real normalization = sqrt(4.0/M_PI) / injectionEnergy.thermalEnergy;
      return normalization * sqrt(x) * exp(-x);
   }

   /** Get injection energy and value of distribution function.
    * @param injectionEnergy Injection energy per amu (in Joules) is written here.
    * @param weight Value of distribution function is written here.
    * @return If true, energy and weight contain sensible values.*/
   bool energyDistribKappaGetEnergy(distrib::InjectionEnergy& injectionEnergy,void* parameters) {
      EnergyDistribKappa& energyKappa = *reinterpret_cast<EnergyDistribKappa*>(parameters);
      
      if (energyKappa.binnedIntervals == true) {
	 // Width of this energy bin:
	 const Real dU = (energyKappa.injectionEnergyMax-energyKappa.injectionEnergyMin) / injectionEnergy.N_particles;
	 const Real dU_super = dU / injectionEnergy.superResolution;

	 // Min,max x-values of this bin:
	 const Real x_min = energyKappa.injectionEnergyMin + injectionEnergy.N*dU +   injectionEnergy.N_super  *dU_super; 
	 const Real x_max = energyKappa.injectionEnergyMin + injectionEnergy.N*dU + (injectionEnergy.N_super+1)*dU_super;

	 // Get distribution function values at min,max energies of this bin:
	 injectionEnergy.energy = x_min * injectionEnergy.thermalEnergy;
	 energyDistribKappaGetDistrib(injectionEnergy,parameters);
	 Real f_min = injectionEnergy.weight;
	 
	 injectionEnergy.energy = 0.5*(x_min + x_max)*injectionEnergy.thermalEnergy;
	 energyDistribKappaGetDistrib(injectionEnergy,parameters);
	 Real f_mid = injectionEnergy.weight;
	 
	 injectionEnergy.energy = x_max * injectionEnergy.thermalEnergy;
	 energyDistribKappaGetDistrib(injectionEnergy,parameters);
	 Real f_max = injectionEnergy.weight;

	 // Get distribution function value at injection energy:
	 Real x = energyKappa.injectionEnergyMin + injectionEnergy.N*dU 
	        + (injectionEnergy.N_super+energyKappa.simClasses->random.uniform())*dU_super;
	 injectionEnergy.energy = x * injectionEnergy.thermalEnergy;

	 // Integrate distribution using Simpson's rule:
	 injectionEnergy.weight = (f_min + 4*f_mid + f_max)*(x_max-x_min)/6.0*injectionEnergy.thermalEnergy;
	 
	 return true;
      }

      // Injection energy has uniform distribution in
      // interval [injectionEnergyMin,injectionEnergyMax]. Interval 
      // is given in units of Maxwell-Boltzmann thermal energy:
      Real x = energyKappa.injectionEnergyMin
	+ (energyKappa.injectionEnergyMax-energyKappa.injectionEnergyMin)
	* energyKappa.simClasses->random.uniform();

      // Calculate physical injection energy (multiply by Boltzmann's constant and temperature):
      injectionEnergy.energy = x * injectionEnergy.thermalEnergy;
      
      // Calculate value of kappa energy distribution at 
      // generated random injection energy:
      return energyDistribKappaGetDistrib(injectionEnergy,parameters);
   }

   bool energyDistribKappaInitialize(Simulation& sim,SimulationClasses& simClasses,
				     ConfigReader& cr,const std::string& regionName,
				     void*& parameters) {
      simClasses.logger << "(ENERGY DISTRIB KAPPA) Starting initialization" << std::endl;
      
      EnergyDistribKappa* energyKappa = new EnergyDistribKappa();
      energyKappa->initialized = true;
      parameters = energyKappa;
      
      energyKappa->simClasses = &simClasses;
      
      // Read injection energy from config file:
      const Real DEF_VALUE = numeric_limits<Real>::infinity();
      cr.add(regionName+".MB_relative_energy_min","Minimum injection energy relative to MB thermal energy per amu (float).",DEF_VALUE);
      cr.add(regionName+".MB_relative_energy_max","Maximum injection energy relative to MB thermal energy per amu (float).",DEF_VALUE);
      cr.add(regionName+".spectral_index","Kappa index of speed distribution (float)",DEF_VALUE);
      cr.add(regionName+".binned_intervals","If 'yes' injection energy is random over binned intervals (string)",string("no"));

      string energyUnits,binnedIntervalsString;
      cr.parse();
      cr.get(regionName+".MB_relative_energy_min",energyKappa->injectionEnergyMin);
      cr.get(regionName+".MB_relative_energy_max",energyKappa->injectionEnergyMax);
      cr.get(regionName+".spectral_index",energyKappa->kappaIndex);
      cr.get(regionName+".binned_intervals",binnedIntervalsString);
      
      if (binnedIntervalsString == "yes") energyKappa->binnedIntervals = true;
      if (energyKappa->injectionEnergyMin == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".injection_energy_min' was not found" << endl;
	 energyKappa->initialized = false;
      }
      if (energyKappa->injectionEnergyMax == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".injection_energy_max' was not found" << endl;
	 energyKappa->initialized = false;
      }
      if (energyKappa->kappaIndex == DEF_VALUE) {
	 simClasses.logger << "\t ERROR: Parameter '" << regionName+".spectral_index' was not found" << endl;
	 energyKappa->initialized = false;
      }

      // Conversion factor from Maxwell-Boltzmann thermal 
      // energy to kappa thermal energy:
      energyKappa->thermalConstant = (2*energyKappa->kappaIndex-3)/(2*energyKappa->kappaIndex);

      // Calculate kappa distribution normalization constant:
      energyKappa->normalizationConstant =
	calculateNormalizationConstant(energyKappa->kappaIndex);

      /*
      // Calculate integral of distribution function over given energy range:
      if (energyKappa->initialized == true) {
	 const uint32_t N = 10000000;
	 const Real T = 1e4;
	 distrib::InjectionEnergy injEnergy;
	 injEnergy.thermalEnergy = T * constants::BOLTZMANN;

	 Real integral = 0.0;
	 const Real relEnergyInterval = energyKappa->injectionEnergyMax - energyKappa->injectionEnergyMin;
	 const Real dU = relEnergyInterval / (N-1);
	 for (uint32_t i=0; i<N; ++i) {
	    injEnergy.energy = (energyKappa->injectionEnergyMin +   i  *dU) * injEnergy.thermalEnergy;
	    energyDistribKappaGetDistrib(injEnergy,energyKappa);
	    Real F_min = injEnergy.weight;
	    injEnergy.energy = (energyKappa->injectionEnergyMax + (i+1)*dU) * injEnergy.thermalEnergy;
	    Real F_max = injEnergy.weight;
	    
	    integral += F_min + 0.5*(F_max-F_min);
	 }
	 energyKappa->integral = integral * dU*injEnergy.thermalEnergy;
      }*/

      // Run tests (if requested):
      if (energyKappa->initialized == true) {
	 if (sim.runTests == true) runTests(sim,simClasses);
      }

      // Write init status and exit:
      simClasses.logger << "\t Initialization completed, status is ";
      if (energyKappa->initialized == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;      
      return energyKappa->initialized;
   }
   
   void runTests(Simulation& sim,SimulationClasses& simClasses) {
      // Only master process runs tests:
      if (sim.mpiRank != sim.MASTER_RANK) return;
      
      // Run tests assuming 1 MK temperature:
      Real T_ion = 1.0e5;
      distrib::InjectionEnergy injEnergy;
      injEnergy.thermalEnergy = constants::BOLTZMANN*T_ion;
      
      // Integrate kappa distributions over [0-100] MB thermal energy range:
      EnergyDistribKappa energyKappa;
      energyKappa.simClasses = &simClasses;
      energyKappa.injectionEnergyMin = 0.0;
      energyKappa.injectionEnergyMax = 10000.0;
      
      // Number of data points used to plot distributions:
      const int N = 50000;
      
      // Store plotted kappa values:
      vector<Real> kappaValues;
      kappaValues.push_back(2.0);
      kappaValues.push_back(3.0);
      kappaValues.push_back(4.0);
      for (size_t k=0; k<kappaValues.size(); ++k) {
	 energyKappa.kappaIndex = kappaValues[k];
	 energyKappa.thermalConstant = (2*kappaValues[k]-3)/(2*kappaValues[k]);
	 const Real kappaThermalEnergy = energyKappa.kappaIndex*energyKappa.thermalConstant*injEnergy.thermalEnergy;

	 // Min,max injection energies are given in units of MB
	 // thermal energy. Compute min,max kappa distribution energies:
	 const Real minKappaEnergy = energyKappa.injectionEnergyMin/(kappaValues[k]*energyKappa.thermalConstant);
	 const Real maxKappaEnergy = energyKappa.injectionEnergyMax/(kappaValues[k]*energyKappa.thermalConstant);
	 
	 // Calculate normalization constant:
	 energyKappa.normalizationConstant =
	   calculateNormalizationConstant(energyKappa.kappaIndex);
	 
	 // Open output file:
	 stringstream ss;
	    ss << "test_energy_distrib_kappa_" << kappaValues[k] << ".txt";
	 fstream out(ss.str(),fstream::out);

	 // Write header:
	 out << "# Thermal energy of this distribution is: ";
	 out << injEnergy.thermalEnergy / constants::CHARGE_ELEMENTARY;
	 out << " eV" << endl;
	 out << "# Energy Rel. to MB (min) : " << energyKappa.injectionEnergyMin << endl;
	 out << "# Energy Rel. to MB (max) : " << energyKappa.injectionEnergyMax << endl;

	 out << "# output values are:" << endl;
	 out << "1: energy / Maxwell-Boltzmann thermal energy" << endl;
	 out << "2: energy (eV)" << endl;
	 out << "3: value of kappa distribution" << endl;
	 out << "4: value of Maxwell-Boltzmann distribution" << endl;
	 
	 const Real relEnergyInterval = energyKappa.injectionEnergyMax - energyKappa.injectionEnergyMin;
	 const Real dEnergy = relEnergyInterval / (N-1);
	 for (int i=0; i<N; ++i) {
	    injEnergy.energy = (energyKappa.injectionEnergyMin + i*dEnergy) * injEnergy.thermalEnergy;
	    energyDistribKappaGetDistrib(injEnergy,&energyKappa);
	    Real value = injEnergy.weight;
	    Real valueMB = getDistributionMB(injEnergy,&energyKappa);

	    out << injEnergy.energy/injEnergy.thermalEnergy << '\t';
	    out << injEnergy.energy/constants::CHARGE_ELEMENTARY << '\t';
	    out << value << '\t';
	    out << valueMB;
	    out << endl;
	 }

	 // Integrate kappa energy distribution using Monte Carlo.
	 // The result should be a value close to unity:
	 const uint32_t N_points = 1000000;
	 Real sum = 0.0;
	 for (uint32_t i=0; i<N_points; ++i) {
	    energyDistribKappaGetEnergy(injEnergy,&energyKappa);
	    sum += injEnergy.weight;
	 }
	 sum /= (1.0*N_points);
	 sum *= (maxKappaEnergy-minKappaEnergy)*kappaThermalEnergy;
	 out << "# Monte-Carlo integral of kappa-distribution is: " << sum << endl;
	 
	 out.close();
      }
   }
   
} // namespace sep
