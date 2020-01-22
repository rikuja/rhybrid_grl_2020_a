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

#ifndef SEP_OPERATOR_SHOCK_ELECTRONS_H
#define SEP_OPERATOR_SHOCK_ELECTRONS_H

#include <vector>
#include <dataoperator.h>

#include "sep_distrib_energy_container.h"
#include "sep_operator_energy_channels.h"

namespace sep {
   
   class OperatorShockElectrons: public DataOperator {
    public:
      OperatorShockElectrons();
      ~OperatorShockElectrons();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool runTests();
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

    private:
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif

      bool testsRun;
      
      void* energyParams;                            /**< Parameters for energy distribution.*/
      sep::finalizeEnergyDistrib finalizeEnergy;     /**< Pointer to energy distribution finalizer.*/
      sep::getEnergyDistribFunction getDistrib;      /**< Pointer to energy distribution.*/
      sep::getEnergyFunction getEnergy;              /**< Pointer to energy distribution function.*/
      std::string energyDistribString;
      std::string energyDistribParamsString;

      Real d_V_par;
      Real V_injection_max;
      Real V_par_min;
      Real V_par_max;
      int32_t xDistribCells;
      int32_t yDistribCells;
      uint32_t N_speedSamples;
      uint32_t N_pitchSamples;

      uint32_t N_shockMeshCells;
      
      std::vector<spacecraft::Instrument> instruments;
      std::vector<std::vector<Real> > instrumentData;
      std::vector<std::vector<Real> > longitudeInstrData;
      std::vector<Real> surfaceAreaSum;
      std::vector<bool> filesCreated;

      bool addConfigFileItems(ConfigReader& cr,const std::string& region);
      int32_t parDistribIndex(const int32_t& i,const int32_t& j);

      bool solveReflection(std::vector<Real>& distributionThermal,
			   std::vector<Real>& distributionReflected,
			   std::vector<Real>& distributionTransmitted,
			   std::vector<Real>& reflectionEfficiency);
      bool solveReflection(Real* centroid,int32_t phiIndex,uint32_t zoneIndex,Real numberDensity,Real temperature,Real R_gas,Real R_magn,
			   Real* V_shock_SIM,Real* shockNormal,Real* V_plasma_SIM,Real* B,
			   std::vector<std::vector<Real> >& distributionThermal,
			   std::vector<std::vector<Real> >& distributionReflected,
			   std::vector<Real>& distributionTransmitted,
			   std::vector<std::vector<Real> >& reflectionEfficiency);
       
      bool writeDistributionMesh(const std::string& meshName);
      bool writeElectronDistribution(const std::vector<Real>& nodeCoords,const std::string& meshName);
   };
   
} // namespace sep

#endif
