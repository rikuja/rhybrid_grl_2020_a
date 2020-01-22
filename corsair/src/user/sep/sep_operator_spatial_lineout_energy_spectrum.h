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

#ifndef SEP_OPERATOR_SPATIAL_LINEOUT_ENERGY_SPECTRUM_H
#define SEP_OPERATOR_SPATIAL_LINEOUT_ENERGY_SPECTRUM_H

#include <vector>

#include <dataoperator.h>

#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_operator_spatial_lineout.h"

namespace sep {
   
   class SpatialLineoutEnergySpectrum: public OperatorSpatialLineout {
    public:
      SpatialLineoutEnergySpectrum();
      ~SpatialLineoutEnergySpectrum();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

    private:
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif
      
      Real d_energy;                              /**< Log10 of energy bin width.*/
      std::vector<Real> energyNodeCoordinates;    /**< Node coordinates of the energy dimension in SI units.*/
      std::vector<Real> energyNodeCoordinatesOut; /**< Node coordinates of the energy dimension in log10 eV.*/
      std::vector<std::string> speciesNames;      /**< Names of accumulated particle species or 'all',
						   * if all species are accumulated to same array.*/
      
      bool addConfigFileItems(ConfigReader& cr);
      bool addParticles(const sep::Species* species,const ParticleListBase* particleList,
			std::vector<pargrid::CellID>& inner,std::vector<pargrid::CellID>& boundary,
			uint8_t lineoutCrd,std::vector<Real>& f);
      void addParticleWeight(pargrid::CellID blockLID,uint8_t lineoutCrd,size_t counter,std::vector<Real>& f,
			     const sep::Species* species,const sep::Particle<Real>& particle);
      Real getCellVolume(pargrid::CellID blockLID,uint8_t lineoutCoordinate,int thirdIndex) const;
      uint32_t getGlobalIndex(uint32_t i,uint32_t j) const;
      bool writeLineoutMesh(size_t line,std::vector<pargrid::CellID>& inner,std::vector<pargrid::CellID>& boundary);
   };
   
} // namespace sep

#endif
