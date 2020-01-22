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

#ifndef SEP_OPERATOR_WAVELENGTH_MESH_PARTICLE_H
#define SEP_OPERATOR_WAVELENGTH_MESH_PARTICLE_H

#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_lagr_definition.h"
#include "sep_lagr_species.h"
#include "sep_operator_spatial_slice.h"
#include "sep_fields_container.h"

namespace sep {
   
   class WavelengthMeshParticles: public SpatialSliceOP {
    public:
      WavelengthMeshParticles();
      ~WavelengthMeshParticles();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);
      
    private:
      bool includeAllSpecies;
      std::vector<std::string> includedSpecies;
      
      bool addConfigFileItems(ConfigReader& cr);
      
      bool writeLagrangian(uint64_t slice,const std::string& speciesName,const LagrangianSpecies* species,
			   const pargrid::ArraySizetype* sizes,LagrangianParticle<Real>** particles);
      bool writeParticles(uint64_t slice,const std::string& speciesName,const Species* species,
			  const pargrid::ArraySizetype* sizes,Particle<Real>** particles);
      
      #if PROFILE_LEVEL > 0
         int lagrangianCopy;
         int particleCopy;
         int totalTime;
      #endif
   };   
   
}

#endif
