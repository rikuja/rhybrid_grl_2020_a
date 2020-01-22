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

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_operator_wavelength_mesh.h"
#include "sep_operator_wavelength_mesh_particles.h"

using namespace std;
   
namespace sep {
         
   extern sep::SimControl simControl;
   extern sep::FieldsContainer fieldsContainer;
   const string PREFIX="WavelengthMeshParticleOperator";
   
   WavelengthMeshParticles::WavelengthMeshParticles(): SpatialSliceOP() {
      #if PROFILE_LEVEL > 0
         lagrangianCopy = -1;
         particleCopy = -1;
         totalTime = -1;
      #endif
      
      finalize();
   }
   
   WavelengthMeshParticles::~WavelengthMeshParticles() {
      finalize();
   }
   
   bool WavelengthMeshParticles::addConfigFileItems(ConfigReader& cr) {
      cr.add(PREFIX+".all_species","If yes, all particle species are included, 'yes' or 'no' (string).",string("yes"));
      cr.addComposed(PREFIX+".species","Names of included particle species");
      return true;
   }
   
   bool WavelengthMeshParticles::finalize() {
      includedSpecies.clear();
      includeAllSpecies = true;
      return true;
   }
   
   std::string WavelengthMeshParticles::getName() const {
      return "WavelengthMeshParticlesOperator";
   }
   
   bool WavelengthMeshParticles::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      // Init base class:
      bool success = true;
      if (SpatialSliceOP::initialize(cr,sim,simClasses) == false) return false;

      addConfigFileItems(cr);
      cr.parse();
      
      string allSpeciesString;
      cr.get(PREFIX+".all_species",allSpeciesString);
      if (allSpeciesString == "no") includeAllSpecies = false;
      
      // Get names of included species, if necessary:
      if (includeAllSpecies == false) {
	 cr.get(PREFIX+".species",includedSpecies);
      }
      return success;
   }
   
   bool WavelengthMeshParticles::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      bool success = true;
      if (baseClassInitialized == false) {
	 simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Base class is not initialized" << endl << write;
	 return false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::start("WLM Particles Total",totalTime);
      #endif
      
      for (size_t species=0; species<particleLists.size(); ++species) {
	 // Check if species is included in plot:
	 const string& speciesName = particleLists[species]->getName();
	 if (includeAllSpecies == false) {
	    bool found = false;
	    for (size_t s=0; s<includedSpecies.size(); ++s) {
	       if (speciesName == includedSpecies[s]) {
		  found = true;
		  break;
	       }
	    }
	    if (found == false) continue;
	 }

	 // Get particle data:
	 pargrid::DataID speciesDataID;
	 if (particleLists[species]->getParticles(speciesDataID) == false) {
	    simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Failed to get species '" << speciesName << "' data" << endl;
	    success = false; continue;
	 }

	 if (particleLists[species]->getSpeciesType() == simControl.particleSpeciesTypename) {
	    // Write particle species:
	    
	    pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	    if (wrapper.valid() == false) {
	       simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Received invalid data from ParGrid for species '" << speciesName << "'" << endl;
	       success = false; continue;
	    }
	    Particle<Real>** particles = wrapper.data();
	    const pargrid::ArraySizetype* const sizes = wrapper.size();

	    // Get species data:
	    const Species* speciesData = reinterpret_cast<const Species*>( particleLists[species]->getSpecies() );
	    if (speciesData == NULL) {
	       simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Failed to get species '" << speciesName << "' struct" << endl;
	       success = false; continue;
	    }

	    for (uint64_t slice=0; slice<getNumberOfSlices(); ++slice) {
	       if (writeParticles(slice,speciesName,speciesData,sizes,particles) == false) {
		  simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Failed to write species '" << speciesName << "' on slice #" << slice << endl;
		  success = false;
	       }
	    }
	 } else if (particleLists[species]->getSpeciesType() == simControl.lagrangianSpeciesTypename) {
	    // Write Lagrangian mesh:
	    
	    pargrid::DataWrapper<LagrangianParticle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<LagrangianParticle<Real> >(speciesDataID);
	    if (wrapper.valid() == false) {
	       simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Received invalid data from ParGrid for Lagrangian species '" << speciesName << "'" << endl;
	       success = false; continue;
	    }
	    LagrangianParticle<Real>** particles = wrapper.data();
	    const pargrid::ArraySizetype* const sizes = wrapper.size();
	    
	    // Get species data:
	    const LagrangianSpecies* speciesData = reinterpret_cast<const LagrangianSpecies*>( particleLists[species]->getSpecies() );
	    if (speciesData == NULL) {
	       simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Failed to get Lagrangian species '" << speciesName << "' struct" << endl;
	       success = false; continue;
	    }
	    
	    for (uint64_t slice=0; slice<getNumberOfSlices(); ++slice) {
	       if (writeLagrangian(slice,speciesName,speciesData,sizes,particles) == false) {
		  simClasses->logger << "(SEP OP WLM PARTICLE) ERROR: Failed to write Lagrangian species '" << speciesName << "' on slice #" << slice << endl;
		  success = false;
	       }
	    }
	 }
	 
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      
      return success;
   }

   bool WavelengthMeshParticles::writeLagrangian(uint64_t slice,const std::string& speciesName,const LagrangianSpecies* species,
						 const pargrid::ArraySizetype* sizes,LagrangianParticle<Real>** particles) {
      const uint8_t slicedCoordinate = getSlicedCoordinate(slice);
      
      stringstream ss;
      ss << "WaveMesh" << slice;
      const string sliceName = ss.str();
      
      vector<Real> particlesOut;
      vector<pargrid::CellID> newBlockGIDs;
      vector<pargrid::CellID> blockLIDs;
      prepareSlice(slice,simControl.N_wavelengthMeshCells);
      getAcceptedBlocks(slice,newBlockGIDs,blockLIDs);
      
      #if PROFILE_LEVEL > 0
         profile::start("Copy Lagrangian To Buffer",lagrangianCopy);
      #endif
      
      // Allocate memory for output buffer:
      size_t N_particles = 0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];
	 N_particles += sizes[blockLID];
      }
      particlesOut.resize(N_particles*3);
      
       size_t counter = 0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];
	       
	 for (pargrid::ArraySizetype p=0; p<sizes[blockLID]; ++p) {
	    // Add particle's coordinates to output buffer:
	    const Real X = particles[blockLID][p].state[lagr::XCRD];
	    const Real Y = particles[blockLID][p].state[lagr::YCRD];
	    const Real Z = particles[blockLID][p].state[lagr::ZCRD];
	    const Real LAMBDA = particles[blockLID][p].state[lagr::LAMBDA];
	    const uint32_t I = static_cast<uint32_t>(X);
	    const uint32_t J = static_cast<uint32_t>(Y);
	    const uint32_t K = static_cast<uint32_t>(Z);
	    //const uint32_t L = static_cast<uint32_t>(LAMBDA);
	    particlesOut[3*counter+0] = sim->x_crds_node[I] + (X-I)*sim->dx_cell[I];
	    particlesOut[3*counter+1] = sim->y_crds_node[J] + (Y-J)*sim->dy_cell[J];
	    particlesOut[3*counter+2] = sim->z_crds_node[K] + (Z-K)*sim->dz_cell[K];

	    // Overwrite sliced coordinate with wavelength:
	    Real lambda = (*simControl.getPhysicalWavelength)(LAMBDA);
	    if (lambda < 0.0) lambda = -log10(-lambda);
	    else lambda = log10(lambda);
	    
	    if (isCylindrical(slice) == true) {
	       particlesOut[3*counter+slicedCoordinate] = particlesOut[3*counter+2];
	       particlesOut[3*counter+2] = lambda;
	    } else {
	       particlesOut[3*counter+slicedCoordinate] = lambda;
	    }
	    ++counter;
	 }
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      
      // Write particles to output file:
      N_particles = particlesOut.size() / 3;
      map<string,string> attribs;
      attribs["name"] = sliceName + "/" + speciesName;
      attribs["type"] = vlsv::mesh::STRING_POINT;
      if (isCylindrical(slice) == true) {
	 attribs["geometry"] = vlsv::getMeshGeometry(vlsv::geometry::CYLINDRICAL);
      } else {
	 attribs["geometry"] = vlsv::getMeshGeometry(vlsv::geometry::CARTESIAN);
      }
      
      if (simClasses->vlsv.writeArray("MESH",attribs,N_particles,3,&(particlesOut[0])) == false) {
	 return false;
      }
      return true;
   }
   
   bool WavelengthMeshParticles::writeParticles(uint64_t slice,const std::string& speciesName,const Species* species,
						const pargrid::ArraySizetype* sizes,Particle<Real>** particles) {
      const uint8_t slicedCoordinate = getSlicedCoordinate(slice);
      
      stringstream ss;
      ss << "WaveMesh" << slice;
      const string sliceName = ss.str();

      Real B[3];
      Real V_wave[3];
      Real dV_wave;

      vector<Real> particlesOut;
      vector<pargrid::CellID> newBlockGIDs;
      vector<pargrid::CellID> blockLIDs;
      prepareSlice(slice,simControl.N_wavelengthMeshCells);
      getAcceptedBlocks(slice,newBlockGIDs,blockLIDs);

      #if PROFILE_LEVEL > 0
         profile::start("Copy Particles To Buffer",particleCopy);
      #endif
      
      // Allocate memory for output buffer:
      size_t N_particles = 0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];
	 N_particles += sizes[blockLID];
      }
      particlesOut.resize(N_particles*3);

      simControl.alfvenSign = +1;
      size_t counter = 0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 const pargrid::CellID blockLID = blockLIDs[b];

	 for (pargrid::ArraySizetype p=0; p<sizes[blockLID]; ++p) {
	    // Add particle's coordinates to output buffer:
	    const Real X = particles[blockLID][p].state[sep::particle::XCRD];
	    const Real Y = particles[blockLID][p].state[sep::particle::YCRD];
	    const Real Z = particles[blockLID][p].state[sep::particle::ZCRD];
	    const uint32_t I = static_cast<uint32_t>(X);
	    const uint32_t J = static_cast<uint32_t>(Y);
	    const uint32_t K = static_cast<uint32_t>(Z);
	    particlesOut[3*counter+0] = sim->x_crds_node[I] + (X-I)*sim->dx_cell[I];
	    particlesOut[3*counter+1] = sim->y_crds_node[J] + (Y-J)*sim->dy_cell[J];
	    particlesOut[3*counter+2] = sim->z_crds_node[K] + (Z-K)*sim->dz_cell[K];

	    // Calculate resonant wavelength:
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,particles[blockLID][p].state,B,V_wave,dV_wave,V_alfven,simControl.alfvenSign);
	    const Real B_mag = vectorMagnitude<3>(B);
	    Real V_wave_par = dotProduct<3>(V_wave,B)/B_mag;

	    // Calculate resonant wavelength in chosen wave rest frame:
	    Real lambda_res = 2*M_PI/species->q_per_m/B_mag * (particles[blockLID][p].state[sep::particle::V_PAR] - V_wave_par);
	    if (lambda_res < 0.0) lambda_res = -log10(-lambda_res);
	    else lambda_res = log10(lambda_res);

	    // Overwrite sliced coordinate with resonant wavelength:
	    if (isCylindrical(slice) == true) {
	       particlesOut[3*counter+slicedCoordinate] = particlesOut[3*counter+2];
	       particlesOut[3*counter+2] = lambda_res;
	    } else {
	       particlesOut[3*counter+slicedCoordinate] = lambda_res;
	    }
	    ++counter;
	 }
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif

      // Write particles to output file:
      N_particles = particlesOut.size() / 3;
      map<string,string> attribs;
      attribs["name"] = sliceName + "/" + speciesName;
      attribs["type"] = vlsv::mesh::STRING_POINT;      
      if (isCylindrical(slice) == true) {
	 attribs["geometry"] = vlsv::getMeshGeometry(vlsv::geometry::CYLINDRICAL);
      } else {
	 attribs["geometry"] = vlsv::getMeshGeometry(vlsv::geometry::CARTESIAN);
      }
      
      if (simClasses->vlsv.writeArray("MESH",attribs,N_particles,3,&(particlesOut[0])) == false) {
	 return false;
      }
      return true;
   }

} // namespace sep
