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

#include "sep_simcontrol.h"
#include "sep_operator_spatial_lineout_energy_spectrum.h"
#include "sep_guiding_center_theory.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   const static string PREFIX="EnergySpectrumLineout";
   const static Real DEF_REAL = numeric_limits<Real>::infinity();
   
   SpatialLineoutEnergySpectrum::SpatialLineoutEnergySpectrum(): OperatorSpatialLineout() { 
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   SpatialLineoutEnergySpectrum::~SpatialLineoutEnergySpectrum() { }
   
   bool SpatialLineoutEnergySpectrum::addConfigFileItems(ConfigReader& cr) {
      cr.add(PREFIX+".bins","Number of energy bins in spectrum (int)",uint32_t(0));
      cr.add(PREFIX+".energy_units","Units in which min/max energies are given eV/keV/MeV/GeV (string)",string(""));
      cr.add(PREFIX+".energy_min","Minimum energy in chosen units (float)",DEF_REAL);
      cr.add(PREFIX+".energy_max","Maximum energy in chosen units (float)",DEF_REAL);
      cr.addComposed(PREFIX+".species","Names of particle species included in energy spectrum (string)");
      return true;
   }

   void SpatialLineoutEnergySpectrum::addParticleWeight(pargrid::CellID blockLID,uint8_t lineoutCrd,size_t counter,std::vector<Real>& f,
							const sep::Species* species,const sep::Particle<Real>& particle) {
      // Get fields at particle position:
      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(blockLID,sim->t,particle.state,E,B,gradB);
      
      // Calculate GC drift velocity:
      Real V_drift[3];
      calculateDriftVelocity(*species,particle,E,B,gradB,V_drift);
      
      // Calculate energy in log(eV) per amu:
      Real energy = particle.state[particle::V_PAR]*particle.state[particle::V_PAR] + vectorMagnitude2<3>(V_drift);
      energy = 0.5*constants::MASS_PROTON*energy;
      energy += particle.state[particle::MU]*vectorMagnitude<3>(B);
      energy = log10(energy / constants::CHARGE_ELEMENTARY);

      pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      int i_block,j_block,k_block;
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

      // Check that energy is not outside array bounds:
      if (energy < energyNodeCoordinatesOut[0] || energy >= energyNodeCoordinatesOut[energyNodeCoordinatesOut.size()-1]) {
	 return;
      }

      // Add particle weight to distribution:
      #warning IMPROVE ME: Energy spectrum only supports NGP for now:
      uint32_t blockWidths[3];
      blockWidths[0] = block::WIDTH_X;
      blockWidths[1] = block::WIDTH_Y;
      blockWidths[2] = block::WIDTH_Z;
      const uint32_t i_cell  = static_cast<uint32_t>(particle.state[lineoutCrd]) % blockWidths[lineoutCrd];
      const int32_t bin = static_cast<uint32_t>((energy - energyNodeCoordinatesOut[0]) / d_energy);
      f[getGlobalIndex(bin,counter*block::WIDTH_X+i_cell)] += particle.state[particle::WEIGHT];
   }

   bool SpatialLineoutEnergySpectrum::addParticles(const sep::Species* species,const ParticleListBase* particleList,
						   std::vector<pargrid::CellID>& inner,std::vector<pargrid::CellID>& boundary,
						   uint8_t lineoutCrd,std::vector<Real>& f) {
      pargrid::DataID speciesID;
      if (particleList->getParticles(speciesID) == false) return false;
      
      typedef Particle<Real> PARTICLE;
      const pargrid::DataWrapper<PARTICLE>& wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(speciesID);

      size_t counter = 0;
      for (size_t b=0; b<inner.size(); ++b) {
	 const pargrid::CellID blockLID = inner[b];
	 const PARTICLE* particles = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    addParticleWeight(blockLID,lineoutCrd,counter,f,species,particles[p]);
	 }
	 ++counter;
      }

      for (size_t b=0; b<boundary.size(); ++b) {
	 const pargrid::CellID blockLID = boundary[b];
	 const PARTICLE* particles = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    addParticleWeight(blockLID,lineoutCrd,counter,f,species,particles[p]);
	 }
	 ++counter;
      }
      
      return true;
   }
   
   bool SpatialLineoutEnergySpectrum::finalize() {
      bool success = true;
      return success;
   }

   Real SpatialLineoutEnergySpectrum::getCellVolume(pargrid::CellID blockLID,uint8_t lineoutCoordinate,int thirdIndex) const {
      Real volume = 0.0;
      switch (lineoutCoordinate) {
       case 0:
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) 
	   volume += simControl.cellVolumes[blockLID*block::SIZE+block::index(thirdIndex,j,k)];
	 break;
       case 1:
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i)
	   volume += simControl.cellVolumes[blockLID*block::SIZE+block::index(i,thirdIndex,k)];
	 break;
       case 2:
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i)
	   volume += simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,thirdIndex)];
	 break;
      }
      return volume;
   }
   
   uint32_t SpatialLineoutEnergySpectrum::getGlobalIndex(uint32_t i,uint32_t j) const {
      return j*(energyNodeCoordinates.size()-1)+i;
   }
   
   std::string SpatialLineoutEnergySpectrum::getName() const {
      return "SpatialLineoutEnergySpectrum";
   }
   
   bool SpatialLineoutEnergySpectrum::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (OperatorSpatialLineout::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Base class failed to initialize" << endl << write;
	 return false;
      }
      
      if (addConfigFileItems(cr) == false) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to add configuration file options" << endl << write;
	 return false;
      }
      
      uint32_t N_energyBins = 0;
      string energyUnitsString;
      Real maxEnergy;
      Real minEnergy;
      vector<string> speciesNamesString;
      cr.parse();
      if (cr.get(PREFIX+".bins",N_energyBins) == false) success = false;
      if (cr.get(PREFIX+".energy_units",energyUnitsString) == false) success = false;
      if (cr.get(PREFIX+".energy_min",minEnergy) == false) success = false;
      if (cr.get(PREFIX+".energy_max",maxEnergy) == false) success = false;
      if (cr.get(PREFIX+".species",speciesNamesString) == false) success = false;
      if (success == false) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to read required configuration file options" << endl << write;
	 return false;
      }

      if (speciesNamesString.size() == 0) {
	 speciesNames.push_back("all");
      } else {
	 for (size_t s=0; s<speciesNamesString.size(); ++s) {
	    if (speciesNamesString[s].size() > 0) speciesNames.push_back(speciesNamesString[s]);
	 }
      }

      if (N_energyBins == 0) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) WARNING: Number of energy bins is zero" << endl << write;
	 return true;
      }
      
      if (minEnergy == DEF_REAL) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Could not read minimum energy." << endl;
	 simClasses.logger << "\t It must be given with option '" << PREFIX+".energy_min'" << endl << write;
	 success = false;
      }
      if (maxEnergy == DEF_REAL) {
	 simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Could not read minimum energy." << endl;
	 simClasses.logger << "\t It must be given with option '" << PREFIX+".energy_min'" << endl << write;
	 success = false;
      }
      
      // Create output energy node coordinates (energies in log units):
      const Real energyUnits = simClasses.constants.getEnergyInSI(energyUnitsString);
      minEnergy = log10(minEnergy * energyUnits / constants::CHARGE_ELEMENTARY);
      maxEnergy = log10(maxEnergy * energyUnits / constants::CHARGE_ELEMENTARY);
      d_energy = (maxEnergy-minEnergy)/(N_energyBins-1);
      for (uint32_t i=0; i<N_energyBins+1; ++i) {
	 energyNodeCoordinatesOut.push_back(minEnergy + i*d_energy);
      }
      
      // Get energy units conversion factor and create energy node coordinates in SI units:
      for (uint32_t i=0; i<energyNodeCoordinatesOut.size(); ++i) {
	 const Real energy = pow(10.0,energyNodeCoordinatesOut[i]*constants::CHARGE_ELEMENTARY);
	 energyNodeCoordinates.push_back(energy);
      }
      
      simClasses.logger << "(SEP OP LINEOUT ENERGY SPECTRUM) Initialization complete, status is ";
      if (success == true) simClasses.logger << "SUCCESS" << endl;
      else simClasses.logger << "FAILURE" << endl;
      simClasses.logger << write;
      
      return success;
   }
   
   bool SpatialLineoutEnergySpectrum::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      bool success = true;
      if (baseClassInitialized == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Base class is not initialized" << endl << write;
	 return false;
      }

      #if PROFILE_LEVEL > 0
         profile::start("Spatial Lineout Energy Spectrum",profTotalTime);
      #endif

      for (size_t line=0; line<getNumberOfLineouts(); ++line) {
	 // Write lineout mesh. We also get local IDs of blocks where we need to
	 // calculate energy spectrum:
	 vector<pargrid::CellID> innerLocalIDs;
	 vector<pargrid::CellID> boundaryLocalIDs;
	 if (writeLineoutMesh(line,innerLocalIDs,boundaryLocalIDs) == false) success = false;

	 for (size_t s=0; s<speciesNames.size(); ++s) {
	 
	    const uint64_t arraySize = (innerLocalIDs.size()+boundaryLocalIDs.size())*(energyNodeCoordinatesOut.size()-1);
	    vector<Real> f(arraySize);
	    for (size_t i=0; i<arraySize; ++i) f[i] = 0.0;
	 
	    const uint8_t lineoutCoordinate = getLineoutCoordinate(line);
	    
	    // If config file doesn't contain any particle species names,
	    // then all species are accumulated to the same variable called 'all'.
	    // Otherwise a separate variable is created for each species.
	    if (speciesNames[s] == "all") {
	       for (size_t list=0; list<particleLists.size(); ++list) {
		  if (particleLists[list]->getSpeciesType() == simControl.lagrangianSpeciesTypename) continue;
		  const sep::Species* species = reinterpret_cast<const sep::Species*>(particleLists[list]->getSpecies());
		  addParticles(species,particleLists[list],innerLocalIDs,boundaryLocalIDs,lineoutCoordinate,f);
	       }
	    } else {
	       const sep::Species* species = NULL;
	       bool found = false;
	       for (size_t list=0; list<particleLists.size(); ++list) {
		  if (particleLists[list]->getSpeciesName() == speciesNames[s]) {
		     species = reinterpret_cast<const sep::Species*>(particleLists[list]->getSpecies());
		     addParticles(species,particleLists[list],innerLocalIDs,boundaryLocalIDs,lineoutCoordinate,f);
		     found = true;
		     break;
		  }
	       }

	       // If species was not found, do not write variable:
	       if (found == false) continue;
	    }

	    size_t counter = 0;
	    for (size_t b=0; b<innerLocalIDs.size(); ++b) {
	       // Calculate block's spatial volume:
	       const pargrid::CellID blockLID = innerLocalIDs[b];
	       
	       int32_t blockWidth = 1;
	       if (lineoutCoordinate == 0) blockWidth = block::WIDTH_X;
	       if (lineoutCoordinate == 1) blockWidth = block::WIDTH_Y;
	       if (lineoutCoordinate == 2) blockWidth = block::WIDTH_Z;
	       
	       for (int j=0; j<blockWidth; ++j) {
		  const Real cellVolume = getCellVolume(blockLID,lineoutCoordinate,j);
		  const size_t N_bins = energyNodeCoordinatesOut.size()-1;
		  for (size_t i=0; i<N_bins; ++i) {
		     const Real U_min = pow(10.0,energyNodeCoordinatesOut[i  ]);
		     const Real U_max = pow(10.0,energyNodeCoordinatesOut[i+1]);
		     const Real dU = U_max-U_min;
		     if (f[(counter*blockWidth+j)*N_bins + i] > 0.0) {
			f[(counter*blockWidth+j)*N_bins + i] = log10(f[(counter*blockWidth+j)*N_bins + i] / (cellVolume*dU));
		     } else {
			f[(counter*blockWidth+j)*N_bins + i] = -100;
		     }
		  }
	       }
	       ++counter;
	    }
	    
	    for (size_t b=0; b<boundaryLocalIDs.size(); ++b) {
	       const pargrid::CellID blockLID = boundaryLocalIDs[b];
	       
	       int32_t blockWidth = 1;
	       if (lineoutCoordinate == 0) blockWidth = block::WIDTH_X;
	       if (lineoutCoordinate == 1) blockWidth = block::WIDTH_Y;
	       if (lineoutCoordinate == 2) blockWidth = block::WIDTH_Z;
	       
	       for (int j=0; j<blockWidth; ++j) {
		  const Real cellVolume = getCellVolume(blockLID,lineoutCoordinate,j);
		  const size_t N_bins = energyNodeCoordinatesOut.size()-1;
		  for (size_t i=0; i<N_bins; ++i) {
		     const Real U_min = pow(10.0,energyNodeCoordinatesOut[i  ]);
		     const Real U_max = pow(10.0,energyNodeCoordinatesOut[i+1]);
		     const Real dU = U_max-U_min;
		     if (f[(counter*blockWidth+j)*N_bins + i] > 0.0) {
			f[(counter*blockWidth+j)*N_bins + i] = log10(f[(counter*blockWidth+j)*N_bins + i] / (cellVolume*dU));
		     } else {
			f[(counter*blockWidth+j)*N_bins + i] = -100;
		     }
		  }  
	       }
	       ++counter;
	    }
	 
	    const string meshName = getLineoutName(line);
	 
	    map<string,string> xmlAttributes;
	    xmlAttributes["name"] = meshName + '/' + "EnergySpectrum/" + speciesNames[s];
	    xmlAttributes["mesh"] = meshName;
	    xmlAttributes["centering"] = "zone";
	    if (simClasses->vlsv.writeArray("VARIABLE",xmlAttributes,arraySize,1,&(f[0])) == false) {
	       simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Could not write energy spectrum on line #" << line << endl << write;
	       success = false;
	    }
	 }
      }

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif

      return success;
   }
   
   bool SpatialLineoutEnergySpectrum::writeLineoutMesh(size_t line,std::vector<pargrid::CellID>& innerLocalIDs,std::vector<pargrid::CellID>& boundaryLocalIDs) {
      bool success = true;
      
      // Get name of lineout mesh and lineout coordinate:
      const string meshName = getLineoutName(line);
      const uint8_t slicedCoordinate = getLineoutCoordinate(line);
      
      // Create mesh bounding box (significant only at master process):
      uint32_t bbox[6];
      bbox[1] = 0;
      bbox[0] = energyNodeCoordinates.size()-1; // Number of mesh blocks in y-direction in mesh bounding box.
      bbox[2] = 1;
      bbox[4] = 0;
      bbox[3] = 1;
      bbox[5] = 1;
      
      // Create array for node coordinate pointers (significant at master process only):
      Real* crds[3];
      Real zCoordinatesDummy[2];
      zCoordinatesDummy[0] = 0.0;
      zCoordinatesDummy[1] = 1.0;
      crds[1] = NULL;
      crds[0] = &(energyNodeCoordinatesOut[0]);
      crds[2] = zCoordinatesDummy;
      
      // Calculate number of nodes in each coordinate direction in 
      // lineout mesh (significant at master process only):
      Real N_nodes[3];
      N_nodes[1] = 0;
      N_nodes[0] = energyNodeCoordinatesOut.size();
      N_nodes[2] = 2;
      switch (slicedCoordinate) {
       case 0:
	 bbox[1] = sim->x_blocks;
	 bbox[4] = block::WIDTH_X;
	 crds[1] = sim->x_crds_node;
	 N_nodes[1] = sim->x_blocks*block::WIDTH_X + 1;
	 break;
       case 1:
	 bbox[1] = sim->y_blocks;
	 bbox[4] = block::WIDTH_Y;
	 crds[1] = sim->y_crds_node;
	 N_nodes[1] = sim->y_blocks*block::WIDTH_Y + 1;
	 break;
       case 2:
	 bbox[1] = sim->z_blocks;
	 bbox[4] = block::WIDTH_Z;
	 crds[1] = sim->z_crds_node;
	 N_nodes[1] = sim->z_blocks*block::WIDTH_Z + 1;
	 break;
       default:
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Invalid lineout coordinate '" << slicedCoordinate << "'" << endl << write;
	 return false;
	 break;
      }
      
      // Only master process writes bounding box and node coordinates. 
      // Set crds and N_nodes to dummy values on other processes:
      uint64_t arraySize = 6;
      uint32_t* ptr = bbox;
      if (sim->mpiRank != sim->MASTER_RANK) {
	 arraySize = 0;
	 ptr = NULL;
	 for (int i=0; i<3; ++i) crds[i]    = NULL;
	 for (int i=0; i<3; ++i) N_nodes[i] = 0;
      }

      // Get blocks that intersect the lineout on this process. This is done here
      // so that if problems come up we can exit before writing any data to file:
      if (getAcceptedBlocks(line,innerLocalIDs,boundaryLocalIDs) == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to get block local IDs" << endl << write;
	 return false;
      }

      // Write mesh bounding box and node coordinates (master process only):
      map<string,string> xmlAttributes;
      xmlAttributes["mesh"] = meshName;
      if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,arraySize,1,ptr) == false) success = false;
      if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttributes,N_nodes[0],1,crds[0]) == false) success = false;
      if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttributes,N_nodes[1],1,crds[1]) == false) success = false;
      if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttributes,N_nodes[2],1,crds[2]) == false) success = false;
      
      // Calculate new global IDs for cells on lineout mesh on this process:
      vector<pargrid::CellID> newGlobalIDs;
      for (size_t j=0; j<innerLocalIDs.size(); ++j) {
	 const pargrid::CellID blockLID = innerLocalIDs[j];
	 #ifndef NDEBUG
	    if (blockLID >= simClasses->pargrid.getNumberOfLocalCells()) {
	       cerr << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: P#" << sim->mpiRank << " invalid local ID #" << blockLID << endl;
	       exit(1);
	    }
         #endif
	 
	 uint32_t i_block,j_block,k_block,blockIndex,blockWidth;
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 
	 switch (slicedCoordinate) {
	  case 0:
	    blockIndex = i_block;
	    blockWidth = block::WIDTH_X;
	    break;
	  case 1:
	    blockIndex = j_block;
	    blockWidth = block::WIDTH_Y;
	    break;
	  case 2:
	    blockIndex = k_block;
	    blockWidth = block::WIDTH_Z;
	    break;
	  default:
	    blockIndex = pargrid::INVALID_CELLID;
	    blockWidth = 1;
	    break;
	 }

	 for (uint32_t k=0; k<blockWidth; ++k) for (size_t i=0; i<energyNodeCoordinates.size()-1; ++i) {
	    const pargrid::CellID newGID = getGlobalIndex(i,blockIndex*blockWidth+k);
	    newGlobalIDs.push_back(newGID);
	 }
      }
      
      for (size_t j=0; j<boundaryLocalIDs.size(); ++j) {
	 const pargrid::CellID blockLID = boundaryLocalIDs[j];
         #ifndef NDEBUG
	    if (blockLID >= simClasses->pargrid.getNumberOfLocalCells()) {
	       cerr << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: P#" << sim->mpiRank << " invalid local ID #" << blockLID << endl;
	       exit(1);
	    }
         #endif
	 
	 uint32_t i_block,j_block,k_block,blockIndex,blockWidth;
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 
	 switch (slicedCoordinate) {
	  case 0:
	    blockIndex = i_block;
	    blockWidth = block::WIDTH_X;
	    break;
	  case 1:
	    blockIndex = j_block;
	    blockWidth = block::WIDTH_Y;
	    break;
	  case 2:
	    blockIndex = k_block;
	    blockWidth = block::WIDTH_Z;
	    break;
	  default:
	    blockIndex = pargrid::INVALID_CELLID;
	    blockWidth = 1;
	    break;
	 }
	 
	 for (uint32_t k=0; k<blockWidth; ++k) for (size_t i=0; i<energyNodeCoordinates.size()-1; ++i) {
	    const pargrid::CellID newGID = getGlobalIndex(i,blockIndex*blockWidth+k);
	    newGlobalIDs.push_back(newGID);
	 }
      }

      // FIXME: N_blocks is local+ghosts:
      const uint64_t N_blocks = newGlobalIDs.size();
      const uint64_t N_ghosts = 0;
      // END FIXME

      // Write mesh global IDs:
      xmlAttributes.clear();
      xmlAttributes["name"] = meshName;
      xmlAttributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
      stringstream sss;
      sss << (pargrid::MPI_processID)sim->mpiProcesses;
      xmlAttributes["domains"] = sss.str();
      xmlAttributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
      xmlAttributes["xlabel"] = "Energy";
      xmlAttributes["xunits"] = "log10(eV/amu)";
      switch (slicedCoordinate) {
       case 0:
	 xmlAttributes["ylabel"] = "x-coordinate";
	 break;
       case 1:
	 xmlAttributes["ylabel"] = "y-coordinate";
	 break;
       case 2:
	 xmlAttributes["ylabel"] = "z-coordinate";
	 break;
      }
      
      if (simClasses->vlsv.writeArray("MESH",xmlAttributes,N_blocks,1,&(newGlobalIDs[0])) == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to write block global IDs" << endl << write;
	 success = false;
      }
      
      // Write domain size:
      xmlAttributes.clear();
      xmlAttributes["mesh"] = meshName;
      int domainSize[2];
      domainSize[0] = N_blocks;
      domainSize[1] = N_ghosts;
      if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,2,domainSize) == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to write domain size" << endl << write;
	 success = false;
      }

      vector<uint32_t> validNeighbors;
      vector<pargrid::MPI_processID> validHosts;
      
      // Write ghost block local IDs:
      if (simClasses->vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttributes,N_ghosts,1,&(validNeighbors[0])) == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to write ghost block local IDs!" << endl << write;
	 success = false;
      }
      
      // Write ghost block domain IDs:
      if (simClasses->vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttributes,N_ghosts,1,&(validHosts[0])) == false) {
	 simClasses->logger << "(SEP OP LINEOUT ENERGY SPECTRUM) ERROR: Failed to write block ghost domain IDs!" << endl << write;
	 success = false;
      }

      return success;
   }
   
} // namespace sep
