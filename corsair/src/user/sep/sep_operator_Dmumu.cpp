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

#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_operator_Dmumu.h"
#include "sep_particle_scatterer.h"
#include "sep_accumulation_stretched.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   
   OperatorDmumu::OperatorDmumu(): SpatialSliceOP() { 
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
      
      N_pitchCells = 20;
   }

   OperatorDmumu::~OperatorDmumu() { }

   bool OperatorDmumu::finalize() {
      bool success = true;
      return success;
   }

   std::string OperatorDmumu::getName() const {
      return "OperatorDmumu";
   }
   
   bool OperatorDmumu::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (SpatialSliceOP::initialize(cr,sim,simClasses) == false) return false;
      
      string prefix = "DifferentialFluxOperator";
      
      vector<string> instrumentNames;
      cr.addComposed(prefix+".instrument_names","Names of energy instruments (string).");
      cr.parse();
      cr.get(prefix+".instrument_names",instrumentNames);

      totalNumberOfArrays=0;
      for (size_t i=0; i<instrumentNames.size(); ++i) {
	 // Skip empty lines
	 if (instrumentNames[i].size() == 0) continue;
	 
	 // Read config file items for instrument and create it:
	 spacecraft::FluxInstrument dummyInstrument;
	 
	 if (dummyInstrument.createInstrument(simClasses,cr,prefix,instrumentNames[i],instruments,instrumentIndices) == false) {
	    simClasses.logger << "(SEP OP DIFF FLUX) ERROR: Failed to read required config file items for ";
	    simClasses.logger << "instrument '" << instrumentNames[i] << "'" << endl << write;
	    initialized = false;
	 } else {
	    totalNumberOfArrays += instruments[instruments.size()-1].minValues.size();
	 }
      }

      return success;
   }
   
   bool OperatorDmumu::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      if (baseClassInitialized == false) {
	 simClasses->logger << "(SEP OP DMUMU) ERROR: Base class is not initialized" << endl << write;
	 return false;
      }

      #if PROFILE_LEVEL > 0
         profile::start("WLM Slices",profTotalTime);
      #endif
      for (size_t slice=0; slice<getNumberOfSlices(); ++slice) {
	 // Tell spatial slice operator that we're processing new slice:
	 prepareSlice(slice,N_pitchCells);

	 // Create a name for the wavelength mesh slice:
	 stringstream ss;
	 ss << "PitchAngleMesh" << slice;
	 meshName = ss.str();

	 // Get sliced coordinate (0-2):
	 const uint8_t slicedCoordinate = getSlicedCoordinate(slice);

	 // Create mesh bounding box:
	 uint32_t bbox[6];
	 bbox[0] = sim->x_blocks;  // Number of mesh blocks in x-direction in mesh bounding box.
	 bbox[1] = sim->y_blocks;  // Number of mesh blocks in y-direction in mesh bounding box.
	 bbox[2] = sim->z_blocks;  // Number of mesh blocks in z-direction in mesh bounding box.
	 bbox[3] = block::WIDTH_X; // Number of cells in each block in x-direction.
	 bbox[4] = block::WIDTH_Y; // Number of cells in each block in y-direction.
	 bbox[5] = block::WIDTH_Z; // Number of cells in each block in z-direction.

	 if (isCylindrical(slice) == true) {
	    bbox[slicedCoordinate+0] = bbox[2];
	    bbox[slicedCoordinate+3] = bbox[5];
	    bbox[2] = N_pitchCells;
	    bbox[5] = 1;
	 } else {
	    bbox[slicedCoordinate+0] = N_pitchCells;
	    bbox[slicedCoordinate+3] = 1;
	 }

	 // Get number of nodes:
	 uint32_t N_nodes[3];
	 N_nodes[0] = sim->x_blocks*block::WIDTH_X + 1;
	 N_nodes[1] = sim->y_blocks*block::WIDTH_Y + 1;
	 N_nodes[2] = sim->z_blocks*block::WIDTH_Z + 1;

	 if (isCylindrical(slice) == true) {
	    N_nodes[slicedCoordinate] = N_nodes[2];
	    N_nodes[2] = N_pitchCells+1;
	 } else {
	    N_nodes[slicedCoordinate] = N_pitchCells+1;
	 }

	 // Create array containing wavelength mesh coordinates:
	 Real* wmesh = new Real[N_pitchCells+1];
	 const Real pitchMin = -1.0;
	 const Real pitchMax = +1.0;
	 const Real d_pitch = (pitchMax-pitchMin)/N_pitchCells;
	 for (size_t i=0; i<N_pitchCells+1; ++i) {
	    wmesh[i] = pitchMin + i*d_pitch;
	 }

	 // Get pointers to mesh node coordinates:
	 Real* crds[3];
	 crds[0] = sim->x_crds_node;
	 crds[1] = sim->y_crds_node;
	 crds[2] = sim->z_crds_node;
	 
	 if (isCylindrical(slice) == true) {
	    crds[slicedCoordinate] = crds[2];
	    crds[2] = wmesh;
	 } else {
	    crds[slicedCoordinate] = wmesh;
	 }

	 uint64_t arraySize = 6;
	 uint32_t* ptr = bbox;
	 if (sim->mpiRank != sim->MASTER_RANK) {
	    arraySize = 0;
	    ptr = NULL;
	    for (int i=0; i<3; ++i) crds[i]    = NULL;
	    for (int i=0; i<3; ++i) N_nodes[i] = 0;
	 }

	 vector<pargrid::CellID> validLocalIDs;
	 vector<pargrid::CellID> newGlobalIDs;
	 getAcceptedBlocks(slice,newGlobalIDs,validLocalIDs);

	 if (newGlobalIDs.size() == 0) {
	    simClasses->logger << "(SEP OP DMUMU) WARNING: Slice has no valid cells, slice: " << slice;
	    simClasses->logger << " sliced coordinate: " << (int)slicedCoordinate << endl;
	    simClasses->logger << "\t slice origin: " << getSliceOrigin(slice) << endl;
	 }

	 // Write mesh bounding box and node coordinates (master process only):
	 map<string,string> xmlAttributes;
	 xmlAttributes["mesh"] = meshName;
	 if (simClasses->vlsv.writeArray("MESH_BBOX",xmlAttributes,arraySize,1,ptr) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttributes,N_nodes[0],1,crds[0]) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttributes,N_nodes[1],1,crds[1]) == false) success = false;
	 if (simClasses->vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttributes,N_nodes[2],1,crds[2]) == false) success = false;

	 #warning IMPROVE ME Add ghost cells to wavelength mesh
	 const uint64_t N_blocks = newGlobalIDs.size(); // This is locals + ghosts
	 const uint64_t N_ghosts = 0;

	 // Write mesh global IDs:
	 xmlAttributes.clear();
	 xmlAttributes["name"] = meshName;
	 xmlAttributes["type"] = vlsv::mesh::STRING_UCD_MULTI;
	 stringstream sss;
	 sss << (pargrid::MPI_processID)sim->mpiProcesses;
	 xmlAttributes["domains"] = sss.str();

	 if (isCylindrical(slice) == true) {
	    xmlAttributes["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 } else {
	    xmlAttributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 }

	 if (simClasses->vlsv.writeArray("MESH",xmlAttributes,N_blocks,1,&(newGlobalIDs[0])) == false) {
	    simClasses->logger << "(OP DMUMU) ERROR: Failed to write block global IDs" << endl << write;
	    success = false;
	 }

	 // Write domain size:
	 xmlAttributes.clear();
	 xmlAttributes["mesh"] = meshName;
	 int domainSize[2];
	 domainSize[0] = N_blocks;
	 domainSize[1] = N_ghosts;
	 if (simClasses->vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttributes,1,2,domainSize) == false) {
	    simClasses->logger << "(OP DMUMU) ERROR: Failed to write domain size" << endl << write;
	    success = false;
	 }

	 vector<uint32_t> validNeighbors;
	 vector<pargrid::MPI_processID> validHosts;

	 // Write ghost block local IDs:
	 if (simClasses->vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttributes,N_ghosts,1,&(validNeighbors[0])) == false) {
	    simClasses->logger << "(OP DMUMU) ERROR: Failed to write ghost block local IDs!" << endl << write;
	    success = false;
	 }

	 // Write ghost block domain IDs:
	 if (simClasses->vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttributes,N_ghosts,1,&(validHosts[0])) == false) {
	    simClasses->logger << "(OP DMUMU) ERROR: Failed to write block ghost domain IDs!" << endl << write;
	    success = false;
	 }

	 delete [] wmesh; wmesh = NULL;
	 
	 for (uint32_t c=0; c<totalNumberOfArrays; ++c) {
	    currentInstrument = instrumentIndices[c].first;
	    currentChannel    = instrumentIndices[c].second;

	    if (writeDmumu(slice,currentInstrument,currentChannel,validLocalIDs,newGlobalIDs) == false) {
	       success = false;
	    }
	 }
      }
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
   bool OperatorDmumu::writeDmumu(int32_t slice,size_t currentInstrument,size_t currentChannel,
				  const std::vector<pargrid::CellID>& validLocalIDs,
				  const std::vector<pargrid::CellID>& newGlobalIDs) {
      bool success = true;

      // Calculate pitch angle diffusion coefficient at speed corresponding
      // to average channel energy:
      const Real avgEnergy = 0.5*(instruments[currentInstrument].minValues[currentChannel]
				  + instruments[currentInstrument].maxValues[currentChannel]);
      const Real channelSpeed = sqrt(2*avgEnergy/constants::MASS_PROTON);

      // Allocate arrays for wave intensity and calculated D_mumu:
      const int32_t NWL = simControl.N_wavelengthMeshCells+2;
      const int32_t WAVE_SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      Real* intensity = new Real[WAVE_SIZE];
      Real* buffer = new Real[validLocalIDs.size()*N_pitchCells];

      simClasses->pargrid.startNeighbourExchange(sim->defaultStencilID,simControl.parAlfvenWaveEnergyDataID);
      simClasses->pargrid.wait(sim->defaultStencilID,simControl.parAlfvenWaveEnergyDataID);

      const int alfvenSign = +1;
      size_t counter=0;
      for (size_t b=0; b<validLocalIDs.size(); ++b) {
	 // Load wave intensity:
	 const pargrid::CellID blockLID = validLocalIDs[b];	 
	 loadIntensity(blockLID,intensity,alfvenSign);

	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 uint32_t i_block,j_block,k_block;
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    // Calculate cell centroid in logical coordinates:
	    Real centroid[4];
	    centroid[0] = i_block*block::WIDTH_X + i + 0.5;
	    centroid[1] = j_block*block::WIDTH_Y + j + 0.5;
	    centroid[2] = k_block*block::WIDTH_Z + k + 0.5;
	    
	    // Get magnetic field magnitude:
	    Real B[3];
	    Real V_wave[3];
	    Real V_alfven,dV_wave;
	    (*simControl.fieldsGetState)(blockLID,sim->t,centroid,B,V_wave,dV_wave,V_alfven,alfvenSign);
	    const Real B_mag = vectorMagnitude<3>(B);
	    Real speed = channelSpeed - dotProduct<3>(V_wave,B)/B_mag;
	    speed = 1.4e7;
	    
	    // Calculate interpolation shape factors at cell centroid:
	    centroid[0] = i + 1.5;
	    centroid[1] = j + 1.5;
	    centroid[2] = k + 1.5;
	    centroid[3] = 1.5;

	    int32_t indices[4];
	    Real shapeFactors[12];
	    switch (simControl.order) {
	     case 0:
	       getShapeFactorsNGP_4D(centroid,indices,shapeFactors);
	       break;
	     case 1:
	       getShapeFactorsCIC_4D(centroid,indices,shapeFactors);
	       break;
	     case 2:
	       getShapeFactorsTSC_4D(centroid,indices,shapeFactors);
	       break;
	     default:
	       cerr << "ERROR: Unknown order of accuracy in sep_operator_mfp.cpp" << endl;
	       exit(1);
	       break;
	    }

	    // Calculate resonant wavelength:
	    const Real omega = constants::CHARGE_ELEMENTARY*B_mag/constants::MASS_PROTON;
	    Real maxResonantWavelength = 2*M_PI/omega * speed;

	    // Calculate diffusion constant:
	    const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);

	    // Calculate pitch angle diffusion coefficients:
	    const Real d_pitch = 2.0 / N_pitchCells;
	    for (uint32_t p=0; p<N_pitchCells; ++p) {
	       const Real pitch = -1.0 + (p+0.5)*d_pitch;

	       Real D_mumu,d_D_mumu;
	       //D_mumu = 1 + (*simControl.getLogicalWavelength)(pitch*maxResonantWavelength);
	       
	       evaluateDiffusionCoefficients(pitch,maxResonantWavelength,diffConstant,centroid,
					     indices,shapeFactors,intensity,D_mumu,d_D_mumu);
	
	       /*
	       centroid[3] = 1 + (*simControl.getLogicalWavelength)(pitch*maxResonantWavelength);
	       if (currentChannel == 0) {
		  getShapeFactorsNGP_4D(centroid,indices,shapeFactors);
		  D_mumu = interpolateScalarNGP_4D(intensity,indices,shapeFactors,NWL);
	       } else if (currentChannel == 1) {
		  getShapeFactorsCIC_4D(centroid,indices,shapeFactors);
		  D_mumu = interpolateScalarCIC_4D(intensity,indices,shapeFactors,NWL);
	       } else {
		  getShapeFactorsTSC_4D(centroid,indices,shapeFactors);
		  D_mumu = interpolateScalarTSC_4D(intensity,indices,shapeFactors,NWL);
	       }
		*/
	       buffer[counter] = D_mumu;
	       ++counter;
	    }
	 }
      }
      
      // Write buffer to output file:
      map<std::string,std::string> attribs;
      attribs["name"] = meshName+'/'+instruments[currentInstrument].name+'/'+instruments[currentInstrument].channelNames[currentChannel];
      attribs["mesh"] = meshName;
      attribs["unit"] = "1/s";
      attribs["centering"] = "zone";
      
      const uint64_t arraySize = validLocalIDs.size()*block::SIZE*N_pitchCells;
      const uint64_t vectorSize = 1;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,buffer) == false) success = false;

      // Deallocate memory and exit:
      delete [] intensity; intensity = NULL;
      delete [] buffer; buffer = NULL;
      return success;
   }
   
} // namespace sep
