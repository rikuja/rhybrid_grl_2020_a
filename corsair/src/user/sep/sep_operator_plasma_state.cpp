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

#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_fields_container.h"
#include "sep_operator_plasma_state.h"
#include "sep_coordinate_transform.h"

using namespace std;

namespace sep {
   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   OperatorPlasmaState::OperatorPlasmaState(): DataOperator() {
      initialized = false;
      
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorPlasmaState::~OperatorPlasmaState() {
      finalize();
   }
   
   bool OperatorPlasmaState::finalize() {
      if (initialized == false) return true;
      initialized = false;
      return true;
   }

   std::string OperatorPlasmaState::getName() const {return "sepOperatorPlasmaState";}
   
   bool OperatorPlasmaState::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      initialized = true;
      
      if (DataOperator::initialize(cr,sim,simClasses) == false) {
	 simClasses.logger << "(SEP OP PLASMA STATE) ERROR: DataOperator failed to initialize!" << endl << write;
	 initialized = false;
      }

      return initialized;
   }
   
   bool OperatorPlasmaState::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
      if (getInitialized() == false) return false;
      bool success = true;

      #if PROFILE_LEVEL > 0
         profile::start("Plasma Variables",profTotalTime);
      #endif
   
      // Get local cell IDs:
      const pargrid::CellID N_localBlocks = simClasses->pargrid.getNumberOfLocalCells();

      // Allocate arrays for E,B vectors:
      Real* B = new Real[3*N_localBlocks*block::SIZE];
      Real* E = new Real[3*N_localBlocks*block::SIZE];
      Real* V = new Real[3*N_localBlocks*block::SIZE];
      Real* V_alfven_par  = new Real[3*N_localBlocks*block::SIZE];
      Real* V_alfven_antipar = new Real[3*N_localBlocks*block::SIZE];
      Real* massDensity = new Real[N_localBlocks*block::SIZE];
      
      // Calculate E & B vectors at local cell centroids:
      size_t index = 0;
      int32_t i_block,j_block,k_block;
      int32_t i_cell,j_cell,k_cell;
      Real centroid[3];
      Real position[3];
      for (pargrid::CellID blockLID=0; blockLID<N_localBlocks; ++blockLID) {
	 // Get block global ID:
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 
	 // Calculate block (i,j,k) indices:
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	    // Calculate cell (i,j,k) indices:
	    i_cell = i_block*block::WIDTH_X + i;
	    j_cell = j_block*block::WIDTH_Y + j;
	    k_cell = k_block*block::WIDTH_Z + k;
	    
	    // Calculate cell centroid in logical coordinates:
	    centroid[0] = i_cell + 0.5;
	    centroid[1] = j_cell + 0.5;
	    centroid[2] = k_cell + 0.5;
	    
	    // Get E,B fields at cell centroid and copy to buffers:
	    //(*simControl.fieldsGetFields)(blockLID,centroid,E_tmp,B_tmp,gradB_tmp);
	    
	    // Get plasma state at cell centroid:
	    PlasmaState plasmaState;
	    (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,centroid,plasmaState);
	    massDensity[index/3] = plasmaState.ionMassDensity;

	    Real V_wave_par[3];
	    Real V_wave_antipar[3];
	    Real dummy[3];
	    Real dV_wave;
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim->t,centroid,dummy,V_wave_par,dV_wave,V_alfven,+1);
	    (*simControl.fieldsGetState)(blockLID,sim->t,centroid,dummy,V_wave_antipar,dV_wave,V_alfven,-1);
	    
	    // Calculate centroid in physical coordinates:
	    position[0] = sim->x_crds_node[i_cell] + 0.5*sim->dx_cell[i_cell];
	    position[1] = sim->y_crds_node[j_cell] + 0.5*sim->dy_cell[j_cell];
	    position[2] = sim->z_crds_node[k_cell] + 0.5*sim->dz_cell[k_cell];
	    
	    // Get vector transformation matrix to Cartesian coordinates:
	    Real matrix[9];
	    getTransformMatrix(position,matrix,simControl.coordinateSystem,sep::CARTESIAN);
	    matrixVectorMultiply<3>(matrix,plasmaState.E,&(E[index+0]));
	    matrixVectorMultiply<3>(matrix,plasmaState.B,&(B[index+0]));
	    matrixVectorMultiply<3>(matrix,plasmaState.V_plasma_SIM,&(V[index+0]));
	    matrixVectorMultiply<3>(matrix,V_wave_par,&(V_alfven_par[index+0]));
	    matrixVectorMultiply<3>(matrix,V_wave_antipar,&(V_alfven_antipar[index+0]));

	    index += 3;
	 }
      }
	    
      // Write arrays to output file:
      std::map<std::string,std::string> attribs;
      attribs["name"] = spatMeshName + "/E";
      attribs["mesh"] = spatMeshName;
      attribs["unit"] = "V/m";
      attribs["centering"] = "zone";

      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
       case sep::CARTESIAN:
	 attribs["geometry"] = vlsv::geometry::STRING_CARTESIAN;
	 break;
       case sep::CYLINDRICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_CYLINDRICAL;
	 break;
       case sep::SPHERICAL:
	 attribs["geometry"] = vlsv::geometry::STRING_SPHERICAL;
	 break;
       default:
	 attribs["geometry"] = vlsv::geometry::STRING_UNKNOWN;
	 break;
      }      
      
      const uint64_t arraySize = N_localBlocks*block::SIZE;
      const uint64_t vectorSize = 3;
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,E) == false) success = false;
      
      attribs["name"] = spatMeshName + "/B";
      attribs["unit"] = "T";
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,B) == false) success = false;
      
      attribs["name"] = spatMeshName + "/V_plasma";
      attribs["unit"] = "m/s";
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,V) == false) success = false;
      
      attribs["name"] = spatMeshName + "/mass_density";
      attribs["unit"] = "kg/m3";
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,massDensity) == false) success = false;
      
      attribs["name"] = spatMeshName + "/V_alfven_antipar";
      attribs["unit"] = "m/s";
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,V_alfven_antipar) == false) success = false;
      
      attribs["name"] = spatMeshName + "/V_alfven_par";
      attribs["unit"] = "m/s";
      if (simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorSize,V_alfven_par) == false) success = false;
            
      // Deallocate memory and exit:
      delete [] E; E = NULL;
      delete [] B; B = NULL;
      delete [] V; V = NULL;
      delete [] massDensity; massDensity = NULL;
      delete [] V_alfven_antipar; V_alfven_antipar = NULL;
      delete [] V_alfven_par; V_alfven_par = NULL;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
      
   
} // namespace sep
