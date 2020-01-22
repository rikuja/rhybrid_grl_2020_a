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
#include <map>
#include <vector>

#include <constants.h>
#include <ucd_mesh.h>
#include "sep_coordinate_transform.h"
#include "sep_simcontrol.h"
#include "sep_operator_shock_mesh.h"
#include "sep_operator_shock_variables.h"
#include "sep_base_class_shock.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;

   const string shockMeshName = "shock";
   
   OperatorShockVariables::OperatorShockVariables(): DataOperator() {
      #if PROFILE_LEVEL > 0
         profTotalTime = -1;
      #endif
   }
   
   OperatorShockVariables::~OperatorShockVariables() { }
   
   bool OperatorShockVariables::finalize() {
      bool success = true;
      return success;
   }
   
   std::string OperatorShockVariables::getName() const {
      return "ShockVariablesWriter";
   }
   
   bool OperatorShockVariables::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;
      if (DataOperator::initialize(cr,sim,simClasses) == false) success = false;
      return success;
   }
   
   bool OperatorShockVariables::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles) {
      bool success = true;
      if (simControl.includeShock == false) return success;
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP SHOCK VARIABLES OP) ERROR: Shock pointer is NULL while includeShock==true" << endl << write;
	 return false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::start("Shock Variables",profTotalTime);
      #endif

      // Allocate shock mesh:
      if (simControl.shock->initializeMesh(sim->t,simControl.N_shockSurfaceRefinements) == false) {
	 #if PROFILE_LEVEL > 0
	    profile::stop();
	 #endif
	 simClasses->logger << "(SEP SHOCK VARIABLES OP) ERROR: Shock mesh initialization failed" << endl << write;
	 return false;
      }

      vector<Real> nodeCoordsLogical;
      simControl.shock->getLogicalNodeCoordinates(nodeCoordsLogical);
      if (writeGasCompressionRatio(nodeCoordsLogical) == false) success = false;
      if (writeLogicalCoordinates(nodeCoordsLogical) == false) success = false;
      if (writeMachNumbers(nodeCoordsLogical) == false) success = false;
      if (writeMagneticField(nodeCoordsLogical) == false) success = false;
      if (writeShockNormals(nodeCoordsLogical) == false) success = false;
      if (writePlasmaVelocityHT(nodeCoordsLogical) == false) success = false;
      if (writePlasmaVelocityShockFrame(nodeCoordsLogical) == false) success = false;
      if (writePlasmaVelocitySNIF(nodeCoordsLogical) == false) success = false;
      if (writeVariableMagneticCompressionRatio(nodeCoordsLogical) == false) success = false;
      if (writeShockPotential(nodeCoordsLogical) == false) success = false;
      if (writeShockVelocity(nodeCoordsLogical) == false) success = false;
      
      // Deallocate shock mesh:
      if (simControl.shock->finalizeMesh() == false) {
	 simClasses->logger << "(SEP SHOCK VARIABLES OP) ERROR: Failed to finalize shock mesh" << endl << write;
	 success = false;
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

   bool OperatorShockVariables::writeGasCompressionRatio(const std::vector<Real>& nodeCoords) {
      bool success = true;
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP OP SHOCK VARS) ERROR: Shock class is NULL" << endl << write;
	 return false;
      }
      if (simControl.shock->isInitialized() == false) {
	 simClasses->logger << "(SEP OP SHOCK VARS) ERROR: Shock class is NULL" << endl << write;
	 return false;
      }

      #if PROFILE_LEVEL > 0
         static int profCounter=-1;
         profile::start("Gas Compr Ratio",profCounter);
      #endif
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> gasCompr(N_nodes);
      vector<Real> magnCompr(N_nodes);
      vector<Real> mirrorRatio(N_nodes);

      // Only master process writes data to output file, all other processes
      // must still call writeVariable below to prevent deadlock:
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;

      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    gasCompr[node] = 1.0;
	    magnCompr[node] = 1.0;
	    mirrorRatio[node] = 1.0;
	    continue;
	 }
	 
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);

	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim->t,&(nodeCoords[3*node+0]),E,B,gradB);
	 Real shockNormal[3];
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);
	 const Real B1_norm = -dotProduct<3>(B,shockNormal);
	 const Real B1_tang = sqrt(vectorMagnitude2<3>(B) - B1_norm*B1_norm);
	 
	 gasCompr[node] = simControl.shock->getGasCompressionRatio(blockLID,sim->t,&(nodeCoords[3*node+0]));
	 magnCompr[node] = simControl.shock->getMagneticCompressionRatio(blockLID,sim->t,&(nodeCoords[3*node+0]));
	 const Real B2_tang2 = magnCompr[node]*magnCompr[node]*B1_tang*B1_tang;
	 mirrorRatio[node] = sqrt(B1_norm*B1_norm+B2_tang2)/sqrt(B1_norm*B1_norm+B1_tang*B1_tang);
      }

exit:      
      if (shockmesh::writeVariable(sim,simClasses,true,gasCompr,1,"gas_compression") == false) success = false;
      if (shockmesh::writeVariable(sim,simClasses,true,magnCompr,1,"magn_compression") == false) success = false;
      if (shockmesh::writeVariable(sim,simClasses,true,mirrorRatio,1,"B2_per_B1") == false) success = false;

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
   bool OperatorShockVariables::writeLogicalCoordinates(const std::vector<Real>& nodeCoords) {
      return shockmesh::writeVariable(sim,simClasses,true,nodeCoords,3,"logical_coords");
   }
   
   bool OperatorShockVariables::writeMachNumbers(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         static int profCounter = -1;
         profile::start("Mach numbers",profCounter);
      #endif
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> alfvenMach(N_nodes);
      vector<Real> sonicMach(N_nodes);
      vector<Real> fastMach(N_nodes);

      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    alfvenMach[node] = 0.0;
	    sonicMach[node] = 0.0;
	    continue;
	 }
	 
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);
	 
	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim->t,&(nodeCoords[3*node+0]),E,B,gradB);
	 
	 Real shockNormal[3];
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);
	 
	 const Real B_norm_mag = dotProduct<3>(B,shockNormal);
	 Real B_tang[3];
	 for (int i=0; i<3; ++i) B_tang[i] = B[i] - B_norm_mag*shockNormal[i];
	 const Real psi_1 = atan(vectorMagnitude<3>(B_tang) / fabs(B_norm_mag));
	 
	 Real V_shock_SIM[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock_SIM);
	 
	 Real V_plasma_SRF[3];
	 for (int i=0; i<3; ++i) V_plasma_SRF[i] = plasmaState.V_plasma_SIM[i] - V_shock_SIM[i];
	 Real V_plasma_norm_mag2 = dotProduct<3>(V_plasma_SRF,shockNormal);
	 V_plasma_norm_mag2 = V_plasma_norm_mag2*V_plasma_norm_mag2;
	 
	 alfvenMach[node] = sqrt(V_plasma_norm_mag2 / plasmaState.alfvenSpeed2);
	 sonicMach[node] = sqrt(V_plasma_norm_mag2 / plasmaState.soundSpeed2);

	 Real fastSpeed2 = plasmaState.alfvenSpeed2 + plasmaState.soundSpeed2;
	 fastSpeed2 = fastSpeed2 + sqrt(fastSpeed2*fastSpeed2 - 4*plasmaState.alfvenSpeed2*plasmaState.soundSpeed2*cos(psi_1)*cos(psi_1));
	 fastSpeed2 = 0.5 * fastSpeed2;
	 fastMach[node] = sqrt(V_plasma_norm_mag2 / fastSpeed2);
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,alfvenMach,1,"Mach_Alfvenic") == false) success = false;
      if (shockmesh::writeVariable(sim,simClasses,true,sonicMach,1,"Mach_sonic") == false) success = false;
      if (shockmesh::writeVariable(sim,simClasses,true,fastMach,1,"Mach_fastMS") == false) success = false;
      
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      
      return success;
   }
   
   bool OperatorShockVariables::writeMagneticField(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> values(N_nodes*3);

      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    for (int j=0; j<3; ++j) values[3*node+j] = 0.0;
	    continue;
	 }

	 // Calculate block global ID:
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    for (int j=0; j<3; ++j) values[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Get fields:
	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim->t,&(nodeCoords[3*node+0]),E,B,gradB);
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),B,&(values[3*node+0]));
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,values,3,"B_shock") == false) success = false;
      return success;
   }

   bool OperatorShockVariables::writePlasmaVelocityHT(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> V_plasma(N_nodes*3);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Calculate block global ID:
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Get fields:
	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim->t,&(nodeCoords[3*node+0]),E,B,gradB);
	 
	 // Get plasma state at particle position:
	 if (simControl.fieldsGetPlasmaState == NULL) {
	    std::cerr << "getPlasmaState is NULL" << std::endl;
	    exit(1);
	 }
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);
      
	 // Get shock velocity and normal:
	 Real V_shock[3];
	 Real shockNormal[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock);
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);

	 // Get transformation velocity:
	 Real V_transform[3];
	 getTransformSimToLocalHT(plasmaState.V_plasma_SIM,B,V_shock,shockNormal,V_transform);

	 Real V_plasma_HT[3];
	 for (int i=0; i<3; ++i) V_plasma_HT[i] = plasmaState.V_plasma_SIM[i] - V_transform[i];
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),V_plasma_HT,&(V_plasma[3*node+0]));
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,V_plasma,3,"V_plasma_HT") == false) success = false;
      return success;
   }
   
   bool OperatorShockVariables::writePlasmaVelocityShockFrame(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> V_plasma(N_nodes*3);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Calculate block global ID:
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Get plasma state at particle position:
	 if (simControl.fieldsGetPlasmaState == NULL) {
	    std::cerr << "getPlasmaState is NULL" << std::endl;
	    exit(1);
	 }
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);
	 
	 // Get local shock velocity in simulation frame:
	 Real V_shock_SIM[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock_SIM);
	 
	 Real V_plasma_SRF[3];
	 for (int i=0; i<3; ++i) V_plasma_SRF[i] = plasmaState.V_plasma_SIM[i] - V_shock_SIM[i];
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),V_plasma_SRF,&(V_plasma[3*node+0]));
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,V_plasma,3,"V_plasma_SRF") == false) success = false;
      return success;
   }
   
   bool OperatorShockVariables::writePlasmaVelocitySNIF(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> V_plasma(N_nodes*3);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	       
	 // Calculate block global ID:
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    for (int j=0; j<3; ++j) V_plasma[3*node+j] = 0.0;
	    continue;
	 }
	 
	 // Get plasma state at particle position:
	 if (simControl.fieldsGetPlasmaState == NULL) {
	    std::cerr << "getPlasmaState is NULL" << std::endl;
	    exit(1);
	 }
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);

	 // Get local shock speed & normal:
	 Real V_shock[3];
	 Real shockNormal[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock);
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);
	 
	 // Get transform velocity:
	 Real V_transform[3];
	 getTransformSimToLocalSNIF(plasmaState.V_plasma_SIM,V_shock,shockNormal,V_transform);

	 Real V_plasma_SNIF[3];
	 for (int i=0; i<3; ++i) V_plasma_SNIF[i] = plasmaState.V_plasma_SIM[i] - V_transform[i];
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),V_plasma_SNIF,&(V_plasma[3*node+0]));
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,V_plasma,3,"V_plasma_SNIF") == false) success = false;
      return success;
   }
	 
   bool OperatorShockVariables::writeShockNormals(const std::vector<Real>& nodeCoords) {
      bool success = true;
      /*
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> normals(N_nodes*3);
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    normals[3*node+0] = 0.0;
	    normals[3*node+1] = 0.0;
	    normals[3*node+2] = 0.0;
	    continue;
	 }

	 // Get shock normal:
	 Real normal[3];
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),normal);
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),normal,&(normals[3*node+0]));
      }

      if (shockmesh::writeVariable(sim,simClasses,true,normals,3,"shock_normal") == false) success = false;
      */
      uint64_t N_surfaces;
      const vector<Real>& faceData = simControl.shock->getFaceData(N_surfaces);
      vector<Real> surfaceAreas(N_surfaces);
      vector<Real> normals(N_surfaces*3);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (uint64_t i=0; i<N_surfaces; ++i) {
	 surfaceAreas[i] = faceData[i*ucdmesh::facedataelement::SIZE + ucdmesh::facedataelement::AREA];
	 normals[i*3+0]  = faceData[i*ucdmesh::facedataelement::SIZE + ucdmesh::facedataelement::NORMAL_X];
	 normals[i*3+1]  = faceData[i*ucdmesh::facedataelement::SIZE + ucdmesh::facedataelement::NORMAL_Y];
	 normals[i*3+2]  = faceData[i*ucdmesh::facedataelement::SIZE + ucdmesh::facedataelement::NORMAL_Z];
      }
      
exit:
      if (shockmesh::writeVariable(sim,simClasses,false,surfaceAreas,1,"surface_area") == false) success = false;
      if (shockmesh::writeVariable(sim,simClasses,false,normals,3,"surface_normal") == false) success = false;
      return success;
   }
   
   bool OperatorShockVariables::writeShockPotential(const std::vector<Real>& nodeCoords) {
      bool success = true;
      #if PROFILE_LEVEL > 0
         static int profID = -1;
         profile::start("Cross-shock potential",profID);
      #endif

      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> potentials(N_nodes);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    potentials[node] = 0.0;
	    continue;
	 }
	 
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    potentials[node] = 0.0;
	    continue;
	 }
	 
	 // Get gas compression ratio:
	 const Real gasCompressionRatio = simControl.shock->getGasCompressionRatio(blockLID,sim->t,&(nodeCoords[3*node+0]));

	 // Get plasma state:
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);

	 // Get local shock velocity in simulation frame:
	 Real V_shock_sim[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock_sim);

	 // Get shock normal:
	 Real shockNormal[3];
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);

	 // Calculate normal component of plasma velocity in local SRF
	 Real V_plasma_SRF[3];
	 for (int i=0; i<3; ++i) V_plasma_SRF[i] = plasmaState.V_plasma_SIM[i] - V_shock_sim[i];
	 const Real V_plasma_SRF_norm = dotProduct<3>(V_plasma_SRF,shockNormal);
	 
	 // Calculate cross-shock potential:
	 potentials[node] = 0.5*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*V_plasma_SRF_norm*V_plasma_SRF_norm;
	 potentials[node] *= (1 - 1/gasCompressionRatio);
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,potentials,1,"node/potential") == false) success = false;
      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }

   bool OperatorShockVariables::writeShockVelocity(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> values(N_nodes*3);
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    for (int i=0; i<3; ++i) values[3*node+i] = 0.0;	    
	    continue;
	 }
	 
	 // Get local shock velocity in simulation frame:
	 Real V_shock[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock);
	 shockmesh::transformToCartesian(sim,&(nodeCoords[3*node+0]),V_shock,&(values[3*node+0]));
      }

exit:
      if (shockmesh::writeVariable(sim,simClasses,true,values,3,"V_shock_SIM") == false) success = false;
      return success;
   }
   
   bool OperatorShockVariables::writeVariableMagneticCompressionRatio(const std::vector<Real>& nodeCoords) {
      bool success = true;
      
      const size_t N_nodes = nodeCoords.size()/3;
      vector<Real> lossConeCosine(N_nodes);
      vector<Real> V_HT(N_nodes*3);
      vector<Real> V_dot_B(N_nodes);
      const Real DUMMY_COSINE = 0.0;
      
      if (sim->mpiRank != sim->MASTER_RANK) goto exit;
      for (size_t node=0; node<N_nodes; ++node) {
	 // Check that node is inside simulation domain. If not, 
	 // write unit value:
	 if (nodeCoords[3*node+0] == numeric_limits<Real>::infinity()) {
	    lossConeCosine[node] = DUMMY_COSINE;
	    V_dot_B[node] = 0.0;
	    continue;
	 }

	 // Calculate block global ID:
	 pargrid::CellID blockLID,blockGID;
	 shockmesh::getBlockIDs(sim,simClasses,&(nodeCoords[3*node+0]),blockLID,blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    lossConeCosine[node] = DUMMY_COSINE;
	    V_dot_B[node] = 0.0;
	    continue;
	 }
	 
	 // Get magnetic field at particle position:
	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim->t,&(nodeCoords[3*node+0]),E,B,gradB);

	 // Get plasma state at particle position:
	 if (simControl.fieldsGetPlasmaState == NULL) {
	    std::cerr << "getPlasmaState is NULL" << std::endl;
	    exit(1);
	 }
	 PlasmaState plasmaState;
	 (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,&(nodeCoords[3*node+0]),plasmaState);

	 // Get local shock velocity in simulation frame:
	 Real V_shock_sim[3];
	 simControl.shock->getLocalShockVelocity(sim->t,&(nodeCoords[3*node+0]),V_shock_sim);
	 
	 // Get shock normal:
	 Real shockNormal[3];
	 simControl.shock->getShockNormal(sim->t,&(nodeCoords[3*node+0]),shockNormal);
	 
	 // Get magnetic compression ratio:
	 const Real R_magn = simControl.shock->getMagneticCompressionRatio(blockLID,sim->t,&(nodeCoords[3*node+0]));
	 
	 // Calculate B2/B1 ratio from R_magn
	 const Real B1_norm_mag = dotProduct<3>(B,shockNormal);
	 const Real B1_mag2 = vectorMagnitude2<3>(B);
	 const Real B1_tang_mag2 = B1_mag2 - B1_norm_mag*B1_norm_mag;
	 const Real B2_mag2 = B1_norm_mag*B1_norm_mag + R_magn*R_magn*B1_tang_mag2;
	 const Real B2_per_B1 = sqrt(B2_mag2 / B1_mag2);
	 
	 // Get transformation from simulation frame to local HT frame:
	 getTransformSimToLocalHT(plasmaState.V_plasma_SIM,B,V_shock_sim,shockNormal,&(V_HT[3*node+0]));

	 // Get loss cone cosine:
	 lossConeCosine[node] = getLossConeCosine(B,shockNormal,B2_per_B1);
	 
	 // Calculate dot product between plasma velocity in HT frame and B:
	 Real V_plasma_HT[3];
	 for (int i=0; i<3; ++i) V_plasma_HT[i] = plasmaState.V_plasma_SIM[i] - V_HT[3*node+i];
	 V_dot_B[node] = fabs(dotProduct<3>(V_plasma_HT,B) / vectorMagnitude<3>(V_plasma_HT) / vectorMagnitude<3>(B));
      }

exit:
      // Write loss cone cosine:
      if (shockmesh::writeVariable(sim,simClasses,true,lossConeCosine,1,"loss_cone_cosine") == false) success = false;
      
      // Write HT transformation vector:
      if (shockmesh::writeVariable(sim,simClasses,true,V_HT,3,"V_dHT") == false) success = false;

      if (shockmesh::writeVariable(sim,simClasses,true,V_dot_B,1,"V_plasma_HT_dot_B") == false) success = false;

      return success;
   }
   
} // namespace sep
