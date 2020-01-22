/** This file is part of Corsair simulation.
 *
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

#include <main.h>
#include <user.h>
#include <gridbuilder.h>
#include <logically_cartesian_builder.h>
#include <dataoperatorcontainer.h>
#include <operator_cellid.h>
#include <operator_mpirank.h>
#include <operator_load.h>
#include <operator_pargrid_array.h>

// Include a struct containing simulation control variables:
#include "example_advection.h"
#include "advection_propagator.h"

using namespace std;

int index(int i,int j,int k) {return k*block::WIDTH_Y*block::WIDTH_X+j*block::WIDTH_X+i;}

/** Propagator function. This function is called once per time step. User should call 
 * here functions that propagate the simulation forward from t=sim.t to t=sim.t+sim.dt.
 * @param sim Struct containing various simulation parameters.
 * @param simClasses Struct containing various simulation classes.
 * @param cr Configugation file reader.
 * @return If true, simulation was propagated successfully.*/
bool propagate(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   bool rvalue = true;
   if (advectionPropagate(sim,simClasses) == false) rvalue = false;
   return rvalue;
}

/** Early initialization function. This function is mainly meant to be used to initialize 
 * dynamically allocated cell data before parallel grid is initialized. For most purposes 
 * userLateInitialization should be used instead.
 * @param sim Struct containing various simulation parameters.
 * @param simClasses Struct containing various simulation classes.
 * @param cr Configugation file reader.
 * @param particleLists Particle lists (if allocated).
 * @return If true, initialization completed successfully.*/
bool userEarlyInitialization(Simulation& sim,SimulationClasses& simClasses,
			     ConfigReader& cr,vector<ParticleListBase*>& particleLists) {return true;}

/** Late initialization function. This function is called after everything else has 
 * been initialized. This function is meant to be used to initialize simulation data etc.
 * @param sim Struct containing various simulation parameters.
 * @param simClasses Struct containing various simulation classes.
 * @param cr Configugation file reader.
 * @param particleLists Particle lists (if allocated).
 * @return If true, initialization completed successfully.*/
bool userLateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			    const ObjectFactories& objectFactories,vector<ParticleListBase*>& particleLists) {
   // ***** ALLOCATE USER DATA IN PARGRID ***** //
   
   // Init struct Advection contents. For the moment transfer IDs 
   // (Advection::avgsTransferID and Advection::fluxTransferID below) must 
   // be guessed:
   Advection::faceNbrStencilID = simClasses.pargrid.invalidStencilID();
   Advection::dataAvgsID       = simClasses.pargrid.invalidDataID();
   Advection::dataFluxID       = simClasses.pargrid.invalidDataID();
   Advection::avgsName         = "avgs";
   Advection::fluxName         = "flux";
   
   // Create a new transfer stencil. Here we tell ParGrid which neighbours are to be 
   // received, i.e. define the neighbour data that is needed to propagate a cell.
   // Send stencil is calculated automatically. Here we only need data from face neighbours.
   vector<pargrid::NeighbourID> nbrTypeIDs;
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0));
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0));
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0));
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0));
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1));
   nbrTypeIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1));
   Advection::faceNbrStencilID = simClasses.pargrid.addStencil(pargrid::localToRemoteUpdates,nbrTypeIDs);
   if (Advection::faceNbrStencilID == simClasses.pargrid.invalidStencilID()) {
      simClasses.logger << "(USER) ERROR: Failed to add face neighbour stencil to ParGrid!" << endl << write;
      return false;
   }

   // Create a parallel data array for function values. Since a patch/block of cells is stored 
   // in each parallel cell, we need (size of block) x (one double) of data per parallel cell:
   Advection::dataAvgsID = simClasses.pargrid.addUserData<Real>(Advection::avgsName,block::SIZE);
   if (Advection::dataAvgsID == simClasses.pargrid.invalidCellID()) {
      simClasses.logger << "(USER) ERROR: Failed to add avgs array to ParGrid!" << endl << write;
      return false;
   }
   
   // Create a parallel data array for fluxes:
   Advection::dataFluxID = simClasses.pargrid.addUserData<Real>(Advection::fluxName,block::SIZE);
   if (Advection::dataFluxID == simClasses.pargrid.invalidCellID()) {
      simClasses.logger << "(USER) ERROR: Failed to add flux array to ParGrid!" << endl << write;
      return false;
   }

   // Tell ParGrid that function values should be transferred using the stencil created above.
   if (simClasses.pargrid.addDataTransfer(Advection::dataAvgsID,Advection::faceNbrStencilID) == false) {
      simClasses.logger << "(USER) ERROR: Failed to add avgs data transfer!" << endl << write;
      return false;
   }
   
   // Tell ParGrid that fluxes should be transferred using the stencil created above.
   if (simClasses.pargrid.addDataTransfer(Advection::dataFluxID,Advection::faceNbrStencilID) == false) {
      simClasses.logger << "(USER) ERROR: Failed to add flux data transfer!" << endl << write;
      return false;
   }
   
   // ***** EXAMPLE OF HOW TO USE CONFIG FILE READER ***** //
   
   // Read advection velocity from config file:
   const Real defaultValue = 0.0;
   cr.add("Advection.velocity_x","Advection velocity x-component in m/s (float)",defaultValue);
   cr.add("Advection.velocity_y","Advection velocity y-component in m/s (float)",defaultValue);
   cr.add("Advection.velocity_z","Advection velocity z-component in m/s (float)",defaultValue);
   cr.parse();
   cr.get("Advection.velocity_x",Advection::Vx);
   cr.get("Advection.velocity_y",Advection::Vy);
   cr.get("Advection.velocity_z",Advection::Vz);

   // ***** INITIAL STATE ***** //
   
   // If simulation was restarted then ParGrid arrays already contain correct data:
   if (sim.restarted == true) return true;
   
   const Real x_max = 30.0e6;
   const Real y_max = 30.0e6;
   const Real z_max = 30.0e6;
   Real* avgs = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Advection::dataAvgsID));
   
   // Iterate over all blocks local to this process:
   const double* blockCrds = getBlockCoordinateArray(sim,simClasses);
   for (size_t block=0; block<simClasses.pargrid.getNumberOfLocalCells(); ++block) {
      Real cellSizes[3];
      getBlockCellSize(simClasses,sim,block,cellSizes);
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	 const Real x = blockCrds[3*block+0] + (i+0.5)*cellSizes[0];
	 const Real y = blockCrds[3*block+1] + (j+0.5)*cellSizes[1];
	 const Real z = blockCrds[3*block+2] + (k+0.5)*cellSizes[2];
	 if (fabs(x) <= x_max && (fabs(y) <= y_max && fabs(z) <= z_max)) {
	    avgs[block*block::SIZE+index(i,j,k)] = 1.0;
	 } else {
	    avgs[block*block::SIZE+index(i,j,k)] = 0.0;
	 }
      }
   }
   
   return true;
}

/** Finalization function that should deallocate or finalize memory that 
 * has been allocated by userEarlyInitialization or userLateInitialization functions.
 * @param sim Struct containing various simulation parameters.
 * @param simClasses Struct containing various simulation classes.
 * @param particleLists Particle lists (if allocated).
 * @return If true, finalization completed successfully.*/
bool userFinalization(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   bool success = true;
   if (simClasses.pargrid.removeUserData(Advection::dataAvgsID) == false) success = false;
   if (simClasses.pargrid.removeUserData(Advection::dataFluxID) == false) success = false;
   return success;
}

/** Add all known or needed GridBuilders to a special container (GridBuilderFactory). The 
 * name of the GridBuilder that is used to construct the simulation mesh is given in configuration file.
 * The GridBuilder name given in config file must match one of the names passed to GridBuilderFactory here.
 * @param gbf Container for GridBuilders.
 * @return If true, all GridBuilders were registered successfully.*/
bool registerObjectMakers(ObjectFactories& objectFactories) {
   bool success = true;
   
   if (objectFactories.gridBuilders.registerMaker("LogicallyCartesian",LCCreator) == false) success = false;
   return success;
}

/** Add all known or needed DataOperators to a special container (DataOperatorContainer). 
 * When simulation data is saved, all registered DataOperators all requested to save their 
 * data to output file. In other words, only data written by DataOperators registered here 
 * are written to output file(s).
 * @param doc Container for DataOperators.
 * @return If true, DataOperators were registered successfully.*/
bool registerDataOperators() {
   bool success = true;
   
   DataOperatorContainer& doc = corsair::getObjectWrapper().dataOperatorContainer;
   if (doc.registerOperator(new MPIRank) == false) success = false;
   if (doc.registerOperator(new CellIDOP) == false) success = false;
   if (doc.registerOperator(new LoadOP) == false) success = false;
   if (doc.registerOperator(new OperatorPargridArray) == false) success = false;
   return success;
}

bool userRunTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists) {return true;}
