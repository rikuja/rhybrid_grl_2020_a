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
#include <cmath>
#include <omp.h>

#include "example_advection.h"
#include "advection_propagator.h"

using namespace std;

// Variables used in profiling
static int totalID = -1;
static int fluxesID = -1;
static int propagID = -1;
static int mpiWaitID = -1;

// Forward declarations 
void fetchData(Real* data,Real* array,SimulationClasses& simClasses,pargrid::CellID blockID);
void saveData(Real* data,Real* array,pargrid::CellID blockID);

/** Propagate cell values forward in time using advection equation. 
 * This example also demonstrates how to use profiling.
 * @param sim Generic simulation control variables.
 * @param simClasses General-user simulation classes.
 * @return If true, cell values were propagated successfully.*/
bool advectionPropagate(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   profile::start("Advection Eqn",totalID);

   // Get pointer to data and flux arrays:
   Real* data = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Advection::dataAvgsID));
   Real* flux = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Advection::dataFluxID));   
   if (data == NULL) {cerr << "ERROR: obtained NULL data array!" << endl; exit(1);}
   if (flux == NULL) {cerr << "ERROR: obtained NULL flux array!" << endl; exit(1);}

   int tid; /**< Thread ID (private).*/

   // Start of parallel construct, creates threads:
   #pragma omp parallel private(tid)
     {
	// Get thread ID:
	tid = omp_get_thread_num();

	// First task is to set boundary conditions. This needs to be done before starting data sync:
	#pragma omp for
	for (size_t block=0; block<simClasses.pargrid.getNumberOfLocalCells(); ++block) {
	   // If all cell's neighbours exist, it is not on the boundary of the simulation domain:
	   if (simClasses.pargrid.getNeighbourFlags(block) == pargrid::ALL_NEIGHBOURS_EXIST) continue;
      
	   // Set boundary value:
	   for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	      data[block*block::SIZE+block::index(i,j,k)] = 0.0;
	   }
	}
   
	// Start remote cell data sync (master thread only):
	if (tid == 0) simClasses.pargrid.startNeighbourExchange(Advection::faceNbrStencilID,Advection::dataAvgsID);
   
	// Calculate inner block fluxes:
	if (tid == 0) profile::start("fluxes",fluxesID);
	const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(Advection::faceNbrStencilID);   
	#pragma omp for
	for (size_t block=0; block<innerBlocks.size(); ++block) calculateFlux(data,flux,sim,simClasses,innerBlocks[block]);
	if (tid == 0) profile::stop();
   
	// Wait for data sync to complete (master thread only):
	if (tid == 0) { 
	   profile::start("MPI waits",mpiWaitID);
	   simClasses.pargrid.wait(Advection::faceNbrStencilID,Advection::dataAvgsID);
	   profile::stop(); 
	}

	// Calculate boundary block fluxes:
	if (tid == 0) profile::start("fluxes",fluxesID);
	const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(Advection::faceNbrStencilID);   
	#pragma omp for
	for (size_t block=0; block<boundaryBlocks.size(); ++block) calculateFlux(data,flux,sim,simClasses,boundaryBlocks[block]);
	if (tid == 0) profile::stop();

	// Propagate:
	if (tid == 0) profile::start("propagation",propagID);   
	#pragma omp for
	for (size_t block=0; block<simClasses.pargrid.getNumberOfLocalCells(); ++block) propagate(data,flux,simClasses,block);
	if (tid == 0) profile::stop();
     } // End parallel construct, threads deleted

   profile::stop();
   return success;
}

/** Calculate fluxes for the given block.
 * @param data Array containing block data values.
 * @param flux Array in which computed fluxes are stored.
 * @param simClasses Struct containing general-use simulation classes.
 * @param blockID Local ID of the block whose fluxes are computed.*/
void calculateFlux(Real* data,Real* flux,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID) {
   // Skip cells on the boundary of simulation domain:
   if (simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   
   // Create a temporary data block:
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size];
   
   // Copy data to temporary block:
   block::fetchValues3D_face(simClasses,blockID,array,data);

   // Calculate flux for each cell in the block:
   Real cellSizes[3];
   getBlockCellSize(simClasses,sim,blockID,cellSizes);
   const Real CONST_X = -0.5*sim.dt/cellSizes[0]*Advection::Vx;
   const Real CONST_Y = -0.5*sim.dt/cellSizes[1]*Advection::Vy;
   const Real CONST_Z = -0.5*sim.dt/cellSizes[2]*Advection::Vz;
   
   for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
      flux[blockID*block::SIZE+block::index(i,j,k)] = CONST_X*(array[block::arrayIndex(i+2,j+1,k+1)] - array[block::arrayIndex(i+0,j+1,k+1)])
	                                            + CONST_Y*(array[block::arrayIndex(i+1,j+2,k+1)] - array[block::arrayIndex(i+1,j+0,k+1)])
	                                            + CONST_Z*(array[block::arrayIndex(i+1,j+1,k+2)] - array[block::arrayIndex(i+1,j+1,k+0)]);
   }
}

/** Propagate given block forward in time.
 * @param data Array containing block data values.
 * @param flux Array containing computed fluxes for each block.
 * @param blockID Local ID of the propagated block.*/
void propagate(Real* data,Real* flux,SimulationClasses& simClasses,pargrid::CellID blockID) {
   // Skip cells on the boundary of simulation domain:
   if (simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;

   Real* const dataTmp = data + blockID*block::SIZE;
   const Real* const fluxTmp = flux + blockID*block::SIZE;
   
   // Propagate each cell in the block:
   for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
      dataTmp[block::index(i,j,k)] += fluxTmp[block::index(i,j,k)];
   }
}
