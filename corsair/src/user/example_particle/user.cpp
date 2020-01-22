/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2015 Finnish Meteorological Institute
 *  Copyright 2016 Arto Sandroos
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
#include <particle_list_skeleton.h>
#include <dataoperatorcontainer.h>
#include <logically_cartesian_builder.h>
#include <operator_mpirank.h>
#include <operator_particle.h>
#include <operator_load.h>
#include <operator_cellid.h>
#include <operator_pargrid_array.h>
#include <operator_timeseries_load.h>
#include <operator_EM_analytic.h>

#include "fields_B_constant.h"
#include "simcontrol.h"
#include "particle_definition.h"
#include "particle_species.h"
#include <test_injector.h>
#include <default_injector.h>
#include "particle_accumulator.h"
#include "particle_propagator_rk2_gc.h"

using namespace std;

const int ORDER = 1;

typedef ConstantB Field;
typedef TestInjector<Field,Species,Prticle<Real> > TEST_INJECTOR;
typedef DefaultInjector<Field,Species,Prticle<Real> > DEFAULT_INJECTOR;
typedef ParticlePropagatorRk2GC<Field,Species,Prticle<Real> > RK2_GC_PROPAGATOR;
typedef Accumulator<Species,Prticle<Real>,ORDER> ACCUM;

//typedef TEST_INJECTOR INJECTOR;
typedef DEFAULT_INJECTOR INJECTOR;
typedef ParticleListSkeleton<Species,Prticle<Real> > PARTICLE_LIST;

bool propagate(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   bool rvalue = true;

   // Propagate all particles:
   for (size_t p=0; p<particleLists.size(); ++p) 
     if (particleLists[p]->propagateBoundaryCellParticles() == false) rvalue = false;
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->propagateInnerCellParticles() == false) rvalue = false;
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->clearAccumulationArrays() == false) rvalue = false;
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->waitParticleSends() == false) rvalue = false;

   // Accumulate particle quantities to simulation mesh:
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->accumulateBoundaryCells() == false) rvalue = false;
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->accumulateInnerCells() == false) rvalue = false;
   
   // Apply boundary conditions:
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->applyBoundaryConditions() == false) rvalue = false;
   
   // Inject new particles:
   for (size_t p=0; p<particleLists.size(); ++p)
     if (particleLists[p]->injectParticles() == false) rvalue = false;
   return rvalue;
}

bool userEarlyInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,vector<ParticleListBase*>& particleLists) {
   return true;
}

bool userLateInitialization(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			    const ObjectFactories& objectFactories,vector<ParticleListBase*>& particleLists) {
   bool success = true;
   // Add propagated particle species to particle list:
   particleLists.push_back(new PARTICLE_LIST);
   
   // Init particle lists added above:
   if (particleLists[particleLists.size()-1]->initialize(sim,simClasses,cr,objectFactories,"proton") == false) success = false;

   // Create a data array where particle density is accumulated:
   SimControl::densityDataID = simClasses.pargrid.addUserData<Real>("density",block::SIZE);
   if (SimControl::densityDataID == pargrid::INVALID_DATAID) {
      simClasses.logger << "(USER) ERROR: Failed to create pargrid static data array!" << endl << write;
      success = false;
   }
   
   // Clear data array to zero values:
   Real* data = simClasses.pargrid.getUserDataStatic<Real>(SimControl::densityDataID);
   for (pargrid::CellID blockLID=0; blockLID<simClasses.pargrid.getNumberOfAllCells(); ++blockLID) {
      for (int i=0; i<block::SIZE; ++i) {
         data[blockLID*block::SIZE+i] = 0.0;
      }
   }
   
   // Create a stencil for transferring density. Here we say that each
   // cell receives updates from each of its remote neighbours:
   vector<pargrid::NeighbourID> nbrIDs;
   for (int k=-1; k<2; ++k) for (int j=-1; j<2; ++j) for (int i=-1; i<2; ++i) {
      if (i == 1 && (j == 1 && k == 1)) continue;
      nbrIDs.push_back(simClasses.pargrid.calcNeighbourTypeID(i,j,k));
   }
   
   // Create a new stencil:
   SimControl::densityStencilID = simClasses.pargrid.addStencil(pargrid::remoteToLocalUpdates,nbrIDs);
   if (SimControl::densityStencilID == pargrid::INVALID_STENCILID) {
      simClasses.logger << "(USER) ERROR: Failed to create a ParGrid stencil!" << endl << write;
      success = false;
   }
   
   // Associate a data transfer between density array and the stencil created above:
   if (simClasses.pargrid.addDataTransfer(SimControl::densityDataID,SimControl::densityStencilID) == false) {
      simClasses.logger << "(USER) ERROR: Failed to add ParGrid data transfer!" << endl << write;
      success = false;
   }

   return success;
}

bool userFinalization(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   return true;
}

bool registerObjectMakers(ObjectFactories& objectFactories) {
   bool success = true;
   
   if (objectFactories.gridBuilders.registerMaker("LogicallyCartesian",LCCreator) == false) success = false;
   if (objectFactories.particleAccumulators.registerMaker("Accumulator",AccumulatorMaker<Species,Prticle<Real>,ORDER>) == false) success = false;
   if (objectFactories.particleBoundaryConditions.registerMaker("BoundaryCondRemove",ParticleBoundaryCondBaseMaker) == false) success = false;
   if (objectFactories.particleInjectors.registerMaker("DefaultInjector",DefInjectorMaker<Field,Species,Prticle<Real> >) == false) success = false;
   if (objectFactories.particlePropagators.registerMaker("RK2GC",RK2GCMaker<Field,Species,Prticle<Real> >) == false) success = false;
   return success;
}

bool registerDataOperators() {
   bool success = true;

   DataOperatorContainer& doc = corsair::getObjectWrapper().dataOperatorContainer;
   if (doc.registerOperator(new MPIRank) == false) success = false;
   if (doc.registerOperator(new CellIDOP) == false) success = false;
   if (doc.registerOperator(new EMAnalyticOperator<Field>) == false) success = false;
   if (doc.registerOperator(new LoadOP) == false) success = false;
   if (doc.registerOperator(new ParticleOperator) == false) success = false;
   if (doc.registerOperator(new LoadTSeriesOP) == false) success = false;
   if (doc.registerOperator(new OperatorPargridArray) == false) success = false;
   
   return success;
}

bool userRunTests(Simulation& sim,SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists) {return true;}
