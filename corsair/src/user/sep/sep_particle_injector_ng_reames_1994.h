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

#ifndef SEP_PARTICLE_INJECTOR_NG_REAMES_1994_H
#define SEP_PARTICLE_INJECTOR_NG_REAMES_1994_H

#include <cstdlib>
#include <climits>
#include <stdint.h>
#include <map>
#include <typeinfo>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include <ucd_mesh.h>
#include <base_class_particle_propagator.h>
#include <main.h>

#include "sep_object_wrapper.h"
#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_particle_definition.h"
#include "sep_fields_container.h"
#include "sep_coordinate_transform.h"
#include "sep_injection_buffer.h"

// TEST
#include "sep_shock_spherical.h"
#include "sep_particle_propagator_coronal_rk2.h"
#include "sep_particle_accelerator.h"
// END TEST

namespace sep {

   extern sep::SimControl simControl;

   struct InjectionRange {
      int32_t indexMin;
      int32_t indexMax;

      InjectionRange() { }
      InjectionRange(const int32_t& indexMin,const int32_t& indexMax): indexMin(indexMin),indexMax(indexMax) { }
   };
   
   template<class SPECIES,class PARTICLE>
   class ParticleInjectorNgReames1994: public ParticleInjectorBase {
    public: 
      ParticleInjectorNgReames1994();
      ~ParticleInjectorNgReames1994();
   
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      Real getMaximumSpeed();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);

    private:
      bool energyPerNucleon;
      bool initialized;
      bool initialInjectionDone;
      Real relativeAbundance;
      SPECIES species;

      uint32_t N_particlesPerCell;

      // TEST
      const ParticleListBase* plist;
      ParticlePropagCoronalRK2<SPECIES,PARTICLE,ParticleAccelerator<SPECIES,PARTICLE> >* propagator;
      ShockSpherical* shock;
      // END TEST

      void* energyParams;                            /**< Parameters for energy distribution.*/
      sep::finalizeEnergyDistrib finalizeEnergy;     /**< Pointer to energy distribution finalizer.*/
      sep::getEnergyFunction getEnergy;              /**< Pointer to energy distribution function.*/
      
      const Real DEF_VALUE;
      
      Real getInjectionPitch() const;
      void injectParticles(pargrid::CellID block,unsigned int* N_particles,pargrid::DataWrapper<PARTICLE>& wrapper,InjectionBuffer<PARTICLE>* injBuffer);
      bool intersects(Real* shockCentroid,Real shockRadius,int32_t i_node,int32_t j_node,int32_t k_node);
      bool intersects(Real* shockCentroid,Real shockRadius,
		      int32_t i1,int32_t j1,int32_t k1,
		      int32_t i_off,int32_t j_off,int32_t k_off);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* PINgReames1994Maker() {return new ParticleInjectorNgReames1994<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorNgReames1994<SPECIES,PARTICLE>::ParticleInjectorNgReames1994(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity())  {
      propagator = NULL;
      shock = NULL;
      finalizeEnergy = NULL;
      getEnergy = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorNgReames1994<SPECIES,PARTICLE>::~ParticleInjectorNgReames1994() {
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".energy_distribution","Name of energy distribution function (string)",std::string(""));
      cr.add(PREFIX+".energy_distribution_parameters","Name of region containing energy distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".energy_per_nucleon","If 'yes' injection energies are per nucleon, defaults to 'no' (string)",std::string("no"));
      cr.add(PREFIX+".macroparticles_per_cell","Number of injected macroparticles per cell (int)",(uint32_t)1);
      cr.add(PREFIX+".injection_propagator_name","Name of propagator that is used in injection (string)",std::string(""));
      cr.add(PREFIX+".relative_abundance","Relative abundance of injected species relative to plasma number density (float)",(Real)1.0);
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      initialInjectionDone = false;

      if (finalizeEnergy != NULL) (*finalizeEnergy)(energyParams);
      finalizeEnergy = NULL;
      getEnergy = NULL;
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   Real ParticleInjectorNgReames1994<SPECIES,PARTICLE>::getInjectionPitch() const {
      return -1.0 + 2*simClasses->random.uniform();
   }

   template<class SPECIES,class PARTICLE> inline
   Real ParticleInjectorNgReames1994<SPECIES,PARTICLE>::getMaximumSpeed() {
      if (sim->t >= simControl.t_setup) {
	 return sqrt(2*1e7*constants::CHARGE_ELEMENTARY/species.mass);
      }
      return 1.0;
   }
   
   /** Called by ParticleListSkeleton.*/
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      bool success = true;

      if (sim->t <= simControl.t_setup) return success;

      // Get injection buffer:
      extern std::map<std::string,InjectionBuffer<PARTICLE> > particleInjectionBuffers;
      particleInjectionBuffers[species.getName()];
      typename std::map<std::string,InjectionBuffer<PARTICLE> >::iterator it
	= particleInjectionBuffers.find(species.getName());
      
      InjectionBuffer<PARTICLE>* injBuffer = &(it->second);
      
      // TEST
      if (propagator == NULL) {
	 typedef ParticleAccelerator<SPECIES,PARTICLE> ACCELERATOR;
	 ParticlePropagatorBase* basePropagator = plist->getPropagator();
	 if (typeid(*basePropagator) != typeid(ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>)) {
	    simClasses->logger << "(SEP PARTICLE INJ NG REAMES) ERROR: Wrong propagator" << std::endl << write;
	    exit(1);
	 }
	 
	 propagator = reinterpret_cast<ParticlePropagCoronalRK2<SPECIES,PARTICLE,ACCELERATOR>*>(basePropagator);
      }
      
      if (shock == NULL) {
	 ShockBaseClass* baseShock = simControl.shock;
	 if (typeid(*baseShock) != typeid(ShockSpherical)) {
	    simClasses->logger << "(SEP PARTICLE INJ NG REAMES) ERROR: Wrong type of shock" << std::endl << write;
	    exit(1);
	 }
	 
	 shock = reinterpret_cast<ShockSpherical*>(simControl.shock);
      }
      // END TEST

      recalculateCellVolumes(*sim,*simClasses);
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);

      #ifndef NDEBUG
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP PARTICLE INJ NG REAMES) ERROR: Shock is NULL" << std::endl << write;
	 return false;
      }
      #endif

      // Get shock centroid in Cartesian coordinates and its radius:
      Real shockCentroid[3];
      //Real shockCentroidSpherical[3];
      Real shockRadius;
      shock->getShockState(sim->t+sim->dt,shockCentroid,shockRadius);
      //transformPositionCartesianToSpherical(shockCentroid,shockCentroidSpherical);
      //getLogicalCoordinates(sim,shockCentroidSpherical,shockCentroid);

      Real t_propag = 0.0;
      const std::vector<pargrid::CellID>& interiorBlocks = simClasses->pargrid.getInteriorCells();
      for (size_t b=0; b<interiorBlocks.size(); ++b) {
      //for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	 const pargrid::CellID blockLID = interiorBlocks[b];
	 int32_t i_block,j_block,k_block;
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

	 int32_t i_cell = i_block*block::WIDTH_X;
	 int32_t j_cell = j_block*block::WIDTH_Y;
	 int32_t k_cell = k_block*block::WIDTH_Z;
	 if (intersects(shockCentroid,shockRadius,i_cell,j_cell,k_cell) == false) continue;

	 //std::cerr << "accepted " << i_cell << ' ' << j_cell << ' ' << k_cell << std::endl;
	 
	 /*
	 if ((i_block+1)*block::WIDTH_X < shockCentroid[0]) continue;
	 
	 
	 Real pos[3];
	 pos[0] = i_block*block::WIDTH_X;
	 pos[1] = j_block*block::WIDTH_Y;
	 pos[2] = k_block*block::WIDTH_Z;

	 bool injectionBlock = false;
	 const int firstRegion = simControl.shock->getShockRegion(sim->t+sim->dt,pos);
	 
	 for (int k=0; k<2; ++k) for (int j=0; j<2; ++j) for (int i=0; i<2; ++i) {
	    pos[0] = i_block*block::WIDTH_X + i*block::WIDTH_X;
	    pos[1] = j_block*block::WIDTH_Y + j*block::WIDTH_Y;
	    pos[2] = k_block*block::WIDTH_Z + k*block::WIDTH_Z;
	    if (simControl.shock->getShockRegion(sim->t+sim->dt,pos) != firstRegion) {
	       injectionBlock = true;
	       break;
	    }
	 }
	 if (injectionBlock == false) continue;*/

	 
	 
	 // Skip blocks that are behind shock:
	 // TODO
	 
	 // Measure block injection time if we are testing for repartitioning:
	 if (sim->countPropagTime == true) t_propag = MPI_Wtime();
	 
	 //injectParticles(blockLID,N_particles,wrapper);
	 injectParticles(blockLID,N_particles,wrapper,injBuffer);

	 // Store block injection time:
	 if (sim->countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }

      return success;
   }

   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorNgReames1994<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,unsigned int* N_particles,
									pargrid::DataWrapper<PARTICLE>& wrapper,
									InjectionBuffer<PARTICLE>* injBuffer) {
      // Get block global ID and calculate bounding box (i,j,k) indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      uint32_t i_block,j_block,k_block;
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);

      //std::cerr << "step " << sim->timestep << " gid#" << blockGID << ' ' << i_block << ' ' << j_block << ' ' << k_block << std::endl;
      
      /*if (blockGID == 813) {
	 std::cerr << "block#" << blockGID << " is injection step " << sim->timestep << std::endl;
	 exit(1);
      }*/
      
      // Make a trial particle that is propagated by one time step.
      // After propagation it is checked if trial particle got 
      // reflected off shock. If yes, it is injected as a macroparticle:
      PARTICLE trial;
      
      distrib::InjectionEnergy injEnergy;
      PlasmaState plasmaState;

      for (int32_t k=0; k<block::WIDTH_Z; ++k) for (int32_t j=0; j<block::WIDTH_Y; ++j) for (int32_t i=0; i<block::WIDTH_X; ++i) {
	 //pargrid::ArraySizetype oldSize = wrapper.size()[blockLID];
	 //pargrid::ArraySizetype counter = 0;
	 
	 // Calculate cell (i,j,k) indices:
	 const uint32_t i_cell = i_block*block::WIDTH_X + i;
	 const uint32_t j_cell = j_block*block::WIDTH_Y + j;
	 const uint32_t k_cell = k_block*block::WIDTH_Z + k;

	 // Get cell volume:
	 const Real volume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	 
	 // Inject particles to cell:
	 //Real weightSum = 0.0;
	 for (uint32_t n=0; n<N_particlesPerCell; ++n) {
	    // Calculate injection position in logical coordinates:
	    trial.state[particle::XCRD] = i_cell + simClasses->random.uniform();
	    trial.state[particle::YCRD] = j_cell + simClasses->random.uniform();
	    trial.state[particle::ZCRD] = k_cell + simClasses->random.uniform();

	    // Do not inject particles to downstream of shock:
	    if (simControl.includeShock == true) {
	       if (simControl.shock->inDownstream(sim->t+sim->dt,trial.state) == true) continue;
	    }

	    // Get plasma state at injection position:
	    (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,trial.state,plasmaState);
	    injEnergy.thermalEnergy = constants::BOLTZMANN*plasmaState.ionTemperature;

	    // Calculate pitch, energy:
	    (*getEnergy)(injEnergy,energyParams);
	    Real pitch = getInjectionPitch();
	    
	    // Calculate parallel speed and magnetic moment:
	    Real speed;
	    if (energyPerNucleon == true) {
	       speed = sqrt(2*injEnergy.energy/constants::MASS_PROTON);
	    } else {
	       speed = sqrt(2*injEnergy.energy/species.mass);
	    }
	    const Real parallelSpeed  = speed * pitch;
	    const Real gyroSpeed      = speed * sqrt(1 - pitch*pitch);
	    const Real B_mag = vectorMagnitude<3>(plasmaState.B);
	    const Real magneticMoment = 0.5*species.mass*gyroSpeed*gyroSpeed / B_mag;

	    trial.state[particle::V_PAR] = parallelSpeed + 
	      dotProduct<3>(plasmaState.V_plasma_SIM,plasmaState.B)/B_mag;
	    trial.state[particle::MU]    = magneticMoment;

	    // Propagate trial particle from t to t+dt and check if it got reflected:
	    PARTICLE tmp = trial;
	    if (propagator->propagateParticle(blockLID,sim->t,tmp) == false) {
	       continue;
	    }

	    /*if ((U_new-U_old)/U_old >= 1.5) {
	       std::cerr << "step: " << sim->timestep << ' ' << i_block << ' ' << j_block << ' ' << k_block;
	       std::cerr << " gain " << (U_new-U_old)/U_old << std::endl;
	    }*/
	    
	    //std::cerr << "step: " << sim->timestep << ' ' << i_block << ' ' << j_block << ' ';
	    //std::cerr << trial.state[particle::V_PAR] << '\t' << tmp.state[particle::V_PAR] << std::endl;

	    // Calculate macroparticle statistical weight == number of real particles 
	    // represented by the macroparticle. injEnergy.weight is the value of 
	    // injection energy distribution function at injection energy, normalized 
	    // to unity. It is multiplied by local plasma number density and spatial 
	    // cell volume and relative abundance:
	    trial.state[particle::WEIGHT] = injEnergy.weight*volume
	      * (plasmaState.ionMassDensity / constants::MASS_PROTON)
		/ N_particlesPerCell * relativeAbundance;
	    
	    tmp.state[particle::WEIGHT] = trial.state[particle::WEIGHT];

	    // Trial particle accepted, inject it to simulation:
	    const Real U_old = 0.5*species.mass*trial.state[particle::V_PAR]*trial.state[particle::V_PAR]
	      + trial.state[particle::MU]*B_mag;
	    const Real U_new = 0.5*species.mass*tmp.state[particle::V_PAR]*tmp.state[particle::V_PAR]
	      + tmp.state[particle::MU]*B_mag;
	    const Real energyGain = (U_new-U_old)/U_old;

	    // TEST
	    // Calculate particle's i,j,k offset indices in block-based mesh:
	    const int32_t I_off = static_cast<uint32_t>(1 + (tmp.state[0]-i_block*block::WIDTH_X) / block::WIDTH_X) - 1;
	    const int32_t J_off = static_cast<uint32_t>(1 + (tmp.state[1]-j_block*block::WIDTH_Y) / block::WIDTH_Y) - 1;
	    const int32_t K_off = static_cast<uint32_t>(1 + (tmp.state[2]-k_block*block::WIDTH_Z) / block::WIDTH_Z) - 1;
	    
	    // Calculate local ID of block where accepted particle ended up:
	    const pargrid::NeighbourID nbrOffset = simClasses->pargrid.calcNeighbourTypeID(I_off,J_off,K_off);
	    const pargrid::CellID newCell = simClasses->pargrid.getCellNeighbourIDs(blockLID)[nbrOffset];
	    if (newCell == pargrid::INVALID_CELLID) continue;
	    // END TEST
	    
	    // If trial particle received a large enough energy gain, split it:
	    const int splitAmount = 100;
	    if (energyGain >= 0.5) {
	       //std::cerr << "INJ energy " << U_new/constants::CHARGE_ELEMENTARY/1e6 << " injector " << injEnergy.energy/constants::CHARGE_ELEMENTARY/1e6 << std::endl;
	       trial.state[particle::WEIGHT] /= splitAmount;
	       tmp.state[particle::WEIGHT] /= splitAmount;
	       //N_particles[blockLID] += splitAmount;
	       //for (int cntr=0; cntr<splitAmount; ++cntr) wrapper.push_back(blockLID,trial);
	       for (int cntr=0; cntr<splitAmount; ++cntr) injBuffer->insert(newCell,tmp);
	    } else if (energyGain >= 0.25) {
	       const int32_t splits = std::min(1,splitAmount/10);
	       tmp.state[particle::WEIGHT] /= splits;
	       for (int cntr=0; cntr<splits; ++cntr) injBuffer->insert(newCell,tmp);
	    } else {
	       //++N_particles[blockLID];
	       //wrapper.push_back(blockLID,trial);
	       injBuffer->insert(newCell,tmp);
	    }
	    //++counter;
	 }
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								   const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP PARTICLE INJ NG REAMES) Starting to init." << std::endl << write;
      
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());
      
      // Parse config file options:
      std::string propagatorName,energyPerNucleonString;
      std::string energyDistribString,energyDistribParamsString;
      std::string pitchDistribString,pitchDistribParamsString;
      addConfigFileItems(cr,regionName);
      cr.parse();
      cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell);
      cr.get(regionName+".injection_propagator_name",propagatorName);
      cr.get(regionName+".energy_distribution",energyDistribString);
      cr.get(regionName+".energy_distribution_parameters",energyDistribParamsString);
      cr.get(regionName+".energy_per_nucleon",energyPerNucleonString);
      cr.get(regionName+".relative_abundance",relativeAbundance);

      energyPerNucleon = false;
      if (energyPerNucleonString == "yes") energyPerNucleon = true;
      
      if (relativeAbundance <= 0.0) {
	 simClasses.logger << "(SEP PARTICLE INJ NG REAMES) ERROR: Relative abundance must be positive" << std::endl << write;
	 initialized = false;
      }

      // Attempt to get energy distribution function:
      sep::initializeEnergyDistrib initializeEnergy;
      sep::getEnergyDistribFunction getDistrib;
      if (sep::getObjectWrapper().energyDistribContainer.getDistribution(energyDistribString,finalizeEnergy,
									 getDistrib,getEnergy,initializeEnergy) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ HOMOG) ERROR: Could not find energy distribution ";
	 simClasses.logger << "function called '";
	 simClasses.logger << energyDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".energy_distribution'" 
	   << std::endl << write;
	 initialized = false;
	 return initialized;
      }
      if (initializeEnergy(sim,simClasses,cr,energyDistribParamsString,energyParams) == false) {
	 simClasses.logger << "(PARTICLE INJ HOMOG) ERROR: Energy distribution function failed ";
	 simClasses.logger << "to initialize" << std::endl << write;
	 initialized = false;
	 finalizeEnergy(energyParams);
      }
      
      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP PARTICLE INJ NG REAMES) ERROR: Parameter '" << regionName+".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
      }
      
      // Write out injection parameters:
      simClasses.logger << "(SEP PARTICLE INJ NG REAMES) Injection parameters are:" << std::endl;
      simClasses.logger << write;
      
      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }

      this->plist = plist;
      return initialized;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::intersects(Real* shockCentroid,Real shockRadius,
								   int32_t i_node,int32_t j_node,int32_t k_node) {

      if (intersects(shockCentroid,shockRadius,i_node  ,j_node  ,k_node  ,1,0,0) == true) return true; // 01
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node  ,k_node  ,0,1,0) == true) return true; // 03
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node  ,k_node  ,0,0,1) == true) return true; // 04
      if (intersects(shockCentroid,shockRadius,i_node+1,j_node  ,k_node  ,0,1,0) == true) return true; // 12
      if (intersects(shockCentroid,shockRadius,i_node+1,j_node  ,k_node  ,0,0,1) == true) return true; // 15
      if (intersects(shockCentroid,shockRadius,i_node+1,j_node+1,k_node  ,0,0,1) == true) return true; // 26
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node+1,k_node  ,1,0,0) == true) return true; // 32
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node+1,k_node  ,0,0,1) == true) return true; // 37
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node  ,k_node+1,1,0,0) == true) return true; // 45
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node  ,k_node+1,0,1,0) == true) return true; // 47
      if (intersects(shockCentroid,shockRadius,i_node+1,j_node  ,k_node+1,0,1,0) == true) return true; // 56
      if (intersects(shockCentroid,shockRadius,i_node  ,j_node+1,k_node+1,1,0,0) == true) return true; // 76
      return false;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorNgReames1994<SPECIES,PARTICLE>::intersects(Real* shockCentroid,Real shockRadius,
								   int32_t i1,int32_t j1,int32_t k1,
								   int32_t i_off,int32_t j_off,int32_t k_off) {
      // Node origin in spherical coords:
      Real spherOrigin[3];
      spherOrigin[0] = sim->x_crds_node[i1];
      spherOrigin[1] = sim->y_crds_node[j1];
      spherOrigin[2] = sim->z_crds_node[k1];
      
      // Direction vector in spherical coords:
      Real spherDir[3];
      spherDir[0] = i_off*sim->dx_cell[i1];
      spherDir[1] = j_off*sim->dy_cell[j1];
      spherDir[2] = k_off*sim->dz_cell[k1];
      
      // Transform direction vector to Cartesian coords:
      Real cartDir[3];
      transformVectorSphericalToCartesian(spherOrigin,spherDir,cartDir);
      
      // Calculate unit vector:
      Real unitCartDir[3];
      for (int i=0; i<3; ++i) unitCartDir[i] = cartDir[i];
      unitVector<3>(unitCartDir);
      
      // Transform node origin to Cartesian coords:
      Real cartOrigin[3];
      transformPositionSphericalToCartesian(spherOrigin,cartOrigin);
      
      // Node origin - shock centroid vector:
      for (int i=0; i<3; ++i) cartOrigin[i] = cartOrigin[i] - shockCentroid[i];
      
      Real l_dot_pos = dotProduct<3>(unitCartDir,cartOrigin);
      Real rooted = l_dot_pos*l_dot_pos - vectorMagnitude2<3>(cartOrigin) + shockRadius*shockRadius;

      // Return immediately if rooted value is negative:
      if (rooted < 0) return false;
      
      // If distance to both intersection points is longer 
      // than edge length, cell does not cross shock surface:
      Real crossing1 = -l_dot_pos + sqrt(rooted);
      Real crossing2 = -l_dot_pos - sqrt(rooted);
      Real edgeLength = vectorMagnitude<3>(cartDir);
      if (fabs(crossing1) > edgeLength && fabs(crossing2) > edgeLength) return false;
      
      return true;
   }

} // namespace sep
   
#endif
