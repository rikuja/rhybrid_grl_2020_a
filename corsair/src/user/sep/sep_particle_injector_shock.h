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

#ifndef SEP_PARTICLE_INJECTOR_SHOCK_H
#define SEP_PARTICLE_INJECTOR_SHOCK_H

#include <cstdlib>
#include <climits>
#include <stdint.h>
#include <map>
#include <set>
#include <vector>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include <ucd_mesh.h>
#include <base_class_particle_propagator.h>
#include <main.h>
#include <ucd_mesh.h>

#include "sep_object_wrapper.h"
#include "sep_namespaces.h"
#include "sep_simcontrol.h"
#include "sep_propagate.h"
#include "sep_particle_definition.h"
#include "sep_coordinate_transform.h"
#include "sep_injection_buffer.h"
#include "sep_injectors_common.h"

// TEST
#include "sep_shock_spherical.h"
#include "sep_particle_propagator_coronal_rk2.h"
#include "sep_particle_accelerator.h"
#include "sep_base_class_particle_splitter.h"
// END TEST

namespace sep {

   extern sep::SimControl simControl;

   /** Particle injector class for SEP simulations. This injector is capable 
    * of injecting particles directly to the upstream of a shock wave. 
    * Injection is done on triangulated shock surface elements.*/
   template<class SPECIES,class PARTICLE>
   class ParticleInjectorShock: public ParticleInjectorBase {
    public: 
      ParticleInjectorShock();
      ~ParticleInjectorShock();
   
      bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
      bool finalize();
      Real getMaximumSpeed();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		      const std::string& regionName,const ParticleListBase* plist);
      bool inject(pargrid::DataID particleDataID,unsigned int* N_particles);

    private:
      bool initialized;
      bool initialInjectionDone;
      bool injectInitial;                               /**< If true, particles are injected everywhere 
							 * to initialize simulation.*/
      Real relativeAbundance;
      SPECIES species;
      particleinjector::frame::Frame injectionFrame;    /**< Frame in which particles are injected, 
						         * affects the energy and pitch of particles in 
							 * simulation frame.*/
      particleinjector::region::Region injectionRegion; /**< Position where particles are injected to.
							 * Default behaviour is to inject particles to the 
							 * shock front.*/
      uint32_t N_particlesPerCell;
      uint32_t inflowBoundaryMask;
      Real densityMultiplier;                           /**< Multiplier applied to particle densities.*/
      bool injectReflectedOnly;                         /**< If true, only particles that reflect off the shock 
							 * to the upstream region are accepted.*/
      Real surfaceNormal_x;
      Real t_injection_min;
      
      void* energyParams;                               /**< Parameters for energy distribution.*/
      sep::finalizeEnergyDistrib finalizeEnergy;        /**< Pointer to energy distribution finalizer.*/
      sep::getEnergyFunction getEnergy;                 /**< Pointer to energy distribution function.*/
      sep::finalizePitchDistrib finalizePitch;          /**< Pointer to pitch distribution finalizer.*/
      sep::getPitchFunction getPitch;                   /**< Pointer to pitch distribution function.*/
      const Real DEF_VALUE;

      void generateParticle(PARTICLE& particle,const Real* V_frame,distrib::InjectionEnergy& injEnergy,
			    const PlasmaState& plasmaState,const Real& volume);
      void getFrameSpeed(pargrid::CellID blockLID,Real t_inj,Real* pos,Real* B,Real* V_frame,Real* V_wave,Real& V_alfven);
      void injectBlock(pargrid::CellID blockLID,pargrid::DataID particleDataID,unsigned int* N_particles,const bool& adjustInjectionRate);
      void injectParticles(pargrid::CellID blockLID,const Real& t_inj,const Real* const logicalCentroid,const Real* const normal,
			   const Real* const V_surface,
			   const Real& area,InjectionBuffer<PARTICLE>* injBuffer,int direction);
      bool injectShockUpstream(pargrid::DataID particleDataID,unsigned int* N_particles,int direction);
      bool injectUpstream(pargrid::DataID particleDataID,unsigned int* N_particles);
   };

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorBase* PIShockMaker() {return new ParticleInjectorShock<SPECIES,PARTICLE>();}
   
   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorShock<SPECIES,PARTICLE>::ParticleInjectorShock(): ParticleInjectorBase(),DEF_VALUE(std::numeric_limits<Real>::infinity())  {
      energyParams = NULL;
      finalizeEnergy = NULL;
      finalizePitch = NULL;
      getEnergy = NULL;
      getPitch = NULL;
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   ParticleInjectorShock<SPECIES,PARTICLE>::~ParticleInjectorShock() {
      finalize();
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& PREFIX) {
      cr.add(PREFIX+".energy_distribution","Name of energy distribution function (string)",std::string(""));
      cr.add(PREFIX+".energy_distribution_parameters","Name of region containing energy distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution","Name of pitch distribution function (string)",std::string(""));
      cr.add(PREFIX+".pitch_distribution_parameters","Name of region containing pitch distribution parameters (string)",std::string(""));
      cr.add(PREFIX+".macroparticles_per_cell","Number of injected macroparticles per cell (int)",(uint32_t)1);
      cr.add(PREFIX+".relative_abundance","Relative abundance of injected species relative to plasma number density (float)",(Real)1.0);
      cr.add(PREFIX+".injection_frame","Injection frame, antiparallel/parallel/simulation (string)",std::string("simulation"));
      cr.add(PREFIX+".injection_region","Injection region upstream/downstream/shock_upstream/shock_downstream (string)",std::string("shock_upstream"));
      cr.add(PREFIX+".inflow_boundary","Inflow boundary (string)",std::string(""));
      cr.add(PREFIX+".density_multiplier","Multiplier applied to particle density (float)",(Real)1.0);
      cr.add(PREFIX+".inject_reflected_only","If 'yes' only reflected particles are accepted (string)",std::string("no"));
      cr.add(PREFIX+".surface_normal_x","If defined, particles are only injected on surfaces where normal has the same sign as this value (float)",(Real)0.0);
      cr.add(PREFIX+".injection_time_min","Minimum injection time in seconds (float)",(Real)0.0);
      return true;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      initialInjectionDone = false;
      injectInitial = false;

      if (finalizeEnergy != NULL) (*finalizeEnergy)(energyParams);
      if (finalizePitch != NULL) (*finalizePitch)();
      finalizeEnergy = NULL;
      getEnergy = NULL;
      finalizePitch = NULL;
      return true;
   }

   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorShock<SPECIES,PARTICLE>::generateParticle(PARTICLE& particle,const Real* V_frame,distrib::InjectionEnergy& injEnergy,
								  const PlasmaState& plasmaState,const Real& volume) {
      const Real B_mag = vectorMagnitude<3>(plasmaState.B);
      Real B[3];
      for (int i=0; i<3; ++i) B[i] = plasmaState.B[i];
      unitVector<3>(B);

      // Calculate pitch, energy:
      (*getEnergy)(injEnergy,energyParams);
      
      //std::cerr << injEnergy.weight << std::endl;
      
      const Real pitch = (*getPitch)();
      const Real V_mag = sqrt(2*injEnergy.energy/species.mass);
      
      // Calculate injection velocity in SIM frame:
      Real V_injection_SIM[3];
      for (int i=0; i<3; ++i) V_injection_SIM[i] = V_mag*pitch*B[i] + V_frame[i];

      // Calculate parallel speed and magnetic moment:
      const Real V_gyro2 = V_mag*V_mag*(1-pitch*pitch);
      particle.state[particle::V_PAR] = dotProduct<3>(V_injection_SIM,B);
      particle.state[particle::MU] = 0.5*species.mass*V_gyro2 / B_mag;

      // Calculate macroparticle weight:
      particle.state[particle::WEIGHT] = densityMultiplier*injEnergy.weight*plasmaState.ionMassDensity/species.mass*volume;
   }

   template<class SPECIES,class PARTICLE> inline
   Real ParticleInjectorShock<SPECIES,PARTICLE>::getMaximumSpeed() {
      if (sim->t < simControl.t_setup) return 1.0;
      if (sim->t < t_injection_min) return 1.0;

      return sqrt(2*1e7*constants::CHARGE_ELEMENTARY/species.mass);
   }

   /** Get injection frame speed.
    * @param blockLID Local ID of block where particles are injected to.
    * @param t_inj Injection time.
    * @param pos Injection position.
    * @param B Magnetic field at injection position (output).
    * @param V_frame Injection frame speed (output).
    * @param V_wave Wave speed at injection position (output).*/
   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorShock<SPECIES,PARTICLE>::getFrameSpeed(pargrid::CellID blockLID,Real t_inj,Real* pos,Real* B,Real* V_frame,Real* V_wave,Real& V_alfven) {
      Real dV_wave;
      PlasmaState plasmaState;
      switch (injectionFrame) {
       case particleinjector::frame::UNKNOWN:
	 simClasses->logger << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown injection frame" << std::endl << write;
	 exit(1);
	 break;
       case particleinjector::frame::SIMULATION:
	 (*simControl.fieldsGetState)(blockLID,t_inj,pos,B,V_wave,dV_wave,V_alfven,+1);
	 for (int i=0; i<3; ++i) V_frame[i] = 0;
	 break;
       case particleinjector::frame::ANTIPARALLEL:
	 (*simControl.fieldsGetState)(blockLID,t_inj,pos,B,V_wave,dV_wave,V_alfven,-1);
	 for (int i=0; i<3; ++i) V_frame[i] = V_wave[i];
	 break;
       case particleinjector::frame::PARALLEL:
	 (*simControl.fieldsGetState)(blockLID,t_inj,pos,B,V_wave,dV_wave,V_alfven,+1);
	 for (int i=0; i<3; ++i) V_frame[i] = V_wave[i];
	 break;
       case particleinjector::frame::PLASMA:
	 (*simControl.fieldsGetState)(blockLID,t_inj,pos,B,V_wave,dV_wave,V_alfven,+1);
	 (*simControl.fieldsGetPlasmaState)(blockLID,t_inj,pos,plasmaState);
	 for (int i=0; i<3; ++i) V_frame[i] = plasmaState.V_plasma_SIM[i];
	 break;
       default:
	 simClasses->logger << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown injection frame" << std::endl << write;
	 exit(1);
	 break;
      }
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::inject(pargrid::DataID particleDataID,unsigned int* N_particles) {
      if (initialized == false) return initialized;
      if (sim->t <= simControl.t_setup) return initialized;
      if (sim->t <= t_injection_min) return initialized;

      // If this is the first time step when particles are 
      // injected, force recalculation of simulation time step:
      static bool injectionStarted = false;
      if (injectionStarted == false) {
	 recalculateTimestep(*sim,*simClasses,corsair::getObjectWrapper().particleLists);
	 injectionStarted = true;
      }
      
      switch (injectionRegion) {
       case particleinjector::region::UNKNOWN:
	 std::cerr << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown or unsupported injection region" << std::endl;
	 exit(1);
	 break;
       case particleinjector::region::UPSTREAM:
	 return injectUpstream(particleDataID,N_particles);
	 break;
       case particleinjector::region::DOWNSTREAM:
	 std::cerr << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown or unsupported injection region" << std::endl;
	 exit(1);
	 break;
       case particleinjector::region::SHOCK_UPSTREAM:
	 return injectShockUpstream(particleDataID,N_particles,-1);
	 break;
       case particleinjector::region::SHOCK_DOWNSTREAM:
	 return injectShockUpstream(particleDataID,N_particles,-1);
	 break;
       default:
	 std::cerr << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown or unsupported injection region" << std::endl;
	 exit(1);
	 break;
      }
      return false;
   }
   
   /** Called by ParticleListSkeleton.*/
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::injectShockUpstream(pargrid::DataID particleDataID,unsigned int* N_particles,int direction) {
      bool success = true;

      const Real t_inj = sim->t + sim->dt;
      if (simControl.shock == NULL) {
	 simClasses->logger << "(SEP PARTICLE INJ SHOCK) ERROR: Pointer to shock is NULL" << std::endl << write;
	 return false;
      }
      
      // Get injection buffer:
      extern std::map<std::string,InjectionBuffer<PARTICLE> > particleInjectionBuffers;
      particleInjectionBuffers[species.getName()];
      typename std::map<std::string,InjectionBuffer<PARTICLE> >::iterator it
	= particleInjectionBuffers.find(species.getName());
            
      InjectionBuffer<PARTICLE>* injBuffer = &(it->second);
      // Get shock surface elements on this process and inject particles:
      std::vector<LocalSurface> localSurfaces;
      uint32_t shockCells;
      if (simControl.shock->getLocalSurfaces(t_inj,localSurfaces,shockCells,simControl.N_shockSurfaceRefinements) == false) {
	 success = false;
      } else {
	 for (size_t s=0; s<localSurfaces.size(); ++s) {
	    if (simControl.shock->getShockRegion(t_inj,localSurfaces[s].position) != 1) {
	       continue;
	    }

	    // Skip injection on surfaces where shock is not compressive:
	    if (simControl.shock->getGasCompressionRatio(localSurfaces[s].localID,t_inj,localSurfaces[s].position) <= 1.05) {
	       continue;
	    }

	    // Get shock velocity and normal vector at surface centroid and inject particles:
	    Real V_shock_SIM[3];
	    simControl.shock->getLocalShockVelocity(t_inj,localSurfaces[s].position,V_shock_SIM);
	    Real shockNormal[3];
	    simControl.shock->getShockNormal(t_inj,localSurfaces[s].position,shockNormal);

	    // Check if surface normal must point to required direction:
	    if (surfaceNormal_x != 0.0) {
	       if (surfaceNormal_x * shockNormal[0] < 0) continue;
	    }
	    
	    injectParticles(localSurfaces[s].localID,t_inj,localSurfaces[s].position,shockNormal,V_shock_SIM,localSurfaces[s].area,injBuffer,direction);
	 }
      }

      return success;
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::injectUpstream(pargrid::DataID particleDataID,unsigned int* N_particles) {
      bool success = true;
      
      #ifndef NDEBUG
         if (simControl.shock == NULL) {
	    simClasses->logger << "(SEP PARTICLE INJ SHOCK) ERROR: Pointer to shock is NULL" << std::endl << write;
	    return false;
	 }
      #endif
      
      // If initial injection has not been done, inject particles 
      // to all cells in shock upstream:
      if (initialInjectionDone == false) {
	 for (pargrid::CellID blockLID=0; blockLID<simClasses->pargrid.getNumberOfLocalCells(); ++blockLID) {
	    injectBlock(blockLID,particleDataID,N_particles,false);
	 }
	 initialInjectionDone = true;
	 return success;
      }

      // Get injection buffer:
      extern std::map<std::string,InjectionBuffer<PARTICLE> > particleInjectionBuffers;
      particleInjectionBuffers[species.getName()];
      typename std::map<std::string,InjectionBuffer<PARTICLE> >::iterator it
	= particleInjectionBuffers.find(species.getName());

      InjectionBuffer<PARTICLE>* injBuffer = &(it->second);
      
      const Real t_inj = sim->t + sim->dt;
      const std::vector<pargrid::CellID>& exteriorCells = simClasses->pargrid.getExteriorCells();
      for (size_t b=0; b<exteriorCells.size(); ++b) {
	 const pargrid::CellID blockLID = exteriorCells[b];
	 const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
	 int32_t indices[3];
	 block::calculateBlockIndices(*sim,blockGID,indices[0],indices[1],indices[2]);

	 const uint32_t existingNbrFlags = simClasses->pargrid.getNeighbourFlags()[blockLID];
	 const uint32_t missingNbrFlags = (existingNbrFlags ^ pargrid::ALL_NEIGHBOURS_EXIST);
	 if ((missingNbrFlags & inflowBoundaryMask) == 0) {
	    continue;
	 }
	 
	 Real normal[3];
	 for (int iii=0; iii<3; ++iii) normal[iii] = 0.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(-1,+0,+0))) > 0) normal[0] = +1.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+1,+0,+0))) > 0) normal[0] = -1.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,-1,+0))) > 0) normal[1] = +1.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+1,+0))) > 0) normal[1] = -1.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+0,-1))) > 0) normal[2] = +1.0;
	 if ((missingNbrFlags & (1 << simClasses->pargrid.calcNeighbourTypeID(+0,+0,+1))) > 0) normal[2] = -1.0;
	 unitVector<3>(normal);

	 Real V_surface[3] = {0,0,0};

	 Real position[3];
	 position[0] = indices[0]*block::WIDTH_X + 0.5*(1 + normal[0]);
	 position[1] = indices[1]*block::WIDTH_Y + 0.5*(1 + normal[1]);
	 position[2] = indices[2]*block::WIDTH_Z + 0.5*(1 + normal[2]);

	 Real area = 0;
	 int32_t I = 0;
	 Real R = 0;
	 switch (simControl.coordinateSystem) {
	  case sep::UNKNOWN:
	    std::cerr << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown coordinate system in injectUpstream" << std::endl;
	    exit(1);
	    break;
	  case sep::CARTESIAN:
	    area = sim->dy_cell[0]*sim->dz_cell[0];
	    break;
	  case sep::CYLINDRICAL:
	    std::cerr << "(SEP PARTICLE INJ SHOCK) ERROR: Cylindrical coordinates not supported in injectUpstream" << std::endl;
	    exit(1);
	    break;
	  case sep::SPHERICAL:
	    I = static_cast<int32_t>(position[0]);
	    R = sim->x_crds_node[I] + (position[0]-I)*sim->dx_cell[I];
	    area = -R*R*(cos(sim->y_crds_node[indices[1]+1])-cos(sim->y_crds_node[indices[1]])) * sim->dz_cell[indices[2]];
	    break;
	 }
	 
	 //injectBlock(blockLID,particleDataID,N_particles,true);
	 injectParticles(blockLID,t_inj,position,normal,V_surface,area,injBuffer,+1.0);
      }

      return success;
   }
   
   /**Inject particles uniformly to the given block.
    * @param blockLID Local ID of block where particles are injected to.
    * @param particleDataID ParGrid DataID of array that contains particles of this species.
    * @param N_particles Array that contains the number of macroparticles per block.
    * @param adjustInjectionRate If true, number of injected macroparticles is modified based
    * on the current simulation time step.
    */
   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorShock<SPECIES,PARTICLE>::injectBlock(pargrid::CellID blockLID,pargrid::DataID particleDataID,
							     unsigned int* N_particles,const bool& adjustInjectionRate) {
      // Measure computation time if we are testing for repartitioning:
      Real t_propag = 0.0;
      if (sim->countPropagTime == true) t_propag = MPI_Wtime();
      
      const Real t_inj = sim->t + sim->dt;
      PARTICLE trial;
      distrib::InjectionEnergy injEnergy;
      PlasmaState plasmaState;
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
      
      // Calculate block i,j,k indices:
      const pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
      uint32_t i_block,j_block,k_block;
      block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
      
      bool checkWaveEnergies = false;
      std::vector<Real> minWaveEnergies;
      if (simControl.applyWaveGrowth == true) {
	 checkWaveEnergies = getMinimumWaveEnergies(blockLID,+1,minWaveEnergies);
      }

      for (int32_t k=0; k<block::WIDTH_Z; ++k) for (int32_t j=0; j<block::WIDTH_Y; ++j) for (int32_t i=0; i<block::WIDTH_X; ++i) {
	 // Calculate cell (i,j,k) indices:
	 const uint32_t i_cell = i_block*block::WIDTH_X + i;
	 const uint32_t j_cell = j_block*block::WIDTH_Y + j;
	 const uint32_t k_cell = k_block*block::WIDTH_Z + k;
	 
	 // Get cell volume:
	 const Real volume = simControl.cellVolumes[blockLID*block::SIZE+block::index(i,j,k)];
	 
	 // Adjust N_particlesPerCell so that ~constant number of macroparticles
	 // is injected per second:
	 uint32_t N_integer;
	 if (adjustInjectionRate == false) {
	    N_integer = N_particlesPerCell*10;
	 } else {
	    // Calculate number of injected macroparticles:
	    Real N_real = N_particlesPerCell * sim->dt;
	    N_integer = std::max((uint32_t)10,static_cast<uint32_t>(floor(N_real)));

	    // Randomly determine if an additional particle is injected:
	    Real fraction = N_real - N_integer;
	    if (simClasses->random.uniform() < fraction) ++N_integer;
	 }

	 // Inject particles to cell:
	 injEnergy.N_particles = N_integer;
	 for (uint32_t n=0; n<N_integer; ++n) {
	    injEnergy.N = n;
	    injEnergy.superResolution = 1;

	    // Calculate random injection position in logical coordinates:
	    trial.state[particle::XCRD] = i_cell + simClasses->random.uniform();
	    trial.state[particle::YCRD] = j_cell + simClasses->random.uniform();
	    trial.state[particle::ZCRD] = k_cell + simClasses->random.uniform();

	    // Do not inject particles to downstream of shock:
	    //if (simControl.shock->getSquaredDistanceToShock(sim->t+sim->dt,trial.state) < 0) continue;
	    if (simControl.shock->getShockRegion(t_inj,trial.state) != 1) continue;

	    // Get injection frame speed:
	    Real V_frame[3] = {0,0,0};
	    Real V_wave[3];
	    Real V_alfven;
	    Real B[3];
	    getFrameSpeed(blockLID,t_inj,trial.state,B,V_frame,V_wave,V_alfven);
	    //const Real V_wave_mag = vectorMagnitude<3>(V_wave);

	    // Get plasma state at injection position:
	    (*simControl.fieldsGetPlasmaState)(blockLID,sim->t,trial.state,plasmaState);
	    injEnergy.thermalEnergy = constants::BOLTZMANN*plasmaState.ionTemperature;

	    /*
	    // Increase superresolution if wave generation is turned on 
	    // and macroparticle weight is too large:
	    if (checkWaveEnergies == true) {
	       const Real omega = species.q_per_m * vectorMagnitude<3>(B);

	       generateParticle(trial,V_frame,injEnergy,plasmaState,volume);
	       origEnergy = injEnergy.energy / constants::CHARGE_ELEMENTARY / 1e6;
	       
	       Real V_mag = sqrt(2*injEnergy.energy/species.mass);

	       Real maxResonantWavelengthPhys = -2*M_PI/omega*V_mag;
	       Real maxResonantWavelengthLogi = (*simControl.getLogicalWavelength)(maxResonantWavelengthPhys);
	       const int32_t l_index = static_cast<int32_t>(maxResonantWavelengthLogi);
	       
	       Real dU = 0.5*densityMultiplier*trial.state[particle::WEIGHT]*species.mass*omega/M_PI*V_wave_mag;
	       if (dU/minWaveEnergies[l_index] > 0.025) {
		  injEnergy.superResolution = std::max(1,static_cast<int32_t>(10*dU/minWaveEnergies[l_index]));
	       }
	    }*/
	    
	    // Generate injEnergy.superResolution particles and add them to particle list:
	    for (int32_t s=0; s<injEnergy.superResolution; ++s) {
	       injEnergy.N_super = s;
	       generateParticle(trial,V_frame,injEnergy,plasmaState,volume);
	       wrapper.push_back(blockLID,trial);
	       ++N_particles[blockLID];
	    }
	 }
      }
      
      // Measure computation time if we are testing for repartitioning:
      if (sim->countPropagTime == true) {
	 t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	 simClasses->pargrid.getCellWeights()[blockLID] += t_propag;
      }
   }

   template<class SPECIES,class PARTICLE> inline
   void ParticleInjectorShock<SPECIES,PARTICLE>::injectParticles(pargrid::CellID blockLID,const Real& t_inj,
								 const Real* const logicalCentroid,
								 const Real* const shockNormal,const Real* const V_surface,
								 const Real& area,
								 InjectionBuffer<PARTICLE>* injBuffer,
								 int direction) {
      // Get current index to injBuffer:
      //const size_t startIndex = injBuffer->particles.size();

      bool checkWaveEnergies = false;
      std::vector<Real> minWaveEnergies;
      if (simControl.applyWaveGrowth == true) {
	 checkWaveEnergies = getMinimumWaveEnergies(blockLID,+1,minWaveEnergies);
      }
      checkWaveEnergies = false;

      // Get plasma state at injection position:
      PlasmaState plasmaState;
      (*simControl.fieldsGetPlasmaState)(blockLID,t_inj,logicalCentroid,plasmaState);

      // Set thermal energy for injection energy distribution:
      distrib::InjectionEnergy injEnergy;
      injEnergy.thermalEnergy = constants::BOLTZMANN*plasmaState.ionTemperature;

      const Real B_mag = vectorMagnitude<3>(plasmaState.B);
      Real B_unit[3];
      for (int i=0; i<3; ++i) B_unit[i] = plasmaState.B[i];
      unitVector<3>(B_unit);

      // Make a dummy particle:
      PARTICLE particle;
      for (int i=0; i<3; ++i) particle.state[i] = logicalCentroid[i];

      // Parallel component of transformation velocity from simulation (SIM)
      // frame to shock normal incidence frame (SNIF):
      Real V_trans_par = 0.0;
      PlasmaParameters plasmaParameters;
      Real d_shock = simControl.shock->getSquaredDistanceToShock(t_inj,particle.state);
      if (injectReflectedOnly == true) {
	 // Get shock velocity in SIM:
	 Real V_shock_SIM[3];
	 simControl.shock->getLocalShockVelocity(t_inj,particle.state,V_shock_SIM);
	 
	 // Transformation velocity:
	 Real V_SIM_2_SNIF[3];
	 getTransformSimToLocalSNIF(plasmaState.V_plasma_SIM,V_shock_SIM,shockNormal,V_SIM_2_SNIF);

	 // Calculate parallel component of V_SIM_2_SNIF:
	 V_trans_par = dotProduct<3>(V_SIM_2_SNIF,B_unit);

	 // Plasma speed in SNIF:
	 Real V_plasma_SNIF[3];
	 for (int i=0; i<3; ++i) V_plasma_SNIF[i] = plasmaState.V_plasma_SIM[i] - V_SIM_2_SNIF[i];

	 // Calculate tangential component of B:
	 Real B_tang[3];
	 const Real B_norm = dotProduct<3>(plasmaState.B,shockNormal);
	 for (int i=0; i<3; ++i) B_tang[i] = plasmaState.B[i] - B_norm*shockNormal[i];

	 // Set plasma parameters for shock accelerator:
	 plasmaParameters.B1_norm        = -B_norm;
	 plasmaParameters.B1_tang        = vectorMagnitude<3>(B_tang);
	 plasmaParameters.V1_plasma_norm = -dotProduct<3>(V_plasma_SNIF,shockNormal);
	 plasmaParameters.V1_plasma_tang = 0.0;
	 plasmaParameters.R_gas          = simControl.shock->getGasCompressionRatio(blockLID,t_inj,particle.state);
	 plasmaParameters.R_magn         = simControl.shock->getMagneticCompressionRatio(blockLID,t_inj,particle.state);
	 plasmaParameters.L_shock        = 1000.0;
      }

      // Simulation time is dynamically adjusted, so we want to inject ~constant 
      // number of macroparticles per second. Adjust injected number of particles:
      Real N_real = N_particlesPerCell * sim->dt;
      uint32_t N_integer = static_cast<uint32_t>(floor(N_real));
      Real fraction = N_real - N_integer;
      if (simClasses->random.uniform() < fraction) ++N_integer;
      
      Real normalization = 1.0;
      normalization = (1.0*N_particlesPerCell)/N_integer;
      /*if (N_integer < 10) {
	 N_integer = 10;
      }*/
      //injEnergy.N_particles = N_integer;
      injEnergy.N_particles = N_particlesPerCell;

      // Get injection frame speed:
      Real V_frame[3] = {0,0,0};
      Real V_wave[3];
      Real V_alfven;
      Real B[3];
      getFrameSpeed(blockLID,t_inj,particle.state,B,V_frame,V_wave,V_alfven);

      // Normalize particle weights to be relative to ambient plasma number density:
      normalization *= densityMultiplier * (plasmaState.ionMassDensity/species.mass) * area * sim->dt;

      for (uint32_t p=0; p<N_integer; ++p) {
	 //injEnergy.N = p;
	 injEnergy.N = N_particlesPerCell*simClasses->random.uniform();
	 injEnergy.superResolution = 1;

	 if (checkWaveEnergies == true) {
	    const Real V_wave_mag = vectorMagnitude<3>(V_wave);
	    const Real omega = species.q_per_m * vectorMagnitude<3>(B);

	    // Get (random) injection energy and pitch:
	    (*getEnergy)(injEnergy,energyParams);
	    const Real V_mag = sqrt(2*injEnergy.energy/species.mass);

	    // Calculate injection velocity in SIM frame:
	    Real V_injection_SIM[3];
	    for (int i=0; i<3; ++i) V_injection_SIM[i] = V_mag*B_unit[i] + V_frame[i];

	    // Injection velocity in shock rest frame:
	    Real V_injection_SRF[3];
	    for (int i=0; i<3; ++i) V_injection_SRF[i] = V_injection_SIM[i] - V_surface[i];
	    const Real V_inj_SRF_n = dotProduct<3>(V_injection_SRF,shockNormal);

	    // Calculate flux-weighted macroparticle weight:
	    particle.state[particle::WEIGHT] = normalization * injEnergy.weight * fabs(V_inj_SRF_n);

	    Real maxResonantWavelengthPhys = -2*M_PI/omega*V_mag;
	    Real maxResonantWavelengthLogi = (*simControl.getLogicalWavelength)(maxResonantWavelengthPhys);
	    const int32_t l_index = static_cast<int32_t>(maxResonantWavelengthLogi);
	    Real dU = 0.5*densityMultiplier*particle.state[particle::WEIGHT]*species.mass*omega/M_PI*V_wave_mag;

	    if (dU/minWaveEnergies[l_index] > 0.025) {
	       injEnergy.superResolution = std::max(1,static_cast<int32_t>(10*dU/minWaveEnergies[l_index]));

	       if (minWaveEnergies[l_index] < 1e9) {
		  pargrid::CellID blockGID = simClasses->pargrid.getGlobalIDs()[blockLID];
		  uint32_t i_block,j_block,k_block;
		  block::calculateBlockIndices(*sim,blockGID,i_block,j_block,k_block);
	       }
	    }
	 }

	 Real multiplier = 1;
	 for (int32_t s=0; s<injEnergy.superResolution; ++s) {
	    injEnergy.N_super = s;

	    // Get (random) injection energy and pitch:
	    (*getEnergy)(injEnergy,energyParams);
	    //std::cerr << N_integer << '\t' << injEnergy.weight << std::endl;
	    Real pitch = (*getPitch)();
	    const Real V_mag = sqrt(2*injEnergy.energy/species.mass);

	    // Calculate injection velocity in SIM frame:
	    Real V_injection_SIM[3];
	    for (int i=0; i<3; ++i) V_injection_SIM[i] = V_mag*pitch*B_unit[i] + V_frame[i];

	    // Injection velocity in shock rest frame:
	    Real V_injection_SRF[3];
	    for (int i=0; i<3; ++i) V_injection_SRF[i] = V_injection_SIM[i] - V_surface[i];
	    const Real V_inj_SRF_n = dotProduct<3>(V_injection_SRF,shockNormal);
	    
	    // Calculate parallel speed and magnetic moment:
	    const Real V_gyro2 = V_mag*V_mag*(1-pitch*pitch);
	    particle.state[particle::V_PAR] = dotProduct<3>(V_injection_SIM,B_unit);
	    particle.state[particle::MU] = 0.5*species.mass*V_gyro2 / B_mag;

	    if (injectReflectedOnly == true) {
	       // Set particle parameters for shock accelerator:
	       ParticleParameters p;
	       if (d_shock >= 0.0) p.state[shockaccelerator::XPOS] = -0.5*plasmaParameters.L_shock;
	       else p.state[shockaccelerator::XPOS]  = +0.5*plasmaParameters.L_shock;
	       p.state[shockaccelerator::YPOS]  = 0.0;
	       p.state[shockaccelerator::ZPOS]  = 0.0;
	       p.state[shockaccelerator::V_PAR] = particle.state[sep::particle::V_PAR] - V_trans_par;
	       p.state[shockaccelerator::MU]    = particle.state[sep::particle::MU];
	       
	       // Solve shock encounter:
	       simControl.shockAccelerator.setSpecies(species);
	       multiplier = simControl.shockAccelerator.getReturnedParticle(plasmaParameters,p,V_alfven,pitch);
	       if (multiplier < 0) continue;

	       // Reject transmitted particles:
	       if (p.state[shockaccelerator::XPOS] > 0.0) continue;

	       // Particle reflected, copy new parallel speed and magnetic moment:
	       particle.state[sep::particle::V_PAR] = p.state[shockaccelerator::V_PAR] + V_trans_par;
	       particle.state[sep::particle::MU]    = p.state[shockaccelerator::MU];
	    } else if (V_inj_SRF_n*direction < 0) {
	       // Reject particle if V_injection_SRF does not point towards the shock:
	       continue;
	    }

	    // Calculate flux-weighted macroparticle weight:
	    particle.state[particle::WEIGHT] = multiplier * normalization * injEnergy.weight * fabs(V_inj_SRF_n);

	    // Insert particle to injection buffer:
	    injBuffer->insert(blockLID,particle);
	 }
      }
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleInjectorShock<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
								   const std::string& regionName,const ParticleListBase* plist) {
      simClasses.logger << "(SEP PARTICLE INJ SHOCK) Starting to init." << std::endl << write;

      // Init base class:
      initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,regionName,plist);
      this->species = *reinterpret_cast<const SPECIES*>(plist->getSpecies());

      // Parse config file options:
      //std::string propagatorName,energyPerNucleonString;
      std::string energyDistribString,energyDistribParamsString;
      std::string pitchDistribString,pitchDistribParamsString;
      std::string injectionFrameString,injectionRegionString,inflowBoundaryString;
      std::string injectReflectedString;
      addConfigFileItems(cr,regionName);
      cr.parse();
      cr.get(regionName+".macroparticles_per_cell",N_particlesPerCell);
      cr.get(regionName+".energy_distribution",energyDistribString);
      cr.get(regionName+".energy_distribution_parameters",energyDistribParamsString);
      cr.get(regionName+".pitch_distribution",pitchDistribString);
      cr.get(regionName+".pitch_distribution_parameters",pitchDistribParamsString);
      cr.get(regionName+".relative_abundance",relativeAbundance);
      cr.get(regionName+".injection_frame",injectionFrameString);
      cr.get(regionName+".injection_region",injectionRegionString);
      cr.get(regionName+".inflow_boundary",inflowBoundaryString);
      cr.get(regionName+".density_multiplier",densityMultiplier);
      cr.get(regionName+".inject_reflected_only",injectReflectedString);
      cr.get(regionName+".surface_normal_x",surfaceNormal_x);
      cr.get(regionName+".injection_time_min",t_injection_min);
/*
      energyPerNucleon = false;
      if (energyPerNucleonString == "yes") energyPerNucleon = true;
*/      
      if (relativeAbundance <= 0.0) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Relative abundance must be positive" << std::endl << write;
	 initialized = false;
      }

      // Determine injection frame and region:
      injectionFrame = getInjectionFrame(injectionFrameString);
      if (injectionFrame == sep::particleinjector::frame::UNKNOWN) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown injection frame" << std::endl << write;
	 initialized = false;
      }
      
      injectionRegion = getInjectionRegion(injectionRegionString);
      if (injectionRegion == sep::particleinjector::region::UNKNOWN) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Unknown injection region" << std::endl << write;
	 initialized = false;
      }
      
      // Check if only reflected particles are accepted:
      if (injectReflectedString == "yes") injectReflectedOnly = true;
      else injectReflectedOnly = false;

      // Attempt to get energy distribution function:
      sep::initializeEnergyDistrib initializeEnergy;
      sep::getEnergyDistribFunction getDistrib;
      if (sep::getObjectWrapper().energyDistribContainer.getDistribution(energyDistribString,finalizeEnergy,
									 getDistrib,getEnergy,initializeEnergy) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Could not find energy distribution ";
	 simClasses.logger << "function called '";
	 simClasses.logger << energyDistribString << "'," << std::endl;
	 simClasses.logger << "\t given with parameter '" << regionName+".energy_distribution'" << std::endl << write;
	 initialized = false;
      }
      if (initialized == true) if (initializeEnergy(sim,simClasses,cr,energyDistribParamsString,energyParams) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Energy distribution function failed ";
	 simClasses.logger << "to initialize" << std::endl << write;
	 initialized = false;
	 finalizeEnergy(energyParams);
      }

      // Attempt to get pitch distribution function and initialize it:
      sep::initializePitchDistrib initializePitch;
      if (sep::getObjectWrapper().pitchDistribContainer.getDistribution(pitchDistribString,finalizePitch,getPitch,initializePitch) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Could not find pitch distribution function called '";
	 simClasses.logger << pitchDistribString << "'," << std::endl << write;
	 simClasses.logger << "\t given with parameter '" << regionName+".pitch_distribution'" << std::endl << write;
	 initialized = false;
      }
      if (initialized == true) if (initializePitch(sim,simClasses,cr,pitchDistribParamsString) == false) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Pitch distribution function failed to initialize" << std::endl << write;
	 initialized = false;
      }

      // Check parsed values:
      if (N_particlesPerCell == 0) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Parameter '" << regionName+".macroparticles_per_cell' was not found" << std::endl << write;
	 initialized = false;
      }
      
      // Determine inflow boundary (if any):
      if (inflowBoundaryString.size() == 0) {
	 simClasses.logger << "(SEP PARTICLE INJ SHOCK) ERROR: Inflow boundary was not specified" << std::endl << write;
	 initialized = false;
      }
      inflowBoundaryMask = 0;
      if (inflowBoundaryString == "-x") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0));
      if (inflowBoundaryString == "+x") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0));
      if (inflowBoundaryString == "-y") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0));
      if (inflowBoundaryString == "+y") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0));
      if (inflowBoundaryString == "-z") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1));
      if (inflowBoundaryString == "+z") inflowBoundaryMask = (inflowBoundaryMask | 1 << simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1));
       
      // Write out injection parameters:
      simClasses.logger << "(SEP PARTICLE INJ SHOCK)  Initialization status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << std::endl << write;
      else simClasses.logger << "FAILURE" << std::endl << write;
      
      // Set some internal variables to correct values if simulation was restarted:
      if (sim.restarted == true) {
	 initialInjectionDone = true;
      }

      //this->plist = plist;
      return initialized;
   }

} // namespace sep
   
#endif
