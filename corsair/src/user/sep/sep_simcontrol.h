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

#ifndef SEP_SIMCONTROL_H
#define SEP_SIMCONTROL_H

#include <definitions.h>
#include <pargrid_definitions.h>
#include <main.h>

#include "sep_fields_container.h"
#include "sep_base_class_shock.h"
#include "sep_shock_accelerator.h"
#include "sep_base_class_wavelength_scaling.h"
#include "sep_mesh_logical.h"
#include "sep_wavelength_mesh_builder.h"
#include "sep_base_class_particle_splitter.h"
#include "sep_distrib_wave_energy_base_class.h"

namespace sep {
   // TODO: Move these into nested namespace
   enum CoordinateSystem {
      UNKNOWN,                      /**< Unknown or unsupported coordinate system.*/
      CARTESIAN,                    /**< Simulation is run using Cartesian coordinates.*/
      CYLINDRICAL,                  /**< Simulation is run using cylindrical coordinates.*/
      SPHERICAL                     /**< Simulation is run using spherical coordinates.*/
   };
   
   struct SimControl {
      SimControl();
      ~SimControl();

      int order;
      
      Real dU_particles;
      Real dU_waves;

      
      const Real EPSILON;
      const Real BLOCK_EPSILON;                          /**< Epsilon added to or substracted from particle 
							  * coordinates when they move across periodic boundaries. 
							  * Used to ensure that coordinates are valid to floating point precision.*/
      CoordinateSystem coordinateSystem;                 /**< Coordinate system in use.*/
      std::string fieldFunctionName;                     /**< Name of EM field function.*/
      Real* cellVolumes;                                 /**< Spatial cell volumes.*/

      // ***** Variables controlling execution of simulation code ***** //
      bool propagateAlfvenWaves;                         /**< If true, Alfven waves are propagated.*/
      bool scatterParticles;                             /**< If true, particles scatter off from Alfven waves. Setting this to
							  * 'true' in log file also sets variables SimControl::scatterAntiparallel
							  * and SimControl::scatterParallel to values 'true' by default.*/
      bool scatterAntiparallel;                          /**< If true, particles scatter off antiparallel waves. Defaults to
							  * the value of SimControl::scatterParticles, but can be turned off 
							  * separately in config file.*/
      bool scatterParallel;                              /**< If true, particles scatter off parallel waves. Defaults to 
							  * the value of SimControl::scatterParticles, but can be turned off 
							  * separately in config file.*/
      bool applyWaveGrowth;                              /**< If true, wave growth factors are applied to Alfven waves.*/
      bool applyAntiparallelWaveGrowth;                  /**< If true, wave growth factors are applied to antiparallel waves. 
							  * Defaults to the value of SimControl::applyWaveGrowth, but can be
							  * separately turned off in config file.*/
      bool applyParallelWaveGrowth;                      /**< If true, wave growth factors are applied to parallel waves. 
							  * Defaults to the value of SimControl::applyWaveGrowth, but can be
							  * separately turned off in config file.*/
      bool includeAntiparWaves;                          /**< If true, antiparallel propagating Alfven waves are
							  * included in simulation.*/
      bool includeParWaves;                              /**< If true, parallel propagating Alfven waves are 
							  * included in simulation.*/
      bool includeShock;                                 /**< If true, a shock wave is included in simulation.*/

      Real courantNumber;                                /**< Courant-Friedrich-Levy number used to calculate simulation time step.*/
      Real R_reference;                                  /**< Reference value of radius, used to scale resonant wavelengths.
							  * Currently this value is set in sep_fields_* initializers (DEPRECATED).*/
      Real t_setup;                                      /**< Particles are not injected into simulation until t >= t_setup.
							  * This allows wave background to settle in.*/
      Real t_shock;                                      /**< Shock is not turned on until t >= t_shock, this allows wave and 
							  * particle background to set up.*/
      Real dt_maximum;                                   /**< Maximum allowed timestep, negative value means that any 
							  * positive dt is accepted.*/
      uint32_t timestepRecalculateInterval;              /**< Interval, in number of timesteps, in which dt is re-evaluated.*/
      int32_t N_lagrangianSpecies;                       /**< Total number of propagated Lagrangian species.*/
      int32_t N_particleSpecies;                         /**< Total number of propagated particle species.*/
      Real dt_particle_coeff;                            /**< Maximum acceptable time step, calculated from particles' 
							  * Courant condition, is multiplied by this factor. This is mainly 
							  * used to take shock drift acceleration into account which can 
							  * increase parallel speeds by a large factor.*/

      ParticleSplitterBase* particleSplitter;              /**< Particle Splitter.*/

      // ***** Variables related to spatial mesh ***** //
      int32_t timestepCellVolumesCalculated;               /**< Timestep when cell volumes in cellVolumes array were calculated.*/
      Real V_frame[3];                                     /**< Constant frame velocity, defaults to 0.*/

      // ***** Variables related to wavelength mesh ***** //
      Real lambda0;                                        /**< Wavelength mesh is linear instead of logarithmic below this wavelength.*/
      int32_t alfvenSign;                                  /**< Sign of the current alfven wave mode used in scattering, 
							    * positive for parallel and negative for antiparallel.*/
      
      getLogicalWavelengthFunction getLogicalWavelength;   /**< Pointer to function that returns the logical wavelength
							    * corresponding to a given physical wavelength.*/
      getPhysicalWavelengthFunction getPhysicalWavelength; /**< Pointer to function that returns the physical wavelength 
							    * corresponding to a given logical wavelength.*/
      getWavelengthScaleFactorFunction getWavelengthScaleFactor;
      
      WaveEnergySpectrumBaseClass* parWaveEnergy;
      WaveEnergySpectrumBaseClass* antiparWaveEnergy;
      
      Real logicalWavelengthCellSize;
      mesh::Type wavelengthMeshType;                     /**< Mesh type (linear, logarithmic, etc.).*/
      uint32_t N_wavelengthMeshCells;                    /**< Number of cells in wavelength meshes.*/
      Real* wavelengthMeshCellSizes;                     /**< Wavelength mesh cell sizes in physical units.*/
      Real* wavelengthMeshNodeCoordinates;               /**< Wavelength mesh node coordinates. Size of array
							  * is N_wavelengthMeshCells+1. Size of cell i is
							  * the difference between node i+1 and i coordinates.*/
      Real maxScatteringDeltaMu;                         /**< Maximum change in pitch in scattering, controls 
							  * the number of scattering substeps.*/
      uint32_t currentMaxScatteringSubsteps;             /**< Current maximum number of scattering substeps taken, 
							  * calculated over all cells on this process.*/
      Real maximumWaveEnergyLossPercentage;              /**< Wave packet can only lose this many percentages of its
							  * current energy during one time step.*/
      Real bohmLimitCoefficient;                         /**< Wave energy spectra are limited to this coefficient times 
							  * Bohm spectrum.*/

      // ***** Variables related to shock waves ***** //
      uint32_t N_shockSurfaceRefinements;                /**< How many times shock surface elements are refined.*/
      bool useShockUpwinding;                            /**< If true, wave packets and particles use modified shape 
							  * factors when near the shock front to upwind/downwind interpolations.*/
      bool useUpwinding;                                 /**< If true, accumulations and interpolations are upwinded or 
							  * downwinded, as appropriate, to reduce diffusion. Difference to 
							  * SimControl::useShockUpwinding is that up/downwinding is now done 
							  * everywhere.*/
      ShockBaseClass* shock;                             /**< Pointer to shock wave class.*/
      ShockAccelerator shockAccelerator;                 /**< Function that calculates particles' energy 
							  * changes in shock crossings.*/
      
      // ***** Pointers to functions that return values of physical parameters ***** //
      finalizeField fieldsFinalize;
      getFields fieldsGetFields;
      getPlasmaStateFunction fieldsGetPlasmaState;
      getStateFunction fieldsGetState;
      initializeField fieldsInitialize;
      
      // ***** ParGrid array names ***** //
      std::string particleWeightName;                     /**<                DEPRECATED                            */
      std::string antiparAlfvenName;                      /**< Name of ParGrid array that contains antiparallel spectral wave energy.
							   * Array's ParGrid ID is stored in SimControl::antiparAlfvenWaveEnergyDataID.*/
      std::string antiparWaveGrowthName;
      std::string parAlfvenName;                          /**< Name of ParGrid array that contains parallel spectral wave energy.
							   * Array's ParGrid ID is stored in SimControl::parAlfvenWaveEnergyDataID.*/
      std::string parWaveGrowthName;
      std::string waveGrowthName;                         /**<                DEPRECATED                            */
      std::string lagrangianSpeciesTypename;              /**< Species typename that identifies that macroparticles
							   * are Alfven wave packets instead of particles.*/
      std::string particleSpeciesTypename;                /**< Species typename that identifies that macroparticles 
							   * are particles instead of Alfven wave packets.*/
      std::string temporaryWaveEnergyName;                /**<                DEPRECATED                            */
      std::string maxAntiparValidWavelengthBinsName;
      std::string maxParValidWavelengthBinsName;
      
      // ***** ParGrid DataID:s for arrays used in SEP simulation ***** //

      pargrid::DataID waveGrowthDataID;               /**< Wave growth factor, DEPRECATED.*/
      pargrid::DataID particleWeightDataID;           /**< Sum of particle weights in 4D (x,y,z,lambda) phase space, DEPRECATED.*/
      
      pargrid::DataID antiparAlfvenDataID;            /**< Antiparallel propagating Alfven waves, particle ID.*/
      pargrid::DataID antiparAlfvenWaveEnergyDataID;  /**< Antiparallel propagating Alfven waves, energy.*/
      pargrid::DataID antiparWaveGrowthDataID;
      pargrid::DataID parAlfvenDataID;                /**< Parallel propagating Alfven waves, particle ID.*/
      pargrid::DataID parAlfvenWaveEnergyDataID;      /**< Parallel propagating Alfven waves, energy.*/
      pargrid::DataID parWaveGrowthDataID;
      pargrid::DataID temporaryWaveEnergyDataID;      /**<                     DEPRECATED                       */
      pargrid::DataID maxAntiparValidWavelengthBinsDataID;
      pargrid::DataID maxParValidWavelengthBinsDataID;
   };

} // namespace sep

#endif
