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

#include <main.h>
#include "sep_register_dataoperators.h"

#include <operator_mpirank.h>
#include <operator_cellid.h>
#include <operator_particle.h>
#include <operator_load.h>
#include "sep_operator_anisotropy.h"
#include "sep_operator_energy_channels.h"
#include "sep_operator_differential_flux.h"
#include "sep_operator_Dmumu.h"
#include "sep_operator_mfp.h"
#include "sep_operator_plasma_state.h"
#include "sep_operator_reduced_scalar.h"
#include "sep_operator_shock_electrons.h"
#include "sep_operator_shock_mesh.h"
#include "sep_operator_shock_variables.h"
#include "sep_operator_spherical_surface_mesh.h"
#include "sep_operator_wavelength_mesh.h"
#include "sep_operator_wavelength_mesh_density.h"
#include "sep_operator_wavelength_mesh_scalar.h"
#include "sep_operator_wavelength_mesh_particles.h"
#include "sep_operator_wavelength_mesh_volume.h"
#include "sep_operator_density.h"
#include "sep_operator_particle_energy.h"
#include "sep_operator_spatial_lineout_energy_spectrum.h"
#include "sep_operator_timeseries_energy.h"
#include "sep_operator_wave_cross_helicity.h"
#include "sep_operator_wave_energy_total.h"
#include "sep_operator_volume.h"

using namespace std;

namespace sep {

   /** Wrapper function for speeding up compilation. Function registers 
    * DataOperators to DataOperatorContainer.
    * @return If true, DataOperators were registered successfully.*/
   bool registerDataOperators() {
      bool success = true;
      DataOperatorContainer& doc = corsair::getObjectWrapper().dataOperatorContainer;
      if (doc.registerOperator(new LoadOP()) == false) success = false;
      if (doc.registerOperator(new MPIRank()) == false) success = false;
      if (doc.registerOperator(new CellIDOP()) == false) success = false;
      if (doc.registerOperator(new ParticleOperator) == false) success = false;
      if (doc.registerOperator(new sep::OperatorAnisotropy(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorEnergyChannels(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorDifferentialFlux(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorDmumu) == false) success = false;
      if (doc.registerOperator(new sep::OperatorMeanFreePath(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorPlasmaState()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorReducedScalar()) == false) success = false;
      if (doc.registerOperator(new sep::WavelengthMeshOP()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorDensity(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorParticleEnergy(2)) == false) success = false;
      if (doc.registerOperator(new sep::OperatorWaveCrossHelicity()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorWaveEnergyTotal()) == false) success = false;
      if (doc.registerOperator(new sep::WavelengthMeshDensityOP()) == false) success = false;
      if (doc.registerOperator(new sep::WavelengthMeshParticles()) == false) success = false;      
      if (doc.registerOperator(new sep::WavelengthMeshScalarOP()) == false) success = false;
      if (doc.registerOperator(new sep::WavelengthMeshVolumeOP()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorTimeseriesEnergy()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorVolume()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorShockElectrons()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorShockMesh()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorShockVariables()) == false) success = false;
      if (doc.registerOperator(new sep::OperatorSphericalSurfaceMesh()) == false) success = false;
      if (doc.registerOperator(new sep::SpatialLineoutEnergySpectrum()) == false) success = false;
      return success;
   }
   
} // namespace sep
