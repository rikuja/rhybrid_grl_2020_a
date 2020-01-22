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

#include <main.h>
#include <object_factory_generic.h>
#include <linear_algebra.h>
#include "sep_object_wrapper.h"
#include "sep_simcontrol.h"
#include "sep_particle_definition.h"
#include "sep_particle_species.h"
#include "sep_lagr_species.h"
#include "sep_particle_scatterer.h"
#include "sep_accumulation_stretched.h"
#include "sep_distrib_wave_energy_base_class.h"
#include "sep_wavelength_mesh_builder.h"

#include "sep_shock_planar.h"

using namespace std;

namespace sep {

   // TEST
   Simulation* sim = NULL;
   // END TEST

   struct CellInfo {
      int i_block;
      int j_block;
      int k_block;
   };
   
   CellInfo cellInfo;
   
   extern sep::SimControl simControl;
   
   //static const Real mu_limit = 0.025;
   static const Real mu_limit = 0.25;
   //static const Real RESONANCE_GAP_CONST_2 = 0.001;

   static Real limitedPitchMin = 0.0;
   static Real limitedPitchMax = 0.0;
   static Real Q_per_M = -1.0;          /**< Charge to mass ratio of particle species currently 
					 * being scattered. Stored here so that evaluateNumberOfSubsteps 
					 * finds it.*/

   typedef sep::Particle<Real> PARTICLE;

   static uint32_t scatteringSubsteps = 1;

   typedef Real (*scatterer)(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			     int32_t* indices,Real* shapeFactors,
			     Real* RESTRICT pos,const Real& mu,const Real& omega,
			     const Real& speed,const Real& B_mag,const Real& minMu,const Real& maxMu);

   // Pointer to scatterer that is used in simulation:
   scatterer particleScatterer = NULL;

   /** This struct is only used when plotting output from 
    *runTestDiffusionCoefficient test.*/
   struct DiffCoeff {
      vector<Real> mu;
      Real** lambda;
      Real** D_mumu_analytic;
      Real** D_mumu_interpolated;
      Real** d_D_mumu_dmu_analytic;
      Real** d_D_mumu_dmu_interpolated;
      Real** d_D_mumu_dmu_diff;
      
      DiffCoeff();
      DiffCoeff(int N_points,int arrays);
      ~DiffCoeff();
   };
   
   DiffCoeff::DiffCoeff() {
      lambda = NULL;
      D_mumu_analytic = NULL;
      D_mumu_interpolated = NULL;
      d_D_mumu_dmu_analytic = NULL;
      d_D_mumu_dmu_interpolated = NULL;
      d_D_mumu_dmu_diff = NULL;
   }
   
   DiffCoeff::DiffCoeff(int N_points,int arrays) {
      Real muMin = -1.0, muMax = +1.0;
      Real d_mu = (muMax-muMin) / (N_points-1);
      
      mu.resize(N_points);
      for (int i=0; i<N_points; ++i) mu[i] = muMin + i*d_mu;
      
      lambda = new Real* [arrays];
      D_mumu_analytic = new Real* [arrays];
      D_mumu_interpolated = new Real* [arrays];
      d_D_mumu_dmu_analytic = new Real* [arrays];
      d_D_mumu_dmu_interpolated = new Real* [arrays];
      d_D_mumu_dmu_diff = new Real* [arrays];
      for (int i=0; i<arrays; ++i) {
	 lambda[i] = new Real[N_points];
	 D_mumu_analytic[i] = new Real[N_points];
	 D_mumu_interpolated[i] = new Real[N_points];
	 d_D_mumu_dmu_analytic[i] = new Real[N_points];
	 d_D_mumu_dmu_interpolated[i] = new Real[N_points];
	 d_D_mumu_dmu_diff[i] = new Real[N_points];
      }
   }
   
   DiffCoeff::~DiffCoeff() {
      delete [] lambda;
      delete [] D_mumu_analytic;
      delete [] D_mumu_interpolated;
      delete [] d_D_mumu_dmu_analytic;
      delete [] d_D_mumu_dmu_interpolated;
      delete [] d_D_mumu_dmu_diff;
   }
   
   #if PROFILE_LEVEL > 0
      int profScatterTotal =-1;
   #endif

   WaveEnergySpectrumBaseClass* getWaveEnergySpectrum(corsair::ObjectWrapper& objectWrapper,int propagationDirection);

   uint32_t evaluateNumberOfSubsteps(Simulation& sim,pargrid::CellID blockLID,
				     const uint32_t* const blockIndices,
				     const wmesh::validDatatype* validBins,
				     const Real* lambdaIntensity,Real maxSpeed);
   
   Real interpolateWaveEnergy(const Real* RESTRICT pos,const Real* RESTRICT waveEnergy);
   Real interpolateWaveEnergyDerivated(const Real* RESTRICT pos,const Real* RESTRICT waveEnergy);
   
   void evaluateDiffusionCoefficients(Real pitch,Real maxResonantLambda,Real* pos,const Real* intensity,Real& D_mumu,Real& d_D_mumu);
   void evaluateDiffusionCoefficientsIsotropic(Real B_mag,Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,
					       int32_t* indices,Real* shapeFactors,
					       const Real* intensity,Real& D_mumu,Real& d_D_mumu);
   void evaluateDiffusionCoefficientsIsotropicAnalytic(Real B_mag,Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,
						       const Real* intensity,Real& D_mumu,Real& d_D_mumu);
   void fetchIntensities(const Real* pos,const Real* lambdaIntensity,Real& I_wave,Real& dI_wave);
   bool loadDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,
				 Real cellVolume,int propagationDirection,Real* RESTRICT waveIntensity);
   bool runTestDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,int propagationDirection);
   bool runTestDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,
				    Real wavelengthConstant,DiffCoeff& diffCoeff,int index,int propagationDirection,Real& wavelength);
   bool runTestScatteringIsotropicity(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,int propagationDirection);
   bool runTestWaveIntensity(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,int propagationDirection);

   Real scatterImplMilstein(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			    Real* RESTRICT pos,const Real& mu,const Real& omega,
			    const Real& speed,const Real& B_mag);
   
   Real scatterMilstein(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			Real* RESTRICT pos,const Real& mu,const Real& omega,
			const Real& speed,const Real& B_mag,const Real& minMu,const Real& maxMu);

   Real scatterPredCorr1(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			 int32_t* indices,Real* shapeFactors,
			 Real* RESTRICT pos,const Real& mu,const Real& omega,
			 const Real& speed,const Real& B_mag,const Real& minMu,const Real& maxMu);
   
   Real scatterPredCorr2(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			 Real* RESTRICT pos,const Real& mu,const Real& omega,
			 const Real& speed,const Real& B_mag,const Real& minMu,const Real& maxMu);
   
   Real scatterIsotropic(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			 int32_t* indices,Real* shapeFactors,
			 Real* RESTRICT pos,const Real& mu,const Real& omega,
			 const Real& speed,const Real& B_mag,const Real& mu_min,const Real& mu_max);
   
   Real scatterIsotropicAnalytic(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
				 int32_t* indices,Real* shapeFactors,
				 Real* RESTRICT pos,const Real& mu,const Real& omega,
				 const Real& speed,const Real& B_mag,const Real& mu_min,const Real& mu_max);

   bool writeIsotropicity(size_t N_particles,size_t N_bins,
			  uint32_t* probabMatrix,const std::string& name,int propagationDirection);
   
   bool writeTransitionProbabilities(size_t N_particles,size_t N_bins,const std::vector<size_t>& bins,
				     uint32_t* probabMatrix,const std::string& name,int propagationDirection);
   
   void accumulateWaveGrowth(const Real& energyChangeFactor,const Real& lambdaMin,const Real& lambdaMax,
			     const int32_t* RESTRICT indices,const Real* RESTRICT shapeFactors,
			     Real* RESTRICT waveGrowthTmp) {
      
      // Calculate average intensity seen by particle:
      const int32_t NWL = simControl.N_wavelengthMeshCells+2;
      const Real d_lambda = lambdaMax-lambdaMin;

      switch (simControl.order) {
       case 0:
	 
	 break;
       case 1:
	 for (int k=0; k<2; ++k) for (int j=0; j<2; ++j) for (int i=0; i<2; ++i) {
	    const Real shapeFactor = shapeFactors[0+i]*shapeFactors[2+j]*shapeFactors[4+k];
	    const int32_t index = block::arrayIndex(indices[0]+i,indices[1]+j,indices[2]+k)*NWL+indices[3];//+l;
	    const Real gamma = energyChangeFactor*d_lambda;
	    waveGrowthTmp[index] += shapeFactor * gamma;
	 }
	 break;
       case 2:
	 for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) {
	    const Real shapeFactor = shapeFactors[1+i_off]*shapeFactors[4+j_off]*shapeFactors[7+k_off];
	    const int32_t index = block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*NWL+indices[3];//+l;
	    const Real gamma = energyChangeFactor*d_lambda;
	    waveGrowthTmp[index] += shapeFactor * gamma;
	 }
	 break;
       default:
	 exit(1);
	 break;
      }
   }

   void accumulateWaveGrowthFactors(const pargrid::CellID& blockLID,SimulationClasses* simClasses,
				    int32_t* RESTRICT indices,const Real* RESTRICT shapeFactors,
				    const Real* pos,Real energyChangeFactor,Real lambda_res_old,
				    Real lambda_res_new,
				    Real* RESTRICT waveGrowthTmp) {
      #ifndef NDEBUG
      bool ok = true;
      const uint32_t MAX_INDEX = simControl.N_wavelengthMeshCells;
      if (lambda_res_old < simControl.wavelengthMeshNodeCoordinates[0]) ok = false;
      if (lambda_res_old > simControl.wavelengthMeshNodeCoordinates[MAX_INDEX]) ok = false;
      if (lambda_res_new < simControl.wavelengthMeshNodeCoordinates[0]) ok = false;
      if (lambda_res_new > simControl.wavelengthMeshNodeCoordinates[MAX_INDEX]) ok = false;
      if (ok == false) {
	 cerr << "(SEP PARTICLE SCATTERER) ERROR: Invalid wavelengths in accumulateWaveGrowthFactors" << endl;
	 cerr << "\t" << lambda_res_old << '\t' << lambda_res_new << endl;
	 exit(1);
      }
      #endif
      
      if (lambda_res_new < lambda_res_old) {
	 energyChangeFactor *= -1.0;
	 const Real tmp = lambda_res_old;
	 lambda_res_old = lambda_res_new;
	 lambda_res_new = tmp;
      }
      
      const uint32_t N = simControl.N_wavelengthMeshCells;
      Real* L_start = lower_bound(simControl.wavelengthMeshNodeCoordinates,simControl.wavelengthMeshNodeCoordinates+N+1,lambda_res_old);
      int L_index = L_start-simControl.wavelengthMeshNodeCoordinates;
      if (*L_start > lambda_res_old) {
	 --L_index;
	 --L_start;
      }
      
      while (simControl.wavelengthMeshNodeCoordinates[L_index+1] < lambda_res_new) {
	 indices[3] = L_index+1;
	 ++L_start;
	 const Real lambdaUpper = *L_start;
	 accumulateWaveGrowth(energyChangeFactor,lambda_res_old,lambdaUpper,indices,shapeFactors,
			      waveGrowthTmp);
	 
	 lambda_res_old = simControl.wavelengthMeshNodeCoordinates[L_index+1];
	 ++L_index;
      }
      
      indices[3] = L_index+1;
      accumulateWaveGrowth(energyChangeFactor,lambda_res_old,lambda_res_new,indices,shapeFactors,
			   waveGrowthTmp);
   }

   #warning DEBUG REMOVE
   void evaluateDiffusionCoefficients(Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,
				      int32_t* indices,Real* shapeFactors,
				      const Real* intensity,Real& D_mumu,Real& d_D_mumu) {
      
      Real L_index;
      int32_t l;
      Real pitchMin,pitchMax,pitchLimited;//,pitchCen;
      //Real lambdaMin,lambdaCen,lambdaMax;
      Real I_neg,I_cen,I_pos;
      const int32_t SIZE = simControl.N_wavelengthMeshCells+2;
      switch (simControl.order) {
       case 0:
	 // Evaluate diffusion coefficient at the middle of bin interval:
	 L_index = (*simControl.getLogicalWavelength)(maxResonantLambda*pitch);
	 pos[3] = 1.0 + L_index;
	 l = static_cast<int32_t>(L_index);
	 pitch = max(-1.0,(*simControl.getPhysicalWavelength)(l+0)/maxResonantLambda);
	 pitch += min(+1.0,(*simControl.getPhysicalWavelength)(l+1)/maxResonantLambda);
	 pitch *= 0.5;

	 evaluateDiffusionCoefficients(pitch,maxResonantLambda,pos,intensity,D_mumu,d_D_mumu);
	 d_D_mumu = 0.0;
	 D_mumu   *= diffConstant;
	 d_D_mumu *= diffConstant;
	 break;
       case 1:
	 pitchLimited = pitch;
	 if (fabs(pitch) < mu_limit) {
	    pitchLimited = limitedPitchMin;
	    if (pitch > 0) pitchLimited = limitedPitchMax;
	 }

	 L_index = (*simControl.getLogicalWavelength)(maxResonantLambda*pitchLimited);
	 l = static_cast<int32_t>(L_index);
	 
	 #ifndef NDEBUG
	 if (l < 0 || l >= simControl.N_wavelengthMeshCells-1) {
	    cerr << "(SEP PARTICLE SCATTERER) ERROR: Invalid logical coordinate" << endl;
	    cerr << "\t wavelength : " << maxResonantLambda*pitchLimited << endl;
	    cerr << "\t L          : " << L_index << endl;
	    cerr << "\t l          : " << l << endl;
	    exit(1);
	 }
	 #endif

	 // Fetch wave intensities from neighboring wavelength cells:
	 I_neg = fetchIntensityCIC_3D(SIZE,l  ,intensity,indices,shapeFactors);
	 I_cen = fetchIntensityCIC_3D(SIZE,l+1,intensity,indices,shapeFactors);
	 I_pos = fetchIntensityCIC_3D(SIZE,l+2,intensity,indices,shapeFactors);

	 // Calculate average intensities at wavelength cell faces:
	 I_neg = 0.5*(I_neg + I_cen);
	 I_pos = 0.5*(I_cen + I_pos);

	 // Calculate diffusion coefficients at wavelength cell faces:
	 pitchMin = max(-1.0,simControl.wavelengthMeshNodeCoordinates[l  ]/maxResonantLambda);
	 pitchMax = min(+1.0,simControl.wavelengthMeshNodeCoordinates[l+1]/maxResonantLambda);
	 I_neg = diffConstant*(1-pitchMin*pitchMin)*fabs(simControl.wavelengthMeshNodeCoordinates[l  ])*I_neg;
	 I_pos = diffConstant*(1-pitchMax*pitchMax)*fabs(simControl.wavelengthMeshNodeCoordinates[l+1])*I_pos;

	 // Interpolate to particle position:
	 d_D_mumu = (I_pos-I_neg)/(pitchMax-pitchMin);
	 D_mumu = max(0.0,I_neg + d_D_mumu * (pitchLimited - pitchMin) + d_D_mumu * (pitch-pitchLimited));
	 break;
       case 2:
	 #warning CIC used instead of TSC in evaluateDiffusionCoefficients

	 pitchLimited = pitch;
	 if (fabs(pitch) < mu_limit) {
	    pitchLimited = limitedPitchMin;
	    if (pitch > 0) pitchLimited = limitedPitchMax;
	 }

	 L_index = (*simControl.getLogicalWavelength)(maxResonantLambda*pitchLimited);
	 l = static_cast<int32_t>(L_index);

	 // Fetch wave intensities from neighboring wavelength cells:
	 I_neg = fetchIntensityTSC_3D(SIZE,l  ,intensity,indices,shapeFactors);
	 I_cen = fetchIntensityTSC_3D(SIZE,l+1,intensity,indices,shapeFactors);
	 I_pos = fetchIntensityTSC_3D(SIZE,l+2,intensity,indices,shapeFactors);

	 // Calculate average intensities at wavelength cell faces:
	 I_neg = 0.5*(I_neg + I_cen);
	 I_pos = 0.5*(I_cen + I_pos);

	 // Calculate diffusion coefficients at wavelength cell faces:
	 pitchMin = max(-1.0,simControl.wavelengthMeshNodeCoordinates[l  ]/maxResonantLambda);
	 pitchMax = min(+1.0,simControl.wavelengthMeshNodeCoordinates[l+1]/maxResonantLambda);
	 I_neg = diffConstant*(1-pitchMin*pitchMin)*fabs(simControl.wavelengthMeshNodeCoordinates[l  ])*I_neg;
	 I_pos = diffConstant*(1-pitchMax*pitchMax)*fabs(simControl.wavelengthMeshNodeCoordinates[l+1])*I_pos;

	 // Interpolate to particle position:
	 d_D_mumu = (I_pos-I_neg)/(pitchMax-pitchMin);
	 D_mumu = max(0.0,I_neg + d_D_mumu * (pitchLimited - pitchMin) + d_D_mumu * (pitch-pitchLimited));

	 break;
      }
   }

   void evaluateDiffusionCoefficientsIsotropic(Real B_mag,Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,
					       int32_t* indices,Real* shapeFactors,
					       const Real* intensity,Real& D_mumu,Real& d_D_mumu) {
      Real L_index;
      int32_t l;
      Real lambdaMin,lambdaMax;
      Real I_neg,I_cen,I_pos;
      
      // Simplified resonance condition, maxResonantLambda 
      // instead of maxResonantLambda*pitchLimited:
      L_index = (*simControl.getLogicalWavelength)(maxResonantLambda);
      l = static_cast<int32_t>(L_index);
      
      const int32_t SIZE = simControl.N_wavelengthMeshCells+2;
      switch (simControl.order) {
       case 0:
	 cerr << "NGP not implemented in scatterIsotropic" << endl; exit(1);
	 break;
       case 1:
	 I_neg = fetchIntensityCIC_3D(SIZE,l  ,intensity,indices,shapeFactors);
	 I_cen = fetchIntensityCIC_3D(SIZE,l+1,intensity,indices,shapeFactors);
	 I_pos = fetchIntensityCIC_3D(SIZE,l+2,intensity,indices,shapeFactors);

	 I_neg = 0.5*(I_neg + I_cen);
	 I_pos = 0.5*(I_cen + I_pos);
	 lambdaMin = (*simControl.getPhysicalWavelength)(l  );
	 lambdaMax = (*simControl.getPhysicalWavelength)(l+1);
	 
	 I_neg = 2*diffConstant*maxResonantLambda*I_neg;
	 I_pos = 2*diffConstant*maxResonantLambda*I_pos;

	 d_D_mumu = (I_pos-I_neg)/(lambdaMax-lambdaMin);
	 D_mumu = max(0.0,I_neg + d_D_mumu * (maxResonantLambda - lambdaMin));
	 break;
       case 2:
	 #warning Using CIC in wavelength direction instead of TSC in scatterIsotropic
	 I_neg = fetchIntensityTSC_3D(SIZE,l  ,intensity,indices,shapeFactors);
	 I_cen = fetchIntensityTSC_3D(SIZE,l+1,intensity,indices,shapeFactors);
	 I_pos = fetchIntensityTSC_3D(SIZE,l+2,intensity,indices,shapeFactors);
	 
	 I_neg = 0.5*(I_neg + I_cen);
	 I_pos = 0.5*(I_cen + I_pos);
	 lambdaMin = (*simControl.getPhysicalWavelength)(l  );
	 lambdaMax = (*simControl.getPhysicalWavelength)(l+1);
	 
	 I_neg = 2*diffConstant*maxResonantLambda*I_neg;
	 I_pos = 2*diffConstant*maxResonantLambda*I_pos;
	 
	 d_D_mumu = (I_pos-I_neg)/(lambdaMax-lambdaMin);
	 D_mumu = max(0.0,I_neg + d_D_mumu * (maxResonantLambda - lambdaMin));
	 break;
      }
   }
   
   void evaluateDiffusionCoefficientsIsotropicAnalytic(Real B_mag,Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,
						       const Real* intensity,Real& D_mumu,Real& d_D_mumu) {
      if (sim == NULL) {
	 sim = &(corsair::getObjectWrapper().sim);
      }

      // Get wave intensity at maximum resonant wavelength:
      Real I_wave = simControl.parWaveEnergy->getIntensity(B_mag,maxResonantLambda);
      if (simControl.shock != NULL) {
	 Real logicalPos[3];
	 logicalPos[0] = cellInfo.i_block*block::WIDTH_X + pos[0] - 1;
	 logicalPos[1] = cellInfo.j_block*block::WIDTH_Y + pos[1] - 1;
	 logicalPos[2] = cellInfo.k_block*block::WIDTH_Z + pos[2] - 1;

	 const Real factor = 20.55;
	 if (simControl.shock->getShockRegion(sim->t,logicalPos) != 1) {
	    I_wave *= factor;
	 }
      }
      
      D_mumu = 2*diffConstant * maxResonantLambda * I_wave;
   }

   void evaluateDiffusionCoefficients(Real pitch,Real maxResonantLambda,Real* pos,
				      const Real* intensity,Real& D_mumu,Real& d_D_mumu) {
      if (fabs(pitch) < mu_limit) {
	 Real limitedPitch = limitedPitchMin;
	 if (pitch > 0.0) limitedPitch = limitedPitchMax;

	 Real resonantLambda = limitedPitch * maxResonantLambda;
	 pos[3] = 1.0 + (*simControl.getLogicalWavelength)(resonantLambda);

	 Real I_wave,dI_wave;
	 fetchIntensities(pos,intensity,I_wave,dI_wave);
	 dI_wave /= fabs(resonantLambda)*simControl.logicalWavelengthCellSize;

	 evaluateDiffusionCoefficients(limitedPitch,maxResonantLambda,I_wave,dI_wave,D_mumu,d_D_mumu);
	 D_mumu = max(0.0,D_mumu+d_D_mumu*(pitch - limitedPitch));
      } else {
	 Real I_wave,dI_wave;
	 fetchIntensities(pos,intensity,I_wave,dI_wave);
	 
	 Real d_lambda;
	 if (simControl.order < 2) {
	    d_lambda = simControl.wavelengthMeshCellSizes[static_cast<int32_t>(pos[3])-1];
	 } else {
	    Real resonantLambda = pitch * maxResonantLambda;
	    d_lambda = fabs(resonantLambda)*simControl.logicalWavelengthCellSize;
	 }
	 dI_wave /= d_lambda;

	 evaluateDiffusionCoefficients(pitch,maxResonantLambda,I_wave,dI_wave,D_mumu,d_D_mumu);
	 D_mumu = max(0.0,D_mumu);
      }
   }

   bool evaluateMasterEquation(const Real* lambdaIntensity,const size_t N_bins,const size_t N_particlesPerBin,
			       const Real& dt_global,const Real& omega,const Real& B_mag,
			       const Real& maxResonantWavelength,const Real& maxSpeed,pargrid::CellID blockLID,
			       corsair::ObjectWrapper& objectWrapper,
			       uint32_t* transitionProbabs,const std::string& name,scatterer scat,
			       int propagationDirection) {
      #if PROFILE_LEVEL > 0
         int profID = -1;
         profile::start("Isotropicity "+name,profID);
      #endif

      // Get array containing min/max valid wavelength bins:
      typedef wmesh::validDatatype validInt;
      validInt* validBins = NULL;
      pargrid::DataID dataID = pargrid::INVALID_DATAID;
      if (propagationDirection < 0) {
	 dataID = simControl.maxAntiparValidWavelengthBinsDataID;
	 validBins = objectWrapper.simClasses.pargrid.getUserDataStatic<validInt>(dataID)
	   + blockLID*block::SIZE*2;
      } else {
	 dataID = simControl.maxParValidWavelengthBinsDataID;
	 validBins = objectWrapper.simClasses.pargrid.getUserDataStatic<validInt>(dataID)
	   + blockLID*block::SIZE*2;
      }

      Real pos[4];
      pos[0] = 0.5;
      pos[1] = 0.5;
      pos[2] = 0.5;
      pos[3] = 0.0;
      
      const int index = block::index((int)pos[0],(int)pos[1],(int)pos[2]);
      const wmesh::validDatatype minValidBin = validBins[index*2+0]+3;
      const wmesh::validDatatype maxValidBin = validBins[index*2+1]  ;
      const Real minWavelength = simControl.wavelengthMeshNodeCoordinates[minValidBin-1];
      const Real maxWavelength = simControl.wavelengthMeshNodeCoordinates[maxValidBin-1];
      const Real minMu = max(-1.0,minWavelength/maxResonantWavelength);
      const Real maxMu = min(+1.0,maxWavelength/maxResonantWavelength);
      
      pos[0] = 1.5;
      pos[1] = 1.5;
      pos[2] = 1.5;
      pos[3] = 0.0;

      uint32_t blockIndices[3];
      pargrid::CellID blockGID = objectWrapper.simClasses.pargrid.getGlobalIDs()[blockLID];
      block::calculateBlockIndices(objectWrapper.sim,blockGID,blockIndices[0],blockIndices[1],blockIndices[2]);
      
      // Evaluate number of scattering substeps for following speeds:
      vector<Real> velocities;
      velocities.push_back(1.4e6); // 10  keV/amu
      velocities.push_back(1.4e7); // 1   MeV/amu
      velocities.push_back(1.4e8); // 100 MeV/amu
      velocities.push_back(1.4e9); // speed of light

      vector<uint32_t> substeps;
      substeps.push_back(evaluateNumberOfSubsteps(objectWrapper.sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[0]));
      substeps.push_back(evaluateNumberOfSubsteps(objectWrapper.sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[1]));
      substeps.push_back(evaluateNumberOfSubsteps(objectWrapper.sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[2]));
      substeps.push_back(evaluateNumberOfSubsteps(objectWrapper.sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[3]));

      // Pick substeps for this particle:
      //scatteringSubsteps = evaluateNumberOfSubsteps(objectWrapper.sim,blockLID,blockIndices,validBins,lambdaIntensity,maxSpeed);
      int bin = max(0,static_cast<int>(1+log10(maxSpeed / 1.4e6)));
      scatteringSubsteps = substeps[bin];
      if (scatteringSubsteps > simControl.currentMaxScatteringSubsteps)
	simControl.currentMaxScatteringSubsteps = scatteringSubsteps;
            
      cout << "EVALUATE MASTER EQUATION" << endl;
      cout << "\t Substeps for pre-evaluated speeds:" << endl;
      for (size_t i=0; i<substeps.size(); ++i) {
	 cout << "\t bin: " << i << "\t speed:\t" << velocities[i] << "\t substeps:\t" << substeps[i] << endl;
      }
      
      cout << "\t N_substeps is " << scatteringSubsteps << endl;
      cout << "\t Scatterer dt is " << dt_global/scatteringSubsteps << endl;
      cout << "\t min,max mu are " << minMu << '\t' << maxMu << endl;
      cout << "\t min,max valid bins are " << minValidBin << '\t' << maxValidBin << endl;
      
      const Real d_mu = 2.0/N_bins;
      for (size_t sourceBin=0; sourceBin<N_bins; ++sourceBin) {
	 for (uint32_t p=0; p<N_particlesPerBin; ++p) {
	    Real mu = -1.0 + (sourceBin + objectWrapper.simClasses.random.uniform()) * d_mu;

	    if (mu >= minMu && mu <= maxMu) {
	       pos[3] = 1.0 + (*simControl.getLogicalWavelength)(maxWavelength*mu);
	       int32_t indices[4];
	       Real shapeFactors[12];
	       
	       switch (simControl.order) {
		case 0:
		  getShapeFactorsNGP_4D(pos,indices,shapeFactors);
		  break;
		case 1:
		  getShapeFactorsCIC_4D(pos,indices,shapeFactors);
		  break;
		case 2:
		  getShapeFactorsTSC_4D(pos,indices,shapeFactors);
		  break;
		default:
		  cerr << "ERROR: Unknown shape factors in in evaluateMasterEquation" << endl;
		  exit(1);
		  break;
	       }
	       
	       mu = (*scat)(objectWrapper.simClasses.random,
			    lambdaIntensity,dt_global,indices,shapeFactors,pos,mu,omega,
			    maxSpeed,B_mag,minMu,maxMu);
	    }

	    const size_t targetBin = (mu+1.0)/d_mu;
	    if (targetBin > N_bins-1) continue;
	    ++transitionProbabs[sourceBin*N_bins+targetBin];
	 }
	 if ((sourceBin+1) % 10 == 0) {
	    cout << "\t" << sourceBin+1 << " / " << N_bins << " bins complete" << endl;
	 }
      }
      
      #if PROFILE_LEVEL > 0
         profile::stop(); //Isotropicity
      #endif
      
      return true;
   }

   /** Evaluate the number of substeps needed when computing pitch angle 
    * scattering for particles in this block.
    * @param sim Struct containing generic simulation parameters.
    * @param blockLID Local ID of block where number of substeps is evaluated.
    * @param blockIndices Block i,j,k indices.
    * @param validBins Array containing min,max wavelength bins that are good for scattering.
    * @param lambdaIntensity Array containing lambda times intensity.
    * @param maxSpeed Particle speed in wave rest frame.*/
   uint32_t evaluateNumberOfSubsteps(Simulation& sim,pargrid::CellID blockLID,const uint32_t* const blockIndices,
				     const wmesh::validDatatype* RESTRICT validBins,
				     const Real* RESTRICT lambdaIntensity,Real maxSpeed) {
      uint32_t substeps = 1;
      
      for (int32_t k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	 // We need correct logical position to get E,B,gradB below:
	 Real pos[4];
	 pos[0] = blockIndices[0]*block::WIDTH_X + 0.5 + i;
	 pos[1] = blockIndices[1]*block::WIDTH_Y + 0.5 + j;
	 pos[2] = blockIndices[2]*block::WIDTH_Z + 0.5 + k;

	 // Get fields at cell center and calculate omega and B_mag for protons:
	 Real E[3];
	 Real B[3];
	 Real gradB[9];
	 (*simControl.fieldsGetFields)(blockLID,sim.t,pos,E,B,gradB);
	 const Real B_mag = vectorMagnitude<3>(B);
	 const Real omega = Q_per_M * B_mag;
	 const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);
	 const Real maxResonantLambda = 2*M_PI*maxSpeed/omega;

	 // In interpolations below we need position as offset to block lower left corner:
	 pos[0] = 1.5 + i;
	 pos[1] = 1.5 + j;
	 pos[2] = 1.5 + k;
	 pos[3] = 1.0;
	 
	 int32_t indices[4];
	 Real shapeFactors[12];
	 switch (simControl.order) {
	  case 0:
	    getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	    break;
	  case 1:
	    getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	    break;
	  case 2:
	    getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	    break;
	 }

	 // Get min,max valid wavelength bins that are good for scattering:
	 const wmesh::validDatatype minValidBin = validBins[block::index(i,j,k)*2+0]+2;
	 const wmesh::validDatatype maxValidBin = validBins[block::index(i,j,k)*2+1]-1;
	 
	 // Calculate min,max logical wavelength range used for scattering:
	 const Real maxWavelength = 2*M_PI*maxSpeed/omega;
	 const Real l_res_min = (*simControl.getLogicalWavelength)(-2*M_PI*maxSpeed/omega);
	 const Real l_res_max = (*simControl.getLogicalWavelength)(+2*M_PI*maxSpeed/omega);
	 int L_index_min = max(static_cast<int>(floor(l_res_min)),(int)minValidBin);
	 int L_index_max = min(static_cast<int>(ceil(l_res_max)),(int)maxValidBin);

	 const Real mu_min = simControl.wavelengthMeshNodeCoordinates[L_index_min] / maxResonantLambda;
	 const Real mu_max = simControl.wavelengthMeshNodeCoordinates[L_index_max] / maxResonantLambda;
	 limitedPitchMin = max(mu_min,-mu_limit);
	 limitedPitchMax = min(mu_max,+mu_limit);

	 // Iterate over L_index_min,L_index_max wavelength bins and 
	 // compute maximum drift and diffusion terms for protons:
	 Real maxDiffusion = -1.0;
	 Real maxDrift     = -1.0;
	 for (int l=L_index_min; l<L_index_max; ++l) {
	    Real diffusionTerm,driftTerm;
	    Real mu = (*simControl.getPhysicalWavelength)(l+0.5)/maxWavelength;
	    if (mu < -1.0) mu = -1.0;
	    if (mu > +1.0) mu = +1.0;

	    evaluateDiffusionCoefficients(mu,maxWavelength,diffConstant,pos,indices,shapeFactors,lambdaIntensity,diffusionTerm,driftTerm);
	    
	    /*
	    Real I_wave,dI_wave;
	    if (mu < 0) {
	       if (mu > limitedPitchMin) {
		  Real d_lambda = -maxWavelength*limitedPitchMin*simControl.logicalWavelengthCellSize;
		  pos[3] = 1 + (*simControl.getLogicalWavelength)(maxWavelength*limitedPitchMin);
		  fetchIntensities(pos,lambdaIntensity,I_wave,dI_wave);
		  I_wave  *= diffConstant;
		  dI_wave *= diffConstant/d_lambda;

		  evaluateDiffusionCoefficients(limitedPitchMin,maxResonantLambda,I_wave,dI_wave,diffusionTerm,driftTerm);
		  diffusionTerm += driftTerm*(mu+limitedPitchMin);
	       } else {
		  Real d_lambda = maxWavelength*fabs(mu)*simControl.logicalWavelengthCellSize;
		  pos[3] = l + 1.5;
		  fetchIntensities(pos,lambdaIntensity,I_wave,dI_wave);
		  I_wave  *= diffConstant;
		  dI_wave *= diffConstant/d_lambda;

		  evaluateDiffusionCoefficients(mu,maxResonantLambda,I_wave,dI_wave,diffusionTerm,driftTerm);
	       }
	    } else {
	       if (mu < limitedPitchMax) {
		  Real d_lambda = maxWavelength*limitedPitchMax*simControl.logicalWavelengthCellSize;
		  pos[3] = 1 + (*simControl.getLogicalWavelength)(maxWavelength*limitedPitchMax);
		  fetchIntensities(pos,lambdaIntensity,I_wave,dI_wave);
		  I_wave  *= diffConstant;
		  dI_wave *= diffConstant/d_lambda;
		  
		  evaluateDiffusionCoefficients(limitedPitchMax,maxResonantLambda,I_wave,dI_wave,diffusionTerm,driftTerm);
		  diffusionTerm += driftTerm*(mu-limitedPitchMax);
	       } else {
		  Real d_lambda = maxWavelength*fabs(mu)*simControl.logicalWavelengthCellSize;
		  pos[3] = l + 1.5;
		  fetchIntensities(pos,lambdaIntensity,I_wave,dI_wave);
		  I_wave  *= diffConstant;
		  dI_wave *= diffConstant/d_lambda;
		  
		  evaluateDiffusionCoefficients(mu,maxResonantLambda,I_wave,dI_wave,diffusionTerm,driftTerm);
	       }
	    }*/

	    if (diffusionTerm > maxDiffusion) {
	       maxDiffusion = diffusionTerm;
	    }
	    if (fabs(driftTerm) > maxDrift) {
	       maxDrift = fabs(driftTerm);
	    }
	 }

	 // Calculate maximum timesteps for drift and diffusion terms:
	 const Real dt_max_diffusion = 0.5*simControl.maxScatteringDeltaMu*simControl.maxScatteringDeltaMu/maxDiffusion;
	 const Real dt_max_drift     = 2.0*simControl.maxScatteringDeltaMu/maxDrift;

	 //cout << "max Dmumu:\t" << maxDiffusion << "(" << dt_max_diffusion << ")\t d_Dmumu:\t" << maxDrift;
	 //cout << "(" << dt_max_drift << ")" << endl;
	 
	 // Calculate required number of scattering substeps:
	 const uint32_t N_substeps = ceil(sim.dt / min(dt_max_diffusion,dt_max_drift));
	 if (N_substeps > substeps) substeps = N_substeps;
      }

      return substeps;
   }
      
   WaveEnergySpectrumBaseClass* getWaveEnergySpectrum(corsair::ObjectWrapper& objectWrapper,int propagationDirection) {
      bool success = true;
      
      // Find particle list that contains Alfven wave packets propagating to given direction,
      // we need to find this to get our hands on the injector below:
      ParticleListBase* particleList = NULL;
      for (size_t i=0; i<objectWrapper.particleLists.size(); ++i) {
	 if (objectWrapper.particleLists[i]->getSpeciesType() != simControl.lagrangianSpeciesTypename) continue;
	 const LagrangianSpecies* species = reinterpret_cast<const LagrangianSpecies*>(objectWrapper.particleLists[i]->getSpecies());
	 if (species->propagationDirection*propagationDirection < 0) continue;
	 particleList = objectWrapper.particleLists[i];
	 break;
      }
      if (particleList == NULL) return NULL;

      // Get wave packet injector name and related config file region:
      string injectorName,injectorRegion;
      particleList->getInjector(injectorName,injectorRegion);

      // Get name and config file region of Alfven wave energy spectrum class:
      string energySpectrumName,energySpectrumParams;
      if (objectWrapper.configReader.get(injectorRegion+".energy_spectrum",energySpectrumName) == false) success = false;
      if (objectWrapper.configReader.get(injectorRegion+".energy_spectrum_parameters",energySpectrumParams) == false) success = false;
      if (success == false) {
	 objectWrapper.simClasses.logger << "(SCATTERER TEST) ERROR: Could not get Alfven wave energy spectrum name" << endl << write;
	 return NULL;
      }

      // Create Alfven wave energy spectrum function for tests:
      WaveEnergySpectrumBaseClass* energySpectrum = sep::getObjectWrapper().waveEnergySpectrumFactory.create(energySpectrumName);
      if (energySpectrum == NULL) {
	 objectWrapper.simClasses.logger << "(SCATTERER TEST) ERROR: Could not create Alfven wave energy spectrum class" << endl << write;
	 return NULL;
      }

      // Init energy spectrum class:
      if (energySpectrum->initialize(objectWrapper.sim,objectWrapper.simClasses,objectWrapper.configReader,energySpectrumParams) == false) {
	 objectWrapper.simClasses.logger << "(SCATTERER TEST) ERROR: Failed to init Alfven wave energy spectrum class" << endl << write;
	 energySpectrum->finalize();
	 delete energySpectrum; energySpectrum = NULL;
      }
      return energySpectrum;
   }

   Real interpolateWaveEnergy(const Real* RESTRICT pos,const Real* RESTRICT waveEnergy) {
      const uint32_t NWL = simControl.N_wavelengthMeshCells+2;
      
      // Calculate shape factors:
      int32_t indices[4];
      Real shapeFactors[12];
      for (int i=0; i<4; ++i) indices[i] = 1;
      for (int i=0; i<12; ++i) shapeFactors[i] = 0.0;
      switch (simControl.order) {
       case 0:
	 getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	 break;
       case 1:
	 getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	 break;
       case 2:
	 getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	 break;
       default:
	 break;
      }

      // Interpolate wave energy to phase-space position pos:
      Real U_wave = 0.0;
      switch (simControl.order) {
       case 0:
	 U_wave = interpolateScalarNGP_4D(waveEnergy,indices,shapeFactors,NWL);
	 break;
       case 1:
	 /*
	 // Clamp shape factor so that R-helicity energy does not flow
	 // into L-helicity side of lambda and vice versa:
	 if (pos[3]-1.5 >= simControl.N_wavelengthMeshCells/2-1 && pos[3]-1 < simControl.N_wavelengthMeshCells/2) {
	    shapeFactors[6] = 1.0;
	    shapeFactors[7] = 0.0;
	 } else if (pos[3]-1 >= simControl.N_wavelengthMeshCells/2 && pos[3]-1 < simControl.N_wavelengthMeshCells/2+0.5) {
	    shapeFactors[6] = 0.0;
	    shapeFactors[7] = 1.0;
	 }*/
	 
	 U_wave = interpolateScalarCIC_4D(waveEnergy,indices,shapeFactors,NWL);
	 break;
       case 2:
	 /*
	 // Clamp shape factor so that R-helicity energy does not flow
	 // into L-helicity side of lambda and vice versa:
	 if (l-1 == simControl.N_wavelengthMeshCells/2-1) {
	    shapeFactors[11] = 0.0;
	 } else if (l-1 == simControl.N_wavelengthMeshCells/2) {
	    shapeFactors[9] = 0.0;
	 }
	 */
	 
	 /*
	 if (indices[3] <= l_min+1 || indices[3] >= l_max-1) return 0.0;	 
	 if (indices[3] == l_min+2) {
	    shapeFactors[10] += shapeFactors[9];
	    shapeFactors[9] = 0.0;
	    //shapeFactors[11] = 0.0;
	 } else if (indices[3] == l_max-2) {
	    shapeFactors[10] += shapeFactors[11];
	    //shapeFactors[9] = 0.0;
	    shapeFactors[11] = 0.0;
	 }*/
	 
	 U_wave = interpolateScalarTSC_4D(waveEnergy,indices,shapeFactors,NWL);
	 break;
       default:
	 break;
      }
      return U_wave;
   }
   
   Real interpolateWaveEnergyDerivated(const Real* RESTRICT pos,const Real* RESTRICT waveEnergy) {
      const int32_t NWL = simControl.N_wavelengthMeshCells+2;
     
      // Calculate shape factors:
      int32_t indices[4];
      Real shapeFactors[12];
      Real U_wave = 0.0;
      for (int i=0; i<4; ++i) indices[i] = 1;
      for (int i=0; i<12; ++i) shapeFactors[i] = 0.0;
      switch (simControl.order) {
       case 0:
	 //getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	 break;
       case 1:
	 getDerivativeShapeFactorsCIC_4D(pos,indices,shapeFactors);
	 if (indices[3] < 1 || indices[3] >= NWL-1) return 0.0;
	 U_wave = interpolateScalarCIC_4D(waveEnergy,indices,shapeFactors,NWL);
	 break;
       case 2:
	 getDerivativeShapeFactorsTSC_4D(pos,indices,shapeFactors);
	 /*
	 // Clamp shape factor so that R-helicity energy does not flow
	 // into L-helicity side of lambda and vice versa:
	 if (l-1 == simControl.N_wavelengthMeshCells/2-1) {
	    shapeFactors[11] = 0.0;
	 } else if (l-1 == simControl.N_wavelengthMeshCells/2) {
	    shapeFactors[9] = 0.0;
	 }
	 */
	 
	 /*
	 if (indices[3] <= l_min+1 || indices[3] >= l_max-1) return 0.0;
	 if (indices[3] == l_min+2) {
	    shapeFactors[10] = -shapeFactors[11];
	    shapeFactors[9] = 0.0;
	 } else if (indices[3] == l_max-2) {
	    shapeFactors[10] = -shapeFactors[9];
	    shapeFactors[11] = 0.0;
	 }*/
	 
	 if (indices[3] < 1 || indices[3] >= NWL-1) return 0.0;
	 U_wave = interpolateScalarTSC_4D(waveEnergy,indices,shapeFactors,NWL);
	 break;
       default:
	 break;
      }
      return U_wave;
   }
   
   Real interpolateWaveEnergyLimited(const PARTICLE& particle,const Real* RESTRICT pos,Real* RESTRICT shapeFactors,
				     Real* RESTRICT shapeFactors2,int32_t* indices,uint32_t NWL,const Real* waveEnergyArray) {
      Real waveEnergy = 0.0;
      switch (simControl.order) {
       case 0:
	 getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	 for (int32_t i=0; i<12; ++i) shapeFactors2[i] = shapeFactors[i]*shapeFactors[i];
	 waveEnergy = interpolateScalarNGP_4D(waveEnergyArray,indices,shapeFactors2,NWL);
	 break;
       case 1:
	 getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	 for (int32_t i=0; i<12; ++i) shapeFactors2[i] = shapeFactors[i]*shapeFactors[i];
	 waveEnergy = interpolateScalarCIC_4D(waveEnergyArray,indices,shapeFactors2,NWL);
	 break;
       case 2:
	 getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	 for (int32_t i=0; i<12; ++i) shapeFactors2[i] = shapeFactors[i]*shapeFactors[i];
	 waveEnergy = interpolateScalarTSC_4D(waveEnergyArray,indices,shapeFactors2,NWL);
	 break;
       default:
	 cerr << "(SEP PARTICLE SCATTERER) ERROR: Unsupported order of accuracy" << endl;
	 exit(1);
	 break;
      }
      return waveEnergy * particle.state[sep::particle::WEIGHT];
   }
   
   bool loadDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,
				 Real cellVolume,int propagationDirection,Real* RESTRICT waveIntensity) {
      // Get global array containing wave energy:
      const Real* RESTRICT waveEnergyGlobal = NULL;
      if (propagationDirection > 0) {
	 waveEnergyGlobal = objectWrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      } else {
	 waveEnergyGlobal = objectWrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      }
      
      if (waveEnergyGlobal == NULL) {
	 objectWrapper.simClasses.logger << "(SEP PARTICLE SCATTERER) ERROR: Failed to get Alfven wave energy array" << endl << write;
	 return false;
      }

      // Load wave energy to temporary array:
      loadScalar4D(&objectWrapper.simClasses,blockLID,waveEnergyGlobal,waveIntensity,
		   simControl.N_wavelengthMeshCells);

      const size_t NWL = simControl.N_wavelengthMeshCells+2;
      
      // Convert wave energy to intensity times fabs(wavelength):
      for (int cell=0; cell<block::SIZE_PLUS_ONE_LAYER; ++cell) {
	 waveIntensity[cell*NWL+0] = 0.0;
	 
	 for (uint32_t l=1; l<NWL-1; ++l) {
	    Real value = waveIntensity[cell*NWL+l];
	    
	    waveIntensity[cell*NWL+l] = value * constants::PERMEABILITY
	      / (cellVolume*simControl.wavelengthMeshCellSizes[l-1]);
	 }

	 waveIntensity[cell*NWL+NWL-1] = 0.0;
      }

      return true;
   }

   bool loadIntensity(pargrid::CellID blockLID,Real* array,int alfvenSign) {
      const int32_t NWL = simControl.N_wavelengthMeshCells+2;
      const int32_t WAVE_SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;

      // Calculate magnetic field magnitude squared, 
      // needed for Bohm limiting:
      const pargrid::CellID blockGID = corsair::getObjectWrapper().simClasses.pargrid.getGlobalIDs()[blockLID];
      int32_t blockIndices[3];
      block::calculateBlockIndices(corsair::getObjectWrapper().sim,blockGID,blockIndices[0],blockIndices[1],blockIndices[2]);

      Real pos[3];
      pos[0] = (blockIndices[0]+0.5)*block::WIDTH_X;
      pos[1] = (blockIndices[1]+0.5)*block::WIDTH_Y;
      pos[2] = (blockIndices[2]+0.5)*block::WIDTH_Z;

      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(corsair::getObjectWrapper().sim.t,blockLID,pos,E,B,gradB);
      Real B_mag2 = vectorMagnitude2<3>(B);
      //cerr << pos[0] << '\t' << blockIndices[0] << '\t' << sqrt(B_mag2) << endl;

      // Load spectral wave energy:
      for (int32_t i=0; i<WAVE_SIZE; ++i) array[i] = 0.0;
      Real* waveEnergy = NULL;
      if (alfvenSign < 0) {
	 waveEnergy = corsair::getObjectWrapper().simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      } else {
	 waveEnergy = corsair::getObjectWrapper().simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      }
      if (waveEnergy == NULL) return false;
      loadScalar4D(&corsair::getObjectWrapper().simClasses,blockLID,waveEnergy,array,simControl.N_wavelengthMeshCells);

      // Load spatial cell volumes:
      Real* tmpCellVolumes = new Real[block::SIZE_PLUS_ONE_LAYER];
      for (int32_t i=0; i<block::SIZE_PLUS_ONE_LAYER; ++i) tmpCellVolumes[i] = numeric_limits<Real>::max();
      block::fetchValues3D(corsair::getObjectWrapper().simClasses,blockLID,tmpCellVolumes,simControl.cellVolumes);

      for (int k=0; k<block::WIDTH_Z+2; ++k) for (int j=0; j<block::WIDTH_Y+2; ++j) for (int i=0; i<block::WIDTH_X+2; ++i) {
	 for (int l=0; l<simControl.N_wavelengthMeshCells/2; ++l) {
	    const Real ONE = 1.0;
	    Real lambdaMin = max(ONE,-simControl.wavelengthMeshNodeCoordinates[l+1]);
	    Real lambdaMax = -simControl.wavelengthMeshNodeCoordinates[l  ];
	    Real maxEnergyDens = simControl.bohmLimitCoefficient*B_mag2/M_PI*log10(lambdaMax/lambdaMin)/constants::PERMEABILITY*2.3;

	    if (array[block::arrayIndex(i,j,k)*NWL+l+1]/tmpCellVolumes[block::arrayIndex(i,j,k)] > maxEnergyDens) {
	       array[block::arrayIndex(i,j,k)*NWL+l+1] = maxEnergyDens*tmpCellVolumes[block::arrayIndex(i,j,k)];
	    }
	 }
	 for (int l=simControl.N_wavelengthMeshCells/2; l<simControl.N_wavelengthMeshCells; ++l) {
	    const Real ONE = 1.0;
	    Real lambdaMin = max(ONE,simControl.wavelengthMeshNodeCoordinates[l]);
	    Real lambdaMax = simControl.wavelengthMeshNodeCoordinates[l+1];
	    Real maxEnergyDens = simControl.bohmLimitCoefficient*B_mag2/M_PI*log10(lambdaMax/lambdaMin)/constants::PERMEABILITY*2.3;

	    if (array[block::arrayIndex(i,j,k)*NWL+l+1]/tmpCellVolumes[block::arrayIndex(i,j,k)] > maxEnergyDens) {
	       array[block::arrayIndex(i,j,k)*NWL+l+1] = maxEnergyDens*tmpCellVolumes[block::arrayIndex(i,j,k)];
	    }
	 }
      }

      // Convert spectral wave energy to intensity:
      for (int k=0; k<block::WIDTH_Z+2; ++k) for (int j=0; j<block::WIDTH_Y+2; ++j) for (int i=0; i<block::WIDTH_X+2; ++i) {
	 const Real spatCellVolume = tmpCellVolumes[block::arrayIndex(i,j,k)];
	 Real d_lambda = simControl.wavelengthMeshCellSizes[0];
	 array[block::arrayIndex(i,j,k)*NWL+0] *= constants::PERMEABILITY / (d_lambda*spatCellVolume);
	 
	 for (int l=1; l<NWL-1; ++l) {
	    d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	    array[block::arrayIndex(i,j,k)*NWL+l] *= constants::PERMEABILITY / (d_lambda*spatCellVolume);
	 }

	 d_lambda = simControl.wavelengthMeshCellSizes[NWL-2];
	 array[block::arrayIndex(i,j,k)*NWL+NWL-1] *= constants::PERMEABILITY / (d_lambda*spatCellVolume);
      }
      delete [] tmpCellVolumes; tmpCellVolumes = NULL;

      return true;
   }

   bool loadWaveIntensity(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real cellVolume,int propagationDirection,
			  Real* RESTRICT waveIntensity) {
      // Get global array containing wave energy:
      const Real* RESTRICT waveEnergyGlobal = NULL;
      if (propagationDirection > 0) {
	 waveEnergyGlobal = objectWrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.parAlfvenWaveEnergyDataID);
      } else {
	 waveEnergyGlobal = objectWrapper.simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparAlfvenWaveEnergyDataID);
      }
      
      // Load wave energy to temporary array:
      loadScalar4D(&objectWrapper.simClasses,blockLID,waveEnergyGlobal,waveIntensity,simControl.N_wavelengthMeshCells);

      // Convert wave energy to intensity:
      const size_t NWL = simControl.N_wavelengthMeshCells+2;
      for (int32_t cell=0; cell<block::SIZE_PLUS_ONE_LAYER; ++cell) {
	 waveIntensity[cell*NWL+0] = 0.0;
	 for (uint32_t l=1; l<NWL-1; ++l) {
	    Real value = waveIntensity[cell*NWL+l];
	    waveIntensity[cell*NWL+l] = value
	      * constants::PERMEABILITY/(cellVolume*simControl.wavelengthMeshCellSizes[l-1]);
	 }
	 waveIntensity[cell*NWL+NWL-1] = 0.0;
      }
      return true;
   }

   bool runScattererTests(Simulation& sim,SimulationClasses& simClasses) {
      bool success = true;

      const pargrid::CellID blockGID = 0;
      const pargrid::CellID blockLID = simClasses.pargrid.getLocalID(blockGID);
      if (blockLID == pargrid::INVALID_CELLID) {
	 cerr << "P#" << sim.mpiRank << " GID#" << blockGID << " not found on this process" << endl;
	 return true;
      }

      corsair::ObjectWrapper& objectWrapper = corsair::getObjectWrapper();
      
      // Calculate block indices:
      uint32_t i_block,j_block,k_block;
      block::calculateBlockIndices(objectWrapper.sim,blockGID,i_block,j_block,k_block);

      // Calculate phase-space cell volume:
      const Real cellVolume = simControl.cellVolumes[blockLID*block::SIZE];
      
      // Get electromagnetic fields at cell centroid:
      Real pos[4];
      Real E[3];
      Real B[3];
      Real gradB[9];
      pos[0] = i_block + 0.5;
      pos[1] = j_block + 0.5;
      pos[2] = k_block + 0.5;
      (*simControl.fieldsGetFields)(blockLID,objectWrapper.sim.t,pos,E,B,gradB);
      const Real B_mag = vectorMagnitude<3>(B);

      // Write parallel and antiparallel pitch angle diffusion coefficients:
      if (runTestDiffusionCoefficient(blockLID,objectWrapper,B_mag,cellVolume,-1) == false) success = false;
      if (runTestDiffusionCoefficient(blockLID,objectWrapper,B_mag,cellVolume,+1) == false) success = false;
      
      // Write parallel and antiparallel wave intensities to output files:
      if (runTestWaveIntensity(blockLID,objectWrapper,B_mag,cellVolume,-1) == false) success = false;
      if (runTestWaveIntensity(blockLID,objectWrapper,B_mag,cellVolume,+1) == false) success = false;
      
      // Run scatterer isotropicity tests:
      if (runTestScatteringIsotropicity(blockLID,objectWrapper,B_mag,cellVolume,-1) == false) success = false;
      if (runTestScatteringIsotropicity(blockLID,objectWrapper,B_mag,cellVolume,+1) == false) success = false;
      
      return success;
   }
   
   bool runTestDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,int propagationDirection) {
      const int N_points = 1000;
      const int N_arrays = 4;
      DiffCoeff diffCoeff(N_points,N_arrays);
      
      Real wavelengthConstants[N_arrays];
      wavelengthConstants[0] = 1.05e-8;
      wavelengthConstants[1] = 1.05e-7;
      wavelengthConstants[2] = 1.05e-6;
      wavelengthConstants[3] = 1.05e-5;

      Real wavelengths[N_arrays];
      
      bool success = true;
      for (int i=0; i<N_arrays; ++i) {
	 if (runTestDiffusionCoefficient(blockLID,objectWrapper,B_mag,cellVolume,wavelengthConstants[i],
					 diffCoeff,i,propagationDirection,wavelengths[i]) == false) success = false;
      }
      if (success == false) return success;

      // Write out results:
      string fname;
      if (propagationDirection > 0) fname = "test_pitch_angle_diff_coeff_par.txt";
      else fname = "test_pitch_angle_diff_coeff_antipar.txt";
      fstream out(fname.c_str(), fstream::out);
      if (out.good() == false) return false;

      // Write file header:
      out << "#Proton energies are:" << endl;
      for (int n=0; n<N_arrays; ++n) {
	 Real omega = constants::CHARGE_ELEMENTARY/constants::MASS_PROTON*B_mag;
	 Real maxLambda = wavelengths[n];
	 Real speed = maxLambda*omega/2/M_PI;
	 Real U = 0.5*constants::MASS_PROTON*speed*speed/constants::CHARGE_ELEMENTARY/1e3;
	 out << "# " << U << " keV" << endl;
      }
      
      out << "#column values are:" << endl;
      out << "#1 : pitch" << endl;
      for (int n=0; n<N_arrays; ++n) {
	 out << "#" << 2+n*3+0 << " : wavelength bin" << endl;
	 out << "#" << 2+n*3+1 << " : D_mumu interpolated for lambda_res/max_lambda " << wavelengthConstants[n] << endl;
	 out << "#" << 2+n*3+2 << " : D_mumu analytic                               " << wavelengthConstants[n] << endl;
      }
      for (int n=0; n<N_arrays; ++n) {
	 out << "#" << 2+N_arrays*3+3*n+0 << " : D_mumu_dmu interpolated for lambda_res/max_lambda " << wavelengthConstants[n] << endl;
	 out << "#" << 2+N_arrays*3+3*n+1 << " : D_mumu_dmu analytic                               " << wavelengthConstants[n] << endl;
	 out << "#" << 2+N_arrays*3+3*n+2 << " : D_mumu_dmu finite difference                      " << wavelengthConstants[n] << endl;
      }
      
      // Write D_mumu and its mu-derivative values:
      for (size_t i=0; i<diffCoeff.mu.size(); ++i) {
	 out << diffCoeff.mu[i] << '\t';
	 for (int j=0; j<N_arrays; ++j) {
	    out << diffCoeff.lambda[j][i] << '\t';
	    out << diffCoeff.D_mumu_interpolated[j][i] << '\t';
	    out << diffCoeff.D_mumu_analytic[j][i] << '\t';
	 }
	 for (int j=0; j<N_arrays; ++j) {
	    out << diffCoeff.d_D_mumu_dmu_interpolated[j][i] << '\t';
	    out << diffCoeff.d_D_mumu_dmu_analytic[j][i] << '\t';
	    out << diffCoeff.d_D_mumu_dmu_diff[j][i] << '\t';
	 }
	 out << endl;
      }
      out.close();

      return success;
   }
   
   bool runTestDiffusionCoefficient(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,
				    Real wavelengthConstant,DiffCoeff& diffCoeff,int index,int propagationDirection,Real& wavelength) {
      // Load wave intensity to temporary array:
      const size_t NWL  = simControl.N_wavelengthMeshCells+2;
      const size_t SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      Real* waveEnergy = new Real[SIZE];
      for (size_t i=0; i<SIZE; ++i) waveEnergy[i] = 0.0;
      if (loadDiffusionCoefficient(blockLID,objectWrapper,cellVolume,propagationDirection,waveEnergy) == false) return false;

      // Get Alfven wave energy spectrum class:
      WaveEnergySpectrumBaseClass* energySpectrum = getWaveEnergySpectrum(objectWrapper,propagationDirection);
      if (energySpectrum == NULL) return false;
      
      // Get maximum wavelengths:
      const Real maxWavelength_L = wavelengthConstant*energySpectrum->getMaximumWavelength(B_mag,+1);
      const Real maxWavelength_R = wavelengthConstant*fabs(energySpectrum->getMaximumWavelength(B_mag,-1));
      wavelength = maxWavelength_L;
      const Real omega = constants::CHARGE_ELEMENTARY*B_mag/constants::MASS_PROTON;
      const Real diffConstant = 0.5*M_PI*omega/B_mag/B_mag;
      const Real maxSpeed_L = 0.5*maxWavelength_L*omega/M_PI;
      const Real maxSpeed_R = 0.5*maxWavelength_R*omega/M_PI;
      const Real maxResonantLambda_L = 2*M_PI*maxSpeed_L/omega;
      const Real maxResonantLambda_R = 2*M_PI*maxSpeed_R/omega;
      Real minMu = max(-1.0,simControl.wavelengthMeshNodeCoordinates[0]     / maxResonantLambda_R);
      Real maxMu = min(+1.0,simControl.wavelengthMeshNodeCoordinates[NWL-2] / maxResonantLambda_L);

      Real speed = maxSpeed_L;
      Real energy = 0.5*constants::MASS_PROTON*speed*speed / constants::CHARGE_ELEMENTARY / 1e3;
      cout << "Evaluating diffusion coefficient for speed:\t" << speed << "\t energy:\t" << energy << " keV \t lambda:\t" << maxResonantLambda_L << endl;
      
      limitedPitchMin = max(-mu_limit,minMu);
      limitedPitchMax = min(+mu_limit,maxMu);

      // Calculate D_mumu due to R helicity waves:
      for (uint32_t i=0; i<diffCoeff.mu.size(); ++i) {
	 const Real mu = diffCoeff.mu[i];
	 Real wavelength = fabs(maxWavelength_R) * mu;
	 if (mu > 0.0) wavelength = maxWavelength_L * mu;
	 
	 Real pos[4];
	 pos[0] = 1.5;
	 pos[1] = 1.5;
	 pos[2] = 1.5;
	 pos[3] = 0.0;
	 Real l = (*simControl.getLogicalWavelength)(wavelength);

	 int32_t indices[4];
	 Real shapeFactors[12];
	 switch (simControl.order) {
	  case 0:
	    getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	    break;
	  case 1:
	    getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	    break;
	  case 2:
	    getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	    break;
	  default:
	    cerr << "Unknown order of accuracy in runTestDiffusionCoefficient" << endl;
	    exit(1);
	    break;
	 }

	 Real D_mumu_interpolated,D_mumu_dmu_interp;
	 evaluateDiffusionCoefficients(mu,maxWavelength_R,diffConstant,pos,indices,shapeFactors,
				       waveEnergy,D_mumu_interpolated,D_mumu_dmu_interp);

	 // Calculate analytic diffusion coefficient:
	 const Real dI_analytic = energySpectrum->getIntensityDerivated(B_mag,wavelength);
	 const Real I_analytic = energySpectrum->getIntensity(B_mag,wavelength);

	 //const Real D_mumu_analytic = I_analytic;
	 const Real D_mumu_analytic     = diffConstant*(1.0-mu*mu)*fabs(wavelength)*I_analytic;
	 Real D_mumu_dmu_analytic = 0.0;
	 if (mu < 0.0) {
	    D_mumu_dmu_analytic =
	      - diffConstant*(1-mu*mu)*2*M_PI*maxSpeed_R/omega*I_analytic
	      - diffConstant*2*mu*fabs(wavelength)*I_analytic
	      + diffConstant*(1-mu*mu)*fabs(wavelength)*2*M_PI*maxSpeed_R/omega*dI_analytic;
	 } else {
	    D_mumu_dmu_analytic =
	      + diffConstant*(1-mu*mu)*2*M_PI*maxSpeed_L/omega*I_analytic
	      - diffConstant*2*mu*fabs(wavelength)*I_analytic
	      + diffConstant*(1-mu*mu)*fabs(wavelength)*2*M_PI*maxSpeed_L/omega*dI_analytic;
	 }

	 diffCoeff.lambda[index][i] = l;
	 diffCoeff.D_mumu_analytic[index][i] = D_mumu_analytic;
	 diffCoeff.D_mumu_interpolated[index][i] = D_mumu_interpolated;
	 diffCoeff.d_D_mumu_dmu_interpolated[index][i] = D_mumu_dmu_interp;
	 diffCoeff.d_D_mumu_dmu_analytic[index][i] = D_mumu_dmu_analytic;
      }
      
      // Calculate derivative of D_mumu using finite differences:
      diffCoeff.d_D_mumu_dmu_diff[index][0] = 
	  diffCoeff.D_mumu_interpolated[index][1]
	- diffCoeff.D_mumu_interpolated[index][0];
      diffCoeff.d_D_mumu_dmu_diff[index][0] /= (diffCoeff.mu[1]-diffCoeff.mu[0]);
      
      for (uint32_t i=1; i<diffCoeff.mu.size()-1; ++i) {
	 diffCoeff.d_D_mumu_dmu_diff[index][i] = 
	     diffCoeff.D_mumu_interpolated[index][i+1] 
	   - diffCoeff.D_mumu_interpolated[index][i-1];
	 diffCoeff.d_D_mumu_dmu_diff[index][i] /= (diffCoeff.mu[i+1]-diffCoeff.mu[i-1]);
      }
      
      diffCoeff.d_D_mumu_dmu_diff[index][diffCoeff.mu.size()-1] = 
	  diffCoeff.D_mumu_interpolated[index][diffCoeff.mu.size()-1]
	- diffCoeff.D_mumu_interpolated[index][diffCoeff.mu.size()-2];
      diffCoeff.d_D_mumu_dmu_diff[index][diffCoeff.mu.size()-1] /=
	(diffCoeff.mu[diffCoeff.mu.size()-1]-diffCoeff.mu[diffCoeff.mu.size()-2]);
      
      // Deallocate memory and exit:
      delete [] waveEnergy; waveEnergy = NULL;
      energySpectrum->finalize();
      delete energySpectrum; energySpectrum = NULL;
      return true;
   }
   
   bool runTestScatteringIsotropicity(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,
				      Real B_mag,Real cellVolume,int propagationDirection) {
      // Load wave intensity into temporary wave array:
      const size_t NWL  = simControl.N_wavelengthMeshCells+2;
      const size_t SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      Real* waveEnergy = new Real[SIZE];
      for (size_t i=0; i<SIZE; ++i) waveEnergy[i] = 0.0;
      
      loadIntensity(blockLID,waveEnergy,propagationDirection);

      /*
      if (loadDiffusionCoefficient(blockLID,objectWrapper,cellVolume,
				   propagationDirection,waveEnergy) == false) {
	 objectWrapper.simClasses.logger << "(SEP PARTICLE SCATTERER) ERROR: Could not run isotropicity test, ";
	 objectWrapper.simClasses.logger << "loadDiffusionCoefficient failed" << endl << write;
	 return false;
      }*/

      // Get Alfven wave energy spectrum class:
      WaveEnergySpectrumBaseClass* energySpectrum 
	= getWaveEnergySpectrum(objectWrapper,propagationDirection);
      if (energySpectrum == NULL) return false;

      // Get maximum wavelengths:
      Real maxWavelength;
      const Real maxWavelength_L = 1.0*energySpectrum->getMaximumWavelength(B_mag,+1);
      const Real maxWavelength_R = 1.0*fabs(energySpectrum->getMaximumWavelength(B_mag,-1));

      const Real omega = constants::CHARGE_ELEMENTARY*B_mag/constants::MASS_PROTON;
      Q_per_M = constants::CHARGE_ELEMENTARY/constants::MASS_PROTON;
      Real maxSpeed;
      Real maxSpeedFactor = 0.0000105;
      const Real maxSpeed_L = maxSpeedFactor*0.5*maxWavelength_L*omega/M_PI;
      const Real maxSpeed_R = maxSpeedFactor*0.5*maxWavelength_R*omega/M_PI;
      if (propagationDirection < 0) {
	 maxSpeed      = maxSpeed_R;
	 maxWavelength = maxSpeedFactor*maxWavelength_R;
      } else {
	 maxSpeed      = maxSpeed_L;
	 maxWavelength = maxSpeedFactor*maxWavelength_L;
      }
      
      cout << "Scattering isotropicity evaluated for speed " << maxSpeed << " energy ";
      cout << 0.5*constants::MASS_PROTON*maxSpeed*maxSpeed/constants::CHARGE_ELEMENTARY/1000;
      cout << " keV" << " max resonant wavelength " << maxWavelength << endl;

      const int32_t N_pitchBins = 160;
      const Real d_mu = 2.0/N_pitchBins;

      const Real dt_global = objectWrapper.sim.dt;
      const uint32_t N_particlesPerBin = 3000;

      vector<size_t> sourceBins;
      sourceBins.push_back(6);
      sourceBins.push_back(N_pitchBins/5-2);
      sourceBins.push_back(N_pitchBins/5-1);
      sourceBins.push_back(N_pitchBins/2-1);
      sourceBins.push_back(5*N_pitchBins/8-9);
      sourceBins.push_back(N_pitchBins-2);
      
      // Write out pitch distribution bin mu values:
      for (size_t i=0; i<sourceBins.size(); ++i) {
	 cout << i << ": bin=" << sourceBins[i] << " mu : ";
	 cout << -1.0 + 0.5*(sourceBins[i]+sourceBins[i]+1)*d_mu << endl;
      }

      vector<scatterer> scatterers;
      vector<string> scattererNames;
      
      if (simControl.wavelengthMeshType == mesh::LINEAR) {
	 scatterers.push_back(scatterPredCorr1);
      } else {
	 //scatterers.push_back(scatterIsotropic);
	 //scatterers.push_back(scatterIsotropicAnalytic);
	 scatterers.push_back(scatterPredCorr1);
      }
      //scattererNames.push_back("Isotropic");
      //scattererNames.push_back("IsotropicAnalytic");
      scattererNames.push_back("PredCorr1");
      //scatterers.push_back(scatterPredCorr2);
      //scattererNames.push_back("PredCorr2");
      //scatterers.push_back(scatterMilstein);
      //scattererNames.push_back("Milstein");
      
      uint32_t* transitionProbabilities = new uint32_t[N_pitchBins*N_pitchBins];
      Real** probabilityDistributions = new Real*[sourceBins.size()];
      for (size_t i=0; i<sourceBins.size(); ++i) probabilityDistributions[i] = new Real[N_pitchBins];
      
      bool success = true;
      for (size_t e=0; e<scatterers.size(); ++e) {
	 for (size_t i=0; i<N_pitchBins*N_pitchBins; ++i) transitionProbabilities[i] = 0;

	 if (evaluateMasterEquation(waveEnergy,N_pitchBins,N_particlesPerBin,dt_global,
				    omega,B_mag,maxWavelength,maxSpeed,blockLID,objectWrapper,
				    transitionProbabilities,scattererNames[e],scatterers[e],
				    propagationDirection) == false) success = false;
	 
	 if (writeIsotropicity(N_particlesPerBin,N_pitchBins,
			       transitionProbabilities,scattererNames[e],
			       propagationDirection) == false)
	   success = false;
	 if (writeTransitionProbabilities(N_particlesPerBin,N_pitchBins,sourceBins,
					  transitionProbabilities,scattererNames[e],
					  propagationDirection) == false)
	   success = false;
      }

      // Deallocate memory:
      delete [] transitionProbabilities;
      transitionProbabilities = NULL;
      for (size_t i=0; i<sourceBins.size(); ++i) {
	 delete [] probabilityDistributions[i];
	 probabilityDistributions[i] = NULL;
      }
      delete [] probabilityDistributions;
      
      // Deallocate analytic wave intensity spectrum class:
      energySpectrum->finalize();
      delete energySpectrum; energySpectrum = NULL;
      
      return success;
   }
   
   bool runTestWaveIntensity(pargrid::CellID blockLID,corsair::ObjectWrapper& objectWrapper,Real B_mag,Real cellVolume,int propagationDirection) {
      // Create temporary wave energy array:
      const size_t NWL  = simControl.N_wavelengthMeshCells+2;
      const size_t SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      Real* waveEnergy = new Real[SIZE];
      for (size_t i=0; i<SIZE; ++i) waveEnergy[i]      = 0.0;
      if (loadWaveIntensity(blockLID,objectWrapper,cellVolume,propagationDirection,waveEnergy) == false) return false;

      // Get Alfven wave energy spectrum class:
      WaveEnergySpectrumBaseClass* energySpectrum = getWaveEnergySpectrum(objectWrapper,propagationDirection);
      if (energySpectrum == NULL) return false;
      
      // Open output file and write header:
      string fname;
      if (propagationDirection > 0) fname = "test_spectral_wave_intensity_par.txt";
      else fname = "test_spectral_wave_intensity_antipar.txt";
      fstream out(fname.c_str(), fstream::out);
      if (out.good() == false) {
	 energySpectrum->finalize();
	 delete energySpectrum; energySpectrum = NULL;
	 return false;
      }
      out << "#lambda(logical) lambda(phys) I_lambda(interp) I_lambda(analytic) dI(interp) dI(analytic)" << endl;

      // Interpolate wave intensity and energy to discrete points, and write interpolated 
      // + analytic values to output file:
      const int resolution = 4;
      Real energy = 0.0;
      for (uint32_t l=1; l<simControl.N_wavelengthMeshCells+1; ++l) {
	 energy += waveEnergy[block::arrayIndex(1,1,1)*NWL+l]*simControl.wavelengthMeshCellSizes[l-1];
	 for (int i=0; i<resolution; ++i) {
	    // Logical position where to interpolate wave energy:
	    Real pos[4];
	    pos[0] = 1.5;
	    pos[1] = 1.5;
	    pos[2] = 1.5;
	    pos[3] = l + 0.5/resolution + (1.0*i)/resolution;

	    const Real wavelength = (*simControl.getPhysicalWavelength)(pos[3]-1);
	    
	    // Calculate spectral intensity:
	    const Real I_lambda = interpolateWaveEnergy(pos,waveEnergy);
	    const Real dI_lambda = interpolateWaveEnergyDerivated(pos,waveEnergy);

	    Real d_lambda = NAN;
	    switch (simControl.order) {
	     case 0:
	       d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	       break;
	     case 1:
	       d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	       break;
	     case 2:
	       d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	       break;
	    }

	    const Real I_analytic = energySpectrum->getIntensity(B_mag,wavelength);
	    const Real dI_analytic = energySpectrum->getIntensityDerivated(B_mag,wavelength);

	    // Calculate abs(wavelength) times intensity:
	    const Real lambda_I = fabs(wavelength)*I_lambda;
	    Real d_lambda_I = I_lambda;
	    if (wavelength < 0.0) d_lambda_I *= -1;
	    d_lambda_I += fabs(wavelength)*dI_lambda/fabs(wavelength);

	    Real lambda_I_analytic = fabs(wavelength)*I_analytic;
	    Real d_lambda_I_analytic = I_analytic;
	    if (wavelength < 0.0) d_lambda_I_analytic *= -1;
	    d_lambda_I_analytic += fabs(wavelength)*dI_analytic;

	    // Write values to output file:
	    out << pos[3]-1 << '\t';
	    out << wavelength << '\t';
	    
	    if (wavelength < 0.0) out << min(-1.0,-log10(-wavelength)) << '\t';
	    else out << max(1.0,log10(wavelength)) << '\t';
	    
	    out << I_lambda << '\t';
	    out << I_analytic << '\t';
	    out << dI_lambda/d_lambda << '\t';
	    out << dI_analytic << '\t';

	    out << lambda_I << '\t';
	    out << d_lambda_I << '\t';
	    out << lambda_I_analytic << '\t';
	    out << d_lambda_I_analytic << '\t';

	    out << endl;
	 }
      }
      
      energy /= constants::PERMEABILITY;
      const Real magneticEnergyDensity = 0.5*B_mag*B_mag/constants::PERMEABILITY;
      out << "# Wave energy density is " << energy << " J/m3" << endl;
      out << "# Relative energy is " << energy/magneticEnergyDensity << endl;
      
      out.close();
      
      // Deallocate memory and exit:
      delete [] waveEnergy; waveEnergy = NULL;
      energySpectrum->finalize();
      delete energySpectrum; energySpectrum = NULL;
      return true;
   }
   
   Real scatterMillstein(SimulationClasses& simClasses,Real dt,Real mu,Real d_mu_max,
			 Real D_mumu,Real d_D_mumu_dmu) {
      Real gaussianRandom = simClasses.random.gaussian();
      Real mu_new = mu + sqrt(2.0*D_mumu*dt) * gaussianRandom 
	          + 0.5*d_D_mumu_dmu*(gaussianRandom*gaussianRandom+1)*dt;      
      int32_t counter=1;

      while (fabs(mu_new) > 1.0 || fabs(mu_new-mu) > d_mu_max) {
	 gaussianRandom = simClasses.random.gaussian();
	 mu_new = mu + sqrt(2.0*D_mumu*dt) * gaussianRandom
	        + 0.5*d_D_mumu_dmu*(gaussianRandom*gaussianRandom+1)*dt;
	 if (counter > 10) {
	    mu_new = mu;
	    break;
	 }
	 ++counter;
      }

      return mu_new;
   }

   /** Fetch lambda times intensity and its derivative with respect to 
    * lambda at the given phase-space position.
    * @param pos Phase-space position.
    * @param lambdaIntensity Array containing lambda times intensity values.
    * @param I_wave Lambda times intensity is written here.
    * @param dI_wave Derivative of lambda times intensity is written here.*/
   void fetchIntensities(const Real* pos,
			 const Real* lambdaIntensity,
			 Real& I_wave,Real& dI_wave) {
      #ifndef NDEBUG
      bool ok = true;
      if (pos[3] >=  1.0 && pos[3] <= simControl.N_wavelengthMeshCells+1) { }
      else ok = false;
      if (pos[0] >= 1.0 || pos[0] <= block::WIDTH_X+1) { } 
      else ok = false;
      if (pos[1] >= 1.0 || pos[1] <= block::WIDTH_Y+1) { } 
      else ok = false;
      if (pos[2] >= 1.0 || pos[2] <= block::WIDTH_Z+1) { } 
      else ok = false;
      
      if (ok == false) {
	 cerr << "(SEP PARTICLE SCATTERER) ERROR: Invalid position in fetchIntensities" << endl;
	 cerr << "\t pos: " << pos[0] << '\t' << pos[1] << '\t' << pos[2] << '\t' << pos[3] << endl;
	 exit(1);
      }
      #endif

      const uint32_t NWL = simControl.N_wavelengthMeshCells+2;
      int32_t indices[4];
      Real shapeFactors[12];
      switch (simControl.order) {
       case 0:
	 getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	 I_wave = interpolateScalarNGP_4D(lambdaIntensity,indices,shapeFactors,NWL);
	 dI_wave = 0.0;
	 break;
       case 1:
	 getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	 I_wave = interpolateScalarCIC_4D(lambdaIntensity,indices,shapeFactors,NWL);
	 shapeFactors[6] = -1.0;
	 shapeFactors[7] = +1.0;
	 dI_wave = interpolateScalarCIC_4D(lambdaIntensity,indices,shapeFactors,NWL);
	 break;
       case 2:
	 getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	 I_wave = interpolateScalarTSC_4D(lambdaIntensity,indices,shapeFactors,NWL);
	 shapeFactors[9]  = pos[3] - indices[3]-1;
	 shapeFactors[11] = pos[3] - indices[3]  ;
	 shapeFactors[10] = -shapeFactors[9]-shapeFactors[11];
	 dI_wave = interpolateScalarTSC_4D(lambdaIntensity,indices,shapeFactors,NWL);
	 break;
       default:
	 I_wave = 0.0;
	 dI_wave = 0.0;
	 break;
      }
      
      #ifndef NDEBUG
      if (I_wave == -numeric_limits<Real>::infinity() || I_wave == +numeric_limits<Real>::infinity()) ok = false;
      if (I_wave < 0.0) ok = false;
      if (I_wave != I_wave) ok = false;
      if (dI_wave == -numeric_limits<Real>::infinity() || dI_wave == +numeric_limits<Real>::infinity()) ok = false;
      if (dI_wave != dI_wave) ok = false;
      if (ok == false) {
	 cerr << "(SEP PARTICLE SCATTERER) ERROR: Invalid position in fetchIntensities" << endl;
	 cerr << "\t I_wave: " << I_wave << "\t dI_wave: " << dI_wave << endl;
	 exit(1);
      }
      #endif
   }
   
   Real scatterIsotropic(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			 int32_t* indices,Real* shapeFactors,
			 Real* RESTRICT pos,const Real& mu,const Real& omega,
			 const Real& speed,const Real& B_mag,const Real& mu_min,const Real& mu_max) {
      const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);
      const Real maxResonantWavelength = 2*M_PI/omega*speed;

      Real dt1 = dt / scatteringSubsteps;
      Real mu_new = mu;
      for (uint32_t i=0; i<scatteringSubsteps; ++i) {
	 Real A0,B0;
	 evaluateDiffusionCoefficientsIsotropic(B_mag,mu_new,maxResonantWavelength,diffConstant,pos,indices,shapeFactors,lambdaIntensity,B0,A0);

	 Real theta = 0.0;
	 do {                
	    theta = random.uniform();
	 } while (theta == 0.0);
	 theta = -1.0 * log(theta);
	 theta = sqrt(2*B0*dt1*theta);
	 
	 Real phi = 2*M_PI*random.uniform();
	 
	 mu_new = mu_new*cos(theta) + sqrt(1-mu_new*mu_new)*sin(theta)*cos(phi);
	 
	 if (mu_new < mu_min) mu_new = 2*mu_min - mu_new;
	 if (mu_new > mu_max) mu_new = 2*mu_max - mu_new;
      }
      return mu_new;
   }
   
   Real scatterIsotropicAnalytic(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
				 int32_t* indices,Real* shapeFactors,
				 Real* RESTRICT pos,const Real& mu,const Real& omega,
				 const Real& speed,const Real& B_mag,const Real& mu_min,const Real& mu_max) {
      const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);
      const Real maxResonantWavelength = 2*M_PI/omega*speed;

      Real dt1 = dt / scatteringSubsteps;
      Real mu_new = mu;
      for (uint32_t i=0; i<scatteringSubsteps; ++i) {
	 Real A0,B0;
	 evaluateDiffusionCoefficientsIsotropicAnalytic(B_mag,mu_new,maxResonantWavelength,diffConstant,pos,lambdaIntensity,B0,A0);
	 
	 Real theta = 0.0;
	 do {
	    theta = random.uniform();
	 } while (theta == 0.0);
	 theta = -1.0 * log(theta);
	 theta = sqrt(2*B0*dt1*theta);
	 
	 Real phi = 2*M_PI*random.uniform();
	 
	 mu_new = mu_new*cos(theta) + sqrt(1-mu_new*mu_new)*sin(theta)*cos(phi);
	 
	 if (mu_new < mu_min) mu_new = 2*mu_min - mu_new;
	 if (mu_new > mu_max) mu_new = 2*mu_max - mu_new;
      }
      return mu_new;
   }

   Real scatterPredCorr1(RandomNumber& random,const Real* RESTRICT lambdaIntensity,const Real& dt,
			 int32_t* indices,Real* shapeFactors,
			 Real* RESTRICT pos,const Real& mu,const Real& omega,
			 const Real& speed,const Real& B_mag,const Real& mu_min,const Real& mu_max) {
      const Real diffConstant = 0.5*M_PI*omega/(B_mag*B_mag);
      const Real eta = 0.5;
      const Real theta = 0.0;
      
      Real dt1 = dt / scatteringSubsteps;
      Real mu_trial = mu;
      Real gaussian;
      const Real maxResonantWavelength = 2*M_PI/omega*speed;
      
      for (uint32_t i=0; i<scatteringSubsteps; ++i) {
	 // Store old position:
	 Real mu_old = mu_trial;

	 // Predictor step:
	 Real A0,B0;
	 evaluateDiffusionCoefficients(mu_old,maxResonantWavelength,diffConstant,pos,indices,shapeFactors,lambdaIntensity,B0,A0);
	 Real A0_trial = (1-eta)*A0;

	 gaussian = random.gaussian();

	 Real delta_mu = A0*dt1 + sqrt(2*B0*dt1)*gaussian;
	 mu_trial = mu_old + delta_mu;
	 #ifndef NDEBUG
	    if (mu_trial != mu_trial) {
	       cerr << "(SEP PARTICLE SCATTERER) ERROR: NAN mu in PredCorr1 predictor step" << endl;
	       cerr << "mu(old,trial)          : " << mu_old << '\t' << mu_trial << endl;
	       cerr << "diff,drift terms       : " << B0 << '\t' << A0 << endl;
	       exit(1);
	    }
	 #endif
	 
	 // Apply reflecting boundary conditions:
	 if (mu_trial < mu_min) mu_trial = 2*mu_min - mu_trial;
	 if (mu_trial > mu_max) mu_trial = 2*mu_max - mu_trial;
	 
	 #ifndef NDEBUG
	 if (mu_trial < mu_min || mu_trial > mu_max) {
	    cerr << "invalid predictor pos " << pos[3] << " substep " << i << '\t' << mu_trial << ' ' << mu_min << ' ' << mu_max << endl;
	    cerr << "\t old mu_trial was " << 2*mu_min - mu_trial << " or " << 2*mu_max - mu_trial << endl;
	    cerr << "\t mu_old " << mu_old << " mode " << simControl.alfvenSign << endl;
	    cerr << "\t a,b " << A0*dt1 << '\t' << sqrt(2*B0*dt1)*gaussian << endl;
	    cerr << "\t d_mu " << delta_mu << "\t gaussian " << gaussian << endl;
	    exit(1);
	 }
         #endif

	 // Corrector step:
	 Real A1,B1;
	 evaluateDiffusionCoefficients(mu_trial,maxResonantWavelength,diffConstant,pos,indices,shapeFactors,lambdaIntensity,B1,A1);
	 const Real A1_trial = (1-eta)*A1;
	 
	 delta_mu = (theta*A1_trial + (1-theta)*A0_trial)*dt1
	          + (eta*sqrt(B1) + (1-eta)*sqrt(B0))*sqrt(2*dt1)*gaussian;
	 mu_trial = mu_old + delta_mu;

	 #ifndef NDEBUG
	    if (mu_trial != mu_trial) {
	       cerr << "(SEP PARTICLE SCATTERER) ERROR: NAN mu in PredCorr1 corrector step" << endl;
	       cerr << "mu(old,trial)          : " << mu_old << '\t' << mu_trial << endl;
	       cerr << "diff,drift terms       : " << B1 << '\t' << A1 << endl;
	       exit(1);
	 }
         #endif
	 
	 // Apply reflecting boundary conditions:
	 if (mu_trial < mu_min) mu_trial = 2*mu_min - mu_trial;
	 if (mu_trial > mu_max) mu_trial = 2*mu_max - mu_trial;

	 #ifndef NDEBUG
	 if (mu_trial < mu_min || mu_trial > mu_max) {
	    cerr << "invalid corrector pos " << pos[3] << " substep " << i << '\t' << mu_trial << ' ' <<mu_min << ' ' << mu_max << endl;
	    cerr << "\t old mu_trial was " << 2*mu_min - mu_trial << " or " << 2*mu_max - mu_trial << endl;
	    cerr << "\t mu_old " << mu_old << " mode " << simControl.alfvenSign << " substeps " << scatteringSubsteps << endl;
	    cerr << "\t a,b " << (theta*A1_trial + (1-theta)*A0_trial)*dt1 << '\t' << (eta*sqrt(B1) + (1-eta)*sqrt(B0))*sqrt(2*dt1)*gaussian << endl;
	    cerr << "\t d_mu " << delta_mu << "\t gaussian " << gaussian << endl;
	    exit(1);
	 }
	 #endif
      }
      
      return mu_trial;
   }

   bool scatterParticles(const std::vector<pargrid::CellID>& blockLIDs,Simulation& sim,SimulationClasses& simClasses,
			 ParticleListBase* particleList,sep::getStateFunction getState) {
      bool success = true;
      
      // Set scatterer function:
      if (particleScatterer == NULL) {
	 string scattererName;
	 corsair::getObjectWrapper().configReader.add("SEP.scatterer","Name of scattering function (string)",string("Anisotropic"));
	 corsair::getObjectWrapper().configReader.parse();
	 corsair::getObjectWrapper().configReader.get("SEP.scatterer",scattererName);
	 
	 if (scattererName == "Anisotropic") {
	    particleScatterer = scatterPredCorr1;
	 } else if (scattererName == "Isotropic") {
	    particleScatterer = scatterIsotropic;
	 } else if (scattererName == "IsotropicAnalytic") {
	    particleScatterer = scatterIsotropicAnalytic;
	 } else {
	    corsair::getObjectWrapper().simClasses.logger << "(SEP PARTICLE SCATTERER) ERROR: Could not create scatterer" << endl << write;
	    return false;
	 }
      }

      // Get spectral wave energy density array:
      //const Real* waveEnergyTemporary = simClasses.pargrid.getUserDataStatic<Real>(simControl.temporaryWaveEnergyDataID);

      // Get wave growth factor array:
      Real* waveGrowthGlobal = NULL;
      if (simControl.alfvenSign < 0) {
	 waveGrowthGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.antiparWaveGrowthDataID);
      } else {
	 waveGrowthGlobal = simClasses.pargrid.getUserDataStatic<Real>(simControl.parWaveGrowthDataID);
      }
      if (waveGrowthGlobal == NULL) return false;
      
      // Get array containing min,max wavelength mesh cell
      // indices that are good for scattering:
      typedef wmesh::validDatatype validInt;
      validInt* validBinsGlobal = NULL;
      if (simControl.alfvenSign < 0) {
	 validBinsGlobal = simClasses.pargrid.getUserDataStatic<validInt>(simControl.maxAntiparValidWavelengthBinsDataID);
      } else {
	 validBinsGlobal = simClasses.pargrid.getUserDataStatic<validInt>(simControl.maxParValidWavelengthBinsDataID);
      }
      #ifndef NDEBUG
         if (validBinsGlobal == NULL) {
	    simClasses.logger << "(SEP PARTICLE SCATTERER) ERROR: validBinsGlobal is NULL" << endl << write;
	    exit(1);
	 }
      #endif

      // Get particle array:
      pargrid::DataID dataID;
      if (particleList->getParticles(dataID) == false) return false;
      pargrid::DataWrapper<PARTICLE> wrapper = simClasses.pargrid.getUserDataDynamic<PARTICLE>(dataID);
      if (wrapper.valid() == false) return false;

      // Get particle species:
      const sep::Species* species = reinterpret_cast<const sep::Species*>(particleList->getSpecies());
      Q_per_M = species->q_per_m;

      const int NWL  = simControl.N_wavelengthMeshCells+2;
      const int SIZE = block::SIZE_PLUS_ONE_LAYER*NWL;
      Real* lambdaIntensity = new Real[SIZE];
      Real* waveGrowthTmp = new Real[SIZE];
      //Real* phaseSpaceVolumes = new Real[SIZE];

      Real t_propag = 0.0;
      for (size_t b=0; b<blockLIDs.size(); ++b) {
	 // Skip blocks with zero particles:
	 if (wrapper.size(blockLIDs[b]) == 0) continue;

	 // Calculate block indices:
	 const pargrid::CellID blockLID = blockLIDs[b];
	 const pargrid::CellID blockGID = simClasses.pargrid.getGlobalIDs()[blockLID];
	 uint32_t blockIndices[3];
	 block::calculateBlockIndices(sim,blockGID,blockIndices[0],blockIndices[1],blockIndices[2]);
	 
	 // TEST
	 //#warning REMOVE ME
	 //cellInfo.i_block = blockIndices[0];
	 //cellInfo.j_block = blockIndices[1];
	 //cellInfo.k_block = blockIndices[2];
	 // END TEST

	 // If interpolations near shock are upwinded/downwinded, classify cells
	 // according to their shock regions:
	 bool shockedBlock = false;
	 int shockRegions[block::WIDTH_X+2];
	 if (simControl.useShockUpwinding == true) {
	    shockedBlock = classifyShockedCells(sim.t,shockRegions,blockIndices);
	 }

	 validInt* validBins = validBinsGlobal + blockLID*block::SIZE*2;

	 // Skip blocks where wave intensity is zero:
	 if (validBins[0] >= validBins[1]) continue;

	 // Skip exterior blocks:
	 const uint32_t nbrFlags = simClasses.pargrid.getNeighbourFlags(blockLID);
	 if (nbrFlags != pargrid::ALL_NEIGHBOURS_EXIST) continue;

	 // Measure block accumulation time if we are testing for repartitioning:
	 if (sim.countPropagTime == true) t_propag = MPI_Wtime();

	 #ifndef NDEBUG
	 bool ok = true;
	 if (validBins[0] < 0 || validBins[0] > NWL/2-1) ok = false;
	 if (validBins[1] < NWL/2 || validBins[1] > NWL-1) ok = false;
	 if (ok == false) {
	    cerr << "GID#" << blockGID << '\t' << blockIndices[0] << ' ' << blockIndices[1] << ' ' << blockIndices[2] << endl;
	    cerr << "\t valid bins invalid: " << validBins[0] << ' ' << validBins[1] << endl;
	    exit(1);
	 }
         #endif

	 // Clear wave growth array:
	 for (int i=0; i<SIZE; ++i) waveGrowthTmp[i] = 0.0;
	 for (int i=0; i<SIZE; ++i) lambdaIntensity[i] = 0.0;

	 loadIntensity(blockLID,lambdaIntensity,simControl.alfvenSign);
	 /*
	 // Calculate mu0 * fabs(wavelength) / phase-space-volume to 
	 // array phaseSpaceVolumes, this is needed in accumulation of 
	 // wave growth factors:
	 Real* tmpCellVolumes = new Real[block::SIZE_PLUS_ONE_LAYER];
	 block::fetchValues3D(simClasses,blockLID,tmpCellVolumes,simControl.cellVolumes);
	 for (int k=0; k<block::WIDTH_Z+2; ++k) for (int j=0; j<block::WIDTH_Y+2; ++j) for (int i=0; i<block::WIDTH_X+2; ++i) {
	    const Real spatCellVolume = tmpCellVolumes[block::arrayIndex(i,j,k)];
	    
	    Real d_lambda = simControl.wavelengthMeshCellSizes[0];
	    phaseSpaceVolumes[block::arrayIndex(i,j,k)*NWL+0] = constants::PERMEABILITY / (d_lambda*spatCellVolume);
	    
	    for (int l=1; l<NWL-1; ++l) {
	       d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	       phaseSpaceVolumes[block::arrayIndex(i,j,k)*NWL+l] = constants::PERMEABILITY / (d_lambda*spatCellVolume);
	    }

	    d_lambda = simControl.wavelengthMeshCellSizes[NWL-2];
	    phaseSpaceVolumes[block::arrayIndex(i,j,k)*NWL+NWL-1] = constants::PERMEABILITY / (d_lambda*spatCellVolume);
	 }
	 delete [] tmpCellVolumes; tmpCellVolumes = NULL;

	 // Load wave energy density to array lambdaIntensity, and convert spectral energy 
	 // density to spectral intensity by multiplying values by mu0 / d_lambda:
	 loadScalar4D(&simClasses,blockLID,waveEnergyTemporary,lambdaIntensity,simControl.N_wavelengthMeshCells);
	 for (int32_t cell=0; cell<block::SIZE_PLUS_ONE_LAYER; ++cell) {
	    for (uint32_t l=1; l<simControl.N_wavelengthMeshCells+1; ++l) {
	       const Real d_lambda = simControl.wavelengthMeshCellSizes[l-1];
	       lambdaIntensity[cell*NWL+l] *= constants::PERMEABILITY/d_lambda;
	    }
	 }*/

	 #ifndef NDEBUG
	 for (int32_t cell=0; cell<block::SIZE_PLUS_ONE_LAYER; ++cell) {
	    for (int32_t l=0; l<NWL; ++l) {
	       if (lambdaIntensity[cell*NWL+l] < 0.0 || lambdaIntensity[cell*NWL+l] != lambdaIntensity[cell*NWL+l]) {
		  cerr << "STEP: " << sim.timestep;
		  cerr << " Invalid value in waveEnergyTmp GID#" << blockGID << " indices: " << blockIndices[0] << ' ';
		  cerr << blockIndices[1] << ' ' << blockIndices[2];
		  cerr << " value: " << lambdaIntensity[cell*NWL+l] << " in bin " << l << endl;
		  success = false;
	       }
	    }
	 }
	 if (success == false) exit(1);
	 #endif

	 // Evaluate number of substeps needed in scattering:
	 #if PROFILE_LEVEL > 0
	    static int profSubsteps=-1;
	    profile::start("Scattering Substeps",profSubsteps);
	 #endif

	 // TEST
	 vector<Real> velocities;
	 velocities.push_back(1.4e6); // 10  keV/amu
	 velocities.push_back(1.4e7); // 1   MeV/amu
	 velocities.push_back(1.4e8); // 100 MeV/amu
	 velocities.push_back(1.4e9); // speed of light

	 vector<uint32_t> substeps;
	 substeps.push_back(evaluateNumberOfSubsteps(sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[0]));
	 substeps.push_back(evaluateNumberOfSubsteps(sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[1]));
	 substeps.push_back(evaluateNumberOfSubsteps(sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[2]));
	 substeps.push_back(evaluateNumberOfSubsteps(sim,blockLID,blockIndices,validBins,lambdaIntensity,velocities[3]));
	 // END TEST

         #if PROFILE_LEVEL > 0
	    profile::stop(); // Scattering substeps
	 #endif
	 
	 Real B[3];
	 Real V_wave[3];
	 Real dV_wave;
	 
/*	 // TEST
	 if (sim.timestep == 1100 && simControl.alfvenSign > 0) {
	    if (blockIndices[0] >= 14 && blockIndices[0] < 25) {
	       int l = simControl.N_wavelengthMeshCells*2/3;
	       Real pos[4];
	       int N = 20;
	       Real dx = 1.0 / N;
	       for (int i=0; i<N; ++i) {
		  Real V_wave[3];
		  Real dV_wave;
		  Real V_alfven;
		  pos[0] = blockIndices[0] + i*dx;
		  pos[1] = 0.5;
		  pos[2] = 0.5;
		  (*simControl.fieldsGetState)(blockLID,sim.t,pos,B,V_wave,dV_wave,V_alfven,+1);
		  int32_t region = simControl.shock->getShockRegion(sim.t,pos);
		  
		  pos[0] = 1 + i*dx;
		  pos[1] = 1.5;
		  pos[2] = 1.5;
		  pos[3] = l;
		  const int32_t SIZE = simControl.N_wavelengthMeshCells+2;
		  int32_t indices[4];
		  Real shapeFactors[12];
		  getShapeFactorsCIC_4D(pos,indices,shapeFactors);
		  
		  if (shockedBlock == true) {
		     calculateUpwindedParticleInterpShapeFactors(shockedBlock,shockRegions,indices,shapeFactors,region,static_cast<int32_t>(pos[0]));
		  }

		  Real I_wave;
		  I_wave = 0.0;
		  for (int32_t k_off=0; k_off<2; ++k_off) for (int32_t j_off=0; j_off<2; ++j_off) for (int32_t i_off=0; i_off<2; ++i_off) {
		     const Real shapeFactor = shapeFactors[0+i_off] * shapeFactors[2+j_off] * shapeFactors[4+k_off];
		     I_wave += lambdaIntensity[block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*SIZE+l] * shapeFactor;
		  }

		  stringstream ss;
		  ss << blockIndices[0]+i*dx << '\t' << I_wave/1e-21 << '\t' << lambdaIntensity[block::arrayIndex(1,1,1)*SIZE+l]/1e-21 << '\t' << V_wave[0]/1000 << '\t' << indices[0] << endl;
		  cerr << ss.str();
	       }
	    }
	 }
	 // END TEST */

	 // Scatter particles:
	 PARTICLE* particles = wrapper.data()[blockLID];
	 for (pargrid::ArraySizetype p=0; p<wrapper.size(blockLID); ++p) {
	    #ifndef NDEBUG
	    bool ok = true;
	    if (particles[p].state[0] < blockIndices[0]*block::WIDTH_X || particles[p].state[0] >= (blockIndices[0]+1)*block::WIDTH_X) ok = false;
	    if (particles[p].state[1] < blockIndices[1]*block::WIDTH_Y || particles[p].state[1] >= (blockIndices[1]+1)*block::WIDTH_Y) ok = false;
	    if (particles[p].state[2] < blockIndices[2]*block::WIDTH_Z || particles[p].state[2] >= (blockIndices[2]+1)*block::WIDTH_Z) ok = false;
	    if (ok == false) {
	       cerr << "(SEP PARTICLE SCATTERER) ERROR: Invalid particle coordinates" << endl;
	       cerr << "\t pos          : " << particles[p].state[0] << '\t' << particles[p].state[1] << '\t' << particles[p].state[2] << endl;
	       cerr << "\t block indices: " << blockIndices[0] << ' ' << blockIndices[1] << ' ' << blockIndices[2] << endl;
	       exit(1);
	    }
	    #endif

	    // Calculate particle's resonant wavelength in logical coordinates:
	    Real V_alfven;
	    (*simControl.fieldsGetState)(blockLID,sim.t,particles[p].state,B,V_wave,dV_wave,V_alfven,simControl.alfvenSign);
	    const Real B_mag = vectorMagnitude<3>(B);
	    const Real waveSpeed = dotProduct<3>(V_wave,B)/B_mag;

	    // Calculate resonant wavelength in selected (antiparallel, parallel) Alfven wave rest frame:
	    Real lambda_res = 2.0*M_PI/species->q_per_m/B_mag * (particles[p].state[sep::particle::V_PAR] - waveSpeed);

	    // Calculate resonant wavelength in logical coordinates:
	    const Real logical_lambda_res = (*simControl.getLogicalWavelength)(lambda_res);

	    Real pos[4];
	    pos[0] = particles[p].state[sep::particle::XCRD] - blockIndices[0]*block::WIDTH_X + 1;
	    pos[1] = particles[p].state[sep::particle::YCRD] - blockIndices[1]*block::WIDTH_Y + 1;
	    pos[2] = particles[p].state[sep::particle::ZCRD] - blockIndices[2]*block::WIDTH_Z + 1;
	    pos[3] = logical_lambda_res + 1;
	    
	    // Transform to wave rest frame:
	    const Real omega = species->q_per_m * B_mag;
	    const Real v_par_dot  = particles[p].state[sep::particle::V_PAR] - waveSpeed;
	    const Real v_gyro_2   = 2.0*particles[p].state[sep::particle::MU]*B_mag / species->mass;
	    const Real v_dot      = sqrt(v_par_dot*v_par_dot + v_gyro_2);
	    const Real mu_dot_old = v_par_dot / v_dot;
	    const Real maxResonantWavelength = 2*M_PI/omega*v_dot;
	    
	    // Approximate number of required scattering substeps:
	    int bin = max(0,static_cast<int>(1+log10(v_dot / 1.4e6)));
	    scatteringSubsteps = substeps[bin];
	    if (scatteringSubsteps > simControl.currentMaxScatteringSubsteps) {
	       simControl.currentMaxScatteringSubsteps = scatteringSubsteps;
	    }

	    // Calculate min,max mu interval where particle can scatter:
	    const wmesh::validDatatype minValidBin = min(NWL-3,validBins[0]+2);
	    const wmesh::validDatatype maxValidBin = max(0    ,validBins[1]-1);
	    const Real minWavelength = simControl.wavelengthMeshNodeCoordinates[minValidBin];
	    const Real maxWavelength = simControl.wavelengthMeshNodeCoordinates[maxValidBin];
	    const Real minMu = max(-1.0,minWavelength/maxResonantWavelength);
	    const Real maxMu = min(+1.0,maxWavelength/maxResonantWavelength);
	    limitedPitchMin = max(minMu,-mu_limit);
	    limitedPitchMax = min(maxMu,+mu_limit);

	    // Skip particles that are outside valid pitch interval:
	    if (mu_dot_old < minMu || mu_dot_old > maxMu) continue;

	    #ifndef NDEBUG
	    if (lambda_res < minWavelength || lambda_res > maxWavelength) ok = false;
	    if (lambda_res < simControl.wavelengthMeshNodeCoordinates[0] || 
		lambda_res > simControl.wavelengthMeshNodeCoordinates[simControl.N_wavelengthMeshCells-1]) ok = false;
	    if (ok == false) {
	       cerr << "(SEP PARTICLE SCATTERER) ERROR: invalid resonant wavelength" << endl;
	       cerr << "\t lambda    :" << lambda_res << endl;
	       cerr << "\t min/max l : " << minWavelength << '\t' << maxWavelength << endl;
	       cerr << "\t mu        : " << mu_dot_old << endl;
	       cerr << "\t min/max mu: " << minMu << '\t' << maxMu << endl;
	       exit(1);
	    }
	    #endif

	    // Calculate shape factors, since particles do not move during scattering 
	    // (except in wavelength) these stay constant:
	    int32_t indices[4];
	    Real shapeFactors[12];
	    switch (simControl.order) {
	     case 0:
	       getShapeFactorsNGP_4D(pos,indices,shapeFactors);
	       break;       
	     case 1:
	       getShapeFactorsCIC_4D(pos,indices,shapeFactors);
	       break;
	     case 2:
	       getShapeFactorsTSC_4D(pos,indices,shapeFactors);
	       break;
	     default:
	       cerr << "ERROR: Unknown shape factors in in scatterParticles" << endl;
	       exit(1);
	       break;
	    }

	    // If upwinding/downwinding is used near the shock, modify shape factors
	    // so that particles in upstream side do not scatter off downstream waves
	    // and vice versa:
	    if (shockedBlock == true) {
	       const int32_t region = simControl.shock->getShockRegion(sim.t,particles[p].state);
	       calculateUpwindedParticleInterpShapeFactors(shockRegions,indices,shapeFactors,region,static_cast<int32_t>(pos[0]));
	    }
	    
	    // Scatter particle:
	    Real mu_dot_new = mu_dot_old;
	    mu_dot_new = (*particleScatterer)(simClasses.random,lambdaIntensity,
					      sim.dt,indices,shapeFactors,pos,mu_dot_old,omega,v_dot,B_mag,minMu,maxMu);

	    // Apply energy (parallel speed & magn. moment) change 
	    // due to scattering:
	    const Real v_par_new = v_dot*mu_dot_new + waveSpeed;
	    particles[p].state[sep::particle::V_PAR]  = v_par_new;
	    particles[p].state[sep::particle::MU]     = 0.5*species->mass*v_dot*v_dot*(1-mu_dot_new*mu_dot_new)/B_mag;

	    #ifndef NDEBUG
	       if (particles[p].state[sep::particle::V_PAR] != particles[p].state[sep::particle::V_PAR]) ok = false;
	       if (particles[p].state[sep::particle::V_PAR] ==  numeric_limits<Real>::infinity()) ok = false;
	       if (particles[p].state[sep::particle::V_PAR] == -numeric_limits<Real>::infinity()) ok = false;
	       if (particles[p].state[sep::particle::MU]    != particles[p].state[sep::particle::MU]) ok = false;
               if (particles[p].state[sep::particle::MU] ==  numeric_limits<Real>::infinity()) ok = false;
	       if (particles[p].state[sep::particle::MU] == -numeric_limits<Real>::infinity()) ok = false;
	       if (ok == false) {
		  cerr << "(SEP PARTICLE SCATTERER) ERROR: NAN parallel speed and/or magnetic moment" << endl;
		  cerr << "\t V_par(new)       : " << particles[p].state[sep::particle::V_PAR] << endl;
		  cerr << "\t magn.moment(new) : " << particles[p].state[sep::particle::MU] << endl;
		  cerr << "\t mu (old,new)     : " << mu_dot_old << '\t' << mu_dot_new << endl;
		  cerr << "\t mu (min,max)     : " << minMu << '\t' << maxMu << endl;
		  exit(1);
	       }

	       const Real particleEnergyChange = species->mass*V_alfven*v_dot*(mu_dot_new-mu_dot_old);
	       simControl.dU_particles += particleEnergyChange*particles[p].state[sep::particle::WEIGHT];
	    #endif

	    // Apply wave energy change. EnergyChangeFactor is the total change 
	    // in wave energy divided by change in resonant particle wavelength:
	    Real energyChangeFactor = -0.5*species->mass*V_alfven*omega/M_PI * particles[p].state[particle::WEIGHT];

	    // Old and new resonant wavelengths:
	    const Real lambda_res_old = lambda_res;
	    const Real lambda_res_new = 2.0*M_PI/omega * v_dot * mu_dot_new;
	    simControl.dU_particles -= energyChangeFactor*(lambda_res_new-lambda_res_old);

	    // If interpolations are upwinded/downwinded near shock, modify 
	    // shape factors so that wave energy changes are accumulated to 
	    // opposite side of shock relative to particle position:
	    if (shockedBlock == true) {
	       indices[0] = static_cast<int32_t>(pos[0]-0.5);
	       calculateUpwindedParticleAccumShapeFactors(shockRegions,indices,shapeFactors,
							  simControl.shock->getShockRegion(sim.t,particles[p].state));
	    }

	    if (simControl.useUpwinding == true) {
	       switch (simControl.order) {
		case 0:
		  break;
		case 1:
		  indices[0] = static_cast<int32_t>(pos[0]);
		  indices[1] = static_cast<int32_t>(pos[1]);
		  indices[2] = static_cast<int32_t>(pos[2]);
		  shapeFactors[0] = 1;
		  shapeFactors[1] = 0;
		  shapeFactors[2] = 1;
		  shapeFactors[3] = 0;
		  shapeFactors[4] = 1;
		  shapeFactors[5] = 0;
		  break;
		case 2:
		  shapeFactors[0] = 0;
		  shapeFactors[1] = 1;
		  shapeFactors[2] = 0;
		  shapeFactors[3] = 0;
		  shapeFactors[4] = 1;
		  shapeFactors[5] = 0;
		  shapeFactors[6] = 0;		  
		  shapeFactors[7] = 1;
		  shapeFactors[8] = 0;
		  break;
		default:
		  cerr << "ERROR: Unknown order of accuracy in scatterParticles" << endl;
		  exit(1);
		  break;
	       }
	    }

	    // Accumulate wave energy growth factors:
	    accumulateWaveGrowthFactors(blockLID,&simClasses,indices,shapeFactors,pos,energyChangeFactor,
					lambda_res_old,lambda_res_new,waveGrowthTmp);
	 }
	 
	 // Add wave energy growth factors to global memory:
	 addValues4D(&simClasses,blockLID,waveGrowthTmp,waveGrowthGlobal,simControl.N_wavelengthMeshCells);

	 // Store block propagation time:
	 if (sim.countPropagTime == true) {
	    t_propag = std::max(0.0,MPI_Wtime() - t_propag);
	    simClasses.pargrid.getCellWeights()[blockLID] += t_propag;
	 }
      }
            
      delete [] lambdaIntensity; lambdaIntensity = NULL;
      delete [] waveGrowthTmp; waveGrowthTmp = NULL;
      return success;
   }
   
   /** Reset values of limitedPitchMin, limitedPitchMax to their default values. 
    * This function is used by some DataOperators that call evaluateDiffusionCoefficient 
    * function.*/
   void resetPitchLimits() {
      limitedPitchMin = -mu_limit;
      limitedPitchMax = +mu_limit;
   }

   bool writeIsotropicity(size_t N_particles,size_t N_bins,
			  uint32_t* probabMatrix,const std::string& name,int propagationDirection) {

      string fname;
      if (propagationDirection < 0) {
	 fname = "test_" + name + "_isotropicity_antipar.txt";
      } else {
	 fname = "test_" + name + "_isotropicity_par.txt";
      }
      fstream out(fname.c_str(), fstream::out);
      if (out.good() == false) return false;
      
      const Real d_mu = 2.0 / N_bins;
      size_t zeroBins = 0;
      Real error = 0;
      Real maxError = -numeric_limits<Real>::max();
      for (size_t j=0; j<N_bins; ++j) {
	 Real isotropicity = 0.0;
	 for (size_t i=0; i<N_bins; ++i) {
	    isotropicity += probabMatrix[i*N_bins+j];
	    if (i == j) continue;
	    
	    Real value1 = probabMatrix[j*N_bins+i];
	    Real value2 = probabMatrix[i*N_bins+j];
	    Real binError = 0.5*(fabs(value1-value2))/N_particles;
	    
	    if (binError > maxError) maxError = binError;	    
	    if (binError < 1e-4) ++zeroBins;
	    error += binError;
	 }
	 isotropicity = (isotropicity - N_particles)/N_particles;

	 out << -1.0 + j*d_mu << '\t' << isotropicity << endl;
	 out << -1.0 + (j+1)*d_mu << '\t' << isotropicity << endl;
      }
      error /= (N_bins*(N_bins-1));
      out.close();

      cout << name << " error: " << error << " max error: " << maxError;
      cout << " % zero bins: " << (1.0*zeroBins)/(N_bins*(N_bins-1)) << endl;
      return true;
   }
   
   bool writeTransitionProbabilities(size_t N_particles,size_t N_bins,const std::vector<size_t>& bins,
				     uint32_t* probabMatrix,const std::string& name,int propagationDirection) {
      string fname;
      if (propagationDirection < 0) {
	 fname = "test_" + name + "_probabs_antipar.txt";
      } else {
	 fname = "test_" + name + "_probabs_par.txt";
      }
      
      fstream out(fname.c_str(), fstream::out);
      if (out.good() == false) return false;
      
      const Real d_mu = 2.0 / N_bins;
      vector<vector<Real> > values(bins.size());
      for (size_t k=0; k<bins.size(); ++k) {
	 values[k].resize(N_bins);
	 
	 const size_t sourceBin = bins[k];
	 out << "# column " << k+1 << " source mu: " << -1.0 + (sourceBin+0.5)*d_mu << endl;
	 for (size_t i=0; i<N_bins; ++i) {
	    values[k][i] = (1.0*probabMatrix[sourceBin*N_bins+i])/N_particles;
	 }
      }
      
      for (size_t i=0; i<N_bins; ++i) {
	 out << -1.0 + i*d_mu << '\t';
	 for (size_t k=0; k<bins.size(); ++k) {
	    out << values[k][i] << '\t';
	 }
	 out << endl;
	 
	 out << -1.0 + (i+1)*d_mu << '\t';
	 for (size_t k=0; k<bins.size(); ++k) {
	    out << values[k][i] << '\t';
	 }
	 out << endl;
      }
      out.close();
      
      return true;
   }
   
} // namespace sep
