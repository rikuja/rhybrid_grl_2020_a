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

#ifndef SEP_PARTICLE_SCATTERER_H
#define SEP_PARTICLE_SCATTERER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <particle_list_base.h>
#include <vector>

#include "sep_fields_container.h"

namespace sep {

   template<typename REAL>
   void evaluateDiffusionCoefficients(const REAL& mu,const Real maxLambda,const REAL& I,const REAL& dI,REAL& diff,REAL& drift);
   
   template<typename REAL>
   REAL fetchIntensityNGP_3D(const int32_t& SIZE,const int32_t& l,const REAL* intensity,const int32_t* indices,const REAL* shapeFactors);
   template<typename REAL>
   REAL fetchIntensityCIC_3D(const int32_t& SIZE,const int32_t& l,const REAL* intensity,const int32_t* indices,const REAL* shapeFactors);
   template<typename REAL>
   REAL fetchIntensityTSC_3D(const int32_t& SIZE,const int32_t& l,const REAL* intensity,const int32_t* indices,const REAL* shapeFactors);

   void evaluateDiffusionCoefficients(Real pitch,Real maxResonantLambda,Real diffConstant,Real* pos,int32_t* indices,
				      Real* shapeFactors,const Real* intensity,Real& D_mumu,Real& d_D_mumu);

   bool loadIntensity(pargrid::CellID blockLID,Real* array,int alfvenSign);
   
   bool runScattererTests(Simulation& sim,SimulationClasses& simClasses);

   Real scatterMillstein(SimulationClasses& simClasses,Real dt,Real mu,Real d_mu_max,
			 Real D_mumu,Real d_D_mumu_dmu);
   
   bool scatterParticles(const std::vector<pargrid::CellID>& blockLIDs,Simulation& sim,SimulationClasses& simClasses,
			 ParticleListBase* particleList,sep::getStateFunction getState);
   
   void resetPitchLimits();

   // ***** TEMPLATE FUNCTION DEFINITIONS ***** //
   
   /** Evaluate drift and diffusion coefficients that are used while solving 
    * pitch angle scattering stochastic differential equation. NOTE: "diffusion coefficient"
    * below equals pi/2 times omega/(B*B).
    * @param mu Particle pitch (cosine of pitch angle) in wave rest frame.
    * @param maxLambda Maximum resonant wavelength in wave rest frame, equals to 
    * 2*pi*speed/omega.
    * @param I Spectral wave intensity times diffusion coefficient at particle position.
    * @param dI Derivative of I with respect to wavelength times diffusion coefficient.
    * @param diff Calculated diffusion coefficient (output).
    * @param drift Calculated drift coefficient (output).*/
   template<typename REAL> inline
   void evaluateDiffusionCoefficients(const REAL& mu,const Real maxLambda,const REAL& I,const REAL& dI,REAL& diff,REAL& drift) {
      const Real fabs_mu = fabs(mu);
      diff  = (1-mu*mu)*maxLambda*fabs_mu*I;
      if (mu < 0.0) {
	 drift = maxLambda * ( -2*mu*fabs_mu*I
			       -(1-mu*mu)*I
			       +(1-mu*mu)*fabs_mu*maxLambda*dI
			     );
      } else {
	 drift = maxLambda * ( -2*mu*fabs_mu*I
			       +(1-mu*mu)*I
			       +(1-mu*mu)*fabs_mu*maxLambda*dI
			     );
      }
      
      #ifndef NDEBUG
      if (diff < 0.0 || diff != diff) {
	 std::cerr << "(SEP PARTICLE SCATTERER) ERROR: Negative diffusion coefficient" << std::endl;
	 std::cerr << "diff,drift terms       : " << diff << '\t' << drift << std::endl;
	 std::cerr << "mu                     : " << mu << std::endl;
	 std::cerr << "I,dI                   : " << I << '\t' << dI << std::endl;
	 std::cerr << "maxLambda              : " << maxLambda << std::endl;
	 exit(1);
      }
      #endif
   }

   template<typename REAL> inline
   REAL fetchIntensityNGP_3D(const int32_t& SIZE,const int32_t& l,const REAL* RESTRICT intensity,
			     const int32_t* RESTRICT indices,const REAL* RESTRICT shapeFactors) {
      return intensity[block::arrayIndex(indices[0],indices[1],indices[2])*SIZE+l];
   }

   template<typename REAL> inline
   REAL fetchIntensityCIC_3D(const int32_t& SIZE,const int32_t& l,const REAL* RESTRICT intensity,
			     const int32_t* RESTRICT indices,const REAL* RESTRICT shapeFactors) {
      REAL I_wave = 0.0;
      for (int32_t k_off=0; k_off<2; ++k_off) for (int32_t j_off=0; j_off<2; ++j_off) for (int32_t i_off=0; i_off<2; ++i_off) {
	 const Real shapeFactor = shapeFactors[0+i_off] * shapeFactors[2+j_off] * shapeFactors[4+k_off];
	 I_wave += intensity[block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*SIZE+l] * shapeFactor;
      }
      return I_wave;
   }
   
   template<typename REAL> inline
   REAL fetchIntensityTSC_3D(const int32_t& SIZE,const int32_t& l,const REAL* RESTRICT intensity,
			     const int32_t* RESTRICT indices,const REAL* RESTRICT shapeFactors) {
      REAL I_wave = 0.0;
      for (int32_t k_off=-1; k_off<2; ++k_off) for (int32_t j_off=-1; j_off<2; ++j_off) for (int32_t i_off=-1; i_off<2; ++i_off) {
	 const Real shapeFactor = shapeFactors[1+i_off] * shapeFactors[4+j_off] * shapeFactors[7+k_off];
	 I_wave += intensity[block::arrayIndex(indices[0]+i_off,indices[1]+j_off,indices[2]+k_off)*SIZE+l] * shapeFactor;
	 #ifndef NDEBUG
	 if (shapeFactor < 0 || shapeFactor > 1) {
	    std::cerr << "Invalid shape factor in fetchIntensityTSC_3D" << std::endl;
	    for (int i=0; i<9; ++i) {
	       std::cerr << '\t' << shapeFactors[i] << std::endl;
	    }
	    exit(1);
	 }
	 #endif
      }
      return I_wave;
   }

} // namespace sep

#endif
