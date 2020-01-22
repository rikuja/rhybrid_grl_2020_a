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

#ifndef SEP_SHOCK_DRIFT_ACCELERATION_H
#define SEP_SHOCK_DRIFT_ACCELERATION_H

#include <cmath>
#include <definitions.h>
#include <linear_algebra.h>
#include <constants.h>

#include "sep_particle_definition.h"

namespace sep {
   template<class PARTICLE>
   void applyShockDriftAcceleration(PARTICLE& particle,Real mass,Real rho_mass,Real* RESTRICT V1_plasma_SIM,Real* RESTRICT B1,
				    Real* RESTRICT V1_shock,Real* RESTRICT shockNormal,
				    Real gasCompressionRatio);

   template<typename REAL>
   REAL getLossConeCosine(const REAL* B1,const REAL* shockNormal,REAL R_magnetic);
   /*
   Real getLossConeCosine(Real rho1_mass,const Real* V1_plasma_sim,const Real* B1,
			  const Real* shockNormal,Real R_magnetic,
			  const Real* V1_trans_SIM_2_HT);
   */
   Real getTransformLocalSNIFToLocalHT(const Real& B1_norm,const Real& B1_tang,const Real V1_norm);
   void getTransformLocalSNIFToLocalHT(const Real* V1_plasma_SNIF,const Real* B1,const Real* shockNormal,Real* V1_trans_SNIF_2_HT);

   void getTransformSimToLocalHT(const Real* V1_plasma_SIM,const Real* B1,
				 const Real* V1_shock,const Real* shockNormal,
				 Real* V1_trans_SIM_2_HT);

   void getTransformSimToLocalSNIF(const Real* V1_plasma_SIM,const Real* V1_shock,const Real* shockNormal,Real* V1_trans_SIM_2_SNIF);

   namespace oblique {
      template<typename REAL> 
      REAL getTangentialMagneticCompressionRatio(REAL cos_psi1,REAL alfvenMach2,REAL gasCompressionRatio);
      bool runTests();
      Real solveGasCompressionRatio(Real cos_psi1,Real sonicMach2,Real alfvenMach2,Real gamma);
   }
   
   namespace perpendicular {
      template<typename REAL>
      REAL getTangentialMagneticCompressionRatio(REAL cos_psi1,REAL alfvenMach2,REAL gasCompressionRatio);
      Real solveGasCompressionRatio(Real cos_psi1,Real sonicMach2,Real alfvenMach2,Real gamma);
   }
   
   /**
    * @param particle Particle hitting the shock.
    * @param mass Mass of particle hitting the shock.
    * @param rho1_mass Ambient plasma mass density in upstream region.
    * @param V1_plasma_SIM Ambient plasma velocity in upstream region in simulation frame.
    * @param B1 Magnetic field in upstream region.
    * @param V1_shock Local shock velocity in upstream region in simulation frame.
    * @param shockNormal Local shock normal.
    * @param gasCompressionRatio Gas compression ratio at impact point.*/
   template<class PARTICLE> inline
   void applyShockDriftAcceleration(PARTICLE& particle,Real mass,Real rho1_mass,Real* RESTRICT V1_plasma_SIM,
				    Real* RESTRICT B1,
				    Real* RESTRICT V1_shock,Real* RESTRICT shockNormal,
				    Real gasCompressionRatio) {
      // Calculate loss cone:
      const Real B1_norm_mag = dotProduct<3>(B1,shockNormal);
      Real B1_tang[3];
      for (int i=0; i<3; ++i) B1_tang[i] = B1[i] - B1_norm_mag * shockNormal[i];
      const Real B1_tang_mag = vectorMagnitude<3>(B1_tang);
      const Real B2_tang_mag = gasCompressionRatio * B1_tang_mag;
      const Real B1_mag2 = vectorMagnitude2<3>(B1);
      const Real B2_mag2 = B1_norm_mag*B1_norm_mag + B2_tang_mag*B2_tang_mag;
      const Real cos_theta_lc = sqrt(std::max(0.0,1.0 - B1_mag2/B2_mag2));

      // Get transformation velocity from simulation frame to local HT frame:
      Real V_transform[3];
      getTransformSimToLocalHT(V1_plasma_SIM,B1,V1_shock,shockNormal,V_transform);
      
      // Calculate parallel component of transformation:
      const Real B1_mag = sqrt(B1_mag2);
      const Real V_trans_par = dotProduct<3>(V_transform,B1) / B1_mag;

      // Calculate Alfvenic Mach number in local HT frame:
      Real V1_plasma_HT[3];
      for (int i=0; i<3; ++i) V1_plasma_HT[i] = V1_plasma_SIM[i] - V_transform[i];
      const Real V1_plasma_HT_mag = vectorMagnitude<3>(V1_plasma_HT);
      const Real alfvenicMach1 = V1_plasma_HT_mag / B1_mag * sqrt(constants::PERMEABILITY * rho1_mass);

      const Real magneticCompressionRatio = 
	  (alfvenicMach1*alfvenicMach1 - 1.0)
	/ (alfvenicMach1*alfvenicMach1 - gasCompressionRatio) 
	* gasCompressionRatio;

      // Calculate GC parallel speed and pitch in local HT frame:
      const Real V1_GC_HT_par = particle.state[sep::particle::V_PAR] - V_trans_par;
      const Real V1_gyro2 = 2.0*particle.state[sep::particle::MU]*B1_mag / mass;
      const Real V1_GC_HT_mag = sqrt(V1_GC_HT_par*V1_GC_HT_par + V1_gyro2);
      const Real cos_theta1_HT = V1_GC_HT_par / V1_GC_HT_mag;
      
      if (B1_norm_mag < 0.0) {
	 if (cos_theta1_HT > cos_theta_lc) {
	    // Particle is in loss cone:
	    
	    //std::cerr << "trans cos_theta: " << cos_theta1_HT << "\t cos_theta_lc: " << cos_theta_lc << std::endl;
	    //std::cerr << "\t V_par: " << V1_GC_HT_par/1000.0 << "\t V_mag: " << V1_GC_HT_mag/1000.0 << std::endl;
	    //std::cerr << "\t B1: " << B1_norm_mag << '\t' << B1_tang_mag << std::endl;
	    //std::cerr << "\t B2: " << B1_norm_mag << '\t' << B2_tang_mag << std::endl;
	    //std::cerr << "\t R: " << gasCompressionRatio << std::endl;
	 } else {
	    // Particle is reflected:
	    
	    //std::cerr << "refl cos_theta: " << cos_theta1_HT << "\t cos_theta_lc: " << cos_theta_lc << std::endl;
	    //std::cerr << "\t V_par: " << V1_GC_HT_par/1000.0 << "\t V_mag: " << V1_GC_HT_mag/1000.0 << std::endl;
	    //std::cerr << "\t B1: " << B1_norm_mag << '\t' << B1_tang_mag << std::endl;
	    //std::cerr << "\t B2: " << B1_norm_mag << '\t' << B2_tang_mag << std::endl;
	    //std::cerr << "\t R: " << gasCompressionRatio << std::endl;
	    //std::cerr << "\t Alfven Mach: " << alfvenicMach1 << "\t R_gas: " << gasCompressionRatio << "\t R_magn: " << magneticCompressionRatio << std::endl;
	    particle.state[sep::particle::V_PAR] = -V1_GC_HT_par + V_trans_par;
	 }
      } else {
	 if (cos_theta1_HT < -cos_theta_lc) {
	    // Particle is in loss cone:
	    
	    //std::cerr << "trans cos_theta: " << cos_theta1_HT << "\t cos_theta_lc: " << -cos_theta_lc << std::endl;
	    //std::cerr << "\t V_par: " << V1_GC_HT_par/1000.0 << "\t V_mag: " << V1_GC_HT_mag/1000.0 << std::endl;
	    //std::cerr << "\t B1: " << B1_norm_mag << '\t' << B1_tang_mag << std::endl;
	    //std::cerr << "\t B2: " << B1_norm_mag << '\t' << B2_tang_mag << std::endl;
	    //std::cerr << "\t R: " << gasCompressionRatio << std::endl;
	 } else {
	    // Particle is reflected:
	    
	    //std::cerr << "refl cos_theta: " << cos_theta1_HT << "\t cos_theta_lc: " << -cos_theta_lc << std::endl;
	    //std::cerr << "\t V_par: " << V1_GC_HT_par/1000.0 << "\t V_mag: " << V1_GC_HT_mag/1000.0 << std::endl;
	    //std::cerr << "\t B1: " << B1_norm_mag << '\t' << B1_tang_mag << std::endl;
	    //std::cerr << "\t B2: " << B1_norm_mag << '\t' << B2_tang_mag << std::endl;
	    //std::cerr << "\t R: " << gasCompressionRatio << std::endl;
	    //std::cerr << "\t Alfven Mach: " << alfvenicMach1 << "\t R_gas: " << gasCompressionRatio << "\t R_magn: " << magneticCompressionRatio << std::endl;
	    particle.state[sep::particle::V_PAR] = -V1_GC_HT_par + V_trans_par;
	 }
      }
   }

   template<typename REAL> inline
   REAL getLossConeCosine(const REAL* RESTRICT B1,const REAL* RESTRICT shockNormal,REAL R_magnetic) {
      // Calculate loss cone cosine:
      const REAL cos_theta_lc = sqrt(std::max(0.0,1.0 - 1/R_magnetic));
      
      // Return loss cone with correct sign:
      const REAL B1_norm_mag = dotProduct<3>(B1,shockNormal);
      if (B1_norm_mag <= 0.0) return cos_theta_lc;
      else return -cos_theta_lc;
   }
   
   template<typename REAL> inline
   REAL oblique::getTangentialMagneticCompressionRatio(REAL cos_psi1,REAL alfvenMach2,REAL gasCompressionRatio) {
      return (alfvenMach2 - cos_psi1*cos_psi1)/(alfvenMach2 - gasCompressionRatio*cos_psi1*cos_psi1)*gasCompressionRatio;
   }
   
   template<typename REAL> inline
   REAL perpendicular::getTangentialMagneticCompressionRatio(REAL cos_psi1,REAL alfvenMach2,REAL gasCompressionRatio) {
      return gasCompressionRatio;
   }
   
} // namespace sep

#endif
