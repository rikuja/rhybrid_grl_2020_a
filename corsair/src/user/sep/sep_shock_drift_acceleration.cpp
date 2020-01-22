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

#include <definitions.h>
#include <linear_algebra.h>
#include <root_finding.h>

#include "sep_shock_drift_acceleration.h"
#include "sep_simcontrol.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   
   static Real cos_psi2 = NAN;
   static Real sin_psi2 = NAN;
   static Real gamma = NAN;
   static Real alfvenMach2 = NAN;
   static Real sonicMach2 = NAN;
   static Real waveRatio;

   int sign(Real value) {
      if (value < 0.0) return -1;
      else if (value == 0.0) return 0;
      return +1;
   }
   
   Real gasCompressionEquationWavesParallel(Real R) {
      Real value = (alfvenMach2-R)*(alfvenMach2-R);
      value *= (2*R - sonicMach2*((gamma+1)-(gamma-1)*R));
      value += waveRatio*waveRatio * sonicMach2 * ((gamma+1)*R*R + (alfvenMach2*(2-gamma) - (gamma+1))*R + gamma*alfvenMach2);
      return -value;
   }

   Real cubicGasCompressionEquationOblique(Real R) {
      // Oblique
      Real value = (alfvenMach2-R*cos_psi2)*(alfvenMach2-R*cos_psi2);
      value *= (((gamma+1)-(gamma-1)*R)*sonicMach2 - 2*R);
      return value - sonicMach2*R*sin_psi2*((gamma+(2-gamma)*R)*alfvenMach2 - ((gamma+1)-(gamma-1)*R)*R*cos_psi2);
      //return value - alfvenMach2*R*sin_psi2*((gamma+(2-gamma)*R)*alfvenMach2 - ((gamma+1)-(gamma-1)*R)*R*cos_psi2)*sonicMach2/alfvenMach2;

      /*
      Real value = -((gamma-1)*sonicMach2 + 2*cos_psi2)*cos_psi2*R*R*R;
      value += (((gamma+1)*sonicMach2 + gamma*sonicMach2*alfvenMach2 + 4*alfvenMach2)*cos_psi2 + (gamma-2)*sonicMach2*alfvenMach2)*R*R;
      value -= (((gamma-1)*sonicMach2+2)*alfvenMach2 + gamma*sonicMach2 + (gamma+2)*sonicMach2*cos_psi2)*alfvenMach2*R;
      value += (gamma+1)*sonicMach2*alfvenMach2*alfvenMach2;
      */
      /*
      Real value = ((gamma+1)*sonicMach2 + gamma*sonicMach2*alfvenMach2 + 4*alfvenMach2)*R;
      value -= ((gamma-1)*sonicMach2 + 2*cos_psi2)*R*R;
      value -= (gamma+2)*sonicMach2*alfvenMach2;
      value *= R*cos_psi2;
      value += ((gamma-2)*sonicMach2*R*R - (((gamma-1)*sonicMach2+2)*alfvenMach2+gamma*sonicMach2)*R + (gamma+1)*sonicMach2*alfvenMach2)*alfvenMach2;
      */

      return value;
   }

   Real cubicGasCompressionEquationPerp(Real R) {
      // Perpendicular:
      const Real beta = 2/gamma*alfvenMach2/sonicMach2;
      Real value = 2*(2-gamma)*R*R;
      value += (2+2*beta + (gamma-1)*beta*sonicMach2)*gamma*R;
      return value - gamma*(gamma+1)*beta*sonicMach2;
   }

   bool oblique::runTests() {
      Real alfvenMach2 = 0;
      //Real alfvenMach2 = 1.705;
      //Real alfvenMach2 = 1.73168;
      
      alfvenMach2 *= alfvenMach2;
      //Real sonicMach2 = 6;
      Real sonicMach2 = 0;

      sep::alfvenMach2 = alfvenMach2;
      sep::sonicMach2 = sonicMach2;
      sep::gamma = 5.0/3.0;
      sep::cos_psi2 = 1.0;
      sep::sin_psi2 = 1.0 - sep::cos_psi2;
      
      fstream out("test_shock_adiabatic_oblique.txt",fstream::out);
      if (out.good() == false) return false;
      
      Real R_min = 0.0;
      Real R_max = 4.0;
      int N = 100;
      Real dR = (R_max-R_min)/(N-1);
      for (int i=0; i<N; ++i) {
	 Real R = R_min + i*dR;
	 out << R << '\t';

	 sep::cos_psi2 = 1.0; 
	 sep::sin_psi2 = 1.0 - sep::cos_psi2;
	 sep::alfvenMach2 = alfvenMach2 * sep::cos_psi2;
	 sep::sonicMach2  = sonicMach2  * sep::cos_psi2;
	 out << cubicGasCompressionEquationOblique(R) << '\t';

	 waveRatio = 1e-3;
	 out << gasCompressionEquationWavesParallel(R) << '\t';

	 /*
	 sep::cos_psi2 = 0.94*0.94;
	 sep::sin_psi2 = 1.0 - sep::cos_psi2;
	 sep::alfvenMach2 = alfvenMach2 * sep::cos_psi2;
	 sep::sonicMach2  = sonicMach2  * sep::cos_psi2;
	 out << cubicGasCompressionEquationOblique(R) << '\t';
	 
	 sep::cos_psi2 = 0.93*0.93;
	 sep::sin_psi2 = 1.0 - sep::cos_psi2;
	 sep::alfvenMach2 = alfvenMach2 * sep::cos_psi2;
	 sep::sonicMach2  = sonicMach2  * sep::cos_psi2;
	 out << cubicGasCompressionEquationOblique(R) << '\t';

	 sep::cos_psi2 = 0.92*0.92;
	 sep::sin_psi2 = 1.0 - sep::cos_psi2;
	 sep::alfvenMach2 = alfvenMach2 * sep::cos_psi2;
	 sep::sonicMach2  = sonicMach2  * sep::cos_psi2;
	 out << cubicGasCompressionEquationOblique(R) << '\t';
	 */
	 out << endl;
      }
      out.close();
      
      out.open("test_shock_gas_compression.txt",fstream::out);
      if (out.good() == false) return false;

      N = 1000;
      Real M_min = 1.0;
      Real M_max = 10.0;
      Real dM = (log10(M_max)-log10(M_min))/(N-1);
      for (int i=0; i<N; ++i) {
	 Real log_alfvenMach = log10(M_min) + i*dM;
	 Real alfvenMach = pow(10.0,log_alfvenMach);
	 Real alfvenMach2 = alfvenMach*alfvenMach;
	 Real cos_psi = 1.0;

	 out << alfvenMach << '\t';
	 out << oblique::solveGasCompressionRatio(cos_psi,sonicMach2,alfvenMach2,5.0/3.0) << '\t';
	 out << endl;
      }
      
      out.close();
      
      return true;
   }
   
   Real oblique::solveGasCompressionRatio(Real cos_psi1,Real sonicMach2,Real alfvenMach2,Real gamma) {
      // Return unit gas compression ratio if a shock cannot physically exist:
      if (alfvenMach2*sonicMach2 <= alfvenMach2 + sonicMach2) return 1.0;

      // Copy parameter values to local variables:
      sep::cos_psi2 = cos_psi1*cos_psi1;
      sep::sin_psi2 = 1.0 - cos_psi2;
      sep::sonicMach2 = sonicMach2;
      sep::alfvenMach2 = alfvenMach2;
      sep::gamma = gamma;

      return rootf::bisection(cubicGasCompressionEquationOblique,0.9,4.0,0.001,0.01,100);
   }
   
   Real perpendicular::solveGasCompressionRatio(Real cos_psi1,Real sonicMach2,Real alfvenMach2,Real gamma) {
      if (alfvenMach2*sonicMach2 <= alfvenMach2 + sonicMach2) return 1.0;
            
      // Copy parameter values to local variables:
      sep::cos_psi2 = cos_psi1*cos_psi1;
      sep::sin_psi2 = 1.0 - cos_psi2;
      sep::sonicMach2 = sonicMach2;
      sep::alfvenMach2 = alfvenMach2;
      sep::gamma = gamma;

      return rootf::bisection(cubicGasCompressionEquationPerp,0.0,4.0,0.01,0.01,100);
   }

   void getTransformLocalSNIFToLocalHT(const Real* RESTRICT V1_plasma_SNIF,const Real* RESTRICT B1,
				       const Real* RESTRICT shockNormal,Real* RESTRICT V1_trans_SNIF_2_HT) {
      // Calculate tangential B:
      Real B1_tang[3];
      const Real B1_norm_mag = dotProduct<3>(B1,shockNormal);
      for (int i=0; i<3; ++i) B1_tang[i] = B1[i] - B1_norm_mag*shockNormal[i];
      
      // Calculate tan of shock obliquity angle in upstream.
      // NOTE: We must ensure that tan_psi1 does not become infinite or nan:
      Real tan_psi1 = vectorMagnitude<3>(B1_tang) / B1_norm_mag;
      if (fabs(B1_norm_mag) < 1.0e-13) {
	 if (B1_norm_mag < 0.0) tan_psi1 = vectorMagnitude<3>(B1_tang) / -1e-13;
	 else tan_psi1 = vectorMagnitude<3>(B1_tang) / 1e-13;
      }
      
      // Calculate unit vector to direction of B_tang:
      Real B1_tang_unit[3];
      for (int i=0; i<3; ++i) B1_tang_unit[i] = B1_tang[i];
      unitVector<3>(B1_tang_unit);
      
      // Calculate transformation velocity from local SNIF to local HT frame:
      const Real V1_plasma_SRF_norm_mag = -dotProduct<3>(V1_plasma_SNIF,shockNormal);
      for (int i=0; i<3; ++i) {
	 const Real speed = V1_plasma_SRF_norm_mag * tan_psi1 * B1_tang_unit[i];
	 V1_trans_SNIF_2_HT[i] = min(speed,constants::SPEED_LIGHT);
	 V1_trans_SNIF_2_HT[i] = max(V1_trans_SNIF_2_HT[i],-constants::SPEED_LIGHT);
      }
   }
   
   Real getTransformLocalSNIFToLocalHT(const Real& B1_norm,const Real& B1_tang,const Real V1_norm) {
      // Calculate tan of shock obliquity angle in upstream.
      // NOTE: We must ensure that tan_psi1 does not become infinite or nan:
      Real tan_psi1 = fabs(B1_tang) / B1_norm;
      if (fabs(B1_norm) < 1.0e-13) {
	 if (B1_norm < 0.0) tan_psi1 = fabs(B1_tang) / -1e-13;
	 else tan_psi1 = fabs(B1_tang) / 1e-13;
      }
      
      Real speed = -V1_norm * tan_psi1;
      speed = min(speed,+constants::SPEED_LIGHT);
      speed = max(speed,-constants::SPEED_LIGHT);
      return speed*B1_tang / sqrt(B1_norm*B1_norm + B1_tang*B1_tang);
   }

   /** Calculate transformation velocity from simulation frame to local de Hofmann-Teller frame.
    * @param V1_plasma_SIM Plasma velocity in simulation frame (upstream region). 
    * @param B1 Magnetic field in simulation frame (upstream region).
    * @param V1_shock Local shock velocity (upstream region).
    * @param shockNormal Local shock normal.
    * @param V1_trans_SIM_2_HT Array where calculated transformation velocity is written to.*/
   void getTransformSimToLocalHT(const Real* RESTRICT V1_plasma_SIM,const Real* RESTRICT B1,
				 const Real* RESTRICT V1_shock,const Real* RESTRICT shockNormal,
				 Real* RESTRICT V1_trans_SIM_2_HT) {
      // NOTE about notation:
      // V1_plasma_SRF means that it is measured in upstream region (1 in V1),
      // it is related to plasma (_plasma_), and it is defined in 
      // local shock rest frame (SRF).
      
      // Calculate transform velocity from simulation frame to local SNIF:
      Real V1_trans_sim_2_SNIF[3];
      getTransformSimToLocalSNIF(V1_plasma_SIM,V1_shock,shockNormal,V1_trans_sim_2_SNIF);

      // *****     PART 2: CALCULATE TRANSFORM VELOCITY FROM     ***** //
      // ***** SIMULATION FRAME TO LOCAL DE HOFFMAN-TELLER FRAME ***** //
      
      // Calculate tangential B:
      Real B1_tang[3];
      const Real B1_norm_mag = dotProduct<3>(B1,shockNormal);
      for (int i=0; i<3; ++i) B1_tang[i] = B1[i] - B1_norm_mag*shockNormal[i];
      
      // Calculate tan of shock obliquity angle in upstream.
      // NOTE: We must ensure that tan_psi1 does not become infinite or nan:
      Real tan_psi1 = vectorMagnitude<3>(B1_tang) / B1_norm_mag;
      if (fabs(B1_norm_mag) < 1.0e-13) {
	 if (B1_norm_mag < 0.0) tan_psi1 = vectorMagnitude<3>(B1_tang) / -1e-13;
	 else tan_psi1 = vectorMagnitude<3>(B1_tang) / 1e-13;
      }
      
      // Calculate unit vector to direction of B_tang:
      Real B1_tang_unit[3];
      for (int i=0; i<3; ++i) B1_tang_unit[i] = B1_tang[i];
      unitVector<3>(B1_tang_unit);

      // Calculate transformation speed from simulation frame
      // to local de Hoffman-Teller frame:
      Real V1_plasma_SNIF[3];
      for (int i=0; i<3; ++i) V1_plasma_SNIF[i] = V1_plasma_SIM[i] - V1_trans_sim_2_SNIF[i];
      Real V1_plasma_SRF_norm_mag = -dotProduct<3>(V1_plasma_SNIF,shockNormal);
      //V1_plasma_SRF_norm_mag = fabs(V1_plasma_SRF_norm_mag);
      
      for (int i=0; i<3; ++i) {
	 const Real speed = V1_trans_sim_2_SNIF[i] + V1_plasma_SRF_norm_mag * tan_psi1 * B1_tang_unit[i];
	 V1_trans_SIM_2_HT[i] = min(speed,constants::SPEED_LIGHT);
	 V1_trans_SIM_2_HT[i] = max(V1_trans_SIM_2_HT[i],-constants::SPEED_LIGHT);
      }
   }

   void getTransformSimToLocalSNIF(const Real* V1_plasma_SIM,const Real* V1_shock,const Real* shockNormal,Real* V1_trans_SIM_2_SNIF) {
      // Calculate plasma velocity in local shock rest frame (SRF):
      Real V1_plasma_SRF[3];
      for (int i=0; i<3; ++i) V1_plasma_SRF[i] = V1_plasma_SIM[i] - V1_shock[i];
      
      // Calculate tangential component of plasma velocity in SRF:
      const Real V1_plasma_SRF_norm_mag = dotProduct<3>(V1_plasma_SRF,shockNormal);
      Real V1_plasma_SRF_tang[3];
      for (int i=0; i<3; ++i) V1_plasma_SRF_tang[i] = V1_plasma_SRF[i] - V1_plasma_SRF_norm_mag*shockNormal[i];
       
      // Calculate transform velocity from simulation frame to local SNIF:
      for (int i=0; i<3; ++i) V1_trans_SIM_2_SNIF[i] = V1_shock[i] + V1_plasma_SRF_tang[i];
   }
   
   
   
} // namespace sep
