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
#include <fstream>

#include <constants.h>
#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_shock_accelerator.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

namespace sep {

   static const string prefix = "Shock.accelerator";
   extern sep::SimControl simControl;
   
   
   ShockAccelerator::ShockAccelerator() { 
      sim = NULL;
      simClasses = NULL;
      phiFactor = NAN;
   }

   ShockAccelerator::~ShockAccelerator() { }

   void ShockAccelerator::borisBunemanAhead(Real* state,const Real* E,const Real* B,Real eforce) {
      const Real fact = 0.5*eforce;
      
      Real v_minus[3];
      v_minus[0] = state[3] + fact * E[0];
      v_minus[1] = state[4] + fact * E[1];
      v_minus[2] = state[5] + fact * E[2];
      
      Real t[3];
      t[0] = fact * B[0];
      t[1] = fact * B[1];
      t[2] = fact * B[2];
      const Real t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
      
      const Real fact2 = 2.0 / (1.0 + t_mag);
      Real s[3];
      s[0] = fact2 * t[0];
      s[1] = fact2 * t[1];
      s[2] = fact2 * t[2];

      Real v_dot[3];
      v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
      v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
      v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
      
      Real v_plus[3];
      v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
      v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
      v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
            
      state[3] = v_plus[0] + fact * E[0];
      state[4] = v_plus[1] + fact * E[1];
      state[5] = v_plus[2] + fact * E[2];
   }
   
   void ShockAccelerator::borisBunemanBehind(Real* state,const Real* E,const Real* B,Real eforce) {
      const Real fact = -0.5*eforce;

      Real v_minus[3];
      v_minus[0] = state[3] + fact * E[0];
      v_minus[1] = state[4] + fact * E[1];
      v_minus[2] = state[5] + fact * E[2];
      
      Real t[3];
      t[0] = fact * B[0];
      t[1] = fact * B[1];
      t[2] = fact * B[2];
      const Real t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
      
      const Real fact2 = 2.0 / (1.0 + t_mag);
      Real s[3];
      s[0] = fact2 * t[0];
      s[1] = fact2 * t[1];
      s[2] = fact2 * t[2];
      
      Real v_dot[3];
      v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
      v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
      v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
      
      Real v_plus[3];
      v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
      v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
      v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
      
      state[3] = v_plus[0] + fact * E[0];
      state[4] = v_plus[1] + fact * E[1];
      state[5] = v_plus[2] + fact * E[2];
   }
   
   void ShockAccelerator::calculateWebbResult(const PlasmaParameters& plasma,
					      const ParticleParameters& particle,
					      Real U_ini,Real& dU,Real& pitch) {
      Real B[3];
      if (particle.state[shockaccelerator::XPOS] < 0.0) {
	 B[0] = plasma.B1_norm;
	 B[1] = 0.0;
	 B[2] = plasma.B1_tang;
      } else {
	 B[0] = plasma.B1_norm;
	 B[1] = 0.0;
	 B[2] = plasma.R_magn * plasma.B1_tang;
      }
      const Real B_mag = vectorMagnitude<3>(B);
      
      Real E[3];
      E[0] = 0.0;
      E[1] = plasma.B1_tang*plasma.V1_plasma_norm - plasma.B1_norm*plasma.V1_plasma_tang;
      E[2] = 0.0;

      Real V_drifts[3];
      crossProduct(E,B,V_drifts);
      for (int i=0; i<3; ++i) V_drifts[i] /= (B_mag*B_mag);
      const Real U_drifts = 0.5*species.mass*vectorMagnitude2<3>(V_drifts);
      
      const Real U_final = 0.5*species.mass*particle.state[shockaccelerator::V_PAR]*particle.state[shockaccelerator::V_PAR]
	                 + particle.state[shockaccelerator::MU]*B_mag
	                 + U_drifts;
      
      const Real V_final = sqrt(particle.state[shockaccelerator::V_PAR]*particle.state[shockaccelerator::V_PAR]
				+ 2/species.mass*B_mag*particle.state[shockaccelerator::MU]);

      dU = (U_final-U_ini)/U_ini;
      pitch = acos(particle.state[shockaccelerator::V_PAR] / V_final) * 180.0/M_PI;
   }

   void ShockAccelerator::getDownstreamValues(const PlasmaParameters& plasma,Real V1_alfven,Real& V2_wave_par,Real* B2,Real& B2_mag,
					      Real* V2_wave,Real& tan_psi2) {
      // Upstream magnetic field:
      Real B1[3];
      B1[0] = plasma.B1_norm;
      B1[1] = 0.0;
      B1[2] = plasma.B1_tang;
      Real B1_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.B1_tang*plasma.B1_tang);
      Real tan_psi1 = B1[2] / B1[0];
      if (fabs(B1[0]) < 1.0e-13) {
	 if (B1[0] < 0.0) tan_psi1 = B1[2] / -1e-13;
	 else tan_psi1 = B1[2] / 1e-13;
      }

      // Downstream magnetic field:
      //Real B2[3];
      B2[0] = plasma.B1_norm;
      B2[1] = 0.0;
      B2[2] = plasma.B1_tang * plasma.R_magn;
      B2_mag = vectorMagnitude<3>(B2);
      tan_psi2 = B2[2] / B2[0];
      if (fabs(B2[0]) < 1.0e-13) {
	 if (B2[0] < 0.0) tan_psi2 = B2[2] / -1e-13;
	 else tan_psi2 = B2[2] / 1e-13;
      }

      // Upstream plasma velocity in SNIF:
      Real V1_plasma_SNIF[3];
      V1_plasma_SNIF[0] = plasma.V1_plasma_norm;
      V1_plasma_SNIF[1] = 0.0;
      V1_plasma_SNIF[2] = 0.0;

      // Downstream plasma velocity in HT:
      Real V2_plasma_SNIF[3];
      V2_plasma_SNIF[0] = plasma.V1_plasma_norm / plasma.R_gas;
      V2_plasma_SNIF[1] = 0.0;
      V2_plasma_SNIF[2] = V2_plasma_SNIF[0]*tan_psi2;
      
      // Downstream plasma velocity in SNIF:
      const Real V_HT1 = -plasma.V1_plasma_norm*tan_psi1;
      V2_plasma_SNIF[2] += V_HT1;
      
      //Real tmp[3];
      //crossProduct(V2_plasma_SNIF,B2,tmp);
      //cerr << "cp is " << vectorMagnitude<3>(tmp) << endl;
      //cerr << "V2_SNIF " << V2_plasma_SNIF[0] << ' ' << V2_plasma_SNIF[1] << ' ' << V2_plasma_SNIF[2] << endl;
      //cerr << "B2    " << B2[0] << ' ' << B2[1] << ' ' << B2[2] << endl;
      //exit(1);

      // Alfven speed in downstream:
      Real V2_alfven = B2_mag/B1_mag/sqrt(plasma.R_gas) * V1_alfven;
      //cerr << "V2_alfven " << V2_alfven << endl;
      /*
      // Shock normal:
      Real shockNormal[3];
      shockNormal[0] = -1.0;
      shockNormal[1] = 0.0;
      shockNormal[2] = 0.0;
      
      // Calculate transformation velocity from SNIF to HT, and 
      // components parallel to magnetic field:
      Real V_SNIF_2_HT[3];
      getTransformLocalSNIFToLocalHT(V1_plasma_SNIF,B1,shockNormal,V_SNIF_2_HT);
      const Real V1_SNIF_2_HT_par = dotProduct<3>(V_SNIF_2_HT,B1) / B1_mag;
      V2_SNIF_2_HT_par = V1_SNIF_2_HT_par * plasma.R_magn/plasma.R_gas;

      // Plasma velocity in downstream HT frame, parallel to B2:
      Real V2_plasma_HT[3];
      for (int i=0; i<3; ++i) V2_plasma_HT[i] = V2_plasma_SNIF[i] - V_SNIF_2_HT[i]*plasma.R_magn/plasma.R_gas;
      */
      // Wave speed in downstream HT frame:
      //Real V2_wave[3];
      //for (int i=0; i<3; ++i) V2_wave[i] = V2_plasma_HT[i] + V2_alfven*B2[i]/B2_mag;
      for (int i=0; i<3; ++i) V2_wave[i] = V2_plasma_SNIF[i];// + V2_alfven*B2[i]/B2_mag;
      V2_wave_par = dotProduct<3>(V2_wave,B2)/B2_mag;

      //Real tmp[3];
      //crossProduct(V2_plasma_HT,B2,tmp);
      
      //cerr << "V2_wave_HT is " << V2_wave[0] << '\t' << V2_wave[1] << '\t' << V2_wave[2] << endl;
      //cerr << "B1 is " << B1[0] << '\t' << B1[1] << '\t' << B1[2] << endl;
      //cerr << "B2 is " << B2[0] << '\t' << B2[1] << '\t' << B2[2] << endl;
      //unitVector<3>(B2);
      //cerr << "B2 is " << B2[0] << '\t' << B2[1] << '\t' << B2[2] << endl;
      
      //V2_wave[0] = V2_plasma_SNIF[0] + V2_alfven*B2[0]/B2_mag;
      //V2_wave[1] = V2_plasma_SNIF[1] + V2_alfven*B2[1]/B2_mag;
      //V2_wave[2] = V2_plasma_SNIF[2] + V2_alfven*B2[2]/B2_mag;
   }

   Real ShockAccelerator::getReturnedParticle(const PlasmaParameters& plasma,ParticleParameters& particle,const Real& V1_alfven,Real& pitch) {
      if (particle.state[shockaccelerator::XPOS] > 0.0) {
	 cerr << "ERROR: Particle initially in downstream" << endl;
	 exit(1);
      }

      // Particle state is given in upstream SNIF, 
      // solve shock encounter:
      solveShockEncounterAnalytic(plasma,particle);
       
      // If particle reflected, return immediately:
      if (particle.state[shockaccelerator::XPOS] < 0.0) {
	 return 1.0;
      }

      Real V1_plasma[3];
      V1_plasma[0] = plasma.V1_plasma_norm;
      V1_plasma[1] = 0;
      V1_plasma[2] = plasma.V1_plasma_tang;
      
      Real B1[3];
      B1[0] = plasma.B1_norm;
      B1[1] = 0;
      B1[2] = fabs(plasma.B1_tang);

      // Calculate upstream electric field and drift:
      Real E[3];
      crossProduct(B1,V1_plasma,E);

      Real V1_E[3];
      crossProduct(E,B1,V1_E);
      Real B1_mag = vectorMagnitude<3>(B1);
      for (int i=0; i<3; ++i) V1_E[i] /= (B1_mag*B1_mag);

      // Particle velocity vector in upstream SNIF (incl. drifts):
      //Real V_US_SNIF[3];
      //for (int i=0; i<3; ++i) V_US_SNIF[i] = particle.state[shockaccelerator::V_PAR]*B1[i]/B1_mag + V1_E[i];

      // Shock normal:
      Real shockNormal[3];
      shockNormal[0] = -1;
      shockNormal[1] = 0;
      shockNormal[2] = 0;
      
      // Transform velocity from upstream SNIF to HT frame:
      Real V_US_SNIF_2_HT[3];
      getTransformLocalSNIFToLocalHT(V1_plasma,B1,shockNormal,V_US_SNIF_2_HT);

      // Upstream plasma velocity in HT:
      Real V1_plasma_HT[3];
      for (int i=0; i<3; ++i) V1_plasma_HT[i] = V1_plasma[i] - V_US_SNIF_2_HT[i];

      // Downstream B:
      Real B2[3];
      B2[0] = B1[0];
      B2[1] = B1[1];
      B2[2] = B1[2] * plasma.R_magn;
      Real B2_mag = vectorMagnitude<3>(B2);

      // Downstream plasma velocity in HT and its component 
      // parallel to downstream magnetic field:
      Real V2_plasma_HT[3];
      V2_plasma_HT[0] = V1_plasma_HT[0] / plasma.R_gas;
      V2_plasma_HT[1] = 0;
      V2_plasma_HT[2] = V1_plasma_HT[2] * plasma.R_magn/plasma.R_gas;      
      Real V2_plasma_HT_par = dotProduct<3>(V2_plasma_HT,B2)/B2_mag;

      // Transform velocity from upstream SNIF to HT frame:
      //Real V_DS_SNIF_2_HT[3];
      //getTransformLocalSNIFToLocalHT(V2_plasma,B2,shockNormal,V_DS_SNIF_2_HT);

      // Component of HT transform velocity parallel to upstream B:
      Real V_SNIF_2_HT_par = dotProduct<3>(V_US_SNIF_2_HT,B1) / B1_mag;

      // Particle parallel speed in HT frame:
      Real V2_par_HT = particle.state[shockaccelerator::V_PAR] - V_SNIF_2_HT_par;

      // Particle parallel speed in downstream plasma rest frame,
      // note that in this coordinate system electric field is still zero:
      Real V2_par_plasma = V2_par_HT - V2_plasma_HT_par;
      Real V_gyro2 = 2*particle.state[shockaccelerator::MU]*B2_mag/species.mass;
      Real V_mag = sqrt(V2_par_plasma*V2_par_plasma + V_gyro2);
      Real pitchOld = V2_par_plasma / V_mag;
      Real pitchLimit = -V2_plasma_HT_par/V_mag;

      // Exit if particle speed is too low to return to shock:
      if (V_mag <= fabs(V2_plasma_HT_par)) return -1;

      // Calculate new random pitch for the particle so that 
      // it can return to shock:
      Real pitchNew;
      if (plasma.B1_norm < 0) {
	 Real P_max = 2*fabs(1+pitchLimit)/((1-pitchLimit)*(1-pitchLimit));
	 Real probab,y;
	 do {
	    pitchNew = pitchLimit + (1-pitchLimit)*simClasses->random.uniform();
	    y = P_max * simClasses->random.uniform();
	    probab = 2*fabs(pitchNew-pitchLimit)/((1-pitchLimit)*(1-pitchLimit));
	 } while (y > probab);
      } else {
	 //pitchNew = -1 + (1+pitchLimit)*simClasses->random.uniform();
	 Real P_max = 2*fabs(1-pitchLimit)/((1-pitchLimit)*(1-pitchLimit));
	 Real probab,y;
	 do {
	    pitchNew = -1 + (1+pitchLimit)*simClasses->random.uniform();
	    y = P_max * simClasses->random.uniform();
	    probab = 2*fabs(pitchNew-pitchLimit)/((1+pitchLimit)*(1+pitchLimit));
	 } while (y > probab);
      }

      // New parallel speed and magnetic moment in HT frame:
      particle.state[shockaccelerator::MU] = 0.5*species.mass*V_mag*V_mag*(1-pitchNew*pitchNew)/B2_mag;
      Real V2_par_HT_new = V_mag*pitchNew + V2_plasma_HT_par;
       
      // Transform to downstream SNIF because solveShockEncounterAnalytic 
      // below expects values in SNIF frame:
      V2_plasma_HT[0] = 0;
      V2_plasma_HT_par = dotProduct<3>(V2_plasma_HT,B2)/B2_mag;
      particle.state[shockaccelerator::V_PAR] = V2_par_HT_new + V_SNIF_2_HT_par - V2_plasma_HT_par;

      // Solve shock encounter:
      solveShockEncounterAnalytic(plasma,particle);

      //#ifndef NDEBUG
         if (particle.state[shockaccelerator::XPOS] > 0) {
	    cerr << "ERROR: particle remains in downstream" << endl;
	    cerr << "pitchLimit " << pitchLimit << endl;
	    cerr << V2_plasma_HT_par << '\t' << V_mag << endl;
	    exit(1);
	 }
      //#endif

      // Returned value is the probability of return multiplied by flux 
      // weighting. Macroparticle statistical weight will be multiplied 
      // by this number.
      /*if (plasma.B1_norm < 0) {
	 return +0.5*(1-pitchLimit) * 2*pitchNew/(1-pitchLimit*pitchLimit);
      } else {
	 return -0.5*(1+pitchLimit) * 2*pitchNew/(1-pitchLimit*pitchLimit);
      }*/
      Real P_ret = (V_mag-V2_plasma_HT_par)/(V_mag+V2_plasma_HT_par);
      return P_ret*P_ret;
      //return +1.0;
   }

   bool ShockAccelerator::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      bool success = true;
      simClasses.logger << "(SEP SHOCK ACCELERATOR) Starting initialization" << endl;
      
      this->sim = &sim;
      this->simClasses = &simClasses;

      // Parse required config file parameters:
      cr.add(prefix+".phi_factor","Fraction of shock ram energy used to calculate cross-shock potential (Real)",(Real)0.08);
      if (cr.parse() == false) exit(1);
      cr.get(prefix+".phi_factor",phiFactor);

      if (phiFactor < 0.0) {
	 simClasses.logger << "\t ERROR: Negative value in phi_factor" << endl;
	 success = false;
      }

      // Write init status and exit:
      simClasses.logger << "\t Initialization done, status is ";
      if (success == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;

      return success;
   }
   
   void ShockAccelerator::propagateBorisBuneman(Real* RESTRICT state,const Real* RESTRICT E,const Real* RESTRICT B,Real eforce,Real dt) {
      // Calculate v-minus velocity
      Real v_minus[3];
      v_minus[0] = state[3] + eforce * E[0];
      v_minus[1] = state[4] + eforce * E[1];
      v_minus[2] = state[5] + eforce * E[2];
      
      // Calculate t and s vectors needed in rotation of velocity
      Real t[3];
      t[0] = eforce * B[0];
      t[1] = eforce * B[1];
      t[2] = eforce * B[2];
      const Real t_mag = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
      
      const Real fact2 = 2.0 / (1.0 + t_mag);
      Real s[3];
      s[0] = fact2 * t[0];
      s[1] = fact2 * t[1];
      s[2] = fact2 * t[2];
      
      // Calculate v' velocity
      Real v_dot[3];
      v_dot[0] = v_minus[0] + v_minus[1]*t[2] - v_minus[2]*t[1];
      v_dot[1] = v_minus[1] + v_minus[2]*t[0] - v_minus[0]*t[2];
      v_dot[2] = v_minus[2] + v_minus[0]*t[1] - v_minus[1]*t[0];
      
      // Calculate v-plus velocity
      Real v_plus[3];
      v_plus[0] = v_minus[0] + v_dot[1]*s[2] - v_dot[2]*s[1];
      v_plus[1] = v_minus[1] + v_dot[2]*s[0] - v_dot[0]*s[2];
      v_plus[2] = v_minus[2] + v_dot[0]*s[1] - v_dot[1]*s[0];
       
      // Write new velocity to particle
      state[3] = v_plus[0] + eforce * E[0];
      state[4] = v_plus[1] + eforce * E[1];
      state[5] = v_plus[2] + eforce * E[2];
      
      // Write new position to particle
      state[0] += state[3] * dt;
      state[1] += state[4] * dt;
      state[2] += state[5] * dt;
   }

   bool ShockAccelerator::runTests() {
      bool success = true;

      #if PROFILE_LEVEL > 0
         static int counterTotal = -1;
         profile::start("ShockAccelerator Test Suite",counterTotal);
      #endif

      //if (testFermiI() == false) success = false;
      //if (testShockReturn() == false) success = false;
      if (testWebb() == false) success = false;
      if (testWriteSampleFields() == false) success = false;

      #if PROFILE_LEVEL > 0
         profile::stop();
      #endif
      return success;
   }
   
   void ShockAccelerator::setSpecies(const sep::Species& species) {
      this->species = species;
   }
   
   int ShockAccelerator::solveShockEncounterAnalytic(const PlasmaParameters& plasma,ParticleParameters& particle) {
      Real B1_mag,B2_mag;
      Real R_magn=1.0, R_gas=1.0;
      Real V1_SNIF_2_HT,V2_SNIF_2_HT;
      if (particle.state[shockaccelerator::XPOS] > 0.0) {
	 // 1: values in downstream region
	 // 2: values in upstream region
	 R_magn = plasma.R_magn;
	 R_gas  = plasma.R_gas;
	 B1_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.R_magn*plasma.R_magn*plasma.B1_tang*plasma.B1_tang);
	 B2_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.B1_tang*plasma.B1_tang);

	 const Real tan_psi1 = plasma.R_magn*plasma.B1_tang/plasma.B1_norm;
	 const Real tan_psi2 = plasma.B1_tang/plasma.B1_norm;

	 const Real Bz1 = plasma.B1_tang*plasma.R_magn;
	 const Real Bz2 = plasma.B1_tang;

	 V1_SNIF_2_HT = -plasma.V1_plasma_norm/plasma.R_gas * tan_psi1 * Bz1/B1_mag;
	 V2_SNIF_2_HT = -plasma.V1_plasma_norm              * tan_psi2 * Bz2/B2_mag;
      } else {
	 // 1: values in upstream region
	 // 2: values in downstream region
	 B1_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.B1_tang*plasma.B1_tang);
	 B2_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.R_magn*plasma.R_magn*plasma.B1_tang*plasma.B1_tang);

	 const Real tan_psi1 = plasma.B1_tang/plasma.B1_norm;
	 const Real tan_psi2 = plasma.R_magn*tan_psi1;
	 const Real Bz1 = plasma.B1_tang;
	 const Real Bz2 = plasma.B1_tang*plasma.R_magn;

	 V1_SNIF_2_HT = -plasma.V1_plasma_norm              * tan_psi1 * Bz1/B1_mag;
	 V2_SNIF_2_HT = -plasma.V1_plasma_norm/plasma.R_gas * tan_psi2 * Bz2/B2_mag;
      }

      // Get parallel component of transformation from shock normal
      // incidence frame to de Hoffmann-Teller frame:
      // NOTE: You can comment out the following code...
      Real V1_plasma[3];
      V1_plasma[0] = plasma.V1_plasma_norm / R_gas;
      V1_plasma[1] = 0.0;
      V1_plasma[2] = plasma.V1_plasma_tang;

      Real B1[3];
      B1[0] = plasma.B1_norm;
      B1[1] = 0.0;
      B1[2] = plasma.B1_tang * R_magn;

      Real shockNormal[3];
      shockNormal[0] = -1.0;
      shockNormal[1] = 0.0;
      shockNormal[2] = 0.0;

      Real V_SNIF_2_HT[3];
      getTransformLocalSNIFToLocalHT(V1_plasma,B1,shockNormal,V_SNIF_2_HT);
      const Real V_SNIF_2_HT_par = dotProduct<3>(V_SNIF_2_HT,B1) / B1_mag;

      // NOTE (CONT'D) ...and replace with this function call.
      // Computational speed is practically the change, however
      // (difference was ~1.5% on Arto's laptop):
      //const Real V_SNIF_2_HT_par = getTransformLocalSNIFToLocalHT(plasma.B1_norm,plasma.B1_tang*R_magn,plasma.V1_plasma_norm/R_gas);

      // Calculate guiding center parallel speed and pitch in HT frame:
      //const Real V_GC_HT_par = particle.state[shockaccelerator::V_PAR] - V_SNIF_2_HT_par;
      const Real V_GC_HT_par = particle.state[shockaccelerator::V_PAR] - V1_SNIF_2_HT;
      const Real V_gyro2 = 2.0*particle.state[shockaccelerator::MU]*B1_mag / species.mass;
      const Real V_GC_HT_mag = sqrt(V_GC_HT_par*V_GC_HT_par + V_gyro2);
      const Real mu_HT = V_GC_HT_par / V_GC_HT_mag;

      // Check if particle encounters shock:
      if (particle.state[shockaccelerator::XPOS] < 0.0) {
	 //if (plasma.B1_norm < 0.0 && mu_HT > 0.0) return 0;
	 //if (plasma.B1_norm > 0.0 && mu_HT < 0.0) return 0;
	 if (plasma.B1_norm*mu_HT < 0.0) {
	    //cerr << "\t SDA no enc " << plasma.B1_norm << '\t' << mu_HT << endl;
	    //exit(1);
	    return -1;
	 }
      } else {
	 //if (plasma.B1_norm < 0.0 && mu_HT < 0.0) return 0;
	 //if (plasma.B1_norm > 0.0 && mu_HT > 0.0) return 0;
	 if (plasma.B1_norm*mu_HT > 0.0) {
	    //cerr << "\t SDA no enc " << mu_HT << " V " << V1_plasma[0] << '\t' << V1_plasma[2] << '\t' << endl;
	    //exit(1);
	    return -1;
	 }
      }

      // Calculate loss cone cosine:
      Real mu_lc = sqrt(max(0.0,1.0 - B1_mag/B2_mag));
      if (plasma.B1_norm < 0.0) mu_lc *= -1.0;

      if (particle.state[shockaccelerator::XPOS] < 0.0) { // Particle initially in upstream:
	 if (plasma.B1_norm <= 0.0) {
	    if (mu_HT <= mu_lc) {
	       // Particle transmitted:
	       const Real mu_HT_new = -sqrt(max(0.0,1.0 - (1.0-mu_HT*mu_HT)*B2_mag/B1_mag));
	       //particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V_SNIF_2_HT_par;
	       particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V2_SNIF_2_HT;
	       particle.state[shockaccelerator::XPOS] = 0.501*plasma.L_shock;
	    } else {
	       // Particle reflected:
	       //particle.state[shockaccelerator::V_PAR] = -V_GC_HT_par + V_SNIF_2_HT_par;
	       particle.state[shockaccelerator::V_PAR] = -V_GC_HT_par + V1_SNIF_2_HT;
	       //cerr << "refl " << mu_HT << '\t' << mu_lc << endl;
	    }
	 } else {
	    if (mu_HT >= mu_lc) {
	       // Particle transmitted:
	       const Real mu_HT_new = +sqrt(max(0.0,1.0 - (1.0-mu_HT*mu_HT)*B2_mag/B1_mag));

	       //particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V_SNIF_2_HT_par;
	       particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V2_SNIF_2_HT;
	       particle.state[shockaccelerator::XPOS] = 0.501*plasma.L_shock;
	    } else {
	       // Particle reflected:
	       //particle.state[shockaccelerator::V_PAR] = -V_GC_HT_par + V_SNIF_2_HT_par;
	       particle.state[shockaccelerator::V_PAR] = -V_GC_HT_par + V1_SNIF_2_HT;
	       //cerr << "refl " << mu_HT << '\t' << mu_lc << endl;
	    }
	 }
      } else {
	 if (plasma.B1_norm < 0.0) { // Particle initially in downstream:
	    const Real mu_HT_new = +sqrt(max(0.0,1.0 - (1.0-mu_HT*mu_HT)*B2_mag/B1_mag));
	    //particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V_SNIF_2_HT_par;
	    particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V2_SNIF_2_HT;
	    particle.state[shockaccelerator::XPOS] = -0.501*plasma.L_shock;
	 } else {
	    const Real mu_HT_new = -sqrt(max(0.0,1.0 - (1.0-mu_HT*mu_HT)*B2_mag/B1_mag));
	    //particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V_SNIF_2_HT_par;
	    particle.state[shockaccelerator::V_PAR] = V_GC_HT_mag*mu_HT_new + V2_SNIF_2_HT;
	    particle.state[shockaccelerator::XPOS] = -0.501*plasma.L_shock;
	 }
      }

      return 1;
   }
   
   int ShockAccelerator::solveShockEncounterGC(const PlasmaParameters& plasma,ParticleParameters& particle) {
      // Calculate non-coplanar electric field (=constant):
      Real E_x = 0.0;
      const Real E_ncp = plasma.B1_tang*plasma.V1_plasma_norm - plasma.B1_norm*plasma.V1_plasma_tang;
      
      // Row 1, column 3 component of grad B tensor:
      const Real B_GRAD_13 = plasma.B1_tang*(plasma.R_magn - 1)/plasma.L_shock;
      
      const Real B1_mag = sqrt(plasma.B1_norm*plasma.B1_norm + plasma.B1_tang*plasma.B1_tang);
      const Real cos_psi1 = fabs(plasma.B1_norm / B1_mag); // fabs needed to make this work with B1_norm < 0
      const Real V_mag = fabs(particle.state[shockaccelerator::V_PAR])*cos_psi1 + fabs(E_ncp)/B1_mag;
      const Real dt_new = plasma.L_shock / V_mag;
      
      const Real dt = 0.0001*dt_new;
      int iterations = 0;
      while (fabs(particle.state[shockaccelerator::XPOS]) < 0.51*plasma.L_shock) {
	 Real R_magn;
	 Real B_grad_13;
	 if (particle.state[shockaccelerator::XPOS] < -0.5*plasma.L_shock) {
	    R_magn = 1.0;
	    B_grad_13 = 0.0;
	    E_x = 0.0;
	 } else if (particle.state[shockaccelerator::XPOS] <= 0.5*plasma.L_shock) {
	    const Real R_gas = 1 + (plasma.R_gas-1)*(particle.state[shockaccelerator::XPOS]+plasma.L_shock*0.5)/plasma.L_shock;
	    const Real V_norm = plasma.V1_plasma_norm / R_gas;
	    E_x = -phiFactor*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*V_norm*V_norm/R_gas*(plasma.R_gas-1)/plasma.L_shock;
	    
	    R_magn = 1 + (plasma.R_magn-1)*(particle.state[shockaccelerator::XPOS]+plasma.L_shock*0.5)/plasma.L_shock;
	    B_grad_13 = B_GRAD_13;
	 } else {
	    R_magn = plasma.R_magn;
	    B_grad_13 = 0.0;
	    E_x = 0.0;
	 }
      
	 // Calculate E,B fields:
	 Real B_x = plasma.B1_norm;
	 Real B_z = plasma.B1_tang * R_magn;
	 Real B_mag2 = B_x*B_x + B_z*B_z;
	 Real B_mag = sqrt(B_mag2);
	 
	 // Calculate x,z components of electric drift:
	 Real V_e_x = +E_ncp / B_mag2 * B_z;
	 //Real V_e_y = -E_x / B_mag2 * B_z;
	 Real V_e_z = -E_ncp / B_mag2 * B_x;

	 // Calculate parallel velocity vector:
	 Real V_par_x = particle.state[shockaccelerator::V_PAR] * B_x / B_mag;
	 Real V_par_z = particle.state[shockaccelerator::V_PAR] * B_z / B_mag;

	 // Calculate parallel acceleration:
	 Real parallelAcceleration = species.q_per_m*E_x*B_x/B_mag - particle.state[shockaccelerator::MU]/species.mass * B_x*B_z/B_mag2 * B_grad_13;
	 parallelAcceleration += V_e_z/B_mag*(V_par_x + V_e_x)*B_grad_13;
	 
	 // Propagate:
	 particle.state[shockaccelerator::XPOS] += (V_par_x + V_e_x)*dt;
	 particle.state[shockaccelerator::ZPOS] += (V_par_z + V_e_z)*dt;
	 particle.state[shockaccelerator::V_PAR] += parallelAcceleration*dt;
	 ++iterations;
      }
      return iterations;
   }
   
   int ShockAccelerator::solveShockEncounterLorentz(const PlasmaParameters& plasma,ParticleParameters& particle) {
      // Create a Lorentz particle that matches the given parameters
      // and assume random gyro phase angle:
      Real B[3];
      B[0] = plasma.B1_norm;
      B[1] = 0.0;
      B[2] = plasma.B1_tang;
      if (particle.state[shockaccelerator::XPOS] > 0.0) B[2] *= plasma.R_magn;
      Real B_mag = vectorMagnitude<3>(B);
      
      // Create an orthonormal basis where w3 points to the 
      // direction of upstream B, w1 and w2 complete the set.
      Real w1[3];
      Real w2[3];
      Real w3[3];
      w3[0] = B[0];
      w3[1] = 0.0;
      w3[2] = B[2];
      unitVector<3>(w3);
      
      for (int i=0; i<3; ++i) w1[i] = 0.0;
      w1[1] = 1.0;
      const Real w1_dot_w3 = dotProduct<3>(w1,w3);
      for (int i=0; i<3; ++i) w1[i] -= w1_dot_w3*w3[i];
      unitVector<3>(w1);
      
      crossProduct(w3,w1,w2);
      unitVector<3>(w2);
      
      // Calculate gyro speed:
      const Real v_gyro = sqrt(2*particle.state[shockaccelerator::MU]*B_mag/species.mass);
      const Real r_gyro = fabs(v_gyro / B_mag / species.q_per_m);

      // Generate random gyro phase:
      const Real phi = 2.0*M_PI*simClasses->random.uniform();
      const Real cos_phi = cos(phi);
      const Real sin_phi = sin(phi);
      //const Real sin_phi = sqrt(1.0 - cos_phi*cos_phi);

      Real state[6];
      state[0] = particle.state[shockaccelerator::XPOS];
      state[1] = particle.state[shockaccelerator::YPOS];
      state[2] = particle.state[shockaccelerator::ZPOS];
      state[3] = particle.state[shockaccelerator::V_PAR]*w3[0] + v_gyro*cos_phi*w1[0] + v_gyro*sin_phi*w2[0];
      state[4] = particle.state[shockaccelerator::V_PAR]*w3[1] + v_gyro*cos_phi*w1[1] + v_gyro*sin_phi*w2[1];
      state[5] = particle.state[shockaccelerator::V_PAR]*w3[2] + v_gyro*cos_phi*w1[2] + v_gyro*sin_phi*w2[2];

      // Offset x-position:
      const Real psi1 = acos(B[0] / B_mag);
      const Real r_gyro_x = r_gyro * cos(0.5*M_PI - psi1);
      if (particle.state[shockaccelerator::XPOS] < 0.0) state[0] -= 2*r_gyro_x;
      else state[0] += 2.0*r_gyro_x;
      
      Real E[3];
      for (int i=0; i<3; ++i) E[i] = 0.0;
      E[1] = plasma.B1_tang*plasma.V1_plasma_norm - plasma.B1_norm*plasma.V1_plasma_tang;
      
      // Add E cross B drift velocity:
      state[3] += E[1]*B[2] / (B_mag*B_mag);
      state[5] -= E[1]*B[0] / (B_mag*B_mag);
      
      const Real omega = species.q_per_m * B_mag;
      const Real dt = 1 / (200 * fabs(omega));
      const Real eforce = 0.5*species.q_per_m*dt;

      const Real x_min = -0.5*plasma.L_shock - 4*r_gyro;
      const Real x_max = +0.5*plasma.L_shock + 4*r_gyro;

      borisBunemanBehind(state,E,B,eforce);
      int iterations = 0;
      while (x_min <= state[0] && state[0] <= x_max) {
	 // Calculate tangential magnetic field:
	 Real R_magn = 1.0;
	 if (state[0] < -0.5*plasma.L_shock) {
	    R_magn = +1.0;
	    E[0] = 0.0;
	 } else if (state[0] > 0.5*plasma.L_shock) {
	    R_magn = plasma.R_magn;
	    E[0] = 0.0;
	 } else {
	    R_magn = 1 + (plasma.R_magn-1)*(state[0]+0.5*plasma.L_shock)/plasma.L_shock;
	    const Real R_gas = 1.0 + (plasma.R_gas-1)*(state[0]+0.5*plasma.L_shock)/plasma.L_shock;
	    const Real V_norm = plasma.V1_plasma_norm / R_gas;
	    E[0] = -phiFactor*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*V_norm*V_norm/R_gas*(plasma.R_gas-1)/plasma.L_shock;
	 }
	 B[2] = R_magn * plasma.B1_tang;

	 propagateBorisBuneman(state,E,B,eforce,dt);
	 ++iterations;
      }
      borisBunemanAhead(state,E,B,eforce);

      // Recalculate V_parallel and magnetic moment from Lorentz particle:
      
      // First remove electric drift velocity:
      B[2] = plasma.B1_tang;
      if (state[0] > 0.0) B[2] *= plasma.R_magn;
      B_mag = vectorMagnitude<3>(B);
      state[3] -= E[1]*B[2] / (B_mag*B_mag);
      state[5] += E[1]*B[0] / (B_mag*B_mag);
      
      // Calculate parallel speed:
      const Real V_parallel = dotProduct<3>(&(state[3]),B) / B_mag;
      particle.state[shockaccelerator::V_PAR] = V_parallel;
      
      // Calculate magnetic moment:
      const Real U_gyro = 0.5*species.mass*(vectorMagnitude2<3>(&(state[3])) 
			- V_parallel*V_parallel);
      particle.state[shockaccelerator::MU] = U_gyro / B_mag;
      
      // Calculate GC position:
      Real r_gyro_vec[3];
      crossProduct(&(state[3]),B,r_gyro_vec);
      for (int i=0; i<3; ++i) r_gyro_vec[i] /= (species.q_per_m*B_mag*B_mag);
      particle.state[shockaccelerator::XPOS] = state[0] - r_gyro_vec[0];
      particle.state[shockaccelerator::YPOS] = state[1] - r_gyro_vec[1];
      particle.state[shockaccelerator::ZPOS] = state[2] - r_gyro_vec[2];
      return iterations;
   }

   bool ShockAccelerator::testShockReturn() {
      bool success = true;

      int N_psi = 150;
      Real cos_psi1_0 = 0.0;
      Real cos_psi1_1 = 1.0;
      Real d_cos_psi = (cos_psi1_1-cos_psi1_0)/(N_psi-1);

      Real U_ini = 1000 * constants::CHARGE_ELEMENTARY * 100.0;
      
      fstream out("test_shock_return.txt",fstream::out);
      if (out.good() == false) return false;
      
      for (int i=0; i<N_psi; ++i) {
	 Real cos_psi1 = cos_psi1_0 + i*d_cos_psi;
         //Real cos_psi1 = 1.0;
	 Real P_return,dU;
	 if (testShockReturn(U_ini,cos_psi1,P_return,dU) == false) {
	    success = false;
	    cerr << "discard " << cos_psi1 << endl;
	    exit(1);
	 } else {
	    out << cos_psi1 << '\t';
	    out << P_return << '\t';
	    out << dU << '\t';
	    out << endl;
	 }
      }

      out.close();
      return success;
   }
   
   bool ShockAccelerator::testFermiI() {
      Real cos_psi1 = 1.0;
      Real U_ini = 1000 * constants::CHARGE_ELEMENTARY * 100.0;
      if (testFermiI(U_ini,cos_psi1) == false) return false;
      return true;
   }

   bool ShockAccelerator::testFermiI(Real U_ini,Real cos_psi1) {
      Real B_mag = 1e-6;
      Real V_shock = -1000 * 1500.0;
      Real R_gas = 3.6;
      
      PlasmaParameters plasma;
      if (testSetupUpstream(R_gas,B_mag,V_shock,cos_psi1,plasma) == false) return false;
      
      Species species;
      species.charge  = constants::CHARGE_ELEMENTARY;
      species.mass    = constants::MASS_PROTON;
      species.q_per_m = species.charge / species.mass;
      setSpecies(species);
      
      ParticleParameters p;
      p.state[shockaccelerator::XPOS] = +0.5*plasma.L_shock;
      p.state[shockaccelerator::YPOS]  = 0.0;
      p.state[shockaccelerator::ZPOS]  = 0.0;

      fstream out("test_Fermi_I.txt",fstream::out);
      if (out.good() == false) return false;
      
      Real U1 = fabs(V_shock);
      Real U2 = fabs(V_shock)/R_gas;
      Real avgGainConstant = 4.0/3.0*(U1-U2);
      
      out << "# Upstream speed     : " << fabs(V_shock) << endl;
      out << "# Downstream speed   : " << fabs(V_shock)/R_gas << endl;
      out << "# Avg. velocity gain : " << avgGainConstant << endl;
      
      Real probability = 1.0;
      Real speed = sqrt(2*U_ini/species.mass);
      Real speed0 = speed;
      int N = 100;
      int N_averages = 500;
      Real totalAvgGain = 0.0;
      Real totalWeightSum = 0.0;

      for (int i=0; i<N; ++i) {
	 Real speed_ini = speed;

	 Real avgEnergyGain = 0.0;
	 Real weightSum = 0.0;
	 Real avgGain = 0.0;
	 Real avgSpeed = 0.0;
	 for (int j=0; j<N_averages; ++j) {
	    // Give particle a random pitch so that it will encounter the shock:
	    Real pitch = simClasses->random.uniform();
	    p.state[shockaccelerator::XPOS]  = -0.5*plasma.L_shock;
	    p.state[shockaccelerator::V_PAR] = speed_ini * pitch - V_shock;
	    p.state[shockaccelerator::MU]    = 0.5*species.mass*speed_ini*speed_ini*(1-pitch*pitch)/B_mag;

	    // Solve shock encounter and get return weight:
	    Real dU_analytic;
	    Real weight = p.state[shockaccelerator::V_PAR];
	    weight *= getReturnedParticle(plasma,p,0.0,dU_analytic);

	    Real V_par = p.state[shockaccelerator::V_PAR] + V_shock;
	    speed = sqrt(V_par*V_par + 2*p.state[shockaccelerator::MU]*B_mag/species.mass);
	    Real U = 0.5*species.mass*V_par*V_par + p.state[shockaccelerator::MU]*B_mag;

	    avgGain += weight * (speed-speed_ini);
	    avgSpeed += weight * speed;
	    weightSum += weight;
	    avgEnergyGain += weight * (U-U_ini);
	 }
	 totalAvgGain += avgGain;
	 avgEnergyGain /= weightSum;
	 avgGain /= weightSum;
	 avgSpeed /= weightSum;
	 totalWeightSum += weightSum;

	 speed = avgSpeed;
	 out << i << '\t' << avgGain << '\t' << speed/speed0 << endl;
      }
      totalAvgGain /= totalWeightSum;
      out << "# Total avg gain " << totalAvgGain << endl;
      out << "# Weight sum     " << totalWeightSum << endl;

      out.close();
      return true;
   }
   
   bool ShockAccelerator::testSetupUpstream(Real R_gas,Real B_mag,Real V_shock,Real cos_psi1,PlasmaParameters& plasma) {
      Real V1_plasma_mag = 0.0;
      Real V1_alfven = 1000 * 700.0;

      Real alfvenMach2 = fabs(V1_plasma_mag-V_shock) / V1_alfven;
      Real R_magn = (alfvenMach2 - cos_psi1*cos_psi1)/(alfvenMach2 - cos_psi1*cos_psi1*R_gas)*R_gas;

      plasma.B1_norm        = B_mag*cos_psi1;
      plasma.B1_tang        = B_mag*sin(acos(cos_psi1));
      plasma.V1_plasma_norm = -V_shock;
      plasma.V1_plasma_tang = 0.0;
      plasma.R_gas          = R_gas;
      plasma.R_magn         = R_magn;
      plasma.L_shock        = 1000.0;
    
      return true;
   }
   
   bool ShockAccelerator::testShockReturn(Real U_ini,Real cos_psi1,Real& P_return_avg,Real& dU_avg) {
      Real B_mag = 1e-6;
      //Real V1_plasma_mag = 0.0;
      Real V_shock = -1000 * 1500.0;
      //Real V1_alfven = 1000 * 700.0;

      //Real R_gas = 3.6;
      //Real alfvenMach2 = fabs(V1_plasma_mag-V_shock) / V1_alfven;
      //Real R_magn = (alfvenMach2 - cos_psi1*cos_psi1)/(alfvenMach2 - cos_psi1*cos_psi1*R_gas)*R_gas;

      // Create plasma parameters for shock accelerator:
      PlasmaParameters plasma;
      //plasma.B1_norm        = B_mag*cos_psi1;
      //plasma.B1_tang        = B_mag*sin(acos(cos_psi1));
      //plasma.V1_plasma_norm = -V_shock;
      //plasma.V1_plasma_tang = 0.0;
      //plasma.R_gas          = R_gas;
      //plasma.R_magn         = R_magn;
      //plasma.L_shock        = 1000.0;
      if (testSetupUpstream(3.6,B_mag,V_shock,cos_psi1,plasma) == false) return false;
      
      Species species;
      species.charge  = constants::CHARGE_ELEMENTARY;
      species.mass    = constants::MASS_PROTON;
      species.q_per_m = species.charge / species.mass;
      setSpecies(species);
      /*
      // Get downstream parameters:
      Real V2_wave_par;
      Real B2[3];
      Real B2_mag;
      Real V2_wave[3];
      Real tan_psi2;
      getDownstreamValues(plasma,V1_alfven,V2_wave_par,B2,B2_mag,V2_wave,tan_psi2);
      */
      ParticleParameters p;
      p.state[shockaccelerator::XPOS] = +0.5*plasma.L_shock;
      p.state[shockaccelerator::YPOS]  = 0.0;
      p.state[shockaccelerator::ZPOS]  = 0.0;

      Real speed = sqrt(2*U_ini/species.mass);
      int N = 10000;
      P_return_avg = 0.0;
      dU_avg = 0.0;
      int N_encounters = 0;
      for (int i=0; i<N; ++i) {
	 // Calculate random pitch:
	 Real pitch = 1.0 - 2.0*corsair::getObjectWrapper().simClasses.random.uniform();

	 // Parallel speed in upstream rest frame is speed*pitch,
	 // transform to shock normal incidence frame:
	 p.state[shockaccelerator::XPOS]  = -0.5*plasma.L_shock;
	 p.state[shockaccelerator::V_PAR] = speed * pitch;// - V_shock*cos_psi1;
	 p.state[shockaccelerator::MU]    = 0.5*species.mass*speed*speed*(1-pitch*pitch)/B_mag;

	 //if (p.state[shockaccelerator::V_PAR]*cos_psi1 < 0) continue;

	 // Solve shock encounter:
	 Real tmp = getReturnedParticle(plasma,p,0.0,cos_psi1);

	 // If particle returned, collect statistics:
	 if (tmp >= 0.0) {
	    if (tmp > 1.0) {
	       cerr << "ERROR: Too large return weight " << tmp << endl;
	       exit(1);
	    }

	    P_return_avg += tmp;
	    Real V_par_2 = p.state[shockaccelerator::V_PAR] ;//+ V_shock*cos_psi1;
	    Real U_fin = 0.5*species.mass*V_par_2*V_par_2 + p.state[shockaccelerator::MU]*B_mag;

	    dU_avg += tmp * (U_fin-U_ini) / U_ini;
	    /*
	    if (cos_psi1 >= 0.77 && cos_psi1 <= 0.778) {
	       cerr << acos(pitch)*180/M_PI << '\t' << tmp << '\t' << tmp * (U_fin-U_ini) / U_ini << endl;
	    }
	    */
	 }
	 ++N_encounters;
      }

      P_return_avg /= N_encounters;
      dU_avg /= N_encounters;
      return true;
   }

   bool ShockAccelerator::testWebb() {
      bool success = true;
      
      #if PROFILE_LEVEL > 0
         static int counterTotal = -1;
         static int counterAnalytic = -1;
         static int counterGC = -1;
         static int counterLorentz = -1;
         static int counterElectrons = -1;
         profile::start("Webb test total",counterTotal);
      #endif
      
      #if PROFILE_LEVEL > 0
         profile::start("GC iterator",counterGC);
      #endif
      /*if (testWebb(400,1e6,60.0,"test_webb_60_00_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(400,1e6,80.0,"test_webb_80_00_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(400,1e6,88.0,"test_webb_88_00_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(400,1e6,89.0,"test_webb_89_00_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(400,1e6,89.92,"test_webb_89_92_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(400,1e6,89.93,"test_webb_89_93_gc.txt",true,100,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      if (testWebb(1500,1e6,39.2,"test_webb_gc.txt",true,300,&ShockAccelerator::solveShockEncounterGC) == false) success = false;
      */
      #if PROFILE_LEVEL > 0
         profile::stop(); // counterGC
      #endif
      
      #if PROFILE_LEVEL > 0
         profile::start("Lorentz iterator",counterLorentz);
      #endif
      /*if (testWebb(400,1e6,60.0,"test_webb_60_00_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(400,1e6,80.0,"test_webb_80_00_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(400,1e6,88.0,"test_webb_88_00_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(400,1e6,89.0,"test_webb_89_00_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(400,1e6,89.92,"test_webb_89_92_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(400,1e6,89.93,"test_webb_89_93_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      if (testWebb(1500,1e5,39.2,"test_webb_lorentz.txt",true,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;
      */
      #if PROFILE_LEVEL > 0
         profile::stop(); // counterLorentz
      #endif

      #if PROFILE_LEVEL > 0
         profile::start("Analytic Solver",counterAnalytic);
      #endif
      /*if (testWebb(400,1e6,60.0,"test_webb_60_00_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(400,1e6,80.0,"test_webb_80_00_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(400,1e6,88.0,"test_webb_88_00_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(400,1e6,89.0,"test_webb_89_00_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(400,1e6,89.92,"test_webb_89_92_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(400,1e6,89.93,"test_webb_89_93_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      if (testWebb(1500,1e5,39.2,"test_webb_analytic.txt",true,300,&ShockAccelerator::solveShockEncounterAnalytic) == false) success = false;
      */
      #if PROFILE_LEVEL > 0
         profile::stop(); // counterAnalytic
         profile::start("Electrons",counterElectrons);
      #endif
      
      if (testWebb(1200.0,1e3,45.0,"test_webb_00_electrons_gc.txt",false,300,&ShockAccelerator::solveShockEncounterGC) == false) success = false;      
      if (testWebb(1200.0,1e3,45.0,"test_webb_00_electrons_lorentz.txt",false,300,&ShockAccelerator::solveShockEncounterLorentz) == false) success = false;

      #if PROFILE_LEVEL > 0
         profile::stop();
         profile::stop(); // counterTotal
      #endif
      return success;
   }
   
   bool ShockAccelerator::testWebb(Real V_plasma,Real energy,Real psi,const std::string& fname,bool ions,uint32_t N_particles,
				   int (ShockAccelerator::*solver)(const PlasmaParameters& plasma,ParticleParameters& particle)) {
      PlasmaParameters params;
      
      // In Webb test gas and magnetic compression ratios are equal.

      const Real B1_mag = 1.0e-8;
      const Real psi1 = psi * M_PI/180.0;
      params.B1_norm = B1_mag * cos(psi1);
      params.B1_tang = B1_mag * sin(psi1);
      params.V1_plasma_norm = V_plasma * 1000.0;
      params.V1_plasma_tang = 0.0;
      params.R_gas = 3.0;
      if (psi > 0.0) {
	 params.R_gas  = min(4.0,sqrt((8 + 9*tan(psi1)*tan(psi1))/(tan(psi1)*tan(psi1))));
      }
      params.R_magn = params.R_gas;
      params.L_shock = 10*100.0;

      if (ions == true) {
	 species.mass = constants::MASS_PROTON;
	 species.charge = constants::CHARGE_ELEMENTARY;
	 species.q_per_m = species.charge / species.mass;
      } else {
	 species.mass = constants::MASS_ELECTRON;
	 species.charge = -constants::CHARGE_ELEMENTARY;
	 species.q_per_m = species.charge / species.mass;
      }

      const Real U_tot = energy * constants::CHARGE_ELEMENTARY;
      const Real E_ncp = params.B1_tang*params.V1_plasma_norm - params.B1_norm*params.V1_plasma_tang;

      const Real mu_min = +1.0;
      const Real mu_max = -1.0;
      const Real d_mu = (mu_max-mu_min)/(N_particles-1);
      
      fstream out(fname.c_str(), fstream::out);
      if (out.good() == false) return false;
      
      const Real X_ini = 0.505*params.L_shock;
      
      const Real B1 = B1_mag;
      const Real U_e1 = 0.5*species.mass*E_ncp*E_ncp/B1/B1;
      const Real U_driftFree1 = U_tot - U_e1;
      ParticleParameters upstream;
      upstream.state[shockaccelerator::XPOS] = -X_ini;
      upstream.state[shockaccelerator::YPOS] = 0.0;
      upstream.state[shockaccelerator::ZPOS] = 0.0;
      upstream.state[shockaccelerator::V_PAR] = sqrt(2*U_driftFree1/species.mass);
      upstream.state[shockaccelerator::MU]    = U_driftFree1/B1;
      if (U_driftFree1 < 0.0) {
	 cerr << "(SEP SHOCK ACCELERATOR) ERROR: Negative energy" << endl;
	 return false;
      }

      const Real B2 = sqrt(params.B1_norm*params.B1_norm + params.R_magn*params.R_magn*params.B1_tang*params.B1_tang);
      const Real U_e2 = 0.5*species.mass*E_ncp*E_ncp/B2/B2;
      const Real U_driftFree2 = U_tot - U_e2;
      ParticleParameters downstream;
      downstream.state[shockaccelerator::XPOS] = +X_ini;
      downstream.state[shockaccelerator::YPOS] = 0.0;
      downstream.state[shockaccelerator::ZPOS] = 0.0;
      downstream.state[shockaccelerator::V_PAR] = sqrt(2*U_driftFree2/species.mass);
      downstream.state[shockaccelerator::MU]    = U_driftFree2/B2;
      if (U_driftFree2 < 0.0) {
	 cerr << "(SEP SHOCK ACCELERATOR) ERROR: Negative energy" << endl;
	 return false;
      }

      out << "#Pitch(ini) Pitch(fin) dU Pitch(fin) dU" << endl;
      for (size_t i=0; i<N_particles; ++i) {
	 const Real mu = mu_min + i*d_mu;
	 
	 // Solve upstream:
	 ParticleParameters particle = upstream;
	 particle.state[shockaccelerator::V_PAR] *= mu;
	 particle.state[shockaccelerator::MU]    *= (1-mu*mu);
	 
	 //solveShockEncounterLorentz(params,particle);
	 (this->*solver)(params,particle);
	 Real dU_up,pitch_up;
	 calculateWebbResult(params,particle,U_tot,dU_up,pitch_up);

	 // Solve downstream:
	 particle = downstream;
	 particle.state[shockaccelerator::V_PAR] *= mu;
	 particle.state[shockaccelerator::MU]    *= (1-mu*mu);

	 //solveShockEncounterLorentz(params,particle);
	 (this->*solver)(params,particle);
	 Real dU_down,pitch_down;
	 calculateWebbResult(params,particle,U_tot,dU_down,pitch_down);
	 
	 // Write out result:
	 out << acos(mu)*180.0/M_PI << '\t' << pitch_up << '\t' << dU_up << '\t' << pitch_down << '\t' << dU_down << endl;
      }
      out.close();
      
      return true;
   }
   
   bool ShockAccelerator::testWriteSampleFields() {
      const Real psi1 = 70.0 * M_PI/180.0;
      const Real Mach2 = 3.0;
      
      PlasmaParameters params;
      params.B1_norm = 1.0e-9;
      params.B1_tang = params.B1_norm * tan(psi1);
      params.V1_plasma_norm = 1000 * 1000.0;
      params.V1_plasma_tang = 0.0;
      params.R_gas  = 3.0;
      params.R_magn = (Mach2 - cos(psi1)*cos(psi1))/(Mach2 - params.R_gas*cos(psi1)*cos(psi1))*params.R_gas;

      const Real B = sqrt(params.B1_norm*params.B1_norm + params.B1_tang*params.B1_tang);
      params.B1_norm = 1.0e-9 * params.B1_norm / B;
      params.B1_tang = 1.0e-9 * params.B1_tang / B;
      
      return writeFields(params,"test_shock_fields.txt");
   }
   
   bool ShockAccelerator::writeFields(const PlasmaParameters& params,const std::string& fileName) {
      // Open output file:
      fstream out(fileName.c_str(), fstream::out);
      if (out.good() == false) return false;
      
      const Real x_min = -1.0;
      const Real x_max = +1.0;
      const int N_points = 100;
      const Real dx = (x_max-x_min)/(N_points-1);
      
      Real R_gas,R_magn;
      for (int i=0; i<N_points; ++i) {
	 const Real x = x_min + i*dx;
	 
	 Real V_tang;
	 if (x < -0.5) {
	    R_magn = 1.0;
	    R_gas = 1.0;
	 } else if (x > +0.5) {
	    R_magn = params.R_magn;
	    R_gas  = params.R_gas;
	 } else {
	    R_magn = 1 + (params.R_magn-1)*(x+0.5);
	    R_gas  = 1 + (params.R_gas-1)*(x+0.5);
	 }

	 const Real B_norm = params.B1_norm;
	 const Real B_tang = R_magn * params.B1_tang;
	 const Real V_norm = params.V1_plasma_norm / R_gas;
	 
	 const Real phi = 0.5*phiFactor*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*(params.V1_plasma_norm*params.V1_plasma_norm - V_norm*V_norm);
	 Real d_phi = 0.0;
	 if (x < -0.5) {
	    
	 } else if (x > +0.5) {
	    
	 } else {
	    d_phi = phiFactor*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY*V_norm*V_norm/R_gas*(params.R_gas-1);
	 }
	 V_tang = V_norm*B_tang/B_norm - params.V1_plasma_norm*params.B1_tang/params.B1_norm;
	 
	 const Real E_x = -d_phi;
	 const Real E_y = B_tang*V_norm - B_norm*V_tang;
	 
	 // Write fields:
	 out << x << '\t';
	 out << V_norm/1000 << '\t' << V_tang/1000 << '\t';
	 out << B_norm/1e-9 << '\t' << B_tang/1e-9 << '\t';
	 out << E_x << '\t' << E_y << '\t';
	 out << phi << '\t' << d_phi << '\t';
	 out << endl;
      }
      
      // Close output file:
      out.close();
      return true;
   }
   
} // namespace sep