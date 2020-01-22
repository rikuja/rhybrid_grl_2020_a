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

#ifndef SEP_PARTICLE_ACCELERATOR_H
#define SEP_PARTICLE_ACCELERATOR_H

#include <cstdlib>
#include <iostream>
#include <cmath>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <linear_algebra.h>

#include "sep_simcontrol.h"
#include "sep_fields_container.h"
#include "sep_particle_definition.h"

namespace sep {

   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class ParticleAccelerator {
    public:
      
      ParticleAccelerator();
      ~ParticleAccelerator();
      
      bool finalize();
      void get(pargrid::CellID blockLID,Real t,const SPECIES& species,Real* state,Real* acceleration);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      
    private:
      bool initialized;

      Simulation* sim;
      SimulationClasses* simClasses;
      SPECIES species;
      
      Real R;
      Real THETA;
      int32_t I;
      int32_t J;
      int32_t K;
   };
   
   template<class SPECIES,class PARTICLE> inline
   ParticleAccelerator<SPECIES,PARTICLE>::ParticleAccelerator() { 
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   ParticleAccelerator<SPECIES,PARTICLE>::~ParticleAccelerator() { 
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool ParticleAccelerator<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      return true;
   }
 
   /** Get guiding center acceleration at given position.
    * state[sep::particle::XCRD] = R,
    * state[sep::particle::YCRD] = PHI,
    * state[sep::particle::ZCRD] = Z.
    */
   template<class SPECIES,class PARTICLE> inline
   void ParticleAccelerator<SPECIES,PARTICLE>::get(pargrid::CellID blockLID,Real t,const SPECIES& species,Real* RESTRICT state,Real* RESTRICT acceleration) {
      // Get fields at particle position:
      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(blockLID,t,state,E,B,gradB);
      
      // Calculate GC electric drift velocity:
      Real V_electric[3];
      crossProduct(E,B,V_electric);
      const Real B_mag2 = vectorMagnitude2<3>(B);
      V_electric[0] /= B_mag2;
      V_electric[1] /= B_mag2;
      V_electric[2] /= B_mag2;
      
      // Calculate GC parallel velocity:
      const Real B_mag = sqrt(B_mag2);
      Real V_GC_par[3];
      V_GC_par[0] = state[sep::particle::V_PAR] * B[0] / B_mag;
      V_GC_par[1] = state[sep::particle::V_PAR] * B[1] / B_mag;
      V_GC_par[2] = state[sep::particle::V_PAR] * B[2] / B_mag;

      // Calculate parallel acceleration due to non-inertial terms:
      Real parallelAcceleration = 
	 V_electric[0] * (  (V_GC_par[0]+V_electric[0]) * gradB[matrixIndex(0,0)]
	   	          + (V_GC_par[1]+V_electric[1]) * gradB[matrixIndex(1,0)]
	   	          + (V_GC_par[2]+V_electric[2]) * gradB[matrixIndex(2,0)])
       + V_electric[1] * (  (V_GC_par[0]+V_electric[0]) * gradB[matrixIndex(0,1)]
		          + (V_GC_par[1]+V_electric[1]) * gradB[matrixIndex(1,1)]
		          + (V_GC_par[2]+V_electric[2]) * gradB[matrixIndex(2,1)])
       + V_electric[2] * (  (V_GC_par[0]+V_electric[0]) * gradB[matrixIndex(0,2)]
		          + (V_GC_par[1]+V_electric[1]) * gradB[matrixIndex(1,2)]
		          + (V_GC_par[2]+V_electric[2]) * gradB[matrixIndex(2,2)]);
      parallelAcceleration /= (-1.0*B_mag);
      
      // Calculate gradient of B:
      Real B_gradient[3];
      B_gradient[0] = (B[0]*gradB[matrixIndex(0,0)] + B[1]*gradB[matrixIndex(0,1)] + B[2]*gradB[matrixIndex(0,2)]) / B_mag;
      B_gradient[1] = (B[0]*gradB[matrixIndex(1,0)] + B[1]*gradB[matrixIndex(1,1)] + B[2]*gradB[matrixIndex(1,2)]) / B_mag;
      B_gradient[2] = (B[0]*gradB[matrixIndex(2,0)] + B[1]*gradB[matrixIndex(2,1)] + B[2]*gradB[matrixIndex(2,2)]) / B_mag;
      
      // Calculate magnitude of gradient force and add it to parallel acceleration:
      const Real gradientForce = dotProduct<3>(B,B_gradient)/B_mag*state[sep::particle::MU]/species.mass;
      parallelAcceleration -= gradientForce;
      
      acceleration[sep::particle::V_PAR] = parallelAcceleration;
      acceleration[sep::particle::MU]    = 0.0;

      // Calculate gradient drift velocity:
      Real V_gradient[3];
      crossProduct(B,B_gradient,V_gradient);
      V_gradient[0] *= (state[sep::particle::MU]/species.charge/B_mag2);
      V_gradient[1] *= (state[sep::particle::MU]/species.charge/B_mag2);
      V_gradient[2] *= (state[sep::particle::MU]/species.charge/B_mag2);
      
      // Calculate curvature drift velocity (B_gradient used as temporary array here):
      const Real curvatureFactor = state[sep::particle::V_PAR]*state[sep::particle::V_PAR]/species.q_per_m/B_mag2/B_mag2;
      B_gradient[0] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,0)] + B[1]*gradB[matrixIndex(1,0)] + B[2]*gradB[matrixIndex(2,0)] );
      B_gradient[1] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,1)] + B[1]*gradB[matrixIndex(1,1)] + B[2]*gradB[matrixIndex(2,1)] );
      B_gradient[2] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,2)] + B[1]*gradB[matrixIndex(1,2)] + B[2]*gradB[matrixIndex(2,2)] );
      Real V_curvature[3];
      crossProduct(B,B_gradient,V_curvature);
      
      // Calculate total drift velocity:
      Real V_drift[3];
      V_drift[0] = V_electric[0] + V_gradient[0] + V_curvature[0];
      V_drift[1] = V_electric[1] + V_gradient[1] + V_curvature[1];
      V_drift[2] = V_electric[2] + V_gradient[2] + V_curvature[2];
      
      // Calculate GC coordinate velocity:
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 acceleration[sep::particle::XCRD] = V_GC_par[0] + V_drift[0];
	 acceleration[sep::particle::YCRD] = V_GC_par[1] + V_drift[1];
	 acceleration[sep::particle::ZCRD] = V_GC_par[2] + V_drift[2];
	 break;
       case sep::CARTESIAN:
	 acceleration[sep::particle::XCRD] = V_GC_par[0] + V_drift[0];
	 acceleration[sep::particle::YCRD] = V_GC_par[1] + V_drift[1];
	 acceleration[sep::particle::ZCRD] = V_GC_par[2] + V_drift[2];
	 break;
       case sep::CYLINDRICAL:
	 I = static_cast<int32_t>(state[sep::particle::XCRD]);
	 R = sim->x_crds_node[I] + (state[sep::particle::XCRD] - I) * sim->dx_cell[I];
	 acceleration[sep::particle::XCRD] = V_GC_par[0] + V_drift[0];
	 acceleration[sep::particle::YCRD] = (V_GC_par[1] + V_drift[1]) / R;
	 acceleration[sep::particle::ZCRD] = V_GC_par[2] + V_drift[2];
	 break;
       case sep::SPHERICAL:
	 I = static_cast<int32_t>(state[sep::particle::XCRD]);
	 if (I < 0) {
	    I = 0;
	 } else if (I >= (int32_t)sim->x_blocks*block::WIDTH_X) {
	    I = sim->x_blocks*block::WIDTH_X-1;
	 } else {

	 }

	 J = static_cast<int32_t>(state[sep::particle::YCRD]);
	 R     = sim->x_crds_node[I] + (state[sep::particle::XCRD] - I) * sim->dx_cell[I];
	 THETA = sim->y_crds_node[J] + (state[sep::particle::YCRD] - J) * sim->dy_cell[J];
	 acceleration[sep::particle::XCRD] =  V_GC_par[0] + V_drift[0];
	 acceleration[sep::particle::YCRD] = (V_GC_par[1] + V_drift[1]) / R;
	 acceleration[sep::particle::ZCRD] = (V_GC_par[2] + V_drift[2]) / ( R*sin(THETA) );
	 break;
       default:
	 acceleration[sep::particle::XCRD] = V_GC_par[0] + V_drift[0];
	 acceleration[sep::particle::YCRD] = V_GC_par[1] + V_drift[1];
	 acceleration[sep::particle::ZCRD] = V_GC_par[2] + V_drift[2];
	 break;
      }
      
      #ifndef NDEBUG
      bool ok = true;
      if (acceleration[sep::particle::XCRD] != acceleration[sep::particle::XCRD]) ok = false;
      if (acceleration[sep::particle::YCRD] != acceleration[sep::particle::YCRD]) ok = false;
      if (acceleration[sep::particle::ZCRD] != acceleration[sep::particle::ZCRD]) ok = false;
      if (ok == false) {
	 std::cerr << "(SEP ACC PARTICLE) ERROR: NAN acceleration" << std::endl;
	 std::cerr << "\t E:   " << E[0] << '\t' << E[1] << '\t' << E[2] << std::endl;
	 std::cerr << "\t B:   " << B[0] << '\t' << B[1] << '\t' << B[2] << std::endl;
	 std::cerr << "\t POS: " << state[0] << '\t' << state[1] << '\t' << state[2] << std::endl;
	 std::cerr << "\t ACC: " << acceleration[0] << '\t' << acceleration[1] << '\t' << acceleration[2] << std::endl;
	 std::cerr << "\t V_PAR: " << state[sep::particle::V_PAR] << " \t MAGN.MOM.: " << state[sep::particle::MU] << std::endl;
	 exit(1);
      }

      I = static_cast<int32_t>(state[sep::particle::XCRD]);
      J = static_cast<int32_t>(state[sep::particle::YCRD]);
      K = static_cast<int32_t>(state[sep::particle::ZCRD]);
      if (acceleration[0]/sim->dx_cell[I]*sim->dt >= 1.0) ok = false;
      if (acceleration[1]/sim->dy_cell[J]*sim->dt >= 1.0) ok = false;
      if (acceleration[2]/sim->dz_cell[K]*sim->dt >= 1.0) ok = false;
      if (ok == false) {
	 std::stringstream ss;
	 ss << "(SEP ACC PARTICLE) ERROR: CFL violated" << std::endl;
	 ss << "\t POS     : " << state[0] << ' ' << state[1] << ' ' << state[2] << std::endl;
	 ss << "\t V_GC_par: " << V_GC_par[0]/sim->dx_cell[I]*sim->dt << ' ' << V_GC_par[1]/sim->dy_cell[J]*sim->dt << ' ' << V_GC_par[2]/sim->dz_cell[K]*sim->dt << ' ' << std::endl;
	 ss << "\t V_drift : " << V_drift[0]/sim->dx_cell[I]*sim->dt << ' ' << V_drift[1]/sim->dy_cell[J]*sim->dt << ' ' << V_drift[2]/sim->dz_cell[K]*sim->dt << ' ' << std::endl;
	 ss << "\t V_par   : " << state[sep::particle::V_PAR] << " dt_max: " << sim->dt << " step " << sim->timestep << std::endl;
	 ss << "\t MU      : " << state[sep::particle::MU] << std::endl;
	 ss << "\t t       : " << sim->t << std::endl;
	 ss << "\t E       : " << E[0] << '\t' << E[1] << '\t' << E[2] << std::endl;
	 ss << "\t B       : " << B[0] << '\t' << B[1] << '\t' << B[2] << std::endl;
	 ss << "\t gradB   : " << gradB[0] << '\t' << gradB[1] << '\t' << gradB[2] << std::endl;
	 ss << "\t         : " << gradB[3] << '\t' << gradB[4] << '\t' << gradB[5] << std::endl;
	 ss << "\t         : " << gradB[6] << '\t' << gradB[7] << '\t' << gradB[8] << std::endl;
	 ss << "\t V_e     : " << V_electric[0] << '\t' << V_electric[1] << '\t' << V_electric[2] << std::endl;
	 ss << "\t V_grad  : " << V_gradient[0] << '\t' << V_gradient[1] << '\t' << V_gradient[2] << std::endl;
	 ss << "\t V_curv  : " << V_curvature[0] << '\t' << V_curvature[1] << '\t' << V_curvature[2] << std::endl;
	 ss << "\t Indices : " << I << ' ' << J << ' ' << K << std::endl;
	 ss << "\t dx      : " << sim->dx_cell[I] << '\t' << sim->dy_cell[J] << '\t' << sim->dz_cell[K] << std::endl;
	 std::cerr << ss.str();
	 exit(1);
      }      
      #endif
   }

   template<class SPECIES,class PARTICLE> inline
   bool ParticleAccelerator<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      initialized = true;
      
      this->sim = &sim;
      this->simClasses = &simClasses;
      return initialized;
   }
   
} // namespace sep
   
#endif
