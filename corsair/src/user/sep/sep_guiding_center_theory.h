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

#ifndef SEP_GUIDING_CENTER_THEORY_H
#define SEP_GUIDING_CENTER_THEORY_H

#include <cmath>
#include <linear_algebra.h>

#include "sep_particle_definition.h"

namespace sep {

   template<class SPECIES,class PARTICLE>
   void calculateDriftVelocity(const SPECIES& species,const PARTICLE& particle,const Real E[3],const Real B[3],const Real gradB[9],Real V_drift[3]);
   
   /** Calculate guiding center drift velocity. Currently electric, gradient,
    * and curvature drifts are included.
    * @param species Particle species, needed for charge and charge-to-mass ratio.
    * @param particle Particle state (parallel speed and magnetic moment).
    * @param E Electric field at particle position.
    * @param B Magnetic field at particle position.
    * @param gradB Gradient B tensor at particle position.
    * @param V_drift Calculated drift velocity is written here.
    */
   template<class SPECIES,class PARTICLE>
   void calculateDriftVelocity(const SPECIES& species,const PARTICLE& particle,const Real E[3],const Real B[3],const Real gradB[9],Real V_drift[3]) {
      // Calculate ExB drift velocity:
      Real V_electric[3];
      crossProduct(E,B,V_electric);
      const Real B_mag2 = vectorMagnitude2<3>(B);
      V_electric[0] /= B_mag2;
      V_electric[1] /= B_mag2;
      V_electric[2] /= B_mag2;
      
      // Calculate gradient of B:
      Real B_gradient[3];
      const Real B_mag = sqrt(B_mag2);
      B_gradient[0] = (B[0]*gradB[matrixIndex(0,0)] + B[1]*gradB[matrixIndex(0,1)] + B[2]*gradB[matrixIndex(0,2)]) / B_mag;
      B_gradient[1] = (B[0]*gradB[matrixIndex(1,0)] + B[1]*gradB[matrixIndex(1,1)] + B[2]*gradB[matrixIndex(1,2)]) / B_mag;
      B_gradient[2] = (B[0]*gradB[matrixIndex(2,0)] + B[1]*gradB[matrixIndex(2,1)] + B[2]*gradB[matrixIndex(2,2)]) / B_mag;

      // Calculate gradient drift velocity:
      Real V_gradient[3];
      crossProduct(B,B_gradient,V_gradient);
      V_gradient[0] *= (particle.state[sep::particle::MU]/species.charge/B_mag2);
      V_gradient[1] *= (particle.state[sep::particle::MU]/species.charge/B_mag2);
      V_gradient[2] *= (particle.state[sep::particle::MU]/species.charge/B_mag2);
      
      // Calculate curvature drift velocity (B_gradient used as temporary array here):
      const Real curvatureFactor = particle.state[sep::particle::V_PAR]*particle.state[sep::particle::V_PAR]/species.q_per_m/B_mag2/B_mag2;
      B_gradient[0] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,0)] + B[1]*gradB[matrixIndex(1,0)] + B[2]*gradB[matrixIndex(2,0)] );
      B_gradient[1] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,1)] + B[1]*gradB[matrixIndex(1,1)] + B[2]*gradB[matrixIndex(2,1)] );
      B_gradient[2] = curvatureFactor * ( B[0]*gradB[matrixIndex(0,2)] + B[1]*gradB[matrixIndex(1,2)] + B[2]*gradB[matrixIndex(2,2)] );
      Real V_curvature[3];
      crossProduct(B,B_gradient,V_curvature);
      
      // Store total drift velocity to V_drift:
      for (int i=0; i<3; ++i) V_drift[i] = V_electric[i] + V_gradient[i] + V_curvature[i];
   }
   
} // namespace sep

#endif
