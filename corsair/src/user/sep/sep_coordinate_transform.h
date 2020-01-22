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

#ifndef SEP_COORDINATE_TRANSFORM_H
#define SEP_COORDINATE_TRANSFORM_H

#include <cmath>

#include <definitions.h>
#include <linear_algebra.h>
#include "sep_simcontrol.h"

namespace sep {

   template<typename REAL>
   void getTransformMatrix(const REAL* position,REAL* matrix,sep::CoordinateSystem original,sep::CoordinateSystem target);

   template<typename REAL> void getTransformMatrixCartesianToCartesian(const REAL* position,REAL* matrix);
   template<typename REAL> void getTransformMatrixCylindricalToCartesian(const REAL* position,REAL* matrix);
   template<typename REAL> void getTransformMatrixSphericalToCartesian(const REAL* position,REAL* matrix);

   template<typename REAL> void transformPositionCartesianToSpherical(const REAL* posIn,REAL* posOut);
   template<typename REAL> void transformPositionSphericalToCartesian(const REAL* posIn,REAL* posOut);
   
   template<typename REAL>
   void transformVector(const REAL* position,const REAL* vectorIn,REAL* vectorOut,sep::CoordinateSystem original,sep::CoordinateSystem target);
   
   template<typename REAL> void transformVectorCartesianToCylindrical(const REAL* position,const REAL* vectorIn,REAL* vectorOut);
   template<typename REAL> void transformVectorCartesianToSpherical(const REAL* position,const REAL* vectorIn,REAL* vectorOut);
   template<typename REAL> void transformVectorCylindricalToCartesian(const REAL* position,const REAL* vectorIn,REAL* vectorOut);
   template<typename REAL> void transformVectorSphericalToCartesian(const REAL* position,const REAL* vectorIn,REAL* vectorOut);

   
   template<typename REAL> inline
   void getTransformMatrix(const REAL* RESTRICT position,REAL* RESTRICT matrix,sep::CoordinateSystem original,sep::CoordinateSystem target) {
      switch (original) {
       case UNKNOWN:
	 std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	 break;
       case CARTESIAN:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    getTransformMatrixCartesianToCartesian(position,matrix);
	    break;
	  case CYLINDRICAL:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  case SPHERICAL:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  default:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       case CYLINDRICAL:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    getTransformMatrixCylindricalToCartesian(position,matrix);
	    break;
	  case CYLINDRICAL:
	    getTransformMatrixCartesianToCartesian(position,matrix);
	    break;
	  case SPHERICAL:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  default:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       case SPHERICAL:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    getTransformMatrixSphericalToCartesian(position,matrix);
	    break;
	  case CYLINDRICAL:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	  case SPHERICAL:
	    getTransformMatrixCartesianToCartesian(position,matrix);
	    break;
	  default:
	    std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       default:
	 std::cerr << "(SEP COORDINATE TRANSFORM) ERROR: Unknown or unsupported coordinate system" << std::endl; exit(1);
	 break;
      }
   }
   
   template<typename REAL> inline 
   void getTransformMatrixCartesianToCartesian(const REAL* RESTRICT position,REAL* RESTRICT matrix) {
      matrix[matrixIndex(0,0)] = 1.0;
      matrix[matrixIndex(0,1)] = 0.0;
      matrix[matrixIndex(0,2)] = 0.0;
      
      matrix[matrixIndex(1,0)] = 0.0;
      matrix[matrixIndex(1,1)] = 1.0;
      matrix[matrixIndex(1,2)] = 0.0;
      
      matrix[matrixIndex(2,0)] = 0.0;
      matrix[matrixIndex(2,1)] = 0.0;
      matrix[matrixIndex(2,2)] = 1.0;
   }
   
   template<typename REAL> 
   inline void getTransformMatrixCylindricalToCartesian(const REAL* RESTRICT position,REAL* RESTRICT matrix) {
      matrix[matrixIndex(0,0)] = cos(position[1]);
      matrix[matrixIndex(0,1)] = -sin(position[1]);
      matrix[matrixIndex(0,2)] = 0.0;
      
      matrix[matrixIndex(1,0)] = sin(position[1]);
      matrix[matrixIndex(1,1)] = cos(position[1]);
      matrix[matrixIndex(1,2)] = 0.0;
      
      matrix[matrixIndex(2,0)] = 0.0;
      matrix[matrixIndex(2,1)] = 0.0;
      matrix[matrixIndex(2,2)] = 1.0;
   }
   
   template<typename REAL> inline 
   void getTransformMatrixSphericalToCartesian(const REAL* RESTRICT position,REAL* RESTRICT matrix) {
      matrix[matrixIndex(0,0)] = sin(position[1])*cos(position[2]);
      matrix[matrixIndex(0,1)] = cos(position[1])*cos(position[2]);
      matrix[matrixIndex(0,2)] = -sin(position[2]);
      
      matrix[matrixIndex(1,0)] = sin(position[1])*sin(position[2]);
      matrix[matrixIndex(1,1)] = cos(position[1])*sin(position[2]);
      matrix[matrixIndex(1,2)] = cos(position[2]);
      
      matrix[matrixIndex(2,0)] = cos(position[1]);
      matrix[matrixIndex(2,1)] = -sin(position[1]);
      matrix[matrixIndex(2,2)] = 0.0;
   }

   template<typename REAL> inline
   void transformPositionCartesianToSpherical(const REAL* RESTRICT posIn,REAL* RESTRICT posOut) {
      posOut[0] = vectorMagnitude<3>(posIn);
      posOut[1] = acos(posIn[2] / posOut[0]);
      posOut[2] = atan(posIn[1] / posIn[0]);
   }
   
   template<typename REAL> inline
   void transformPositionSphericalToCartesian(const REAL* RESTRICT posIn,REAL* RESTRICT posOut) {
      posOut[0] = posIn[0] * sin(posIn[1]) * cos(posIn[2]);
      posOut[1] = posIn[0] * sin(posIn[1]) * sin(posIn[2]);
      posOut[2] = posIn[0] * cos(posIn[1]);
   }

   template<typename REAL> inline
   void transformVector(const REAL* RESTRICT position,const REAL* RESTRICT vectorIn,REAL* RESTRICT vectorOut,
			sep::CoordinateSystem original,sep::CoordinateSystem target) {
      switch (original) {
       case UNKNOWN:
	 std::cerr << "transformVector ERROR: Unknown source coordinate system" << std::endl; exit(1);
	 break;
       case CARTESIAN:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    for (int i=0; i<3; ++i) vectorOut[i] = vectorIn[i];
	    break;
	  case CYLINDRICAL:
	    transformVectorCartesianToCylindrical(position,vectorIn,vectorOut);
	    break;
	  case SPHERICAL:
	    transformVectorCartesianToSpherical(position,vectorIn,vectorOut);
	    break;
	  default:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       case CYLINDRICAL:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    transformVectorCylindricalToCartesian(position,vectorIn,vectorOut);
	    break;
	  case CYLINDRICAL:
	    for (int i=0; i<3; ++i) vectorOut[i] = vectorIn[i];
	    break;
	  case SPHERICAL:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	  default:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       case SPHERICAL:
	 switch (target) {
	  case UNKNOWN:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	  case CARTESIAN:
	    transformVectorSphericalToCartesian(position,vectorIn,vectorOut);
	    break;
	  case CYLINDRICAL:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	  case SPHERICAL:
	    for (int i=0; i<3; ++i) vectorOut[i] = vectorIn[i];
	    break;
	  default:
	    std::cerr << "transformVector ERROR: Unknown target coordinate system" << std::endl; exit(1);
	    break;
	 }
	 break;
       default:
	 std::cerr << "transformVector ERROR: Unknown source coordinate system" << std::endl; exit(1);
	 break;
      }
   }

   template<typename REAL> inline
   void transformVectorCartesianToCylindrical(const REAL* RESTRICT position,const REAL* RESTRICT vectorIn,REAL* RESTRICT vectorOut) {
      const Real rho = sqrt(position[0]*position[0]+position[1]*position[1]);
      
      const Real cos_phi = position[0] / rho;
      const Real sin_phi = position[1] / rho;
      
      vectorOut[0] =  vectorIn[0]*cos_phi + vectorIn[1]*sin_phi;
      vectorOut[1] = -vectorIn[0]*sin_phi + vectorIn[1]*cos_phi;
      vectorOut[2] =  vectorIn[2];
   }
      
   template<typename REAL> inline
   void transformVectorCartesianToSpherical(const REAL* RESTRICT position,const REAL* RESTRICT vectorIn,REAL* RESTRICT vectorOut) {
      const Real rho2 = position[0]*position[0]+position[1]*position[1];
      const Real r    = sqrt(rho2 + position[2]*position[2]);
      const Real rho  = sqrt(rho2);
      
      const Real cos_theta = position[2] / r;
      const Real sin_theta = rho / r;
      const Real cos_phi   = position[0] / rho;
      const Real sin_phi   = position[1] / rho;
      
      vectorOut[0] = vectorIn[0]*sin_theta*cos_phi + vectorIn[1]*sin_theta*sin_phi + vectorIn[2]*cos_theta;
      vectorOut[1] = vectorIn[0]*cos_theta*cos_phi + vectorIn[1]*cos_theta*sin_phi - vectorIn[2]*sin_theta;
      vectorOut[2] = -vectorIn[0]*sin_phi          + vectorIn[1]*cos_phi;
   }
   
   template<typename REAL> inline
   void transformVectorCylindricalToCartesian(const REAL* RESTRICT position,const REAL* RESTRICT vectorIn,REAL* RESTRICT vectorOut) {
      vectorOut[0] = vectorIn[0]*cos(position[1]) - vectorIn[1]*sin(position[1]);
      vectorOut[1] = vectorIn[0]*sin(position[1]) + vectorIn[1]*cos(position[1]);
      vectorOut[2] = vectorIn[2];
   }
   
   template<typename REAL> inline
   void transformVectorSphericalToCartesian(const REAL* RESTRICT position,const REAL* RESTRICT vectorIn,REAL* RESTRICT vectorOut) {
      vectorOut[0] = vectorIn[0]*sin(position[1])*cos(position[2]) + vectorIn[1]*cos(position[1])*cos(position[2]) - vectorIn[2]*sin(position[2]);
      vectorOut[1] = vectorIn[0]*sin(position[1])*sin(position[2]) + vectorIn[1]*cos(position[1])*sin(position[2]) + vectorIn[2]*cos(position[2]);
      vectorOut[2] = vectorIn[0]*cos(position[1])                  - vectorIn[1]*sin(position[1]);
   }
   
} // namespace sep

#endif
