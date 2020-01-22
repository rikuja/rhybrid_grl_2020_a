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

#ifndef SEP_ACCUMULATION_STRETCHED_H
#define SEP_ACCUMULATION_STRETCHED_H

#include <stdint.h>
#include <simulationclasses.h>

#include "sep_simcontrol.h"

namespace sep {

   extern sep::SimControl simControl;
   
   template<typename REAL>
   void accumScalarCentroidLogisticNGP_3D(REAL* array,const REAL* const pos,REAL value);
   template<typename REAL>
   void accumScalarCentroidLogisticCIC_3D(REAL* array,const REAL* const pos,REAL value);
   template<typename REAL>
   void accumScalarCentroidLogisticTSC_3D(REAL* array,const REAL* const pos,REAL value);
   template<typename REAL>
   bool accumScalarCentroidLogisticNGP_4D(REAL* array,const REAL* const pos,REAL value,int SIZE);
   template<typename REAL>
   bool accumScalarCentroidLogisticCIC_4D(REAL* array,const REAL* const pos,REAL value,int SIZE);
   template<typename REAL>
   bool accumScalarCentroidLogisticTSC_4D(REAL* array,const REAL* const pos,REAL value,int SIZE);

   template<typename REAL>
   void accumMultiplyScalarCentroidLogicalTSC_4D(REAL* RESTRICT array,const REAL* RESTRICT array2,
						 const int32_t* RESTRICT const indices,
						 const REAL* RESTRICT const sf,REAL value,int SIZE);

   template<typename REAL>
   void accumScalarCentroidLogicalNGP_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE);
   template<typename REAL>
   void accumScalarCentroidLogicalCIC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE);
   template<typename REAL>
   void accumScalarCentroidLogicalTSC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE);
   
   // TEST
   template<typename REAL>
   void accumScalarCentroidLogicalShockTSC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,
					      const REAL* RESTRICT const sf,const int32_t* cellRegions,
					      int32_t region,REAL value,int SIZE);
   // END TEST
   
   template<typename REAL>
   void addValues4D(SimulationClasses* simClasses,pargrid::CellID blockLID,const REAL* const source,REAL* dest,uint32_t SIZE);

   template<typename REAL>
   void getDerivativeShapeFactorsCIC_4D(const REAL* const pos,int32_t* indices,REAL* shapeFactors);
   template<typename REAL>
   void getDerivativeShapeFactorsTSC_4D(const REAL* const pos,int32_t* indices,REAL* shapeFactors);

   template<typename REAL>
   void getShapeFactorsNGP_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors);
   template<typename REAL>
   void getShapeFactorsCIC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors);
   template<typename REAL>
   void getShapeFactorsTSC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors);

   template<typename REAL>
   REAL interpolateScalarNGP_4D(const REAL* array,const int32_t* const indices,const REAL* const sf,int SIZE);
   template<typename REAL> 
   REAL interpolateScalarCIC_4D(const REAL* array,const int32_t* const indices,const REAL* const sf,int SIZE);
   template<typename REAL> 
   REAL interpolateScalarTSC_4D(const REAL* array,const int32_t* const indices,const REAL* const sf,int SIZE);
   
   template<typename REAL>
   REAL interpolateScalarDerivLogicalTSC_4D(const REAL* array,const REAL* const pos,int SIZE);

   template<typename REAL>
   REAL interpolateScalarLogicalNGP_4D(const REAL* array,const REAL* const pos,int SIZE);
   template<typename REAL>
   REAL interpolateScalarLogicalCIC_4D(const REAL* array,const REAL* const pos,int SIZE);
   template<typename REAL>
   REAL interpolateScalarLogicalTSC_4D(const REAL* array,const REAL* const pos,int SIZE);
   
   template<typename REAL>
   void loadScalar4D(SimulationClasses* simClasses,pargrid::CellID blockLID,const REAL* const source,REAL* dest,uint32_t SIZE);
   
   template<typename REAL> inline
   void accumScalarCentroidLogisticNGP_3D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value) {
      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);

      #ifndef NDEBUG
         if (i < 1 || i > block::WIDTH_X) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
         if (j < 1 || j > block::WIDTH_Y) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
         if (k < 1 || k > block::WIDTH_Z) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
      #endif
      
      array[block::arrayIndex(i,j,k)] += value;
   }
   
   template<typename REAL> inline
   void accumScalarCentroidLogisticCIC_3D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const int32_t i = static_cast<int32_t>(X);
      const int32_t j = static_cast<int32_t>(Y);
      const int32_t k = static_cast<int32_t>(Z);
      const REAL W_X = X - i;
      const REAL W_Y = Y - j;
      const REAL W_Z = Z - k;
      
      #ifndef NDEBUG
         if (i < 0 || i > block::WIDTH_X) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
         if (j < 0 || j > block::WIDTH_Y) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
         if (k < 0 || k > block::WIDTH_Z) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds" << std::endl; exit(1);}
      #endif
      
      array[block::arrayIndex(i  ,j  ,k  )] += (1-W_X)*(1-W_Y)*(1-W_Z)*value;
      array[block::arrayIndex(i+1,j  ,k  )] +=   W_X  *(1-W_Y)*(1-W_Z)*value;
      array[block::arrayIndex(i  ,j+1,k  )] += (1-W_X)*  W_Y  *(1-W_Z)*value;
      array[block::arrayIndex(i+1,j+1,k  )] +=   W_X  *  W_Y  *(1-W_Z)*value;
      array[block::arrayIndex(i  ,j  ,k+1)] += (1-W_X)*(1-W_Y)*  W_Z  *value;
      array[block::arrayIndex(i+1,j  ,k+1)] +=   W_X  *(1-W_Y)*  W_Z  *value;
      array[block::arrayIndex(i  ,j+1,k+1)] += (1-W_X)*  W_Y  *  W_Z  *value;
      array[block::arrayIndex(i+1,j+1,k+1)] +=   W_X  *  W_Y  *  W_Z  *value;
   }
   
   template<typename REAL> inline
   void accumScalarCentroidLogisticTSC_3D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value) {
      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);

      #ifndef NDEBUG
         if (i < 1 || i > block::WIDTH_X) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds, x=" << pos[0] << std::endl; exit(1);}
         if (j < 1 || j > block::WIDTH_Y) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds, y=" << pos[1] << std::endl; exit(1);}
         if (k < 1 || k > block::WIDTH_Z) {std::cerr << i << ' ' << j << ' ' << k << " out of bounds, z=" << pos[2] << std::endl; exit(1);}
      #endif
      
      const REAL W_X_N = 0.5*(i+1-pos[0])*(i+1-pos[0]);
      const REAL W_X_P = 0.5*(pos[0] -i )*(pos[0] -i );
      const REAL W_X_C = 1-W_X_N-W_X_P;
      const REAL W_Y_N = 0.5*(j+1-pos[1])*(j+1-pos[1]);
      const REAL W_Y_P = 0.5*(pos[1] -j )*(pos[1] -j );
      const REAL W_Y_C = 1-W_Y_N-W_Y_P;
      const REAL W_Z_N = 0.5*(k+1-pos[2])*(k+1-pos[2]);
      const REAL W_Z_P = 0.5*(pos[2] -k )*(pos[2] -k );
      const REAL W_Z_C = 1-W_Z_N-W_Z_P;
      array[block::arrayIndex(i-1,j-1,k-1)] += W_X_N * W_Y_N * W_Z_N * value;
      array[block::arrayIndex(i  ,j-1,k-1)] += W_X_C * W_Y_N * W_Z_N * value;
      array[block::arrayIndex(i+1,j-1,k-1)] += W_X_P * W_Y_N * W_Z_N * value;
      array[block::arrayIndex(i-1,j  ,k-1)] += W_X_N * W_Y_C * W_Z_N * value;
      array[block::arrayIndex(i  ,j  ,k-1)] += W_X_C * W_Y_C * W_Z_N * value;
      array[block::arrayIndex(i+1,j  ,k-1)] += W_X_P * W_Y_C * W_Z_N * value;
      array[block::arrayIndex(i-1,j+1,k-1)] += W_X_N * W_Y_P * W_Z_N * value;
      array[block::arrayIndex(i  ,j+1,k-1)] += W_X_C * W_Y_P * W_Z_N * value;
      array[block::arrayIndex(i+1,j+1,k-1)] += W_X_P * W_Y_P * W_Z_N * value;
      
      array[block::arrayIndex(i-1,j-1,k  )] += W_X_N * W_Y_N * W_Z_C * value;
      array[block::arrayIndex(i  ,j-1,k  )] += W_X_C * W_Y_N * W_Z_C * value;
      array[block::arrayIndex(i+1,j-1,k  )] += W_X_P * W_Y_N * W_Z_C * value;
      array[block::arrayIndex(i-1,j  ,k  )] += W_X_N * W_Y_C * W_Z_C * value;
      array[block::arrayIndex(i  ,j  ,k  )] += W_X_C * W_Y_C * W_Z_C * value;
      array[block::arrayIndex(i+1,j  ,k  )] += W_X_P * W_Y_C * W_Z_C * value;
      array[block::arrayIndex(i-1,j+1,k  )] += W_X_N * W_Y_P * W_Z_C * value;
      array[block::arrayIndex(i  ,j+1,k  )] += W_X_C * W_Y_P * W_Z_C * value;
      array[block::arrayIndex(i+1,j+1,k  )] += W_X_P * W_Y_P * W_Z_C * value;
      
      array[block::arrayIndex(i-1,j-1,k+1)] += W_X_N * W_Y_N * W_Z_P * value;
      array[block::arrayIndex(i  ,j-1,k+1)] += W_X_C * W_Y_N * W_Z_P * value;
      array[block::arrayIndex(i+1,j-1,k+1)] += W_X_P * W_Y_N * W_Z_P * value;
      array[block::arrayIndex(i-1,j  ,k+1)] += W_X_N * W_Y_C * W_Z_P * value;
      array[block::arrayIndex(i  ,j  ,k+1)] += W_X_C * W_Y_C * W_Z_P * value;
      array[block::arrayIndex(i+1,j  ,k+1)] += W_X_P * W_Y_C * W_Z_P * value;
      array[block::arrayIndex(i-1,j+1,k+1)] += W_X_N * W_Y_P * W_Z_P * value;
      array[block::arrayIndex(i  ,j+1,k+1)] += W_X_C * W_Y_P * W_Z_P * value;
      array[block::arrayIndex(i+1,j+1,k+1)] += W_X_P * W_Y_P * W_Z_P * value;
   }

   template<typename REAL> inline
   bool accumScalarCentroidLogisticNGP_4D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value,int SIZE) {
      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);
      const int32_t l = static_cast<int32_t>(pos[3]);
      
      array[block::arrayIndex(i,j,k)*SIZE+l] += value;
      return true;
   }
   
   template<typename REAL> inline
   bool accumScalarCentroidLogisticCIC_4D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value,int SIZE) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const REAL L = pos[3] - 0.5;
      const int32_t i = static_cast<int32_t>(X);
      const int32_t j = static_cast<int32_t>(Y);
      const int32_t k = static_cast<int32_t>(Z);
      const int32_t l = static_cast<int32_t>(L);
      const REAL W_X = X - i;
      const REAL W_Y = Y - j;
      const REAL W_Z = Z - k;
      const REAL W_L = L - l;

      array[block::arrayIndex(i  ,j  ,k  )*SIZE+l  ] += (1-W_X)*(1-W_Y)*(1-W_Z)*(1-W_L)*value;
      array[block::arrayIndex(i  ,j  ,k  )*SIZE+l+1] += (1-W_X)*(1-W_Y)*(1-W_Z)*  W_L  *value;
      array[block::arrayIndex(i+1,j  ,k  )*SIZE+l  ] +=   W_X  *(1-W_Y)*(1-W_Z)*(1-W_L)*value;
      array[block::arrayIndex(i+1,j  ,k  )*SIZE+l+1] +=   W_X  *(1-W_Y)*(1-W_Z)*  W_L  *value;
      
      array[block::arrayIndex(i  ,j+1,k  )*SIZE+l  ] += (1-W_X)*  W_Y  *(1-W_Z)*(1-W_L)*value;
      array[block::arrayIndex(i  ,j+1,k  )*SIZE+l+1] += (1-W_X)*  W_Y  *(1-W_Z)*  W_L  *value;
      array[block::arrayIndex(i+1,j+1,k  )*SIZE+l  ] +=   W_X  *  W_Y  *(1-W_Z)*(1-W_L)*value;
      array[block::arrayIndex(i+1,j+1,k  )*SIZE+l+1] +=   W_X  *  W_Y  *(1-W_Z)*  W_L  *value;
      
      array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l  ] += (1-W_X)*(1-W_Y)*  W_Z  *(1-W_L)*value;
      array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l+1] += (1-W_X)*(1-W_Y)*  W_Z  *  W_L  *value;
      array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l  ] +=   W_X  *(1-W_Y)*  W_Z  *(1-W_L)*value;
      array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l+1] +=   W_X  *(1-W_Y)*  W_Z  *  W_L  *value;
      
      array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l  ] += (1-W_X)*  W_Y  *  W_Z  *(1-W_L)*value;
      array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l+1] += (1-W_X)*  W_Y  *  W_Z  *  W_L  *value;
      array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l  ] +=   W_X  *  W_Y  *  W_Z  *(1-W_L)*value;
      array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l+1] +=   W_X  *  W_Y  *  W_Z  *  W_L  *value;
      return true;
   }

   template<typename REAL> inline
   bool accumScalarCentroidLogisticTSC_4D(REAL* RESTRICT array,const REAL* RESTRICT const pos,REAL value,int SIZE) {
      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);
      const int32_t l = static_cast<int32_t>(pos[3]);
      
      const REAL W_X_N = 0.5*(i+1-pos[0])*(i+1-pos[0]);
      const REAL W_X_P = 0.5*(pos[0] -i )*(pos[0] -i );
      const REAL W_X_C = 1-W_X_N-W_X_P;
      
      const REAL W_Y_N = 0.5*(j+1-pos[1])*(j+1-pos[1]);
      const REAL W_Y_P = 0.5*(pos[1] -j )*(pos[1] -j );
      const REAL W_Y_C = 1-W_Y_N-W_Y_P;
      
      const REAL W_Z_N = 0.5*(k+1-pos[2])*(k+1-pos[2]);
      const REAL W_Z_P = 0.5*(pos[2] -k )*(pos[2] -k );
      const REAL W_Z_C = 1-W_Z_N-W_Z_P;
      
      const REAL W_L_N = 0.5*(l+1-pos[3])*(l+1-pos[3]);
      const REAL W_L_P = 0.5*(pos[3] -l )*(pos[3] -l );
      const REAL W_L_C = 1-W_L_N-W_L_P;
      

      #ifndef NDEBUG
      // Check that shape factors are valid:
      bool ok = true;
      if (W_X_N < 0.0 || W_X_N > 1.0) ok = false;
      if (W_X_P < 0.0 || W_X_P > 1.0) ok = false;
      if (W_Y_N < 0.0 || W_Y_N > 1.0) ok = false;
      if (W_Y_P < 0.0 || W_Y_P > 1.0) ok = false;
      if (W_Z_N < 0.0 || W_Z_N > 1.0) ok = false;
      if (W_Z_P < 0.0 || W_Z_P > 1.0) ok = false;
      if (W_L_N < 0.0 || W_L_N > 1.0) ok = false;
      if (W_L_P < 0.0 || W_L_P > 1.0) ok = false;
      if (ok == false) {
	 std::cerr << "WEIGHTS: " << W_X_N << '\t' << W_X_P << '\t' << W_Y_N << '\t' << W_Y_P;
	 std::cerr << '\t' << W_Z_N << '\t' << W_Z_P << '\t' << W_L_N << '\t' << W_L_P << std::endl;
      }
      
      // Check that indices are within bounds:
      if (i < 1 || i > block::WIDTH_X) ok = false;
      if (j < 1 || j > block::WIDTH_Y) ok = false;
      if (k < 1 || k > block::WIDTH_Z) ok = false;
      if (l < 1 || l > SIZE-1) ok = false;
      if (ok == false) {
	 std::cerr << "Indices " << i << ' ' << j << ' ' << k << ' ' << l << " OUT OF BOUNDS!" << std::endl;
      }
      #endif
      // K - 1, J - 1
      array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l-1] += W_X_N * W_Y_N * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l  ] += W_X_N * W_Y_N * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l+1] += W_X_N * W_Y_N * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l-1] += W_X_C * W_Y_N * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l  ] += W_X_C * W_Y_N * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l+1] += W_X_C * W_Y_N * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l-1] += W_X_P * W_Y_N * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l  ] += W_X_P * W_Y_N * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l+1] += W_X_P * W_Y_N * W_Z_N * W_L_P * value;
      
      // J + 0
      array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l-1] += W_X_N * W_Y_C * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l  ] += W_X_N * W_Y_C * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l+1] += W_X_N * W_Y_C * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l-1] += W_X_C * W_Y_C * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l  ] += W_X_C * W_Y_C * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l+1] += W_X_C * W_Y_C * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l-1] += W_X_P * W_Y_C * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l  ] += W_X_P * W_Y_C * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l+1] += W_X_P * W_Y_C * W_Z_N * W_L_P * value;
      
      // J + 1
      array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l-1] += W_X_N * W_Y_P * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l  ] += W_X_N * W_Y_P * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l+1] += W_X_N * W_Y_P * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l-1] += W_X_C * W_Y_P * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l  ] += W_X_C * W_Y_P * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l+1] += W_X_C * W_Y_P * W_Z_N * W_L_P * value;
      
      array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l-1] += W_X_P * W_Y_P * W_Z_N * W_L_N * value;
      array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l  ] += W_X_P * W_Y_P * W_Z_N * W_L_C * value;
      array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l+1] += W_X_P * W_Y_P * W_Z_N * W_L_P * value;
      
      // K + 0, J - 1
      array[block::arrayIndex(i-1,j-1,k  )*SIZE+l-1] += W_X_N * W_Y_N * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i-1,j-1,k  )*SIZE+l  ] += W_X_N * W_Y_N * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i-1,j-1,k  )*SIZE+l+1] += W_X_N * W_Y_N * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i  ,j-1,k  )*SIZE+l-1] += W_X_C * W_Y_N * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i  ,j-1,k  )*SIZE+l  ] += W_X_C * W_Y_N * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i  ,j-1,k  )*SIZE+l+1] += W_X_C * W_Y_N * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i+1,j-1,k  )*SIZE+l-1] += W_X_P * W_Y_N * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i+1,j-1,k  )*SIZE+l  ] += W_X_P * W_Y_N * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i+1,j-1,k  )*SIZE+l+1] += W_X_P * W_Y_N * W_Z_C * W_L_P * value;
      
      // J + 0
      array[block::arrayIndex(i-1,j  ,k  )*SIZE+l-1] += W_X_N * W_Y_C * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i-1,j  ,k  )*SIZE+l  ] += W_X_N * W_Y_C * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i-1,j  ,k  )*SIZE+l+1] += W_X_N * W_Y_C * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i  ,j  ,k  )*SIZE+l-1] += W_X_C * W_Y_C * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i  ,j  ,k  )*SIZE+l  ] += W_X_C * W_Y_C * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i  ,j  ,k  )*SIZE+l+1] += W_X_C * W_Y_C * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i+1,j  ,k  )*SIZE+l-1] += W_X_P * W_Y_C * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i+1,j  ,k  )*SIZE+l  ] += W_X_P * W_Y_C * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i+1,j  ,k  )*SIZE+l+1] += W_X_P * W_Y_C * W_Z_C * W_L_P * value;
      
      // J + 1
      array[block::arrayIndex(i-1,j+1,k  )*SIZE+l-1] += W_X_N * W_Y_P * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i-1,j+1,k  )*SIZE+l  ] += W_X_N * W_Y_P * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i-1,j+1,k  )*SIZE+l+1] += W_X_N * W_Y_P * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i  ,j+1,k  )*SIZE+l-1] += W_X_C * W_Y_P * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i  ,j+1,k  )*SIZE+l  ] += W_X_C * W_Y_P * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i  ,j+1,k  )*SIZE+l+1] += W_X_C * W_Y_P * W_Z_C * W_L_P * value;
      
      array[block::arrayIndex(i+1,j+1,k  )*SIZE+l-1] += W_X_P * W_Y_P * W_Z_C * W_L_N * value;
      array[block::arrayIndex(i+1,j+1,k  )*SIZE+l  ] += W_X_P * W_Y_P * W_Z_C * W_L_C * value;
      array[block::arrayIndex(i+1,j+1,k  )*SIZE+l+1] += W_X_P * W_Y_P * W_Z_C * W_L_P * value;
      
      // K + 1, J - 1
      array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l-1] += W_X_N * W_Y_N * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l  ] += W_X_N * W_Y_N * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l+1] += W_X_N * W_Y_N * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l-1] += W_X_C * W_Y_N * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l  ] += W_X_C * W_Y_N * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l+1] += W_X_C * W_Y_N * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l-1] += W_X_P * W_Y_N * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l  ] += W_X_P * W_Y_N * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l+1] += W_X_P * W_Y_N * W_Z_P * W_L_P * value;
      
      // J + 0
      array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l-1] += W_X_N * W_Y_C * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l  ] += W_X_N * W_Y_C * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l+1] += W_X_N * W_Y_C * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l-1] += W_X_C * W_Y_C * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l  ] += W_X_C * W_Y_C * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l+1] += W_X_C * W_Y_C * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l-1] += W_X_P * W_Y_C * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l  ] += W_X_P * W_Y_C * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l+1] += W_X_P * W_Y_C * W_Z_P * W_L_P * value;
      
      // J + 1
      array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l-1] += W_X_N * W_Y_P * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l  ] += W_X_N * W_Y_P * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l+1] += W_X_N * W_Y_P * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l-1] += W_X_C * W_Y_P * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l  ] += W_X_C * W_Y_P * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l+1] += W_X_C * W_Y_P * W_Z_P * W_L_P * value;
      
      array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l-1] += W_X_P * W_Y_P * W_Z_P * W_L_N * value;
      array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l  ] += W_X_P * W_Y_P * W_Z_P * W_L_C * value;
      array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l+1] += W_X_P * W_Y_P * W_Z_P * W_L_P * value;
      
      #ifndef NDEBUG
         return ok;
      #else
         return true;
      #endif
   }
   
   template<typename REAL> inline
   void accumScalarCentroidLogicalNGP_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
            
      array[block::arrayIndex(i,j,k)*SIZE+l] += sf[0]*value;
   }

   template<typename REAL> inline
   void accumScalarCentroidLogicalCIC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
      
      for (int k_off=0; k_off<2; ++k_off) for (int j_off=0; j_off<2; ++j_off) for (int i_off=0; i_off<2; ++i_off) for (int l_off=0; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off] * sf[2+j_off] * sf[4+k_off] * sf[6+l_off];
	 array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] += shapeFactor * value;
      }
   }

   template<typename REAL> inline
   void accumScalarCentroidLogicalTSC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,REAL value,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
      
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) for (int l_off=-1; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off+1] * sf[3+j_off+1] * sf[6+k_off+1] * sf[9+l_off+1];
	 #ifndef NDEBUG
	    const int32_t index = block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off;
	    if (index >= block::SIZE_PLUS_ONE_LAYER*SIZE) {
	       std::cerr << "(SEP ACCUM STRETCHED) ERROR: Invalid index in accumScalarCentroidLogicalTSC_4D" << std::endl;
	       std::cerr << "\t" << indices[0] << ' ' << indices[1] << ' ' << indices[2] << ' ' << indices[3] << std::endl;
	       exit(1);
	    }
	    array[index] += shapeFactor * value;
	 #else
	     array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] += shapeFactor * value;
	 #endif
      }
   }

   // TEST
   template<typename REAL> inline
   void accumScalarCentroidLogicalShockTSC_4D(REAL* RESTRICT array,const int32_t* RESTRICT const indices,
					      const REAL* RESTRICT const sf,const int32_t* cellRegions,
					      int32_t region,REAL value,int SIZE) {
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];           

      Real remainder = 0.0;
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) for (int l_off=-1; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off+1] * sf[3+j_off+1] * sf[6+k_off+1] * sf[9+l_off+1];
	 
	 if (cellRegions[block::arrayIndex(i+i_off,j+j_off,k+k_off)] != region) {
	    remainder += shapeFactor;
	    continue;
	 }
	 
	 array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] += shapeFactor * value;
      }

      array[block::arrayIndex(i,j+j,k+k)*SIZE+l] += remainder * value;
   }
   // END TEST
   
   template<typename REAL> inline
   void accumMultiplyScalarCentroidLogicalTSC_4D(REAL* RESTRICT array,const REAL* RESTRICT array2,const int32_t* RESTRICT const indices,
						 const REAL* RESTRICT const sf,REAL value,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
      
      //Real sum = 0.0;
      
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) for (int l_off=-1; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off+1] * sf[3+j_off+1] * sf[6+k_off+1] * sf[9+l_off+1];
	 Real addedValue = shapeFactor * value * array2[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off];
	 array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] += addedValue;
	 //sum += addedValue;
      }
      
      //std::cerr << "\t accumulated total: " << sum << std::endl;
   }
   
   template<typename REAL>
   void addValues4D(SimulationClasses* simClasses,pargrid::CellID blockLID,REAL* RESTRICT source,REAL* RESTRICT dest,uint32_t SIZE) {
      REAL* RESTRICT destTmp = NULL;
      const int BLOCK_DATA_SIZE = block::SIZE*SIZE;
      const int SRC_SIZE = SIZE+2;
      
      pargrid::CellID nbrLID = simClasses->pargrid.invalidCellID();
      const pargrid::CellID* const RESTRICT nbrs = simClasses->pargrid.getCellNeighbourIDs(blockLID);
      
      // ***** ADD VALUES TO THIS BLOCK ***** //
      destTmp = dest + blockLID*BLOCK_DATA_SIZE;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	destTmp[block::index(i,j,k)*SIZE+l] += source[block::arrayIndex(i+1,j+1,k+1)*SRC_SIZE+l+1];

      // ***** ADD VALUES TO 6 FACE NEIGHBOURS ***** //
      // -x neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,j,k)*SIZE+l] += source[block::arrayIndex(0,j+1,k+1)*SRC_SIZE+l+1];
      }
      // +x neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,j,k)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)*SRC_SIZE+l+1];
      }
      // -y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,block::WIDTH_Y-1,k)*SIZE+l] += source[block::arrayIndex(i+1,0,k+1)*SRC_SIZE+l+1];
      }      
      // +y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,0,k)*SIZE+l] += source[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)*SRC_SIZE+l+1];
      }
      // -z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,j,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(i+1,j+1,0)*SRC_SIZE+l+1];
      }
      // +z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,j,0)*SIZE+l] += source[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      // ***** ADD VALUES TO 12 EDGE NEIGHBOURS ***** //
      // Add values to -x,-y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k)*SIZE+l] += source[block::arrayIndex(0,0,k+1)*SRC_SIZE+l+1];
      }
      // Add values to +x,-y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,block::WIDTH_Y-1,k)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,0,k+1)*SRC_SIZE+l+1];
      }
      // Add values to -x,+y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,0,k)*SIZE+l] += source[block::arrayIndex(0,block::WIDTH_Y+1,k+1)*SRC_SIZE+l+1];
      }
      // Add values to +x,+y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,0,k)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)*SRC_SIZE+l+1];
      }
      // Add values to -x,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(0,j+1,0)*SRC_SIZE+l+1];
      }
      // Add values to +x,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,j,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,j+1,0)*SRC_SIZE+l+1];
      }
      // Add values to -x,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,j,0)*SIZE+l] += source[block::arrayIndex(0,j+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      // Add values to +x,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,j,0)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      // Add values to -y,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(i+1,0,0)*SRC_SIZE+l+1];
      }
      // Add values to +y,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,0,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(i+1,block::WIDTH_Y+1,0)*SRC_SIZE+l+1];
      }
      // Add values to -y,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,block::WIDTH_Y-1,0)*SIZE+l] += source[block::arrayIndex(i+1,0,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      // Add values to +y,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(i,0,0)*SIZE+l] += source[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }

      // ***** ADD VALUES TO CORNER NEIGHBOURS ***** //
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(0,0,0)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,0,0)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(0,block::WIDTH_Y+1,0)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,0,block::WIDTH_Z-1)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0)*SIZE+l] += source[block::arrayIndex(0,0,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,block::WIDTH_Y-1,0)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(block::WIDTH_X-1,0,0)*SIZE+l] += source[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 destTmp = dest + nbrLID*BLOCK_DATA_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   destTmp[block::index(0,0,0)*SIZE+l] += source[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*SRC_SIZE+l+1];
      }
   }

   template<typename REAL> inline
   void getDerivativeShapeFactorsCIC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const REAL L = pos[3] - 0.5;
      
      indices[0] = static_cast<int32_t>(X);
      indices[1] = static_cast<int32_t>(Y);
      indices[2] = static_cast<int32_t>(Z);
      indices[3] = static_cast<int32_t>(L);
      
      const REAL W_X = X - indices[0];
      const REAL W_Y = Y - indices[1];
      const REAL W_Z = Z - indices[2];
      //const REAL W_L = L - indices[3];
      
      shapeFactors[0] = 1-W_X;
      shapeFactors[1] =  W_X;
      shapeFactors[2] = 1-W_Y;
      shapeFactors[3] =  W_Y;
      shapeFactors[4] = 1-W_Z;
      shapeFactors[5] =  W_Z;
      shapeFactors[6] = -1.0;
      shapeFactors[7] =  1.0;
   }
      
   template<typename REAL> inline
   void getDerivativeShapeFactorsTSC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors) {
      indices[0] = static_cast<int32_t>(pos[0]);
      indices[1] = static_cast<int32_t>(pos[1]);
      indices[2] = static_cast<int32_t>(pos[2]);
      indices[3] = static_cast<int32_t>(pos[3]);
      
      // X Shape Factors:
      shapeFactors[0]  = 0.5*(indices[0]+1-pos[0])*(indices[0]+1-pos[0]);
      shapeFactors[2]  = 0.5*(pos[0] -indices[0] )*(pos[0] -indices[0] );
      shapeFactors[1]  = 1-shapeFactors[0]-shapeFactors[2];
      // Y Shape Factors:
      shapeFactors[3]  = 0.5*(indices[1]+1-pos[1])*(indices[1]+1-pos[1]);
      shapeFactors[5]  = 0.5*(pos[1] - indices[1] )*(pos[1] - indices[1] );
      shapeFactors[4]  = 1-shapeFactors[3]-shapeFactors[5];
      // Z Shape Factors:
      shapeFactors[6]  = 0.5*(indices[2]+1-pos[2])*(indices[2]+1-pos[2]);
      shapeFactors[8]  = 0.5*(pos[2] -indices[2] )*(pos[2] -indices[2] );
      shapeFactors[7]  = 1-shapeFactors[6]-shapeFactors[8];
      // L Shape Factors:
      shapeFactors[9]  = pos[3] - indices[3]-1;
      shapeFactors[11] = pos[3] - indices[3];
      shapeFactors[10] = -shapeFactors[9]-shapeFactors[11];
   }
      
   template<typename REAL> inline
   void getShapeFactorsNGP_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors) {
      indices[0] = static_cast<int32_t>(pos[0]);
      indices[1] = static_cast<int32_t>(pos[1]);
      indices[2] = static_cast<int32_t>(pos[2]);
      indices[3] = static_cast<int32_t>(pos[3]);
      shapeFactors[0] = 1.0;
   }
   
   template<typename REAL> inline
   void getShapeFactorsCIC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const REAL L = pos[3] - 0.5;
      
      indices[0] = static_cast<int32_t>(X);
      indices[1] = static_cast<int32_t>(Y);
      indices[2] = static_cast<int32_t>(Z);
      indices[3] = static_cast<int32_t>(L);

      const REAL W_X = X - indices[0];
      const REAL W_Y = Y - indices[1];
      const REAL W_Z = Z - indices[2];
      const REAL W_L = L - indices[3];

      shapeFactors[0] = 1-W_X;
      shapeFactors[1] =  W_X;
      shapeFactors[2] = 1-W_Y;
      shapeFactors[3] =  W_Y;
      shapeFactors[4] = 1-W_Z;
      shapeFactors[5] =  W_Z;
      shapeFactors[6] = 1-W_L;
      shapeFactors[7] =  W_L;
   }
   
   template<typename REAL> inline
   void getShapeFactorsTSC_4D(const REAL* pos,int32_t* indices,REAL* shapeFactors) {
      indices[0] = static_cast<int32_t>(pos[0]);
      indices[1] = static_cast<int32_t>(pos[1]);
      indices[2] = static_cast<int32_t>(pos[2]);
      indices[3] = static_cast<int32_t>(pos[3]);

      #ifndef NDEBUG
      bool ok = true;
      if (indices[0] < 0 || indices[0] > block::WIDTH_X) ok = false;
      if (indices[1] < 0 || indices[1] > block::WIDTH_Y) ok = false;
      if (indices[2] < 0 || indices[2] > block::WIDTH_Z) ok = false;
      if (indices[3] < 0 || indices[3] > (int32_t)simControl.N_wavelengthMeshCells) ok = false;
      if (ok == false) {
	 std::cerr << "(SEP ACCUM STRETCHED) ERROR: Invalid indices in getShapeFactorsTSC_4D" << std::endl;
	 std::cerr << "\t" << indices[0] << ' ' << indices[1] << ' ' << indices[2] << ' ' << indices[3] << std::endl;
	 exit(1);
      }
      #endif
      
      shapeFactors[0]  = 0.5*(indices[0]+1-pos[0])*(indices[0]+1-pos[0]);
      shapeFactors[2]  = 0.5*(pos[0] -indices[0] )*(pos[0] -indices[0] );
      shapeFactors[1]  = 1-shapeFactors[0]-shapeFactors[2];
      
      shapeFactors[3]  = 0.5*(indices[1]+1-pos[1])*(indices[1]+1-pos[1]);
      shapeFactors[5]  = 0.5*(pos[1] - indices[1] )*(pos[1] - indices[1] );
      shapeFactors[4]  = 1-shapeFactors[3]-shapeFactors[5];
      
      shapeFactors[6]  = 0.5*(indices[2]+1-pos[2])*(indices[2]+1-pos[2]);
      shapeFactors[8]  = 0.5*(pos[2] -indices[2] )*(pos[2] -indices[2] );
      shapeFactors[7]  = 1-shapeFactors[6]-shapeFactors[8];
      
      shapeFactors[9]  = 0.5*(indices[3]+1-pos[3])*(indices[3]+1-pos[3]);
      shapeFactors[11] = 0.5*(pos[3] -indices[3] )*(pos[3] -indices[3] );
      shapeFactors[10] = 1-shapeFactors[9]-shapeFactors[11];
   }

   template<typename REAL> inline
   REAL interpolateScalarNGP_4D(const REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
      return array[block::arrayIndex(i,j,k)*SIZE+l];
   }
   
   template<typename REAL> inline
   REAL interpolateScalarCIC_4D(const REAL* RESTRICT array,const int32_t* const indices,const REAL* RESTRICT const sf,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];
      
      REAL value = 0.0;
      for (int k_off=0; k_off<2; ++k_off) for (int j_off=0; j_off<2; ++j_off) for (int i_off=0; i_off<2; ++i_off) for (int l_off=0; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off] * sf[2+j_off] * sf[4+k_off] * sf[6+l_off];
	 value += array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] * shapeFactor;
      }
      return value;
   }

   template<typename REAL> inline
   REAL interpolateScalarTSC_4D(const REAL* RESTRICT array,const int32_t* RESTRICT const indices,const REAL* RESTRICT const sf,int SIZE) {
      // Create aliases for indices:
      const int32_t i = indices[0];
      const int32_t j = indices[1];
      const int32_t k = indices[2];
      const int32_t l = indices[3];

      #ifndef NDEBUG
      bool ok=true;
      if (indices[0] < 1 || indices[0] > block::WIDTH_X) ok = false;
      if (indices[1] < 1 || indices[1] > block::WIDTH_Y) ok = false;
      if (indices[2] < 1 || indices[2] > block::WIDTH_Z) ok = false;
      if (indices[3] < 1 || indices[3] > SIZE-2) ok = false;
      if (ok == false) {
	 std::cerr << "(SEP ACCUM STRETCHED) ERROR: Invalid indices in interpolateScalarTSC_4D" << std::endl;
	 std::cerr << "\t" << indices[0] << ' ' << indices[1] << ' ' << indices[2] << ' ' << indices[3] << std::endl;
	 exit(1);
      }
      #endif
      
      REAL value = 0.0;
      for (int k_off=-1; k_off<2; ++k_off) for (int j_off=-1; j_off<2; ++j_off) for (int i_off=-1; i_off<2; ++i_off) for (int l_off=-1; l_off<2; ++l_off) {
	 const Real shapeFactor = sf[0+i_off+1] * sf[3+j_off+1] * sf[6+k_off+1] * sf[9+l_off+1];
	 value += array[block::arrayIndex(i+i_off,j+j_off,k+k_off)*SIZE+l+l_off] * shapeFactor;
      }

      return value;
   }

   template<typename REAL> inline
   REAL interpolateScalarLogicalNGP_4D_DerivL(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      return 0.0;
   }
   
   template<typename REAL> inline
   REAL interpolateScalarLogicalCIC_4D_DerivL(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const REAL L = pos[3] - 0.5;
      
      int32_t indices[4];
      indices[0] = static_cast<int32_t>(X);
      indices[1] = static_cast<int32_t>(Y);
      indices[2] = static_cast<int32_t>(Z);
      indices[3] = static_cast<int32_t>(L);

      const REAL W_X = X - indices[0];
      const REAL W_Y = Y - indices[1];
      const REAL W_Z = Z - indices[2];
      const REAL W_L = L - indices[3];
      
      Real shapeFactors[8];
      shapeFactors[0] = 1-W_X;
      shapeFactors[1] = W_X;
      shapeFactors[2] = 1-W_Y;
      shapeFactors[3] = W_Y;
      shapeFactors[4] = 1-W_Z;
      shapeFactors[5] = W_Z;
      shapeFactors[6] = -1.0;
      shapeFactors[7] = 1.0;
      
      return interpolateScalarCIC_4D(array,indices,shapeFactors,SIZE);
   }
   
   template<typename REAL> inline
   REAL interpolateScalarLogicalTSC_4D_DerivL(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      int32_t indices[4];
      indices[0] = static_cast<int32_t>(pos[0]);
      indices[1] = static_cast<int32_t>(pos[1]);
      indices[2] = static_cast<int32_t>(pos[2]);
      indices[3] = static_cast<int32_t>(pos[3]);
      
      Real shapeFactors[12];
      // X Shape Factors:
      shapeFactors[0]  = 0.5*(indices[0]+1-pos[0])*(indices[0]+1-pos[0]);
      shapeFactors[2]  = 0.5*(pos[0] -indices[0] )*(pos[0] -indices[0] );
      shapeFactors[1]  = 1-shapeFactors[0]-shapeFactors[2];
      // Y Shape Factors:
      shapeFactors[3]  = 0.5*(indices[1]+1-pos[1])*(indices[1]+1-pos[1]);
      shapeFactors[5]  = 0.5*(pos[1] - indices[1] )*(pos[1] - indices[1] );
      shapeFactors[4]  = 1-shapeFactors[3]-shapeFactors[5];
      // Z Shape Factors:
      shapeFactors[6]  = 0.5*(indices[2]+1-pos[2])*(indices[2]+1-pos[2]);
      shapeFactors[8]  = 0.5*(pos[2] -indices[2] )*(pos[2] -indices[2] );
      shapeFactors[7]  = 1-shapeFactors[6]-shapeFactors[8];
      // L Shape Factors:
      shapeFactors[9]  = pos[3] - indices[3]-1;
      shapeFactors[11] = pos[3] - indices[3];
      shapeFactors[10] = -shapeFactors[9]-shapeFactors[11];

      return interpolateScalarTSC_4D(array,indices,shapeFactors,SIZE);
   }

   template<typename REAL> inline
   REAL interpolateScalarLogicalNGP_4D(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);
      const int32_t l = static_cast<int32_t>(pos[3]);
    
      return array[block::arrayIndex(i  ,j  ,k  )*SIZE+l  ];
   }
   
   template<typename REAL> inline
   REAL interpolateScalarLogicalCIC_4D(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      const REAL X = pos[0] - 0.5;
      const REAL Y = pos[1] - 0.5;
      const REAL Z = pos[2] - 0.5;
      const REAL L = pos[3] - 0.5;
      
      const int32_t i = static_cast<int32_t>(X);
      const int32_t j = static_cast<int32_t>(Y);
      const int32_t k = static_cast<int32_t>(Z);
      const int32_t l = static_cast<int32_t>(L);
    
      const REAL W_X = X - i;
      const REAL W_Y = Y - j;
      const REAL W_Z = Z - k;
      const REAL W_L = L - l;
      
      REAL value = 0.0;
      value += array[block::arrayIndex(i  ,j  ,k  )*SIZE+l  ] * (1-W_X)*(1-W_Y)*(1-W_Z)*(1-W_L);
      value += array[block::arrayIndex(i  ,j  ,k  )*SIZE+l+1] * (1-W_X)*(1-W_Y)*(1-W_Z)*  W_L  ;
      value += array[block::arrayIndex(i+1,j  ,k  )*SIZE+l  ] *   W_X  *(1-W_Y)*(1-W_Z)*(1-W_L);
      value += array[block::arrayIndex(i+1,j  ,k  )*SIZE+l+1] *   W_X  *(1-W_Y)*(1-W_Z)*  W_L  ;
      value += array[block::arrayIndex(i  ,j+1,k  )*SIZE+l  ] * (1-W_X)*  W_Y  *(1-W_Z)*(1-W_L);
      value += array[block::arrayIndex(i  ,j+1,k  )*SIZE+l+1] * (1-W_X)*  W_Y  *(1-W_Z)*  W_L  ;
      value += array[block::arrayIndex(i+1,j+1,k  )*SIZE+l  ] *   W_X  *  W_Y  *(1-W_Z)*(1-W_L);
      value += array[block::arrayIndex(i+1,j+1,k  )*SIZE+l+1] *   W_X  *  W_Y  *(1-W_Z)*  W_L  ;
      
      value += array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l  ] * (1-W_X)*(1-W_Y)*  W_Z  *(1-W_L);
      value += array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l+1] * (1-W_X)*(1-W_Y)*  W_Z  *  W_L  ;
      value += array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l  ] *   W_X  *(1-W_Y)*  W_Z  *(1-W_L);
      value += array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l+1] *   W_X  *(1-W_Y)*  W_Z  *  W_L  ;
      value += array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l  ] * (1-W_X)*  W_Y  *  W_Z  *(1-W_L);
      value += array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l+1] * (1-W_X)*  W_Y  *  W_Z  *  W_L  ;
      value += array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l  ] *   W_X  *  W_Y  *  W_Z  *(1-W_L);
      value += array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l+1] *   W_X  *  W_Y  *  W_Z  *  W_L  ;
      
      return value;
   }
      
   template<typename REAL> inline
   REAL interpolateScalarLogicalTSC_4D(const REAL* RESTRICT array,const REAL* RESTRICT const pos,int SIZE) {
      int32_t indices[4];
      indices[0] = static_cast<int32_t>(pos[0]);
      indices[1] = static_cast<int32_t>(pos[1]);
      indices[2] = static_cast<int32_t>(pos[2]);
      indices[3] = static_cast<int32_t>(pos[3]);
      
      REAL shapeFactors[12];
      shapeFactors[0]  = 0.5*(indices[0]+1-pos[0])*(indices[0]+1-pos[0]);
      shapeFactors[2]  = 0.5*(pos[0] -indices[0] )*(pos[0] -indices[0] );
      shapeFactors[1]  = 1-shapeFactors[0]-shapeFactors[2];
      
      shapeFactors[3]  = 0.5*(indices[1]+1-pos[1])*(indices[1]+1-pos[1]);
      shapeFactors[5]  = 0.5*(pos[1] - indices[1] )*(pos[1] - indices[1] );
      shapeFactors[4]  = 1-shapeFactors[3]-shapeFactors[5];
      
      shapeFactors[6]  = 0.5*(indices[2]+1-pos[2])*(indices[2]+1-pos[2]);
      shapeFactors[8]  = 0.5*(pos[2] -indices[2] )*(pos[2] -indices[2] );
      shapeFactors[7]  = 1-shapeFactors[6]-shapeFactors[8];
      
      shapeFactors[9]  = 0.5*(indices[3]+1-pos[3])*(indices[3]+1-pos[3]);
      shapeFactors[11] = 0.5*(pos[3] -indices[3] )*(pos[3] -indices[3] );
      shapeFactors[10] = 1-shapeFactors[9]-shapeFactors[11];
      
      return interpolateScalarTSC_4D(array,indices,shapeFactors,SIZE);

      const int32_t i = static_cast<int32_t>(pos[0]);
      const int32_t j = static_cast<int32_t>(pos[1]);
      const int32_t k = static_cast<int32_t>(pos[2]);
      const int32_t l = static_cast<int32_t>(pos[3]);
      
      const REAL W_X_N = 0.5*(i+1-pos[0])*(i+1-pos[0]);
      const REAL W_X_P = 0.5*(pos[0] -i )*(pos[0] -i );
      const REAL W_X_C = 1-W_X_N-W_X_P;
      
      const REAL W_Y_N = 0.5*(j+1-pos[1])*(j+1-pos[1]);
      const REAL W_Y_P = 0.5*(pos[1] -j )*(pos[1] -j );
      const REAL W_Y_C = 1-W_Y_N-W_Y_P;
      
      const REAL W_Z_N = 0.5*(k+1-pos[2])*(k+1-pos[2]);
      const REAL W_Z_P = 0.5*(pos[2] -k )*(pos[2] -k );
      const REAL W_Z_C = 1-W_Z_N-W_Z_P;
      
      const REAL W_L_N = 0.5*(l+1-pos[3])*(l+1-pos[3]);
      const REAL W_L_P = 0.5*(pos[3] -l )*(pos[3] -l );
      const REAL W_L_C = 1-W_L_N-W_L_P;
      
      REAL value = 0.0;

      // K - 1, J - 1
      value += array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l-1] * W_X_N * W_Y_N * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l  ] * W_X_N * W_Y_N * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i-1,j-1,k-1)*SIZE+l+1] * W_X_N * W_Y_N * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l-1] * W_X_C * W_Y_N * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l  ] * W_X_C * W_Y_N * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i  ,j-1,k-1)*SIZE+l+1] * W_X_C * W_Y_N * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l-1] * W_X_P * W_Y_N * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l  ] * W_X_P * W_Y_N * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i+1,j-1,k-1)*SIZE+l+1] * W_X_P * W_Y_N * W_Z_N * W_L_P;

      // J + 0
      value += array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l-1] * W_X_N * W_Y_C * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l  ] * W_X_N * W_Y_C * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i-1,j  ,k-1)*SIZE+l+1] * W_X_N * W_Y_C * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l-1] * W_X_C * W_Y_C * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l  ] * W_X_C * W_Y_C * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i  ,j  ,k-1)*SIZE+l+1] * W_X_C * W_Y_C * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l-1] * W_X_P * W_Y_C * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l  ] * W_X_P * W_Y_C * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i+1,j  ,k-1)*SIZE+l+1] * W_X_P * W_Y_C * W_Z_N * W_L_P;

      // J + 1
      value += array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l-1] * W_X_N * W_Y_P * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l  ] * W_X_N * W_Y_P * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i-1,j+1,k-1)*SIZE+l+1] * W_X_N * W_Y_P * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l-1] * W_X_C * W_Y_P * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l  ] * W_X_C * W_Y_P * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i  ,j+1,k-1)*SIZE+l+1] * W_X_C * W_Y_P * W_Z_N * W_L_P;
      
      value += array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l-1] * W_X_P * W_Y_P * W_Z_N * W_L_N;
      value += array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l  ] * W_X_P * W_Y_P * W_Z_N * W_L_C;
      value += array[block::arrayIndex(i+1,j+1,k-1)*SIZE+l+1] * W_X_P * W_Y_P * W_Z_N * W_L_P;
      
      // K + 0, J - 1
      value += array[block::arrayIndex(i-1,j-1,k  )*SIZE+l-1] * W_X_N * W_Y_N * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i-1,j-1,k  )*SIZE+l  ] * W_X_N * W_Y_N * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i-1,j-1,k  )*SIZE+l+1] * W_X_N * W_Y_N * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i  ,j-1,k  )*SIZE+l-1] * W_X_C * W_Y_N * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i  ,j-1,k  )*SIZE+l  ] * W_X_C * W_Y_N * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i  ,j-1,k  )*SIZE+l+1] * W_X_C * W_Y_N * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i+1,j-1,k  )*SIZE+l-1] * W_X_P * W_Y_N * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i+1,j-1,k  )*SIZE+l  ] * W_X_P * W_Y_N * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i+1,j-1,k  )*SIZE+l+1] * W_X_P * W_Y_N * W_Z_C * W_L_P;
      
      // J + 0
      value += array[block::arrayIndex(i-1,j  ,k  )*SIZE+l-1] * W_X_N * W_Y_C * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i-1,j  ,k  )*SIZE+l  ] * W_X_N * W_Y_C * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i-1,j  ,k  )*SIZE+l+1] * W_X_N * W_Y_C * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i  ,j  ,k  )*SIZE+l-1] * W_X_C * W_Y_C * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i  ,j  ,k  )*SIZE+l  ] * W_X_C * W_Y_C * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i  ,j  ,k  )*SIZE+l+1] * W_X_C * W_Y_C * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i+1,j  ,k  )*SIZE+l-1] * W_X_P * W_Y_C * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i+1,j  ,k  )*SIZE+l  ] * W_X_P * W_Y_C * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i+1,j  ,k  )*SIZE+l+1] * W_X_P * W_Y_C * W_Z_C * W_L_P;
      
      // J + 1
      value += array[block::arrayIndex(i-1,j+1,k  )*SIZE+l-1] * W_X_N * W_Y_P * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i-1,j+1,k  )*SIZE+l  ] * W_X_N * W_Y_P * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i-1,j+1,k  )*SIZE+l+1] * W_X_N * W_Y_P * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i  ,j+1,k  )*SIZE+l-1] * W_X_C * W_Y_P * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i  ,j+1,k  )*SIZE+l  ] * W_X_C * W_Y_P * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i  ,j+1,k  )*SIZE+l+1] * W_X_C * W_Y_P * W_Z_C * W_L_P;
      
      value += array[block::arrayIndex(i+1,j+1,k  )*SIZE+l-1] * W_X_P * W_Y_P * W_Z_C * W_L_N;
      value += array[block::arrayIndex(i+1,j+1,k  )*SIZE+l  ] * W_X_P * W_Y_P * W_Z_C * W_L_C;
      value += array[block::arrayIndex(i+1,j+1,k  )*SIZE+l+1] * W_X_P * W_Y_P * W_Z_C * W_L_P;
      
      // K + 1, J - 1
      value += array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l-1] * W_X_N * W_Y_N * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l  ] * W_X_N * W_Y_N * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i-1,j-1,k+1)*SIZE+l+1] * W_X_N * W_Y_N * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l-1] * W_X_C * W_Y_N * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l  ] * W_X_C * W_Y_N * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i  ,j-1,k+1)*SIZE+l+1] * W_X_C * W_Y_N * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l-1] * W_X_P * W_Y_N * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l  ] * W_X_P * W_Y_N * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i+1,j-1,k+1)*SIZE+l+1] * W_X_P * W_Y_N * W_Z_P * W_L_P;
      
      // J + 0
      value += array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l-1] * W_X_N * W_Y_C * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l  ] * W_X_N * W_Y_C * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i-1,j  ,k+1)*SIZE+l+1] * W_X_N * W_Y_C * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l-1] * W_X_C * W_Y_C * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l  ] * W_X_C * W_Y_C * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i  ,j  ,k+1)*SIZE+l+1] * W_X_C * W_Y_C * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l-1] * W_X_P * W_Y_C * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l  ] * W_X_P * W_Y_C * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i+1,j  ,k+1)*SIZE+l+1] * W_X_P * W_Y_C * W_Z_P * W_L_P;
      
      // J + 1
      value += array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l-1] * W_X_N * W_Y_P * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l  ] * W_X_N * W_Y_P * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i-1,j+1,k+1)*SIZE+l+1] * W_X_N * W_Y_P * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l-1] * W_X_C * W_Y_P * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l  ] * W_X_C * W_Y_P * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i  ,j+1,k+1)*SIZE+l+1] * W_X_C * W_Y_P * W_Z_P * W_L_P;
      
      value += array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l-1] * W_X_P * W_Y_P * W_Z_P * W_L_N;
      value += array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l  ] * W_X_P * W_Y_P * W_Z_P * W_L_C;
      value += array[block::arrayIndex(i+1,j+1,k+1)*SIZE+l+1] * W_X_P * W_Y_P * W_Z_P * W_L_P;
      
      return value;
   }
   
   template<typename REAL> inline
   void loadScalar4D(SimulationClasses* simClasses,pargrid::CellID blockLID,const REAL* RESTRICT source,REAL* RESTRICT dest,uint32_t SIZE) {
      const size_t SRC_SIZE = block::SIZE*SIZE;
      const size_t DEST_SIZE = SIZE+2;
      
      pargrid::CellID nbrLID = simClasses->pargrid.invalidCellID();
      const pargrid::CellID* const RESTRICT nbrs = simClasses->pargrid.getCellNeighbourIDs(blockLID);
      const Real* srcTmp = NULL;

      // ***** LOAD VALUES FROM THIS BLOCK ***** //
      srcTmp = source + blockLID*SRC_SIZE;
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	dest[block::arrayIndex(i+1,j+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,j,k)*SIZE+l];

      // ***** LOAD VALUES FROM 6 FACE NEIGHBOURS ***** //
      // -x neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,j+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,j,k)*SIZE+l];
      }
      // +x neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,j,k)*SIZE+l];
      }
      // -y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,0,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,block::WIDTH_Y-1,k)*SIZE+l];
      }
      // +y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,0,k)*SIZE+l];
      }
      // -z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,j+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(i,j,block::WIDTH_Z-1)*SIZE+l];
      }
      // +z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,j,0)*SIZE+l];
      }

      // ***** LOAD VALUES FROM 12 EDGE NEIGHBOURS ***** //
      // -x,-y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,0,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k)*SIZE+l];
      }
      // +x,-y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,0,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,block::WIDTH_Y-1,k)*SIZE+l]; 
      }
      // -x,+y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,block::WIDTH_Y+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,0,k)*SIZE+l]; 
      }
      // +x,+y neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,+0)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int k=0; k<block::WIDTH_Z; ++k) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,0,k)*SIZE+l];
      }
      // -x,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,j+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1)*SIZE+l]; 
      }
      // +x,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(0,j,block::WIDTH_Z-1)*SIZE+l];
      }
      // -x,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,j+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,j,0)*SIZE+l];
      }
      // +x,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+0,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int j=0; j<block::WIDTH_Y; ++j) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,j,0)*SIZE+l];
      }
      // -y,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,0,0)*DEST_SIZE+l+1] = srcTmp[block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l];
      }
      // +y,-z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(i,0,block::WIDTH_Z-1)*SIZE+l];
      }
      // -y,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,0,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,block::WIDTH_Y-1,0)*SIZE+l];
      }
      // +y,+z neighbour:
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+0,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (int i=0; i<block::WIDTH_X; ++i) for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(i,0,0)*SIZE+l];
      }

      // ***** LOAD VALUES FROM 8 CORNER NEIGHBOURS ***** //
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,0,0)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,0,0)*DEST_SIZE+l+1] = srcTmp[block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,block::WIDTH_Y+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,-1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)*DEST_SIZE+l+1] = srcTmp[block::index(0,0,block::WIDTH_Z-1)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,0,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,-1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,block::WIDTH_Y-1,0)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(-1,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(block::WIDTH_X-1,0,0)*SIZE+l];
      }
      
      nbrLID = nbrs[simClasses->pargrid.calcNeighbourTypeID(+1,+1,+1)];
      if (nbrLID != simClasses->pargrid.invalidCellID()) {
	 srcTmp = source + nbrLID*SRC_SIZE;
	 for (uint32_t l=0; l<SIZE; ++l)
	   dest[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*DEST_SIZE+l+1] = srcTmp[block::index(0,0,0)*SIZE+l];
      }

   }
   
}

#endif
