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

#ifndef SEP_LAGR_DEFINITION_H
#define SEP_LAGR_DEFINITION_H

#include <mpi.h>
#include <mpiconversion.h>

// First three data fields are reserved for coordinates:

namespace sep {

   namespace lagr {
      enum STATE {XCRD,YCRD,ZCRD,LAMBDA,ENERGY,DATASIZE};
   }

   template<typename REAL>
   struct LagrangianParticle {
      REAL state[lagr::DATASIZE];
      
      static void getDatatype(MPI_Datatype& datatype);
   };

   template<typename REAL> inline
   void LagrangianParticle<REAL>::getDatatype(MPI_Datatype& datatype) {
      MPI_Type_contiguous(lagr::DATASIZE,MPI_Type<REAL>(),&datatype);
   }
   
} // namespace sep

#endif
