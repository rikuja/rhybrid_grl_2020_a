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

#ifndef SEP_PARTICLE_DEFINITION_H
#define SEP_PARTICLE_DEFINITION_H

#include <mpi.h>
#include <mpiconversion.h>

// First three data fields are reserved for coordinates:

namespace sep {

   namespace particle {
      enum STATE {
	 XCRD,             /**< Guiding center x coordinate.*/
	 YCRD,             /**< Guiding center y coordinate.*/
	 ZCRD,             /**< Guiding center z coordinate.*/
	 V_PAR,            /**< Guiding center parallel speed.*/
	 MU,               /**< Magnetic moment.*/
	 WEIGHT,           /**< Number of real particles represented by this macroparticle.*/
	 DATASIZE
      };
   }

   template<typename REAL>
   struct Particle {
      REAL state[particle::DATASIZE];
      
      static void getDatatype(MPI_Datatype& datatype);
   };

   template<typename REAL> inline
   void Particle<REAL>::getDatatype(MPI_Datatype& datatype) {
      MPI_Type_contiguous(particle::DATASIZE,MPI_Type<REAL>(),&datatype);
   }
   
} // namespace sep

#endif
