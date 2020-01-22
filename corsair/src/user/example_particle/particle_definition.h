/** This file is part of Corsair simulation.
 *
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

#ifndef PARTICLE_DEFINITION_H
#define PARTICLE_DEFINITION_H

#include <mpi.h>

// First three data fields are reserved for coordinates:

namespace particle {
   enum STATE {XPOS,YPOS,ZPOS,VPAR,MU,WEIGHT,REFLECTED,SIZE};
}

template<typename REAL>
struct Prticle {
   REAL state[particle::SIZE];
   
   REAL& operator[](const int& element);
   static void getDatatype(MPI_Datatype& datatype);
};

template<typename REAL> inline
void Prticle<REAL>::getDatatype(MPI_Datatype& datatype) {
   MPI_Type_contiguous(particle::SIZE,MPI_Type<REAL>(),&datatype);
}

template<typename REAL> inline
REAL& Prticle<REAL>::operator[](const int& element) {return state[element];}

#endif
