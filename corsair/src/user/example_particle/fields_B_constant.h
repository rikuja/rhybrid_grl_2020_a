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

#ifndef FIELDS_B_CONSTANT_H
#define FIELDS_B_CONSTANT_H

#include <cmath>
#include "configreader.h"
#include "simulation.h"
#include "simulationclasses.h"
#include "linear_algebra.h"
#include "particle_definition.h"

class ConstantB {
 public:
   ConstantB();
   ~ConstantB();
   
   static bool finalize();
   template<typename REAL,class SPECIES> static void getAcceleration(pargrid::CellID cellID,const REAL& t,const REAL& dt,const SPECIES& species,const REAL* state,REAL* acc);
   template<typename REAL> static void getFields(pargrid::CellID cellID,const REAL& t,const REAL* state,REAL E[3],REAL B[3]);
   template<typename REAL> static void getFields(pargrid::CellID cellID,const REAL& t,const REAL* state,REAL E[3],REAL B[3],REAL dB[9]);
   template<typename REAL,class SPECIES> static REAL getMaximumTimestep(pargrid::CellID cellID,REAL t,REAL dt_suggested,const SPECIES& species,const REAL* state);
   static bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   
 private:
     
   static Real B[3];       /**< Value of constant magnetic field.*/
};

template<typename REAL,class SPECIES>
void ConstantB::getAcceleration(pargrid::CellID cellID,const REAL& t,const REAL& dt,const SPECIES& species,const REAL* state,REAL* acc) {
   const REAL B_mag = vectorMagnitude<3>(B);
   acc[particle::XPOS] = state[particle::VPAR]*B[particle::XPOS]/B_mag;
   acc[particle::YPOS] = state[particle::VPAR]*B[particle::YPOS]/B_mag;
   acc[particle::ZPOS] = state[particle::VPAR]*B[particle::ZPOS]/B_mag;
   acc[particle::VPAR] = 0.0;
}

template<typename REAL> 
void ConstantB::getFields(pargrid::CellID cellID,const REAL& t,const REAL* state,REAL E[3],REAL B[3]) {
   for (int i=0; i<3; ++i) E[i] = 0.0;
   for (int i=0; i<3; ++i) B[i] = ConstantB::B[i];
}

template<typename REAL>
void ConstantB::getFields(pargrid::CellID cellID,const REAL& t,const REAL* state,REAL E[3],REAL B[3],REAL dB[9]) {
   for (int i=0; i<3; ++i) E[i] = 0.0;
   for (int i=0; i<3; ++i) B[i] = ConstantB::B[i];
   for (int i=0; i<9; ++i) dB[i] = 0.0;
}

template<typename REAL,class SPECIES>
REAL ConstantB::getMaximumTimestep(pargrid::CellID cellID,REAL t,REAL dt_suggested,const SPECIES& species,const REAL* state) {
   return dt_suggested;
}

#endif
