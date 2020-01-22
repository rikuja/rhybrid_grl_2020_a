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

#include <cmath>
#include "configreader.h"
#include "fields_B_constant.h"
#include "linear_algebra.h"

using namespace std;

Real ConstantB::B[3];

ConstantB::ConstantB() { }

ConstantB::~ConstantB() { }

bool ConstantB::finalize() {return true;}

bool ConstantB::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
   bool success = true;
   
   // Define parameters read from configuration file(s) and parse:
   const string configName = "ConstantB";
   const Real DEFVALUE = NAN;
   Real B_mag = NAN;
   cr.add(configName+".direction_x","Vector to direction of B, x-component (float).",DEFVALUE);
   cr.add(configName+".direction_y","Vector to direction of B, y-component (float).",DEFVALUE);
   cr.add(configName+".direction_z","Vector to direction of B, z-component (float).",DEFVALUE);
   cr.add(configName+".magnitude","Magnitude of B in nT (float).",DEFVALUE);
   cr.parse();
   cr.get(configName+".direction_x",B[0]);
   cr.get(configName+".direction_y",B[1]);
   cr.get(configName+".direction_z",B[2]);
   cr.get(configName+".magnitude",B_mag);
   
   // Check input values for sanity:
   if (B[0] != B[0]) {
      simClasses.logger << "(CONSTANT_B) ERROR: x-component of vector to direction of B was not given with parameter '" << configName << ".direction_x' !" << endl << write;
      success = false;
   }
   if (B[1] != B[1]) {
      simClasses.logger << "(CONSTANT_B) ERROR: y-component of vector to direction of B was not given with parameter '" << configName << ".direction_y' !" << endl << write;
      success = false;
   }
   if (B[2] != B[2]) {
      simClasses.logger << "(CONSTANT_B) ERROR: z-component of vector to direction of B was not given with parameter '" << configName << ".direction_z' !" << endl << write;
      success = false;
   }
   if (B_mag != B_mag) {
      simClasses.logger << "(CONSTANT_B) ERROR: Magnitude of B was not given with parameter '" << configName << ".magnitude' !" << endl << write;
      success = false;
   }
   
   const Real magnitude = vectorMagnitude<3>(B);
   for (int i=0; i<3; ++i) B[i] = B[i]*B_mag/magnitude;
   
   return success;
}
