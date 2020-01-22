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

#include <climits>

#include "particle_species.h"

using namespace std;

bool Species::finalize() {
   return true;   
}

std::string Species::getName() const {
   return name;
}

const std::string& Species::getSpeciesType() const {
   return speciesType;
}

bool Species::readParameters(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& name) {
   bool success = true;
   this->name = name;
   
   // Read species' parameters from config file:
   string q_unit,m_unit;
   cr.add(name+".mass_unit","Unit in which 'mass' is given (string).",string(""));
   cr.add(name+".charge_unit","Unit which 'charge' is given (string).",string(""));
   cr.add(name+".mass","Mass in mass units (float).",numeric_limits<Real>::infinity());
   cr.add(name+".charge","Charge in charge units (float).",numeric_limits<Real>::infinity());
   cr.add(name+".type","What type of species this is (string)?",string(""));
   cr.parse();
   cr.get(name+".mass_unit",m_unit);
   cr.get(name+".charge_unit",q_unit);
   cr.get(name+".mass",m);
   cr.get(name+".charge",q);
   cr.get(name+".type",speciesType);
   
   // Check input parameters for sanity:
   Real charge = simClasses.constants.get(q_unit);
   if (charge == simClasses.constants.notFound()) {
      simClasses.logger << "(SPECIES) ERROR: illegal charge unit '" << q_unit << "' !" << endl << write;
      success = false;
   }
   Real mass = simClasses.constants.get(m_unit);
   if (mass == simClasses.constants.notFound()) {
      simClasses.logger << "(SPECIES) ERROR: illegal mass unit '" << m_unit << "' !" << endl << write;
      success = false;
   }   
   if (q == numeric_limits<Real>::infinity()) {
      simClasses.logger << "(SPECIES) ERROR: Charge was not specified with parameter '" << name+".charge' !" << endl << write;
      success = false;
   }
   if (m == numeric_limits<Real>::infinity()) {
      simClasses.logger << "(SPECIES) ERROR: Mass was not specified with parameter '" << name+".mass' !" << endl << write;
      success = false;
   }
   if (speciesType == "") {
      simClasses.logger << "(SPECIES) ERROR: Species type was not specified with parameter '" << name+".type' !" << endl << write;
      success = false;
   }
   
   q *= charge;
   m *= mass;
   q_per_m = q/m;
   
   return success;
}

