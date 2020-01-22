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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>

#include "constants.h"

using namespace std;

Constants::Constants() {
   values["BOLTZMANN"] = constants::BOLTZMANN;
   values["CHARGE_ELEMENTARY"] = constants::CHARGE_ELEMENTARY;
   values["DIST_ASTRONOMICAL_UNIT"] = constants::DIST_ASTRONOMICAL_UNIT;
   values["DIST_EARTH_RADIUS"] = constants::DIST_EARTH_RADIUS;
   values["DIST_METER"] = constants::DIST_METER;
   values["DIST_SOLAR_RADIUS"] = constants::DIST_SOLAR_RADIUS;
   values["DIST_VENUS_RADIUS"] = constants::DIST_VENUS_RADIUS;
   values["GRAVITY"] = constants::GRAVITY;
   values["MASS_ELECTRON"] = constants::MASS_ELECTRON;
   values["MASS_PROTON"] = constants::MASS_PROTON;
   values["MASS_SOLAR"] = constants::MASS_SOLAR;
   values["OMEGA_SOLAR_EQUATOR"] = constants::OMEGA_SOLAR_EQUATOR;
   values["PERMEABILITY"] = constants::PERMEABILITY;
   values["PERMITTIVITY"] = constants::PERMITTIVITY;
   values["SPEED_LIGHT"] = constants::SPEED_LIGHT;
}

double Constants::get(const string& name) const {
   map<string,double>::const_iterator it = values.find(name);
   if (it == values.end()) return numeric_limits<double>::infinity();
   return it->second;
}

double Constants::getDistanceInSI(const string& name) const {
   if (name == "m") return 1.0;
   else if (name == "km") return 1000.0;
   else if (name == "AU") return constants::DIST_ASTRONOMICAL_UNIT;
   else if (name == "RE") return constants::DIST_EARTH_RADIUS;
   else if (name == "RS") return constants::DIST_SOLAR_RADIUS;
   else if (name == "RV") return constants::DIST_VENUS_RADIUS;
   else return numeric_limits<double>::infinity();
}

double Constants::getEnergyInSI(const string& name) const {
   if (name == "eV") return constants::CHARGE_ELEMENTARY;
   else if (name == "keV") return 1e3*constants::CHARGE_ELEMENTARY;
   else if (name == "MeV") return 1e6*constants::CHARGE_ELEMENTARY;
   else if (name == "GeV") return 1e9*constants::CHARGE_ELEMENTARY;
   else return numeric_limits<double>::infinity();
}

double Constants::notFound() const {return numeric_limits<double>::infinity();}
