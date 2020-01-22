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

#include "example_advection.h"

using namespace std;

// Init static variables:
pargrid::StencilID Advection::faceNbrStencilID;
pargrid::DataID Advection::dataAvgsID;
pargrid::DataID Advection::dataFluxID;
std::string Advection::avgsName;
std::string Advection::fluxName;
Real Advection::Vx;
Real Advection::Vy;
Real Advection::Vz;

