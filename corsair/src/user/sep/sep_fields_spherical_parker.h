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

#ifndef SEP_FIELDS_SPHERICAL_PARKER_H
#define SEP_FIELDS_SPHERICAL_PARKER_H

#include <configreader.h>
#include <simulation.h>
#include <simulationclasses.h>

namespace sep {

   struct SphericalParker {
      bool initialized;
      const double* blockCoordinates;
      Real B_earth_AU2;               /**< (Signed) magnitude of radial B at 1 AU, multiplied by AU squared.*/
      Real V_sw_r;                    /**< Magnitude of radial solar wind velocity.*/   
      Real rho_mass_AU2;              /**< Solar wind mass density at 1 AU, multiplied by AU squared.*/
      Simulation* sim;
      SimulationClasses* simClasses;
      
      SphericalParker();
   };
   
   bool sphericalParkerFieldFinalize();
   void sphericalParkerFieldGetFields(pargrid::CellID cellID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]);
   void sphericalParkerFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* position,PlasmaState& plasmaState);
   void sphericalParkerFieldGetState(pargrid::CellID cellID,Real t,const Real* position,Real* B,Real* V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign);
   bool sphericalParkerFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   
} // namespace sep
   
#endif
