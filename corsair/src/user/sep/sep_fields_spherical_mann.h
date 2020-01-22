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

#ifndef SEP_FIELDS_SPHERICAL_MANN_H
#define SEP_FIELDS_SPHERICAL_MANN_H

#include <configreader.h>
#include <simulation.h>
#include <simulationclasses.h>

#include "sep_fields_container.h"

namespace sep {

   struct SphericalMann {
      bool initialized;
      const double* blockCoordinates;

      int resolution;
      int N_points;                   /**< Size of arrays V_alfven, V_radial.*/
      Real* V_alfven;                 /**< Alfven speed.*/
      Real* V_radial;                 /**< Solar wind speed.*/
      Real r_min;                     /**< Minimum value of radial coordinate, used in interpolations.*/
      Real r_max;                     /**< Maximum value of radial coordinate, used in interpolations.*/
      Real dr;                        /**< Cell size in radial direction.*/
      
      // ***** INPUT PARAMETERS ***** //
      Real numberDensity;             /**< Solar wind number density at reference radius.*/
      Real radialMagneticField;       /**< Solar wind radial magnetic field at reference radius.*/
      Real radialSpeed;               /**< Solar wind radial velocity at reference radius.*/
      Real referenceRadius;           /**< Radius at which parameters were defined, defaults to 1 AU.*/
      Real temperature;               /**< Solar wind temperature (isothermal).*/
      Real meanMolecularWeight;       /**< Mean molecular weight.*/
      Real ionPolytropicIndex;        /**< Ion polytropic index.*/
      
      // ***** DERIVED PARAMETERS ***** //
      Real mannModelConstant;
      Real criticalRadius;
      Real criticalSpeed;
      Real B_radial_R2;
      
      Simulation* sim;
      SimulationClasses* simClasses;
      
      SphericalMann();
      ~SphericalMann();
      void interpolate(Real radius,Real& V_sw,Real& V_alfven,Real& dV_alfven,int alfvenSign);
   };

   bool sphericalMannFieldFinalize();
   void sphericalMannFieldGetFields(pargrid::CellID cellID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]);
   void sphericalMannFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* position,PlasmaState& plasmaState);
   void sphericalMannFieldGetState(pargrid::CellID cellID,Real t,const Real* position,Real* B,Real* V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign);
   bool sphericalMannFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   
} // namespace sep
   
#endif
