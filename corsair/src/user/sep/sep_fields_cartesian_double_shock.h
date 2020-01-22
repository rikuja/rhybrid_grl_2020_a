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

#ifndef SEP_FIELDS_CARTESIAN_DOUBLE_SHOCK_H
#define SEP_FIELDS_CARTESIAN_DOUBLE_SHOCK_H

#include <configreader.h>
#include <simulation.h>
#include <simulationclasses.h>

namespace sep {

   struct CartesianDoubleShock {
      bool initialized;                      /**< If true, all variables have been initialized correctly.*/
      Real ionNumberDensity;         /**< Undisturbed plasma ion number density.*/
      Real ionPolytropicIndex;       /**< Ion polytropic index.*/
      Real massIon;                  /**< Mass of plasma ion species.*/
      Simulation* sim;
      SimulationClasses* simClasses;

      Real ionPressure1;                     /**< Plasma pressure upstream of leading shock.*/
      Real ionPressure2;                     /**< Plasma pressure in region between shocks.*/
      Real ionPressure3;                     /**< Plasma pressure downstream of trailing shock.*/
      Real ionTemperature1;                  /**< Ion temperature upstream of leading shock.*/
      Real ionTemperature2;                  /**< Ion temperature in region between shocks.*/
      Real ionTemperature3;                  /**< Ion temperature downstream of trailing shock.*/

      Real leadingShockCentroid[3];          /**< Leading shock centroid at current simulation time.*/
      Real leadingShockGasCompressionRatio;
      Real leadingShockMagnCompressionRatio;
      Real leadingShockInitialCentroid[3];   /**< Leading shock centroid position at t=0.*/
      Real leadingShockNormal[3];            /**< Leading shock normal.*/
      Real leadingShockVelocity[3];          /**< Leading shock velocity vector.*/

      Real trailingShockCentroid[3];         /**< Trailing shock centroid at current simulation time.*/
      Real trailingShockGasCompressionRatio; 
      Real trailingShockMagnCompressionRatio;
      Real trailingShockInitialCentroid[3];  /**< Trailing shock centroid position at t=0.*/
      Real trailingShockNormal[3];           /**< Trailing shock normal.*/
      Real trailingShockVelocity[3];         /**< Trailing shock velocity vector.*/
      Real relativeShockCentroid[3];         /**< Trailing shock's position wrt leading shock.*/
      Real relativeShockVelocity[3];         /**< Trailing shock's velocity relative to leading shock.*/

      Real doubleGasCompressionRatio;        /**< Product of leading and trailing shock gas compression ratios.*/
      Real doubleMagnCompressionRatio;       /**< Product of leading and trailing shock magnetic compression ratios.*/

      Real B1[3];                            /**< Magnetic field in simulation frame in upstream region of leading shock.*/
      Real B2[3];                            /**< Magnetic field in simulation frame in region between leading and trailing shocks.*/
      Real B3[3];                            /**< Magnetic field in simulation frame in downstream region of trailing shock.*/

      Real E1_SIM[3];                        /**< Electric field in simulation frame in upstream region of leading shock.*/
      Real E2_SIM[3];                        /**< Electric field in simulation frame in region between leading and trailing shocks.*/
      Real E3_SIM[3];                        /**< Electric field in simulation frame in dowstream region of trailing shock.*/

      Real V1_alfven[3];
      Real V2_alfven[3];
      Real V3_alfven[3];
      Real V_alfven_surface;
      
      Real V1_plasma_SIM[3];                 /**< Plasma velocity in simulation frame in upstream region of leading shock.*/
      Real V2_plasma_SIM[3];                 /**< Plasma velocity in simulation frame in region between leading and trailing shocks.*/
      Real V3_plasma_SIM[3];                 /**< Plasma velocity in simulation frame in dowstream region of trailing shock.*/

      Real referenceDistance;
      
      CartesianDoubleShock();
   };
   
   bool cartesianDoubleShockFieldFinalize();
   void cartesianDoubleShockFieldGetFields(pargrid::CellID cellID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]);
   void cartesianDoubleShockFieldGetPlasmaState(pargrid::CellID cellID,Real t,const Real* position,PlasmaState& plasmaState);
   void cartesianDoubleShockFieldGetState(pargrid::CellID cellID,Real t,const Real* position,Real* B,Real* V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign);
   bool cartesianDoubleShockFieldInitialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

} // namespace sep
   
#endif
