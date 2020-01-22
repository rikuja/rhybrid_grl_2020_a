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

#ifndef SEP_FIELDS_CONTAINER_H
#define SEP_FIELDS_CONTAINER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_definitions.h>

namespace sep {

   struct PlasmaState {
      Real B[3];
      Real E[3];
      Real alfvenMachNumber;                 /**< Alfvenic Mach number, uninitialized in constructor.*/
      Real alfvenSpeed2;                     /**< Alfven speed squared, uninitialized in constructor.*/
      Real electronNumberDensity;
      Real electronTemperature;
      Real ionMassDensity;
      Real ionPolytropicIndex;
      Real ionTemperature;                   /**< Ion temperature.*/
      Real soundSpeed2;                      /**< Sound speed squared.*/
      Real V_plasma_SIM[3];
      
      PlasmaState();
   };
   
   typedef bool (*finalizeField)();
   typedef void (*getFields)(pargrid::CellID cellID,Real t,const Real* position,Real E[3],Real B[3],Real gradB[9]);
   typedef void (*getPlasmaStateFunction)(pargrid::CellID cellID,Real t,const Real* position,PlasmaState& plasmaState);
   typedef void (*getStateFunction)(pargrid::CellID cellID,Real t,const Real* position,Real* B,Real* V_wave,Real& dV_wave,Real& V_alfven,int alfvenSign);
   typedef bool (*initializeField)(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
   
   struct FieldFunction {
      finalizeField finalize;
      getFields getFieldsFunction;
      getPlasmaStateFunction getPlasmaState;
      getStateFunction getState;
      initializeField initialize;
      
      FieldFunction();
      FieldFunction(finalizeField finalize,getFields getFieldsFunction,getPlasmaStateFunction getPlasma,
		    getStateFunction getState,initializeField initialize);
   };
   
   class FieldsContainer {
    public:
      FieldsContainer();
      ~FieldsContainer();
      bool getField(const std::string& name,finalizeField& finalize,getFields& getFieldsFunction,
		    getPlasmaStateFunction& getPlasma,getStateFunction& getState,initializeField& initialize);
      bool registerField(const std::string& name,finalizeField finalize,getFields getFieldsFunction,
			 getPlasmaStateFunction getPlasma,getStateFunction getState,initializeField initialize);
      
    private:
      std::map<std::string,FieldFunction> registeredFields;
   };

} // namespace sep
   
#endif
