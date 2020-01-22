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

#include <cstdlib>
#include <iostream>

#include "sep_fields_container.h"

using namespace std;

namespace sep {

   FieldFunction::FieldFunction() {
      finalize= NULL;
      getFieldsFunction = NULL;
      getPlasmaState = NULL;
      getState = NULL;
      initialize = NULL;
   }
   
   FieldFunction::FieldFunction(finalizeField finalize,getFields getFieldsFunction,getPlasmaStateFunction getPlasma,
				getStateFunction getState,initializeField initialize): 
     finalize(finalize), getFieldsFunction(getFieldsFunction), getPlasmaState(getPlasma), getState(getState), initialize(initialize) {
      
   }
   
   FieldsContainer::FieldsContainer() { }
   
   FieldsContainer::~FieldsContainer() { 
      for (map<string,FieldFunction>::iterator it=registeredFields.begin(); it!=registeredFields.end(); ++it) { 
	 (*it->second.finalize)();
      }
   }
   
   bool FieldsContainer::getField(const std::string& name,finalizeField& finalize,getFields& getFieldsFunction,
				  getPlasmaStateFunction& getPlasmaState,getStateFunction& getState,initializeField& initialize) {
      map<string,FieldFunction>::iterator it = registeredFields.find(name);
      if (it == registeredFields.end()) return false;
      
      finalize          = it->second.finalize;
      getFieldsFunction = it->second.getFieldsFunction;
      getPlasmaState    = it->second.getPlasmaState;
      getState          = it->second.getState;
      initialize        = it->second.initialize;
      return true;
   }
   
   bool FieldsContainer::registerField(const std::string& name,finalizeField finalize,getFields getFieldsFunction,
				       getPlasmaStateFunction getPlasmaState,getStateFunction getState,initializeField initialize) {
      if (finalize == NULL) return false;
      if (getFieldsFunction == NULL) return false;
      if (getPlasmaState == NULL) return false;
      if (getState == NULL) return false;
      if (initialize == NULL) return false;
      if (registeredFields.find(name) != registeredFields.end()) return false;
      registeredFields[name] = FieldFunction(finalize,getFieldsFunction,getPlasmaState,getState,initialize);
      return true;
   }
   
   PlasmaState::PlasmaState() {
      // Default value for all parameters is infinity:
      const Real DEF_VALUE = numeric_limits<Real>::infinity();

      for (int i=0; i<3; ++i) B[i] = DEF_VALUE;
      for (int i=0; i<3; ++i) E[i] = DEF_VALUE;
      electronNumberDensity = DEF_VALUE;
      electronTemperature = DEF_VALUE;
      ionPolytropicIndex = DEF_VALUE;
      ionMassDensity = DEF_VALUE;
      for (int i=0; i<3; ++i) V_plasma_SIM[i] = DEF_VALUE;
   }
   
} // namespace sep
