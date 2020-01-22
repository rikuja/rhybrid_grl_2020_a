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

#ifndef SEP_LAGR_ACCELERATOR_H
#define SEP_LAGR_ACCELERATOR_H

#include <cstdlib>
#include <iostream>
#include <cmath>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <linear_algebra.h>
#include <constants.h>

#include "sep_simcontrol.h"
#include "sep_fields_container.h"
#include "sep_lagr_definition.h"

namespace sep {

   extern sep::FieldsContainer fieldsContainer;
   extern sep::SimControl simControl;
   
   template<class SPECIES,class PARTICLE>
   class LagrangianAccelerator {
    public:
      
      LagrangianAccelerator();
      ~LagrangianAccelerator();
      
      bool finalize();
      void get(pargrid::CellID blockLID,Real t,int direction,Real* state,Real* acceleration);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      
    private:
      bool initialized;

      Simulation* sim;
      SimulationClasses* simClasses;
      
      Real R;
      Real THETA;
      uint32_t I;
      uint32_t J;
   };
   
   template<class SPECIES,class PARTICLE> inline
   LagrangianAccelerator<SPECIES,PARTICLE>::LagrangianAccelerator() { 
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   LagrangianAccelerator<SPECIES,PARTICLE>::~LagrangianAccelerator() { 
      finalize();
   }
   
   template<class SPECIES,class PARTICLE> inline
   bool LagrangianAccelerator<SPECIES,PARTICLE>::finalize() {
      initialized = false;
      sim = NULL;
      simClasses = NULL;
      return true;
   }
 
   /** Get guiding center acceleration at given position.
    * state[sep::particle::XCRD] = R,
    * state[sep::particle::YCRD] = PHI,
    * state[sep::particle::ZCRD] = Z.
    */
   template<class SPECIES,class PARTICLE> inline
   void LagrangianAccelerator<SPECIES,PARTICLE>::get(pargrid::CellID blockLID,Real t,int direction,Real* RESTRICT state,Real* RESTRICT acceleration) {
      // Get Alfven speed, its parallel derivative, and magnetic field at wave packet position
      // NOTE: Alfven speed includes plasma speed:
      Real B_unit[3];
      Real V_wave[3];
      Real dV_wave;
      Real V_alfven;
      (*simControl.fieldsGetState)(blockLID,t,state,B_unit,V_wave,dV_wave,V_alfven,direction);
      unitVector<3>(B_unit);
      
      // Calculate coordinate velocity:
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 acceleration[sep::lagr::XCRD] = V_wave[0];
	 acceleration[sep::lagr::YCRD] = V_wave[1];
	 acceleration[sep::lagr::ZCRD] = V_wave[2];
	 break;
       case sep::CARTESIAN:
	 acceleration[sep::lagr::XCRD] = V_wave[0];
	 acceleration[sep::lagr::YCRD] = V_wave[1];
	 acceleration[sep::lagr::ZCRD] = V_wave[2];
	 break;
       case sep::CYLINDRICAL:
	 I = static_cast<uint32_t>(state[sep::particle::XCRD]);
	 R = sim->x_crds_node[I] + (state[sep::particle::XCRD] - I) * sim->dx_cell[I];
	 acceleration[sep::lagr::XCRD] = V_wave[0];
	 acceleration[sep::lagr::YCRD] = V_wave[1] / R;
	 acceleration[sep::lagr::ZCRD] = V_wave[2];
	 break;
       case sep::SPHERICAL:
	 // Guard against the case that wave packet ends up outside simulation domain:
	 I = static_cast<uint32_t>(state[sep::lagr::XCRD]);
	 if (I < 0) {
	    I = 0;
	 } else if (I >= sim->x_blocks*block::WIDTH_X) {
	    I = sim->x_blocks*block::WIDTH_X-1;
	 } else {

	 }

	 J = static_cast<uint32_t>(state[sep::lagr::YCRD]);
	 R     = sim->x_crds_node[I] + (state[sep::particle::XCRD] - I) * sim->dx_cell[I];
	 THETA = sim->y_crds_node[J] + (state[sep::lagr::YCRD] - J) * sim->dy_cell[J];
	 acceleration[sep::lagr::XCRD] = V_wave[0];
	 acceleration[sep::lagr::YCRD] = V_wave[1] / R;
	 acceleration[sep::lagr::ZCRD] = V_wave[2] / ( R*sin(THETA) );
	 break;
       default:
	 acceleration[sep::lagr::XCRD] = V_wave[0];
	 acceleration[sep::lagr::YCRD] = V_wave[1];
	 acceleration[sep::lagr::ZCRD] = V_wave[2];
	 break;
      }

      switch (simControl.wavelengthMeshType) {
       case sep::mesh::UNKNOWN:
	 std::cerr << "(SEP LAGR ACCELERATOR) ERROR: Wavelength mesh type is 'UNKNOWN'" << std::endl;
	 std::cerr << "\t This is likely due to programming error" << std::endl;
	 exit(1);
	 break;
       case sep::mesh::LINEAR:
	 acceleration[sep::lagr::LAMBDA] = state[sep::lagr::LAMBDA] * dV_wave;
	 break;
       case sep::mesh::LOGARITHMIC:
	 acceleration[sep::lagr::LAMBDA] = dV_wave;
	 break;
       default:
	 acceleration[sep::lagr::LAMBDA] = 0.0;
	 std::cerr << "(SEP LAGR ACCELERATOR) ERROR: Unknown wavelength mesh type" << std::endl;
	 exit(1);
	 break;
      }
      
      acceleration[sep::lagr::ENERGY] = 0.0;
   }

   template<class SPECIES,class PARTICLE> inline
   bool LagrangianAccelerator<SPECIES,PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      initialized = true;
      
      this->sim = &sim;
      this->simClasses = &simClasses;
      return initialized;
   }
   
} // namespace sep
   
#endif
