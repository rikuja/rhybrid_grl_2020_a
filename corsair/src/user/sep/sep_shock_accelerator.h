/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#ifndef SEP_SHOCK_ACCELERATOR_H
#define SEP_SHOCK_ACCELERATOR_H

#include <cstring>
#include <definitions.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_particle_species.h"

namespace sep {

   namespace shockaccelerator {
      enum Elements {
	 XPOS,
	 YPOS,
	 ZPOS,
	 V_PAR,
	 MU,
	 DATASIZE
      };
   }
   
   struct ParticleParameters {
      Real state[sep::shockaccelerator::DATASIZE];
   };
   
   struct PlasmaParameters {
      Real B1_norm;
      Real B1_tang;
      Real V1_plasma_norm;
      Real V1_plasma_tang;
      Real R_gas;
      Real R_magn;
      Real L_shock;
   };

   class ShockAccelerator {
    public:
      ShockAccelerator();
      ~ShockAccelerator();

      Real getReturnedParticle(const PlasmaParameters& plasma,ParticleParameters& particle,const Real& V1_alfven,Real& pitch);
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);
      
      bool runTests();
      void setSpecies(const sep::Species& species);
      int solveShockEncounterAnalytic(const PlasmaParameters& plasma,ParticleParameters& particle);
      int solveShockEncounterGC(const PlasmaParameters& plasma,ParticleParameters& particle);
      int solveShockEncounterLorentz(const PlasmaParameters& plasma,ParticleParameters& particle);

    private:
      Real phiFactor;
      Simulation* sim;
      SimulationClasses* simClasses;
      sep::Species species;
      
      void borisBunemanAhead(Real* state,const Real* E,const Real* B,Real eforce);
      void borisBunemanBehind(Real* state,const Real* E,const Real* B,Real eforce);
      void calculateWebbResult(const PlasmaParameters& plasma,const ParticleParameters& particle,
			       Real U_ini,Real& dU,Real& pitch);
      void getDownstreamValues(const PlasmaParameters& plasma,Real V1_alfven,Real& V2_wave_par,Real* B2,Real& B2_mag,
			       Real* V2_wave,Real& tan_psi2);
      void propagateBorisBuneman(Real* state,const Real* E,const Real* B,Real eforce,Real dt);
      bool testFermiI();
      bool testFermiI(Real U_ini,Real cos_psi1);
      bool testSetupUpstream(Real R_gas,Real B_mag,Real V_shock,Real cos_psi1,PlasmaParameters& plasma);
      bool testShockReturn();
      bool testShockReturn(Real U_ini,Real psi1,Real& P_return_avg,Real& dU_avg);
      bool testWebb();
      bool testWebb(Real V_plasma,Real energy,Real psi,const std::string& fname,bool ions,uint32_t N_particles,
		    int (ShockAccelerator::*solver)(const PlasmaParameters& plasma,ParticleParameters& particle));
      bool testWriteSampleFields();
      bool writeFields(const PlasmaParameters& params,const std::string& fileName);
   };
} // namespace sep

#endif