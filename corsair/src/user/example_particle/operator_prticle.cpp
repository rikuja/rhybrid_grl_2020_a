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

#include "operator_prticle.h"
#include "simcontrol.h"

using namespace std;

PrticleOperator::PrticleOperator(): DataOperator() {
   #ifdef PROFILE
      totalTimeID = -1;
   #endif
}

PrticleOperator::~PrticleOperator() {finalize();}

bool PrticleOperator::finalize() {return true;}

std::string PrticleOperator::getName() const {
   return "particles";
}

bool PrticleOperator::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleList>& particleLists) {
   bool success = true;
   this->sim = &sim;
   this->simClasses = &simClasses;
   return success;
}

bool PrticleOperator::writeData(const std::string& spatMeshName,const std::vector<ParticleList>& particles) {
   bool success = true;
   #ifdef PROFILE
      profile::start("Particles",totalTimeID);
   #endif

   if (SimControl::plist != NULL) {
      if (SimControl::plist->writeParticles(spatMeshName) == false) success = false;
   }
   
   #ifdef PROFILE
      profile::stop();
   #endif
   return success;
}

