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

#ifndef SEP_OBJECT_WRAPPER_H
#define SEP_OBJECT_WRAPPER_H

#include <object_factory_generic.h>

#include "sep_fields_container.h"
#include "sep_distrib_energy_container.h"
#include "sep_distrib_pitch_container.h"
#include "sep_base_class_shock.h"
#include "sep_distrib_wave_energy_base_class.h"
#include "sep_base_class_particle_splitter.h"

namespace sep {

   struct ObjectWrapper {
      ObjectWrapper();
      ~ObjectWrapper();

      DistribEnergyContainer energyDistribContainer;
      DistribPitchContainer pitchDistribContainer;
      FieldsContainer fieldsContainer;

      ObjectFactoryGeneric<ShockBaseClass> shockFactory;
      ObjectFactoryGeneric<ParticleSplitterBase> splitterFactory;
      ObjectFactoryGeneric<WaveEnergySpectrumBaseClass> waveEnergySpectrumFactory;
      
    private:
      ObjectWrapper(const sep::ObjectWrapper& ow);
   };
   
   sep::ObjectWrapper& getObjectWrapper();
   
} // namespace sep

#endif
