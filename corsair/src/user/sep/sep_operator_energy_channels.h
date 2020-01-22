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

#ifndef SEP_OPERATOR_ENERGY_CHANNELS_H
#define SEP_OPERATOR_ENERGY_CHANNELS_H

#include <stdint.h>
#include <map>
#include <vector>

#include <sep_operator_accumulation_base.h>

namespace sep {

   namespace spacecraft {
      struct Instrument {
	 std::vector<Real> binWidths;                  /**< Energy bin widths in eV.*/
	 std::vector<std::string> channelNames;        /**< Names of energy channel.*/
	 bool divideByBinWidth;                        /**< If true, measured counts are divided by channel bin
							* width. Physical units will be 1/m^3/eV instead of 1/m^3.*/
	 bool energyPerNucleon;                        /**< If true, energy limits in minValues, maxValues are
							* per nucleon (energy * proton mass / particle mass).*/
	 std::vector<Real> maxValues;                  /**< Maximum energy for each channel.*/
	 std::vector<Real> minValues;                  /**< Minimum energy for each channel.*/
	 std::string name;                             /**< Name of the Instrument.*/
	 std::vector<std::string> speciesNames;        /**< Names of accumulated particle species.*/
      };
      
      bool createInstrument(SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegion,
			    const std::string& name,Instrument& instrument);
   }

   class OperatorEnergyChannels: public OperatorAccumulationBase {
    public:
      OperatorEnergyChannels(int32_t order=1);
      ~OperatorEnergyChannels();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   
    private:

      std::vector<spacecraft::Instrument> instruments;
      std::vector<std::pair<size_t,size_t> > instrumentIndices;

      uint32_t channelIndex;                             /**< Index of channel within the current Instrument whose data is 
							  * being calculated. This indexes vectors defined in spacecraft::Instrument.*/
      uint32_t currentChannel;                           /**< Index of channel that is currently accumulated.*/
      uint32_t currentInstrument;                        /**< Index of instrument that is currently accumulated.*/
      uint32_t instrumentIndex;                          /**< Index of Instrument whose data is currently being calculated,
							  * this is an index to OperatorEnergyChannels::instruments vector.*/
      uint32_t totalNumberOfArrays;                      /**< Total number of arrays written to output files,
							  * sum of number of instruments times number of 
							  * channels in each instrument.*/

      void accumulateBlock(pargrid::CellID blockLID,Real* accumArray,const std::vector<ParticleListBase*>& particleLists);
      bool addConfigFileItems(ConfigReader& cr);
      bool createInstrument(SimulationClasses& simClasses,ConfigReader& cr,const std::string& name);
      uint32_t getNumberOfArrays() const;
      std::string getOutputName(uint32_t arrayIndex) const;
      std::string getOutputUnits(uint32_t arrayIndex) const;
      bool postProcessData(Real* array);
      bool setAccumulatedArray(uint32_t arrayIndex);
   };

} // namespace sep
   
#endif
