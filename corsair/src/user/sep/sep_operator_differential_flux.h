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

#ifndef SEP_OPERATOR_FLUX_H
#define SEP_OPERATOR_FLUX_H

#include <stdint.h>
#include <map>
#include <vector>

#include <sep_operator_accumulation_base.h>

namespace sep {

   namespace spacecraft {
      struct FluxInstrument {
	 std::vector<Real> binWidths;                  /**< Flux bin widths in Joules.*/
	 std::vector<std::string> channelNames;        /**< Names of energy channels.*/
	 bool divideByBinWidth;                        /**< If true, measured counts are divided by channel bin
							* width. Physical units will be 1/m^3/eV instead of 1/m^3.*/
	 bool energyPerNucleon;                        /**< If true, measured counts are divided by mass (in amu).*/
	 std::vector<Real> maxValues;                  /**< Maximum energy for each channel.*/
	 std::vector<Real> minValues;                  /**< Minimum energy for each channel.*/
	 std::string name;                             /**< Name of the Instrument.*/
	 std::vector<std::string> speciesNames;        /**< Names of accumulated particle species.*/

	 bool createInstrument(SimulationClasses& simClasses,ConfigReader& cr,
			       const std::string& prefix,const std::string& name,
			       std::vector<spacecraft::FluxInstrument>& instruments,
			       std::vector<std::pair<size_t,size_t> >& instrumentIndices);
      };
   }

   class OperatorDifferentialFlux: public OperatorAccumulationBase {
    public:
      OperatorDifferentialFlux(int32_t order=1);
      virtual ~OperatorDifferentialFlux();
      
      virtual bool finalize();
      virtual std::string getName() const;
      virtual bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   
    protected:

      std::vector<spacecraft::FluxInstrument> instruments;
      std::vector<std::pair<size_t,size_t> > instrumentIndices;

      Real normal[3];
      uint32_t channelIndex;                             /**< Index of channel within the current Instrument whose data is 
							  * being calculated. This indexes vectors defined in spacecraft::Instrument.*/
      uint32_t currentChannel;                           /**< Index of channel that is currently accumulated.*/
      uint32_t currentInstrument;                        /**< Index of instrument that is currently accumulated.*/
      uint32_t instrumentIndex;                          /**< Index of Instrument whose data is currently being calculated,
							  * this is an index to OperatorDifferentialFlux::instruments vector.*/
      uint32_t totalNumberOfArrays;                      /**< Total number of arrays written to output files,
							  * sum of number of instruments times number of 
							  * channels in each instrument.*/

      virtual void accumulateBlock(pargrid::CellID blockLID,Real* accumArray,const std::vector<ParticleListBase*>& particleLists);
      bool addConfigFileItems(ConfigReader& cr);
      uint32_t getNumberOfArrays() const;
      virtual std::string getOutputName(uint32_t arrayIndex) const;
      virtual std::string getOutputUnits(uint32_t arrayIndex) const;
      virtual bool preProcessData(Real* array);
      virtual bool postProcessData(Real* array);
      bool setAccumulatedArray(uint32_t arrayIndex);
   };

} // namespace sep
   
#endif
