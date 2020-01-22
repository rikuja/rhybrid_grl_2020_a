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

#ifndef SEP_OPERATOR_ACCUMULATION_BASE_H
#define SEP_OPERATOR_ACCUMULATION_BASE_H

#include <stdint.h>
#include <map>
#include <vector>

#include <dataoperator.h>

namespace sep {
   
   class OperatorAccumulationBase: public DataOperator {
    public:
      OperatorAccumulationBase(uint32_t vectorSize,uint32_t order);
      virtual ~OperatorAccumulationBase();
      
      bool finalize();
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);
   
    protected:
      
      #if PROFILE_LEVEL > 0
         std::string profileName;
         int accumulationID;
         int bufferCopyID;
         int densityOperatorID;
         int mpiOverheadID;
         int mpiWaitID;
         int mpiWaitSendsID;
         void setProfileName(const std::string& name);
      #endif

      uint32_t getOrder() const;
      uint32_t getVectorSize() const;
      bool writeData(uint32_t arrayIndex,const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists);

      // Get a unique name of DataOperator, used to register this 
      // DataOperator to DataOperatorContainer. Note: this is defined 
      // in src/include/dataoperator.h file.
      //virtual std::string getName() const =0;
      
      virtual void accumulateBlock(pargrid::CellID blockLID,Real* accumArray,const std::vector<ParticleListBase*>& particleLists) =0;
      virtual uint32_t getNumberOfArrays() const =0;
      virtual std::string getOutputName(uint32_t arrayIndex) const =0;
      virtual std::string getOutputUnits(uint32_t arrayIndex) const =0;
      virtual bool preProcessData(Real* accumArray);
      virtual bool postProcessData(Real* accumArray);
      virtual bool setAccumulatedArray(uint32_t arrayIndex) =0;

    protected:
      Real* accumArray;
      uint32_t currentArray;     /**< Index of array currently being accumulated.*/
      uint32_t N_arrays;         /**< Number of arrays written to file.*/
      const uint32_t order;      /**< Order of accuracy of particle shape factors used in accumulation.*/
      const uint32_t vectorSize; /**< Number of elements in each data vector (3 for 3D vectors).*/
   };

} // namespace sep
   
#endif
