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

#ifndef SEP_OPERATOR_SHOCK_MESH_H
#define SEP_OPERATOR_SHOCK_MESH_H

#include <vector>
#include <dataoperator.h>

namespace sep {

   namespace shockmesh {
      void getBlockIDs(Simulation* sim,SimulationClasses* simClasses,const Real* nodeCoords,
		       pargrid::CellID& blockLID,pargrid::CellID& blockGID);
      void transformToCartesian(Simulation* sim,const Real* pos,const Real* vectorIn,Real* vectorout);
      bool writeVariable(Simulation* sim,SimulationClasses* simClasses,bool nodeCentered,
			 const std::vector<Real>& values,size_t vectorSize,const std::string& name);
   }
   
   class OperatorShockMesh: public DataOperator {
    public:
      OperatorShockMesh();
      ~OperatorShockMesh();
      
      bool finalize();
      std::string getName() const;
      bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
      bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);

    private:
      #if PROFILE_LEVEL > 0
         int profTotalTime;
      #endif

      //Real countReflectedElectrons(Real numberDensity,Real temperature,Real mu_lc,Real V_HT_par);
      //void getBlockIDs(const Real* nodeCoords,pargrid::CellID& blockLID,pargrid::CellID& blockGID) const;
      //void transformToCartesian(const Real* pos,const Real* vectorIn,Real* vectorout) const;
      //bool writeElectronEmission(const std::vector<Real>& nodeCoords,const std::string& meshName);
      //bool writeEmissionMesh(const std::string& meshName);
      //bool writeVariable(bool nodeCentered,const std::vector<Real>& values,size_t vectorSize,const std::string& name);
      bool writeGasCompressionRatio(const std::vector<Real>& nodeCoords);
      bool writeLogicalCoordinates(const std::vector<Real>& nodeCoords);
      bool writeMachNumbers(const std::vector<Real>& nodeCoords);
      bool writeMagneticField(const std::vector<Real>& nodeCoords);
      bool writePlasmaVelocityHT(const std::vector<Real>& nodeCoords);
      bool writePlasmaVelocityShockFrame(const std::vector<Real>& nodeCoords);
      bool writePlasmaVelocitySNIF(const std::vector<Real>& nodeCoords);
      bool writeShockNormals(const std::vector<Real>& nodeCoords);
      bool writeShockPotential(const std::vector<Real>& nodeCoords);
      bool writeShockVelocity(const std::vector<Real>& nodeCoords);
      bool writeVariableMagneticCompressionRatio(const std::vector<Real>& nodeCoords);
      /*
      uint32_t xFreqBins;
      uint32_t yFreqBins;
      std::vector<Real> emission;
      Real freqMax;
      Real freqMin;
      Real dx;
       */
   };
   
} // namespace sep

#endif
