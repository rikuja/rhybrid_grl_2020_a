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

#ifndef SEP_SHOCK_PLANAR_H
#define SEP_SHOCK_PLANAR_H

#include <vector>

#include <vlsv_common.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_base_class_shock.h"

namespace sep {

   namespace doubleshock {
      enum DoubleShockRegions {
	 DUMMY,
	 UPSTREAM,                /**< Region upstream to leading shock.*/
	 IN_BETWEEN,              /**< Region between leading and trailing shocks.*/
	 DOWNSTREAM               /**< Region downstream to trailing shock.*/
      };
   }
   
   class ShockPlanar: public ShockBaseClass {
    public:
      ShockPlanar();
      ~ShockPlanar();

      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

      Real getAllowedDistance() const;
      Real getDistanceToShock(Real t,const Real* pos);
      Real getGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos);
      void getLocalShockVelocity(Real t,const Real* pos,Real* shockVelocity);
      void getLogicalNodeCoordinates(std::vector<Real>& nodeCoords);
      Real getMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos);
      void getShockNormal(Real t,const Real* pos,Real* normal);
      void getShockNormal(Real t,const Real* pos,const Real* shockCentroid,Real* normal);
      int32_t getShockRegion(Real t,const Real* pos);
      void getShockState(Real t,Real* centroid,Real& radius);
      Real getSquaredDistanceToShock(Real t,const Real* pos);
      bool inDownstream(Real t,const Real* pos);
      
      bool finalizeMesh();
      std::vector<uint32_t>& getCellConnectivity();
      std::vector<uint32_t>& getCellOffsets();
      std::vector<Real>& getFaceData(uint64_t& N_surfaces);
      std::vector<Real>& getNodeCoordinates();
      uint32_t getNumberOfCells() const;
      uint32_t getNumberOfNodes() const;
      bool initializeMesh(Real t,int refinements);
      
    private:
      Real leadingCentroid[3];                /**< Leading shock centroid position at current simulation time.*/
      Real leadingInitialCentroid[3];         /**< Leading shock centroid position at t=0.*/
      Real leadingNormal[3];                  /**< Leading shock normal in simulation frame.*/
      Real leadingVelocity_SIM[3];            /**< Leading shock velocity in simulation frame.*/
      
      Real trailingCentroid[3];               /**< Trailing shock centroid position at current simulation time.*/
      Real trailingInitialCentroid[3];        /**< Trailing shock centroid position at t=0.*/
      Real trailingNormal[3];                 /**< Trailing shock normal in simulation frame.*/
      Real trailingVelocity_SIM[3];           /**< Trailing shock velocity in simulation frame.*/
      Real trailingRelativeCentroid[3];       /**< Trailing shock's relative position wrt leading shock at t=0.*/
      Real relativeVelocity_SIM[3];           /**< Velocity of trailing shock relative to leading shock.*/

      std::vector<uint32_t> cellConnectivity;
      std::vector<uint32_t> cellOffsets;
      std::vector<Real> faceData;
      std::vector<Real> nodeCoordinates;
      Real leading_R_gas;                     /**< Gas compression ratio, only used for 'constant' shock model.*/
      Real leading_R_magnetic;                /**< Magnetic compression ratio, only used for 'constant' shock model.*/
      Real trailing_R_gas;                    /**< Gas compression ratio, only used for 'constant' shock model.*/
      Real trailing_R_magnetic;               /**< Magnetic compression ratio, only used for 'constant' shock model.*/
      
      bool addConfigFileItems(ConfigReader& cr);
      Real calculateGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* position);
      Real calculateMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* position);
      void calculateShockCentroids(Real t);
      bool createShockMesh(Real t,int refinements);
   };

   ShockBaseClass* PlanarShockMaker();
   
} // namespace sep

#endif
