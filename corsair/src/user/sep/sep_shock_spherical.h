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

#ifndef SEP_SHOCK_SPHERICAL_H
#define SEP_SHOCK_SPHERICAL_H

#include <vector>

#include <vlsv_common.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include "sep_base_class_shock.h"

namespace sep {

   namespace sphericalshock {
      enum Regions {
	 DUMMY,
	 UPSTREAM,
	 DOWNSTREAM
      };
   }
   
   class ShockSpherical: public ShockBaseClass {
    public:
      ShockSpherical();
      ~ShockSpherical();

      bool finalize();
      bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

      Real getAllowedDistance() const;
      Real getDistanceToShock(Real t,const Real* pos);
      Real getDistanceToShock(Real t,const Real* pos,const Real* shockCentroid);
      Real getGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos);
      void getLocalShockVelocity(Real t,const Real* pos,Real* shockVelocity);
      void getLogicalNodeCoordinates(std::vector<Real>& nodeCoords);
      Real getMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos);
      void getShockNormal(Real t,const Real* pos,Real* normal);
      void getShockNormal(Real t,const Real* pos,const Real* shockCentroid,Real* normal);
      int32_t getShockRegion(Real t,const Real* pos);
      void getShockState(Real t,Real* centroid,Real& radius);
      Real getSquaredDistanceToShock(Real t,const Real* pos);
      Real getSquaredDistanceToShock(Real t,const Real* pos,const Real* shockCentroid);
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
      std::vector<uint32_t> cellConnectivity;
      std::vector<uint32_t> cellOffsets;
      Real centroid[3];                       /**< Shock centroid position at time t=0 in Cartesian coordinates.*/
      std::vector<Real> faceData;
      Real initialRadius;
      std::vector<Real> nodeCoordinates;
      Real accelerationCentroid;              /**< Centroid acceleration between t_shock < t < t_shock+t_cruising.*/
      Real R_gas;                             /**< Gas compression ratio, only used for 'constant' shock model.*/
      Real R_magnetic;                        /**< Magnetic compression ratio, only used for 'constant' shock model.*/
      Real speed;                             /**< Shock centroid speed.*/
      Real t_cruising;                        /**< Time it takes for CME to reach cruising speed.*/
      Real V_direction[3];
      Real V_cartDirection[3];                /**< Direction of shock centroid propagation in Cartesian coordinates.*/
      Real V_expansion;                       /**< Shock lateral expansion speed.*/

      bool addConfigFileItems(ConfigReader& cr);
      Real calculateGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* position);
      Real calculateMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* position);
      void calculateRelativePosition(const Real* position,const Real* cartesianCentroid,Real* relPos);
      void calculateShockCentroid(Real t,Real* currentCentroid);
      Real calculateShockRadius(Real t);
      Real calculateSquaredShockRadius(Real t);
      bool createShockMesh(Real t,int refinements);
      void writeValues();

      template<typename REAL> void cylindricalToCylindrical(const Real* pos,const Real* vectorIn,Real* vectorOut) const;
      template<typename REAL> void sphericalToSpherical(const Real* pos,const Real* vectorIn,Real* vectorOut) const;
   };

   ShockBaseClass* SphericalShockMaker();
   
   template<typename REAL> inline
   void ShockSpherical::cylindricalToCylindrical(const Real* RESTRICT pos,const Real* RESTRICT vectorIn,Real* RESTRICT vectorOut) const {
      const uint32_t J = static_cast<uint32_t>(pos[1]);
      const Real phi = sim->y_crds_node[J] + (pos[1]-J)*sim->dy_cell[J];
      
      const Real cos_phi = cos(phi);
      const Real sin_phi = sin(phi);
      
      vectorOut[0] =  vectorIn[0]*cos_phi + vectorIn[1]*sin_phi;
      vectorOut[1] = -vectorIn[0]*sin_phi + vectorIn[1]*cos_phi;
      vectorOut[2] =  vectorIn[2];
   }

   template<typename REAL> inline
   void ShockSpherical::sphericalToSpherical(const Real* RESTRICT pos,const Real* RESTRICT vectorIn,Real* RESTRICT vectorOut) const {
      const uint32_t I = static_cast<uint32_t>(pos[0]);
      const uint32_t J = static_cast<uint32_t>(pos[1]);
      const uint32_t K = static_cast<uint32_t>(pos[2]);
      Real physPos[3];
      physPos[0] = sim->x_crds_node[I] + (pos[0]-I)*sim->dx_cell[I];
      physPos[1] = sim->y_crds_node[J] + (pos[1]-J)*sim->dy_cell[J];
      physPos[2] = sim->z_crds_node[K] + (pos[2]-K)*sim->dz_cell[K];
      
      const Real cos_theta = cos(physPos[1]);
      const Real sin_theta = sin(physPos[1]);
      const Real cos_phi   = cos(physPos[2]);
      const Real sin_phi   = sin(physPos[2]);
      
      vectorOut[0] = vectorIn[0]*sin_theta*cos_phi + vectorIn[1]*sin_theta*sin_phi + vectorIn[2]*cos_theta;
      vectorOut[1] = vectorIn[0]*cos_theta*cos_phi + vectorIn[1]*cos_theta*sin_phi - vectorIn[2]*sin_theta;
      vectorOut[2] = -vectorIn[0]*sin_phi          + vectorIn[1]*cos_phi;
   }
   
} // namespace sep

#endif
