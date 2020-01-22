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

#ifndef SEP_BASE_CLASS_SHOCK_H
#define SEP_BASE_CLASS_SHOCK_H

#include <vector>

#include <vlsv_common.h>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

namespace sep {

   namespace gascompression {
      enum Model {
	 CONSTANT,
	 OBLIQUE,
	 PERPENDICULAR
      };
   }

   /** Definition of a shock surface element whose centroid is 
    * inside a block hosted on this process.*/
   struct LocalSurface {
      pargrid::CellID localID;  /**< Block local ID.*/
      pargrid::CellID globalID; /**< Block global ID.*/
      Real position[3];         /**< Surface centroid (logical) coordinates.*/
      Real area;                /**< Surface area.*/
      uint32_t zoneIndex;
   };

   class ShockBaseClass {
    public:
      ShockBaseClass();
      virtual ~ShockBaseClass();

      virtual bool finalize();
      virtual bool getLocalSurfaces(Real t,std::vector<LocalSurface>& localSurfaces,uint32_t& N_shockMeshCells,int refinements);
      virtual bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr);

      virtual Real getAllowedDistance() const =0;
      virtual Real getDistanceToShock(Real t,const Real* pos) =0;
      virtual Real getGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) =0;
      virtual void getLocalShockVelocity(Real t,const Real* pos,Real* shockVelocity) =0;
      virtual void getLogicalNodeCoordinates(std::vector<Real>& nodeCoords) =0;
      virtual Real getMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) =0;
      virtual void getShockNormal(Real t,const Real* pos,Real* normal) =0;
      virtual int32_t getShockRegion(Real t,const Real* pos) =0;
      virtual Real getSquaredDistanceToShock(Real t,const Real* pos) =0;
      virtual bool inDownstream(Real t,const Real* pos) =0;
      bool isInitialized() const;

      //***** FUNCTIONS RELATED TO PLOTTING AND *****//
      //***** GENERATION OF SHOCK MESH          *****//
      virtual bool finalizeMesh() = 0;
      virtual std::vector<uint32_t>& getCellConnectivity() =0;
      virtual std::vector<uint32_t>& getCellOffsets() =0;
      virtual std::vector<Real>& getFaceData(uint64_t& N_surfaces) =0;
      virtual std::vector<Real>& getNodeCoordinates() =0;
      virtual uint32_t getNumberOfCells() const =0;
      virtual uint32_t getNumberOfNodes() const =0;
      virtual bool initializeMesh(Real t,int refinements) = 0;

    protected:
      gascompression::Model gasCompressionModel;
      bool initialized;
      Simulation* sim;
      SimulationClasses* simClasses;
   };
   
   void calculateUpwindedParticleAccumShapeFactors(int32_t* cellRegions,const int32_t* indices,Real* shapeFactors,const int32_t& particleRegion);
   void calculateUpwindedParticleInterpShapeFactors(int32_t* cellRegions,int32_t* indices,Real* shapeFactors,const int32_t& particleRegion,const int32_t& i);
   void calculateUpwindedWaveShapeFactors(int32_t* cellRegions,const int32_t* indices,Real* shapeFactors,const int32_t& particleRegion);
   bool classifyShockedCells(const Real& t,int32_t* cellRegions,const uint32_t* blockIndices);
   
} // namespace sep

#endif
