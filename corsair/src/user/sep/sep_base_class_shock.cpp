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

#include <ucd_mesh.h>

#include "sep_simcontrol.h"
#include "sep_base_class_shock.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   
   ShockBaseClass::ShockBaseClass() { 
      sim = NULL;
      simClasses = NULL;
      initialized = false;
   }
   
   ShockBaseClass::~ShockBaseClass() { 
      ShockBaseClass::finalize();
   }
   
   bool ShockBaseClass::finalize() {
      sim = NULL;
      simClasses = NULL;
      initialized = false;
      return true;
   }
   
   bool ShockBaseClass::getLocalSurfaces(Real t,std::vector<LocalSurface>& localSurfaces,uint32_t& N_shockCells,int refinements) {
      // Create shock mesh:
      if (initializeMesh(t,refinements) == false) {
	 simClasses->logger << "(SEP SHOCK BASE) ERROR: Failed to create shock mesh" << endl << write;
	 return false;
      }

      uint64_t N_surfaces = 0;
      const vector<Real>& faceData = getFaceData(N_surfaces);
      const vector<Real>& nodeCoordinates = getNodeCoordinates();
      const vector<uint32_t>& cellConnectivity = getCellConnectivity();
      N_shockCells = N_surfaces;
      
      // Iterate over all shock surface elements. If the element is on a cell hosted on 
      // this process (determined by cell centroid), inject particles to it:
      size_t connIndex = 0;
      for (size_t s=0; s<N_surfaces; ++s) {
	 const size_t areaIndex = s * ucdmesh::facedataelement::SIZE;
	 const size_t connSize = cellConnectivity[connIndex+1]+2;
	 
	 const size_t node1 = cellConnectivity[connIndex+2];
	 const size_t node2 = cellConnectivity[connIndex+3];
	 const size_t node3 = cellConnectivity[connIndex+4];
	 
	 Real centroid[3];
	 centroid[0] = (nodeCoordinates[3*node1+0] + nodeCoordinates[3*node2+0] + nodeCoordinates[3*node3+0])/3;
	 centroid[1] = (nodeCoordinates[3*node1+1] + nodeCoordinates[3*node2+1] + nodeCoordinates[3*node3+1])/3;
	 centroid[2] = (nodeCoordinates[3*node1+2] + nodeCoordinates[3*node2+2] + nodeCoordinates[3*node3+2])/3;
	 
	 Real normal[3];
	 normal[0] = faceData[areaIndex + ucdmesh::facedataelement::NORMAL_X];
	 normal[1] = faceData[areaIndex + ucdmesh::facedataelement::NORMAL_Y];
	 normal[2] = faceData[areaIndex + ucdmesh::facedataelement::NORMAL_Z];
	 
	 // Convert centroid position to logical coordinates and check that it 
	 // is inside simulation domain (coords are not inf):
	 Real logical[3];
	 Real xy,r,theta,phi;
	 Real spherCentroid[3];
	 switch (simControl.coordinateSystem) {
	  case sep::UNKNOWN:
	    finalizeMesh();
	    return false;
	    break;
	  case sep::CARTESIAN:
	    getLogicalCoordinates(sim,centroid,logical);
	    if (logical[0] == numeric_limits<Real>::infinity()) {
	       connIndex += connSize;
	       continue;
	    }
	    break;
	  case sep::CYLINDRICAL:
	    finalizeMesh();
	    return false;
	    break;
	  case sep::SPHERICAL:
	    // Need to convert to spherical coordinates first:
	    xy = centroid[0]*centroid[0] + centroid[1]*centroid[1];
	    r = xy + centroid[2]*centroid[2];
	    r = sqrt(r);
	    xy = sqrt(xy);
	    
	    theta = acos(centroid[2]/r);
	    phi = acos(centroid[0] / xy);
	    if (centroid[1] < 0.0) phi = -phi;
	    
	    spherCentroid[0] = r;
	    spherCentroid[1] = theta;
	    spherCentroid[2] = phi;
	    
	    // Check logical coords:
	    getLogicalCoordinates(sim,spherCentroid,logical);
	    if (logical[0] == std::numeric_limits<Real>::infinity()) {
	       connIndex += connSize;
	       continue;
	    }
	    break;
	  default:
	    finalizeMesh();
	    return false;
	    break;
	 }
	 
	 // Check that injection position is in upstream side:
	 Real d_shock = getSquaredDistanceToShock(t,logical);
	 bool acceptSurface = true;
	 while (d_shock <= 1.0) {
	    // Move injection point to the direction of shock normal:
      	    for (int i=0; i<3; ++i) centroid[i] += (0.1+1.001*fabs(d_shock))*normal[i];

	    switch (simControl.coordinateSystem) {
	     case sep::UNKNOWN:
	       finalizeMesh();
	       return false;
	       break;
	     case sep::CARTESIAN:
	       getLogicalCoordinates(sim,centroid,logical);
	       if (logical[0] == std::numeric_limits<Real>::infinity()) {
		  acceptSurface = false;
	       }
	       break;
	     case sep::CYLINDRICAL:
	       finalizeMesh();
	       return false;
	       break;
	     case sep::SPHERICAL:
	       xy = centroid[0]*centroid[0] + centroid[1]*centroid[1];
	       r = sqrt(xy + centroid[2]*centroid[2]);
	       xy = sqrt(xy);
	       theta = acos(centroid[2]/r);
	       phi = acos(centroid[0] / xy);
	       if (centroid[1] < 0.0) phi = -phi;
	       
	       spherCentroid[0] = r;
	       spherCentroid[1] = theta;
	       spherCentroid[2] = phi;
	       
	       getLogicalCoordinates(sim,spherCentroid,logical);
	       if (logical[0] == numeric_limits<Real>::infinity()) acceptSurface = false;
	       break;
	     default:
	       finalizeMesh();
	       return false;
	       break;
	    }
	    
	    if (acceptSurface == false) break;
	    d_shock = simControl.shock->getSquaredDistanceToShock(t,logical);
	 }

	 if (acceptSurface == false) {
	    connIndex += connSize;
	    continue;
	 }
	 
	 // Calculate cell indices:
	 int32_t i_block = static_cast<int32_t>(logical[0]) / block::WIDTH_X;
	 int32_t j_block = static_cast<int32_t>(logical[1]) / block::WIDTH_Y;
	 int32_t k_block = static_cast<int32_t>(logical[2]) / block::WIDTH_Z;
	               
	 // Calculate block global index and check if this process has it:
	 pargrid::CellID blockGID = block::calculateGlobalIndex(*sim,i_block,j_block,k_block);
	 pargrid::CellID blockLID = simClasses->pargrid.getLocalID(blockGID);
	 if (blockLID == pargrid::INVALID_CELLID) {
	    connIndex += connSize;
	    continue;
	 }
	 if (simClasses->pargrid.getHosts()[blockLID] != sim->mpiRank) {
	    connIndex += connSize;
	    continue;
	 }

	 // Surface accepted, add to output vector:
	 LocalSurface surface;
	 surface.localID = blockLID;
	 surface.globalID = blockGID;
	 for (int i=0; i<3; ++i) surface.position[i] = logical[i];
	 surface.area = faceData[areaIndex+ucdmesh::facedataelement::AREA];
	 surface.zoneIndex = s;
	 localSurfaces.push_back(surface);

	 connIndex += connSize;
      }

      if (finalizeMesh() == false) return false;
      return true;
   }

   bool ShockBaseClass::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      this->sim = &sim;
      this->simClasses = &simClasses;
      return true;
   }
   
   bool ShockBaseClass::isInitialized() const {
      return initialized;
   }

   /** Modify particle accumulation shape factors to be upwinded or downwinded near the shock front, 
    * if shock sharpening is used. These shape factors are only used for wave energy change accumulations.
    * @param cellRegions Array of size block::WIDTH_X+2 containing spatial cell shock regions.
    * @param indices Temporary accumulation block indices of the cell in which particle resides.
    * @param shapeFactors Array containing shape factors.
    * @param particleRegion Particle shock region.*/
   void calculateUpwindedParticleAccumShapeFactors(int32_t* cellRegions,const int32_t* indices,
						   Real* shapeFactors,const int32_t& particleRegion) {
      
      switch (simControl.order) {
       case 0:
	 std::cerr << "ERROR: NGP shape factors not implemented in calculateUpwindedParticleAccumShapeFactors" << std::endl;
	 exit(1);
	 break;
       case 1:
	 if (abs(cellRegions[indices[0]+0]) != particleRegion) {
	    shapeFactors[0] = 0;
	    shapeFactors[1] = 1;
	 }
	 if (abs(cellRegions[indices[0]+1]) != particleRegion) {
	    shapeFactors[0] = 1;
	    shapeFactors[1] = 0;
	 }
	 break;
       case 2:
	 std::cerr << "ERROR: TSC shape factors not implemented in calculateUpwindedparticleAccumShapeFactors" << std::endl;
	 exit(1);
	 break;
      }
   }

   /** Modify particle interpolation shape factors to be upwinded or downwinded near the shock front, 
    * if shock sharpening is used. These shape factors are only used to interpolate wave energy to particle 
    * phase-space position.
    * @param cellRegions Array of size block::WIDTH_X+2 containing spatial cell shock regions.
    * @param indices Temporary accumulation block indices of the cell in which particle resides.
    * @param shapeFactors Array containing shape factors.
    * @param particleRegion Particle shock region.*/
   void calculateUpwindedParticleInterpShapeFactors(int32_t* cellRegions,int32_t* indices,
						    Real* shapeFactors,const int32_t& particleRegion,const int32_t& i) {
      
      switch (simControl.order) {
       case 0:
	 std::cerr << "ERROR: NGP shape factors not implemented in calculateUpwindedParticleInterpShapeFactors" << std::endl;
	 exit(1);
	 break;
       case 1:
	 if (cellRegions[i] < 0) {
	    if (cellRegions[0] == particleRegion) {
	       indices[0] = 0;
	       shapeFactors[0] = 1;
	       shapeFactors[1] = 0;
	    } else {
	       indices[0] = 1;
	       shapeFactors[0] = 0;
	       shapeFactors[1] = 1;
	    }
	 } else {
	    if (cellRegions[indices[0]+0] < 0) {
	       shapeFactors[0] = 0;
	       shapeFactors[1] = 1;
	    }
	    if (cellRegions[indices[0]+1] < 0) {
	       shapeFactors[0] = 1;
	       shapeFactors[1] = 0;
	    }
	 }
	 break;
       case 2:
	 std::cerr << "ERROR: TSC shape factors not implemented in calculateUpwindedParticleInterpShapeFactors" << std::endl;
	 exit(1);
	 break;
      }
   }
   
   /** Modify wave packet shape factors to be upwinded or downwinded near the shock front, if 
    * shock sharpening is used. Modified shape factors are used for wave energy accumulations and 
    * and wave energy change interpolations. If the cell in which wave packet resides has 
    * been flagged to be near the shock front, one on the shape factors is set to zero value to 
    * prevent wave energy from being accumulated to opposite side of shock.
    * @param cellRegions Array of size block::WIDTH_X+2 containing cell shock regions.
    * @param indices Array of size 3 containing indices of the cell within a temporary accumulation 
    * @param (or interpolation block) in which wave packet resides.
    * @param shapeFactors Array containing shape factors.
    * @param particleRegion Shock region of the wave packet.*/
   void calculateUpwindedWaveShapeFactors(int32_t* cellRegions,const int32_t* indices,
					  Real* shapeFactors,const int32_t& particleRegion) {
      
      switch (simControl.order) {
       case 0:
	 std::cerr << "ERROR: NGP shape factors not implemented in calculateUpwindedWaveShapeFactors" << std::endl;
	 exit(1);
	 break;
       case 1:
	 if (abs(cellRegions[indices[0]+0]) != particleRegion) {
	    shapeFactors[0] = 0;
	 }
	 if (abs(cellRegions[indices[0]+1]) != particleRegion) {
	    shapeFactors[1] = 0;
	 }
	 break;
       case 2:
	 std::cerr << "ERROR: TSC shape factors not implemented in calculateUpwindedWaveShapeFactors" << std::endl;
	 exit(1);
	 break;
      }
   }

   /** Classify spatial cells according to their shock regions (upstream, downstream, ...) in 1D.
    * This function is mostly used in conjunction with the temporary blocks used in accumulations 
    * and interpolations, i.e., the size of cellRegions array has to be block::WIDTH_X+2.
    * For each cell, shock regions are calculated for its +/- x-faces (nodes). If both faces have the 
    * same shock region, the that will be the cell's region. If face regions differ, the region of cell 
    * will be the region of cell centroid.
    * @param t Simulation time.
    * @param cellRegions Array of size block::WIDTH_X+2 where cell regions are written to.
    * @param blockIndices Array of size 3 that contains temporary block's indices.
    * @return If true, shock front is inside this block.*/
   bool classifyShockedCells(const Real& t,int32_t* cellRegions,const uint32_t* blockIndices) {
      bool shockedBlock = false;
      
      Real pos[3];
      pos[0] = blockIndices[0]*block::WIDTH_X - 1.0;
      pos[1] = blockIndices[1]*block::WIDTH_Y + 0.5;
      pos[2] = blockIndices[2]*block::WIDTH_Z + 0.5;

      int i=0;
      int prevRegion = simControl.shock->getShockRegion(t,pos);
      int nextRegion;
      ++i;
      do {
	 pos[0] = blockIndices[0]*block::WIDTH_X+i-1.0;
	 nextRegion = simControl.shock->getShockRegion(t,pos);

	 if (prevRegion == nextRegion) {
	    cellRegions[i-1] = nextRegion;
	 } else {
	    pos[0] = blockIndices[0]*block::WIDTH_X+i-1.5;
	    const Real centroidRegion = simControl.shock->getShockRegion(t,pos);
	    cellRegions[i-1] = -1*centroidRegion;
	    shockedBlock = true;
	 }
	 ++i;
      } while (i <= block::WIDTH_X+2);
      
      return shockedBlock;
   }
   
} // namespace sep

