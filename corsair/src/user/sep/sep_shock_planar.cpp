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

#include <linear_algebra.h>
#include <triangulated_plane.h>

#include "sep_simcontrol.h"
#include "sep_shock_planar.h"
#include "sep_mesh_logical.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

static const Real DEF_REAL = numeric_limits<Real>::infinity();
static const string PREFIX = "Shock.planar.double";

namespace sep {

   extern sep::SimControl simControl;

   ShockBaseClass* PlanarShockMaker() {return new ShockPlanar();}
   
   ShockPlanar::ShockPlanar() { 

   }
   
   ShockPlanar::~ShockPlanar() { 
      finalize();
   }
   
   bool ShockPlanar::addConfigFileItems(ConfigReader& cr) {      
      cr.add(PREFIX+".leading_shock.speed","Speed of leading shock in simulation frame (float).",DEF_REAL);
      cr.add(PREFIX+".leading_shock.centroid_x","Leading shock centroid position at t=0, x-component (float).",DEF_REAL);
      cr.add(PREFIX+".leading_shock.centroid_y","Leading shock centroid position at t=0, y-component (float).",DEF_REAL);
      cr.add(PREFIX+".leading_shock.centroid_z","Leading shock centroid position at t=0, z-component (float).",DEF_REAL);
      cr.add(PREFIX+".leading_shock.gas_compression_ratio","Leading shock gas compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      cr.add(PREFIX+".leading_shock.magnetic_compression_ratio","Leading shock magnetic compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      
      cr.add(PREFIX+".trailing_shock.speed","Speed of trailing shock in simulation frame (float).",DEF_REAL);
      cr.add(PREFIX+".trailing_shock.centroid_x","Trailing shock centroid position at t=0, x-component (float).",DEF_REAL);
      cr.add(PREFIX+".trailing_shock.centroid_y","Trailing shock centroid position at t=0, y-component (float).",DEF_REAL);
      cr.add(PREFIX+".trailing_shock.centroid_z","Trailing shock centroid position at t=0, z-component (float).",DEF_REAL);
      cr.add(PREFIX+".trailing_shock.gas_compression_ratio","Trailing shock gas compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      cr.add(PREFIX+".trailing_shock.magnetic_compression_ratio","Trailing shock magnetic compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      
      cr.add(PREFIX+".length_units","Units in which positions and radii are given, defaults to m (m/km/AU/RE/RS/RV)",string("m"));
      cr.add(PREFIX+".shock_model","Which shock model is used, 'constant', 'oblique' or 'perpendicular' (string).",string(""));
      return true;
   }

   Real ShockPlanar::calculateGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* RESTRICT pos) {
      // Get shock velocity:
      Real V_shock_SIM[3];
      getLocalShockVelocity(t,pos,V_shock_SIM);

      // Get shock normal:
      Real normal[3];
      getShockNormal(t,pos,normal);

      // Get fields at given position and calculate cosine of local obliquity angle:
      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(blockLID,t,pos,E,B,gradB);

      const Real B_norm_mag = dotProduct<3>(B,normal);
      const Real cos_psi1 = B_norm_mag / vectorMagnitude<3>(B);
      
      // Get plasma state and calculate Mach numbers:
      PlasmaState plasmaState;
      (*simControl.fieldsGetPlasmaState)(blockLID,t,pos,plasmaState);
      Real V_plasma_SRF[3];
      for (int i=0; i<3; ++i) V_plasma_SRF[i] = plasmaState.V_plasma_SIM[i] - V_shock_SIM[i];
      const Real V_plasma_norm = dotProduct<3>(V_plasma_SRF,normal);
      if (V_plasma_norm > 0) return 1.0;
      const Real V_plasma_norm2 = V_plasma_norm*V_plasma_norm;
      const Real alfvenMach2 = V_plasma_norm2 / plasmaState.alfvenSpeed2;
      const Real sonicMach2 = V_plasma_norm2 / plasmaState.soundSpeed2;
      
      // Solve gas compression ratio:
      calculateShockCentroids(t);
      switch (gasCompressionModel) {
       case gascompression::CONSTANT:
	 if (simControl.coordinateSystem == sep::CARTESIAN) {
	    if (pos[0] <= leadingCentroid[0]) return leading_R_gas;
	    else return trailing_R_gas;
	 } else {
	    if (pos[0] >= leadingCentroid[0]) return leading_R_gas;
	    else return trailing_R_gas;
	 }
	 break;
       case gascompression::OBLIQUE:
	 return sep::oblique::solveGasCompressionRatio(cos_psi1,sonicMach2,alfvenMach2,plasmaState.ionPolytropicIndex);
	 break;
       case gascompression::PERPENDICULAR:
	 return sep::perpendicular::solveGasCompressionRatio(cos_psi1,sonicMach2,alfvenMach2,plasmaState.ionPolytropicIndex);
	 break;
       default:
	 return -1.0;
	 break;
      }
   }
   
   Real ShockPlanar::calculateMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) {
      // Get shock velocity:
      Real V_shock_SIM[3];
      getLocalShockVelocity(t,pos,V_shock_SIM);
               
      // Get shock normal:
      Real normal[3];
      getShockNormal(t,pos,normal);

      // Get fields at given position and calculate cosine of local obliquity angle:
      Real E[3];
      Real B[3];
      Real gradB[9];
      (*simControl.fieldsGetFields)(blockLID,t,pos,E,B,gradB);
      
      const Real B_norm_mag = dotProduct<3>(B,normal);
      const Real cos_psi1 = B_norm_mag / vectorMagnitude<3>(B);
      
      // Get plasma state and calculate Mach numbers:
      PlasmaState plasmaState;
      (*simControl.fieldsGetPlasmaState)(blockLID,t,pos,plasmaState);
      Real V_plasma_SRF[3];
      for (int i=0; i<3; ++i) V_plasma_SRF[i] = plasmaState.V_plasma_SIM[i] - V_shock_SIM[i];
      const Real V_plasma_norm = dotProduct<3>(V_plasma_SRF,normal);
      if (V_plasma_norm > 0) return 1.0;
      const Real V_plasma_norm2 = V_plasma_norm*V_plasma_norm;
      const Real alfvenMach2 = V_plasma_norm2 / plasmaState.alfvenSpeed2;
      const Real sonicMach2 = V_plasma_norm2 / plasmaState.soundSpeed2;
      
      // Solve gas compression ratio and compute magnetic enhancement:
      Real R_gas = -1.0;
      Real R_magn = 1.0;
      calculateShockCentroids(t);
      switch (gasCompressionModel) {
       case gascompression::CONSTANT:
	 if (simControl.coordinateSystem == sep::CARTESIAN) {
	    if (pos[0] <= leadingCentroid[0]) return leading_R_magnetic;
	    else return trailing_R_magnetic;
	 } else {
	    if (pos[0] >= leadingCentroid[0]) return leading_R_magnetic;
	    else return trailing_R_magnetic;
	 }
	 break;
       case gascompression::OBLIQUE:
	 R_gas = sep::oblique::solveGasCompressionRatio(cos_psi1,sonicMach2,alfvenMach2,plasmaState.ionPolytropicIndex);
	 R_magn = (alfvenMach2 - cos_psi1*cos_psi1)/(alfvenMach2 - cos_psi1*cos_psi1*R_gas)*R_gas;
	 break;
       case gascompression::PERPENDICULAR:
	 R_gas = sep::perpendicular::solveGasCompressionRatio(cos_psi1,sonicMach2,alfvenMach2,plasmaState.ionPolytropicIndex);
	 R_magn = R_gas;
	 break;
       default:
	 R_gas = -1.0;
	 break;
      }

      return max(1.0,R_magn);
   }

   void ShockPlanar::calculateShockCentroids(Real t) {
      const Real dt = max((Real)0.0,t - simControl.t_setup);
      Real temporaryLeadingCentroid[3];
      Real temporaryTrailingCentroid[3];

      for (int i=0; i<3; ++i) {
	 temporaryLeadingCentroid[i] = leadingInitialCentroid[i]
	   + (leadingVelocity_SIM[i] - simControl.V_frame[i]) * dt;
      }
      for (int i=0; i<3; ++i) {
	 temporaryTrailingCentroid[i] = temporaryLeadingCentroid[i]
	   + trailingRelativeCentroid[i] + relativeVelocity_SIM[i] * t;
      }

      // Convert centroids into logical coordinates:
      getLogicalCoordinates(sim,temporaryLeadingCentroid,leadingCentroid);
      getLogicalCoordinates(sim,temporaryTrailingCentroid,trailingCentroid);
   }
   
   bool ShockPlanar::createShockMesh(Real t,int refinements) {
      #if PROFILE_LEVEL > 0
         static int profileMeshGeneration = -1;
         profile::start("Shock Mesh Generation",profileMeshGeneration);
      #endif

      // Clear previous mesh:
	{
	   vector<Real> dummy1;
	   vector<uint32_t> dummy2;
	   vector<uint32_t> dummy3;
	   vector<Real> dummy4;
	   nodeCoordinates.swap(dummy1);
	   cellConnectivity.swap(dummy2);
	   cellOffsets.swap(dummy3);
	   faceData.swap(dummy4);
	}
      
      size_t N_connections = 0;
      
      // The "padding" here ensures that all nodes in shock mesh 
      // are inside simulation domain. Node coordinates are used to 
      // find simulation cells where nodes are, and this ensures that 
      // cells are correctly identified:
      const Real Y_PADDING = (sim->y_crds_node[sim->y_blocks*block::WIDTH_Y] - sim->y_crds_node[0]) / 1e6;
      const Real Z_PADDING = (sim->z_crds_node[sim->z_blocks*block::WIDTH_Z] - sim->z_crds_node[0]) / 1e6;
      
      Real y_min = sim->y_crds_node[0] + Y_PADDING;
      Real y_max = sim->y_crds_node[sim->y_blocks*block::WIDTH_Y] - Y_PADDING;
      Real z_min = sim->z_crds_node[0] + Z_PADDING;
      Real z_max = sim->z_crds_node[sim->z_blocks*block::WIDTH_Z] - Z_PADDING;

      // Calculate leading and trailing shock centroid x-coordinates, 
      // taking into account that shocks are not allowed to move until 
      // current timestep exceeds simControl.setupTimesteps:
      const Real dt = max((Real)0.0,t - simControl.t_setup);
      Real x_leading,x_trailing;
      x_leading  = leadingInitialCentroid[0]  + leadingVelocity_SIM[0]  * dt;
      x_trailing = x_leading + trailingRelativeCentroid[0] + relativeVelocity_SIM[0] * sim->t;

      vector<Real> corners;
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 cerr << "(SEP SHOCK PLANAR) Unknown coordinate system" << endl;
	 exit(1);
	 break;
       case sep::CARTESIAN:
	 if (x_leading > sim->x_crds_node[0] && x_leading < sim->x_crds_node[sim->x_blocks*block::WIDTH_X]) {
	    corners.push_back(x_leading);
	    corners.push_back(y_min);
	    corners.push_back(z_min);
	    corners.push_back(x_leading);
	    corners.push_back(y_min);
	    corners.push_back(z_max);
	    corners.push_back(x_leading);
	    corners.push_back(y_max);
	    corners.push_back(z_max);
	    corners.push_back(x_leading);
	    corners.push_back(y_max);
	    corners.push_back(z_min);
	 }
	 if (x_trailing > sim->x_crds_node[0] && x_trailing < sim->x_crds_node[sim->x_blocks*block::WIDTH_X]) {
	    corners.push_back(x_trailing);
	    corners.push_back(y_min);
	    corners.push_back(z_min);
	    corners.push_back(x_trailing);
	    corners.push_back(y_min);
	    corners.push_back(z_max);
	    corners.push_back(x_trailing);
	    corners.push_back(y_max);
	    corners.push_back(z_max);
	    corners.push_back(x_trailing);
	    corners.push_back(y_max);
	    corners.push_back(z_min);
	 }
	 
	 break;
       case sep::CYLINDRICAL:
	 cerr << "(SEP SHOCK PLANAR) cylindrical coordinates not implemented" << endl;
	 exit(1);
	 break;
       case sep::SPHERICAL:
	 if (x_leading > sim->x_crds_node[0] && x_leading < sim->x_crds_node[sim->x_blocks*block::WIDTH_X]) {
	    corners.push_back(x_leading * sin(y_min) * cos(z_min) );
	    corners.push_back(x_leading * sin(y_min) * sin(z_min) );
	    corners.push_back(x_leading * cos(y_min) );
	    corners.push_back(x_leading * sin(y_max) * cos(z_min) );
	    corners.push_back(x_leading * sin(y_max) * sin(z_min) );
	    corners.push_back(x_leading * cos(y_max) );
	    corners.push_back(x_leading * sin(y_max) * cos(z_max) );
	    corners.push_back(x_leading * sin(y_max) * sin(z_max) );
	    corners.push_back(x_leading * cos(y_max) );
	    corners.push_back(x_leading * sin(y_min) * cos(z_max) );
	    corners.push_back(x_leading * sin(y_min) * sin(z_max) );
	    corners.push_back(x_leading * cos(y_min) );
	 }
	 if (x_trailing > sim->x_crds_node[0] && x_trailing < sim->x_crds_node[sim->x_blocks*block::WIDTH_X]) {
	    corners.push_back(x_trailing * sin(y_min) * cos(z_min) );
	    corners.push_back(x_trailing * sin(y_min) * sin(z_min) );
	    corners.push_back(x_trailing * cos(y_min) );
	    corners.push_back(x_trailing * sin(y_max) * cos(z_min) );
	    corners.push_back(x_trailing * sin(y_max) * sin(z_min) );
	    corners.push_back(x_trailing * cos(y_max) );
	    corners.push_back(x_trailing * sin(y_max) * cos(z_max) );
	    corners.push_back(x_trailing * sin(y_max) * sin(z_max) );
	    corners.push_back(x_trailing * cos(y_max) );
	    corners.push_back(x_trailing * sin(y_min) * cos(z_max) );
	    corners.push_back(x_trailing * sin(y_min) * sin(z_max) );
	    corners.push_back(x_trailing * cos(y_min) );
	 }
	 break;
       default:
	 cerr << "(SEP SHOCK PLANAR) Unknown coordinate system" << endl;
	 exit(1);
	 break;
      }
      
      // Exit if both shocks are outside simulation domain:
      if (corners.size() == 0) {
	 #if PROFILE_LEVEL > 0
	    profile::stop(); // profileMeshGeneration
         #endif
	 return true;
      }

      // Generate triangulated planar surface(s):
      vector<ucdmesh::Node> nodes;
      vector<ucdmesh::TriangularFace> faces;
      const bool result = ucdmesh::createPlane(corners,refinements,nodes,faces,N_connections);
      if (result == false) {
	 simClasses->logger << "(SEP SHOCK DOUBLE) ERROR: ucdmesh::createPlane failed" << endl << write;
	 #if PROFILE_LEVEL > 0
	    profile::stop(); // profileMeshGeneration
	 #endif
	 return result;
      }
      
      // Get face data (areas and unit vectors):
      ucdmesh::getPlanarFaceData(nodes,faces,faceData);

      // Temporary solution
      #warning FIXME inefficient method of shock mesh generation for plotting
      nodeCoordinates.resize(nodes.size()*3);
      
      switch (simControl.coordinateSystem) {
       case sep::UNKNOWN:
	 cerr << "(SEP SHOCK DOUBLE) ERROR: Unknown mesh geometry" << endl;
	 exit(1);
	 break;
       case sep::CARTESIAN:
	 for (size_t i=0; i<nodes.size(); ++i) {
	    nodeCoordinates[3*i+0] = nodes[i].x;
	    nodeCoordinates[3*i+1] = nodes[i].y;
	    nodeCoordinates[3*i+2] = nodes[i].z;
	 }
	 break;
       case sep::CYLINDRICAL:
	 cerr << "(SEP SHOCK DOUBLE) ERROR: Cylindrical double shock not implemented" << endl;
	 exit(1);
	 break;
       case sep::SPHERICAL:
	 for (size_t i=0; i<nodes.size(); ++i) {
	    nodeCoordinates[3*i+0] = nodes[i].x;
	    nodeCoordinates[3*i+1] = nodes[i].y;
	    nodeCoordinates[3*i+2] = nodes[i].z;
	 }
	 break;
       default:
	 for (size_t i=0; i<nodes.size(); ++i) {
	    nodeCoordinates[3*i+0] = nodes[i].x;
	    nodeCoordinates[3*i+1] = nodes[i].y;
	    nodeCoordinates[3*i+2] = nodes[i].z;
	 }
	 break;
      }

      cellConnectivity.resize(faces.size()*5);
      for (size_t i=0; i<faces.size(); ++i) {
	 cellConnectivity[5*i+0] = faces[i].cellType;
	 cellConnectivity[5*i+1] = faces[i].N_nodes;
	 cellConnectivity[5*i+2] = faces[i].node1;
	 cellConnectivity[5*i+3] = faces[i].node2;
	 cellConnectivity[5*i+4] = faces[i].node3;
      }
      
      // For now, cell offsets only needs to be of correct size:
      cellOffsets.resize(N_connections);
      
      #if PROFILE_LEVEL > 0
         profile::stop(); // profileMeshGeneration
      #endif
      
      return result;
   }
   
   bool ShockPlanar::finalize() {
      return true;
   }
   
   bool ShockPlanar::finalizeMesh() {
      vector<Real> dummy1;
      vector<uint32_t> dummy2;
      vector<uint32_t> dummy3;
      vector<Real> dummy4;
      nodeCoordinates.swap(dummy1);
      cellConnectivity.swap(dummy2);
      cellOffsets.swap(dummy3);
      faceData.swap(dummy4);
      return true;
   }

   Real ShockPlanar::getAllowedDistance() const {
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 return 1000.0;
      } else {
	 return 1000.0;
      }
   }

   vector<uint32_t>& ShockPlanar::getCellConnectivity() {return cellConnectivity;}
   
   vector<uint32_t>& ShockPlanar::getCellOffsets() {return cellOffsets;}

   Real ShockPlanar::getDistanceToShock(Real t,const Real* RESTRICT pos) {
      cerr << "getDistanceToShock called" << endl;
      exit(1);
      return NAN;
   }

   vector<Real>& ShockPlanar::getFaceData(uint64_t& N_surfaces) {
      N_surfaces = faceData.size() / ucdmesh::facedataelement::SIZE;
      return faceData;
   }

   Real ShockPlanar::getGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* RESTRICT pos) {
      return calculateGasCompressionRatio(blockLID,t,pos);
   }

   void ShockPlanar::getLocalShockVelocity(Real t,const Real* RESTRICT pos,Real* RESTRICT shockVelocity) {
      // Calculate shock centroids:
      calculateShockCentroids(t);
      
      // Return the velocity of shock closer to position pos:
      const Real distanceLeading  = fabs(leadingCentroid[0] - pos[0]);
      const Real distanceTrailing = fabs(trailingCentroid[0] - pos[0]);

      if (distanceLeading < distanceTrailing) {
	 for (int i=0; i<3; ++i) shockVelocity[i] = leadingVelocity_SIM[i] - simControl.V_frame[i];
      } else {
	 for (int i=0; i<3; ++i) shockVelocity[i] = trailingVelocity_SIM[i] - simControl.V_frame[i];
      }

   }
   
   void ShockPlanar::getLogicalNodeCoordinates(std::vector<Real>& nodeCoords) {
      nodeCoords.clear();
      nodeCoords.resize(nodeCoordinates.size());
      
      Real tmp[3];
      const size_t N_nodes = nodeCoordinates.size()/3;
      for (size_t node=0; node<N_nodes; ++node) {
	 // Transform Cartesian node coordinates into simulation coordinate system:
	 Real r,theta,phi,xy;
	 switch (simControl.coordinateSystem) {
	  case UNKNOWN:
	    simClasses->logger << "(SEP SHOCK DOUBLE) ERROR: Unknown coordinate system in getLogicalNodeCoordinates" << endl << write;
	    exit(1);
	    break;
	  case CARTESIAN:
	    tmp[0] = nodeCoordinates[3*node+0];
	    tmp[1] = nodeCoordinates[3*node+1];
	    tmp[2] = nodeCoordinates[3*node+2];
	    break;
	  case CYLINDRICAL:
	    break;
	  case SPHERICAL:
	    xy = nodeCoordinates[3*node+0]*nodeCoordinates[3*node+0]
	       + nodeCoordinates[3*node+1]*nodeCoordinates[3*node+1];
	    r = xy + nodeCoordinates[3*node+2]*nodeCoordinates[3*node+2];
	    r = sqrt(r);
	    xy = sqrt(xy);

	    theta = acos(nodeCoordinates[3*node+2]/r);
	    phi = acos(nodeCoordinates[3*node+0] / xy);
	    if (nodeCoordinates[3*node+1] < 0.0) phi = -phi;

	    tmp[0] = r;
	    tmp[1] = theta;
	    tmp[2] = phi;
	    break;
	  default:
	    simClasses->logger << "(SEP SHOCK DOUBLE) ERROR: Unknown coordinate system in getLogicalNodeCoordinates" << endl << write;
	    exit(1);
	    break;
	 }

	 getLogicalCoordinates(sim,tmp,&(nodeCoords[3*node+0]));
      }
   }
   
   Real ShockPlanar::getMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) {
      return calculateMagneticCompressionRatio(blockLID,t,pos);
   }

   vector<Real>& ShockPlanar::getNodeCoordinates() {return nodeCoordinates;}
   
   uint32_t ShockPlanar::getNumberOfCells() const {return cellConnectivity.size()/5;}
   
   uint32_t ShockPlanar::getNumberOfNodes() const {return nodeCoordinates.size()/3;}

   void ShockPlanar::getShockNormal(Real t,const Real* RESTRICT pos,Real* RESTRICT normal) {
      for (int i=0; i<3; ++i) normal[i] = this->leadingNormal[i];
   }

   int32_t ShockPlanar::getShockRegion(Real t,const Real* pos) {
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 // Calculate shock centroids:
	 calculateShockCentroids(t);
	 
	 if (pos[0] < leadingCentroid[0]) return doubleshock::UPSTREAM;
	 else if (pos[0] < trailingCentroid[0]) return doubleshock::IN_BETWEEN;
	 else return doubleshock::DOWNSTREAM;
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 // Calculate radial coordinate in SIM frame:
	 int32_t I = static_cast<uint32_t>(pos[0]);
	 if (I < 0) {
	    I = 0;
	 } else if (I >= (int32_t)sim->x_blocks*block::WIDTH_X) {
	    I = sim->x_blocks*block::WIDTH_X-1;
	 } else {

	 }
	 Real R_SIM = sim->x_crds_node[I] + (pos[0]-I)*sim->dx_cell[I];
	 
	 // Calculate shock radial coordinates in SIM frame:
	 Real dt = max((Real)0.0,t-simControl.t_setup);
	 Real leadingCentroid = leadingInitialCentroid[0] + leadingVelocity_SIM[0]*dt;
	 Real trailingCentroid = leadingCentroid + trailingRelativeCentroid[0] + relativeVelocity_SIM[0]*t;

	 if (R_SIM > leadingCentroid) return doubleshock::UPSTREAM;
	 else if (R_SIM > trailingCentroid) {
	    return doubleshock::IN_BETWEEN;
	 }
	 else {
	    return doubleshock::DOWNSTREAM;
	 }
      } else {
	 exit(1);
      }
   }

   void ShockPlanar::getShockState(Real t,Real* centroid,Real& radius) {
      radius = numeric_limits<Real>::infinity();
   }

   Real ShockPlanar::getSquaredDistanceToShock(Real t,const Real* RESTRICT pos) {
      // Return signed distance to closer of the two shocks:
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 // Calculate shock centroids:
	 Real dt = max((Real)0.0,t-simControl.t_setup);
	 Real R_leading  = leadingInitialCentroid[0] + (leadingVelocity_SIM[0] - simControl.V_frame[0])*dt;
	 Real R_trailing = R_leading + trailingRelativeCentroid[0] + relativeVelocity_SIM[0]*t;
	 
	 int32_t I = static_cast<uint32_t>(pos[0]);
	 if (I < 0) {
	    I = 0;
	 } else if (I >= (int32_t)sim->x_blocks*block::WIDTH_X) {
	    I = sim->x_blocks*block::WIDTH_X-1;
	 } else {
	    
	 }
	 Real R_SIM = sim->x_crds_node[I] + (pos[0]-I)*sim->dx_cell[I];
	 
	 if (R_SIM < R_leading) {
	    return R_leading - R_SIM;
	 } else if (fabs(R_SIM-R_leading) < fabs(R_trailing-R_SIM)) {
	    return R_leading - R_SIM;
	 } else {
	    return R_trailing - R_SIM;
	 }
	 
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 // Calculate radial coordinate in SIM frame:
	 int32_t I = static_cast<uint32_t>(pos[0]);
	 if (I < 0) {
	    I = 0;
	 } else if (I >= (int32_t)sim->x_blocks*block::WIDTH_X) {
	    I = sim->x_blocks*block::WIDTH_X-1;
	 } else {
	    
	 }
	 Real R_SIM = sim->x_crds_node[I] + (pos[0]-I)*sim->dx_cell[I];

	 // Calculate shock radial coordinates in SIM frame:
	 Real dt = max((Real)0.0,t-simControl.t_setup);
	 Real R_leading  = leadingInitialCentroid[0] + leadingVelocity_SIM[0]*dt;
	 Real R_trailing = R_leading + trailingRelativeCentroid[0] + relativeVelocity_SIM[0]*t;

	 if (R_SIM > R_leading) {
	    return R_SIM - R_leading;
	 } else if (R_leading-R_SIM < R_SIM-R_trailing) {
	    return R_SIM - R_leading;
	 } else {
	    return R_SIM - R_trailing;
	 }
      } else {
	 exit(1);
      }
   }

   bool ShockPlanar::inDownstream(Real t,const Real* RESTRICT pos) {
      cerr << "(SEP SHOCK DOUBLE) ERROR: inDownstream called!" << endl;
      exit(1);
   }

   bool ShockPlanar::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Init base class:
      initialized = ShockBaseClass::initialize(sim,simClasses,cr);

      // Read config file:
      if (addConfigFileItems(cr) == false) initialized = false;
      cr.parse();

      string lengthUnitsString,shockModelString;
      cr.get(PREFIX+".leading_shock.speed",leadingVelocity_SIM[0]);
      cr.get(PREFIX+".leading_shock.centroid_x",leadingInitialCentroid[0]);
      cr.get(PREFIX+".leading_shock.centroid_y",leadingInitialCentroid[1]);
      cr.get(PREFIX+".leading_shock.centroid_z",leadingInitialCentroid[2]);
      
      cr.get(PREFIX+".trailing_shock.speed",trailingVelocity_SIM[0]);
      cr.get(PREFIX+".trailing_shock.centroid_x",trailingInitialCentroid[0]);
      cr.get(PREFIX+".trailing_shock.centroid_y",trailingInitialCentroid[1]);
      cr.get(PREFIX+".trailing_shock.centroid_z",trailingInitialCentroid[2]);
      
      cr.get(PREFIX+".length_units",lengthUnitsString);
      cr.get(PREFIX+".shock_model",shockModelString);

      // Sanity check on input parameters:
      const double lengthUnits = simClasses.constants.getDistanceInSI(lengthUnitsString);
      if (lengthUnits == numeric_limits<double>::infinity()) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid length units" << endl << write;
	 initialized = false;
      }

      if (configreader::checkInputValues<3>(leadingVelocity_SIM,DEF_REAL) == false) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid leading shock velocity vector" << endl;
	 simClasses.logger << "\t values: " << trailingVelocity_SIM[0] << '\t' << trailingVelocity_SIM[1] << '\t';
	 simClasses.logger << trailingVelocity_SIM[2] << endl << write;
	 initialized = false;
      }
      
      if (configreader::checkInputValues<3>(leadingInitialCentroid,DEF_REAL) == false) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid leading shock centroid position" << endl;
	 simClasses.logger << "\t values: " << leadingInitialCentroid[0] << '\t' << leadingInitialCentroid[1] << '\t';
	 simClasses.logger << leadingInitialCentroid[2] << endl << write;
	 initialized = false;
      }
      
      if (configreader::checkInputValues<3>(trailingInitialCentroid,DEF_REAL) == false) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid trailing shock centroid position" << endl;
	 simClasses.logger << "\t values: " << trailingInitialCentroid[0] << '\t' << trailingInitialCentroid[1] << '\t'; 
	 simClasses.logger << trailingInitialCentroid[2] << endl << write;
	 initialized = false;
      }

      if (configreader::checkInputValues<3>(trailingVelocity_SIM,DEF_REAL) == false) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid trailing shock velocity vector" << endl;
	 simClasses.logger << "\t values: " << trailingVelocity_SIM[0] << '\t' << trailingVelocity_SIM[1] << '\t';
	 simClasses.logger << trailingVelocity_SIM[2] << endl << write;
	 initialized = false;
      }

      // Get gas compression ratio model:
      if (shockModelString == "constant") {
	 gasCompressionModel = gascompression::CONSTANT;
	 cr.get(PREFIX+".leading_shock.gas_compression_ratio",leading_R_gas);
	 cr.get(PREFIX+".leading_shock.magnetic_compression_ratio",leading_R_magnetic);
	 cr.get(PREFIX+".trailing_shock.gas_compression_ratio",trailing_R_gas);
	 cr.get(PREFIX+".trailing_shock.magnetic_compression_ratio",trailing_R_magnetic);
      } else if (shockModelString == "oblique") {
	 gasCompressionModel = gascompression::OBLIQUE;
      } else if (shockModelString == "perpendicular") {
	 gasCompressionModel = gascompression::PERPENDICULAR;
      } else {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Invalid shock model '" << shockModelString << "'" << endl;
	 simClasses.logger << "                           Should be 'perpendicular' or 'oblique'" << endl << write;
	 initialized = false;
      }

      // Simulation is run in leading shock's rest frame, and 
      // shock normals are restricted to point to -x direction:
      for (int i=1; i<3; ++i) leadingVelocity_SIM[i] = 0.0;
      for (int i=1; i<3; ++i) trailingVelocity_SIM[i] = 0.0;
      for (int i=0; i<3; ++i) leadingNormal[i] = 0.0;
      for (int i=0; i<3; ++i) trailingNormal[i] = 0.0;
      
      if (simControl.coordinateSystem == sep::CARTESIAN) {
	 leadingNormal[0] = -1.0;
	 trailingNormal[0] = -1.0;
	 trailingVelocity_SIM[0] *= -1.0;
	 simClasses.logger << "(SEP SHOCK DOUBLE) NOTE: Current implementation forces shock normals to point to -x direction," << endl;
	 simClasses.logger << "                         and runs in leading shock's normal incidence frame." << endl << write;
	 
	 // Transform centroid positions into simulation units:
	 for (int i=0; i<3; ++i) leadingInitialCentroid[i] *= lengthUnits;
	 for (int i=0; i<3; ++i) trailingInitialCentroid[i] *= lengthUnits;	 
      } else if (simControl.coordinateSystem == sep::SPHERICAL) {
	 leadingNormal[0] = +1.0;
	 trailingNormal[0] = +1.0;
	 leadingInitialCentroid[0]  *= lengthUnits;
	 trailingInitialCentroid[0] *= lengthUnits;
      }

      for (int i=0; i<3; ++i) trailingRelativeCentroid[i] 
	= trailingInitialCentroid[i] - leadingInitialCentroid[i];
      for (int i=0; i<3; ++i) relativeVelocity_SIM[i] 
	= trailingVelocity_SIM[i] - leadingVelocity_SIM[i];

      // Check that we're using Cartesian coordinates:
      if (simControl.coordinateSystem != sep::CARTESIAN && simControl.coordinateSystem != sep::SPHERICAL) {
	 simClasses.logger << "(SEP SHOCK DOUBLE) ERROR: Current implementation only works in Cartesian and Spherical geometries." << endl << write;
	 initialized = false;
      }

      simClasses.logger << "(SEP SHOCK DOUBLE) Initialization completed, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << endl << write;
      else simClasses.logger << "FAILURE" << endl << write;
      return initialized;
   }
   
   bool ShockPlanar::initializeMesh(Real t,int refinements) {
      return createShockMesh(t,refinements);
   }
   
} // namespace sep

