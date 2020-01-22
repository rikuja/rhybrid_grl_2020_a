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
#include <triangulated_sphere.h>

#include "sep_simcontrol.h"
#include "sep_shock_spherical.h"
#include "sep_mesh_logical.h"
#include "sep_shock_drift_acceleration.h"

using namespace std;

static const Real DEF_REAL = numeric_limits<Real>::infinity();
static const string PREFIX = "Shock.spherical";

namespace sep {

   extern sep::SimControl simControl;

   ShockBaseClass* SphericalShockMaker() {return new ShockSpherical();}
   
   ShockSpherical::ShockSpherical() { 

   }
   
   ShockSpherical::~ShockSpherical() { 
      finalize();
   }
   
   bool ShockSpherical::addConfigFileItems(ConfigReader& cr) {      
      cr.add(PREFIX+".centroid_x","Shock origin at time t=0, x-coordinate (float).",DEF_REAL);
      cr.add(PREFIX+".centroid_y","Shock origin at time t=0, y-coordinate (float).",DEF_REAL);
      cr.add(PREFIX+".centroid_z","Shock origin at time t=0, z-coordinate (float).",DEF_REAL);
      cr.add(PREFIX+".length_units","Units in which positions and radii are given, defaults to m (m/km/AU/RE/RS/RV)",string("m"));
      cr.add(PREFIX+".initial_radius","Shock radius at time t=0 (float).",DEF_REAL);
      cr.add(PREFIX+".expansion_speed","Shock superexpansion speed in m/s (float).",DEF_REAL);
      cr.add(PREFIX+".propagation_speed","Shock propagation speed in m/s, defaults to zero (float).",(Real)0.0);
      cr.add(PREFIX+".propagation_direction_x","Direction of shock propagation, x-component (float).",(Real)1.0);
      cr.add(PREFIX+".propagation_direction_y","Direction of shock propagation, y-component (float).",(Real)0.0);
      cr.add(PREFIX+".propagation_direction_z","Direction of shock propagation, z-component (float).",(Real)0.0);
      cr.add(PREFIX+".shock_model","Which shock model is used, 'oblique' or 'perpendicular' (string).",string(""));
      cr.add(PREFIX+".gas_compression_ratio","Gas compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      cr.add(PREFIX+".magnetic_compression_ratio","Magnetic compression ratio, only used for 'constant' shock model (float).",(Real)1.0);
      cr.add(PREFIX+".cruising_time","Time is takes for CME to reach cruising speed (float).",DEF_REAL);
      return true;
   }

   Real ShockSpherical::calculateGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* RESTRICT pos) {
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
      switch (gasCompressionModel) {
       case gascompression::CONSTANT:
	 return R_gas;
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

   Real ShockSpherical::calculateMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) {
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
      switch (gasCompressionModel) {
       case gascompression::CONSTANT:
	 return R_magnetic;
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

      /*if (R_magn < 0.0) {
	 cerr << "SHOCK: " << alfvenMach2 << '\t' << R_gas << '\t' << R_magn << '\t' << acos(cos_psi1) << endl;
	 cerr << "\t\t" << sonicMach2 << endl;
      }*/

      return max(1.0,fabs(R_magn));
      //return max(1.0,R_magn);
      /*
      const Real B1_mag2 = vectorMagnitude2<3>(B);
      const Real B1_tang_mag2 = B1_mag2 - B_norm_mag*B_norm_mag;
      const Real B2_mag2 = B_norm_mag*B_norm_mag + R_magn*R_magn*B1_tang_mag2;
      return sqrt(B2_mag2 / B1_mag2);
       */
   }

   void ShockSpherical::calculateRelativePosition(const Real* RESTRICT pos,const Real* RESTRICT cartesianCentroid,Real* RESTRICT relPos) {
      const int32_t I = static_cast<int32_t>(pos[0]);
      const int32_t J = static_cast<int32_t>(pos[1]);
      const int32_t K = static_cast<int32_t>(pos[2]);
      const Real X = sim->x_crds_node[I] + (pos[0]-I)*sim->dx_cell[I];
      const Real Y = sim->y_crds_node[J] + (pos[1]-J)*sim->dy_cell[J];
      const Real Z = sim->z_crds_node[K] + (pos[2]-K)*sim->dz_cell[K];
      
      switch (simControl.coordinateSystem) {
       case UNKNOWN:
	 relPos[0] = X - cartesianCentroid[0];
	 relPos[1] = Y - cartesianCentroid[1];
	 relPos[2] = Z - cartesianCentroid[2];
	 simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getSquaredDistanceToShock" << endl << write;
	 exit(1);
	 break;
       case CARTESIAN:
	 relPos[0] = X - cartesianCentroid[0];
	 relPos[1] = Y - cartesianCentroid[1];
	 relPos[2] = Z - cartesianCentroid[2];
	 break;
       case CYLINDRICAL:
	 relPos[0] = X*cos(Y) - cartesianCentroid[0];
	 relPos[1] = X*sin(Y) - cartesianCentroid[1];
	 relPos[2] = Z        - cartesianCentroid[2];
	 break;
       case SPHERICAL:
	 // Transform to Cartesian coordinates:
	 relPos[0] = X*sin(Y)*cos(Z) - cartesianCentroid[0];
	 relPos[1] = X*sin(Y)*sin(Z) - cartesianCentroid[1];
	 relPos[2] = X*cos(Y)        - cartesianCentroid[2];
	 break;
       default:
	 relPos[0] = X - cartesianCentroid[0];
	 relPos[1] = Y - cartesianCentroid[1];
	 relPos[2] = Z - cartesianCentroid[2];
	 simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getSquaredDistanceToShock" << endl << write;
	 exit(1);
	 break;
      }
   }

   void ShockSpherical::calculateShockCentroid(Real t,Real* RESTRICT cartCentroid) {
      // Calculate current distance centroid has travelled
      // to direction of V_cartDirection:
      Real R_centroid = 0.0;
      if (t < simControl.t_shock) {
	 
      } else if (t < t_cruising) {
	 const Real dt = t - simControl.t_shock;
	 R_centroid += accelerationCentroid*dt*dt/2 - accelerationCentroid*dt*dt*dt/(t_cruising - simControl.t_shock)/6;
      } else {
	 const Real dt = t_cruising - simControl.t_shock;
	 R_centroid += accelerationCentroid*dt*dt/3 + accelerationCentroid*dt*(t-t_cruising)/2;
      }

      cartCentroid[0] = centroid[0] + R_centroid*V_cartDirection[0];
      cartCentroid[1] = centroid[1] + R_centroid*V_cartDirection[1];
      cartCentroid[2] = centroid[2] + R_centroid*V_cartDirection[2];
   }
   
   Real ShockSpherical::calculateShockRadius(Real t) {
      // Get shock centroid:
      Real cartCentroid[3];
      calculateShockCentroid(t,cartCentroid);
      Real R_centroid   = vectorMagnitude<3>(cartCentroid);
      Real R_centroid_0 = vectorMagnitude<3>(centroid);

      // Calculate expanded radius (non-geometric expansion only):
      Real expandedRadius = initialRadius;
      if (t < simControl.t_shock) {
	 
      } else if (t < t_cruising) {
	 const Real dt = t - simControl.t_shock;
	 //expandedRadius += V_expansion*dt - V_expansion*dt*dt/(t_cruising-simControl.t_shock)/2;
	 expandedRadius += V_expansion*dt - V_expansion*dt*dt/(t_cruising-simControl.t_shock)
	   + V_expansion*dt*dt*dt/(t_cruising-simControl.t_shock)/(t_cruising-simControl.t_shock)/3;
      } else {
	 const Real dt = t_cruising-simControl.t_shock;
	 //expandedRadius += V_expansion*dt/2;
	 expandedRadius += V_expansion*dt/3;
      }
      
      // Scale radius to take geometric expansion into account:
      return expandedRadius * R_centroid / R_centroid_0;
   }
   
   Real ShockSpherical::calculateSquaredShockRadius(Real t) {
      const Real radius = calculateShockRadius(t);
      return radius*radius;
   }

   bool ShockSpherical::createShockMesh(Real t,int refinements) {
      // Clear previous mesh:
	{
	   vector<Real> dummy1;
	   vector<uint32_t> dummy2;
	   nodeCoordinates.swap(dummy1);
	   cellConnectivity.swap(dummy2);
	}
      
      const Real radius = calculateShockRadius(t);
      size_t N_connections = 0;
      Real currentCentroid[3];
      calculateShockCentroid(t,currentCentroid);
      
      vector<ucdmesh::Node> nodes;
      vector<ucdmesh::TriangularFace> faces;
      const bool result = ucdmesh::createSphere(ucdmesh::spheregenerator::icosahedron,refinements,currentCentroid,
						radius,nodes,faces,N_connections);
      if (result == false) {
	 simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: ucdmesh::createSphere failed" << endl;
	 simClasses->logger << "\t centroid: " << currentCentroid[0] << '\t' << currentCentroid[1] << '\t' << currentCentroid[2] << endl;
	 simClasses->logger << "\t refinements: " << refinements << "\t radius: " << radius << endl;
	 simClasses->logger << write;
	 return result;
      }
      
      // Get face data (areas and unit vectors):
      ucdmesh::getFaceData(nodes,faces,faceData);
      
      // Temporary solution
      #warning FIXME inefficient method of shock mesh generation for plotting
      nodeCoordinates.resize(nodes.size()*3);
      for (size_t i=0; i<nodes.size(); ++i) {
	 nodeCoordinates[3*i+0] = nodes[i].x;
	 nodeCoordinates[3*i+1] = nodes[i].y;
	 nodeCoordinates[3*i+2] = nodes[i].z;
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
      return result;
   }
   
   bool ShockSpherical::finalize() {
      return true;
   }
   
   bool ShockSpherical::finalizeMesh() {
      vector<Real> dummy1;
      vector<uint32_t> dummy2;
      vector<uint32_t> dummy3;
      nodeCoordinates.swap(dummy1);
      cellConnectivity.swap(dummy2);
      cellOffsets.swap(dummy3);
      return true;
   }

   Real ShockSpherical::getAllowedDistance() const {return 10000.0*10000.0;}
   
   vector<uint32_t>& ShockSpherical::getCellConnectivity() {return cellConnectivity;}
   
   vector<uint32_t>& ShockSpherical::getCellOffsets() {return cellOffsets;}

   Real ShockSpherical::getDistanceToShock(Real t,const Real* RESTRICT pos) {
      return sqrt(getSquaredDistanceToShock(t,pos));
   }

   vector<Real>& ShockSpherical::getFaceData(uint64_t& N_surfaces) {
      N_surfaces = faceData.size() / ucdmesh::facedataelement::SIZE;
      return faceData;
   }

   Real ShockSpherical::getGasCompressionRatio(pargrid::CellID blockLID,Real t,const Real* RESTRICT pos) {
      return calculateGasCompressionRatio(blockLID,t,pos);
   }

   void ShockSpherical::getLocalShockVelocity(Real t,const Real* RESTRICT pos,Real* RESTRICT shockVelocity) {
      // Calculate shock centroid position:
      Real cartesianCentroid[3];
      calculateShockCentroid(t,cartesianCentroid);
      Real R_centroid   = vectorMagnitude<3>(cartesianCentroid);
      Real R_centroid_0 = vectorMagnitude<3>(centroid);
      
      // Calculate current lateral expansion speed, 
      // centroid speed, and expanded radius:
      Real radialSpeed = 0.0;
      Real expandedRadius = initialRadius;
      Real centroidSpeed = 0.0;
      if (t < simControl.t_shock) { 
	 
      } else if (t < t_cruising) {
	 const Real dt = t - simControl.t_shock;
	 const Real dt_cruise = t_cruising - simControl.t_shock;
	 //radialSpeed = V_expansion*(1 - dt / dt_cruise);
	 radialSpeed = V_expansion*(1 - 2*dt/dt_cruise + dt*dt/dt_cruise/dt_cruise);
	 centroidSpeed = accelerationCentroid*dt*(1 - dt/dt_cruise/2);
	 expandedRadius += V_expansion*dt*(1 - dt/dt_cruise + dt*dt/dt_cruise/dt_cruise/3);
      } else {
	 const Real dt_cruise = t_cruising-simControl.t_shock;
	 centroidSpeed = accelerationCentroid*dt_cruise/2;
	 expandedRadius += V_expansion*dt_cruise/3;
      }
      radialSpeed     = radialSpeed*(R_centroid/R_centroid_0) + expandedRadius*centroidSpeed/R_centroid_0;

      // Calculate shock normal at given position relative to shock origin:
      Real normal[3];
      getShockNormal(t,pos,cartesianCentroid,normal);
      
      Real V_centroid[3];
      sphericalToSpherical<Real>(pos,V_cartDirection,V_centroid);
      for (int i=0; i<3; ++i) shockVelocity[i] = centroidSpeed*V_centroid[i] + radialSpeed*normal[i];
   }
   
   void ShockSpherical::getLogicalNodeCoordinates(std::vector<Real>& nodeCoords) {
      nodeCoords.clear();
      nodeCoords.resize(nodeCoordinates.size());
      
      Real tmp[3];
      const size_t N_nodes = nodeCoordinates.size()/3;
      for (size_t node=0; node<N_nodes; ++node) {
	 // Transform Cartesian node coordinates into simulation coordinate system:
	 Real r,theta,phi,xy;
	 switch (simControl.coordinateSystem) {
	  case UNKNOWN:
	    simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getLogicalNodeCoordinates" << endl << write;
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
	    simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getLogicalNodeCoordinates" << endl << write;
	    exit(1);
	    break;
	 }

	 getLogicalCoordinates(sim,tmp,&(nodeCoords[3*node+0]));
      }
   }
   
   Real ShockSpherical::getMagneticCompressionRatio(pargrid::CellID blockLID,Real t,const Real* pos) {
      return calculateMagneticCompressionRatio(blockLID,t,pos);
   }

   vector<Real>& ShockSpherical::getNodeCoordinates() {return nodeCoordinates;}
   
   uint32_t ShockSpherical::getNumberOfCells() const {return cellConnectivity.size()/5;}
   
   uint32_t ShockSpherical::getNumberOfNodes() const {return nodeCoordinates.size()/3;}

   void ShockSpherical::getShockNormal(Real t,const Real* RESTRICT pos,Real* RESTRICT normal) {
      Real cartCentroid[3];
      calculateShockCentroid(t,cartCentroid);
      getShockNormal(t,pos,cartCentroid,normal);
   }
   
   void ShockSpherical::getShockNormal(Real t,const Real* RESTRICT pos,const Real* RESTRICT cartesianCentroid,Real* RESTRICT normal) {
      // Transform to simulation frame (if necessary):
      Real tmpNormal[3];
      switch (simControl.coordinateSystem) {
       case UNKNOWN:
	 simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getShockNormal" << endl << write;
	 exit(1);
	 break;
       case CARTESIAN:
	 calculateRelativePosition(pos,cartesianCentroid,normal);
	 unitVector<3>(normal);
	 break;
       case CYLINDRICAL:
	 calculateRelativePosition(pos,cartesianCentroid,tmpNormal);
	 unitVector<3>(tmpNormal);
	 cylindricalToCylindrical<Real>(pos,tmpNormal,normal);
	 break;
       case SPHERICAL:
	 calculateRelativePosition(pos,cartesianCentroid,tmpNormal);
	 unitVector<3>(tmpNormal);
	 sphericalToSpherical<Real>(pos,tmpNormal,normal);
	 break;
       default:
	 simClasses->logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system in getShockNormal" << endl << write;
	 exit(1);
	 break;
      }
	 
      return;
   }

   int32_t ShockSpherical::getShockRegion(Real t,const Real* RESTRICT pos) {
      Real cartCentroid[3];
      calculateShockCentroid(t,cartCentroid);
      Real relPos[3];
      calculateRelativePosition(pos,cartCentroid,relPos);
      if (vectorMagnitude2<3>(relPos) - calculateSquaredShockRadius(t) >= 0) return sphericalshock::UPSTREAM;
      else return sphericalshock::DOWNSTREAM;
   }
   
   void ShockSpherical::getShockState(Real t,Real* centroid,Real& radius) {
      calculateShockCentroid(t,centroid);
      radius = calculateShockRadius(t);
   }

   Real ShockSpherical::getSquaredDistanceToShock(Real t,const Real* RESTRICT pos) {
      Real cartCentroid[3];
      calculateShockCentroid(t,cartCentroid);
      Real relPos[3];
      calculateRelativePosition(pos,cartCentroid,relPos);
      
      Real radius = vectorMagnitude<3>(relPos);
      Real shockRadius = calculateShockRadius(t);
      return radius - shockRadius;
   }

   bool ShockSpherical::inDownstream(Real t,const Real* RESTRICT pos) {
      return getSquaredDistanceToShock(t,pos) < 0.0;
   }
   
   bool ShockSpherical::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      // Init base class:
      initialized = ShockBaseClass::initialize(sim,simClasses,cr);

      // Read config file:
      if (addConfigFileItems(cr) == false) initialized = false;
      cr.parse();
      
      string lengthUnitsString,shockModelString;
      cr.get(PREFIX+".centroid_x",centroid[0]);
      cr.get(PREFIX+".centroid_y",centroid[1]);
      cr.get(PREFIX+".centroid_z",centroid[2]);
      cr.get(PREFIX+".length_units",lengthUnitsString);
      cr.get(PREFIX+".initial_radius",initialRadius);
      cr.get(PREFIX+".expansion_speed",V_expansion);
      cr.get(PREFIX+".propagation_speed",speed);
      cr.get(PREFIX+".propagation_direction_x",V_direction[0]);
      cr.get(PREFIX+".propagation_direction_y",V_direction[1]);      
      cr.get(PREFIX+".propagation_direction_z",V_direction[2]);
      cr.get(PREFIX+".shock_model",shockModelString);
      cr.get(PREFIX+".cruising_time",t_cruising);
      
      // Sanity check on input parameters:
      const double lengthUnits = simClasses.constants.getDistanceInSI(lengthUnitsString);
      if (lengthUnits == numeric_limits<double>::infinity()) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid length units" << endl << write;
	 initialized = false;
      }
      
      if (t_cruising == DEF_REAL) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Cruising time was not given in config file" << endl << write;
	 initialized = false;
      }

      if (centroid[0] == DEF_REAL) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid initial x-position" << endl << write;
	 initialized = false;
      }
      if (centroid[1] == DEF_REAL) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid initial y-position" << endl << write;
	 initialized = false;
      }
      if (centroid[2] == DEF_REAL) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid initial z-position" << endl << write;
	 initialized = false;
      }
      if (V_expansion == DEF_REAL || V_expansion < 0.0) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid radial expansion speed" << endl << write;
	 initialized = false;
      }
      /*if (speed < 0.0) {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Shock propagation speed cannot be negative" << endl << write;
	 initialized = false;
      }*/

      if (shockModelString == "constant") {
	 gasCompressionModel = gascompression::CONSTANT;
	 cr.get(PREFIX+".gas_compression_ratio",R_gas);
	 cr.get(PREFIX+".magnetic_compression_ratio",R_magnetic);
      } else if (shockModelString == "oblique") {
	 gasCompressionModel = gascompression::OBLIQUE;
      } else if (shockModelString == "perpendicular") {
	 gasCompressionModel = gascompression::PERPENDICULAR;
      } else {
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid shock model '" << shockModelString << "'" << endl;
	 simClasses.logger << "                             Should be 'constant', 'perpendicular' or 'oblique'" << endl << write;
	 initialized = false;
      }
      
      initialRadius *= lengthUnits;

      Real r,phi,theta;
      switch (simControl.coordinateSystem) {
       case sep::CARTESIAN:
	 V_cartDirection[0] = V_direction[0];
	 V_cartDirection[1] = V_direction[1];
	 V_cartDirection[2] = V_direction[2];
	 break;
       case sep::CYLINDRICAL:
	 r   = centroid[0];
	 phi = centroid[1];
	 if (r < 0.0) {
	    simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid centroid r-coordinate" << endl << write;
	    initialized = false;
	 }
	 
	 // Transform direction of propagation vector to Cartesian coordinates.
	 // Note that transformation depends on position:
	 V_cartDirection[0] = V_direction[0]*cos(phi) - V_direction[1]*sin(phi);
	 V_cartDirection[1] = V_direction[0]*sin(phi) + V_direction[1]*cos(phi);
	 V_cartDirection[2] = V_direction[2];
	 
	 // Transform centroid position to Cartesian coordinates:
	 centroid[0] = r*cos(phi);
	 centroid[1] = r*sin(phi);
	 break;
       case sep::SPHERICAL:
	 r     = centroid[0];
	 theta = centroid[1];
	 phi   = centroid[2];
	 if (r < 0.0) {
	    simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid centroid r-coordinate" << endl << write;
	    initialized = false;
	 }
	 if (theta < 0.0 || theta > M_PI) {
	    simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Invalid centroid theta-coordinate" << endl << write;
	    initialized = false;
	 }
	 
	 // Transform direction of propagation vector to Cartesian coordinates.
	 // Note that transformation depends on position:
	 V_cartDirection[0] = V_direction[0]*sin(theta)*cos(phi) + V_direction[1]*cos(theta)*cos(phi) - V_direction[2]*sin(phi);
	 V_cartDirection[1] = V_direction[0]*sin(theta)*sin(phi) + V_direction[1]*cos(theta)*sin(phi) + V_direction[2]*cos(phi);
	 V_cartDirection[2] = V_direction[0]*cos(theta) - V_direction[1]*sin(theta);

	 // Transform centroid position to Cartesian coordinates:
	 centroid[0] = r*sin(theta)*cos(phi);
	 centroid[1] = r*sin(theta)*sin(phi);
	 centroid[2] = r*cos(theta);
	 break;
       default:
	 simClasses.logger << "(SEP SHOCK SPHERICAL) ERROR: Unknown coordinate system" << endl << write;
	 initialized = false;
	 break;
      }
      unitVector<3>(V_cartDirection);
      for (int i=0; i<3; ++i) centroid[i] *= lengthUnits;
      
      // Write init status to log file:
      simClasses.logger << "(SEP SHOCK SPHERICAL) Initialization complete, status is ";
      if (initialized == true) simClasses.logger << "SUCCESS" << endl;
      else simClasses.logger << "FAILURE" << endl;

      const Real R0 = vectorMagnitude<3>(centroid);

      const Real expandedRadius = initialRadius + V_expansion*t_cruising/3;
      speed = speed / (1.0 + expandedRadius / R0);
      //speed = speed / (1.0 + initialRadius / R0);
      accelerationCentroid = 2*speed / t_cruising;
      t_cruising += simControl.t_shock;

      // Write values of calculated parameters to log file:
      simClasses.logger << "\t Centroid initial position     : " << centroid[0] << '\t' << centroid[1] << '\t' << centroid[2] << endl;
      simClasses.logger << "\t Expansion speed               : " << initialRadius/R0*speed/1000.0 << " km/s" << endl;
      simClasses.logger << "\t Lateral expansion speed       : " << V_expansion/1000.0 << " km/s" << endl;
      simClasses.logger << "\t Velocity direction of centroid: " << V_cartDirection[0] << ' ' << V_cartDirection[1] << ' ' << V_cartDirection[2] << endl;
      simClasses.logger << "\t Centroid acceleration:        : " << accelerationCentroid/1000.0 << " km/s2" << endl;
      simClasses.logger << "\t Initial lateral acceleration  : " << V_expansion/(t_cruising-simControl.t_shock)/1000.0 << " km/s2" << endl;
      simClasses.logger << write;
      
      if (sim.mpiRank == sim.MASTER_RANK) writeValues();
      
      return initialized;
   }
   
   bool ShockSpherical::initializeMesh(Real t,int refinements) {
      return createShockMesh(t,refinements);
   }
   
   void ShockSpherical::writeValues() {
      int N = 1000;
      Real t_min = simControl.t_shock;
      Real t_max = simControl.t_shock + 2*(t_cruising-simControl.t_shock);
      Real dt = (t_max-t_min)/(N-1);
      
      const Real dt_cruise = t_cruising-simControl.t_shock;

      // Open output file and write header:
      fstream out("test_shock_stats.txt",fstream::out);
      out << "# Column data are:" << endl;
      out << "#1: Time (minutes)" << endl;
      out << "#2: Leading edge position (solar radii)" << endl;
      out << "#3: Leading edge speed (km/s)" << endl;
      out << "#4: Leading edge acceleration (km/s^2)" << endl;
      
      for (int i=0; i<N; ++i) {
	 Real t = t_min + i*dt;
	 
	 Real pos[3];
	 Real radius = calculateShockRadius(t);
	 calculateShockCentroid(t,pos);
	 const Real R_centroid   = vectorMagnitude<3>(pos);
	 const Real R_centroid_0 = vectorMagnitude<3>(centroid);
	 
	 Real centroidSpeed = accelerationCentroid*(t-simControl.t_shock) 
	   - accelerationCentroid*(t-simControl.t_shock)*(t-simControl.t_shock)/dt_cruise/2;
	 
	 Real R_leadingEdge = (radius + R_centroid) / constants::DIST_SOLAR_RADIUS;

	 Real V_leadingEdge = 0.0;
	 Real A_leadingEdge = 0.0;
	 Real centroidAcceleration = 0.0;
	 Real expansionAcceleration = 0.0;
	 Real expansionSpeed = 0.0;
	 if (t <= t_cruising) {
	    centroidAcceleration = accelerationCentroid*(1 - (t-simControl.t_shock)/dt_cruise);
	    
	    expansionAcceleration = -2*V_expansion/dt_cruise;
	    expansionAcceleration *= (1 - (t-simControl.t_shock)/dt_cruise);
	    expansionSpeed = V_expansion - 2*V_expansion*(t-simControl.t_shock)/dt_cruise
	      + V_expansion*(t-simControl.t_shock)*(t-simControl.t_shock)/dt_cruise/dt_cruise;

	    Real expandedRadius = initialRadius 
	      + V_expansion*(t-simControl.t_shock)*(1 - (t-simControl.t_shock)/dt_cruise 
						    + (t-simControl.t_shock)*(t-simControl.t_shock)/dt_cruise/dt_cruise/3);
	    V_leadingEdge = centroidSpeed*(1 + expandedRadius/R_centroid_0) 
	      + expansionSpeed*R_centroid/R_centroid_0;
	    
	    A_leadingEdge = centroidAcceleration + expansionAcceleration*R_centroid/R_centroid_0
	      + 2*expansionSpeed*centroidSpeed/R_centroid_0
	      + radius*centroidAcceleration/R_centroid_0;
	 } else {
	    centroidSpeed = accelerationCentroid*dt_cruise/2;
	    V_leadingEdge = centroidSpeed + (initialRadius + V_expansion*dt_cruise/3)*centroidSpeed/R_centroid_0;
	 }
	 
	 out << i*dt/60.0 << '\t';
	 out << R_leadingEdge << '\t';
	 out << V_leadingEdge/1000.0 << '\t';
	 out << A_leadingEdge/1000.0 << '\t';
	 out << endl;
      }
      
      out.close();
   }
   
} // namespace sep

