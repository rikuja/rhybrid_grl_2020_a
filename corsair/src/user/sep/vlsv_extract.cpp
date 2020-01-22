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

#include <cstdlib>
#include <iostream>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include <vlsv_reader.h>

using namespace std;

uint64_t calculateAcceptedIndex(const double& crd,const double* nodeCoords,const size_t& size);
uint64_t calculateGlobalIndex(uint64_t i,uint64_t j,uint64_t k,uint64_t isize,uint64_t jsize,uint64_t ksize);
void calculateIndices(uint64_t globalID,uint64_t& i,uint64_t& j,uint64_t& k,uint64_t isize,uint64_t jsize,uint64_t ksize);

struct MeshInfo {
   vector<uint64_t> acceptedGlobalIDs;
   vector<uint64_t> boundingBox;
   
   //vector<size_t> globalIDs;
   vector<size_t> indices;
   double* x_crds;
   double* y_crds;
   double* z_crds;
   
   MeshInfo() {
      x_crds = NULL;
      y_crds = NULL;
      z_crds = NULL;
   }
   
   ~MeshInfo() {
      delete [] x_crds;
      delete [] y_crds;
      delete [] z_crds;
   }
};

bool calculateAcceptedCells(vlsv::Reader& vlsvReader,const map<string,string>& arguments,MeshInfo& meshInfo) {
   double xcrd=NAN,ycrd=NAN,zcrd=NAN;
   uint64_t i_min,i_max;
   uint64_t j_min,j_max;
   uint64_t k_min,k_max;

   // Check that x-coordinate was specified in command line:
   map<string,string>::const_iterator it = arguments.find("x");
   if (it == arguments.end()) {
      cerr << "ERROR: x-coordinate was not specified in command line with --x=<value>" << endl;
      return false;
   } else {
      if (it->second == "*") xcrd = numeric_limits<double>::infinity();
      else xcrd = atof(it->second.c_str());
   }
   
   // Check that y-coordinate was specified in command line:
   it = arguments.find("y");
   if (it == arguments.end()) {
      cerr << "ERROR: y-coordinate was not specified in command line with --y=<value>" << endl;
      return false;
   } else {
      if (it->second == "*") ycrd = numeric_limits<double>::infinity();
      else ycrd = atof(it->second.c_str());
   }
   
   // Check that z-coordinate was specified in command line:
   it = arguments.find("z");
   if (it == arguments.end()) {
      cerr << "ERROR: z-coordinate was not specified in command line with --z=<value>" << endl;
      return false;
   } else {
      if (it->second == "*") zcrd = numeric_limits<double>::infinity();
      else zcrd = atof(it->second.c_str());
   }
   
   // Calculate i index limits for accepted cells:
   if (xcrd != numeric_limits<double>::infinity()) {
      i_min = calculateAcceptedIndex(xcrd,meshInfo.x_crds,meshInfo.boundingBox[0]*meshInfo.boundingBox[3]);
      i_max = i_min+1;
   } else {
      i_min = 0;
      i_max = meshInfo.boundingBox[0]*meshInfo.boundingBox[3];
   }

   // Calculate j index limits for accepted cells:
   if (ycrd != numeric_limits<double>::infinity()) {
      j_min = calculateAcceptedIndex(ycrd,meshInfo.y_crds,meshInfo.boundingBox[1]*meshInfo.boundingBox[4]);
      j_max = j_min+1;
   } else {
      j_min = 0;
      j_max = meshInfo.boundingBox[1]*meshInfo.boundingBox[4];
   }
   
   // Calculate k index limits for accepted cells:
   if (zcrd != numeric_limits<double>::infinity()) {
      k_min = calculateAcceptedIndex(zcrd,meshInfo.z_crds,meshInfo.boundingBox[2]*meshInfo.boundingBox[5]);
      k_max = k_min+1;
   } else {
      k_min = 0;
      k_max = meshInfo.boundingBox[2]*meshInfo.boundingBox[5];
   }

   // Check that given x,y,z coordinate values are inside simulation box:
   bool success = true;
   if (i_min == numeric_limits<uint64_t>::max()) success = false;
   if (i_max == numeric_limits<uint64_t>::max()) success = false;
   if (success == false) {
      cerr << "ERROR: Given x-coordinate '" << xcrd << "' is out of bounds, min/max values are: ";
      cerr << meshInfo.x_crds[0] << '\t' << meshInfo.x_crds[meshInfo.boundingBox[0]*meshInfo.boundingBox[3]] << endl;
      return false;
   }
   
   if (j_min == numeric_limits<uint64_t>::max()) success = false; 
   if (j_max == numeric_limits<uint64_t>::max()) success = false;
   if (success == false) {
      cerr << "ERROR: Given y-coordinate '" << ycrd << "' is out of bounds, min/max values are: ";
      cerr << meshInfo.y_crds[0] << '\t' << meshInfo.y_crds[meshInfo.boundingBox[1]*meshInfo.boundingBox[4]] << endl;
      return false;
   }
   
   if (k_min == numeric_limits<uint64_t>::max()) success = false;
   if (k_max == numeric_limits<uint64_t>::max()) success = false;
   if (success == false) {
      cerr << "ERROR: Given z-coordinate '" << zcrd << "' is out of bounds, min/max values are: ";
      cerr << meshInfo.z_crds[0] << '\t' << meshInfo.z_crds[meshInfo.boundingBox[2]*meshInfo.boundingBox[5]] << endl;
      return false;
   }
   
   // Store accepted cell global indices:
   uint64_t i_size = meshInfo.boundingBox[0]*meshInfo.boundingBox[3];
   uint64_t j_size = meshInfo.boundingBox[1]*meshInfo.boundingBox[4];
   uint64_t k_size = meshInfo.boundingBox[2]*meshInfo.boundingBox[5];
   for (uint64_t k=k_min; k<k_max; ++k) for (uint64_t j=j_min; j<j_max; ++j) for (uint64_t i=i_min; i<i_max; ++i) {
      meshInfo.acceptedGlobalIDs.push_back(calculateGlobalIndex(i,j,k,i_size,j_size,k_size));
   }

   // Read mesh info:
   vlsv::datatype::type datatype;
   uint64_t arraySize,vectorSize,byteSize;
   size_t* bbox = NULL;
   list<pair<string,string> > args;
   map<string,string>::const_iterator meshName = arguments.find("mesh");
   args.push_back(make_pair(string("name"),meshName->second));
   
   if (vlsvReader.getArrayInfo("MESH",args,arraySize,vectorSize,datatype,byteSize) == false) {
      cerr << "ERROR: Failed to read mesh info" << endl;
      return false;
   }
   
   // Read cell global IDs:
   uint64_t* globalIDs;
   if (vlsvReader.read("MESH",args,0,arraySize,globalIDs,true) == false) {
      cerr << "ERROR: Failed to read cell global indices" << endl;
      delete [] globalIDs;
      return false;
   }
   
   sort(meshInfo.acceptedGlobalIDs.begin(),meshInfo.acceptedGlobalIDs.end());
   meshInfo.indices.resize(meshInfo.acceptedGlobalIDs.size());
   vector<uint64_t>::iterator iter;
   uint64_t counter = 0;
   for (uint64_t i=0; i<arraySize; ++i) {
      iter = find(meshInfo.acceptedGlobalIDs.begin(),meshInfo.acceptedGlobalIDs.end(),globalIDs[i]);
      if (iter != meshInfo.acceptedGlobalIDs.end()) {
	 size_t position = &(*iter) - &(meshInfo.acceptedGlobalIDs[0]);
	 meshInfo.indices[position] = i;
      }
      ++counter;
   }

   delete [] globalIDs;
   return true;
}

uint64_t calculateAcceptedIndex(const double& crd,const double* nodeCoords,const size_t& size) {
   for (size_t i=0; i<size; ++i) {
      if (crd >= nodeCoords[i] && crd <= nodeCoords[i+1]) return i;
   }
   return numeric_limits<uint64_t>::max();
}

uint64_t calculateGlobalIndex(uint64_t i,uint64_t j,uint64_t k,uint64_t isize,uint64_t jsize,uint64_t ksize) {
   return k*jsize*isize + j*isize + i;
}

void calculateIndices(uint64_t globalID,uint64_t& i,uint64_t& j,uint64_t& k,uint64_t isize,uint64_t jsize,uint64_t ksize) {
   k = globalID / (isize*jsize);
   globalID -= k*isize*jsize;
   j = globalID / isize;
   i = globalID - j*isize;
}

bool parseArguments(int argn,char** args,map<string,string>& arguments) {
   
   for (int i=1; i<argn; ++i) {
      // Each argument-value pair starts with "--"
      if (args[i][0] != '-') continue;
      if (args[i][1] != '-') continue;
      
      string s(args[i]);
      s = s.substr(2,string::npos);
      
      // Check that argument has a value:
      if (s.find_first_of("=",0) == string::npos) continue;
      
      // Parse argument name and check that its size is nonzero:
      string argument = s.substr(0,s.find_first_of("=",0));
      if (argument.size() == 0) continue;
      
      // Parse value and add argument,value pair to map:
      string value = s.substr(s.find_first_of("=",0)+1,string::npos);
      arguments[argument] = value;
   }
   
   return true;
}

bool readMesh(vlsv::Reader& vlsvReader,const map<string,string>& arguments,MeshInfo& meshInfo) {
   // Read names of meshes in file:
   set<string> meshNames;
   if (vlsvReader.getUniqueAttributeValues("MESH","name",meshNames) == false) return false;

   // Check that mesh name was given in command line:
   map<string,string>::const_iterator meshName = arguments.find("mesh");
   if (meshName == arguments.end()) {
      cerr << "ERROR: Mesh name was not specified in command line with --mesh=<name>" << endl;
      cerr << "       File contains following meshes" << endl;
      for (set<string>::iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
	 cout << '\t' << *it << endl;
      }
      return false;
   }

   // Check that given mesh exists in file:
   if (meshNames.find(meshName->second) == meshNames.end()) {
      cerr << "ERROR: File does not contain a mesh called '" << meshName->second << "'" << endl;
      cerr << "       Existing meshes are:" << endl;
      for (set<string>::iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
	 cout << '\t' << *it << endl;
      }
      return false;
   }

   // Read mesh bounding box:
   vlsv::datatype::type datatype;
   uint64_t arraySize,vectorSize,byteSize;
   size_t* bbox = NULL;
   list<pair<string,string> > args;
   args.push_back(make_pair(string("mesh"),meshName->second));
   
   if (vlsvReader.getArrayInfo("MESH_BBOX",args,arraySize,vectorSize,datatype,byteSize) == false) {
      cerr << "ERROR: Failed to read bounding box info" << endl;
      return false;
   }
   
   if (vlsvReader.read("MESH_BBOX",args,0,arraySize,bbox,true) == false) {
      cerr << "ERROR: Failed to read mesh bounding box" << endl;
      delete [] bbox;
      return false;
   }
   
   for (uint64_t i=0; i<arraySize; ++i) meshInfo.boundingBox.push_back(bbox[i]);
   delete [] bbox;
   
   // Read node coordinates:
   if (vlsvReader.read("MESH_NODE_CRDS_X",args,0,meshInfo.boundingBox[0]*meshInfo.boundingBox[3]+1,meshInfo.x_crds,true) == false) {
      cerr << "ERROR: Failed to read x node coordinates" << endl;
      return false;
   }
   if (vlsvReader.read("MESH_NODE_CRDS_Y",args,0,meshInfo.boundingBox[1]*meshInfo.boundingBox[4]+1,meshInfo.y_crds,true) == false) {
      cerr << "ERROR: Failed to read y node coordinates" << endl;
      return false;
   }
   if (vlsvReader.read("MESH_NODE_CRDS_Z",args,0,meshInfo.boundingBox[2]*meshInfo.boundingBox[5]+1,meshInfo.z_crds,true) == false) {
      cerr << "ERROR: Failed to read z node coordinates" << endl;
      return false;
   }

   return true;
}

bool writeData(vlsv::Reader& vlsvReader,const map<string,string>& arguments,MeshInfo& meshInfo) {
   bool success = true;
   
   // Check that a variable name was given in command line:
   map<string,string>::const_iterator it = arguments.find("variable");
   if (it == arguments.end()) {
      cerr << endl << "ERROR: Variable name was not specified on command line with --variable=<name>" << endl;
      success = false;
   }
   string variableName;
   if (it != arguments.end()) variableName = it->second;
   
   // Check that variable exists in mesh:
   it = arguments.find("mesh");
   string meshName = it->second;
   list<pair<string,string> > attribsIn;
   map<string,string> attribsOut;
   attribsIn.push_back(make_pair("mesh",meshName));
   attribsIn.push_back(make_pair("name",variableName));
   if (vlsvReader.getArrayAttributes("VARIABLE",attribsIn,attribsOut) == false) {
      cerr << "Existing variable names are:" << endl;

      set<string> output;
      vlsvReader.getUniqueAttributeValues("VARIABLE","name",output);

      for (set<string>::iterator it=output.begin(); it!=output.end(); ++it) {
	 attribsIn.clear();
	 attribsOut.clear();
	 attribsIn.push_back(make_pair("name",*it));
	 vlsvReader.getArrayAttributes("VARIABLE",attribsIn,attribsOut);

	 if (attribsOut.find("mesh") != attribsOut.end()) {
	    if (attribsOut["mesh"] == meshName) 
	      cerr << '\t' << *it << endl;
	 }
      }
      cerr << endl;
      success = false;
   }
   if (success == false) return success;
   
   attribsIn.clear();
   attribsIn.push_back(make_pair("mesh",meshName));
   attribsIn.push_back(make_pair("name",variableName));
   double value;
   double* ptr = &value;
   uint64_t i_index,j_index,k_index;
   uint64_t i_size = meshInfo.boundingBox[0]*meshInfo.boundingBox[3];
   uint64_t j_size = meshInfo.boundingBox[1]*meshInfo.boundingBox[4];
   uint64_t k_size = meshInfo.boundingBox[2]*meshInfo.boundingBox[5];

   cout.precision(14);
   for (uint64_t i=0; i<meshInfo.indices.size(); ++i) {
      vlsvReader.read("VARIABLE",attribsIn,meshInfo.indices[i],1,ptr,false);
      calculateIndices(meshInfo.acceptedGlobalIDs[i],i_index,j_index,k_index,i_size,j_size,k_size);
      
      cout << meshInfo.x_crds[i_index] << '\t';
      cout << meshInfo.y_crds[j_index] << '\t';
      cout << meshInfo.z_crds[k_index] << '\t';      
      cout << value << '\t';
      cout << endl;
   }

   return true;
}

int main(int argn,char* args[]) {
   // Parse argument,value pairs from command line:
   map<string,string> arguments;
   if (parseArguments(argn,args,arguments) == false) return 1;

   // Open VLSV file for reading:
   vlsv::Reader vlsvReader;
   if (arguments.find("file") == arguments.end()) {
      cerr << "ERROR: Input file name was not given with --file=<file name>" << endl;
      return 1;
   }
   if (vlsvReader.open(arguments["file"]) == false) {
      cerr << "ERROR: Failed to open file '" << arguments["file"] << "' for reading" << endl;
      return 1;
   }

   MeshInfo meshInfo;
   if (readMesh(vlsvReader,arguments,meshInfo) == false) return 1;
   
   if (calculateAcceptedCells(vlsvReader,arguments,meshInfo) == false) return 1;

   if (writeData(vlsvReader,arguments,meshInfo) == false) return 1;
   
   vlsvReader.close();
   return 0;
}




