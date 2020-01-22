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

#include <cstdlib>
#include <iostream>

#include "sep_simcontrol.h"
#include "sep_wavelength_mesh_builder.h"
#include "sep_mesh_logical.h"

using namespace std;

namespace sep {

   extern sep::SimControl simControl;
   static const string PREFIX = "WavelengthMesh";
   static const Real DEF_VALUE = numeric_limits<Real>::infinity();
   static Real cellSize;
   
   namespace wmesh {
      namespace linear {
	 Real getLogicalWavelength(const Real& wavelength) {
	    return (wavelength - simControl.wavelengthMeshNodeCoordinates[0]) / cellSize;
	 }
	 
	 Real getPhysicalWavelength(const Real& l) {
	    return simControl.wavelengthMeshNodeCoordinates[0] + l*cellSize;
	 }

	 Real getWavelengthScaleFactor(const Real& wavelength) {
	    return 1.0 / cellSize;
	 }
	 
      } // namespace linear

      namespace logarithmic {
	 
	 /** Get logical wavelength.
	  * @param wavelength Wavelength in physical units.
	  * @return Wavelength in logical units.*/
	 Real getLogicalWavelength(const Real& wavelength) {
	    const uint32_t N = simControl.N_wavelengthMeshCells;

	    if (wavelength < 0) {
	       if (wavelength > -simControl.lambda0) {
		  // Linear mesh in first cell:
		  return N/2 + wavelength/simControl.lambda0;
	       } else {
		  // Logarithmic mesh elsewhere:
		  return N/2-1 - log(-wavelength/simControl.lambda0) / cellSize;
	       }
	    } else {
	       if (wavelength < simControl.lambda0) {
		  // Linear mesh in first cell:
		  return N/2 + wavelength/simControl.lambda0;
	       } else {
		  // Logarithmic mesh elsewhere:
		  return N/2+1 + log(wavelength/simControl.lambda0) / cellSize;
	       }
	    }
	 }

	 /** Get wavelength in physical units.
	  * @param l Wavelength in logical units.
	  * @return Wavelength in physical units.*/
	 Real getPhysicalWavelength(const Real& l) {
	    const uint32_t N = simControl.N_wavelengthMeshCells;
	    
	    if (l < N/2) {
	       if (l >= N/2-1) {
		  // Linear mesh in first cell:
		  return simControl.lambda0 * (l-N/2);
	       } else {
		  // Logarithmic mesh elsewhere:
		  return -simControl.lambda0 * exp((N/2-1-l)*cellSize);
	       }
	    } else {
	       if (l < N/2+1) {
		  // Linear mesh in first cell:
		  return simControl.lambda0 * (l-N/2);
	       } else {
		  // Logarithmic mesh elsewhere:
		  return simControl.lambda0 * exp((l-N/2-1)*cellSize);
	       }
	    }
	 }

	 Real getWavelengthScaleFactor(const Real& wavelength) {
	    const Real d_log_lambda = cellSize;
	    return 1.0 / (wavelength*d_log_lambda);
	 }
	 
      } // namespace logarithmic
      
   } // namespace wmesh
   
   /** Add required configuration file parameters to configuration file reader.
    * @param cr Configuration file reader.*/
   void addConfigParameters(ConfigReader& cr) {
      cr.add(PREFIX+".lambda_min","Min. wavelength in mesh (float).",DEF_VALUE);
      cr.add(PREFIX+".lambda_max","Max. wavelength in mesh (float).",DEF_VALUE);
      cr.add(PREFIX+".size","Number of wavelength cells (int).",(uint32_t)0);
      cr.add(PREFIX+".mesh_type","Mesh type 'linear/logarithmic' (string).",string("linear"));
   }
   
   bool createLogarithmicMesh(Real minLambda,Real maxLambda) {
      bool success = true;

      //cellSize = 2*(log(maxLambda / simControl.lambda0)) / simControl.N_wavelengthMeshCells;
      cellSize = log(maxLambda/simControl.lambda0) / (simControl.N_wavelengthMeshCells/2-1);
      simControl.logicalWavelengthCellSize = cellSize;
      
      const uint32_t NWL = simControl.N_wavelengthMeshCells;
      simControl.wavelengthMeshNodeCoordinates = new Real[NWL+1];

      minLambda = -log(fabs(minLambda));
      maxLambda = +log(fabs(maxLambda));

      for (uint32_t i=0; i<simControl.N_wavelengthMeshCells/2; ++i) {
	 simControl.wavelengthMeshNodeCoordinates[i] = minLambda + i*cellSize;
	 simControl.wavelengthMeshNodeCoordinates[i] =
	   -exp(-simControl.wavelengthMeshNodeCoordinates[i]);
      }
      simControl.wavelengthMeshNodeCoordinates[NWL/2] = 0.0;
      for (uint32_t i=simControl.N_wavelengthMeshCells/2+1; i<simControl.N_wavelengthMeshCells+1; ++i) {
	 simControl.wavelengthMeshNodeCoordinates[i] = log(simControl.lambda0) + (i - NWL/2-1)*cellSize;
	 simControl.wavelengthMeshNodeCoordinates[i] =
	   exp(simControl.wavelengthMeshNodeCoordinates[i]);
      }
      
      simControl.wavelengthMeshCellSizes = new Real[simControl.N_wavelengthMeshCells];
      for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	 simControl.wavelengthMeshCellSizes[l] = simControl.wavelengthMeshNodeCoordinates[l+1]
	   - simControl.wavelengthMeshNodeCoordinates[l];
      }

      simControl.getLogicalWavelength = wmesh::logarithmic::getLogicalWavelength;
      simControl.getPhysicalWavelength = wmesh::logarithmic::getPhysicalWavelength;
      simControl.getWavelengthScaleFactor = wmesh::logarithmic::getWavelengthScaleFactor;

      return success;
   }
   
   /** Create wavelength mesh with uniform cell spacing.
    * @param minLambda Minimum value for wavelength.
    * @param maxLambda Maximum value for wavelength.
    * @return If true, mesh was created successfully.*/
   bool createUniformMesh(Real minLambda,Real maxLambda) {
      bool success = true;

      cellSize = (maxLambda-minLambda) / simControl.N_wavelengthMeshCells;
      simControl.logicalWavelengthCellSize = cellSize;
      
      simControl.wavelengthMeshNodeCoordinates = new Real[simControl.N_wavelengthMeshCells+1];
      for (uint32_t i=0; i<simControl.N_wavelengthMeshCells+1; ++i)
	simControl.wavelengthMeshNodeCoordinates[i] = minLambda + i*cellSize;
      simControl.wavelengthMeshNodeCoordinates[simControl.N_wavelengthMeshCells/2] = 0.0;
      
      simControl.wavelengthMeshCellSizes = new Real[simControl.N_wavelengthMeshCells];
      for (uint32_t l=0; l<simControl.N_wavelengthMeshCells; ++l) {
	 simControl.wavelengthMeshCellSizes[l] = simControl.wavelengthMeshNodeCoordinates[l+1]
	   - simControl.wavelengthMeshNodeCoordinates[l];
      }

      simControl.getLogicalWavelength = wmesh::linear::getLogicalWavelength;
      simControl.getPhysicalWavelength = wmesh::linear::getPhysicalWavelength;
      simControl.getWavelengthScaleFactor = wmesh::linear::getWavelengthScaleFactor;
      
      return success;
   }
   
   /** Read the value of configuration file parameter 
    * PREFIX+".cell_spacing" and return the corresponding enumerated
    * value defined in sep::wavecellspacing.
    * @param cr Configuration file reader.
    * @return Method which is used to create wavelength mesh, or 
    * sep::wavecellspacing::UNKNOWN if configuration file 
    * had erroneous method.*/
   sep::mesh::Type parseCellSpacingMethod(ConfigReader& cr) {
      string cellSpacingString;
      cr.get(PREFIX+".mesh_type",cellSpacingString);
      
      if (cellSpacingString == "linear") return mesh::Type::LINEAR;
      else if (cellSpacingString == "logarithmic") return mesh::Type::LOGARITHMIC;
      else return mesh::Type::UNKNOWN;
   }
   
   /** Read all required parameters from configuration file
    * and create the wavelength mesh. Specifically, this method
    * creates the array sep::simControl.wavelengthMeshNodeCoordinates, which
    * contains wavelength mesh node coordinates. Its
    * size is equal to sep::simControl.N_wavelengthMeshCells+1.
    * @param sim Generic simulation variables.
    * @param simClasses Generic simulation classes.
    * @param cr Configuration file reader.
    * @return If true, wavelength mesh was created correctly.*/
   bool buildWaveMesh(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr) {
      bool success = true;
      
      // Define required config file parameters:
      addConfigParameters(cr);
      
      // Parse values from config file:
      cr.parse();      
      Real minLambda, maxLambda;
      string cellSpacingString;
      cr.get(PREFIX+".lambda_min",minLambda);
      cr.get(PREFIX+".lambda_max",maxLambda);
      cr.get(PREFIX+".size",simControl.N_wavelengthMeshCells);
      simControl.wavelengthMeshType = parseCellSpacingMethod(cr);
      
      // Check config parameters for errors:
      if (minLambda == DEF_VALUE) {
	 simClasses.logger << "(SEP WAVEMESH BUILDER) ERROR: Parameter '";
	 simClasses.logger << PREFIX+".lambda_min' was not found." << endl << write;
	 success = false;
      }
      if (maxLambda == DEF_VALUE) {
	 simClasses.logger << "(SEP WAVEMESH BUILDER) ERROR: Parameter '"; 
	 simClasses.logger << PREFIX+".lambda_max' was not found." << endl << write;
	 success = false;
      }
      if (simControl.N_wavelengthMeshCells == 0) {
	 simClasses.logger << "(SEP WAVEMESH BUILDER) ERROR: Parameter '";
	 simClasses.logger << PREFIX+".size' was not found or has zero value." << endl << write;
	 success = false;
      }
      if (simControl.N_wavelengthMeshCells % 2 != 0) {
	 simClasses.logger << "(SEP WAVEMESH BUILDER) WARNING: Wavelength mesh has odd number of cells, increasing by one to make it even." << endl << write;
	 simControl.N_wavelengthMeshCells += 1;
      }
      if (simControl.wavelengthMeshType == sep::mesh::Type::UNKNOWN) {
	 simClasses.logger << "(SEP WAVEMESH BUILDER) ERROR: Parameter '";
	 simClasses.logger << PREFIX+".mesh_type' has unsupported value." << endl << write;
	 success = false;
      }
      
      // Exit if errors were detected:
      if (success == false) return success;
      
      // Create wavelength mesh with desired cell size method:
      switch (simControl.wavelengthMeshType) {
       case (sep::mesh::Type::UNKNOWN):
	 success = false;
	 break;
       case (sep::mesh::Type::LINEAR):
	 if (createUniformMesh(minLambda,maxLambda) == false) success = false;
	 break;
       case (sep::mesh::Type::LOGARITHMIC):
	 if (createLogarithmicMesh(minLambda,maxLambda) == false) success = false;
	 break;
       default:
	 success = false;
	 break;
      }
      
      return success;
   }

   uint32_t waveIndex(uint32_t i,uint32_t j,uint32_t k,uint32_t l) {
      return block::index(i,j,k)*simControl.N_wavelengthMeshCells + l;
   }
   
} // namespace sep
