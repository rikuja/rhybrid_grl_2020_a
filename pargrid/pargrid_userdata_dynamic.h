/** This file is part of ParGrid parallel grid.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PARGRID_USERDATA_DYNAMIC_H
#define PARGRID_USERDATA_DYNAMIC_H

#include <set>
#include <mpi.h>

#include "pargrid_definitions.h"

namespace pargrid {

   typedef unsigned int ArraySizetype; /**< Datatype used for array sizes and capacities.*/
   
   // ********************************************** //
   // ***** CLASS USERDATADYNAMIC DECLARATIONS ***** //
   // ********************************************** //
   
   template<class PARGRID>
   class UserDataDynamic {
    public:
      UserDataDynamic();
      UserDataDynamic(const UserDataDynamic& udd);
      ~UserDataDynamic();
      
      void finalize();
      char** getArrayPointer() const;
      ArraySizetype* getCapacityPointer() const;
      const std::string& getDatatype() const;
      void getDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype);
      ArraySizetype getElementByteSize() const;
      const std::string& getName() const;
      CellID getNumberOfCells() const;
      void getSizeDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype);
      ArraySizetype* getSizePointer() const;
      void initialize(PARGRID* pargrid,const std::string& name,CellID N_cells,unsigned int elementByteSize,const std::string& datatype);
      void reallocate(CellID start,CellID end);
      void swap(UserDataDynamic& udd);
      
    private:
      char** arrays;                 /**< Arrays containing data for each cell.*/
      MPI_Datatype basicDatatype;    /**< Derived MPI datatype that transfers a single element.*/
      ArraySizetype* capacities;     /**< Array capacities for each cell, in number of elements.*/
      std::string datatype;          /**< Basic datatype stored in the array, either "int", "uint", or "float".*/
      ArraySizetype elementByteSize; /**< Byte size of array element.*/
      std::string name;              /**< Name of the user data array.*/
      CellID N_cells;                /**< Number of cells.*/
      PARGRID* pargrid;              /**< Pointer to ParGrid.*/
      ArraySizetype* sizes;          /**< Current amount of elements in each cell.*/
   };
   
   // ********************************************* //
   // ***** CLASS USERDATADYNAMIC DEFINITIONS ***** //
   // ********************************************* //
   
   template<class PARGRID> inline
   UserDataDynamic<PARGRID>::UserDataDynamic() {
      arrays = NULL;
      basicDatatype = MPI_DATATYPE_NULL;
      capacities = NULL;
      N_cells = 0;
      sizes = NULL;
      finalize();
   }
   
   template<class PARGRID> inline
   UserDataDynamic<PARGRID>::UserDataDynamic(const UserDataDynamic& udd) {
      finalize();
      elementByteSize = udd.elementByteSize;
      name            = udd.name;
      N_cells         = udd.N_cells;
      pargrid         = udd.pargrid;
      datatype        = udd.datatype;
      
      // Copy contents of array sizes:
      sizes = new ArraySizetype[N_cells];
      for (CellID cell=0; cell<N_cells; ++cell) sizes[cell] = udd.sizes[cell];

      capacities = new ArraySizetype[N_cells];
      for (CellID cell=0; cell<N_cells; ++cell) capacities[cell] = udd.capacities[cell];
      
      // Allocate data arrays and copy contents:
      arrays = new char* [N_cells];
      for (CellID cell=0; cell<N_cells; ++cell) {
	      arrays[cell] = new char[elementByteSize*capacities[cell]];
	      for (ArraySizetype i=0; i<elementByteSize*sizes[cell]; ++i)
	      arrays[cell][i] = udd.arrays[cell][i];
      }
      
      // Create a derived MPI datatype that transfers a single array element:
      MPI_Type_contiguous(elementByteSize,MPI_Type<char>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);
   }
   
   template<class PARGRID> inline
   UserDataDynamic<PARGRID>::~UserDataDynamic() {finalize();}
   
   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::finalize() {
      for (DataID cell=0; cell<N_cells; ++cell) {
         delete [] arrays[cell];
	      arrays[cell] = NULL;
      }
      delete [] arrays; arrays = NULL;
      delete [] capacities; capacities = NULL;
      delete [] sizes; sizes = NULL;

      datatype = "";
      elementByteSize = 0;
      name = "";
      N_cells = 0;
      pargrid = NULL;
      
      if (basicDatatype != MPI_DATATYPE_NULL) MPI_Type_free(&basicDatatype);
   }

   template<class PARGRID> inline
   char** UserDataDynamic<PARGRID>::getArrayPointer() const {return arrays;}
   
   template<class PARGRID> inline
   ArraySizetype* UserDataDynamic<PARGRID>::getCapacityPointer() const {return capacities;}

   template<class PARGRID> inline
   const std::string& UserDataDynamic<PARGRID>::getDatatype() const {return datatype;}
   
   /** Get a derived MPI datatype that sends or receives given cell data
    * with a single transfer. Note that array sizes must have been transferred, 
    * i.e. the amount of data sent or received, and reallocate function called, 
    * before getDatatype is called. Otherwise this function may reference 
    * invalid memory areas. This function commits the newly created datatype.
    * @param globalIDs List of transferred cells.
    * @param datatype Variable in which the derived datatype is to be written.*/   
   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::getDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype) {
      // Calculate displacements relative to MPI_BOTTOM (in bytes). 
      // Array blockLengths needs to be 'int' array here because of MPI:
      int* blockLengths           = new int[globalIDs.size()];
      MPI_Aint* displacements     = new MPI_Aint[globalIDs.size()];
      CellID counter = 0;
      for (std::set<CellID>::const_iterator globalID=globalIDs.begin(); globalID!=globalIDs.end(); ++globalID) {
	      const CellID localID = pargrid->getLocalID(*globalID);
	      MPI_Get_address(arrays[localID],displacements+counter);
	      blockLengths[counter] = sizes[localID];
	      ++counter;
      }
      // Create a derived MPI datatype for sending or receiving all requested
      // data with a single send and commit the created datatype:
      MPI_Type_hindexed(globalIDs.size(),blockLengths,displacements,basicDatatype,&datatype);
      MPI_Type_commit(&datatype);
      delete [] blockLengths; blockLengths = NULL;
      delete [] displacements; displacements = NULL;
   }

   template<class PARGRID> inline
   ArraySizetype UserDataDynamic<PARGRID>::getElementByteSize() const {return elementByteSize;}
   
   template<class PARGRID> inline
   const std::string& UserDataDynamic<PARGRID>::getName() const {return name;}
   
   template<class PARGRID> inline
   CellID UserDataDynamic<PARGRID>::getNumberOfCells() const {return N_cells;}
   
   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::getSizeDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype) {
      // Calculate displacements relative to start of array sizes:
      int* displacements = new int[globalIDs.size()];
      CellID counter = 0;
      for (std::set<CellID>::const_iterator globalID=globalIDs.begin(); globalID!=globalIDs.end(); ++globalID) {
	      const CellID localID = pargrid->getLocalID(*globalID);
	      displacements[counter] = localID;
	      ++counter;
      }
      
      // Create a derived datatype for transferring all array sizes
      // with a single send and commit the created datatype:
      MPI_Type_create_indexed_block(globalIDs.size(),1,displacements,MPI_Type<ArraySizetype>(),&datatype);
      MPI_Type_commit(&datatype);
      delete [] displacements; displacements = NULL;
   }

   template<class PARGRID> inline
   ArraySizetype* UserDataDynamic<PARGRID>::getSizePointer() const {return sizes;}
   
   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::initialize(PARGRID* pargrid,const std::string& name,CellID N_cells,
					                              unsigned int elementByteSize,const std::string& datatype) {
      this->pargrid         = pargrid;
      this->name            = name;
      this->N_cells         = N_cells;
      this->elementByteSize = elementByteSize;
      this->datatype        = datatype;
      
      sizes = new ArraySizetype[N_cells];
      for (CellID i=0; i<N_cells; ++i) sizes[i]  = 0;
      
      capacities = new ArraySizetype[N_cells];
      for (CellID i=0; i<N_cells; ++i) capacities[i] = sizes[i];
      
      arrays = new char* [N_cells];
      for (CellID i=0; i<N_cells; ++i) arrays[i] = NULL;
      
      MPI_Type_contiguous(elementByteSize,MPI_Type<char>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);
   }

   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::reallocate(CellID start,CellID end) {
      // Iterate over the given range of cells and check that 
      // array capacities are at least as large as the requested sizes:
      for (CellID cell=start; cell<end; ++cell) {
	      if (sizes[cell] > capacities[cell]) {
	         delete [] arrays[cell];
	         capacities[cell] = sizes[cell];
	         arrays[cell] = new char[elementByteSize*capacities[cell]];
	      }
      }
   }

   /** Swap the contents of two dynamic user data arrays.
    * @param udd Dynamic user data array whose contents are to be swapped with.*/
   template<class PARGRID> inline
   void UserDataDynamic<PARGRID>::swap(UserDataDynamic& udd) {
      // Make dummy copies of member variables:
      char** dummyArrays                       = arrays;
      ArraySizetype* dummyCapacities           = capacities;
      const std::string dummyDatatype          = datatype;
      const ArraySizetype dummyElementByteSize = elementByteSize;
      const std::string dummyName              = name;
      const CellID dummy_N_cells               = N_cells;
      ArraySizetype* dummySizes                = sizes;

      // Set the contents of this container:
      arrays = udd.arrays;
      capacities = udd.capacities;
      datatype = udd.datatype;
      elementByteSize = udd.elementByteSize;
      name = udd.name;
      N_cells = udd.N_cells;
      sizes = udd.sizes;
      MPI_Type_free(&basicDatatype);
      MPI_Type_contiguous(elementByteSize,MPI_Type<char>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);

      // Set the contents of the other container:
      udd.arrays = dummyArrays;
      udd.capacities = dummyCapacities;
      udd.datatype = dummyDatatype;
      udd.elementByteSize = dummyElementByteSize;
      udd.name = dummyName;
      udd.N_cells = dummy_N_cells;
      udd.sizes = dummySizes;
      MPI_Type_free(&udd.basicDatatype);
      MPI_Type_contiguous(udd.elementByteSize,MPI_Type<char>(),&udd.basicDatatype);
      MPI_Type_commit(&udd.basicDatatype);
   }
   
} // namespace pargrid

#endif
