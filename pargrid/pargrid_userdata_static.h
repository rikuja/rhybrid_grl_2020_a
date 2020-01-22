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

#ifndef PARGRID_USERDATA_STATIC_H
#define PARGRID_USERDATA_STATIC_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <mpi.h>
#include "pargrid_definitions.h"

namespace pargrid {
   // ********************************************* //
   // ***** CLASS STATICUSERDATA DECLARATIONS ***** //
   // ********************************************* //
   
   /** Wrapper for a user-defined parallel data array in ParGrid.
    * The purpose of this struct is to include all information that is 
    * needed to create, delete, and copy parallel arrays. Some helper 
    * functions are also provided.*/
   template<class PARGRID> 
   struct UserDataStatic {
    public:
      UserDataStatic();
      UserDataStatic(const UserDataStatic& udw);
      ~UserDataStatic();
      
      void copy(const UserDataStatic& udw,CellID newCell,CellID oldCell);
      bool finalize();
      void* getAddress();
      void getDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype);
      unsigned int getElementSize() const {return N_elements*byteSize;}
      const std::string& getName() const;
      void initialize(PARGRID* pargrid,const std::string& name,size_t N_cells,unsigned int elements,
		      unsigned int byteSize,const std::string& datatype);
      void swap(UserDataStatic& uds);
      
      char* array;                   /**< Pointer to array containing the data.*/
      std::string name;              /**< Name of the array.*/
      unsigned int byteSize;         /**< Byte size of a single array element, i.e. sizeof(double) for a double array.*/
      std::string datatype;          /**< Datatype stored in array, "int", "uint", or "float".*/
      unsigned int N_cells;          /**< Number of cells (=elements) in the array, equal to the number of 
				      *                                       * cells (local+remote) on this process.*/
      unsigned int N_elements;       /**< How many values are reserved per cell.*/
      PARGRID* pargrid;              /**< Pointer to ParGrid class that owns this UserDataStatic.*/
    private:
      MPI_Datatype basicDatatype;    /**< MPI derived datatype that transfers data on a single cell.*/
   };
   
   // ******************************************** //
   // ***** CLASS STATICUSERDATA DEFINITIONS ***** //
   // ******************************************** //
   
   /** Default constructor. Calls UserDataStatic::finalize().*/
   template<class PARGRID> inline
   UserDataStatic<PARGRID>::UserDataStatic() {
      array=NULL; 
      basicDatatype=MPI_DATATYPE_NULL; 
      finalize();
   }

   /** Copy-constructor. Calls UserDataStatic::initialize and then copies the data.
    * @param udw UserDataStatic whose copy is to be made.*/
   template<class PARGRID> inline
   UserDataStatic<PARGRID>::UserDataStatic(const UserDataStatic& udw) {
      initialize(const_cast<PARGRID*>(udw.pargrid),udw.name,udw.N_cells,udw.N_elements,udw.byteSize,udw.datatype);
      for (size_t i=0; i<N_cells*N_elements*byteSize; ++i) array[i] = udw.array[i];
   }

   /** Destructor. Calls UserDataStatic::finalize().*/
   template<class PARGRID> inline
   UserDataStatic<PARGRID>::~UserDataStatic() {finalize();}
   
   /** Copy data from given UserDataStatic.
    * @param udw UserDataStatic whose data is to be copied.
    * @param newCell Target cell in which data is written.
    * @param oldCell Source cell in udw in which data is read.*/
   template<class PARGRID> inline
   void UserDataStatic<PARGRID>::copy(const UserDataStatic& udw,CellID newCell,CellID oldCell) {
      for (unsigned int i=0; i<N_elements*byteSize; ++i)
	array[newCell*N_elements*byteSize+i] = udw.array[oldCell*N_elements*byteSize+i];
   }
   
   /** Deallocates UserDataStatic::array and sets internal variables to dummy values.
    * @return If true, finalization completed successfully.*/
   template<class PARGRID> inline
   bool UserDataStatic<PARGRID>::finalize() {
      byteSize = 0;
      name = "";
      datatype = "";
      pargrid = NULL;
      N_elements = 0;
      N_cells = 0;
      delete [] array; array = NULL;
      return true;
   }
   
   /** Get a pointer to UserDataStatic::array.
    * @return Pointer to user data array. NULL value is returned for uninitialized array.*/
   template<class PARGRID> inline
   void* UserDataStatic<PARGRID>::getAddress() {return array;}
   
   /** Get a derived MPI datatype that will transfer all the given cells with a single send (or receive).
    * This function calls MPI_Type_commit for the newly created datatype, MPI_Type_free must be called elsewhere.
    * @param globalIDs Global IDs of the cells.
    * @param datatype Variable in which the derived datatype is to be written.*/
   template<class PARGRID> inline
   void UserDataStatic<PARGRID>::getDatatype(const std::set<CellID>& globalIDs,MPI_Datatype& datatype) {
      // For each global ID, get the corresponding local ID (=index into array) and 
      // store its offset relative to memory address MPI_BOTTOM:
      int* displacements = new int[globalIDs.size()];
      size_t counter = 0;
      for (std::set<CellID>::const_iterator cellGID=globalIDs.begin(); cellGID!=globalIDs.end(); ++cellGID) {
	 const CellID cellLID = pargrid->getLocalID(*cellGID);
	 displacements[counter] = cellLID;
	 ++counter;
      }

      // Create a derived MPI datatype and commit it:
      if (MPI_Type_create_indexed_block(globalIDs.size(),1,displacements,basicDatatype,&datatype) != MPI_SUCCESS) {
	 std::cerr << "(USERDATAWRAPPER) FATAL ERROR: Failed to create MPI datatype!" << std::endl;
	 exit(1);
      }
      MPI_Type_commit(&datatype);
      delete [] displacements; displacements = NULL;
   }

   template<class PARGRID> inline
   const std::string& UserDataStatic<PARGRID>::getName() const {return name;}
   
   /** Initialize UserDataStatic. Memory for UserDataStatic::array is allocated here.
    * @param pargrid Pointer to ParGrid class that called this function.
    * @param name Name of the user data array.
    * @param N_cells Number of cells in the array, should be equal to the total number of cells 
    * on this process (local + remote).
    * @param N_elements How many values of byte size byteSize are allocated per cell. This makes it 
    * possible to allocate, say, five doubles per cell.
    * @param byteSize Byte size of single value, i.e. sizeof(double) for doubles.
    * @param datatype Basic datatype stored in the array, either "int", "uint", or "float".*/
   template<class PARGRID> inline
   void UserDataStatic<PARGRID>::initialize(PARGRID* pargrid,const std::string& name,size_t N_cells,
					    unsigned int N_elements,unsigned int byteSize,const std::string& datatype) {
      this->pargrid = pargrid;
      this->name = name;
      this->byteSize = byteSize;
      this->N_elements = N_elements;
      this->N_cells = N_cells;
      this->datatype = datatype;
      array = new char[N_cells*N_elements*byteSize];
      
      // Create an MPI datatype that transfers a single array element:
      MPI_Type_contiguous(N_elements*byteSize,MPI_Type<char>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);
   }
   
   /** Swap the contents of two UserDataStatic containers.
    * @param uds UserDataStatic container to be swapped with.*/
   template<class PARGRID> inline
   void UserDataStatic<PARGRID>::swap(UserDataStatic& uds) {
      // Make dummy copies of member variables:
      char* dummyArray = array;
      const std::string dummyName = name;
      const unsigned int dummyByteSize = byteSize;
      const std::string dummyDatatype = datatype;
      const unsigned int dummy_N_cells = N_cells;
      const unsigned int dummy_N_elements = N_elements;

      // Set the contents of this UserDataStatic:
      array = uds.array;
      name = uds.name;
      byteSize = uds.byteSize;
      datatype = uds.datatype;
      N_cells = uds.N_cells;
      N_elements = uds.N_elements;
      MPI_Type_free(&basicDatatype);
      MPI_Type_contiguous(N_elements*byteSize,MPI_Type<char>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);
      
      // Set the contents ot the other UserDataStatic:
      uds.array = dummyArray;
      uds.name = dummyName;
      uds.byteSize = dummyByteSize;
      uds.datatype = dummyDatatype;
      uds.N_cells = dummy_N_cells;
      uds.N_elements = dummy_N_elements;
      MPI_Type_free(&uds.basicDatatype);
      MPI_Type_contiguous(uds.N_elements*uds.byteSize,MPI_Type<char>(),&uds.basicDatatype);
      MPI_Type_commit(&uds.basicDatatype);
   }

} // namespace pargrid

#endif
