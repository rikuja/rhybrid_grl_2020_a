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

#ifndef PARGRID_DATAWRAPPER_H
#define PARGRID_DATAWRAPPER_H

#include <cstdlib>
#include <iostream>
#include "pargrid_userdata_dynamic.h"

namespace pargrid {
   /** Wrapper for class pargrid::UserDataDynamic. This class defines 
    * the interface through which end user can read and modify the 
    * contents of dynamic data arrays in ParGrid.*/
   template<typename T>
   class DataWrapper {
    public:
      DataWrapper(char** arrays,ArraySizetype* capacities,CellID N_cells,ArraySizetype* sizes,uint64_t elementByteSize);
      const ArraySizetype* capacity() const;
      ArraySizetype capacity(CellID cell) const;
      T** data() const;
      void push_back(CellID cell,const T& element);
      void recapacitate(CellID cell,ArraySizetype newCapacity);
      void reserve(CellID cell,ArraySizetype newCapacity);
      void resize(CellID cell,ArraySizetype newSize);
      const ArraySizetype* size() const;
      ArraySizetype size(CellID cell) const;
      bool valid() const;
      
    private:
      T** arrays;                   /**< Pointers to arrays containing user data for each cell.*/
      ArraySizetype* capacities;    /**< Current capacity of each cell, counted in array elements having byte size elementByteSize.*/
      uint64_t elementByteSize;     /**< Byte size of an array element.*/
      CellID N_cells;               /**< Number of cells in arrays, capacities, and sizes.*/
      ArraySizetype* sizes;         /**< Current size of each cell, counted in array elements having byte size elementByteSiz.*/
      
      /** Private default constructor to prevent the creation 
       * of DataWrappers that point to an invalid ParGrid dynamic data array.*/
      DataWrapper();
   };

   /** Create a new DataWrapper for given ParGrid dynamic data array.
    * @param arrays Arrays containing the data for each cell.
    * @param capacities Array containing current capacity for each cell, counted in array elements having byte size elementByteSize.
    * @param N_cells Number of cells in the mesh local to this process.
    * @param sizes Array containing current size for each cell, counted in array elements having byte size elementByteSize.
    * @param elementByteSize Byte size of an array element.*/
   template<typename T> inline
   DataWrapper<T>::DataWrapper(char** arrays,ArraySizetype* capacities,CellID N_cells,ArraySizetype* sizes,uint64_t elementByteSize) {
      this->arrays          = reinterpret_cast<T**>(arrays);
      this->capacities      = capacities;
      this->N_cells         = N_cells;
      this->sizes           = sizes;
      this->elementByteSize = elementByteSize;
   }

   /** Get the array containing current capacities of cells.
    * @return Pointer to capacity array.*/
   template<typename T> inline
   const ArraySizetype* DataWrapper<T>::capacity() const {return capacities;}
   
   /** Get the current capacity of the given cell.
    * @param cell Local ID of the cell.
    * @return Current capacity of the cell.*/
   template<typename T> inline
   ArraySizetype DataWrapper<T>::capacity(CellID cell) const {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: cell #" << cell << " out of bounds in capacity!" << std::endl;
	    exit(1);
	 }
      #endif
      return capacities[cell];
   }
   
   /** Get arrays containg user data. Number of elements in each cell 
    * can be obtained from the array returned by size() member function.
    * @return Pointers to arrays containing dynamic user data.*/
   template<typename T> inline
   //T** DataWrapper<T>::data() const {return reinterpret_cast<T**>(arrays);}
   T** DataWrapper<T>::data() const {return arrays;}

   /** Insert a new element to the given cell. This function increases the 
    * capacity of the cell by a factor of two if the current capacity is too 
    * small to hold the new element.
    * @param cell Local ID of the cell.
    * @param element Element to be inserted.*/
   template<typename T> inline
   void DataWrapper<T>::push_back(CellID cell,const T& element) {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: Given cell ID #" << cell << " out of bounds in push_back!" << std::endl;
	    exit(1);
	 }         
      #endif

      // If array does not have sufficient capacity for a new element,
      // increase the capacity, create a new array, and copy contents 
      // from the old array before inserting the new element:
      if (sizes[cell]+1 > capacities[cell]) {
	 capacities[cell] = 2*(sizes[cell]+1);
	 T* tmp = new T[capacities[cell]];
	 for (ArraySizetype i=0; i<sizes[cell]; ++i) tmp[i] = arrays[cell][i];
	 delete [] arrays[cell];
	 arrays[cell] = tmp;
      }
      
      #ifndef NDEBUG
         if (arrays[cell] == NULL) {
	    std::cerr << "(PARGRID) ERROR: DataWrapper array[" << cell << "] is NULL!" << std::endl;
	    exit(1);
	 }
      #endif
      
      // Insert given element and increase array size:
      arrays[cell][sizes[cell]] = element;
      ++sizes[cell];
   }

   /** Force the capacity of the given cell to be newCapacity. If newCapacity is 
    * less than the current capacity, a new array with smaller capacity is allocated 
    * and elements [0,newCapacity-1] are copied from the old array into the new one. 
    * The difference between this function and reserve is that this function will 
    * also reduce array capacity. Note that this function will not call constructors 
    * or destructors for inserted or removed data elements.
    * @param cell Local ID of the cell.
    * @param newCapacity New capacity for the cell.*/
   template<typename T> inline
   void DataWrapper<T>::recapacitate(CellID cell,ArraySizetype newCapacity) {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: cell #" << cell << " out of bounds in reserve!" << std::endl;
	    exit(1);
	 }
      #endif
      
      if (newCapacity < capacities[cell]) {
	 // If requested capacity is smaller than the current one, 
	 // create a new array with reduces capacity and copy 
	 // some of the contents from the old array:
	 capacities[cell] = newCapacity;
	 T* tmp = new T[newCapacity];
	 for (ArraySizetype i=0; i<newCapacity; ++i) tmp[i] = arrays[cell][i];
	 delete [] arrays[cell];
	 arrays[cell] = tmp;
	 sizes[cell] = newCapacity;
      } else {
	 // If requested capacity is larger than the current one,
	 // just call reserve:
	 reserve(cell,newCapacity);
      }
   }
   
   /** Make a request that given cell is able to hold at least 
    * newCapacity number of elements. This function will increase 
    * the current capacity of the given cell if the requested capacity 
    * is larger than the current one. Note that this function will not 
    * reduce current capacity.
    * @param cell Local ID of the cell.
    * @param newCapacity Requested new capacity.
    * @see recapacitate.*/
   template<typename T> inline
   void DataWrapper<T>::reserve(CellID cell,ArraySizetype newCapacity) {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: cell #" << cell << " out of bounds in reserve!" << std::endl;
	    exit(1);
	 }
      #endif
      
      // If requested capacity is larger than the current capacity,
      // create a new array with requested capacity and copy
      // contents from the old one:
      if (newCapacity > capacities[cell]) {
	 capacities[cell] = newCapacity;
	 T* tmp = new T[newCapacity];
	 for (ArraySizetype i=0; i<sizes[cell]; ++i) tmp[i] = arrays[cell][i];
	 delete [] arrays[cell];
	 arrays[cell] = tmp;
      }
   }
   
   /** Resize, i.e. insert or remove elements to/from, the given cell. If 
    * newSize is larger than the current capacity, a new array with capacity equal 
    * to newSize is allocated. Note that this function will not call 
    * constructors or destructors for the inserted or removed elements.
    * @param cell Local ID of the cell.
    * @param newSize New size of the cell.*/
   template<typename T> inline
   void DataWrapper<T>::resize(CellID cell,ArraySizetype newSize) {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: Given cell ID #" << cell << " out of bounds in resize!" << std::endl;
	    exit(1);
	 }
      #endif
      
      if (newSize <= capacities[cell]) {
	 // Just reduce array size:
	 sizes[cell] = newSize;
      } else {
	 // Old capacity is too small, increase it and 
	 // copy contents from old array into a new one:
	 capacities[cell] = newSize;
	 char* ptr = reinterpret_cast<char*>(arrays[cell]);
	 char* tmp = new char[newSize*elementByteSize];
	 for (uint64_t i=0; i<sizes[cell]*elementByteSize; ++i) tmp[i] = ptr[i];
	 delete [] arrays[cell];
	 arrays[cell] = reinterpret_cast<T*>(tmp);
	 
	 // Resize array:
	 sizes[cell] = newSize;
      }
   }
   
   /** Get array containing the current size of each cell.
    * @return Pointer to array containing the size of each cell.*/
   template<typename T> inline
   const ArraySizetype* DataWrapper<T>::size() const {return sizes;}

   /** Get current size of given cell.
    * @param cell Local ID of the cell.
    * @return Current number of elements in the cell.*/
   template<typename T> inline
   ArraySizetype DataWrapper<T>::size(CellID cell) const {
      #ifndef NDEBUG
         if (cell >= N_cells) {
	    std::cerr << "(PARGRID DATAWRAPPER) ERROR: Cell #" << cell << " out of bounds in size!" << std::endl;
	    exit(1);
	 }
      #endif
      return sizes[cell];
   }
   
   /** If true, this DataWrapper points to a valid dynamic ParGrid data array.
    * DataWrappers become invalid, e.g. after repartitioning.
    * @return If true, this DataWrapper points to a valid dynamic
    * data array and can be used safely.*/
   template<typename T> inline
   bool DataWrapper<T>::valid() const {
      if (arrays == NULL) return false;
      return true;
   }
   
} // namespace pargrid

#endif
