/** This file is part of ParGrid parallel grid.
 * 
 *  Copyright 2011-2013 Finnish Meteorological Institute
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

#ifndef PARGRID_BUFFERS_H
#define PARGRID_BUFFERS_H

#include <vector>
#include "pargrid_definitions.h"

namespace pargrid {

   namespace buffermetadata {
      /** Buffer class allocates extra elements in std::vector<uint32_t> Buffer<T>::blockSizes 
       * that are used to exchange metadata between sender and receiver. This enumeration
       * defines what those extra elements are.*/
      enum ExtraBlocks {
         N_ELEMENTS_TOTAL,                   /**< Total number of elements in buffer.*/
	 N_ELEMENTS_FIRST,                   /**< Number of elements copied with first message.*/
	 SIZE                                /**< Total number of extra metadata blocks in blockSizes array.*/
      };
   }

   /** Class for copying buffered data between source and destination.
    * Basic element stored in Buffer is given by template parameter T. 
    * It must be byte-copyable, i.e. cannot contain dynamically allocated data, 
    * and call to sizeof(T) must return its byte size.
    * Elements are grouped into one or more blocks, number of elements in 
    * each block is stored in array returned by method getBlockSizeArray.
    * For example, typically number of blocks is the number of cells 
    * (in the simulation mesh) whose data is copied, and block size could 
    * be the number of particles in each cell.
    * 
    * Note that the actual data copy is performed by one of the CopyProtocol 
    * classes defined in namespace pargrid, e.g. pargrid::CopyProtocolMPI.*/
   template<typename T>
   class Buffer {
    public:
      Buffer();

      uint32_t* getBlockSizes();
      std::vector<T>& getBuffer();
      void* getBufferPointer();
      size_t getBufferSize() const;
      size_t getElementByteSize() const;
      size_t getNumberOfBlocks() const;
      size_t getNumberOfCopiedElements() const;
      bool resize(uint32_t N_blocks);
      bool setBufferSize(size_t newSize);
      void setState(bool copyOngoing,size_t N_copiedElements);

    private:
      std::vector<uint32_t> blockSizes;               /**< Number of elements in each block, e.g. 
						       * number of particles in block.
						       * numberArray[N-2] = number of elements in first send,
						       * numberArray[N-1] = total number of incoming elements,
						       * where N is the size of numberArray.*/
      std::vector<T> buffer;                          /**< Receive buffer.*/
      size_t bufferCapacity;                          /**< Current capacity of buffer, in number of elements.*/
      bool copyOngoing;                               /**< If true, buffer contents are currently being copied.*/
      std::vector<T> dummyBuffer;                     /**< Empty buffer, it is returned by member functions if
						       * buffer contents cannot be accessed at the time.*/
      size_t N_copiedElements;                        /**< Number of elements that were actually copied. This 
						       * is equal to or less than the current buffer size.*/
   };
   
   /** Constructor for Buffer.*/
   template<typename T> inline
   Buffer<T>::Buffer() {
      copyOngoing = false;
      N_copiedElements = 0;
   }

   /** Get pointer to array containing number of elements in each message block.
    * If buffer contents are currently being copied then a NULL pointer is returned.
    * @return Array containing current number of elements in each message block.*/
   template<typename T> inline
   uint32_t* Buffer<T>::getBlockSizes() {
      // If buffer is being copied its contents cannot be read or written:
      if (copyOngoing == true) return NULL;
      return &(this->blockSizes[0]);
   }
   
   /** Get buffer where copied elements are stored. If buffer contents 
    * are currently being copied then a vector of size zero is returned.
    * @return Vector containing elements currently in buffer.*/
   template<typename T> inline
   std::vector<T>& Buffer<T>::getBuffer() {
      if (copyOngoing == true) return dummyBuffer;
      return buffer;
   }

   /** Get a pointer to array containing elements currently in buffer.
    * This function is used by CopyProtocol classes.
    * @return Pointer to array used by std::vector<T> buffer.*/
   template<typename T> inline
   void* Buffer<T>::getBufferPointer() {return &(buffer[0]);}

   /** Get number of elements in buffer. This function is used by CopyProtocol classes.
    * @return Size of std::vector<T> buffer.*/
   template<typename T> inline
   size_t Buffer<T>::getBufferSize() const {return buffer.size();}
   
   /** Get byte size of element stored in buffer. This function is used by CopyProtocol classes.
    * @return Byte size of element, as given by sizeof function.*/
   template<typename T> inline
   size_t Buffer<T>::getElementByteSize() const {return sizeof(T);}
   
   /** Get number of message blocks. This is the size of array 
    * returned by Buffer::getBlockSizes().
    * @return Number of message blocks.*/
   template<typename T> inline
   size_t Buffer<T>::getNumberOfBlocks() const {return blockSizes.size() - buffermetadata::SIZE;}

   /** Get number of elements that were actually copied by CopyProtocol class.
    * This is equal to or less than the size of vector returned by Buffer::getBuffer().
    * @return Number of elements copied.*/
   template<typename T> inline
   size_t Buffer<T>::getNumberOfCopiedElements() const {return N_copiedElements;}
   
   /** Change the amount of message blocks used by Buffer. This function
    * fails if buffer contents are currently being copied.
    * @param N_blocks New amount of message blocks.
    * @return If true, amount of message blocks was successfully changed.*/
   template<typename T> inline
   bool Buffer<T>::resize(uint32_t N_blocks) {
      // Buffer cannot be changed if its contents are being copied:
      if (copyOngoing == true) return false;
      
      // Reallocate blockSizes array:
	{
	   std::vector<uint32_t> dummy;
	   blockSizes.swap(dummy);
	}
      blockSizes.resize(N_blocks + buffermetadata::SIZE);
      N_copiedElements = 0;
      return true;
   }

   /** Change the size of buffer, calls std::vector<T> Buffer<T>::buffer.resize function.
    * This function is used by CopyProtocol classes. This function will fail if 
    * buffer contents are currently being copied.
    * @param newSize New buffer size.
    * @return If true, buffer size was successfully changed.*/
   template<typename T> inline
   bool Buffer<T>::setBufferSize(size_t newSize) {
      if (copyOngoing == true) return false;
      buffer.resize(newSize);
      return true;
   }

   /** Set buffer state. This function is used by CopyProtocol classes.
    * @param copyOngoing If true, buffer contents are being copied.
    * @param N_copiedElements Number of elements that were copied.*/
   template<typename T> inline
   void Buffer<T>::setState(bool copyOngoing,size_t N_copiedElements) {
      this->copyOngoing = copyOngoing;
      this->N_copiedElements = N_copiedElements;
   }
   
} // namespace pargrid

#endif
