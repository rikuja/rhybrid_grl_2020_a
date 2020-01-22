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

#ifndef PARGRID_COPY_PROTOCOL_MPI_H
#define PARGRID_COPY_PROTOCOL_MPI_H

#include <cmath>
#include "mpiconversion.h"

namespace pargrid {

   namespace mpiprotocol {
      const double BUFFER_INCREMENT_FACTOR = 1.2;     /**< If receiving buffer is too small, it is resized to 
						       * total number of copied elements times this factor.*/
      const size_t INITIAL_BUFFER_SIZE = 1;           /**< Initial size of receiving buffer.*/
      
      /** MPI tags used by CopyProtocolMPI class.*/
      enum tags {
	 BLOCK_SIZES,                                 /**< Message used to copy Buffer<T>::blockSizes vector.*/
	 BUFFER_1ST_COPY,                             /**< First message used to copy Buffer<T>::buffer data.*/
	 BUFFER_2ND_COPY,                             /**< Second message used to copy Buffer<T>::buffer data.*/
	 SIZE                                         /**< Total amount of MPI tags used by CopyProtocolMPI class.*/
      };
   }

   /** Class CopyProtocolMPI copies data between two pargrid::Buffer classes over MPI.
    * Class given as template parameter is the buffer class this CopyProtocolMPI should use.*/
   template<class BUFFER>
   class CopyProtocolMPI {
    public:
      CopyProtocolMPI();
      CopyProtocolMPI(const CopyProtocolMPI& cp);
      ~CopyProtocolMPI();

      bool clear();
      bool set(BUFFER* buffer,bool sender,MPI_Comm comm,pargrid::MPI_processID);
      bool start();
      bool wait(const std::string& name = "");
      
    protected:
      BUFFER* buffer;                          /**< BUFFER that is used in data copy.*/
      MPI_Comm comm;                           /**< MPI communicator that should be used.*/
      MPI_Datatype datatype;                   /**< Derived MPI datatype used to copy buffer contents.*/
      pargrid::MPI_processID myRank;           /**< MPI rank of this process in communicator comm.*/
      size_t partnerBufferSize;                /**< Size of partner's buffer, only significant at sender.*/
      pargrid::MPI_processID partnerRank;      /**< MPI rank of CopyProtocolMPI that this class is paired with.*/
      MPI_Request requests[mpiprotocol::SIZE]; /**< MPI requests used in buffer copy.*/
      bool sender;                             /**< If true, this class is sending buffer contents instead of receiving.*/
      bool setCalled;                          /**< If true, set member function has been successfully called.*/
      bool transferStarted;                    /**< If true, buffer contents are currently being copied.*/
      
      void invalidate();
      
      pargrid::MPI_processID test; // TEST
   };

   /** Constructor for class CopyProtocolMPI. Calls invalidate().
    * You need to call the following member functions before CopyProtocolMPI 
    * can be used: CopyProtocolMPI<BUFFER>::set.
    * @see CopyProtocolMPI<BUFFER>::set.*/
   template<class BUFFER> inline
   CopyProtocolMPI<BUFFER>::CopyProtocolMPI() {
      invalidate();
   }
   
   /** Copy-constructor for class CopyProtocolMPI. Calls invalidate().
    * This needs to be done in order to set MPI-related member variables to 
    * correct NULL states. Note that copy-constructor does not actually make a 
    * copy of the given class.
    * @param cp CopyProtocolMPI class to be copied.*/
   template<class BUFFER> inline
   CopyProtocolMPI<BUFFER>::CopyProtocolMPI(const CopyProtocolMPI<BUFFER>& cp) {
      invalidate();
   }

   /** Destructor for class CopyProtocolMPI(). Calls invalidate() to 
    * cleanly deallocate MPI-related variables.*/
   template<class BUFFER> inline
   CopyProtocolMPI<BUFFER>::~CopyProtocolMPI() {
      transferStarted = false;
      clear();
   }
   
   /** Reset class state to invalid. This function will deallocate all MPI-related variables.
    * This function will fail if CopyProtocolMPI is currently copying buffer contents.
    * @return If true, CopyProtocolMPI was successfully set to invalid state.*/
   template<class BUFFER> inline
   bool CopyProtocolMPI<BUFFER>::clear() {
      if (transferStarted == true) return false;
      
      buffer = NULL;
      //if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
      if (datatype != MPI_DATATYPE_NULL) MPI_Type_free(&datatype);
      partnerRank = -1;
      for (int i=0; i<mpiprotocol::SIZE; ++i) requests[i] = MPI_REQUEST_NULL;
      setCalled = false;
      transferStarted = false;

      partnerBufferSize = pargrid::mpiprotocol::INITIAL_BUFFER_SIZE;
      return true;
   }

   /** Invalidate CopyProtocolMPI class. Difference between this function and 
    * CopyProtocolMPI<BUFFER>::clear is that this function will not deallocate 
    * MPI-related variables.*/
   template<class BUFFER> inline
   void CopyProtocolMPI<BUFFER>::invalidate() {
      buffer = NULL;
      comm = MPI_COMM_NULL;
      datatype = MPI_DATATYPE_NULL;
      partnerBufferSize = pargrid::mpiprotocol::INITIAL_BUFFER_SIZE;
      partnerRank = -1;
      for (int i=0; i<mpiprotocol::SIZE; ++i) requests[i] = MPI_REQUEST_NULL;
      setCalled = false;
      transferStarted = false;
   }
   
   /** Pair CopyProtocolMPI with another instance of class and set the buffer to be copied.
    * You need to call this function on sending and receiving process to prevent
    * CopyProtocolMPI<BUFFER>::start and CopyProtocolMPI<BUFFER>::wait from deadlocking.
    * This function will fail if CopyProtocolMPI is currently copying buffer contents.
    * @param buffer Pointer to buffer used in data copy.
    * @param sender If true, then this CopyProtocolMPI is sending data instead of receiving.
    * @param comm MPI communicator to use.
    * @param partnerRank MPI rank of partner CopyProtocolMPI class within the given communicator.
    * @return If true, CopyProtocolMPI was successfully paired with another instance of class.*/
   template<class BUFFER> inline
   bool CopyProtocolMPI<BUFFER>::set(BUFFER* buffer,bool sender,MPI_Comm comm,pargrid::MPI_processID partnerRank) {
      // Cannot set new state if a buffer copy is in progress:
      if (transferStarted == true) return false;
      
      // If set has already been called, clear contents before continuing:
      if (setCalled == true) clear();

      bool success = true;
      if (buffer == NULL) success = false;
      this->buffer = buffer;
      this->sender = sender;

      // Make a copy of the given communicator:
      //MPI_Comm_dup(comm,&(this->comm));
      this->comm = comm;
      MPI_Comm_rank(comm,&myRank);

      // Check that partnerRank is valid:
      int commSize;
      MPI_Comm_size(comm,&commSize);
      this->partnerRank = partnerRank;
      if (partnerRank < 0 || partnerRank >= commSize) success = false;
      if (myRank == partnerRank) success = false;

      // If errors occurred reset communication protocol to invalid state and exit:
      if (success == false) {
	 clear();
	 return success;
      }
      
      // Create MPI datatype for sending a single element in buffer:
      const size_t elementByteSize = buffer->getElementByteSize();      
      MPI_Type_contiguous(elementByteSize,MPI_BYTE,&datatype);
      MPI_Type_commit(&datatype);
      
      // If this protocol is receiving data, resize buffer:
      buffer->setState(false,0);
      if (sender == false) {
	 if (buffer->setBufferSize(pargrid::mpiprotocol::INITIAL_BUFFER_SIZE) == false) success = false;
      }
      partnerBufferSize = pargrid::mpiprotocol::INITIAL_BUFFER_SIZE;
      
      setCalled = true;
      return success;
   }

   /** Start data copy between two buffers. This function will fail if CopyProtocolMPI<BUFFER>::set 
    * has not been called, or if CopyProtocolMPI is already copying data. Note that buffer 
    * contents are invalid until CopyProtocolMPI<BUFFER>::wait has been successfully called.
    * @return If true, data copy was successfully started.
    * @see CopyProtocolMPI<BUFFER>::wait.*/
   template<class BUFFER> inline
   bool CopyProtocolMPI<BUFFER>::start() {
      if (setCalled == false) return false;
      if (transferStarted == true) return false;
      
      const size_t N          = buffer->getNumberOfBlocks();
      uint32_t* blockSizes    = buffer->getBlockSizes();
      char* bufferPointer     = reinterpret_cast<char*>(buffer->getBufferPointer());
      const size_t bufferSize = buffer->getBufferSize();
      
      if (sender == true) {
	 // Receiver may not have sufficient capacity to receive all elements.
	 // Compute how many elements we can send with the first message:
	 const size_t N_total = bufferSize;
	 size_t N_firstMessage = 0;
	 size_t N_secondMessage = 0;
	 if (N_total > partnerBufferSize) {
	    N_firstMessage  = partnerBufferSize;
	    N_secondMessage = N_total - N_firstMessage;
	 } else {
	    N_firstMessage  = N_total;
	    N_secondMessage = 0;
	 }
	 
	 // Write total number of copied elements, and number of
	 // elements copied with first message, to end of array blockSizes:
	 blockSizes[N+pargrid::buffermetadata::N_ELEMENTS_TOTAL] = N_total;
	 blockSizes[N+pargrid::buffermetadata::N_ELEMENTS_FIRST] = N_firstMessage;

	 // Send blockSizes array and first batch of buffer elements:
	 MPI_Isend(blockSizes,N+pargrid::buffermetadata::SIZE,MPI_Type<uint32_t>(),partnerRank,pargrid::mpiprotocol::BLOCK_SIZES,comm,&(requests[0]));
	 MPI_Isend(bufferPointer,N_firstMessage,datatype,partnerRank,pargrid::mpiprotocol::BUFFER_1ST_COPY,comm,&(requests[1]));

	 // If receiving buffer does not have capacity to receive all elements with 
	 // a single message, send the rest of elements with a second message and 
	 // increase perceived partner buffer capacity:
	 if (N_secondMessage > 0) {
	    const size_t elementByteSize = buffer->getElementByteSize();
	    MPI_Isend(bufferPointer+elementByteSize*N_firstMessage,N_secondMessage,datatype,partnerRank,pargrid::mpiprotocol::BUFFER_2ND_COPY,comm,&(requests[2]));
	    partnerBufferSize = static_cast<size_t>(floor(N_total*pargrid::mpiprotocol::BUFFER_INCREMENT_FACTOR));
	 }	 
      } else {
	 // Post receives for blockSizes array, and at most bufferSize elements:
	 MPI_Irecv(blockSizes,N+pargrid::buffermetadata::SIZE,MPI_Type<uint32_t>(),partnerRank,pargrid::mpiprotocol::BLOCK_SIZES,comm,&(requests[0]));
	 MPI_Irecv(bufferPointer,bufferSize,datatype,partnerRank,pargrid::mpiprotocol::BUFFER_1ST_COPY,comm,&(requests[1]));
      }
      
      // Lock buffer:
      buffer->setState(true,0);
      transferStarted = true;
      return true;
   }

   /** Wait until buffer contents have been copied. This function will fail 
    * if CopyProtocolMPI<BUFFER>::start has not been called prior to calling
    * this function.
    * @return If true, buffer contents were successfully copied.*/
   template<class BUFFER> inline
   bool CopyProtocolMPI<BUFFER>::wait(const std::string& name) {
      if (transferStarted == false) return false;

      // Wait for MPI copies to finish:
      #ifndef NDEBUG
         int flag = false;
         uint64_t waitTime      = 10000;
         uint64_t waitedTime    = 0;
         uint64_t maxWaitedTime = 5000000000;
         timespec timeSpec;
         timeSpec.tv_sec  = 0;
         timeSpec.tv_nsec = waitTime;
         do {
	    MPI_Testall(3,requests,&flag,MPI_STATUSES_IGNORE);
	    nanosleep(&timeSpec,NULL);
	    waitedTime += waitTime;
	    if (flag == true) break;
	    if (waitedTime >= maxWaitedTime) {
	       std::cerr << "Killing execution in ParGrid::CopyProtocolMPI with name '" << name << "'" << std::endl;
	       exit(1);
	    }
	 } while(true);
      #else
         MPI_Waitall(3,requests,MPI_STATUSES_IGNORE);
      #endif
      
      // Unlock buffer:
      buffer->setState(false,0); // FIXME

      if (sender == false) {
	 // Check if we were able to receive all buffer elements with the first message.
	 // If not, increase buffer capacity and post receive for rest of data:
	 const size_t N          = buffer->getNumberOfBlocks();
	 uint32_t* blockSizes    = buffer->getBlockSizes();
	 const size_t bufferSize = buffer->getBufferSize();

	 // Total number of copied elements:
	 const size_t N_total        = blockSizes[N+pargrid::buffermetadata::N_ELEMENTS_TOTAL];
	 // Number of elements copied with first message:
	 const size_t N_firstMessage = blockSizes[N+pargrid::buffermetadata::N_ELEMENTS_FIRST];

	 if (N_total > bufferSize) {
	    // Increase buffer size:
	    const size_t newSize = static_cast<size_t>(floor(N_total*pargrid::mpiprotocol::BUFFER_INCREMENT_FACTOR));
	    buffer->setBufferSize(newSize);
	    
	    // Receive elements that did not fit into first message:
	    const size_t N_secondMessage = N_total - N_firstMessage;
	    char* bufferPointer = reinterpret_cast<char*>(buffer->getBufferPointer());
	    const size_t elementByteSize = buffer->getElementByteSize();
	    MPI_Irecv(bufferPointer+elementByteSize*N_firstMessage,N_secondMessage,datatype,partnerRank,pargrid::mpiprotocol::BUFFER_2ND_COPY,comm,&(requests[2]));
	    MPI_Waitall(3,requests,MPI_STATUSES_IGNORE);
	 }

	 buffer->setState(false,N_total);
      }
      transferStarted = false;
      return true;
   }
      
} // namespace pargrid

#endif
