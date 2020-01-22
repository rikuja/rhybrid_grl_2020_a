/** This file is part of ParGrid parallel grid.
 * 
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

#ifndef PARGRID_STENCIL_H
#define PARGRID_STENCIL_H

#include <time.h>
#include <map>
#include <set>
#include <vector>
#include <mpi.h>
#include "pargrid_definitions.h"

namespace pargrid {
   
   // ************************************** //
   // ***** CLASS STENCIL DECLARATIONS ***** //
   // ************************************** //
   
   /** Class used to calculate send and receive lists for user-defined stencils in ParGrid. 
    * This class does most of the dirty work related to transferring data with MPI.*/
   template<class PARGRID,class C>
   struct Stencil {
    public:
      Stencil();
      ~Stencil();
      
      bool addUserDataTransfer(DataID userDataID,bool isDynamic);
      void clear();
      const std::vector<CellID>& getBoundaryCells() const;
      const std::vector<CellID>& getInnerCells() const;
      const std::map<MPI_processID,std::set<CellID> >& getRecvs() const;
      const std::map<MPI_processID,std::set<CellID> >& getSends() const;
      bool getRemoteUpdates(DataID userDataID,unsigned int*& offsets,char*& buffer) const;
      bool initialize(PARGRID& pargrid,StencilType stencilType,const std::vector<NeighbourID>& receives);
      bool removeTransfer(DataID userDataID);
      bool startTransfer(DataID userDataID);
      bool update();
      bool wait(DataID userDataID);
      bool wait(DataID userDataID,const std::string& name);
      
    private:
      /** Wrapper for MPI_Datatype. This was created to make copying and deletion of
       * MPI_Datatype work correctly. Calls MPI_Type_dup in appropriate places.*/
      struct TypeWrapper {
	 MPI_Datatype type;                               /**< MPI datatype.*/
	 
	 TypeWrapper();
	 TypeWrapper(const TypeWrapper& tw);
	 ~TypeWrapper();
	 TypeWrapper& operator=(const TypeWrapper& tw);
      };
      
      /** MPI Datatype cache for data transfer.*/
      struct TypeCache {
	 std::vector<TypeWrapper> recvs; /**< MPI datatypes for receiving data.*/
	 std::vector<TypeWrapper> sends; /**< MPI datatypes for sending data.*/
      };
    
      /** Information on data transfer.*/
      struct TypeInfo {
	 int N_receives;          /**< Total number of messages received during this transfer.*/
	 int N_sends;             /**< Total number of messages sent during this transfer.*/
	 bool typeVolatile;       /**< If true, MPI Datatypes need to be recalculated every time 
				   * this transfer is started.*/
	 DataID userDataID;       /**< ID of user-defined array. Valid if this transfer is associated with a user-defined array.*/
	 bool started;            /**< If true, this transfer has started and MPI requests are valid.*/
	 MPI_Request* requests;   /**< MPI requests associated with this transfer.*/
      };
      
      /** Receive buffer for StencilType::remoteToLocalUpdates. These arrays 
       * are only allocated if needed. RecvBuffer is only used for 
       * user-defined ParGrid data arrays.*/
      struct RecvBuffer {
	 RecvBuffer();
	 RecvBuffer(const RecvBuffer& rbuffer);
	 ~RecvBuffer();

	 int offsetsSize;
	 int bufferSize;
	 int elementSize;
	 unsigned int* offsets;
	 char* buffer;
      };
      
      std::vector<CellID> boundaryCells;                           /**< List of boundary cells of this stencil.*/
      std::vector<CellID> innerCells;                              /**< List of inner cells of this stencil.*/
      bool initialized;                                            /**< If true, Stencil has initialized successfully and is ready for use.*/
      std::vector<NeighbourID> receivedNbrTypeIDs;                 /**< Neighbour type IDs indicating which cells to receive data from.*/
      std::vector<NeighbourID> sentNbrTypeIDs;                     /**< Neighbour type IDs indicating which cells to send data.*/
      
      std::map<DataID,std::map<MPI_processID,TypeCache> > typeCachesUser;
      std::map<DataID,TypeInfo> typeInfoUser;
      
      PARGRID* parGrid;                                            /**< Pointer to parallel grid.*/
      StencilType stencilType;
      std::map<MPI_processID,std::set<CellID> > recvs;             /**< List of received cells (ordered by global ID) from each neighbour process.*/
      std::map<MPI_processID,std::set<CellID> > sends;             /**< List of cells sent (ordered by global ID) to each neighbour process.*/
      std::map<CellID,std::set<MPI_processID> > recvCounts;        /**< For each local cell, identified by global ID, ranks of remote 
								    * processes whom to receive an update from. This map is only 
								    * used for stencilType remoteToLocalUpdates.*/
      std::map<DataID,RecvBuffer*> recvBuffers;                    /**< Receive buffers for StencilType::remoteToLocalUpdates. These 
								    * are allocated only when needed.*/
      
      bool calcLocalUpdateSendsAndReceives();
      bool calcTypeCacheUser(DataID userDataID);
      bool calcTypeCacheDynamic(DataID userDataID,TypeInfo& info);
      bool calcTypeCacheStatic(DataID userDataID,TypeInfo& info);
      bool startStaticTransfer(DataID userDataID);
      bool waitStatic(DataID userDataID);
      bool waitStatic(DataID userDataID,const std::string& name);
   };
         
   // ************************************************** //
   // ***** CLASS STENCIL::TYPEWRAPPER DEFINITIONS ***** //
   // ************************************************** //
   
   /** Default constructor for TypeWrapper. Sets type to MPI_DATATYPE_NULL.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::TypeWrapper::TypeWrapper(): type(MPI_DATATYPE_NULL) { }
   
   /** Copy-constructor for TypeWrapper. Makes a copy of MPI datatype with 
    * MPI_Type_dup unless the given datatype is MPI_DATATYPE_NULL.
    * @param tw TypeWrapper to be copied.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::TypeWrapper::TypeWrapper(const Stencil<PARGRID,C>::TypeWrapper& tw) {
      if (tw.type != MPI_DATATYPE_NULL) MPI_Type_dup(tw.type,&type);
      else type = MPI_DATATYPE_NULL;
   }
   
   /** Destructor for TypeWrapper. Frees the MPI datatype with MPI_Type_free 
    * unless the datatype is MPI_DATATYPE_NULL.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::TypeWrapper::~TypeWrapper() {
      if (type != MPI_DATATYPE_NULL) MPI_Type_free(&type);
   }
   
   /** Assignment operator for TypeWrapper. Makes a copy of the given MPI datatype 
    * with MPI_Type_dup, unless the given datatype is MPI_DATATYPE_NULL. Also frees 
    * the current datatype with MPI_DATATYPE_NULL if necessary.
    * @param tw TypeWrapper to make copy of.
    * @return Reference to this TypeWrapper.*/
   template<class PARGRID,class C> inline
   typename Stencil<PARGRID,C>::TypeWrapper& Stencil<PARGRID,C>::TypeWrapper::operator=(const Stencil<PARGRID,C>::TypeWrapper& tw) {
      if (type != MPI_DATATYPE_NULL) MPI_Type_free(&type);
      if (tw.type == MPI_DATATYPE_NULL) type = MPI_DATATYPE_NULL;
      else MPI_Type_dup(tw.type,&type);
      return *this;
   }

   // ************************************************* //
   // ***** CLASS STENCIL::RECVBUFFER DEFINITIONS ***** //
   // ************************************************* //
   
   /** Default constructor for RecvBuffer. Sets array pointers to NULL.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::RecvBuffer::RecvBuffer(): offsetsSize(0),bufferSize(0),elementSize(0),offsets(NULL),buffer(NULL) { }
   
   /** RecvBuffer copy-constructor. Allocates new offsets and buffer arrays.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::RecvBuffer::RecvBuffer(const RecvBuffer& rbuffer) {
      offsetsSize = rbuffer.offsetsSize;
      bufferSize = rbuffer.offsetsSize;
      elementSize = rbuffer.elementSize;
      offsets = new unsigned int[offsetsSize];
      buffer  = new char[bufferSize];
      for (size_t i=0; i<offsetsSize; ++i) offsets[i] = rbuffer.offsets[i];
      for (size_t i=0; i<bufferSize; ++i) buffer[i] = rbuffer.buffer[i];
   }
   
   /** Destructor for RecvBuffer. Calls delete for arrays.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::RecvBuffer::~RecvBuffer() {
      delete [] offsets; offsets = NULL;
      delete [] buffer; buffer = NULL;
   }

   // ************************************* //
   // ***** CLASS STENCIL DEFINITIONS ***** //
   // ************************************* //
   
   /** Default constructor.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::Stencil() { }
   
   /** Destructor. Frees all MPI datatypes and deallocates arrays.*/
   template<class PARGRID,class C> inline
   Stencil<PARGRID,C>::~Stencil() {
      // Delete MPI Requests (static data):
      for (typename std::map<DataID,TypeInfo>::iterator i=typeInfoUser.begin(); i!=typeInfoUser.end(); ++i) {
	 delete [] i->second.requests; i->second.requests = NULL;
      }
      // Delete receive buffers:
      for (typename std::map<DataID,RecvBuffer*>::iterator it=recvBuffers.begin(); it!=recvBuffers.end(); ++it) {
	 delete it->second;
	 it->second = NULL;
      }
   }

   /** Add a transfer for user-defined ParGrid data array.
    * @param userDataID ID number of the user data array.
    * @param transfer ID number of the transfer.
    * @param isDynamic If true, the associated MPI datatypes should be recalculated 
    * every time transfers are started.
    * @return If true, transfer was added to Stencil successfully.*/
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::addUserDataTransfer(DataID userDataID,bool isDynamic) {
      if (initialized == false) return false;
      if (typeCachesUser.find(userDataID) != typeCachesUser.end()) {
	 return false;
      }
      typeCachesUser[userDataID];
      typeInfoUser[userDataID].typeVolatile = isDynamic;
      typeInfoUser[userDataID].requests     = NULL;
      typeInfoUser[userDataID].started      = false;
      typeInfoUser[userDataID].userDataID   = userDataID;
      calcTypeCacheUser(userDataID);
      return true;
   }
   
   /** Calculate send and receive lists for this Stencil.
    * @return If true, send and receive lists were calculated successfully.*/
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::calcLocalUpdateSendsAndReceives() {
      if (initialized == false) return false;
      bool success = true;
      clear();
	
      // Iterate over all local cells' neighbours. If the neighbour is not in 
      // sentNbrTypeIDs or receivedNbrTypeIDs, skip it.
      // If the neighbour is remote, add the local cell into sends and the neighbour into recvs. 
      // All local cells with one or more remote neighbours are inserted into boundaryCells.
      // Cells with zero remote neighbours are inserted into innerCells instead.
      const std::vector<CellID>& globalIDs = parGrid->getGlobalIDs();
      for (CellID i=0; i<parGrid->getNumberOfLocalCells(); ++i) {
	 unsigned int N_remoteNeighbours = 0;
	 CellID* nbrIDs = parGrid->getCellNeighbourIDs(i);
	 for (size_t nbr=0; nbr<N_neighbours; ++nbr) {
	    // Check that neighbour exists and that it is not local:
	    const CellID nbrLocalID = nbrIDs[nbr];
	    if (nbrLocalID == parGrid->invalid()) continue;
	    const MPI_processID nbrHost = parGrid->getHosts()[nbrLocalID];
	    if (nbrHost == parGrid->getRank()) continue;
	    const CellID nbrGlobalID = globalIDs[nbrLocalID];
	    
	    // If neighbour type ID is in sentNbrTypeIDs, add a send:
	    if (std::find(sentNbrTypeIDs.begin(),sentNbrTypeIDs.end(),nbr) != sentNbrTypeIDs.end()) {
	       if (stencilType == localToRemoteUpdates) sends[nbrHost].insert(globalIDs[i]);
	       else sends[nbrHost].insert(nbrGlobalID);
	    }
	    
	    // If neighbour type ID is in receivedNbrTypeIDs, add a receive:
	    if (std::find(receivedNbrTypeIDs.begin(),receivedNbrTypeIDs.end(),nbr) != receivedNbrTypeIDs.end()) {
	       if (stencilType == localToRemoteUpdates) {
		  recvs[nbrHost].insert(nbrGlobalID);
	       } else {
		  recvs[nbrHost].insert(globalIDs[i]);
		  recvCounts[globalIDs[i]].insert(nbrHost);
	       }
	       ++N_remoteNeighbours;
	    }
	 }
	 // Add local cell either into innerCells or boundaryCells:
	 if (N_remoteNeighbours == 0) innerCells.push_back(i);
	 else boundaryCells.push_back(i);
      }
      return success;
   }

   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::calcTypeCacheUser(DataID userDataID) {
      typename std::map<DataID,TypeInfo>::iterator info = typeInfoUser.find(userDataID);
      if (info == typeInfoUser.end()) return false;
      return calcTypeCacheStatic(userDataID,info->second);
   }

   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::calcTypeCacheDynamic(DataID userDataID,TypeInfo& info) {
      bool success = true;
      return success;
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::calcTypeCacheStatic(DataID userDataID,TypeInfo& info) {
      // Free old datatypes:
      typename std::map<DataID,std::map<MPI_processID,TypeCache> >::iterator it=typeCachesUser.find(userDataID);
      for (typename std::map<MPI_processID,TypeCache>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
	 jt->second.sends.clear();
	 jt->second.recvs.clear();
      }
      
      // Delele old MPI_Requests:
      delete [] info.requests; info.requests = NULL;
      
      // Remove old process entries from typeCache:
      it->second.clear();
      
      // Insert entry for each neighbouring process:
      for (std::map<MPI_processID,std::set<CellID> >::const_iterator i=sends.begin(); i!=sends.end(); ++i) {
	 (it->second)[i->first];
      }
      for (std::map<MPI_processID,std::set<CellID> >::const_iterator i=recvs.begin(); i!=recvs.end(); ++i) {
	 (it->second)[i->first];
      }
      
      // Insert entry to recvBuffers if one does not already exist:
      if (stencilType == remoteToLocalUpdates) {
	 // Allocate offset array for each local boundary cell and 
	 // calculate offsets into receive buffer. New RecvBuffer is not 
	 // allocated if one already exists:
	 std::pair<typename std::map<DataID,RecvBuffer*>::iterator,bool> result
	   = recvBuffers.insert(std::make_pair(userDataID,(RecvBuffer*)NULL));
	 
	 // Allocate arrays if a new RecvBuffer was created:
	 if (result.second == true) {
	    result.first->second = new RecvBuffer();
	    result.first->second->offsetsSize = recvCounts.size()+1;
	    result.first->second->offsets = new unsigned int[recvCounts.size()+1];
	    unsigned int* offsets = result.first->second->offsets;
	    offsets[0] = 0;
	    for (size_t cell=0; cell<boundaryCells.size(); ++cell) {
	       offsets[cell+1] = offsets[cell] + recvCounts[parGrid->getGlobalIDs()[boundaryCells[cell]]].size();
	    }

	    // Allocate memory for receive buffer:
	    const size_t N_receivedCells = offsets[recvCounts.size()];
	    const size_t elementByteSize = parGrid->getUserDataElementSize(userDataID);
	    if (elementByteSize == 0) {
	       std::cerr << "(PARGRID) Stencil::calcTypeCache ERROR: User data array element size is zero!" << std::endl;
	       exit(1);
	    }
	    result.first->second->bufferSize = N_receivedCells*elementByteSize;
	    result.first->second->buffer = new char[N_receivedCells*elementByteSize];
	    result.first->second->elementSize = elementByteSize;
	 }
      }
      
      // Create MPI datatypes for sending and receiving data:
      info.N_receives = 0;
      info.N_sends    = 0;
      std::map<CellID,int> receivesPosted;
      for (typename std::map<MPI_processID,TypeCache>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
	 // Allocate arrays for MPI datatypes:
	 size_t N_recvs = 0;
	 std::map<MPI_processID,std::set<pargrid::CellID> >::const_iterator tmp = recvs.find(jt->first);
	 if (tmp != recvs.end()) N_recvs = tmp->second.size();
	 size_t N_sends = 0;
	 tmp = sends.find(jt->first);
	 if (tmp != sends.end()) N_sends = tmp->second.size();
	 
	 // Create MPI datatype for receiving all user data, associated with given 
	 // transfer ID, at once from process jt->first:
	 if (N_recvs > 0) {
	    // Allocate arrays for creating derived MPI datatype:
	    MPI_Datatype basicType;
	    RecvBuffer* rbuffer = NULL;
	    int* disps = NULL;
	    jt->second.recvs.push_back(TypeWrapper());
	    
	    switch (stencilType) {
	     case localToRemoteUpdates:
	       parGrid->getUserDatatype(userDataID,recvs[jt->first],jt->second.recvs.back().type,false);
	       break;
	     case remoteToLocalUpdates:
	       // Create MPI derived datatype that transfers user data of single cell:
	       typename std::map<DataID,RecvBuffer*>::iterator tmp = recvBuffers.find(userDataID);
               #ifndef NDEBUG
	          if (tmp == recvBuffers.end()) {
		     std::cerr << "(PARGRID) Stencil::calcTypeCache ERROR: Could not find RecvBuffer!" << std::endl;
		     exit(1);
		  }
               #endif
	       rbuffer = tmp->second;
	       MPI_Type_contiguous(rbuffer->elementSize,MPI_Type<char>(),&basicType);
	       MPI_Type_commit(&basicType);
	       
	       // Allocate memory for displacements array, and calculate displacements 
	       // into receive buffer for cells received from process jt->first:
	       disps = new int[recvs[jt->first].size()];
	       size_t counter = 0;
	       for (std::set<CellID>::const_iterator i=recvs[jt->first].begin(); i!=recvs[jt->first].end(); ++i) {
		  const CellID localID      = parGrid->getLocalID(*i);
		  std::vector<CellID>::const_iterator ptr = std::find(boundaryCells.begin(),boundaryCells.end(),localID);
                  #ifndef NDEBUG
		     if (ptr == boundaryCells.end()) {
			std::cerr << "(PARGRID) Stencil::calcTypeCache ERROR: could not find offset into recv buffer!" << std::endl;
			exit(1);
		     }
                  #endif
		  const unsigned int offset = (rbuffer->offsets)[ptr - boundaryCells.begin()];
		  
		  disps[counter] = offset + receivesPosted[*i];
		  ++counter;
		  ++receivesPosted[*i];
	       }
	       
	       #ifndef NDEBUG
	          for (size_t index=0; index<recvs[jt->first].size(); ++index) {
		     if (disps[index] >= rbuffer->bufferSize) {
			std::cerr << "(PARGRID) Stencil::calcTypeCache ERROR: calculated displacement into buffer too large!" << std::endl;
			exit(1);
		     }
		  }
               #endif
	       
	       // Create MPI datatype for receiving all data at once from process jt->second:
	       MPI_Type_create_indexed_block(recvs[jt->first].size(),1,disps,basicType,&(jt->second.recvs.back().type));
	       MPI_Type_commit(&(jt->second.recvs.back().type));
	       break;
	    }
	    ++info.N_receives;
	    delete [] disps; disps = NULL;
	 }
	  
	 // Create MPI datatype for sending all user data, associated with given
	 // transfer ID, at once to process jt->first:
	 if (N_sends > 0) {
	    jt->second.sends.push_back(TypeWrapper());
	    switch (stencilType) {
	     case localToRemoteUpdates:
	       parGrid->getUserDatatype(userDataID,sends[jt->first],jt->second.sends.back().type,false);
	       break;
	     case remoteToLocalUpdates:
	       parGrid->getUserDatatype(userDataID,sends[jt->first],jt->second.sends.back().type,true);
	       break;
	    }
	    ++info.N_sends;
	 }
      }

      // Allocate enough MPI requests:
      info.requests = new MPI_Request[info.N_receives+info.N_sends];
      return true;
   }
      
   /** Clear some of Stencil internal variables.*/
   template<class PARGRID,class C> inline
   void Stencil<PARGRID,C>::clear() {
      boundaryCells.clear();
      innerCells.clear();
      recvs.clear();
      sends.clear();
      recvCounts.clear();
   }
   
   /** Get boundary cell list for this Stencil.
    * @return Vector containing boundary cell local IDs.*/
   template<class PARGRID,class C> inline
   const std::vector<CellID>& Stencil<PARGRID,C>::getBoundaryCells() const {return boundaryCells;}
   
   /** Get inner cell list for this Stencil.
    * @return Vector containing inner cell local IDs.*/
   template<class PARGRID,class C> inline
   const std::vector<CellID>& Stencil<PARGRID,C>::getInnerCells() const {return innerCells;}

   /** Get list(s) of cells to receive from each neighbour process, ordered by their global IDs.
    * @return List(s) of cells to receive.*/
   template<class PARGRID,class C> inline
   const std::map<MPI_processID,std::set<CellID> >& Stencil<PARGRID,C>::getRecvs() const {return recvs;}
   
   /** Get offset and buffer array for given user-defined ParGrid data array. This 
    * function should only be called for StencilType::remoteToLocalUpdates.
    * @param userDataID ID number of the user-defined data array.
    * @param offsets Address of offsets array is copied to this variable.
    * @param buffer Address of receive buffer array is copied to this variable.
    * @return If true, offsets and buffer variables contain valid values.*/
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::getRemoteUpdates(DataID userDataID,unsigned int*& offsets,char*& buffer) const {
      if (stencilType != remoteToLocalUpdates) return false;
      typename std::map<DataID,RecvBuffer*>::const_iterator it = recvBuffers.find(userDataID);
      if (it == recvBuffers.end()) return false;
      offsets = it->second->offsets;
      buffer  = it->second->buffer;
      return true;
   }

   /** Get list(s) of cells sent to each neighbour process, ordered by their global IDs.
    * @return List of cells to send to each neighbour process.*/
   template<class PARGRID,class C> inline
   const std::map<MPI_processID,std::set<CellID> >& Stencil<PARGRID,C>::getSends() const {return sends;}
   
   /** Initialize Stencil.
    * @param parGrid Pointer to ParGrid.
    * @param stencilType Type of the Stencil.
    * @param receives List of received neighbours' type IDs.
    * @return If true, Stencil initialized successfully.*/
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::initialize(PARGRID& parGrid,StencilType stencilType,const std::vector<NeighbourID>& receives) {
      initialized = false;
      this->parGrid = &parGrid;
      this->stencilType = stencilType;
      this->receivedNbrTypeIDs = receives;
      sentNbrTypeIDs.reserve(receives.size());
	
      // Send stencil is the inverse of receive stencil:
      int i_off = 0;
      int j_off = 0;
      int k_off = 0;
      for (size_t i=0; i<receives.size(); ++i) {
	 parGrid.calcNeighbourOffsets(receives[i],i_off,j_off,k_off);
	 i_off *= -1;
	 j_off *= -1;
	 k_off *= -1;
	 sentNbrTypeIDs.push_back(parGrid.calcNeighbourTypeID(i_off,j_off,k_off));
      }
      
      // Sort vectors containing sent and received neighbour types:
      std::sort(sentNbrTypeIDs.begin(),sentNbrTypeIDs.end());
      std::sort(receivedNbrTypeIDs.begin(),receivedNbrTypeIDs.end());

      // Calculate requested stencil:
      initialized = true;
      if (update() == false) initialized = false;
      return initialized;
   }
   
   /** Remove transfer with given identifier from Stencil.
    * @param ID ID number of the user data array.
    * @return If true, transfer was removed successfully.*/
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::removeTransfer(DataID userDataID) {
      // Check that transfer exists:
      if (initialized == false) return false;
      #ifndef NDEBUG
         if (typeCachesUser.find(userDataID) == typeCachesUser.end()) return false;
         if (typeInfoUser.find(userDataID) == typeInfoUser.end()) return false;
      #endif

      // Erase transfer:
      typeCachesUser.erase(userDataID);
      typeInfoUser.erase(userDataID);
      
      // Erase RecvBuffer associated with this array:
      typename std::map<DataID,RecvBuffer*>::iterator jt = recvBuffers.find(userDataID);
      if (jt != recvBuffers.end()) {
	 delete jt->second;
	 jt->second = NULL;
	 recvBuffers.erase(jt);
      }
      return true;
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::startTransfer(DataID userDataID) {
      return startStaticTransfer(userDataID);
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::startStaticTransfer(DataID userDataID) {
      typename std::map<DataID,std::map<MPI_processID,TypeCache> >::iterator it = typeCachesUser.find(userDataID);
      typename std::map<DataID,TypeInfo>::iterator info = typeInfoUser.find(userDataID);
      if (it == typeCachesUser.end()) return false;
      if (info->second.started == true) return false;
      if (info->second.typeVolatile == true) calcTypeCacheUser(userDataID);
      
      // Post sends and receives:
      unsigned int counter = 0;
      MPI_Request* requests = info->second.requests;
      for (typename std::map<MPI_processID,TypeCache>::iterator proc=it->second.begin(); proc!=it->second.end(); ++proc) {
	 // Get send and receive buffers -- these are different for StencilType::remoteToLocalUpdates
	 void* rcvBuffer = NULL;
	 void* sndBuffer = NULL;
	 if (stencilType == localToRemoteUpdates) {
	    rcvBuffer = parGrid->getUserData(userDataID);
	    sndBuffer = rcvBuffer;
	 } else {
	    sndBuffer = parGrid->getUserData(userDataID);
	    typename std::map<DataID,RecvBuffer*>::iterator it = recvBuffers.find(userDataID);
	    if (it != recvBuffers.end()) rcvBuffer = it->second->buffer;
	 }
	 
	 #ifndef NDEBUG
	    if (rcvBuffer == NULL) {
	       std::cerr << "(STENCIL) ERROR: User data ID #" << userDataID << " NULL recv buffer!" << std::endl; exit(1);
	    }
	    if (sndBuffer == NULL) {
	       std::cerr << "(STENCIL) ERROR: User data ID #" << userDataID << " NULL send buffer!" << std::endl; exit(1);
	    }
         #endif
	    
	 for (size_t i=0; i<proc->second.recvs.size(); ++i) {
	    MPI_Irecv(rcvBuffer,1,proc->second.recvs[i].type,proc->first,proc->first,parGrid->getComm(),requests+counter);
	    ++counter;
	 }
	 for (size_t i=0; i<proc->second.sends.size(); ++i) {
	    MPI_Isend(sndBuffer,1,proc->second.sends[i].type,proc->first,parGrid->getRank(),parGrid->getComm(),requests+counter);
	    ++counter;
	 }
      }
      info->second.started = true;
      return true;
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::update() {
      if (initialized == false) return false;
	
      // Deallocate RecvBuffers:
      for (typename std::map<DataID,RecvBuffer*>::iterator it=recvBuffers.begin(); it!=recvBuffers.end(); ++it) {
	 delete it->second;
	 it->second = NULL;
      }
      recvBuffers.clear();
      
      // Recalculate new send & recv lists:
      bool success = calcLocalUpdateSendsAndReceives();
      
      // Recalculate MPI datatype caches (static data):
      for (typename std::map<DataID,std::map<MPI_processID,TypeCache> >::iterator it=typeCachesUser.begin(); it!=typeCachesUser.end(); ++it)
	calcTypeCacheUser(it->first);
      
      return success;
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::wait(DataID userDataID) {
      return waitStatic(userDataID);
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::wait(DataID userDataID,const std::string& name) {
      return waitStatic(userDataID,name);
   }

   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::waitStatic(DataID userDataID) {
      typename std::map<DataID,TypeInfo>::iterator info = typeInfoUser.find(userDataID);
      if (info == typeInfoUser.end()) return false;
      if (info->second.started == false) return false;
      MPI_Waitall(info->second.N_receives+info->second.N_sends,info->second.requests,MPI_STATUSES_IGNORE);
      info->second.started = false;
      return true;
   }
   
   template<class PARGRID,class C> inline
   bool Stencil<PARGRID,C>::waitStatic(DataID userDataID,const std::string& name) {
      typename std::map<DataID,TypeInfo>::iterator info = typeInfoUser.find(userDataID);
      if (info == typeInfoUser.end()) return false;
      if (info->second.started == false) return false;
      
      int flag = false;
      uint64_t waitTime      = 10000;
      uint64_t waitedTime    = 0;
      uint64_t maxWaitedTime = 5000000000;
      timespec timeSpec;
      timeSpec.tv_sec  = 0;
      timeSpec.tv_nsec = waitTime;
      do {
	 MPI_Testall(info->second.N_receives+info->second.N_sends,info->second.requests,&flag,MPI_STATUSES_IGNORE);
	 nanosleep(&timeSpec,NULL);
	 waitedTime += waitTime;
	 if (waitedTime >= maxWaitedTime) {
	    std::cerr << "Killing execution in ParGrid::Stencil::waitStatic with name '" << name << "'" << std::endl;
	    exit(1);
	 }
      } while (flag == false);

      info->second.started = false;
      return true;
   }
      
} // namespace pargrid

#endif