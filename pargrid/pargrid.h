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

#ifndef PARGRID_H
#define PARGRID_H

// Includes for standard headers:

#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <sstream>

// Includes for external package headers:

#include <mpi.h>
#include <zoltan_cpp.h>

#ifdef PROFILE
   #include "profiler.h"
#endif

// Includes for ParGrid package headers:
#include "mpiconversion.h"
#include "pargrid_definitions.h"
#include "pargrid_zoltan_callbacks.h"
#include "pargrid_userdata_static.h"
#include "pargrid_userdata_dynamic.h"
#include "pargrid_stencil.h"
#include "pargrid_datawrapper.h"

namespace pargrid {

   template<class C> class ParGrid {
    public:
      ParGrid();
      ~ParGrid();

      // ***************************************** //
      // ***** INTERFACE AIMED FOR END USERS ***** //
      // ***************************************** //
      
      bool addDataTransfer(DataID userDataID,StencilID stencilID);
      StencilID addStencil(pargrid::StencilType stencilType,const std::vector<NeighbourID>& recvNbrTypeIDs);
      template<typename T> DataID addUserData(const std::string& name,unsigned int N_elements,bool isDynamic=false);
      DataID addUserData(const std::string& name,uint64_t N_elements,const std::string& datatype,uint64_t dataSize,bool isDynamic=false);
      void calcNeighbourOffsets(NeighbourID nbrTypeID,int& i_off,int& j_off,int& k_off) const;
      NeighbourID calcNeighbourTypeID(int i_off,int j_off,int k_off) const;
      bool checkPartitioningStatus(int& counter) const;
      bool checkSuccess(const bool& value) const;
      const std::vector<CellID>& getBoundaryCells(StencilID stencilID) const;
      CellID* getCellNeighbourIDs(CellID cellID);
      std::vector<CellWeight>& getCellWeights();
      const std::vector<CellID>& getExteriorCells();
      const std::vector<CellID>& getGlobalIDs() const;
      const std::vector<MPI_processID>& getHosts() const;
      const std::vector<CellID>& getInnerCells(StencilID stencilID) const;
      const std::vector<CellID>& getInteriorCells();
      CellID getLocalID(CellID globalID) const;
      const std::vector<CellID>& getLocalIDs() const;
      const std::vector<uint32_t>& getNeighbourFlags() const;
      uint32_t getNeighbourFlags(CellID cellID) const;
      const std::set<MPI_processID>& getNeighbourProcesses() const;
      CellID getNumberOfAllCells() const;
      CellID getNumberOfLocalCells() const;
      bool getRemoteNeighbours(CellID cellID,const std::vector<NeighbourID>& nbrTypeIDs,std::vector<CellID>& nbrIDs);
      template<typename T> bool getRemoteUpdates(StencilID stencilID,DataID userDataID,unsigned int*& offsets,T*& buffer) const;
      DataID getUserDataID(const std::string& name,bool& isDynamic) const;
      template<typename T> T* getUserDataStatic(DataID userDataID);
      template<typename T> T* getUserDataStatic(const std::string& name);
      std::size_t getUserDataStaticElements(DataID userDataID) const;
      template<typename T> DataWrapper<T> getUserDataDynamic(DataID userDataID);
      template<typename T> DataWrapper<T> getUserDataDynamic(const std::string& name);
      CellID invalid() const;
      CellID invalidCellID() const;
      DataID invalidDataID() const;
      StencilID invalidStencilID() const;
      bool localCellExists(CellID cellID);
      bool removeDataTransfer(DataID userDataID,StencilID stencilID);
      bool removeStencil(StencilID stencilID);
      bool removeUserData(DataID userDataID);
      bool startNeighbourExchange(StencilID stencilID,DataID userDataID);
      bool wait(StencilID stencilID,DataID userDataID);
      bool wait(StencilID stencilID,DataID userDataID,const std::string& name);
      
      // **************************************************** //
      // ***** INTERFACE FOR MAIN PROGRAM USING PARGRID ***** //
      // **************************************************** //
      
      bool addCell(CellID cellID,const std::vector<CellID>& nbrIDs,const std::vector<NeighbourID>& nbrTypes);
      bool addCellFinished();
      bool balanceLoad();
      void barrier() const;
      void clearCellWeights();
      bool finalize();
      MPI_Comm getComm() const;
      void getDynamicUserDataInfo(std::vector<DataID>& dataIDs,std::vector<std::string>& names,
				  std::vector<ArraySizetype*>& sizeArrays,std::vector<std::string>& datatypes,
				  std::vector<unsigned int>& byteSizes,std::vector<const char**>& pointers) const;
      bool getInitialized() const;
      std::size_t getNumberOfReceives(StencilID stencilID,MPI_processID hostID) const;
      std::size_t getNumberOfSends(StencilID stencilID,MPI_processID hostsID) const;
      MPI_processID getProcesses() const;
      MPI_processID getRank() const;
      const std::map<MPI_processID,std::set<CellID> >& getReceives(StencilID stencilID) const;
      const std::map<MPI_processID,std::set<CellID> >& getSends(StencilID stencilID) const;
      void getStaticUserDataInfo(std::vector<DataID>& dataIDs,std::vector<std::string>& names,std::vector<unsigned int>& elements,
				 std::vector<std::string>& datatypes,std::vector<unsigned int>& byteSizes,std::vector<const char*>& pointers) const;
      char* getUserData(DataID userDataID);
      char* getUserData(const std::string& name);
      unsigned int getUserDataElementSize(DataID userDataID) const;
      bool getUserDatatype(DataID transferID,const std::set<CellID>& globalIDs,MPI_Datatype& datatype,bool reverseStencil);
      bool initialize(MPI_Comm comm,const std::vector<std::map<InputParameter,std::string> >& parameters);
      bool setPartitioningMode(PartitioningMode pm);

      // ***** ZOLTAN CALLBACK FUNCTION DECLARATIONS ***** //

      void getAllCellCoordinates(int N_globalEntries,int N_localEntries,int N_cellIDs,ZOLTAN_ID_PTR globalID,
				 ZOLTAN_ID_PTR localID,int N_coords,double* geometryData,int* rcode);
      void getCellCoordinates(int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
			      ZOLTAN_ID_PTR localID,double* geometryData,int* rcode);
      void getAllCellEdges(int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			   ZOLTAN_ID_PTR localIDs,int* N_edges,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			   int N_weights,CellWeight* edgeWeights,int* rcode);
      void getCellEdges(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			int N_weights,CellWeight* weight,int* rcode);
      void getHierarchicalParameters(int level,Zoltan_Struct* zs,int* rcode);
      int getHierarchicalPartNumber(int level,int* rcode);
      void getHyperedges(int N_globalIDs,int N_vtxedges,int N_pins,int format,
			 ZOLTAN_ID_PTR vtxedge_GID,int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode);
      void getHyperedgeWeights(int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
			       ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,CellWeight* edgeWeights,int* rcode);
      void getLocalCellList(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
			    ZOLTAN_ID_PTR localIDs,int N_weights,CellWeight* cellWeights,int* rcode);
      int getMeshDimension(int* rcode);
      void getNumberOfAllEdges(int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			      ZOLTAN_ID_PTR localIDs,int* N_edges,int* rcode);
      int getNumberOfEdges(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			   ZOLTAN_ID_PTR localID,int* rcode);
      int getNumberOfHierarchicalLevels(int* rcode);
      void getNumberOfHyperedges(int* N_lists,int* N_pins,int* format,int* rcode);
      void getNumberOfHyperedgeWeights(int* N_edges,int* rcode);
      int getNumberOfLocalCells(int* rcode);

    private:
      ParGrid(const ParGrid<C>& p);

      bool recalculateExteriorCells;                                      /**< If true, exteriorCells needs to be recalculated when 
									   * it is requested.*/
      bool recalculateInteriorCells;                                      /**< If true, interiorCells needs to be recalculated 
									   * when it is requested.*/
      std::vector<CellID> exteriorCells;                                  /**< List of exterior cells on this process, i.e. cells that lie on 
									   * the boundary of the simulation domain. These cells have one or more 
									   * missing neighbours. This list is recalculated only when needed.*/
      std::vector<CellID> interiorCells;                                  /**< List of interior cells on this process, i.e. cells with 
									   * zero missing neighbours. This list is recalculated only when needed.*/      
      CellWeight cellWeight;                                              /**< Cell weight scale, used to calculate cell weights for Zoltan.*/
      std::vector<CellWeight> cellWeights;                                /**< Computational load of each local cell.*/
      bool cellWeightsUsed;                                               /**< If true, cell weights are calculated.*/
      std::vector<CellID> globalIDs;                                      /**< Global IDs of cells stored in this process.*/
      std::vector<CellID> localIDs;                                       /**< Local IDs of cells stored in this process. For remote cells 
									   * this vector contains the local ID of that cell in the remote process.*/
      CellID N_localCells;                                                /**< Number of local cells on this process.*/
      CellID N_totalCells;                                                /**< Total number of cells (local+buffered) on this process.*/
      std::map<CellID,CellID> global2LocalMap;                            /**< Global-to-local ID mapping.*/
      std::vector<MPI_processID> hosts;                                   /**< Host process for each cell in vector cells.*/
      std::vector<CellID> cellNeighbours;                                 /**< Local IDs of neighbours of cells, only valid for local cells
									   as remote cells (that are cached on this process) neighbours do 
									   not have local copies.*/
      std::vector<uint32_t> neighbourFlags;                               /**< Neighbour existence flags for local cells.*/
      std::map<DataID,UserDataDynamic<ParGrid<C> >*> userDataDynamic;     /**< Arrays containing dynamically allocated user data.*/
      std::map<DataID,UserDataStatic<ParGrid<C> >*> userDataStatic;       /**< Arrays containing static-size user data.*/
      std::set<unsigned int> userDataHoles;                               /**< List of holes in userData.*/

      MPI_Comm comm;                                                      /**< MPI communicator used by ParGrid.*/
      CellWeight edgeWeight;                                              /**< Edge weight scale, used to calculate edge weights for Zoltan.*/
      bool edgeWeightsUsed;                                               /**< If true, edge weights are calculated.*/
      bool initialized;                                                   /**< If true, ParGrid initialized correctly and is ready for use.*/
      std::vector<std::list<std::pair<std::string,std::string> > > 
	                                         loadBalancingParameters; /**< LB parameters for Zoltan for each hierarchical level.*/
      MPI_processID myrank;                                               /**< MPI rank of this process in communicator comm.*/
      std::set<MPI_processID> nbrProcesses;                               /**< MPI ranks of neighbour processes. A process is considered 
									   * to be a neighbour if it has any of this processes local 
									   * cells' neighbours. Calculated in balanceLoad.*/
      MPI_processID N_processes;                                          /**< Number of MPI processes in communicator comm.*/
      int partitioningCounter;                                            /**< Counter that is increased every time the mesh is repartitioned.*/
      std::vector<MPI_Request> recvRequests;                              /**< List of MPI requests that are used for receiving data.*/
      std::vector<MPI_Request> sendRequests;                              /**< List of MPI requests that are used for sending data.*/
      std::map<StencilID,Stencil<ParGrid<C>,C > > stencils;                 /**< Transfer stencils, identified by their ID numbers.*/
      Zoltan* zoltan;                                                     /**< Pointer to Zoltan.*/
      
      #ifdef PROFILE
         int profZoltanLB;
         int profParGridLB;
         int profMPI;
         int profTotalLB;
         int profStaticData;
         int profDynamicData;
         int profStencilRecalc;
         int profZoltanCB;
      #endif
      
      bool checkInternalStructures() const;
      void invalidate();
      void startMetadataRepartitioning(std::vector<MPI_Request>& nbrRecvRequests,
				       std::map<MPI_processID,std::vector<int> >& importDisplacements,
				       std::vector<MPI_Request>& nbrSendRequests,
				       std::map<MPI_processID,std::vector<int> >& exportDisplacements,
				       std::vector<CellID>& newCellNeighbours);
      bool repartitionUserDataStatic(size_t N_cells,CellID newLocalsBegin,int N_import,int* importProcesses,
				     const std::map<MPI_processID,std::vector<int> >& importDisplacements,
				     int N_export,int* exportProcesses,ZOLTAN_ID_PTR exportLocalIDs,
				     const std::map<MPI_processID,std::vector<int> >& exportDisplacements);

      bool repartitionUserDataDynamic(size_t N_cells,CellID newLocalsBegin,int N_import,int* importProcesses,
				      const std::map<MPI_processID,std::vector<int> >& importDisplacements,
				      int N_export,int* exportProcesses,ZOLTAN_ID_PTR exportLocalIDs,
				      const std::map<MPI_processID,std::vector<int> >& exportDisplacements);
      bool syncCellHosts();
      bool syncLocalIDs();
      void waitMetadataRepartitioning(std::vector<MPI_Request>& recvRequests,std::vector<MPI_Request>& sendRequests);
   };


   template<typename T> inline
   std::string getDatatype() {return "unknown";}

   template<> inline
   std::string getDatatype<bool>() {return "int";}
   template<> inline
   std::string getDatatype<char>() {return "int";}
     
   template<> inline
   std::string getDatatype<int8_t>() {return "int";}
   template<> inline
   std::string getDatatype<int16_t>() {return "int";}
   template<> inline
   std::string getDatatype<int32_t>() {return "int";}
   template<> inline
   std::string getDatatype<int64_t>() {return "int";}
      
   template<> inline
   std::string getDatatype<uint8_t>() {return "uint";}
   template<> inline
   std::string getDatatype<uint16_t>() {return "uint";}
   template<> inline
   std::string getDatatype<uint32_t>() {return "uint";}
   template<> inline
   std::string getDatatype<uint64_t>() {return "uint";}
      
   template<> inline
   std::string getDatatype<float>() {return "float";}
   template<> inline
   std::string getDatatype<double>() {return "float";}
   template<> inline
   std::string getDatatype<long double>() {return "float";}

   
   /** Constructor for ParGrid. Initializes Zoltan. Note that MPI_Init must 
    * have been called prior to calling ParGrid constructor.
    * @param hierarchicalLevels Number of hierarchical partitioning levels to use.
    * @param parameters Load balancing parameters for all hierarchical levels.
    * The parameters for each level are given in a list, whose contents are pairs 
    * formed from parameter types and their string values. These lists themselves 
    * are packed into a list, whose first item (list) is used for hierarchical level
    * 0, second item for hierarchical level 1, and so forth.*/   
   template<class C> inline
   ParGrid<C>::ParGrid(): initialized(false) {
      recalculateInteriorCells = true;
      recalculateExteriorCells = true;
      N_localCells = 0;
      N_totalCells = 0;
      partitioningCounter = 0;
      
      #ifdef PROFILE
         profZoltanLB  = -1;
         profZoltanCB  = -1;
         profParGridLB = -1;
         profMPI       = -1;
         profTotalLB   = -1;
         profStaticData  = -1;
         profDynamicData = -1;
         profStencilRecalc = -1;
      #endif
   }
   
   /** Destructor for ParGrid. Deallocates are user data in cells, i.e. destructor 
    * is called for ParCell::userData for each existing cell. Note that ParGrid 
    * destructor does not call MPI_Finalize.*/
   template<class C> inline
   ParGrid<C>::~ParGrid() { }
   
   /** Add a new cell to ParGrid on this process. This function should only 
    * be called when adding initial cells to the grid, i.e. after calling 
    * ParGrid::initialize but before calling initialLoadBalance.
    * @return If true, the cell was inserted successfully.*/
   template<class C> inline
   bool ParGrid<C>::addCell(CellID cellID,const std::vector<CellID>& nbrIDs,const std::vector<NeighbourID>& nbrTypes) {
      if (getInitialized() == false) return false;
      bool success = true;
      
      // Check that the cell doesn't already exist:
      if (global2LocalMap.find(cellID) != global2LocalMap.end()) return false;
      global2LocalMap[cellID] = N_localCells;
      hosts.push_back(getRank());
      globalIDs.push_back(cellID);
      
      // Copy cell's neighbours and increase reference count to cell's neighbours:
      const size_t offset = cellNeighbours.size();
      const size_t flagOffset = neighbourFlags.size();
      cellNeighbours.insert(cellNeighbours.end(),N_neighbours,invalid());
      neighbourFlags.insert(neighbourFlags.end(),0);
            
      uint32_t nbrFlag = (1 << calcNeighbourTypeID(0,0,0));
      for (size_t n=0; n<nbrIDs.size(); ++n) {
	 //if (nbrIDs[n] == invalid()) success = false;
	 if (nbrIDs[n] == invalid()) continue;
	 if (nbrTypes[n] > calcNeighbourTypeID(1,1,1)) success = false;
	 
	 cellNeighbours[offset+nbrTypes[n]]     = nbrIDs[n];
	 nbrFlag = (nbrFlag | (1 << nbrTypes[n]));
      }
      neighbourFlags[flagOffset] = nbrFlag;
      cellNeighbours[offset+calcNeighbourTypeID(0,0,0)] = cellID;
      
      ++N_localCells;
      ++N_totalCells;
      return success;
   }

   /** Tell ParGrid that all new cells have been added. 
    * This function contains massive amount of MPI transfer, as 
    * each process must give every other process information on the 
    * cells it has.
    * @return If true, information was shared successfully with other processes.*/
   template<class C> inline
   bool ParGrid<C>::addCellFinished() {
      if (getInitialized() == false) return false;
      
      // Insert remote cells and convert global IDs in neighbour lists into local IDs:
      for (size_t i=0; i<N_localCells; ++i) {
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrGID = cellNeighbours[i*N_neighbours+n];
	    if (nbrGID == invalid()) continue;

	    std::map<CellID,CellID>::const_iterator it = global2LocalMap.find(nbrGID);
	    if (it != global2LocalMap.end()) {
	       // Neighbour is a local cell, replace global ID with local ID:
	       const CellID nbrLID = it->second;
	       cellNeighbours[i*N_neighbours+n] = nbrLID;
	    } else {
	       // Neighbour is a remote cell. Insert a local copy of the neighbour, 
	       // and replace global ID with the new local ID:
	       ++N_totalCells;
	       hosts.push_back(MPI_PROC_NULL);
	       globalIDs.push_back(nbrGID);
	       global2LocalMap[nbrGID] = N_totalCells-1;
	       cellNeighbours[i*N_neighbours+n] = N_totalCells-1;
	    }
	 }
      }
      
      // Update remote neighbour hosts:
      if (syncCellHosts() == false) {
	 std::cerr << "PARGRID FATAL ERROR: sync cell hosts failed!" << std::endl;
	 return false;
      }
      
      // Calculate neighbour processes:
      nbrProcesses.clear();
      for (size_t i=N_localCells; i<hosts.size(); ++i) {
	 nbrProcesses.insert(hosts[i]);
      }

      #ifndef NDEBUG
         // Check that data from user is ok:
         int successSum = 0;
         int mySuccess = 0;
         if (checkInternalStructures() == false) ++mySuccess;      
         MPI_Allreduce(&mySuccess,&successSum,1,MPI_Type<int>(),MPI_SUM,comm);
         if (successSum > 0) return false;
      #endif
      
      // Add default stencil (all neighbours exchange data):
      std::vector<NeighbourID> nbrTypeIDs(27);
      for (size_t i=0; i<27; ++i) {
	 if (i == 13) continue;
	 nbrTypeIDs[i] = i;
      }
      stencils[pargrid::DEFAULT_STENCIL].initialize(*this,localToRemoteUpdates,nbrTypeIDs);

      // Invalidate all internal variables that depend on cell partitioning:
      invalidate();
      return true;
   }

   /** Create a new Stencil for MPI transfers.
    * @param stencilType Type of Stencil, one of the values defined in enum StencilType.
    * @param recvNbrTypeIDs Neighbour type IDs that identify the neighbours whom to receive data from.
    * @return If greater than zero the Stencil was created successfully.*/
   template<class C> inline
   StencilID ParGrid<C>::addStencil(pargrid::StencilType stencilType,const std::vector<NeighbourID>& recvNbrTypeIDs) {
      int currentSize = stencils.size();
      if (stencils[currentSize].initialize(*this,stencilType,recvNbrTypeIDs) == false) {
	 stencils.erase(currentSize);
	 currentSize = -1;
      }
      return currentSize;
   }
   
   /** Add a parallel data array to ParGrid. Basic datatype, e.g. int or double, is 
    * given as template parameter T. You still need to associate the array with a transfer 
    * if you want to synchronize neighbour data.
    * @param name Unique name for the array. Array insertion will fail if an array with 
    * the same name already exists.
    * @param N_elements Number of basic datatypes per parallel cell. For example, in 
    * order to get an array of size five per cell, pass value five here. Note: this parameter 
    * is ignored for dynamically allocated data (isDynamic==true).
    * @param isDynamic If true, user data is dynamically allocated (=number of element per cell varies)
    * instead of being static (=same number of elements per cell).
    * @return ID number of the new data array or invalidDataID() if array was not created.
    * @see addUserDataTransfer.*/
   template<class C> template<typename T> inline
   DataID ParGrid<C>::addUserData(const std::string& name,unsigned int N_elements,bool isDynamic) {
      // Check that a user data array (static or dynamic) with the given name doesn't already exist.
      // If it does return the DataID of the existing array:
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator 
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second->getName() == name) return it->first;
      }
      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::const_iterator
	   it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second->getName() == name) return it->first;
      }
  
      // Assign an ID number to the array:
      DataID userDataID = invalidDataID();
      if (userDataHoles.size() > 0) {
	 userDataID = *userDataHoles.begin();
	 userDataHoles.erase(userDataID);
      } else {
	 userDataID = userDataStatic.size() + userDataDynamic.size();
      }

      // Initialize the new user data array:
      if (isDynamic == false) {
	 if (userDataStatic.find(userDataID) != userDataStatic.end()) {
	    std::cerr << "(PARGRID) ERROR: user data array with ID " << userDataID << " already exists" << std::endl;
	    exit(1);
	 }
	 
	 userDataStatic[userDataID] = new UserDataStatic<ParGrid<C> >();
	 userDataStatic[userDataID]->initialize(this,name,N_totalCells,N_elements,sizeof(T),pargrid::getDatatype<T>());
      } else {
	 if (userDataDynamic.find(userDataID) != userDataDynamic.end()) {
	    std::cerr << "(PARGRID) ERROR: user data array with ID " << userDataID << " already exists" << std::endl;
	    exit(1);
	 }

	 userDataDynamic[userDataID] = new UserDataDynamic<ParGrid<C> >();
	 userDataDynamic[userDataID]->initialize(this,name,N_totalCells,sizeof(T),pargrid::getDatatype<T>());
      }
      return userDataID;
   }
   
   template<class C> inline
   DataID ParGrid<C>::addUserData(const std::string& name,uint64_t N_elements,const std::string& datatype,uint64_t dataSize,bool isDynamic) {
      // Check that a user data array (static or dynamic) with the given name doesn't already exist.
      // If it does return the DataID of the existing array:
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second->getName() == name) return it->first;
      }
      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::const_iterator
	   it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second->getName() == name) return it->first;
      }
      
      // Assign a unique DataID to the array:
      DataID userDataID = invalidDataID();
      if (userDataHoles.size() > 0) {
	 userDataID = *userDataHoles.begin();
	 userDataHoles.erase(userDataID);
      } else {
	 userDataID = userDataStatic.size() + userDataDynamic.size();
      }
      
      // Initialize the new user data array:
      if (isDynamic == false) {
	 userDataStatic[userDataID] = new UserDataStatic<ParGrid<C> >();
	 userDataStatic[userDataID]->initialize(this,name,N_totalCells,N_elements,dataSize,datatype);
      } else {
	 userDataDynamic[userDataID] = new UserDataDynamic<ParGrid<C> >();
	 userDataDynamic[userDataID]->initialize(this,name,N_totalCells,dataSize,datatype);
      }
      return userDataID;
   }
   
   /** Add a new MPI transfer for a user-defined data array.
    * @param userDataID ID number of the user data array, the value returned by addUserData function.
    * @param stencilID ID number of the transfer stencil. For a user-defined stencil this is the value 
    * returned by addStencil function. For the default stencil a value 0.
    * @param transferID A unique ID number for the transfer.
    * @return If true, a new transfer was added successfully.*/
   template<class C> inline
   bool ParGrid<C>::addDataTransfer(DataID userDataID,StencilID stencilID) {
      if (getInitialized() == false) return false;
      if (userDataStatic.find(userDataID) != userDataStatic.end()) {
	 typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator it = stencils.find(stencilID);
	 if (it == stencils.end()) return false;
	 return it->second.addUserDataTransfer(userDataID,false);
      } else {
	 std::cerr << "(PARGRID) ERROR: Could not find user data with ID#" << userDataID << " in addDataTransfer!" << std::endl;
	 return false;
      }
   }
   
   /** Balance load based on criteria given in ParGrid::initialize. Load balancing 
    * invalidates lists of cells etc. stored outside of ParGrid.
    * @return If true, mesh was repartitioned successfully.*/
   template<class C> inline
   bool ParGrid<C>::balanceLoad() {
      #ifdef PROFILE
         profile::start("ParGrid load balance",profTotalLB);
      #endif
      bool success = true;
      
      // ********************************************************
      // ***** STEP 1: REQUEST NEW PARTITIONING FROM ZOLTAN *****
      // ********************************************************
      
      // Request load balance from Zoltan, and get cells which should be imported and exported.
      // NOTE that import/export lists may contain cells that already are on this process, at
      // least with RANDOM load balancing method!
      #ifdef PROFILE
         profile::start("Zoltan LB",profZoltanLB);
      #endif
      int changes,N_globalIDs,N_localIDs,N_import,N_export;
      int* importProcesses = NULL;
      int* importParts     = NULL;
      int* exportProcesses = NULL;
      int* exportParts     = NULL;
      ZOLTAN_ID_PTR importGlobalIDs = NULL;
      ZOLTAN_ID_PTR importLocalIDs  = NULL;
      ZOLTAN_ID_PTR exportGlobalIDs = NULL;
      ZOLTAN_ID_PTR exportLocalIDs  = NULL;
      if (zoltan->LB_Partition(changes,N_globalIDs,N_localIDs,N_import,importGlobalIDs,importLocalIDs,
			       importProcesses,importParts,N_export,exportGlobalIDs,exportLocalIDs,
			       exportProcesses,exportParts) != ZOLTAN_OK) {
	 std::cerr << "ParGrid FATAL ERROR: Zoltan failed on load balancing!" << std::endl << std::flush;
	 zoltan->LB_Free_Part(&importGlobalIDs,&importLocalIDs,&importProcesses,&importParts);
	 zoltan->LB_Free_Part(&exportGlobalIDs,&exportLocalIDs,&exportProcesses,&exportParts);
	 success = false;
	 #ifdef PROFILE
	    profile::stop();
	    profile::stop();
	 #endif
	 return success;
      }
      #ifdef PROFILE
         profile::stop();
         profile::start("ParGrid LB",profParGridLB);
      #endif

      // Allocate memory for new cell array and start to receive imported cells.
      // Note. N_newLocalCells is correct even if N_import and N_export contain 
      // transfers from process A to process A.
      const double totalToLocalRatio = (1.0*N_totalCells) / (N_localCells+1.0e-20);
      const CellID N_newLocalCells = N_localCells + N_import - N_export;
      const CellID newCapacity = static_cast<CellID>(std::ceil((totalToLocalRatio+0.1)*N_newLocalCells));

      // Calculate index into newCells where imported cells from other processes are 
      // inserted. N_realExports != N_export because export list may contain 
      // exports from process A to process A, i.e. not all cells are migrated.
      int N_realExports = 0;
      for (int i=0; i<N_export; ++i) if (exportProcesses[i] != getRank()) ++N_realExports;
      const CellID newLocalsBegin = N_localCells - N_realExports;

      // Calculate displacement arrays for imported and exported static-size data. 
      // "Static-size" here means that all cells store the same amount of data, i.e. 
      // the data can be transferred with indexed blocks. Import/export displacement 
      // arrays created here are used to transfer static user data arrays and 
      // cell neighbour IDs:
      std::map<MPI_processID,std::vector<int> > exportDisplacements;
      std::map<MPI_processID,std::vector<int> > importDisplacements;
      size_t counter = newLocalsBegin;
      for (int i=0; i<N_import; ++i) {
	 if (importProcesses[i] == getRank()) continue;
	 importDisplacements[importProcesses[i]].push_back(counter);
	 ++counter;
      }
      for (int i=0; i<N_export; ++i) {
	 if (exportProcesses[i] == getRank()) continue;
	 exportDisplacements[exportProcesses[i]].push_back(exportLocalIDs[i]);
      }
      
      //std::vector<ParCell<C> > newCells;
      std::vector<MPI_processID> newHosts;
      std::map<CellID,CellID> newGlobal2LocalMap;
      std::vector<CellID> newGlobalIDs;
      std::vector<CellID> newCellNeighbours;
      std::vector<uint32_t> newNeighbourFlags;      
      newHosts.reserve(newCapacity);
      newGlobalIDs.reserve(newCapacity);
      newCellNeighbours.resize(N_neighbours*N_newLocalCells);
      newNeighbourFlags.reserve(N_newLocalCells);

      // Count the number of cells exported per neighbouring process, and number of cells 
      // imported per neighbouring process. All data needs to be sent (and received) with a 
      // single MPI call per remote neighbour, otherwise unexpected message buffers may run out!
      // Note: importsPerProcess and exportsPerProcess may contain transfers from process A to process A.
      std::map<MPI_processID,unsigned int> exportsPerProcess;
      std::map<MPI_processID,unsigned int> importsPerProcess;
      for (int i=0; i<N_import; ++i) ++importsPerProcess[importProcesses[i]];
      for (int i=0; i<N_export; ++i) ++exportsPerProcess[exportProcesses[i]];

      // Erase imports/exports from process A to process A:
      exportsPerProcess.erase(getRank());
      importsPerProcess.erase(getRank());

      // Swap all cells' neighbour local IDs to global IDs:
      for (CellID c=0; c<N_localCells; ++c) {
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrLID = cellNeighbours[c*N_neighbours+n];
	    if (nbrLID == invalid()) continue;
	    cellNeighbours[c*N_neighbours+n] = globalIDs[nbrLID];
	 }
      }       

      std::vector<MPI_Request> nbrRecvRequests;
      std::vector<MPI_Request> nbrSendRequests;
      startMetadataRepartitioning(nbrRecvRequests,importDisplacements,nbrSendRequests,exportDisplacements,newCellNeighbours);

      // Set new host for all exported cells, it should not matter 
      // if cells are exported from process A to process A here:
      for (int c=0; c<N_export; ++c) {
	 hosts[exportLocalIDs[c]] = exportProcesses[c];
      }

      // Copy non-exported local cells to newCells:
      counter = 0;
      for (CellID c=0; c<N_localCells; ++c) {
	 if (hosts[c] != getRank()) continue;
	 newGlobalIDs.push_back(globalIDs[c]);
	 newHosts.push_back(getRank());
	 newGlobal2LocalMap[globalIDs[c]] = counter;
	 
	 for (size_t n=0; n<N_neighbours; ++n) {
	    newCellNeighbours[counter*N_neighbours+n] = cellNeighbours[c*N_neighbours+n];
	 }
	 newNeighbourFlags[counter] = neighbourFlags[c];
	 ++counter;
      }

      #ifndef NDEBUG
         if (counter*N_neighbours > newCellNeighbours.size()) {
	    std::cerr << "(PARGRID) ERROR: newCellNeighbours array out of bounds!" << std::endl;
	    exit(1);
	 }
      #endif
      
      // Set host, global ID, and global2Local entries for imported cells:
      counter = 0;
      for (int i=0; i<N_import; ++i) {
	 if (importProcesses[i] == getRank()) continue;
	 
	 newHosts.push_back(getRank());
	 newGlobalIDs.push_back(importGlobalIDs[i]);
	 newGlobal2LocalMap[importGlobalIDs[i]] = newLocalsBegin + counter;
	 ++counter;
      }

      // ************************************************************ //
      // ***** INSERT REMOTE CELLS AND UPDATE REMOTE CELL HOSTS ***** //
      // ************************************************************ //

      waitMetadataRepartitioning(nbrRecvRequests,nbrSendRequests);
      
      // Iterate over new cell list and replace neighbour GIDs with LIDs.
      // Also calculate neighbourFlags for all local cells:
      for (CellID c=0; c<N_newLocalCells; ++c) {
	 uint32_t nbrFlag = (1 << calcNeighbourTypeID(0,0,0));
	 
	 std::map<CellID,CellID>::const_iterator it;
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrGID = newCellNeighbours[c*N_neighbours+n];
	    if (nbrGID == invalid()) continue;
	    
	    CellID nbrLID = invalid();
	    it = newGlobal2LocalMap.find(nbrGID);
	    if (it != newGlobal2LocalMap.end()) {
	       nbrLID = it->second;
	    } else {
	       // Insert remote cell to newCells:
	       nbrLID = newHosts.size();
	       newGlobalIDs.push_back(nbrGID);
	       newGlobal2LocalMap[nbrGID] = nbrLID;
	       newHosts.push_back(MPI_PROC_NULL);
	       
	       // If remote cell already exists on this process, copy host ID from hosts:
	       it = global2LocalMap.find(nbrGID);
	       if (it != global2LocalMap.end()) {
		  const CellID nbrOldLID = it->second;
		  newHosts[nbrLID] = hosts[nbrOldLID];
	       }
	    }

	    newCellNeighbours[c*N_neighbours+n] = nbrLID;
	    nbrFlag = (nbrFlag | (1 << n));
	 }
	 newNeighbourFlags[c] = nbrFlag;
      }

      // Repartition user-defined data arrays:
      #ifdef PROFILE
         profile::start("static user data",profStaticData);
      #endif
      repartitionUserDataStatic(newHosts.size(),newLocalsBegin,N_import,importProcesses,importDisplacements,
				N_export,exportProcesses,exportLocalIDs,exportDisplacements);
      #ifdef PROFILE
         profile::stop();
         profile::start("dynamic user data",profDynamicData);
      #endif

      repartitionUserDataDynamic(newHosts.size(),newLocalsBegin,N_import,importProcesses,importDisplacements,
				 N_export,exportProcesses,exportLocalIDs,exportDisplacements);
      #ifdef PROFILE
         profile::stop();
      #endif
      
      // First host update pass:
      // Send list of exported cells and their new hosts to every remote process:
      int* neighbourChanges                      = new int[nbrProcesses.size()];              // Number of exports nbrs. process has
      ZOLTAN_ID_TYPE** neighbourMigratingCellIDs = new ZOLTAN_ID_TYPE* [nbrProcesses.size()]; // Global IDs of cells nbr. process is exporting
      MPI_processID** neighbourMigratingHosts    = new MPI_processID* [nbrProcesses.size()];  // New hosts for cells nbr. process exports
      for (size_t i=0; i<nbrProcesses.size(); ++i) {
	 neighbourMigratingCellIDs[i] = NULL;
	 neighbourMigratingHosts[i]   = NULL;
      }

      // Send the number of cells this process is exporting to all neighbouring process,
      // and receive the number of exported cells per neighbouring process:
      if (recvRequests.size() < 2*nbrProcesses.size()) recvRequests.resize(2*nbrProcesses.size());
      if (sendRequests.size() < 2*nbrProcesses.size()) sendRequests.resize(2*nbrProcesses.size());
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=nbrProcesses.begin(); it!=nbrProcesses.end(); ++it) {
	 MPI_Irecv(&(neighbourChanges[counter]),1,MPI_INT,*it,   *it,comm,&(recvRequests[counter]));
	 MPI_Isend(&N_export,                   1,MPI_INT,*it,myrank,comm,&(sendRequests[counter]));
	 ++counter;
      }
      
      // Wait for information to arrive from neighbours:
      #ifdef PROFILE
         profile::start("MPI Waits",profMPI);
      #endif
      MPI_Waitall(nbrProcesses.size(),&(recvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(nbrProcesses.size(),&(sendRequests[0]),MPI_STATUSES_IGNORE);
      #ifdef PROFILE
         profile::stop();
      #endif
      
      // Allocate arrays for receiving migrating cell IDs and new hosts
      // from neighbouring processes. Exchange data with neighbours:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=nbrProcesses.begin(); it!=nbrProcesses.end(); ++it) {
	 neighbourMigratingCellIDs[counter] = new ZOLTAN_ID_TYPE[neighbourChanges[counter]];
	 neighbourMigratingHosts[counter]   = new MPI_processID[neighbourChanges[counter]];
	 
	 MPI_Irecv(neighbourMigratingCellIDs[counter],neighbourChanges[counter],MPI_Type<ZOLTAN_ID_TYPE>(),*it,*it,comm,&(recvRequests[2*counter+0]));
	 MPI_Irecv(neighbourMigratingHosts[counter]  ,neighbourChanges[counter],MPI_Type<MPI_processID>() ,*it,*it,comm,&(recvRequests[2*counter+1]));
	 MPI_Isend(exportGlobalIDs,N_export,MPI_Type<ZOLTAN_ID_TYPE>(),*it,myrank,comm,&(sendRequests[2*counter+0]));
	 MPI_Isend(exportProcesses,N_export,MPI_Type<MPI_processID>() ,*it,myrank,comm,&(sendRequests[2*counter+1]));
	 ++counter;
      }
      #ifndef NDEBUG
         if (2*counter > recvRequests.size() || 2*counter > sendRequests.size()) {
	    std::cerr << "(PARGRID) ERROR: 2*counter " << 2*counter << " recvRequests.size() " << recvRequests.size();
	    std::cerr << " sendRequests.size() " << sendRequests.size() << std::endl;
	    exit(1);
	 }
      #endif

      #ifdef PROFILE
         profile::start("MPI Waits",profMPI);
      #endif
      MPI_Waitall(2*nbrProcesses.size(),&(recvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(2*nbrProcesses.size(),&(sendRequests[0]),MPI_STATUSES_IGNORE);
      #ifdef PROFILE
         profile::stop();
      #endif
      
      // Update hosts, this should work even if export list from process A has exports to process A:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=nbrProcesses.begin(); it!=nbrProcesses.end(); ++it) {
	 for (int i=0; i<neighbourChanges[counter]; ++i) {
	    // Update entry in hosts (needed for second host update pass below):
	    std::map<CellID,CellID>::const_iterator it = global2LocalMap.find(neighbourMigratingCellIDs[counter][i]);
	    if (it != global2LocalMap.end()) hosts[it->second] = neighbourMigratingHosts[counter][i];
	    
	    // Update entry in newHosts:
	    it = newGlobal2LocalMap.find(neighbourMigratingCellIDs[counter][i]);
	    if (it == newGlobal2LocalMap.end()) continue;
	    
	    // Skip cells imported to this process:
	    const CellID nbrLID = it->second;
	    if (nbrLID < N_newLocalCells) continue;	    
	    newHosts[nbrLID] = neighbourMigratingHosts[counter][i];
	 }
	 ++counter;
      }
      
      // Deallocate memory:
      delete [] neighbourChanges; neighbourChanges = NULL;
      for (size_t i=0; i<nbrProcesses.size(); ++i) {
	 delete neighbourMigratingCellIDs[i]; neighbourMigratingCellIDs[i] = NULL;
	 delete neighbourMigratingHosts[i]; neighbourMigratingHosts[i] = NULL;
      }
      delete [] neighbourMigratingCellIDs; neighbourMigratingCellIDs = NULL;
      delete [] neighbourMigratingHosts; neighbourMigratingHosts = NULL;
      
      // Calculate unique import / export processes:
      std::set<MPI_processID> importProcs;
      std::set<MPI_processID> exportProcs;
      for (int i=0; i<N_import; ++i) importProcs.insert(importProcesses[i]);
      for (int i=0; i<N_export; ++i) exportProcs.insert(exportProcesses[i]);
      importProcs.erase(getRank());
      exportProcs.erase(getRank());
      
      // Second host update pass.
      // Go over exported cells' neighbours. If neighbour is not in exported cell's 
      // new host, add an entry to hostUpdates:
      std::map<MPI_processID,std::set<std::pair<CellID,MPI_processID> > > hostUpdates;
      for (int i=0; i<N_export; ++i) {
	 if (exportProcesses[i] == getRank()) continue;

	 const CellID exportLID = exportLocalIDs[i];
	 const MPI_processID newHost = hosts[exportLID];
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrGID = cellNeighbours[exportLID*N_neighbours+n];
	    if (nbrGID == invalid()) continue;
	    std::map<CellID,CellID>::const_iterator it=global2LocalMap.find(nbrGID);
	    #ifndef NDEBUG
	       if (it == global2LocalMap.end()) {
		  std::cerr << "(PARGRID) ERROR: Nbr cell with GID " << nbrGID << " was not found from global2LocalMap" << std::endl;
		  exit(1);
	       }
	    #endif
	    const MPI_processID nbrHost = hosts[it->second];
	    if (nbrHost != newHost) hostUpdates[newHost].insert(std::make_pair(nbrGID,nbrHost));
	 }
      }
      
      if (recvRequests.size() < 2*importProcs.size()) recvRequests.resize(2*importProcs.size());
      if (sendRequests.size() < 2*exportProcs.size()) sendRequests.resize(2*exportProcs.size());
      
      // Recv hostUpdates list from each process in importProcs:
      counter = 0;
      size_t* incomingUpdates = new size_t[importProcs.size()];
      for (std::set<MPI_processID>::const_iterator it=importProcs.begin(); it!=importProcs.end(); ++it) {
	 MPI_Irecv(incomingUpdates+counter,1,MPI_Type<size_t>(),*it,*it,comm,&(recvRequests[counter]));
	 ++counter;
      }
      
      counter = 0;
      size_t* outgoingUpdates = new size_t[exportProcs.size()];
      for (std::set<MPI_processID>::const_iterator it=exportProcs.begin(); it!=exportProcs.end(); ++it) {
	 outgoingUpdates[counter] = hostUpdates[*it].size();
	 MPI_Isend(outgoingUpdates+counter,1,MPI_Type<size_t>(),*it,myrank,comm,&(sendRequests[counter]));
	 ++counter;
      }
      
      // Wait for number of host updates sends and recvs to complete:
      #ifdef PROFILE
         profile::start("MPI Waits",profMPI);
      #endif
      MPI_Waitall(importProcs.size(),&(recvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(exportProcs.size(),&(sendRequests[0]),MPI_STATUSES_IGNORE);
      #ifdef PROFILE
         profile::stop();
      #endif
      
      // Allocate buffers for sending and receiving host updates:
      CellID** incomingCellIDs      = new CellID* [importProcs.size()];
      CellID** outgoingCellIDs      = new CellID* [exportProcs.size()];
      MPI_processID** incomingHosts = new MPI_processID* [importProcs.size()];
      MPI_processID** outgoingHosts = new MPI_processID* [exportProcs.size()];
      for (size_t i=0; i<importProcs.size(); ++i) {
	 incomingCellIDs[i] = NULL;
	 incomingHosts[i] = NULL;
      }
      for (size_t i=0; i<exportProcs.size(); ++i) {
	 outgoingCellIDs[i] = NULL;
	 outgoingHosts[i] = NULL;
      }

      // Copy data to send buffers:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=exportProcs.begin(); it!=exportProcs.end(); ++it) {
	 outgoingCellIDs[counter] = new CellID[outgoingUpdates[counter]];
	 outgoingHosts[counter]   = new MPI_processID[outgoingUpdates[counter]];
	 
	 unsigned int i = 0;
	 for (std::set<std::pair<CellID,MPI_processID> >::const_iterator jt=hostUpdates[*it].begin(); jt!=hostUpdates[*it].end(); ++jt) {
	    outgoingCellIDs[counter][i] = jt->first;
	    outgoingHosts[counter][i]   = jt->second;
	    ++i;
	 }
	 ++counter;
      }
      
      // Allocate buffers for incoming cell updates:
      for (size_t i=0; i<importProcs.size(); ++i) {
	 incomingCellIDs[i] = new CellID[incomingUpdates[i]];
	 incomingHosts[i]   = new MPI_processID[incomingUpdates[i]];
      }
      
      // Receive host updates from neighbouring processes:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=importProcs.begin(); it!=importProcs.end(); ++it) {
	 MPI_Irecv(incomingCellIDs[counter],incomingUpdates[counter],MPI_Type<CellID>()       ,*it,*it,comm,&(recvRequests[2*counter+0]));
	 MPI_Irecv(incomingHosts[counter]  ,incomingUpdates[counter],MPI_Type<MPI_processID>(),*it,*it,comm,&(recvRequests[2*counter+1]));
	 ++counter;
      }
      
      // Send this process' update list to neighbouring processes:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=exportProcs.begin(); it!=exportProcs.end(); ++it) {
	 MPI_Isend(outgoingCellIDs[counter],outgoingUpdates[counter],MPI_Type<CellID>()       ,*it,myrank,comm,&(sendRequests[2*counter+0]));
	 MPI_Isend(outgoingHosts[counter]  ,outgoingUpdates[counter],MPI_Type<MPI_processID>(),*it,myrank,comm,&(sendRequests[2*counter+1]));
	 ++counter;
      }
      
      // Wait for host updates to complete:
      #ifdef PROFILE
         profile::start("MPI Waits",profMPI);
      #endif
      MPI_Waitall(2*importProcs.size(),&(recvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(2*exportProcs.size(),&(sendRequests[0]),MPI_STATUSES_IGNORE);
      #ifdef PROFILE
         profile::stop();
      #endif

      // Update hosts based on received information:
      counter = 0;
      for (std::set<MPI_processID>::const_iterator it=importProcs.begin(); it!=importProcs.end(); ++it) {
	 for (size_t i=0; i<incomingUpdates[counter]; ++i) {
	    std::map<CellID,CellID>::const_iterator it = newGlobal2LocalMap.find(incomingCellIDs[counter][i]);
	    if (it == newGlobal2LocalMap.end()) continue;
	    const CellID LID = it->second;
	    newHosts[LID] = incomingHosts[counter][i];
	 }
	 ++counter;
      }

      // Deallocate arrays used in host updates:
      delete [] incomingUpdates; incomingUpdates = NULL;
      delete [] outgoingUpdates; outgoingUpdates = NULL;
      for (size_t i=0; i<importProcs.size(); ++i) {
	 delete [] incomingCellIDs[i]; incomingCellIDs[i] = NULL;
	 delete [] incomingHosts[i]; incomingHosts[i] = NULL;
      }
      delete [] incomingHosts; incomingHosts = NULL;
      delete [] incomingCellIDs; incomingCellIDs = NULL;
      for (size_t i=0; i<exportProcs.size(); ++i) {
	 delete [] outgoingCellIDs[i]; outgoingCellIDs[i] = NULL;
	 delete [] outgoingHosts[i]; outgoingHosts[i] = NULL;
      }
      delete [] outgoingHosts; outgoingHosts = NULL;
      delete [] outgoingCellIDs; outgoingCellIDs = NULL;

      // Deallocate Zoltan arrays:
      zoltan->LB_Free_Part(&importGlobalIDs,&importLocalIDs,&importProcesses,&importParts);
      zoltan->LB_Free_Part(&exportGlobalIDs,&exportLocalIDs,&exportProcesses,&exportParts);

      // Swap cell data arrays:
      hosts.swap(newHosts);
      globalIDs.swap(newGlobalIDs);
      global2LocalMap.swap(newGlobal2LocalMap);
      N_localCells = N_newLocalCells;
      N_totalCells = hosts.size();
      cellNeighbours.swap(newCellNeighbours);
      neighbourFlags.swap(newNeighbourFlags);

      // Recalculate and/or invalidate all other internal data that depends on partitioning:
      #ifdef PROFILE
         profile::start("Stencil Recalculation",profStencilRecalc);
      #endif
      nbrProcesses.clear();
      for (CellID c=N_localCells; c<hosts.size(); ++c) {
	 nbrProcesses.insert(hosts[c]);
      }
      invalidate();
      #ifdef PROFILE
         profile::stop();
      #endif
      
      #ifndef NDEBUG
         checkInternalStructures();
      #endif

      ++partitioningCounter;
      
      #ifdef PROFILE
         profile::stop();
         profile::stop();
      #endif
      return success;
   }
   
   /** Synchronize MPI processes in the communicator ParGrid is using.*/
   template<class C> inline
   void ParGrid<C>::barrier() const {
      MPI_Barrier(comm);
   }

   /** Calculate neighbour (i,j,k) offset from the given neighbour type ID number.
    * @param nbrTypeID Neighbour type ID number.
    * @param i_off Variable in which calculated i-offset is written.
    * @param j_off Variable in which calculated j-offset is written.
    * @param k_off Variable in which calculated k-offset is written.*/
   template<class C> inline
   void ParGrid<C>::calcNeighbourOffsets(NeighbourID nbrTypeID,int& i_off,int& j_off,int& k_off) const {
      int tmp = nbrTypeID;
      k_off = tmp / 9;
      tmp -= k_off*9;
      j_off = tmp / 3;
      tmp -= j_off*3;
      i_off = tmp;
      --i_off;
      --j_off;
      --k_off;
   }
   
   /** Calculate neighbour type ID corresponding to the given cell offsets.
    * Valid offset values are [-1,0,+1].
    * @param i_off Cell offset in first coordinate direction.
    * @param j_off Cell offset in second coordinate direction.
    * @param k_off Cell offset in third coordinate direction.
    * @return Calculated neighbour type ID.*/
   template<class C> inline
   NeighbourID ParGrid<C>::calcNeighbourTypeID(int i_off,int j_off,int k_off) const {
      return (k_off+1)*9 + (j_off+1)*3 + i_off+1;
   }

   /** Debugging function, checks ParGrid internal structures for correctness.
    * @return If true, everything is ok.*/
   template<class C> inline
   bool ParGrid<C>::checkInternalStructures() const {
      if (getInitialized() == false) return false;
      bool success = true;
      std::map<CellID,int> tmpReferences;
      std::set<CellID> tmpRemoteCells;
      
      if (N_totalCells != hosts.size() || N_totalCells != globalIDs.size()) {
	 std::cerr << "N_totalCells has incorrect value " << N_totalCells << " should be ";
	 std::cerr << hosts.size() << '\t' << globalIDs.size() << std::endl;
	 success = false;
      }
      
      // Count neighbour references, and collect global IDs of remote neighbours,
      // Also check neighbourFlags for correctness:
      for (size_t cell=0; cell<N_localCells; ++cell) {
	 for (size_t i=0; i<N_neighbours; ++i) {
	    const CellID nbrLID = cellNeighbours[cell*N_neighbours+i];
	    if (nbrLID == invalid()) {
	       if (((neighbourFlags[cell] >> i) & 1) != 0) {
		  std::cerr << "P#" << myrank << " LID#" << cell << " GID#" << globalIDs[cell] << " nbrFlag is one for non-existing nbr type " << i << std::endl;
	       }
	       continue;
	    }
	    if (((neighbourFlags[cell] >> i ) & 1) != 1) {
	       std::cerr << "P#" << myrank << " LID#" << cell << " GID#" << globalIDs[cell] << " nbrFlag is zero for existing nbr type " << i;
	       std::cerr << " nbr LID#" << nbrLID << " GID#" << globalIDs[nbrLID];
	       std::cerr << std::endl;
	    }
		
	    if (nbrLID == cell) continue;
	    ++tmpReferences[nbrLID];
	    if (nbrLID >= N_localCells) tmpRemoteCells.insert(nbrLID);
	 }
      }

      // Assume that localCells is correct. Check that remoteCells is correct.
      for (std::set<CellID>::const_iterator it=tmpRemoteCells.begin(); it!=tmpRemoteCells.end(); ++it) {
	 const CellID localID = *it;
	 if (localID < N_localCells || localID >= N_totalCells) {
	    std::cerr << "P#" << myrank << " remote cell LID#" << *it << " has invalid LID!" << std::endl;
	    success = false;
	 }
      }

      // Check that remoteCells does not contain unnecessary entries:
      for (size_t i=N_localCells; i<N_totalCells; ++i) {
	 std::set<CellID>::const_iterator jt=tmpRemoteCells.find(i);
	 if (jt == tmpRemoteCells.end()) {
	    std::cerr << "P#" << myrank << " unnecessary remote cell entry in cells, LID#" << i << std::endl;
	    success = false;
	 }
      }
      
      // Check that localCells and remoteCells do not contain duplicate cells:
      std::set<CellID> tmpGlobalIDs;
      for (size_t i=0; i<N_totalCells; ++i) {
	 std::set<CellID>::const_iterator it = tmpGlobalIDs.find(globalIDs[i]);
	 if (it != tmpGlobalIDs.end()) {
	    std::cerr << "P#" << myrank << " LID#" << i << " GID#" << globalIDs[i] << " duplicate entry" << std::endl;
	    success = false;
	 }
	 tmpGlobalIDs.insert(globalIDs[i]);
      }
      
      // Check that all hosts are valid:
      for (size_t i=0; i<N_totalCells; ++i) {
	 if (i < N_localCells) {
	    // Check that all local cells have this process as their host:
	    if (hosts[i] != getRank()) {
	       std::cerr << "P#" << myrank << " LID#" << i << " GID#" << globalIDs[i] << " host " << hosts[i] << " should be " << getRank() << std::endl;
	       success = false;
	    }
	 } else {
	    // Remote cells must have sensible host value:
	    if (hosts[i] >= getProcesses()) {
	       std::cerr << "P#" << myrank << " LID#" << i << " GID#" << globalIDs[i] << " host " << hosts[i] << std::endl;
	       success = false;
	    }
	 }
	 if (i >= N_localCells) continue;
	 
	 // Check that all neighbour hosts have reasonable values:
	 for (size_t n=0; n<N_neighbours; ++n) {
	    if (cellNeighbours[i*N_neighbours+n] == invalid()) continue;
	 }
      }
      
      // Check that all remote cells have a correct host. Each process sends everyone else
      // a list of cells it owns. Hosts can be checked based on this information:
      for (MPI_processID p = 0; p < getProcesses(); ++p) {
	 CellID N_cells;
	 if (p == getRank()) {
	    // Tell everyone how many cells this process has:
	    N_cells = N_localCells;
	    MPI_Bcast(&N_cells,1,MPI_Type<CellID>(),p,comm);
	    
	    // Create an array containing cell IDs and broadcast:
	    CellID* buffer = const_cast<CellID*>(&(globalIDs[0]));
	    MPI_Bcast(buffer,N_cells,MPI_Type<CellID>(),p,comm);
	 } else {
	    // Receive number of cells process p is sending:
	    MPI_Bcast(&N_cells,1,MPI_Type<CellID>(),p,comm);
	    // Allocate array for receiving cell IDs and receive:
	    CellID* remoteCells = new CellID[N_cells];
	    MPI_Bcast(remoteCells,N_cells,MPI_Type<CellID>(),p,comm);
	    
	    for (CellID c=0; c<N_cells; ++c) {
	       // Check that received cells are not local to this process:
	       for (size_t i=0; i<N_localCells; ++i) {
		  if (globalIDs[i] == remoteCells[c]) {
		     std::cerr << "P#" << myrank << " GID#" << remoteCells[c] << " from P#" << p << " is local to this process!" << std::endl;
		     success = false;
		  }
	       }
	       // Check that remote cell has the correct host:
	       std::map<CellID,CellID>::const_iterator it = global2LocalMap.find(remoteCells[c]);
	       if (it != global2LocalMap.end()) {
		  const CellID localID = it->second;
		  if (localID < N_localCells) {
		     std::cerr << "P#" << myrank << " GID#" << remoteCells[c] << " from P#" << p << " is local to this process!" << std::endl;
		     success = false;
		  }
		  if (hosts[localID] != p) {
		     std::cerr << "P#" << myrank << " GID#" << remoteCells[c] << " from P#" << p << " has wrong host " << hosts[localID] << std::endl;
		     success = false;
		  }
	       }
	    }
	    delete [] remoteCells; remoteCells = NULL;
	 }
      }	    
      return success;
   }

   /** Check if partitioning has changed since last time checkPartitioningStatus 
    * was called. User must define a counter elsewhere and pass it as a parameter when 
    * calling this function. ParGrid then compares the value of that counter to internal 
    * value and if they differ, grid has been repartitioned.
    * @param counter User-defined partitioning counter. Initialize this to negative value.
    * @return If true, partitioning has changed since last time this function was called.*/
   template<class C> inline
   bool ParGrid<C>::checkPartitioningStatus(int& counter) const {
      #ifndef NDEBUG
         if (counter > partitioningCounter) {
	    std::cerr << "(PARGRID) ERROR: User-given counter value '" << counter << "' is invalid in checkPartitioningStatus!" << std::endl;
	 }
      #endif
      bool rvalue = false;
      if (counter < partitioningCounter) rvalue = true;
      counter = partitioningCounter;
      return rvalue;
   }

   /** Helper function for checking success or failure of simulation. Each process
    * must call this function simultaneously with their own success status as 
    * parameter. Value 'true' is returned if and only if all processes called 
    * this function with status 'true'.
    * @param status Success status of this process.
    * @return If true, all processes have success status true.*/
   template<class C> inline
   bool ParGrid<C>::checkSuccess(const bool& status) const {
      int32_t mySuccess = 0;
      int32_t globalSuccess = 0;
      if (status == false) mySuccess = 1;
      MPI_Reduce(&mySuccess,&globalSuccess,1,MPI_Type<int32_t>(),MPI_SUM,0,comm);
      MPI_Bcast(&globalSuccess,1,MPI_Type<int32_t>(),0,comm);
      if (globalSuccess == 0) return true;
      return false;
   }
   
   /** Set contents of cell weight array to default value.*/
   template<class C> inline
   void ParGrid<C>::clearCellWeights() {
      for (size_t i=0; i<cellWeights.size(); ++i) cellWeights[i] = DEFAULT_CELL_WEIGHT;
   }
   
   /** Finalize ParGrid. After this function returns ParGrid cannot be used 
    * without re-initialisation.
    * @return If true, ParGrid finalized successfully.*/
   template<class C> inline
   bool ParGrid<C>::finalize() {
      if (initialized == false) return false;
      initialized = false;
      stencils.clear();
      delete zoltan; zoltan = NULL;
      
      // Deallocate user data:
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second == NULL) continue;
	 it->second->finalize();
	 delete it->second; it->second = NULL;
      }
      
      // Deallocate dynamic user data:
      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::iterator
	   it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second == NULL) continue;
	 it->second->finalize();
	 delete it->second; it->second = NULL;
      }
      return true;
   }
   
   /** Get a list of boundary cells in the given Stencil. Boundary cells are cells that have 
    * one or more remote neighbours. Typically these cells cannot be propagated until 
    * neighbour data has been synchronized.
    * @param stencilID ID number of the Stencil.
    * @return Vector containing local IDs of exterior cells.*/
   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getBoundaryCells(StencilID stencilID) const {
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator it = stencils.find(stencilID);
      #ifndef NDEBUG
         if (it == stencils.end()) {
	    std::cerr << "(PARGRID) ERROR: Non-existing stencil " << stencilID << " requested in getBoundaryCells!" << std::endl;
	    exit(1);
	 }
      #endif
      return it->second.getBoundaryCells();
   }
   
   /** Get cell's neighbours. Non-existing neighbours have their global IDs 
    * set to value ParGrid::invalid(). This function may return a pointer to 
    * invalid memory area if an invalid local ID is passed.
    * @param localID Local ID of cell whose neighbours are requested.
    * @param Reference to vector containing neihbours' global IDs. Size 
    * of vector is given by pargrid::N_neighbours. Vector can be indexed with ParGrid::calcNeighbourTypeID.*/
   template<class C> inline
   CellID* ParGrid<C>::getCellNeighbourIDs(CellID localID) {
      #ifndef NDEBUG
         if (localID >= N_localCells) {
	    std::cerr << "(PARGRID) ERROR: getCellNeighbourIDs local ID#" << localID << " is too large!" << std::endl;
	    exit(1);
	 }
      #endif
      return &(cellNeighbours[localID*N_neighbours]);
   }
   
   /** Get array containing cell weights that are used in load balancing.
    * @return Vector containing cell weights, indexed with local IDs.*/
   template<class C> inline
   std::vector<CellWeight>& ParGrid<C>::getCellWeights() {return cellWeights;}

   /** Get MPI communicator that ParGrid uses.
    * @return Communicator used by ParGrid.*/
   template<class C> inline
   MPI_Comm ParGrid<C>::getComm() const {return comm;}

   /** Get information on data stored in dynamic arrays.
    * @param dataIDs Vector in which Data IDs of allocated arrays are copied.
    * @param names Names of allocated arrays are copied here.
    * @param sizeArrays Pointers to size arrays for each (dynamic) array are copied here.
    * @param datatypes String representations of data stored in each array are copied here.
    * @param byteSizes Byte sizes of data stored in each array are copied here.
    * @param pointers Pointers to actual data stored in each array are copied here.*/
   template<class C> inline
   void ParGrid<C>::getDynamicUserDataInfo(std::vector<DataID>& dataIDs,std::vector<std::string>& names,
					   std::vector<ArraySizetype*>& sizeArrays,std::vector<std::string>& datatypes,
					   std::vector<unsigned int>& byteSizes,std::vector<const char**>& pointers) const {
      dataIDs.clear();
      names.clear();
      sizeArrays.clear();
      datatypes.clear();
      byteSizes.clear();
      pointers.clear();
      
      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::const_iterator
	   it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second == NULL) continue;
	 dataIDs.push_back(it->first);
	 names.push_back(it->second->getName());
	 sizeArrays.push_back(it->second->getSizePointer());
	 datatypes.push_back(it->second->getDatatype());
	 byteSizes.push_back(it->second->getElementByteSize());
	 const char** ptr = const_cast<const char**>(it->second->getArrayPointer());
	 pointers.push_back(ptr);
      }
   }
   
   /** Get list of exterior cells. Exterior cells are cells that have at least one missing 
    * neighbour, i.e. these cells are situated at the boundary of the simulation domain.
    * @return Vector containing local IDs of exterior cells.*/
   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getExteriorCells() {
      if (recalculateExteriorCells == true) {
	 exteriorCells.clear();
	 for (size_t i=0; i<N_localCells; ++i) {
	    if (neighbourFlags[i] != ALL_NEIGHBOURS_EXIST) exteriorCells.push_back(i);
	 }
	 recalculateExteriorCells = false;
      }
      return exteriorCells;
   }

   /** Get global IDs of cells stored on this process. The returned vector 
    * contains both local and remote cell global IDs.
    * @return Global IDs of all cells on this process.*/
   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getGlobalIDs() const {return globalIDs;}
   
   /** Get host process IDs of all cells (local + remote) stored on this process.
    * @return Hosts for all cells that this process has a copy of.*/
   template<class C> inline 
   const std::vector<MPI_processID>& ParGrid<C>::getHosts() const {return hosts;}
   
   /** Query if ParGrid has initialized correctly.
    * The value returned by this function is set in ParGrid::initialize.
    * @return If true, ParGrid is ready for use.*/
   template<class C> inline
   bool ParGrid<C>::getInitialized() const {return initialized;}

   /** Get list of inner cells in the given Stencil. Inner cells are cells with zero 
    * remote neighbours. This function kills program if an invalid stencilID is used.
    * @param stencilID ID number of the Stencil.
    * @return Vector containing local IDs of inner cells.*/
   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getInnerCells(StencilID stencilID) const {
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator it = stencils.find(stencilID);
      if (it == stencils.end()) {
	 std::cerr << "(PARGRID) ERROR: Non-existing stencil " << stencilID << " requested in getInnerCells!" << std::endl;
	 exit(1);
      }
      return it->second.getInnerCells();
   }

   /** Get list of interior cells. Interior cells are cells whose all neighbours exist, i.e. 
    * they are not situated at the boundaries of the simulation domain.
    * @return Vector containing local IDs of interior cells.*/
   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getInteriorCells() {
      if (recalculateInteriorCells == true) {
	 const unsigned int ALL_EXIST = 134217728 - 1; // This value is 2^27 - 1, i.e. integer with first 27 bits flipped
	 interiorCells.clear();
	 for (size_t i=0; i<N_localCells; ++i) {
	    if (neighbourFlags[i] == ALL_EXIST) interiorCells.push_back(i);
	 }
	 recalculateInteriorCells = false;
      }
      return interiorCells;
   }

   /** Get the local ID of a cell with given global ID. If the specified cell was found on this 
    * process, either as a local cell or as copy of a remote cell, the returned local ID is valid.
    * @param globalID Global ID of the cell.
    * @return Local ID of the cell if the cell was found, or numeric_limits<pargrid::CellID>::max() otherwise.*/
   template<class C> inline
   CellID ParGrid<C>::getLocalID(CellID globalID) const {
      std::map<CellID,CellID>::const_iterator it = global2LocalMap.find(globalID);
      #ifndef NDEBUG
         if (it == global2LocalMap.end()) {
	    std::cerr << "(PARGRID) ERROR: In getLocalID, cell with global ID#" << globalID << " does not exist on P#" << myrank << std::endl;
	 }
      #endif
      if (it == global2LocalMap.end()) return invalidCellID();
      return it->second;
   }

   template<class C> inline
   const std::vector<CellID>& ParGrid<C>::getLocalIDs() const {return localIDs;}
   
   /** Get vector containing neighbour existence flags of local cells. The array 
    * is indexed with local IDs. The size of array is given by getNumberOfLocalCells().
    * @return Vector containing neighbour existence flags.*/
   template<class C> inline
   const std::vector<uint32_t>& ParGrid<C>::getNeighbourFlags() const {return neighbourFlags;}
   
   /** Get neighbour existence flags of given cell. This function may crash the simulation 
    * or return completely wrong information if an invalid local ID is used.
    * @param localID Local ID of the cell.
    * @return Neighbour flags of the given cell.*/
   template<class C> inline
   uint32_t ParGrid<C>::getNeighbourFlags(CellID localID) const {
      #ifndef NDEBUG
         if (localID >= N_localCells) {
	    std::cerr << "(PARGRID) ERROR: Local ID#" << localID << " too large in getNeighbourFlags!" << std::endl;
	 }
      #endif
      return neighbourFlags[localID];
   }
   
   /** Get a list of neighbour processes. A process is considered to be a neighbour 
    * if it has one or more this process' local cells' neighbours.
    * @return List of neighbour process IDs.*/
   template<class C> inline
   const std::set<MPI_processID>& ParGrid<C>::getNeighbourProcesses() const {return nbrProcesses;}
   
   /** Get the total number of cells on this process. This includes remote cells buffered on this process.
    * @return Total number of cells hosted and buffered on this process.*/
   template<class C> inline
   CellID ParGrid<C>::getNumberOfAllCells() const {return N_totalCells;}
   
   /** Get the number of cells on this process.
    * @return Number of local cells.*/
   template<class C> inline
   CellID ParGrid<C>::getNumberOfLocalCells() const {return N_localCells;}
   
   template<class C> inline
   std::size_t ParGrid<C>::getNumberOfReceives(StencilID stencilID,MPI_processID hostID) const {
      // Attempt to find the requested Stencil:
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator stencil = stencils.find(stencilID);
      if (stencil == stencils.end()) return 0;
      
      // Get received cells from Stencil:
      const std::map<MPI_processID,std::set<CellID> >& recvs = stencil->second.getRecvs();
      
      // Check that the Stencil receives cells from the given MPI process, and if it does,
      // return the number of received cells:
      typename std::map<MPI_processID,std::set<CellID> >::const_iterator host = recvs.find(hostID);
      if (host == recvs.end()) return 0;
      return host->second.size();
   }
   
   template<class C> inline
   std::size_t ParGrid<C>::getNumberOfSends(StencilID stencilID,MPI_processID hostID) const {
      // Attempt to find the requested Stencil:
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator stencil = stencils.find(stencilID);
      if (stencil == stencils.end()) return 0;
      
      // Get sent cells from Stencil:
      const std::map<MPI_processID,std::set<CellID> >& sends = stencil->second.getSends();
      
      // Check that the Stencil sends cells to the given MPI process, and if it does,
      // return the number of sent cells:
      typename std::map<MPI_processID,std::set<CellID> >::const_iterator host = sends.find(hostID);
      if (host == sends.end()) return 0;
      return host->second.size();
   }
   
   /** Get the number of MPI processes in the communicator used by ParGrid.
    * The value returned by this function is set in ParGrid::initialize.
    * @return Number of MPI processes in communicator comm.*/
   template<class C> inline
   MPI_processID ParGrid<C>::getProcesses() const {return N_processes;}
   
   /** Get the rank of this process in the MPI communicator used by ParGrid.
    * The value returned by this function is set in ParGrid::initialize.
    * @return MPI rank of this process in communicator comm.*/
   template<class C> inline
   MPI_processID ParGrid<C>::getRank() const {return myrank;}

   /** Get lists of cells received from remote processes when using the given Stencil.
    * @param stencilID ID of the Stencil.
    * @return List of cells received from each neighbour process.*/
   template<class C> inline
   const std::map<MPI_processID,std::set<CellID> >& ParGrid<C>::getReceives(StencilID stencilID) const {
      #ifndef NDEBUG
         if (getInitialized() == false) exit(1);
      #endif
      
      // Attemp to find the given Stencil:
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator stencil = stencils.find(stencilID);
      #ifndef NDEBUG
         if (stencil == stencils.end()) {
	    std::cerr << "(PARGRID) ERROR: Could not find Stencil with ID#" << stencilID << " in getReceives!" << std::endl;
	    exit(1);
	 }
      #endif
      return stencil->second.getRecvs();
   }
   
   /** Get the remote neighbours of a given local cell.
    * @param localID Local ID of the cell.
    * @param nbrTypeIDs Searched neighours.
    * @param nbrIDs Global IDs of searched remote neighbours are written here.
    * @param hosts MPI host processes of searched remote neighbours are written here.*/
   template<class C> inline
   bool ParGrid<C>::getRemoteNeighbours(CellID localID,const std::vector<NeighbourID>& nbrTypeIDs,std::vector<CellID>& nbrIDs) {
      nbrIDs.clear();
      hosts.clear();
      if (localID >= N_localCells) return false;
      
      // Iterate over given neighbour type IDs and check if this cell has 
      // those neighbours, and if those neighbours are remote:
      for (size_t n=0; n<nbrTypeIDs.size(); ++n) {
	 const CellID nbrLID = cellNeighbours[localID*N_neighbours + n];
	 if (cellNeighbours[nbrLID] == INVALID_CELLID) continue;
	 if (hosts[nbrLID] == getRank()) continue;
	 nbrIDs.push_back(nbrLID);
      }	 
      return true;
   }

   /** Get array containing remote updates received from neighbouring processes.
    * @param stencilID ID of the Stencil that was used to transfer updates.
    * @param userDataID ID of the user data array.
    * @param offsets Address of array containing offsets to array buffer is copied here.
    * @param buffer Address of array containing remote updates is copied here.
    * @return If true, offsets and buffer arrays contain valid pointers.*/
   template<class C> template<typename T> inline
   bool ParGrid<C>::getRemoteUpdates(StencilID stencilID,DataID userDataID,unsigned int*& offsets,T*& buffer) const {
      offsets = NULL;
      buffer = NULL;
      #ifndef NDEBUG
         if (initialized == false) return false;
         if (userDataStatic.find(userDataID) == userDataStatic.end()) return false;
      #endif
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator it=stencils.find(stencilID);
      if (it == stencils.end()) return false;
      
      char* ptr = NULL;
      const bool rvalue = it->second.getRemoteUpdates(userDataID,offsets,ptr);
      buffer = reinterpret_cast<T*>(ptr);
      return rvalue;
   }
   
   /** Get lists of cells sent to remote processes when using the given Stencil.
    * @param stencilID ID of the Stencil.
    * @return List of cells sent to each neighbour process.*/
   template<class C> inline
   const std::map<MPI_processID,std::set<CellID> >& ParGrid<C>::getSends(StencilID stencilID) const {
      #ifndef NDEBUG
         if (getInitialized() == false) exit(1);
      #endif
	
      // Attemp to find the given Stencil:
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::const_iterator stencil = stencils.find(stencilID);
      #ifndef NDEBUG
         if (stencil == stencils.end()) {
	    std::cerr << "(PARGRID) ERROR: Could not find Stencil with ID#" << stencilID << " in getSends!" << std::endl;
	    exit(1);
	 }
      #endif
      return stencil->second.getSends();
   }

   /** Get information on static user data arrays stored in ParGrid.
    * @param dataIDs Unique ParGrid data ID numbers.
    * @param names Unique array names.
    * @param elements Number of elements per cell in each array.
    * @param byteSizes Byte sizes of basic datatypes in each array.
    * @param pointers Pointer to each array.*/
   template<class C> inline
   void ParGrid<C>::getStaticUserDataInfo(std::vector<DataID>& dataIDs,std::vector<std::string>& names,std::vector<unsigned int>& elements,
					  std::vector<std::string>& datatypes,
					  std::vector<unsigned int>& byteSizes,std::vector<const char*>& pointers) const {
      dataIDs.clear();
      names.clear();
      elements.clear();
      datatypes.clear();					     
      byteSizes.clear();
      pointers.clear();
      
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second == NULL) continue;
	 dataIDs.push_back(it->first);
	 names.push_back(it->second->name);
	 elements.push_back(it->second->N_elements);
	 datatypes.push_back(it->second->datatype);
	 byteSizes.push_back(it->second->byteSize);
	 pointers.push_back(it->second->array);
      }
   }
   
   /** Get pointer to user-defined data array.
    * @param userDataID ID number of the array, as returned by addUserData.
    * @return Pointer to array or NULL if an array with given ID does not exist.
    * @see addUserData.*/
   template<class C> inline
   char* ParGrid<C>::getUserData(DataID userDataID) {
      if (userDataStatic.find(userDataID) == userDataStatic.end()) return NULL;
      return userDataStatic[userDataID]->array;
   }
   
   /** Get pointer to user-defined data array.
    * @param name Unique name of the array. This name was given in addUserData.
    * @return Pointer to array or NULL if an array with given name does not exist.
    * @see addUserData.*/
   template<class C> inline
   char* ParGrid<C>::getUserData(const std::string& name) {
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second->getName() == name) return it->second->array;
      }      
      return NULL;
   }
   
   /** Get number of elements in static user data array. Number of elements is the 
    * value of N_elements parameter given in addUserData function.
    * @param userDataID DataID of the array.
    * @return Number of elements in given array. Zero value is returned if the array does not exist.
    * @see addUserData.*/
   template<class C> inline
   std::size_t ParGrid<C>::getUserDataStaticElements(DataID userDataID) const {
      typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator staticIt = userDataStatic.find(userDataID);
      if (staticIt != userDataStatic.end()) return staticIt->second->N_elements;
      return 0;
   }
   
   /** Get DataID of given user data array.
    * @param name Name of the static or dynamic user data array.
    * @param isDynamic If user data array is dynamic, this variable will be set to value true.
    * @return DataID of the array, or pargrid::INVALID_DATAID, if a user data array
    * with given name does not exist.*/
   template<class C> inline
   DataID ParGrid<C>::getUserDataID(const std::string& name,bool& isDynamic) const {
      // Search static arrays for given name:
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second->getName() == name) {
	    isDynamic = false;
	    return it->first;
	 }
      }
      
      // Search dynamic arrays for given name:
      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::const_iterator 
	   it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second->getName() == name) {
	    isDynamic = true;
	    return it->first;
	 }
      }
      
      return INVALID_DATAID;
   }
   
   /** Get pointer to user-defined data array.
    * @param userDataID ID number of the array, as returned by addUserData.
    * @return Pointer to array or NULL if an array with given ID does not exist.
    * @see addUserData.*/
   template<class C> template<typename T> inline
   T* ParGrid<C>::getUserDataStatic(DataID userDataID) {
      typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::iterator it = userDataStatic.find(userDataID);
      if (it == userDataStatic.end()) return NULL;
      return reinterpret_cast<T*>(it->second->array);
   }
   
   /** Get pointer to user-defined data array.
    * @param name Unique name of the array. This name was given in addUserData.
    * @return Pointer to array or NULL if an array with given name does not exist.
    * @see addUserData.*/
   template<class C> template<typename T> inline
   T* ParGrid<C>::getUserDataStatic(const std::string& name) {
      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator
	   it=userDataStatic.begin(); it!=userDataStatic.end(); ++it) {
	 if (it->second->getName() == name) return reinterpret_cast<T*>(it->second->array);
      }
      return NULL;
   }
   
   /** Get user array containing dynamically allocated data.
    * @param userDataID Data ID of the array.
    * @return Wrapper to array.*/
   template<class C> template<typename T> inline
   DataWrapper<T> ParGrid<C>::getUserDataDynamic(DataID userDataID) {
      typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::iterator it = userDataDynamic.find(userDataID);
      
      // If dynamic data array with given ID does not exist,
      // return an invalid DataWrapper:
      if (it == userDataDynamic.end()) {
	 return DataWrapper<T>(NULL,NULL,0,NULL,0);
      }
      
      // Return DataWrapper associated with requested dynamic data array:
      return DataWrapper<T>(it->second->getArrayPointer(),
			    it->second->getCapacityPointer(),
			    it->second->getNumberOfCells(),
			    it->second->getSizePointer(),
			    it->second->getElementByteSize());
   }
   
   /** Get user array containing dynamically allocated data.
    * @param name Name of the array.
    * @return Wrapper to array.*/
   template<class C> template<typename T> inline
   DataWrapper<T> ParGrid<C>::getUserDataDynamic(const std::string& name) {
      // Iterate over all dynamic user data and search for the given name:
      typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::iterator it;
      for (it=userDataDynamic.begin(); it!=userDataDynamic.end(); ++it) {
	 if (it->second->getName() == name) {
	    return DataWrapper<T>(it->second->getArrayPointer(),
				  it->second->getCapacityPointer(),
				  it->second->getNumberOfCells(),
				  it->second->getSizePointer(),
				  it->second->getElementByteSize());
	 }
      }
      // Dynamic user data with given name was not found:
      return DataWrapper<T>(NULL,NULL,0,NULL,0);
   }

   /** Get the byte size of an array element in given user-defined array. 
    * The value is equal to number of elements per cell times byte size of basic datatype.
    * @param userDataID ID number of the user data array.
    * @return Zero if an array with given ID number does not exist, otherwise byte size of array element.*/
   template<class C> inline
   unsigned int ParGrid<C>::getUserDataElementSize(DataID userDataID) const {
      typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator 
	it=userDataStatic.find(userDataID);
      #ifndef NDEBUG
         if (it == userDataStatic.end()) return 0;
      #endif
      return it->second->getElementSize();
   }

   /** Get an MPI derived datatype that transfers all given cells of a user data array. 
    * MPI_Type_commit is called for the created datatype, MPI_Type_free must be called elsewhere.
    * @param userDataID ID number of the user array.
    * @param globalIDs Global IDs of transferred cells.
    * @param datatype MPI_Type_commit is called for this variable.
    * @param reverseStencil FIX ME
    * @return If true, user array with given ID number exists and an MPI datatype was created successfully.*/
   template<class C> inline
   bool ParGrid<C>::getUserDatatype(DataID userDataID,const std::set<CellID>& globalIDs,MPI_Datatype& datatype,bool reverseStencil) {
      typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::const_iterator it = userDataStatic.find(userDataID);
      #ifndef NDEBUG
         if (it == userDataStatic.end()) return false;
      #endif
      it->second->getDatatype(globalIDs,datatype);
      return true;
   }
   
   /** Initialize ParGrid and Zoltan. Note that MPI_Init must
    * have been called prior to calling this function.
    * @param comm MPI communicator that ParGrid should use.
    * @param parameters Load balancing parameters for all hierarchical levels.
    * The parameters for each hierarchical level are given in a map, whose contents are pairs
    * formed from parameter types and their string values. These maps themselves
    * are packed into a vector, whose first item (map) is used for hierarchical level
    * 0, second item for hierarchical level 1, and so forth. Zoltan is set to use 
    * hierarchical partitioning if vector size is greater than one, otherwise the 
    * load balancing method given in the first element is used.
    * @return If true, ParGrid initialized correctly.*/
   template<class C> inline
   bool ParGrid<C>::initialize(MPI_Comm comm,const std::vector<std::map<InputParameter,std::string> >& parameters) {
      zoltan = NULL;
      MPI_Comm_dup(comm,&(this->comm));
      
      // Get the number MPI of processes, and the rank of this MPI process, in given communicator:
      MPI_Comm_size(comm,&N_processes);
      MPI_Comm_rank(comm,&myrank);
      
      // Check that parameters vector is not empty:
      if (parameters.size() == 0) {
	 std::cerr << "PARGRID ERROR: parameters vector in constructor is empty!" << std::endl;
	 return false;
      }
      
      std::map<InputParameter,std::string> zoltanParameters;
      zoltanParameters[imbalanceTolerance]    = "IMBALANCE_TOL";
      zoltanParameters[loadBalancingMethod]   = "LB_METHOD";
      zoltanParameters[processesPerPartition] = "PROCS_PER_PART";
      
      // Parse user-defined load balancing parameters into a string,string container
      // which is more convenient to use with Zoltan:
      unsigned int level = 0;
      loadBalancingParameters.clear();
      loadBalancingParameters.resize(parameters.size());
      for (std::vector<std::map<InputParameter,std::string> >::const_iterator it=parameters.begin(); it!=parameters.end(); ++it) {
	 for (std::map<InputParameter,std::string>::const_iterator jt=it->begin(); jt!=it->end(); ++jt) {
	    std::map<InputParameter,std::string>::const_iterator kt=zoltanParameters.find(jt->first);
	    if (kt != zoltanParameters.end())
	      loadBalancingParameters[level].push_back(make_pair(kt->second,jt->second));
	 }
	 ++level;
      }

      // Determine if cell weights are calculated and passed to Zoltan:
      std::string objWeightDim;
      std::map<InputParameter,std::string>::const_iterator it = parameters.front().find(cellWeightScale);
      if (it != parameters.front().end()) {
	 cellWeightsUsed = true;
	 cellWeight      = atof(it->second.c_str());
	 objWeightDim    = "1";
      } else {
	 cellWeightsUsed = false;
	 cellWeight      = 0.0;
	 objWeightDim    = "0";
      }
      
      // Determine if edge weights are calculated and passed to Zoltan:
      std::string edgeWeightDim;
      it = parameters.front().find(edgeWeightScale);
      if (it != parameters.front().end()) {
	 edgeWeightsUsed = true;
	 edgeWeight      = atof(it->second.c_str());
	 edgeWeightDim   = "1";
      } else {
	 edgeWeightsUsed = false;
	 edgeWeight      = 0.0;
	 edgeWeightDim   = "0";
      }

      if (cellWeightsUsed == false && edgeWeightsUsed == false) {
	 std::cerr << "PARGRID ERROR: You must specify cell weight scale, edge weight scale, or both in ";
	 std::cerr << "the first element in parameters vector in constructor!" << std::endl;
	 return false;
      }

      // Create a new Zoltan object and set some initial parameters:
      zoltan = new Zoltan(comm);
      zoltan->Set_Param("NUM_GID_ENTRIES","1");
      zoltan->Set_Param("NUM_LID_ENTRIES","1");
      zoltan->Set_Param("RETURN_LISTS","ALL");
      zoltan->Set_Param("OBJ_WEIGHT_DIM",objWeightDim.c_str());
      zoltan->Set_Param("EDGE_WEIGHT_DIM",edgeWeightDim.c_str());
      zoltan->Set_Param("DEBUG_LEVEL","0");
      zoltan->Set_Param("REMAP","1");
      //zoltan->Set_Param("PHG_CUT_OBJECTIVE","CONNECTIVITY");
      //zoltan->Set_Param("CHECK_HYPERGRAPH","1");
      
      // Check if hierarchical partitioning should be enabled:
      if (parameters.size() > 1) {
	 
	 std::stringstream ss;
	 int counter = 0;
	 for (size_t i=0; i<loadBalancingParameters.size(); ++i) {
	    for (std::list<std::pair<std::string,std::string> >::const_iterator it=loadBalancingParameters[i].begin();
		 it != loadBalancingParameters[i].end(); ++it) {
	       if (it->first == "PROCS_PER_PART") {
		  if (counter > 0) ss << ',';
		  ss << it->second;
		  ++counter;
	       }
	    }
	 }
	 
	 zoltan->Set_Param("LB_METHOD","HIER");
	 zoltan->Set_Param("HIER_ASSIST","1");
	 //zoltan->Set_Param("HIER_CHECKS","0");
	 //zoltan->Set_Param("HIER_DEBUG_LEVEL","0");
	 zoltan->Set_Param("TOPOLOGY",ss.str());
	 
	 // These are set because of a Zoltan bug:
	 zoltan->Set_Param("EDGE_WEIGHT_DIM","0");
	 edgeWeightsUsed = false;
	 edgeWeightDim = "0";

	 if (getRank() == 0) std::cerr << "Enabling hierarchical partitioning" << std::endl;
      } else {
	 for (std::list<std::pair<std::string,std::string> >::const_iterator it=loadBalancingParameters[0].begin();
	      it != loadBalancingParameters[0].end(); ++it) {
	    if (it->first == "PROCS_PER_PART") continue;
	    zoltan->Set_Param(it->first,it->second);
	 }
      }

      // Register Zoltan callback functions:
      zoltan->Set_Num_Obj_Fn(&cb_getNumberOfLocalCells<ParGrid<C> >,this);               // All methods
      zoltan->Set_Obj_List_Fn(&cb_getLocalCellList<ParGrid<C> >,this);
      zoltan->Set_Num_Geom_Fn(&cb_getMeshDimension<ParGrid<C> >,this);                   // Geometrical
      zoltan->Set_Geom_Multi_Fn(&cb_getAllCellCoordinates<ParGrid<C> >,this);
      zoltan->Set_Num_Edges_Multi_Fn(&cb_getNumberOfAllEdges<ParGrid<C> >,this);         // Graph
      zoltan->Set_Edge_List_Multi_Fn(&cb_getAllCellEdges<ParGrid<C> >,this);
      zoltan->Set_HG_Size_CS_Fn(&cb_getNumberOfHyperedges<ParGrid<C> >,this);            // Hypergraph
      zoltan->Set_HG_CS_Fn(&cb_getHyperedges<ParGrid<C> >,this);
      zoltan->Set_HG_Size_Edge_Wts_Fn(cb_getNumberOfHyperedgeWeights<ParGrid<C> >,this);
      zoltan->Set_HG_Edge_Wts_Fn(cb_getHyperedgeWeights<ParGrid<C> >,this);
      //zoltan->Set_Hier_Num_Levels_Fn(cb_getNumberOfHierarchicalLevels<C>,this);        // Hierarchical
      //zoltan->Set_Hier_Part_Fn(cb_getHierarchicalPartNumber<C>,this);
      //zoltan->Set_Hier_Method_Fn(&cb_getHierarchicalParameters<C>,this);

      initialized = true;
      return initialized;
   }
   
   /** Return invalid cell ID. Cell IDs obtained elsewhere may be 
    * tested against this value to see if they are valid. DEPRECATED.
    * @return Invalid cell global ID.*/
   template<class C> inline
   CellID ParGrid<C>::invalid() const {return std::numeric_limits<CellID>::max();}
   
   /** Return an invalid cell ID.
    * @return Invalid cell (local or global) ID.*/
   template<class C> inline
   CellID ParGrid<C>::invalidCellID() const {return std::numeric_limits<CellID>::max();}
   
   /** Return an invalid user data ID.
    * @return Invalid user data ID.*/
   template<class C> inline
   DataID ParGrid<C>::invalidDataID() const {return std::numeric_limits<DataID>::max();}
   
   /** Return an invalid stencil ID.
    * @return Invalid stencil ID.*/
   template<class C> inline
   StencilID ParGrid<C>::invalidStencilID() const {return std::numeric_limits<StencilID>::max();}

   /** Invalidate ParGrid internal variables. This function gets called after mesh has been repartitioned.*/
   template<class C> inline
   void ParGrid<C>::invalidate() {
      recalculateInteriorCells = true;
      recalculateExteriorCells = true;
      for (typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator it=stencils.begin(); it!=stencils.end(); ++it) {
	 it->second.update();
      }
      
      // Reset size of vector cellWeights to N_localCells:
	{
	   std::vector<CellWeight> newCellWeights;
	   cellWeights.swap(newCellWeights);
	}
      if (cellWeightsUsed == true) {
	 cellWeights.resize(N_localCells);
	 for (size_t i=0; i<cellWeights.size(); ++i) cellWeights[i] = DEFAULT_CELL_WEIGHT;
      }
      
      // Sync local IDs:
      if (syncLocalIDs() == false) {
	 std::cerr << "(PARGRID) ERROR: Failed to sync localIDs vector!" << std::endl;
	 exit(1);
      }
   }
   
   /** Check if a cell with given global ID exists on this process.
    * @param localID Local ID of the searched cell.
    * @return If true, the cell exists on this process.*/
   template<class C> inline
   bool ParGrid<C>::localCellExists(CellID localID) {
      if (localID >= N_localCells) return false;
      return true;
   }

   /** Remove a data transfer that has been previously associated with the given Stencil.
    * @param userDataID ID of the user data array.
    * @param stencilID ID of the Stencil.
    * @return If true, the transfer was removed successfully.*/
   template<class C> inline
   bool ParGrid<C>::removeDataTransfer(DataID userDataID,StencilID stencilID) {
      if (getInitialized() == false) return false;
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator it = stencils.find(stencilID);
      if (it == stencils.end()) return false;
      return it->second.removeTransfer(userDataID);
   }
   
   /** Remove a transfer stencil from ParGrid.
    * @param stencilID ID of the Stencil.
    * @return If true, the Stencil was removed successfully.*/
   template<class C> inline
   bool ParGrid<C>::removeStencil(StencilID stencilID) {
      if (stencilID == pargrid::DEFAULT_STENCIL) return false;
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator it = stencils.find(stencilID);
      if (it == stencils.end()) return false;
      stencils.erase(it);
      return true;
   }
   
   /** Remove a user data array that has been previously allocated with addUserData.
    * @param userDataID ID number of the user data array, as returned by addUserData.
    * @return If true, the user data array was removed successfully.
    * @see addUserData.*/
   template<class C> inline
   bool ParGrid<C>::removeUserData(DataID userDataID) {
      typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::iterator user=userDataStatic.find(userDataID);
      if (user != userDataStatic.end()) {
	 // Remove transfers associated with this user data array. Note that 
	 // we do not know transferID(s) associated with this userDataID:
	 for (typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator it=stencils.begin(); it!=stencils.end(); ++it) {
	    it->second.removeTransfer(userDataID);
	 }
	 
	 userDataHoles.insert(userDataID);
	 delete user->second; user->second = NULL;
	 userDataStatic.erase(userDataID);
	 return true;
      }
      
      typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::iterator userDyn=userDataDynamic.find(userDataID);
      if (userDyn != userDataDynamic.end()) {
	 // Remove transfers associated with this dynamic user data array:
	 // TODO
	 
	 // Delete dynamic user data array:
	 userDataHoles.insert(userDataID);
	 delete userDyn->second; userDyn->second = NULL;
	 userDataDynamic.erase(userDataID);
	 return true;
      }
      return false;
   }
   
   /** Repartition user-defined data arrays. This function is called by balanceLoad.
    * @param N_cells Number of needed parallel cells in new mesh partitioning.
    * @param newLocalsBegin Index into arrays from which new (imported) cell data begin.
    * @param N_import Number of imported cells.
    * @param importProcesses MPI rank of process that imports cell i to this process.
    * @param importsPerProcess Number of cells imported from each process.
    * @param N_export Number of exported cells.
    * @param exportProcesses MPI rank of process that receives cell i from this process.
    * @param exportLocalIDs Local ID of exported cell i.
    * @param exportGlobalIDs Global ID of exported cell i.
    * @return If true, user data was repartitioned successfully.*/
   template<class C> inline
   bool ParGrid<C>::repartitionUserDataStatic(size_t N_cells,CellID newLocalsBegin,int N_import,int* importProcesses,
					      const std::map<MPI_processID,std::vector<int> >& importDisplacements,
					      int N_export,int* exportProcesses,ZOLTAN_ID_PTR exportLocalIDs,
					      const std::map<MPI_processID,std::vector<int> >& exportDisplacements) {
      // Skip if there are no static user data arrays:
      if (userDataStatic.size() == 0) return true;

      for (typename std::map<DataID,UserDataStatic<ParGrid<C> >*>::iterator
	   userData=userDataStatic.begin(); userData!=userDataStatic.end(); ++userData) {
	 // Allocate a temporary container for the new user data array:
	 UserDataStatic<ParGrid<C> > newUserData;
	 newUserData.initialize(this,userData->second->name,N_cells,userData->second->N_elements,
				userData->second->byteSize,userData->second->datatype);

	 // Create an MPI datatype that transfers a single user data array element:
	 const int byteSize = userData->second->N_elements * userData->second->byteSize;

	 // If static user data array contains no data, just swap arrays and skip MPI:
	 if (byteSize == 0) {
	    userData->second->swap(newUserData);
	    newUserData.finalize();
	    continue;
	 }

	 MPI_Datatype basicDatatype;
	 MPI_Type_contiguous(byteSize,MPI_Type<char>(),&basicDatatype);
	 MPI_Type_commit(&basicDatatype);

	 std::vector<MPI_Request> recvRequests(importDisplacements.size());
	 std::vector<MPI_Request> sendRequests(exportDisplacements.size());
	 std::vector<MPI_Datatype> recvDatatypes(importDisplacements.size());
	 std::vector<MPI_Datatype> sendDatatypes(exportDisplacements.size());

	 // Iterate over all processes sending cells to this process:
	 size_t recvRequestCounter = 0;
	 for (typename std::map<MPI_processID,std::vector<int> >::const_iterator it=importDisplacements.begin(); it!=importDisplacements.end(); ++it) {
	    // Create an MPI datatype that transfers all incoming cells:
	    int* displacements = const_cast<int*>(it->second.data());
	    MPI_Type_create_indexed_block(it->second.size(),1,displacements,basicDatatype,&(recvDatatypes[recvRequestCounter]));
	    MPI_Type_commit(&(recvDatatypes[recvRequestCounter]));

	    // Post a receive:
	    const MPI_processID source = it->first;
	    const int tag              = it->first;
	    char* buffer               = newUserData.array;
	    MPI_Irecv(buffer,1,recvDatatypes[recvRequestCounter],source,tag,comm,&(recvRequests[recvRequestCounter]));
	    ++recvRequestCounter;
	 }
	 
	 // Iterate over all processes whom to send cells:
	 size_t sendRequestCounter = 0;
	 for (typename std::map<MPI_processID,std::vector<int> >::const_iterator it=exportDisplacements.begin(); it!=exportDisplacements.end(); ++it) {
	    // Create an MPI datatype that transfers all outgoing cells:
	    int* displacements = const_cast<int*>(it->second.data());
	    MPI_Type_create_indexed_block(it->second.size(),1,displacements,basicDatatype,&(sendDatatypes[sendRequestCounter]));
	    MPI_Type_commit(&(sendDatatypes[sendRequestCounter]));
	    
	    // Post a send:
	    const MPI_processID dest = it->first;
	    const int tag            = getRank();
	    char* buffer             = userData->second->array;
	    MPI_Isend(buffer,1,sendDatatypes[sendRequestCounter],dest,tag,comm,&(sendRequests[sendRequestCounter]));
	    ++sendRequestCounter;
	 }

	 // Copy the cells that remain on this process:
	 CellID counter = 0;
	 for (CellID cell=0; cell<N_localCells; ++cell) {
	    if (hosts[cell] != getRank()) continue;
	    newUserData.copy(*(userData->second),counter,cell);
	    ++counter;
	 }

	 // Wait for MPI transfers:
	 MPI_Waitall(recvRequestCounter,&(recvRequests[0]),MPI_STATUSES_IGNORE);
	 MPI_Waitall(sendRequestCounter,&(sendRequests[0]),MPI_STATUSES_IGNORE);

	 // Swap user data containers and deallocate memory:
	 userData->second->swap(newUserData);
	 newUserData.finalize();

	 for (size_t i=0; i<recvRequestCounter; ++i) MPI_Type_free(&(recvDatatypes[i]));
	 for (size_t i=0; i<sendRequestCounter; ++i) MPI_Type_free(&(sendDatatypes[i]));
	 MPI_Type_free(&basicDatatype);
      }

      return true;
   }

   /**
    * @param N_cells Total number of cells (local + remote) in new mesh partition.
    * @param newLocalsBegin Start index to data arrays where imported cells are copied.
    * @param N_import Total number of cells imported to this process. 
    * @param importProcesses For each imported cell, the MPI process ID who sends the data to this process.
    * Note that importProcess may equal the MPI process ID of this process.
    * @param importDisplacements For each MPI process sending cells to this process, a list of displacements 
    * (starting from newLocalsBegin) where the cell data is to be copied.
    * @param N_export Total number of cells this process exports to other processes.
    * @param exportProcesses For each exported cell, the MPI process ID of the process receiving the cell.
    * @param exportLocalIDs The local index of each exported cell.
    * @param exportDisplacements For each MPI process receiving cells from this process, a list of cell
    * local indices of cells exported to that proces.
    * @return If true, dynamic user data was successfully exported.*/
   template<class C> inline
   bool ParGrid<C>::repartitionUserDataDynamic(size_t N_cells,CellID newLocalsBegin,int N_import,int* importProcesses,
					       const std::map<MPI_processID,std::vector<int> >& importDisplacements,
					       int N_export,int* exportProcesses,ZOLTAN_ID_PTR exportLocalIDs,
					       const std::map<MPI_processID,std::vector<int> >& exportDisplacements) {
      // Return if there are no dynamic user data arrays:
      if (userDataDynamic.size() == 0) return true;

      bool success = true;

      for (typename std::map<DataID,UserDataDynamic<ParGrid<C> >*>::iterator
	   userData=userDataDynamic.begin(); userData!=userDataDynamic.end(); ++userData) {
	 // Create a new dynamic user data container:
	 UserDataDynamic<ParGrid<C> > newUserData;
	 newUserData.initialize(this,userData->second->getName(),N_cells,userData->second->getElementByteSize(),userData->second->getDatatype());
	 
	 if (userData->second->getElementByteSize() == 0) {
	    userData->second->swap(newUserData);
	    newUserData.finalize();
	    continue;
	 }
	 
	 std::vector<MPI_Request> recvRequests(importDisplacements.size());
	 std::vector<MPI_Request> sendRequests(exportDisplacements.size());

	 std::vector<MPI_Datatype> recvDatatypes(importDisplacements.size());
	 std::vector<MPI_Datatype> sendDatatypes(exportDisplacements.size());

	 // Iterate over all importing processes and receive the size array:
	 size_t recvRequestCounter = 0;
	 for (std::map<MPI_processID,std::vector<int> >::const_iterator it=importDisplacements.begin(); it!=importDisplacements.end(); ++it) {
	    int* displacements = const_cast<int*>(it->second.data());
	    MPI_Type_create_indexed_block(it->second.size(),1,displacements,MPI_Type<ArraySizetype>(),&(recvDatatypes[recvRequestCounter]));
	    MPI_Type_commit(&(recvDatatypes[recvRequestCounter]));
	    
	    const MPI_processID source = it->first;
	    const int tag              = it->first;
	    void* buffer               = newUserData.getSizePointer();
	    MPI_Irecv(buffer,1,recvDatatypes[recvRequestCounter],source,tag,comm,&(recvRequests[recvRequestCounter]));
	    ++recvRequestCounter;
	 }

	 size_t sendRequestCounter = 0;
	 for (std::map<MPI_processID,std::vector<int> >::const_iterator it=exportDisplacements.begin(); it!=exportDisplacements.end(); ++it) {
	    int* displacements = const_cast<int*>(it->second.data());
	    MPI_Type_create_indexed_block(it->second.size(),1,displacements,MPI_Type<ArraySizetype>(),&(sendDatatypes[sendRequestCounter]));
	    MPI_Type_commit(&(sendDatatypes[sendRequestCounter]));

	    const MPI_processID source = it->first;
	    const int tag              = getRank();
	    void* buffer               = userData->second->getSizePointer();
	    MPI_Isend(buffer,1,sendDatatypes[sendRequestCounter],source,tag,comm,&(sendRequests[sendRequestCounter]));
	    ++sendRequestCounter;
	 }

	 // Copy sizes entries of cells remaining on this process:
	 const ArraySizetype* oldSizes = userData->second->getSizePointer();
	 ArraySizetype* newSizes       = newUserData.getSizePointer();
	 size_t counter = 0;
	 for (CellID cell=0; cell<N_localCells; ++cell) {
	    if (hosts[cell] != getRank()) continue;
	    newSizes[counter] = oldSizes[cell];
	    ++counter;
	 }

	 // Wait for MPI transfers of size-arrays to complete:
	 MPI_Waitall(recvRequestCounter,&(recvRequests[0]),MPI_STATUSES_IGNORE);
	 MPI_Waitall(sendRequestCounter,&(sendRequests[0]),MPI_STATUSES_IGNORE);

	 for (size_t i=0; i<recvRequestCounter; ++i) MPI_Type_free(&(recvDatatypes[i]));
	 for (size_t i=0; i<sendRequestCounter; ++i) MPI_Type_free(&(sendDatatypes[i]));

	 // Allocate memory for incoming elements in new user data container:
	 newUserData.reallocate(0,N_cells);

	 // Create an MPI datatype that transfers a single element:
	 MPI_Datatype basicDatatype;
	 MPI_Type_contiguous(userData->second->getElementByteSize(),MPI_Type<char>(),&basicDatatype);
	 MPI_Type_commit(&basicDatatype);

	 recvRequestCounter = 0;
	 for (std::map<MPI_processID,std::vector<int> >::const_iterator it=importDisplacements.begin(); it!=importDisplacements.end(); ++it) {
	    // Iterate over all cells received from process it->first, 
	    // and store the buffer addresses to vector offsets:
	    std::vector<MPI_Aint> offsets;
	    std::vector<int> blockLengths;
	    for (size_t i=0; i<it->second.size(); ++i) {
	       // Get the local ID of the received cell:
	       const CellID cellLID = it->second[i];
	       
	       // Get pointer to start of the data:
	       void* ptr = newUserData.getArrayPointer()[cellLID];

	       // Calculate offset relative to MPI_BOTTOM:
	       MPI_Aint address;
	       if (ptr == NULL) address = 0;
	       else MPI_Get_address(ptr,&address);

	       offsets.push_back(address);
	       blockLengths.push_back(newUserData.getSizePointer()[cellLID]);
	    }

	    // Post a receive for cells imported from process it->first:
	    const int amount           = it->second.size();
	    const MPI_processID source = it->first;
	    const int tag              = it->first;

	    MPI_Type_hindexed(amount,&(blockLengths[0]),&(offsets[0]),basicDatatype,&(recvDatatypes[recvRequestCounter]));
	    MPI_Type_commit(&(recvDatatypes[recvRequestCounter]));
	    MPI_Irecv(MPI_BOTTOM,1,recvDatatypes[recvRequestCounter],source,tag,comm,&(recvRequests[recvRequestCounter]));
	    ++recvRequestCounter;
	 }
	 
	 sendRequestCounter = 0;
	 for (std::map<MPI_processID,std::vector<int> >::const_iterator it=exportDisplacements.begin(); it!=exportDisplacements.end(); ++it) {
	    // Iterate over all cells exported to process it->first,
	    // and store the buffer addresses to vector offsets:
	    std::vector<MPI_Aint> offsets;
	    std::vector<int> blockLengths;
	    for (size_t i=0; i<it->second.size(); ++i) {
	       // Get the local ID of the exported cell:
	       const CellID cellLID = it->second[i];
	       
	       // Get pointer to start of the data:
	       void* ptr = userData->second->getArrayPointer()[cellLID];

	       // Calculate offset relative to MPI_BOTTOM:
	       MPI_Aint address;
	       if (ptr == NULL) address = 0;
	       else MPI_Get_address(ptr,&address);

	       offsets.push_back(address);
	       blockLengths.push_back(userData->second->getSizePointer()[cellLID]);
	    }
	    
	    // Post a send for cells exported to process it->first:
	    const int amount         = it->second.size();
	    const MPI_processID dest = it->first;	    
	    const int tag            = getRank();

	    MPI_Type_hindexed(amount,&(blockLengths[0]),&(offsets[0]),basicDatatype,&(sendDatatypes[sendRequestCounter]));
	    MPI_Type_commit(&(sendDatatypes[sendRequestCounter]));
	    MPI_Isend(MPI_BOTTOM,1,sendDatatypes[sendRequestCounter],dest,tag,comm,&(sendRequests[sendRequestCounter]));
	    ++sendRequestCounter;
	 }

	 // Copy cells remaining on this process:
	 const ArraySizetype elementByteSize = userData->second->getElementByteSize();
	 char** oldArray                     = userData->second->getArrayPointer();
	 char** newArray                     = newUserData.getArrayPointer();
	 counter = 0;
	 for (CellID cell=0; cell<N_localCells; ++cell) {
	    if (hosts[cell] != getRank()) continue;
	    for (size_t i=0; i<userData->second->getSizePointer()[cell]*elementByteSize; ++i)
	      newArray[counter][i] = oldArray[cell][i];
	    ++counter;
	 }

	 // Wait for MPI transfers to complete:
	 MPI_Waitall(recvRequestCounter,&(recvRequests[0]),MPI_STATUSES_IGNORE);
	 MPI_Waitall(sendRequestCounter,&(sendRequests[0]),MPI_STATUSES_IGNORE);

	 // Swap arrays:
	 userData->second->swap(newUserData);
	 newUserData.finalize();

	 for (size_t i=0; i<recvRequestCounter; ++i) MPI_Type_free(&(recvDatatypes[i]));
	 for (size_t i=0; i<sendRequestCounter; ++i) MPI_Type_free(&(sendDatatypes[i]));
	 MPI_Type_free(&basicDatatype);
      }
      return success;
   }

   /** Set partitioning mode. This only has effect if GRAPH or HYPERGRAPH 
    * partitioners is used. Valid values are partition (static load balancing),
    * repartition (dynamic load balancing), and refine (fast improvement).
    * @param pm New partitioning mode. These values are passed to Zoltan. Note 
    * that one needs to call ParGrid::balanceLoad before the new mode has any effect.
    * @return If true, new mode was successfully taken into use.*/
   template<class C> inline
   bool ParGrid<C>::setPartitioningMode(PartitioningMode pm) {
      bool rvalue = true;
      switch (pm) {
       case partition:
	 zoltan->Set_Param("LB_APPROACH","PARTITION");
	 break;
       case repartition:
	 zoltan->Set_Param("LB_APPROACH","REPARTITION");
	 break;
       case refine:
	 zoltan->Set_Param("LB_APPROACH","REFINE");
	 break;
       default:
	 rvalue = false;
	 break;
      }
      return rvalue;
   }
   
   template<class C> inline
   void ParGrid<C>::startMetadataRepartitioning(std::vector<MPI_Request>& nbrRecvRequests,
						std::map<MPI_processID,std::vector<int> >& importDisplacements,
						std::vector<MPI_Request>& nbrSendRequests,
						std::map<MPI_processID,std::vector<int> >& exportDisplacements,
						std::vector<CellID>& newCellNeighbours) {
      // Create an MPI datatype for transferring all neighbour IDs of a single cell:
      MPI_Datatype basicDatatype;
      MPI_Datatype datatype;
      MPI_Type_contiguous(N_neighbours,MPI_Type<CellID>(),&basicDatatype);
      MPI_Type_commit(&basicDatatype);

      // Allocate enough MPI requests:
      nbrRecvRequests.resize(importDisplacements.size());
      nbrSendRequests.resize(exportDisplacements.size());

      // Receive all neighbour IDs:
      size_t counter = 0;
      for (std::map<MPI_processID,std::vector<int> >::iterator
	   it=importDisplacements.begin(); it!=importDisplacements.end(); ++it) {
	 // Create a datatype for receiving all neighbour IDs with a single receive:
	 MPI_Type_create_indexed_block(it->second.size(),1,&(it->second[0]),basicDatatype,&datatype);
	 MPI_Type_commit(&datatype);

	 // Post a receive:
	 const MPI_processID source = it->first;
	 const int tag              = it->first;
	 void* buffer               = &(newCellNeighbours[0]);

	 MPI_Irecv(buffer,1,datatype,source,tag,comm,&(nbrRecvRequests[counter]));
	 MPI_Type_free(&datatype);
	 ++counter;
      }

      // Send cell neighbour IDs:
      counter = 0;
      for (std::map<MPI_processID,std::vector<int> >::iterator
	   it=exportDisplacements.begin(); it!=exportDisplacements.end(); ++it) {
	 // Create a datatype for sending all neighbour IDs with a single send:
	 MPI_Type_create_indexed_block(it->second.size(),1,&(it->second[0]),basicDatatype,&datatype);
	 MPI_Type_commit(&datatype);

	 // Post a send:
	 const MPI_processID dest = it->first;
	 const int tag            = getRank();
	 void* buffer             = &(cellNeighbours[0]);

	 MPI_Isend(buffer,1,datatype,dest,tag,comm,&(nbrSendRequests[counter]));
	 MPI_Type_free(&datatype);
	 ++counter;
      }

      MPI_Type_free(&basicDatatype);
   }
   
   /** Start remote cell data synchronization. 
    * @param stencilID ID of the stencil.
    * @param transferID ID of the data transfer.
    * @return If true, synchronization started successfully.*/
   template<class C> inline
   bool ParGrid<C>::startNeighbourExchange(StencilID stencilID,DataID userDataID) {
      if (getInitialized() == false) return false;
      if (stencils.find(stencilID) == stencils.end()) return false;
      return stencils[stencilID].startTransfer(userDataID);
   }

   /** Every process tells every other process the global IDs of the cells it hosts. 
    * Based on this information processes know who hosts the remote neighbours of 
    * their local cells. This function should only be called once at the start of simulation.
    * @return If true, process information was exchanged successfully.*/
   template<class C> inline
   bool ParGrid<C>::syncCellHosts() {
      if (getInitialized() == false) return false;
      
      // Each process sends every other process a list of cells it owns. 
      // We can then figure out remote neighbour hosts from these lists:
      CellID N_cells;
      for (MPI_processID p = 0; p < getProcesses(); ++p) {
	 if (p == getRank()) {
	    // It is this processes turn to send a list of local cells. 
	    // First tell how many cells this process has:
	    N_cells = N_localCells;
	    MPI_Bcast(&N_cells,1,MPI_Type<CellID>(),p,comm);
	    MPI_Bcast(&(globalIDs[0]),N_cells,MPI_Type<CellID>(),p,comm);
	 } else {
	    // Receive a list of global cell IDs from another process:
	    MPI_Bcast(&N_cells,1,MPI_Type<CellID>(),p,comm);

	    CellID* remoteCells = new CellID[N_cells];
	    MPI_Bcast(remoteCells,N_cells,MPI_Type<CellID>(),p,comm);

	    // Go through the received list and check if any of 
	    // this processes remote neighbours are on that list:
	    for (CellID c=0; c<N_cells; ++c) {
	       std::map<CellID,CellID>::const_iterator it = global2LocalMap.find(remoteCells[c]);
	       if (it == global2LocalMap.end()) continue;
	       const CellID localID = it->second;
	       #ifndef NDEBUG
	          // DEBUG: Check that obtained cell is not local to this process:
	          if (localID < N_localCells) {
		     std::cerr << "(PARGRID) ERROR: P#" << getRank() << " remote cell GID#" << remoteCells[c] << " is local to this process!" << std::endl;
		     exit(1);
	          }
	          // DEBUG: Check that obtained local ID does not exceed array size:
	          if (localID >= N_totalCells) {
		     std::cerr << "(PARGRID) ERROR: P#" << getRank() << " remote cell GID#" << remoteCells[c] << " local ID# ";
		     std::cerr << localID << " above array bounds!" << std::endl;
		     exit(1);
		  }
	          // DEBUG: Check that globalIDs entry has correct value:
	          if (globalIDs[localID] != remoteCells[c]) {
		     std::cerr << "(PARGRID) ERROR: P#" << getRank() << " remote cell GID#" << remoteCells[c] << " has LID#" << localID;
		     std::cerr << " but globalIDs[localID] has value " << globalIDs[localID] << std::endl;
		  }
	       #endif
	       hosts[localID] = p;
	    }
	    delete [] remoteCells; remoteCells = NULL;
	 }
      }
      return true;
   }
   
   /** Recalculate the contents of vector localIDs and sync data with neighbour processes.
    * This function is called after mesh has been (re)partitioned.
    * @return If true, localIDs vector was re-constructed and synced successfully.*/
   template<class C> inline
   bool ParGrid<C>::syncLocalIDs() {
      #ifndef NDEBUG
         if (getInitialized() == false) return false;
      #endif
      
      // Reallocate vector localIDs:
	{
	   std::vector<CellID> dummy;
	   localIDs.swap(dummy);
	}
      localIDs.resize(N_totalCells);
      for (CellID c=0; c<N_localCells; ++c) localIDs[c] = c;
      
      // Get lists of cells to send and receive, ordered by global IDs:
      const std::map<MPI_processID,std::set<CellID> >& recvs = stencils[DEFAULT_STENCIL].getRecvs();
      const std::map<MPI_processID,std::set<CellID> >& sends = stencils[DEFAULT_STENCIL].getSends();
      
      // Make sure enough MPI Requests are available:
      if (recvRequests.size() < recvs.size()) recvRequests.resize(recvs.size());
      if (sendRequests.size() < sends.size()) sendRequests.resize(sends.size());
      
      // Calculate receive displacements and post receive(s):
      int recvCounter = 0;
      std::vector<int> displacements;
      for (typename std::map<MPI_processID,std::set<CellID> >::const_iterator proc=recvs.begin(); proc!=recvs.end(); ++proc) {
	 // Calculate displacements:
	 displacements.clear();
	 for (std::set<CellID>::const_iterator cell=proc->second.begin(); cell!=proc->second.end(); ++cell) {
	    const CellID localID = global2LocalMap[*cell];
	    displacements.push_back(localID);
	 }
	 
	 // Create derived MPI datatype for receiving all cell IDs from 
	 // neighbour process proc and post a receive:
	 MPI_Datatype datatype;
	 MPI_Type_create_indexed_block(displacements.size(),1,&(displacements[0]),MPI_Type<CellID>(),&datatype);
	 MPI_Type_commit(&datatype);
	 MPI_Irecv(&(localIDs[0]),1,datatype,proc->first,proc->first,comm,&(recvRequests[recvCounter]));
	 MPI_Type_free(&datatype);
	 ++recvCounter;
      }
      
      // Calculate send displacements and post send(s):
      int sendCounter = 0;
      for (typename std::map<MPI_processID,std::set<CellID> >::const_iterator proc=sends.begin(); proc!=sends.end(); ++proc) {
	 // Calculate displacements:
	 displacements.clear();
	 for (std::set<CellID>::const_iterator cell=proc->second.begin(); cell!=proc->second.end(); ++cell) {
	    const CellID localID = global2LocalMap[*cell];
	    displacements.push_back(localID);
	 }
	 
	 // Create derived MPI datatype for sending all cell IDs to 
	 // neighbour process proc and post a send:
	 MPI_Datatype datatype;
	 MPI_Type_create_indexed_block(displacements.size(),1,&(displacements[0]),MPI_Type<CellID>(),&datatype);
	 MPI_Type_commit(&datatype);
	 MPI_Isend(&(localIDs[0]),1,datatype,proc->first,getRank(),comm,&(sendRequests[sendCounter]));
	 MPI_Type_free(&datatype);
	 
	 ++sendCounter;
      }
      
      // Wait for transfers to complete:
      MPI_Waitall(recvCounter,&(recvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(sendCounter,&(sendRequests[0]),MPI_STATUSES_IGNORE);
      return true;
   }

   /** Wait for remote cell data synchronization to complete.
    * @param stencilID ID of the stencil.
    * @param transferID ID of the transfer.
    * @return If true, data synchronization completed successfully.*/
   template<class C> inline
   bool ParGrid<C>::wait(StencilID stencilID,DataID userDataID) {
      if (getInitialized() == false) return false;
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator sten = stencils.find(stencilID);
      if (sten == stencils.end()) return false;
      return sten->second.wait(userDataID);
   }

   /** Wait for remote cell data synchronization to complete. If -DNDEBUG compiler flag 
    * is not defined, this function will
    * kill simulation execution if any process waits too long for data syncs to 
    * complete, printing error message with the given label. This is intented to 
    * be used for debbugging purposes. If -DNDEBUG compiler flag is defined, 
    * this function behaves similar to the other ParGrid::wait.
    * @param stencilID ID of the stencil.
    * @param transferID ID of the transfer.
    * @param name Label for the wait.
    * @return If true, data synchronization completed successfully.*/
   template<class C> inline
   bool ParGrid<C>::wait(StencilID stencilID,DataID userDataID,const std::string& name) {
      if (getInitialized() == false) return false;
      typename std::map<StencilID,Stencil<ParGrid<C>,C > >::iterator sten = stencils.find(stencilID);
      if (sten == stencils.end()) return false;
      #ifndef NDEBUG
         return sten->second.wait(userDataID,name);
      #else
         return sten->second.wait(userDataID);
      #endif
   }

   template<class C> inline
   void ParGrid<C>::waitMetadataRepartitioning(std::vector<MPI_Request>& nbrRecvRequests,
					       std::vector<MPI_Request>& nbrSendRequests) {
      // Wait for transfers to complete:
      MPI_Waitall(nbrRecvRequests.size(),&(nbrRecvRequests[0]),MPI_STATUSES_IGNORE);
      MPI_Waitall(nbrSendRequests.size(),&(nbrSendRequests[0]),MPI_STATUSES_IGNORE);
   }

   // ***************************************************** //
   // ***** DEFINITIONS FOR ZOLTAN CALLBACK FUNCTIONS ***** //
   // ***************************************************** //
   
   /** Definition of Zoltan callback function ZOLTAN_GEOM_MULTI_FN. This function is required by geometry-based 
    * load balancing (BLOCK,RCB,RIB,HSFC,REFTREE) functions. The purpose of this function is to tell Zoltan 
    * the physical x/y/z coordinates of all cells local to this process.
    * @param N_globalIDs Size of global ID.
    * @param N_localIDs Size of local ID.
    * @param N_cellIDs Number of cells whose coordinates are requested.
    * @param globalIDs Global IDs of cells whose coordinates are requested.
    * @param localIDs Local IDs of cells whose coordinates are requested.
    * @param N_coords Dimensionality of the mesh.
    * @param geometryData Array in which cell coordinates are to be written.
    * @param rcode Return code, upon successful exit value ZOLTAN_OK is written here.*/
   template<class C> inline
   void ParGrid<C>::getAllCellCoordinates(int N_globalIDs,int N_localIDs,int N_cellIDs,ZOLTAN_ID_PTR globalIDs,
					  ZOLTAN_ID_PTR localIDs,int N_coords,double* geometryData,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      const double* coordinates = reinterpret_cast<double*>(getUserData(C::getCoordinateDataID()));
      for (int i=0; i<N_cellIDs; ++i) {
	 geometryData[3*i+0] = coordinates[3*localIDs[i]+0];
	 geometryData[3*i+1] = coordinates[3*localIDs[i]+1];
	 geometryData[3*i+2] = coordinates[3*localIDs[i]+2];
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
   
   /** Definition of Zoltan callback function ZOLTAN_GEOM_FN. This function is required by
    * geometry-based load balancing (BLOCK,RCB,RIB,HSFC,Reftree). The purpose of this function is to tell
    * Zoltan the physical x/y/z coordinates of a given cell on this process.
    * @param N_globalEntries The size of array globalID.
    * @param N_localEntries The size of array localID.
    * @param globalID The global ID of the cell whose coordinates are queried.
    * @param localID The local ID of the cell whose coordinates are queried.
    * @param geometryData Array where the coordinates should be written. Zoltan has reserved this array, its size
    * is determined by the getGridDimensions function.
    * @param rcode The return code. Upon success this should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getCellCoordinates(int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
				       ZOLTAN_ID_PTR localID,CellCoordinate* geometryData,int* rcode) {
      // Check that the given cell exists on this process:
      #ifdef DEBUG_PARGRID
         if (localID[0] >= N_localCells) {
	    std::cerr << "ParGrid ERROR: getCellCoordinates queried non-existing cell #" << localID[0] << std::endl;
	    *rcode = ZOLTAN_FATAL;
	    return;
	 }
      #endif

      // Get cell coordinates from user:
      const double* coordinates = reinterpret_cast<double*>(getUserData(C::getCoordinateDataID()));
      geometryData[0] = coordinates[3*localID[0]+0];
      geometryData[1] = coordinates[3*localID[0]+1];
      geometryData[2] = coordinates[3*localID[0]+2];
      *rcode = ZOLTAN_OK;
   }

   /** Definition for Zoltan callback function ZOLTAN_EDGE_LIST_MULTI_FN. This function is required
    * for graph-based load balancing (GRAPH). The purpose of this function is to tell Zoltan
    * the global IDs of each neighbour of all cells local to this process, as well as the ranks of the 
    * MPI processes who own the neighbours and the weights of the edges (if edge weights are used).
    * @param N_globalIDs Size of global ID.
    * @param N_localIDs Size of local ID.
    * @param N_cells
    * @param globalIDs Global IDs cells whose neighbours are queried.
    * @param localIDs Local IDs cells whose neighbours are queried.
    * @param nbrGlobalIDs Array in which the global IDs of neighbours of all given cells are to be written.
    * @param nbrHosts For each neighbour, the rank of the MPI process who owns the cell.
    * @param N_weights The size of edge weight.
    * @param weight Array in which the weights of each edge are to be written.
    * @param rcode Return code, upon successful exit value ZOLTAN_OK is written here.*/
   template<class C> inline
   void ParGrid<C>::getAllCellEdges(int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
				    ZOLTAN_ID_PTR localIDs,int* N_edges,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
				    int N_weights,CellWeight* edgeWeights,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      size_t counter = 0;
      if (N_weights == 0) {
	 // Edge weights are not calculated
	 for (int i=0; i<N_cells; ++i) {
	    // Copy cell's neighbour information:
	    for (size_t n=0; n<N_neighbours; ++n) {
	       const CellID nbrLID = cellNeighbours[i*N_neighbours+n];
	       if (nbrLID == invalid()) continue;
	       
	       // Copy neighbour global ID and host process ID:
	       nbrGlobalIDs[counter] = this->globalIDs[nbrLID];
	       nbrHosts[counter]     = hosts[nbrLID];
	       ++counter;
	    }
	 }
      } else {
	 // Edge weights are calculated
	 for (int i=0; i<N_cells; ++i) {
	    // Copy cell's neighbour information:
	    for (size_t n=0; n<N_neighbours; ++n) {
	       const CellID nbrLID = cellNeighbours[localIDs[i]*N_neighbours+n];
	       if (nbrLID == invalid()) continue;
	       
	       // Copy neighbour global ID and host process ID:
	       nbrGlobalIDs[counter] = this->globalIDs[nbrLID];
	       nbrHosts[counter]     = hosts[nbrLID];
	       edgeWeights[counter]  = edgeWeight;
	       ++counter;
	    }
	 }	 
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
   
   /** Definition for Zoltan callback function ZOLTAN_EDGE_LIST_FN. This function is required
    * for graph-based load balancing (GRAPH). The purpose is to give the global IDs of each neighbour of
    * a given cell, as well as the ranks of the MPI processes which have the neighbouring cells.
    * @param N_globalIDs The size of array globalID.
    * @param N_localIDs The size of array localID.
    * @param globalID The global ID of the cell whose neighbours are queried.
    * @param localID The local ID of the cell whose neighbours are queried.
    * @param nbrGlobalIDs An array where the global IDs of the neighbours are written. Note that
    * Zoltan has already allocated this array based on a call to getNumberOfEdges function.
    * @param nbrHosts For each neighbour, the rank of the MPI process which has the cell.
    * @param N_weights The size of array weight.
    * @param weight The weight of each edge.
    * @param rcode The return code. Upon success should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getCellEdges(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
				   ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
				   int N_weights,CellWeight* weight,int* rcode) {
      // Count global IDs of cell's neighbours into Zoltan structures and calculate 
      // edge weight (if in use):
      int counter = 0;
      if (edgeWeightsUsed == true) {
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrLID = cellNeighbours[localID[0]*N_neighbours+n];
	    if (nbrLID == invalid()) continue;
	    nbrGlobalIDs[counter] = globalIDs[nbrLID];
	    nbrHosts[counter]     = hosts[nbrLID];
	    weight[counter]       = edgeWeight;
	    ++counter;
	 }
      } else {
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrLID = cellNeighbours[localID[0]*N_neighbours+n];
	    if (nbrLID == invalid()) continue;
	    nbrGlobalIDs[counter] = globalIDs[nbrLID];
	    nbrHosts[counter]     = hosts[nbrLID];
	    ++counter;
	 }
      }
      *rcode = ZOLTAN_OK;
   }
   
   template<class C> inline
   void ParGrid<C>::getHierarchicalParameters(int level,Zoltan_Struct* zs,int* rcode) {
      #ifndef DEBUG
         // Sanity check on input parameters:
         if (level < 0 || level >= (int)loadBalancingParameters.size()) {*rcode = ZOLTAN_FATAL; return;}
      #endif
      
      // Copy user-defined load balancing parameters to Zoltan structure:
      for (std::list<std::pair<std::string,std::string> >::const_iterator it=loadBalancingParameters[level].begin();
	   it!=loadBalancingParameters[level].end(); ++it) {
	 if (it->first == "PROCS_PER_PART") continue;
	 Zoltan_Set_Param(zs,it->first.c_str(),it->second.c_str());
      }
   }
   
   template<class C> inline
   int ParGrid<C>::getHierarchicalPartNumber(int level,int* rcode) {
      #ifndef NDEBUG
         if (level < 0 || level >= loadBalancingParameters.size()) {*rcode = ZOLTAN_FATAL; return -1;}
      #endif
      
      MPI_processID rank = getRank();
      for (int i=0; i<level; ++i) {
	 bool found = false;
	 int procsPerPart = 0;
	 for (std::list<std::pair<std::string,std::string> >::const_iterator it=loadBalancingParameters[i].begin();
	      it != loadBalancingParameters[i].end(); ++it) {
	    if (it->first == "PROCS_PER_PART") {
	       procsPerPart = atoi(it->second.c_str());
	       found = true;
	    }
	 }
	 
	 if (found == false) {
	    *rcode = ZOLTAN_FATAL;
	    return -1;
	 }
	 
	 if (i == level) {
	    rank = rank / procsPerPart;
	 } else {
	    rank = rank % procsPerPart;
	 }
      }
      
      *rcode = ZOLTAN_OK;
      return rank;
   }
   
   /** Definition for Zoltan callback function ZOLTAN_HG_CS_FN. This function is required for
    * hypergraph-based load balancing (HYPERGRAPH). The purpose is to give Zoltan the hypergraph in a compressed format.
    * @param N_globalIDs The size of globalID.
    * @param N_vtxedges The number of entries that need to be written to vtxedge_GID.
    * @param N_pins The number of pins that need to be written to pin_GID.
    * @param format The format that is used to represent the hypergraph, either ZOLTAN_COMPRESSED_EDGE or ZOLTAN_COMPRESSED_VERTEX.
    * @param vtxedge_GID An array where the hypergraph global IDs are written into.
    * @param vtxedge_ptr An array where, for each hyperedge, an index into pin_GID is given from where the pins for that
    * hyperedge are given.
    * @param pin_GID An array where the pins are written to.
    * @param rcode The return code. Upon success should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getHyperedges(int N_globalIDs,int N_vtxedges,int N_pins,int format,ZOLTAN_ID_PTR vtxedge_GID,
				  int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      // Check that correct hyperedge format is requested:
      if (format != ZOLTAN_COMPRESSED_VERTEX) {
	 *rcode = ZOLTAN_FATAL;
	 return;
      }
      // ------------ TARKISTA TOIMIIKO TÄMÄ ---------------- //
      // ONKO pinGID[pinCounter] OK ?????
      int pinCounter = 0;

      // Create list of hyperedges and pins:
      for (CellID i=0; i<N_localCells; ++i) {
	 vtxedge_GID[i]      = globalIDs[i];
	 vtxedge_ptr[i]      = pinCounter;
	 pin_GID[pinCounter] = globalIDs[i];
	 
	 // Add pin to this cell and to every existing neighbour:
	 for (size_t n=0; n<N_neighbours; ++n) {
	    const CellID nbrLID = cellNeighbours[i*N_neighbours+n];
	    if (nbrLID == invalid()) continue;
	    pin_GID[pinCounter] = globalIDs[nbrLID];
	    ++pinCounter;
	 }
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }

   /** Definition for Zoltan callback function ZOLTAN_HG_EDGE_WTS_FN. This is an optional function
    * for hypergraph-based load balancing (HYPEREDGE). The purpose is to tell Zoltan the weight of each hyperedge.
    * @param N_globalIDs The size of edgeGlobalID entry.
    * @param N_localIDs The size of edgeLocalID entry.
    * @param N_edges The number of hyperedge weights that need to be written to edgeWeights.
    * @param N_weights Number of weights per hyperedge that need to be written.
    * @param edgeGlobalID An array where the global IDs of each weight-supplying hyperedge are to be written.
    * @param edgeLocalID An array where the local IDs of each weight-supplying hyperedge are to be written.
    * This array can be left empty.
    * @param edgeWeights An array where the hyperedge weights are written into.
    * @param rcode The return code. Upon success should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getHyperedgeWeights(int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
					ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,CellWeight* edgeWeights,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      unsigned int counter = 0;
      
      if (edgeWeightsUsed == true) {
	 for (CellID i=0; i<N_localCells; ++i) {
	    edgeGlobalID[counter] = globalIDs[i];
	    edgeWeights[counter]  = edgeWeight;
	    ++counter;
	 }
      } else {
	 for (CellID i=0; i<N_localCells; ++i) {
	    edgeGlobalID[counter] = globalIDs[i];
	    ++counter;
	 }
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }

   /** Definition for Zoltan callback function ZOLTAN_OBJ_LIST_FN. This function
    * is required to use Zoltan. The purpose is to tell Zoltan the global and local
    * IDs of the cells assigned to this process, as well as their weights.
    * @param N_globalIDs The number of array entries used to describe one global ID.
    * @param N_localIDs The number of array entries used to describe one local ID.
    * @param globalIDs An array which is to be filled with the global IDs of the cells
    * currently assigned to this process. This array has been allocated by Zoltan.
    * @param localIDs An array which is to be filled with the local IDs of the cells
    * currently assigned to this process. This array has been allocated by Zoltan.
    * @param N_weights
    * @param cellWeights An array which is to be filled with the cell weights.
    * @param rcode The return code. Upon success this should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getLocalCellList(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
				     ZOLTAN_ID_PTR localIDs,int N_weights,CellWeight* cellWeights,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      #ifndef NDEBUG
         if (N_globalIDs != 1 || N_localIDs != 1) {
	    std::cerr << "(PARGRID) ERROR: Incorrect number of global/local IDs!" << std::endl;
	    exit(1);
	 }
         for (CellID i=0; i<N_localCells; ++i) {
	    if (this->cellWeights[i] < 0.0) {
	       std::cerr << "(PARGRID) ERROR: Negative cell weight in LID#" << globalIDs[i];
	       std::cerr << "\t GID#" << localIDs[i] << std::endl;
	    }
	 }      
      #endif
      
      if (N_weights == 1) {
	 // Iterate over all local cells, and get the cell weights from user. This 
	 // allows support for variable cell weights.
	 for (CellID i=0; i<N_localCells; ++i) {
	    globalIDs[i]   = this->globalIDs[i];
	    localIDs[i]    = i;
	    cellWeights[i] = this->cellWeights[i];
	 }
      } else {
	 // Iterate over all local cells and just copy global IDs to Zoltan structures:
	 for (CellID i=0; i<N_localCells; ++i) {
	    globalIDs[i] = this->globalIDs[i];
	    localIDs[i]  = i;
	 }
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
 
   /** Definition for Zoltan callback function ZOLTAN_NUM_GEOM_FN. This function is
    * required for geometry-based load balancing (BLOCK,RCB,RIB,HSFC,Reftree).
    * The purpose is 
    * to tell Zoltan the dimensionality of the grid. ParGrid always uses 
    * three-dimensional mesh internally.
    * @param parGridPtr A pointer to ParGrid.
    * @param rcode The return code. Upon success this should be ZOLTAN_OK.
    * @return The number of physical dimensions. Returns a value three.*/
   template<class C> inline
   int ParGrid<C>::getMeshDimension(int* rcode) {
      *rcode = ZOLTAN_OK;
      return 3;
   }

   /** Definition of Zoltan callback function ZOLTAN_NUM_EDGES_MULTI_FN. This function is required 
    * for graph-based load balancing (GRAPH). The purpose of this function is to tell Zoltan how 
    * many edges each cell local to this process has.
    * @param N_globalIDs Size of global ID.
    * @param N_localIDs Size of local ID.
    * @param N_cells Number of cells whose number of edges are requested.
    * @param globalIDs Global IDs of cells whose number of edges are requested.
    * @param localIDs Local IDs of cells whose number of edges are requested.
    * @param N_edges Array in which the number of edges each cell has are to be written.
    * @param rcode Return code, upon successful exit a value ZOLTAN_OK is written here.*/
   template<class C> inline
   void ParGrid<C>::getNumberOfAllEdges(int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
				       ZOLTAN_ID_PTR localIDs,int* N_edges,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      for (int i=0; i<N_cells; ++i) {
	 const CellID localID = localIDs[i];
	 int edgeSum = 0;
	 for (size_t n=0; n<N_neighbours; ++n) {
	    if (cellNeighbours[localID*N_neighbours+n] == invalid()) continue;
	    ++edgeSum;
	 }
	 N_edges[i] = edgeSum;
      }
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
   
   /** Definition of Zoltan callback function ZOLTAN_NUM_EDGES_FN. This function is required
    * for graph-based load balancing (GRAPH). The purpose is to tell how many edges a given cell has, i.e.
    * how many neighbours it has to share data with.
    * @param N_globalIDs The size of array globalID.
    * @param N_localIDs The size of array localID.
    * @param globalID The global ID of the cell whose edges are queried.
    * @param localID The local ID of the cell whose edges are queried.
    * @param rcode The return code. Upon success this should be ZOLTAN_OK.
    * @return The number of edges the cell has. For three-dimensional box grid this is between 3 and 6,
    * depending on if the cell is on the edge of the simulation volume.*/
   template<class C> inline
   int ParGrid<C>::getNumberOfEdges(int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
				    ZOLTAN_ID_PTR localID,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      // Count the number of neighbours the cell has:
      const CellID LID = localID[0];
      int edgeSum = 0;
      for (size_t n=0; n<N_neighbours; ++n) {
	 if (cellNeighbours[LID*N_neighbours+n] == invalid()) continue;
	 ++edgeSum;
      }
      
      // Return the number of edges:
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif      
      return edgeSum;
   }
   
   template<class C>
   int ParGrid<C>::getNumberOfHierarchicalLevels(int* rcode) {
      *rcode = ZOLTAN_OK;
      return loadBalancingParameters.size();
   }

   /** Definition for Zoltan callback function ZOLTAN_HG_SIZE_CS_FN. This function is required
    * for hypergraph-based load balancing (HYPERGRAPH). The purpose is to tell Zoltan which hypergraph format
    * is used (ZOLTAN_COMPRESSED_EDGE or ZOLTAN_COMPRESSED_VERTEX), how many hyperedges and
    * vertices there will be, and how many pins.
    * @param N_lists The total number of vertices or hyperedges (depending on the format)
    * is written to this variable.
    * @param N_pins The total number of pins (connections between vertices and hyperedges) is 
    * written to this variable.
    * @param format The chosen hyperedge storage format is written to this variable.
    * @param rcode The return code. Upon success should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getNumberOfHyperedges(int* N_lists,int* N_pins,int* format,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif

      *format = ZOLTAN_COMPRESSED_VERTEX;
      
      // Each local cell is a vertex:
      *N_lists = N_localCells;

      // Calculate the total number of pins:
      int totalNumberOfPins = 0;
      for (CellID i=0; i<N_localCells; ++i) {
	 for (size_t n=0; n<N_neighbours; ++n) {
	    if (cellNeighbours[i*N_neighbours+n] == invalid()) continue;
	    ++totalNumberOfPins;
	 }
      }
      *N_pins = totalNumberOfPins;
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
   
   /** Definition for Zoltan callback function ZOLTAN_HG_SIZE_EDGE_WTS_FN. This is an optional function
    * for hypergraph-based load balancing (HYPEREDGE). The purpose is to tell Zoltan how many hyperedges will have
    * a weight factor. Here we give a weight to each hyperedge.
    * @param parGridPtr A pointer to ParGrid.
    * @param N_edges A parameter where the number of weight-supplying hyperedges is written into.
    * @param rcode The return code. Upon success should be ZOLTAN_OK.*/
   template<class C> inline
   void ParGrid<C>::getNumberOfHyperedgeWeights(int* N_edges,int* rcode) {
      #ifdef PROFILE
         profile::start("ParGrid callbacks",profZoltanCB);
      #endif
      
      *N_edges = N_localCells;
      *rcode = ZOLTAN_OK;
      
      #ifdef PROFILE
         profile::stop();
      #endif
   }
   
   /** Definition for Zoltan callback function ZOLTAN_NUM_OBJ_FN. This function
    * is required to use Zoltan. The purpose is to tell Zoltan how many cells
    * are currently assigned to this process.
    * @param rcode The return code. Upon success this should be ZOLTAN_OK.
    * @return The number of cells assigned to this process.*/
   template<class C> inline
   int ParGrid<C>::getNumberOfLocalCells(int* rcode) {
      // Get the size of localCells container:
      *rcode = ZOLTAN_OK;
      return N_localCells;
   }
   
} // namespace pargrid

#endif
