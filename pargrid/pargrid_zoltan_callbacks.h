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

#ifndef PARGRID_ZOLTAN_CALLBACKS_H
#define PARGRID_ZOLTAN_CALLBACKS_H

#include <zoltan_cpp.h>

namespace pargrid {
   // ************************************************* //
   // ***** ZOLTAN CALLBACK FUNCTION DECLARATIONS ***** //
   // ************************************************* //
   template<class PARGRID>
   void cb_getAllCellCoordinates(void* pargrid,int N_globalEntries,int N_localEntries,int N_cellIDs,ZOLTAN_ID_PTR globalID,
				 ZOLTAN_ID_PTR localID,int N_coords,double* geometryData,int* rcode);
   
   template<class PARGRID>
   void cb_getCellCoordinates(void* pargrid,int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
			      ZOLTAN_ID_PTR localID,double* geometryData,int* rcode);
   
   template<class PARGRID>
   void cb_getAllCellEdges(void* pargrid,int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			   ZOLTAN_ID_PTR localIDs,int* N_edges,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			   int N_weights,CellWeight* edgeWeights,int* rcode);
   
   template<class PARGRID>
   void cb_getCellEdges(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			int N_weights,CellWeight* weight,int* rcode);
   
   template<class PARGRID>
   void cb_getHierarchicalParameters(void* pargrid,int level,Zoltan_Struct* zs,int* rcode);
   
   template<class PARGRID>
   int cb_getHierarchicalPartNumber(void* pargrid,int level,int* rcode);
   
   template<class PARGRID>
   void cb_getHyperedges(void* pargrid,int N_globalIDs,int N_vtxedges,int N_pins,int format,
			 ZOLTAN_ID_PTR vtxedge_GID,int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode);
   
   template<class PARGRID>
   void cb_getHyperedgeWeights(void* pargrid,int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
			       ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,CellWeight* edgeWeights,int* rcode);

   template<class PARGRID>
   int cb_getMeshDimension(void* pargrid,int* rcode);

   template<class PARGRID>
   void cb_getLocalCellList(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
			    ZOLTAN_ID_PTR localIDs,int N_weights,CellWeight* cellWeights,int* rcode);
   
   template<class PARGRID>
   void cb_getNumberOfAllEdges(void* pargrid,int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			       ZOLTAN_ID_PTR localIDs,int* N_edges,int* rcode);

   template<class PARGRID>
   int cb_getNumberOfEdges(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			   ZOLTAN_ID_PTR localID,int* rcode);
   
   template<class PARGRID>
   int cb_getNumberOfHierarchicalLevels(void* pargrid,int* rcode);
   
   template<class PARGRID>
   void cb_getNumberOfHyperedges(void* pargrid,int* N_lists,int* N_pins,int* format,int* rcode);
   
   template<class PARGRID>
   void cb_getNumberOfHyperedgeWeights(void* pargrid,int* N_edges,int* rcode);
   
   template<class PARGRID>
   int cb_getNumberOfLocalCells(void* pargrid,int* rcode);
   
   // ************************************************ //
   // ***** ZOLTAN CALLBACK FUNCTION DEFINITIONS ***** //
   // ************************************************ //
   
   template<class PARGRID> inline
   void cb_getAllCellCoordinates(void* pargrid,int N_globalEntries,int N_localEntries,int N_cellIDs,ZOLTAN_ID_PTR globalID,
				 ZOLTAN_ID_PTR localID,int N_coords,double* geometryData,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getAllCellCoordinates(N_globalEntries,N_localEntries,N_cellIDs,globalID,localID,N_coords,geometryData,rcode);
   }
   
   template<class PARGRID> inline
   void cb_getCellCoordinates(void* pargrid,int N_globalEntries,int N_localEntries,ZOLTAN_ID_PTR globalID,
			      ZOLTAN_ID_PTR localID,CellCoordinate* geometryData,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getCellCoordinates(N_globalEntries,N_localEntries,globalID,localID,geometryData,rcode);
   }

   template<class PARGRID> inline
   void cb_getAllCellEdges(void* pargrid,int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			   ZOLTAN_ID_PTR localIDs,int* N_edges,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			   int N_weights,CellWeight* edgeWeights,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getAllCellEdges(N_globalIDs,N_localIDs,N_cells,globalIDs,localIDs,N_edges,
			   nbrGlobalIDs,nbrHosts,N_weights,edgeWeights,rcode);
   }

   template<class PARGRID> inline
   void cb_getCellEdges(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,
			ZOLTAN_ID_PTR localID,ZOLTAN_ID_PTR nbrGlobalIDs,int* nbrHosts,
			int N_weights,CellWeight* weight,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getCellEdges(N_globalIDs,N_localIDs,globalID,localID,nbrGlobalIDs,nbrHosts,N_weights,weight,rcode);
   }

   template<class PARGRID> inline
   void cb_getHierarchicalParameters(void* pargrid,int level,Zoltan_Struct* zs,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getHierarchicalParameters(level,zs,rcode);
   }

   template<class PARGRID> inline
   int cb_getHierarchicalPartNumber(void* pargrid,int level,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      return ptr->getHierarchicalPartNumber(level,rcode);
   }

   template<class PARGRID> inline
   void cb_getHyperedges(void* pargrid,int N_globalIDs,int N_vtxedges,int N_pins,int format,ZOLTAN_ID_PTR vtxedge_GID,
			 int* vtxedge_ptr,ZOLTAN_ID_PTR pin_GID,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getHyperedges(N_globalIDs,N_vtxedges,N_pins,format,vtxedge_GID,vtxedge_ptr,pin_GID,rcode);
   }

   template<class PARGRID> inline
   void cb_getHyperedgeWeights(void* pargrid,int N_globalIDs,int N_localIDs,int N_edges,int N_weights,
			       ZOLTAN_ID_PTR edgeGlobalID,ZOLTAN_ID_PTR edgeLocalID,CellWeight* edgeWeights,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getHyperedgeWeights(N_globalIDs,N_localIDs,N_edges,N_weights,edgeGlobalID,edgeLocalID,edgeWeights,rcode);
   }

   template<class PARGRID> inline
   void cb_getLocalCellList(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalIDs,
			    ZOLTAN_ID_PTR localIDs,int N_weights,CellWeight* cellWeights,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getLocalCellList(N_globalIDs,N_localIDs,globalIDs,localIDs,N_weights,cellWeights,rcode);
   }

   template<class PARGRID> inline
   int cb_getMeshDimension(void* pargrid,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      return ptr->getMeshDimension(rcode);
   }

   template<class PARGRID> inline
   void cb_getNumberOfAllEdges(void* pargrid,int N_globalIDs,int N_localIDs,int N_cells,ZOLTAN_ID_PTR globalIDs,
			       ZOLTAN_ID_PTR localIDs,int* N_edges,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getNumberOfAllEdges(N_globalIDs,N_localIDs,N_cells,globalIDs,localIDs,N_edges,rcode);
   }

   template<class PARGRID> inline
   int cb_getNumberOfEdges(void* pargrid,int N_globalIDs,int N_localIDs,ZOLTAN_ID_PTR globalID,ZOLTAN_ID_PTR localID,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      return ptr->getNumberOfEdges(N_globalIDs,N_localIDs,globalID,localID,rcode);
   }

   template<class PARGRID> inline
   int cb_getNumberOfHierarchicalLevels(void* pargrid,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      return ptr->getNumberOfHierarchicalLevels(rcode);
   }

   template<class PARGRID> inline
   void cb_getNumberOfHyperedges(void* pargrid,int* N_lists,int* N_pins,int* format,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getNumberOfHyperedges(N_lists,N_pins,format,rcode);
   }

   template<class PARGRID> inline
   void cb_getNumberOfHyperedgeWeights(void* pargrid,int* N_edges,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      ptr->getNumberOfHyperedgeWeights(N_edges,rcode);
   }

   template<class PARGRID> inline
   int cb_getNumberOfLocalCells(void* pargrid,int* rcode) {
      PARGRID* ptr = reinterpret_cast<PARGRID*>(pargrid);
      return ptr->getNumberOfLocalCells(rcode);
   }
   
} // namespace pargrid

#endif