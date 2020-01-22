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

#ifndef PARGRID_DEFINITIONS_H
#define PARGRID_DEFINITIONS_H

#include <cstdlib>
#include <stdint.h>
#include <limits>
#include <mpi.h>

namespace MPITypes {
   typedef int blocklength;                     /**< Datatype used by MPI in derived datatype blocklengths.*/
   typedef MPI_Aint displacement;               /**< Datatype used by MPI in derived datatype displacements.*/
   typedef int rank;                            /**< Datatype used by MPI to enumerate processes.*/
} // namespace MPI_Type

namespace pargrid {
   typedef float CellWeight;                        /**< Datatype Zoltan uses for cell and edge weights.*/
   typedef double CellCoordinate;                   /**< Datatype Zoltan uses for cell coordinates.*/
   typedef unsigned int CellID;                     /**< Datatype used for global cell IDs.*/
   typedef int MPI_processID;                       /**< Datatype MPI uses for process IDs.*/
   typedef unsigned char NeighbourID;               /**< Datatype ParGrid uses for neighbour type IDs.*/
   typedef unsigned int StencilID;                  /**< Datatype ParGrid uses for stencil IDs.*/
   typedef unsigned int DataID;                     /**< Datatype ParGrid uses for user-defined array IDs.*/
   
   enum InputParameter {
      cellWeightScale,                              /**< Value used for cell weights.*/
      edgeWeightScale,                              /**< Value (or its scale) used to calculate edge weights.*/
      imbalanceTolerance,                           /**< Load imbalance tolerance value.*/
      loadBalancingMethod,                          /**< Load balancing method to use.*/
      processesPerPartition                         /**< Processes per partition (hierarchical partitioning only).*/
   };
   
   enum StencilType {
      localToRemoteUpdates,
      remoteToLocalUpdates
   };
   
   /** Different partitioning modes that ParGrid supports.
    * These are directly tied to the partitioning modes supported by Zoltan.*/
   enum PartitioningMode {
      partition,                                    /**< New partitioning should always be computed from scratch.*/
      repartition,                                  /**< While repartitioning cells the existing partitioning
						     * should be taken into account to minimize data transfers.*/
      refine                                        /**< Partitioning calls should only slightly improve 
						     * the existing partitioning.*/
   };
   
   const size_t N_neighbours = 27;                      /**< Number of neighbours reserved each parallel cell has (including the cell itself).*/
   const uint32_t ALL_NEIGHBOURS_EXIST = 134217728 - 1; /**< If a cell's neighbour flags field equals this value all its 
							 * neighbours exist, i.e. the cell is an inner cell.*/
   const CellWeight DEFAULT_CELL_WEIGHT = 1.0e-06;      /**< Default cell weight (measures computational load).*/
   const StencilID DEFAULT_STENCIL = 0;                 /**< ID of the default transfer Stencil. This Stencil always exists.*/

   const CellID INVALID_CELLID       = std::numeric_limits<pargrid::CellID>::max();
   const DataID INVALID_DATAID       = std::numeric_limits<pargrid::DataID>::max();
   const StencilID INVALID_STENCILID = std::numeric_limits<pargrid::StencilID>::max();

   template<typename T> inline
   T getInvalidID() {return std::numeric_limits<T>::max();}
   
} // namespace pargrid
   
#endif