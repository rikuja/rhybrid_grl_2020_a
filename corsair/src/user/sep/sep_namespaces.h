/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#ifndef SEP_NAMESPACES_H
#define SEP_NAMESPACES_H

namespace sep {

   namespace particleinjector {
      namespace frame {
	 enum Frame {
	    UNKNOWN,
	    SIMULATION,                              /**< Injection in simulation frame.*/
	    ANTIPARALLEL,                            /**< Injection in antiparallel Alfven wave rest frame.*/
	    PARALLEL,                                /**< Injection in parallel Alfven rest frame.*/
	    PLASMA,                                  /**< Injection in plasma rest frame.*/
	    SIZE
	 };
      }
      
      namespace region {
	 enum Region {
	    UNKNOWN,
	    UPSTREAM,
	    DOWNSTREAM,
	    SHOCK_UPSTREAM,
	    SHOCK_DOWNSTREAM,
	    SIZE
	 };
      }
   }
   
   namespace lagrinjector {
      enum InflowBoundary {
	 X_NEG,                                   /**< Wave packets always injected to -x boundary.*/
	 X_POS,                                   /**< Wave packets always injected to +x boundary.*/
	 VELOCITY_DIR                             /**< Wave packets injected when wave velocity vector points
						   * towards simulation domain.*/
      };
   }
   
} // namespace sep

#endif
