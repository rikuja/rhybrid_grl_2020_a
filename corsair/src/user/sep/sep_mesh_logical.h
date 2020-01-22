/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2013 Finnish Meteorological Institute
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

#ifndef SEP_MESH_LOGICAL_H
#define SEP_MESH_LOGICAL_H

#include <definitions.h>
#include <simulation.h>

namespace sep {

   namespace mesh {
      enum Type {
	 UNKNOWN,
	 LINEAR,                 /**< Mesh is linear.*/
	 LOGARITHMIC             /**< Mesh is logarithmic.*/
      };
   }
   
   void getLogicalCoordinates(Simulation* sim,const Real* physical,Real* logical);

   Real lambdaScalingInvRadius(const Real* pos);
   Real lambdaScalingNone(const Real* pos);
   Real lambdaScalingRadius(const Real* pos);
   Real lambdaScalingRadius2(const Real* pos);
   
} // namespace sep

#endif
