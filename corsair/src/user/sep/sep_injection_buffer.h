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

#ifndef SEP_INJECTION_BUFFER_H
#define SEP_INJECTION_BUFFER_H

#include <vector>

namespace sep {
   template<class C>
   struct InjectionBuffer {
    public:
      void clear();
      void insert(const pargrid::CellID& cellID,const C& particle);
      
      std::vector<pargrid::CellID> cellIDs;
      std::vector<C> particles;
   };
   
   template<class C> inline
   void InjectionBuffer<C>::clear() {
      cellIDs.clear();
      particles.clear();
   }
   
   template<class C> inline
   void InjectionBuffer<C>::insert(const pargrid::CellID& cellID,const C& particle) {
      cellIDs.push_back(cellID);
      particles.push_back(particle);
   }
}
   
#endif
