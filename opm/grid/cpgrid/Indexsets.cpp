//===========================================================================
//
// File: Indexsets.hpp
//
// Created: Fri May 29 23:30:01 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
Copyright 2009, 2010, 2022 Equinor ASA.

This file is part of The Open Porous Media project  (OPM).

OPM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OPM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Entity.hpp"
#include "Indexsets.hpp"
#include "CpGridData.hpp"

namespace Dune
{
namespace cpgrid
{
IndexSet::IndexType IndexSet::subIndex(const cpgrid::Entity<0>& e, int i, grid_size_t cc) const
{
    switch(cc) {
    case 0: return index(e.subEntity<0>(i));
    case 1: return index(e.subEntity<1>(i));
    case 2: return index(e.subEntity<2>(i));
    case 3: return index(e.subEntity<3>(i));
    default: OPM_THROW(std::runtime_error,
                       "Codimension " + std::to_string(cc) + " not supported.");
    }
}

IdSet::IdType IdSet::subId(const cpgrid::Entity<0>& e, int i, int cc) const
{
    switch (cc) {
    case 0: return id(e.subEntity<0>(i));
    case 1: return id(e.subEntity<1>(i));
    case 2: return id(e.subEntity<2>(i));
    case 3: return id(e.subEntity<3>(i));
    default: OPM_THROW(std::runtime_error,
                       "Cannot get subId of codimension " + std::to_string(cc));
    }
    return -1;
}

LevelGlobalIdSet::IdType LevelGlobalIdSet::subId(const cpgrid::Entity<0>& e, int i, int cc) const
{
    assert(view_ == e.pgrid_);

    switch (cc) {
    case 0: return id(e.subEntity<0>(i));
        //case 1: return id(*e.subEntity<1>(i));
        //case 2: return id(*e.subEntity<2>(i));
    case 3: return id(e.subEntity<3>(i));
    default: OPM_THROW(std::runtime_error,
                       "Cannot get subId of codimension " + std::to_string(cc));
    }
    return -1;
}

GlobalIdSet::IdType GlobalIdSet::subId(const cpgrid::Entity<0>& e, int i, int cc) const
{
    return levelIdSet(e.pgrid_).subId(e, i, cc);
}

void GlobalIdSet::insertIdSet(const CpGridData& view)
{
    idSets_.insert(std::make_pair(&view,view.global_id_set_));
}
GlobalIdSet::GlobalIdSet(const CpGridData& view)
{
    idSets_.insert(std::make_pair(&view,view.global_id_set_));
}
} // end namespace cpgrid
} // end namespace Dune
