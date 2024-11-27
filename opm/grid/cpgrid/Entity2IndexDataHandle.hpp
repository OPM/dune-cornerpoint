//===========================================================================
//
// File: Entity2IndexDataHandle.hpp
//
// Created: Mon Nov 4 2013
//
// Author(s): Markus Blatt <markus@dr-blatt.de>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2013 Statoil ASA.
  Copyright 2013 Dr. Markus Blatt - HPC-Simulation-Software & Services

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
#ifndef OPM_ENTITY2INDEXDATAHANDLE_HEADER
#define OPM_ENTITY2INDEXDATAHANDLE_HEADER

#include <dune/common/version.hh>

#include <cstddef>

namespace Dune
{
namespace cpgrid
{
template<int codim> class Entity;
class CpGridData;

/// \brief Wrapper that turns a data handle suitable for dune-grid into one based on
/// integers instead of entities.
///
/// \tparam DataHandle The type of the data handle to wrap. Has to adhere to the interface
/// of Dune::DataHandleIf
///  \tparam codim The codimension to use when mapping indices to Entities.
template<class DataHandle, int codim>
class Entity2IndexDataHandle
{
public:
    typedef typename DataHandle::DataType DataType;

    Entity2IndexDataHandle(const CpGridData& grid, DataHandle& data)
        : fromGrid_(grid), toGrid_(grid), data_(data)
    {}
    Entity2IndexDataHandle(const CpGridData& fromGrid, const CpGridData& toGrid, DataHandle& data)
        : fromGrid_(fromGrid), toGrid_(toGrid), data_(data)
    {}

#if DUNE_VERSION_GTE(DUNE_COMMON, 2, 8)
    bool fixedSize()
    {
        return data_.fixedSize(3, codim);
    }
#else
    bool fixedsize()
    {
        return data_.fixedSize(3, codim);
    }
#endif
    std::size_t size(std::size_t i)
    {
        return data_.size(Entity<codim>(fromGrid_, i, true));
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        data_.gather(buffer, Entity<codim>(fromGrid_, i, true));
    }
    template<class B>
    void scatter(B& buffer, std::size_t i, std::size_t s)
    {
        data_.scatter(buffer, Entity<codim>(toGrid_, i, true), s);
    }
private:
    const CpGridData& fromGrid_;
    const CpGridData& toGrid_;
    DataHandle& data_;

};
} // end namespace cpgrid
} // end namespace Dune
#endif
