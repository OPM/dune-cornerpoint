//===========================================================================
//
// File: Intersection.hpp
//
// Created: Tue Jun  9 11:17:13 2009
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

#ifndef OPM_INTERSECTION_HEADER
#define OPM_INTERSECTION_HEADER




#include <dune/grid/common/gridenums.hh>

#include <opm/common/ErrorMacros.hpp>

// The next statement is a layering violation: we only #include
// preprocess.h to get at its "enum face_tag" definition.  Enum
// face_tag is needed in method Intersection::boundaryId().  This hack
// is in dire need of a better solution!
#include <opm/grid/cpgpreprocess/preprocess.h>

#include "Geometry.hpp"
#include "OrientedEntityTable.hpp"
namespace Dune
{
    namespace cpgrid
    {
    template<int>
    class Entity;
    class CpGridData;

        /// @brief
        /// @todo Doc me!
        /// @tparam
        class Intersection
        {
        public:
            /// @brief
            /// @todo Doc me!
            enum { dimension = 3 };
            enum { dimensionworld = 3 };
            /// @brief
            /// @todo Doc me!
            typedef cpgrid::Entity<0> Entity;
            typedef cpgrid::Geometry<2,3> Geometry;
            typedef cpgrid::Geometry<2,3> LocalGeometry;
            typedef double ctype;
            typedef FieldVector<ctype, 2> LocalCoordinate;
            typedef FieldVector<ctype, 3> GlobalCoordinate;

            /// @brief
            /// @todo Doc me!
            /// @param
            Intersection()
                : pgrid_(0),
                  index_(-1),
                  subindex_(-1),
                  faces_of_cell_(),
                  nbcell_(-1), // Init to self, which is invalid.
                  is_on_boundary_(false)
            {
            }
            /// @brief
            /// @todo Doc me!
            /// @param
            Intersection(const CpGridData& grid, const EntityRep<0>& cell, int subindex, bool update_now = true);

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            bool operator==(const Intersection& other) const
            {
                return subindex_ == other.subindex_  &&  index_ == other.index_  &&  pgrid_ == other.pgrid_;
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            bool operator!=(const Intersection& other) const
            {
                return !operator==(other);
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            bool boundary() const
            {
                return is_on_boundary_;
            }

            /// Returns the boundary id of this intersection.
            int boundaryId() const;


            /// Returns the boundary segment index of this intersection.
            int boundarySegmentIndex() const;

            /// @brief
            /// @todo Doc me!
            /// @return
            bool neighbor() const
            {
                return !boundary() && nbcell_!=std::numeric_limits<int>::max();
            }

            /// @brief
            /// @todo Doc me!
            /// @return
            Entity inside() const;

            /// @brief
            /// @todo Doc me!
            /// @return
            Entity outside() const;

            /// @brief
            /// @todo Doc me!
            /// @return
            bool conforming() const
            {
                return boundary(); // I.e. we are assuming all nonconforming interior.
            }

            // Geometrical information about this intersection in
            // local coordinates of the inside() entity.
            /// @brief
            /// @todo Doc me!
            /// @return
            const LocalGeometry& geometryInInside() const
            {
                OPM_THROW(std::runtime_error, "This intersection class does not support geometryInInside().");
            }

            // Geometrical information about this intersection in
            // local coordinates of the outside() entity.
            /// @brief
            /// @todo Doc me!
            /// @return
            const LocalGeometry& geometryInOutside() const
            {
                if (boundary()) {
                    OPM_THROW(std::runtime_error, "Cannot access geometryInOutside(), intersection is at a boundary.");
                }
                OPM_THROW(std::runtime_error, "This intersection class does not support geometryInOutside().");
            }

            /// @brief
            /// @todo Doc me!
            /// @return
            Geometry geometry() const;

            /// @brief
            /// @todo Doc me!
            /// @return
            GeometryType type() const
            {
                return geometry().type();
            }

            /// \brief Local index of codim 1 entity in the inside() entity
            /// where intersection is contained in.
            ///
            /// Note that CpGrid does not support codim-1 entities and return
            /// the index of the face of the underlying cartesian grid assuming that
            /// all faces are perpendicular to the axis:
            /// For faces perpendicular to x-axis: 0 (normal points contrary to axis) 1 (else)
            /// For faces perpendicular to x-axis: 2 (normal points contrary to axis) 3 (else)
            /// For faces perpendicular to x-axis: 4 (normal points contrary to axis) 5 (else)
            /// For faces not connecting neighboring elements: -1
            int indexInInside() const;

            /// \brief Local index of codim 1 entity in outside() entity
            /// where intersection is contained in.
            ///
            /// Note that CpGrid does not support codim-1 entities and return
            /// the index of the face of the underlying cartesian grid assuming that
            /// all faces are perpendicular to the axis:
            /// For faces perpendicular to x-axis: 0 (normal points contrary to axis) 1 (else)
            /// For faces perpendicular to x-axis: 2 (normal points contrary to axis) 3 (else)
            /// For faces perpendicular to x-axis: 4 (normal points contrary to axis) 5 (else)
            /// For faces not connecting neighboring elements: -1
            int indexInOutside() const
            {
                int in_inside = indexInInside();
                if (in_inside == -1) {
                    // NNC face, return -1 here as well.
                    return -1;
                }
                return in_inside + ((in_inside % 2) ? -1 : 1);
            }

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            FieldVector<ctype, 3> outerNormal(const FieldVector<ctype, 2>&) const;

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            FieldVector<ctype, 3> integrationOuterNormal(const FieldVector<ctype, 2>& unused) const;

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            FieldVector<ctype, 3> unitOuterNormal(const FieldVector<ctype, 2>&) const;

            /// @brief
            /// @todo Doc me!
            /// @param
            /// @return
            FieldVector<ctype, 3> centerUnitOuterNormal() const;

            int id() const
            {
                const EntityRep<1>& face = faces_of_cell_[subindex_];
                return face.index();
            }

        protected:
            const CpGridData* pgrid_;
            int index_;
            int subindex_;
            OrientedEntityTable<0,1>::row_type faces_of_cell_;
            int nbcell_;
            bool is_on_boundary_;

            void increment();

            void update();

            void setAtEnd()
            {
                subindex_ = faces_of_cell_.size();
            }

            bool isAtEnd() const
            {
                return subindex_ == faces_of_cell_.size();
            }

            int nbcell() const
            {
                if (is_on_boundary_) {
                    OPM_THROW(std::runtime_error, "There is no outside cell, intersection is at boundary.");
                }
                if(nbcell_==std::numeric_limits<int>::max())
                    OPM_THROW(std::runtime_error, "There is no outside cell, intersection is at processor boundary.");
                return nbcell_;
            }
        };





        class IntersectionIterator : public Intersection
        {
        public:
            typedef cpgrid::Intersection Intersection;

            IntersectionIterator()
                : Intersection()
            {
            }

            IntersectionIterator(const CpGridData& grid, const EntityRep<0>& cell, bool at_end)
                : Intersection(grid, cell, 0, !at_end)
            {
                if (at_end) {
                    Intersection::setAtEnd();
                } else {
                    Intersection::update();
                }
            }

            IntersectionIterator& operator++()
            {
                Intersection::increment();
                return *this;
            }

            const Intersection* operator->() const
            {
                assert(!Intersection::isAtEnd());
                return this;
            }

            const Intersection& operator*() const
            {
                assert(!Intersection::isAtEnd());
                return *this;
            }

        };





    } // namespace cpgrid
} // namespace Dune

#endif // OPM_INTERSECTION_HEADER
