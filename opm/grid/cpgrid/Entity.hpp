//===========================================================================
//
// File: Entity.hpp
//
// Created: Fri May 29 20:26:48 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Brd Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPM_ENTITY_HEADER
#define OPM_ENTITY_HEADER

#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>

#include "PartitionTypeIndicator.hpp"
#include "EntityRep.hpp"

namespace Dune
{
    namespace cpgrid
    {

       template<int,int> class Geometry;
       template<int,PartitionIteratorType> class Iterator;
       class IntersectionIterator;
       class HierarchicIterator;
       class CpGridData;
       class LevelGlobalIdSet;


        /// @brief
        /// @todo Doc me!
        /// @tparam
        template <int codim>
        class Entity : public EntityRep<codim>
        {
            friend class LevelGlobalIdSet;
            friend class GlobalIdSet;

        public:
        /// @brief
        /// @todo Doc me!
            enum { codimension = codim };
            enum { dimension = 3 };
            enum { mydimension = dimension - codimension };
            enum { dimensionworld = 3 };

            // the official DUNE names
            typedef Entity    EntitySeed;

            /// @brief
            /// @todo Doc me!
            /// @tparam
            template <int cd>
            struct Codim
            {
                typedef cpgrid::Entity<cd> Entity;
            };


            typedef cpgrid::Geometry<3-codim,3> Geometry;
            typedef Geometry LocalGeometry;

            typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
            typedef cpgrid::IntersectionIterator LevelIntersectionIterator;
            typedef cpgrid::HierarchicIterator HierarchicIterator;

            typedef double ctype;

            /// Constructor taking a grid and an integer entity representation.
            /// This constructor should probably be removed, since it exposes
            /// details of the implementation of \see EntityRep, see comment in
            /// EntityRep<>::EntityRep(int).
//             Entity(const CpGridData& grid, int entityrep)
//              : EntityRep<codim>(entityrep), pgrid_(&grid)
//          {
//          }

            /// Constructor creating empty entity
            Entity()
                : EntityRep<codim>(), pgrid_( 0 )
            {
            }

            /// Constructor taking a grid and an entity representation.
            Entity(const CpGridData& grid, EntityRep<codim> entityrep)
                : EntityRep<codim>(entityrep), pgrid_(&grid)
            {
            }

            /// Constructor taking a grid, entity index and orientation.
            Entity(const CpGridData& grid, int index_arg, bool orientation_arg)
                : EntityRep<codim>(index_arg, orientation_arg), pgrid_(&grid)
            {
            }

            /// Equality.
            bool operator==(const Entity& other) const
            {
                return EntityRep<codim>::operator==(other)  &&  pgrid_ == other.pgrid_;
            }

            /// Inequality.
            bool operator!=(const Entity& other) const
            {
                return !operator==(other);
            }

            /// Return an entity seed.
            /// For CpGrid, EntitySeed and EntityPtr are the same class.
            EntitySeed seed() const
            {
                return EntitySeed( impl() );
            }

            /// Returns the geometry of the entity (does not depend on its orientation).
            Geometry geometry() const;

            /// We do not support refinement, so level() is always 0.
            int level() const
            {
                return 0;
            }

            /// The entity is always on the leaf grid, since we have no refinement.
            bool isLeaf() const
            {
                return true;
            }

            /// Refinement is not defined for CpGrid.
            bool isRegular() const
            {
                return true;
            }

            /// For now, the grid is serial and the only partitionType() is InteriorEntity.
            PartitionType partitionType() const;

            /// Using the cube type for all entities now (cells and vertices).
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(3 - codim);
            }

            /// The count of subentities of codimension cc
            unsigned int subEntities ( const unsigned int cc ) const;
            
            /// Obtain subentity.
            template <int cc>
            typename Codim<cc>::Entity subEntity(int i) const;

            /// Start iterator for the cell-cell intersections of this entity.
            inline LevelIntersectionIterator ilevelbegin() const;

            /// End iterator for the cell-cell intersections of this entity.
            inline LevelIntersectionIterator ilevelend() const;

            /// Start iterator for the cell-cell intersections of this entity.
            inline LeafIntersectionIterator ileafbegin() const;

            /// End iterator for the cell-cell intersections of this entity.
            inline LeafIntersectionIterator ileafend() const;

            /// Dummy first child iterator.
            HierarchicIterator hbegin(int) const;

            /// Dummy beyond last child iterator.
            HierarchicIterator hend(int) const;

            /// \brief Returns true, if the entity has been created during the last call to adapt(). Dummy.
            bool isNew() const
            {
                return false;
            }

            /// \brief Returns true, if entity might disappear during the next call to adapt(). Dummy.
            bool mightVanish() const
            {
                return false;
            }

            /// No hierarchy in this grid.
            bool hasFather() const
            {
                return false;
            }


            /// Dummy, returning this.
            Entity father() const
            {
                return *this;
            }


            /// Dummy, returning default geometry.
            LocalGeometry geometryInFather() const
            {
                return LocalGeometry();
            }

            /// Returns true if any of my intersections are on the boundary.
            /// Implementation note:
            /// This is a slow, computed, function. Could be speeded
            /// up by putting boundary info in the CpGrid class.
            bool hasBoundaryIntersections() const;

            // Mimic Dune entity wrapper

            const Entity& impl() const
            {
                return *this;
            }

            Entity& impl()
            {
                return *this;
            }

            /// isValid method for EntitySeed
            /// \return return true if seed is pointing to a valid entity
            bool isValid () const;

        protected:
            const CpGridData* pgrid_;
        };

    } // namespace cpgrid
} // namespace Dune

// now we include the Iterators.hh We need to do this here because for hbegin/hend the compiler
// needs to know the size of hierarchicIterator
#include "Iterators.hpp"
#include "Entity.hpp"
#include "Intersection.hpp"
namespace Dune
{
  namespace cpgrid
  {
  template<int codim>
  typename Entity<codim>::LevelIntersectionIterator Entity<codim>::ilevelbegin() const
  {
                static_assert(codim == 0, "");
                return LevelIntersectionIterator(*pgrid_, *this, false);
            }

  template<int codim>
  typename Entity<codim>::LevelIntersectionIterator Entity<codim>::ilevelend() const
  {
      static_assert(codim == 0, "");
      return LevelIntersectionIterator(*pgrid_, *this, true);
  }

  template<int codim>
  typename Entity<codim>::LeafIntersectionIterator Entity<codim>::ileafbegin() const
  {
      static_assert(codim == 0, "");
      return LeafIntersectionIterator(*pgrid_, *this, false);
  }

  template<int codim>
  typename Entity<codim>::LeafIntersectionIterator Entity<codim>::ileafend() const
  {
      static_assert(codim == 0, "");
      return LeafIntersectionIterator(*pgrid_, *this, true);
  }


    template<int codim>
    HierarchicIterator Entity<codim>::hbegin(int) const
    {
        return HierarchicIterator(*pgrid_);
    }

    /// Dummy beyond last child iterator.
    template<int codim>
    HierarchicIterator Entity<codim>::hend(int) const
    {
        return HierarchicIterator(*pgrid_);
    }

    template <int codim>
    PartitionType Entity<codim>::partitionType() const
    {
        return pgrid_->partition_type_indicator_->getPartitionType(*this);
    }
    } // namespace cpgrid
} // namespace Dune

#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Intersection.hpp>

namespace Dune {
namespace cpgrid {

template<int codim>
unsigned int Entity<codim>::subEntities ( const unsigned int cc ) const
{
    if (cc == codim) {
        return 1;
    } else if ( codim == 0 ){
        if ( cc == 3 ) {
            return 8;
        }
    }
    return 0;
}

template <int codim>
typename Entity<codim>::Geometry Entity<codim>::geometry() const
{
    return pgrid_->geomVector<codim>()[*this];
}

template <int codim>
template <int cc>
typename Entity<codim>::template Codim<cc>::Entity Entity<codim>::subEntity(int i) const
{
    static_assert(codim == 0, "");
    if (cc == 0) {
        assert(i == 0);
        typename Codim<cc>::Entity se(*pgrid_, EntityRep<codim>::index(), EntityRep<codim>::orientation());
        return se;
    } else if (cc == 3) {
        assert(i >= 0 && i < 8);
        int corner_index = pgrid_->cell_to_point_[EntityRep<codim>::index()][i];
        typename Codim<cc>::Entity se(*pgrid_, corner_index, true);
        return se;
    } else {
        OPM_THROW(std::runtime_error,
                  "No subentity exists of codimension " + std::to_string(cc));
    }
}

template <int codim>
bool Entity<codim>::hasBoundaryIntersections() const
{
    // Copied implementation from EntityDefaultImplementation,
    // except for not checking LevelIntersectionIterators.
    typedef LeafIntersectionIterator Iter;
    Iter end = ileafend();
    for (Iter it = ileafbegin(); it != end; ++it) {
        if (it->boundary()) return true;
    }
    return false;
}

template <int codim>
bool Entity<codim>::isValid() const
{
    return pgrid_ ? EntityRep<codim>::index() < pgrid_->size(codim) : false;
}

}}


#endif // OPM_ENTITY_HEADER
