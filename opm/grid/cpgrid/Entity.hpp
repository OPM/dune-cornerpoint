//===========================================================================
//
// File: Entity.hpp
//
// Created: Fri May 29 20:26:48 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Brd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
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

#include <opm/grid/UnstructuredGrid.h>
#include "PartitionTypeIndicator.hpp"
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>

// To be able to test local and global ids of vertices
void refinePatch_and_check(Dune::CpGrid&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::string>&);

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
    friend class HierarchicIterator;
    friend class CpGridData;
    friend void ::refinePatch_and_check(Dune::CpGrid&,
                                        const std::vector<std::array<int,3>>&,
                                        const std::vector<std::array<int,3>>&,
                                        const std::vector<std::array<int,3>>&,
                                        const std::vector<std::string>&);

public:
    /// @brief
    /// @todo Doc me!
    static constexpr int codimension = codim;
    static constexpr int dimension = 3;
    static constexpr int mydimension = dimension - codimension;
    static constexpr int dimensionworld = 3;

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

    /// Constructor taking a grid, entity index, and orientation.
    Entity(const CpGridData& grid, int index_arg, bool orientation_arg)
        : EntityRep<codim>(index_arg, orientation_arg), pgrid_(&grid)
    {
    }

    /// Constructor taking a entity index, and orientation.
    Entity(int index_arg, bool orientation_arg)
        : EntityRep<codim>(index_arg, orientation_arg), pgrid_()
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

    /// @brief Return an entity seed (light-weight entity).
    ///        EntitySeed objects are used to obtain an Entity back when combined with the corresponding grid.
    ///        For CpGrid, EntitySeed and EntityPtr are the same class.
    EntitySeed seed() const
    {
        return EntitySeed( impl() );
    }

    /// @brief Return the geometry of the entity (does not depend on its orientation).
    const Geometry& geometry() const;

    /// @brief Return the level of the entity in the grid hierarchy. Level = 0 represents the coarsest grid.
    int level() const;

    /// @brief Check if the entity is in the leafview.
    ///
    ///        @TODO: Modify the definition to cover serial and parallel cases.
    ///        Serial: an element is a leaf <-> hbegin and hend return the same iterator
    ///        Parallel: true <-> the element is a leaf entity of the global refinement hierarchy.
    bool isLeaf() const;

    /// Refinement is not defined for CpGrid.
    bool isRegular() const
    {
        return true;
    }

    /// @brief In serial run, the only partitionType() is InteriorEntity.
    ///        Only needed when distributed_data_ is not empty.
    ///        In parallel run, the entity each cell inherits the partition type of its origin
    ///        (either parent cell or equivalent cell in level 0).
    ///        Only needed when distributed_data_ is not empty.
    PartitionType partitionType() const;

    /// @brief Return marker object (GeometryType object) representing the reference element of the entity.
    ///        Currently, cube type for all entities (cells and vertices).
    GeometryType type() const
    {
        return Dune::GeometryTypes::cube(3 - codim);
    }

    /// @brief Return the number of all subentities of the entity of a given codimension cc.
    grid_size_t subEntities ( const grid_size_t cc ) const;

    /// @brief Obtain subentity.
    ///        Example: If cc = 3 and i = 5, it returns the 5th corner/vertex of the entity.
    template <int cc>
    typename Codim<cc>::Entity subEntity(int i) const;

    /// Start level-iterator for the cell-cell intersections of this entity.
    inline LevelIntersectionIterator ilevelbegin() const;

    /// End level-iterator for the cell-cell intersections of this entity.
    inline LevelIntersectionIterator ilevelend() const;

    /// Start leaf-iterator for the cell-cell intersections of this entity.
    inline LeafIntersectionIterator ileafbegin() const;

    /// End leaf-iterator for the cell-cell intersections of this entity.
    inline LeafIntersectionIterator ileafend() const;


    /// @brief Iterator begin over the children. [If requested, also over descendants more than one generation away.]
    HierarchicIterator hbegin(int) const;

    /// @brief Iterator end over the children/beyond last child iterator.
    HierarchicIterator hend(int) const;

    /// \brief Returns true, if the entity has been created during the last call to adapt(). Dummy.
    bool isNew() const;

    /// \brief Returns true, if entity might disappear during the next call to adapt().
    ///
    ///        Currently, returns true if the element has been marked to be refined. Therefore, it "vanishes"
    ///        in the sense that it gets replaced by its children.
    ///        To be modified: When supporting coarsening, mightVanish should return true for all the "siblings"
    ///        refined cells of one refined cell that got marked to be coarsened.
    bool mightVanish() const;

    /// @brief ONLY FOR CELLS (Entity<0>)
    ///        Check if the entity comes from an LGR, i.e., it has been created via refinement from coarser level.
    ///
    ///        @TODO: When distributed_data_ is not empty, check whether the father element exists on the
    ///        local process, which can be used to test whether it is safe to call father.
    bool hasFather() const;

    /// @brief  ONLY FOR CELLS (Entity<0>). Get the father Entity (in the level-grid the father cell was born), in case entity.hasFather() is true.
    ///
    /// @return father-entity
    Entity<0> father() const;

    /// @brief Return LocalGeometry representing the embedding of the entity into its father (when hasFather() is true).
    ///        Map from the entity's reference element into the reference element of its father.
    ///        Currently, LGR is built via refinement of a block-shaped patch from the coarse grid. So the LocalGeometry
    ///        of an entity coming from the LGR is one of the refined cells of the unit cube, with suitable amount of cells
    ///        in each direction.
    Dune::cpgrid::Geometry<3,3> geometryInFather() const;

    /// Returns true if any of my intersections are on the boundary.
    /// Implementation note:
    /// This is a slow, computed, function. Could be speeded
    /// up by putting boundary info in the CpGrid class.
    bool hasBoundaryIntersections() const;

    // Mimic Dune entity wrapper
    /// @brief Access the actual implementation class behind Entity interface class.
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

    /// \brief Returns (1) parent entity in the level-grid the parent cell was born, if the entity was born in any LGR.
    ///                (2) if the entity has no father, and is a leaf-grid-view entity, it returns the equivalent element on the level
    ///                    that it was born. Namely, for coarse cells on the leaf never involved in any refinement process, we get the
    ///                    equivalent entity in the GLOBAL grid (level 0). Notice that if it's a refined cell on the leaf, it does have a
    ///                    father, in that case, this method returns the father entity.
    ///                (3) Otherwise, returns itself. Notice that if the element is a refined one, this method returns the father() entity.
    ///                    For a coarse cell never involved in any refinement, if it's also a leaf-grid-cell, it returns the equivalent
    ///                    cell from level 0. Otherwise, it returns itself, which means that its grid is level 0 and the entity itself is already
    ///                    its "origin".
    Entity<0> getOrigin() const;

    /// \brief To be invoked only for leaf-grid-view entities. Get equivalent element on the level grid the leaf-entity was born.
    Entity<0> getLevelElem() const;

    /// \brief Get equivalent element on the level grid where the entity was born, if grid = leaf-grid-view. Otherwise, return itself.
    Entity<0> getEquivLevelElem() const;

    /// \brief Get Cartesian Index in the level grid view where the Entity was born.
    int getLevelCartesianIdx() const;

protected:
    const CpGridData* pgrid_;
};

} // namespace cpgrid
} // namespace Dune

// now we include the Iterators.hh We need to do this here because for hbegin/hend the compiler
// needs to know the size of hierarchicIterator
#include "Iterators.hpp"
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
HierarchicIterator Entity<codim>::hbegin(int maxLevel) const
{
    // Creates iterator with first child as target if there is one. Otherwise empty stack and target.
    return HierarchicIterator(*this, maxLevel);
}

/// Dummy beyond last child iterator.
template<int codim>
HierarchicIterator Entity<codim>::hend(int maxLevel) const
{
    // Creates iterator with empty stack and target.
    return HierarchicIterator(maxLevel);
}

template <int codim>
PartitionType Entity<codim>::partitionType() const
{
    return pgrid_->partition_type_indicator_->getPartitionType(*this);
}

} // namespace cpgrid
} // namespace Dune


#include <opm/grid/cpgrid/CpGridData.hpp>

namespace Dune {
namespace cpgrid {

template<int codim>
grid_size_t Entity<codim>::subEntities ( const grid_size_t cc ) const
{
    if (cc == codim) {
        return 1;
    }
    else if ( codim == 0 ){ // Cell/element/Entity<0>
        if ( cc == 3 ) { // Get number of corners of the element.
            return 8;
        }
    }
    return 0;
}

template <int codim>
const typename Entity<codim>::Geometry& Entity<codim>::geometry() const
{
    return pgrid_->geomVector<codim>()[*this];
}

template <int codim>
template <int cc>
typename Entity<codim>::template Codim<cc>::Entity Entity<codim>::subEntity(int i) const
{
    static_assert(codim == 0, "");
    if (cc == 0) { // Cell/element/Entity<0>
        assert(i == 0);
        typename Codim<cc>::Entity se(*pgrid_, this->index(), this->orientation());
        return se;
    } else if (cc == 3) { // Corner/Entity<3>
        assert(i >= 0 && i < 8);
        int corner_index = pgrid_->cell_to_point_[this->index()][i];
        typename Codim<cc>::Entity se(*pgrid_, corner_index, true);
        return se;
    }
    else {
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
    return pgrid_ ?  EntityRep<codim>::index() < pgrid_->size(codim) : false;
}


// level() It simply returns the level of the entity in the grid hierarchy.
template <int codim>
int Entity<codim>::level() const
{
    // Parallel and LGRs:
    // If the grid has been distributed and - after that - LGRs have been added, then level_data_ptr_
    // points at distributed_data_ which has size > 1.
    //
    // Serial and LGRs:
    // If the grid is not distributed and there LGRs have been added, level_data_ptr_ points at data_ which
    // has size > 1.
    //
    // - leaf_to_level_cells_ is non-empty only on the leaf grid view of a grid that has been refined.
    // - level_ is set equal to zero when instantiating a pgrid, and rewriten when it corresponds to a refined level grid.
    return pgrid_->leaf_to_level_cells_.empty()? pgrid_->level_ : pgrid_->leaf_to_level_cells_[this-> index()][0];
}

// isLeaf()
// - if distributed_data_ is empty: an element is a leaf <-> hbegin and hend return the same iterator. Then,
//    *cells from level 0 (coarse grid) that are not parents, are Leaf.
//    *cells from any LGR are Leaf, since they do not children (nested refinement not supported).
// - if distrubuted_data_ is NOT empty: there may be children on a different process. Therefore,
// isLeaf() returns true <-> the element is a leaf entity of the global refinement hierarchy. Equivalently,
// it can be checked whether parent_to_children_cells_ is empty.

template<int codim>
bool Entity<codim>::isLeaf() const
{
    if (pgrid_ -> parent_to_children_cells_.empty()){ // LGR cells
        return true;
    }
    else {
        return (std::get<0>((pgrid_ -> parent_to_children_cells_)[this-> index()]) == -1);  // Cells from GLOBAL, not involved in any LGR
    }
}

template<int codim>
bool Entity<codim>::isNew() const
{
    // WIP
    // For an Entity "to be new" means that
    // - it has a father
    // - it has been created in the last call of adapt()
    // To determine that, we track the level of the Entity and
    // check if it was marked for refinement in the previos state
    // of current_view_data_
    return (isLeaf() && hasFather());
}


template<int codim>
bool Entity<codim>::mightVanish() const
{
    const auto refinementMark = pgrid_ -> getMark(*this);
    return (refinementMark == 1);
}

template<int codim>
bool Entity<codim>::hasFather() const
{
    if ((pgrid_ -> child_to_parent_cells_.empty()) || (pgrid_ -> child_to_parent_cells_[this->index()][0] == -1)){
        return false;
    }
    else{
        return true;
    }
}

template<int codim>
Entity<0> Entity<codim>::father() const
{
    if (this->hasFather()){
        const int& coarser_level = pgrid_ -> child_to_parent_cells_[this->index()][0];
        const int& parent_cell_index = pgrid_ -> child_to_parent_cells_[this->index()][1];
        return Entity<0>( *((*(pgrid_ -> level_data_ptr_))[coarser_level].get()), parent_cell_index, true);
    }
    else{
        OPM_THROW(std::logic_error, "Entity has no father.");
    }
}

template<int codim>
Dune::cpgrid::Geometry<3,3> Dune::cpgrid::Entity<codim>::geometryInFather() const
{
    if (!(this->hasFather())){
        OPM_THROW(std::logic_error, "Entity has no father.");
    }

    // Indices of corners in entity's geometry in father reference element.
    static constexpr std::array<int,8> in_father_reference_elem_corner_indices = {0,1,2,3,4,5,6,7};
    // 'static': The returned object Geometry<3,3> stores a pointer to in_father_reference_elem_corner_indices. Therefore,
    // this variable is declared static to prolongate its lifetime beyond this function (static storage duration).

    auto idx_in_parent_cell = pgrid_ -> cell_to_idxInParentCell_[this->index()];
    if (idx_in_parent_cell !=-1) {
        const auto& cells_per_dim =  (*(pgrid_ -> level_data_ptr_))[this->level()] -> cells_per_dim_;
        const auto& auxArr = pgrid_ -> getReferenceRefinedCorners(idx_in_parent_cell, cells_per_dim);
        FieldVector<double, 3> corners_in_father_reference_elem_temp[8] =
            { auxArr[0], auxArr[1], auxArr[2], auxArr[3], auxArr[4], auxArr[5], auxArr[6], auxArr[7]};
        auto in_father_reference_elem_corners = std::make_shared<EntityVariable<cpgrid::Geometry<0, 3>, 3>>();
        EntityVariableBase<cpgrid::Geometry<0, 3>>& mutable_in_father_reference_elem_corners = *in_father_reference_elem_corners;
        // Assign the corners. Make use of the fact that pointers behave like iterators.
        mutable_in_father_reference_elem_corners.assign(corners_in_father_reference_elem_temp,
                                                        corners_in_father_reference_elem_temp + 8);
        // Compute the center of the 'local-entity'.
        Dune::FieldVector<double, 3> center_in_father_reference_elem = {0., 0.,0.};
        for (int corn = 0; corn < 8; ++corn) {
            for (int c = 0; c < 3; ++c)
            {
                center_in_father_reference_elem[c] += corners_in_father_reference_elem_temp[corn][c]/8.;
            }
        }
        // Compute the volume of the 'local-entity'.
        double volume_in_father_reference_elem = double(1)/(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
        // Construct (and return) the Geometry<3,3> of 'child-cell in the reference element of its father (unit cube)'.
        return Dune::cpgrid::Geometry<3,3>(center_in_father_reference_elem, volume_in_father_reference_elem,
                                           in_father_reference_elem_corners, in_father_reference_elem_corner_indices.data());
    }
    else {
        OPM_THROW(std::logic_error, "Entity has no father.");
    }
}

template<int codim>
Dune::cpgrid::Entity<0> Dune::cpgrid::Entity<codim>::getOrigin() const
{
    if (hasFather())
    {
        return this->father();
    }
    if (!(pgrid_ -> leaf_to_level_cells_.empty())) // entity on the LeafGridView
    {  // leaf_to_level_cells_ [leaf idx] = { level where entity was born, cell idx in that level}
        const int& levelElem = pgrid_->leaf_to_level_cells_[this->index()][0];
        const int& levelElemIdx = pgrid_->leaf_to_level_cells_[this->index()][1];
        return Dune::cpgrid::Entity<0>( *((*(pgrid_ -> level_data_ptr_))[levelElem].get()), levelElemIdx, true);
    }
    else
    {
        return *this;
    }
}

template<int codim>
Dune::cpgrid::Entity<0> Dune::cpgrid::Entity<codim>::getLevelElem() const
{
    // Check that the element belongs to the leaf grid view
    // This is needed to get the index of the element in the level it was born.
    // leaf_to_level_cells_ [leaf idx] = {level where the entity was born, equivalent cell idx in that level}
    if (!(pgrid_ -> leaf_to_level_cells_.empty())) // entity on the LeafGridView
    {
        const int& entityLevelIdx = pgrid_->leaf_to_level_cells_[this->index()][1];
        return Dune::cpgrid::Entity<0>( *((*(pgrid_ -> level_data_ptr_))[this->level()].get()), entityLevelIdx, true);
    }
    else {
        throw std::invalid_argument("The entity provided does not belong to the leaf grid view. ");
    }
}

template<int codim>
Dune::cpgrid::Entity<0> Dune::cpgrid::Entity<codim>::getEquivLevelElem() const
{
    // Check if the element belongs to the leaf grid view
    // This is needed to get the index of the element in the level it was born.
    // leaf_to_level_cells_ [leaf idx] = {level where the entity was born, equivalent cell idx in that level}
    if (!(pgrid_ -> leaf_to_level_cells_.empty())) // entity on the LeafGridView
    {
        const int& entityLevelIdx = pgrid_->leaf_to_level_cells_[this->index()][1];
        return Dune::cpgrid::Entity<0>( *((*(pgrid_ -> level_data_ptr_))[this->level()].get()), entityLevelIdx, true);
    }
    else {
        return *this;
    }
}

template<int codim>
int Dune::cpgrid::Entity<codim>::getLevelCartesianIdx() const
{
    const auto& level_data = (*(pgrid_ -> level_data_ptr_))[level()].get();
    // getLevelElem() throws when the entity does not belong to the leaf grid view.
    return level_data -> global_cell_[getLevelElem().index()];
}

} // namespace cpgrid
} // namespace Dune


#endif // OPM_ENTITY_HEADER
