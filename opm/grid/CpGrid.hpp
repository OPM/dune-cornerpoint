//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri May 29 20:26:36 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2014, 2022-2023 Equinor ASA.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulartion-Software & Services
  Copyright 2015       NTNU

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

#ifndef OPM_CPGRID_HEADER
#define OPM_CPGRID_HEADER

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>  // Not really needed it seems, but alas.

#include <dune/grid/common/grid.hh>
#include <opm/grid/cpgrid/CpGridDataTraits.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/OrientedEntityTable.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h> //  Not really needed it seems, but alas.
#include "common/GridEnums.hpp"
#include <opm/grid/utility/OpmWellType.hpp>

#include <set>

namespace Opm
{
struct NNCdata;
class EclipseGrid;
class EclipseState;
}

namespace Dune
{
    class CpGrid;

    namespace cpgrid
    {
    class CpGridData;
    template <int> class Entity;
    template<int,int> class Geometry;
    class HierarchicIterator;
    class IntersectionIterator;
    template<int, PartitionIteratorType> class Iterator;
    class LevelGlobalIdSet;
    class GlobalIdSet;
    class Intersection;
    class IntersectionIterator;
    class IndexSet;
    class IdSet;

    }
}

namespace Dune
{

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridTraits
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridTraits
    {
        /// \brief The type that implements the grid.
        typedef CpGrid Grid;

        /// \brief The type of the intersection at the leafs of the grid.
        typedef cpgrid::Intersection LeafIntersection;
        /// \brief The type of the intersection at the levels of the grid.
        typedef cpgrid::Intersection LevelIntersection;
        /// \brief The type of the intersection iterator at the leafs of the grid.
        typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
        /// \brief The type of the intersection iterator at the levels of the grid.
        typedef cpgrid::IntersectionIterator LevelIntersectionIterator;

        /// \brief The type of the  hierarchic iterator.
        typedef cpgrid::HierarchicIterator HierarchicIterator;

        /// \brief Traits associated with a specific codim.
        /// \tparam cd The codimension.
        template <int cd>
        struct Codim
        {
            /// \brief The type of the geometry associated with the entity.
            /// IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
             typedef cpgrid::Geometry<3-cd, 3> Geometry;
            //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> Geometry;
            /// \brief The type of the local geometry associated with the entity.
             typedef cpgrid::Geometry<3-cd, 3> LocalGeometry;
            //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> LocalGeometry;
            /// \brief The type of the entity.
            typedef cpgrid::Entity<cd> Entity;

            /// \brief The type of the iterator over all level entities of this codim.
            typedef cpgrid::Iterator<cd, All_Partition> LevelIterator;

            /// \brief The type of the iterator over all leaf entities of this codim.
            typedef cpgrid::Iterator<cd, All_Partition> LeafIterator;

            /// \brief The type of the entity pointer for entities of this codim.
            typedef cpgrid::Entity<cd> EntitySeed;

            /// \brief Traits associated with a specific grid partition type.
            /// \tparam pitype The type of the grid partition.
            template <PartitionIteratorType pitype>
            struct Partition
            {
                /// \brief The type of the iterator over the level entities of this codim on this partition.
                typedef cpgrid::Iterator<cd, pitype> LevelIterator;
                /// \brief The type of the iterator over the leaf entities of this codim on this partition.
                typedef cpgrid::Iterator<cd, pitype> LeafIterator;
            };
        };

        /// \brief Traits associated with a specific grid partition type.
        /// \tparam pitype The type of the grid partition.
        template <PartitionIteratorType pitype>
        struct Partition
        {
            /// \brief The type of the level grid view associated with this partition type.
            typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid> > LevelGridView;
            /// \brief The type of the leaf grid view associated with this partition type.
            typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid> > LeafGridView;

        };

        /// \brief The type of the level grid view associated with this partition type.
        typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid>> LevelGridView;
        /// \brief The type of the leaf grid view associated with this partition type.
        typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid>> LeafGridView;

        /// \brief The type of the level index set.
        typedef cpgrid::IndexSet LevelIndexSet;
        /// \brief The type of the leaf index set.
        typedef cpgrid::IndexSet LeafIndexSet;
        /// \brief The type of the global id set.
        typedef cpgrid::GlobalIdSet GlobalIdSet;
        /// \brief The type of the local id set.
        typedef GlobalIdSet LocalIdSet;

        /// \brief The type of the collective communication.
        using Communication = cpgrid::CpGridDataTraits::Communication;
        using CollectiveCommunication = cpgrid::CpGridDataTraits::CollectiveCommunication;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridFamily
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridFamily
    {
        typedef CpGridTraits Traits;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGrid
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> provides \ref Dune::Grid </em>]
    class CpGrid
        : public GridDefaultImplementation<3, 3, double, CpGridFamily>
    {
        friend class cpgrid::CpGridData;
        friend class cpgrid::Entity<0>;
        friend class cpgrid::Entity<1>;
        friend class cpgrid::Entity<2>;
        friend class cpgrid::Entity<3>;
        template<int dim>
        friend cpgrid::Entity<dim> createEntity(const CpGrid&,int,bool);

    public:

        // --- Typedefs ---


        /// Family typedef, why is this not defined by Grid<>?
        typedef CpGridFamily GridFamily;


        // --- Methods ---


        /// Default constructor
        CpGrid();

        CpGrid(MPIHelper::MPICommunicator comm);

#if HAVE_ECL_INPUT
        /// Read the Eclipse grid format ('grdecl').
        ///
        /// \return Vector of global indices to the cells which have
        ///         been removed in the grid processing due to small pore volume. Function only returns
        ///         indices on rank 0, the vector is empty of other ranks.
        /// \param ecl_grid the high-level object from opm-parser which represents the simulation's grid
        ///        In a parallel run this may be a nullptr on all ranks but rank zero.
        /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
        ///        aquifer information when ecl_state is available. NNC and aquifer connection
        ///        information will also be updated during the function call when available and necessary.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        /// \param pinchActive Force specific pinch behaviour. If true a face will connect two vertical cells, that are
        ///           topological connected, even if there are cells with zero volume between them. If false these
        ///           cells will not be connected despite their faces coinciding.
        std::vector<std::size_t>
        processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                             Opm::EclipseState* ecl_state,
                             bool periodic_extension, bool turn_normals, bool clip_z,
                             bool pinchActive);

        /// Read the Eclipse grid format ('grdecl').
        ///
        /// Pinch behaviour is determind from the parameter ecl_grid. If ecl_grid is a nullptr or PINCH was specified for
        /// the grid, then a face will connect two vertical cells, that are topological connected, even if there are
        /// cells with zero volume between them, Otherwise these cells will not be connected despite their faces coinciding.
        ///
        /// \return Vector of global indices to the cells which have
        ///         been removed in the grid processing due to small pore volume. Function only returns
        ///         indices on rank 0, the vector is empty of other ranks.
        /// \param ecl_grid the high-level object from opm-parser which represents the simulation's grid
        ///        In a parallel run this may be a nullptr on all ranks but rank zero.
        /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
        ///        aquifer information when ecl_state is available. NNC and aquifer connection
        ///        information will also be updated during the function call when available and necessary.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        std::vector<std::size_t>
        processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                             Opm::EclipseState* ecl_state,
                             bool periodic_extension, bool turn_normals = false, bool clip_z = false);

#endif

        /// Read the Eclipse grid format ('grdecl').
        /// \param input_data the data in grdecl format, declared in preprocess.h.
        /// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
        void processEclipseFormat(const grdecl& input_data, bool remove_ij_boundary, bool turn_normals = false);

        //@}

        /// \name Cartesian grid extensions.
        ///
        /// A cornerpoint grid can be seen as a degenerated and distorted cartesian
        /// grid. Therefore it provides mappings from cells to the underlying cartesian
        /// index.
        //@{
        /// Create a cartesian grid.
        /// \param dims the number of cells in each cartesian direction.
        /// \param cellsize the size of each cell in each dimension.
        /// \param shift The origin of the grid, i.e. the corner of the cell with index (0,0,0) where
        ///              the left, bottom, and top face of that cell intersect, is at the coordinate
        ///              origin per default. This parameter shifts that corner to lie at
        ///              (shift[0]*cellsize[0], ..., shift[2]*cellsize[2]).
        void createCartesian(const std::array<int, 3>& dims,
                             const std::array<double, 3>& cellsize,
                             const std::array<int, 3>& shift = {0,0,0});

        /// The logical cartesian size of the global grid.
        /// This function is not part of the Dune grid interface,
        /// and should be used with caution.
        const std::array<int, 3>& logicalCartesianSize() const;

        /// Retrieve mapping from internal ("compressed") active grid
        /// cells to external ("uncompressed") cells.  Specifically,
        /// @code globalCell()[i] @endcode is the linearized Cartesian
        /// index of grid cell @code i @endcode.  This method should
        /// only be used by classes which really need it, such as
        /// those dealing with permeability fields from the input deck
        /// from whence the current CpGrid was constructed.
        const std::vector<int>& globalCell() const;

        /// @brief Returns either data_ or distributed_data_(if non empty).
        const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& currentData() const;

        /// @brief Returns either data_ or distributed_data_(if non empty).
        std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& currentData();

        /// @brief
        ///    Extract Cartesian index triplet (i,j,k) of an active cell.
        ///
        /// @param [in] c
        ///    Active cell index.
        ///
        /// @param [out] ijk  Cartesian index triplet
        void getIJK(const int c, std::array<int,3>& ijk) const;
        //@}

        /// Is the grid currently using unique boundary ids?
        /// \return true if each boundary intersection has a unique id
        ///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
        bool uniqueBoundaryIds() const;

        /// Set whether we want to have unique boundary ids.
        /// \param uids if true, each boundary intersection will have a unique boundary id.
        void setUniqueBoundaryIds(bool uids);
       

        // --- Dune interface below ---

        /// \name The DUNE grid interface implementation
        // \@{
        /// \brief Get the grid name.
        ///
        /// It's the same as the class name.
        /// What did you expect, something funny?
        std::string name() const;
        
        /// Return maximum level defined in this grid. Levels are 0 and 1,  maxlevel = 1 (not counting leafview), 0 = the coarsest level.
        int maxLevel() const;

        /// Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const;
        /// one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const;
        
        /// Iterator to first entity of given codim on level and PartitionIteratorType
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const;
        /// one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const;

        /// Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const;
        /// one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const;

        /// Iterator to first leaf entity of given codim and PartitionIteratorType
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const;
        /// one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const;

        /// \brief Number of grid entities per level and codim
        int size (int level, int codim) const;

        /// number of leaf entities per codim in this process
        int size (int codim) const;

        /// number of entities per level and geometry type in this process
        int size (int level, GeometryType type) const;

        /// number of leaf entities per geometry type in this process
        int size (GeometryType type) const;

        /// \brief Access to the GlobalIdSet
        const Traits::GlobalIdSet& globalIdSet() const;

        /// \brief Access to the LocalIdSet
        const Traits::LocalIdSet& localIdSet() const;
        
        /// \brief Access to the LevelIndexSets
        const Traits::LevelIndexSet& levelIndexSet(int level) const;

        /// \brief Access to the LeafIndexSet
        const Traits::LeafIndexSet& leafIndexSet() const;

        /// \brief  Refine the grid refCount times using the default refinement rule.
        ///         This behaves like marking all elements for refinement and then calling preAdapt, adapt and postAdapt.
        ///         The state after globalRefine is comparable to the state after postAdapt.
        void globalRefine (int refCount);

        const std::vector<Dune::GeometryType>& geomTypes(const int) const;

        /// given an EntitySeed (or EntityPointer) return an entity object
        template <int codim>
        cpgrid::Entity<codim> entity(const cpgrid::Entity<codim>& seed) const;

        /// @brief Create a grid out of a coarse one and (at most) 2 refinements(LGRs) of selected block-shaped disjoint patches
        ///        of cells from that coarse grid.
        ///
        /// Level0 refers to the coarse grid, assumed to be this-> data_[0]. Level1 and level2 refer to the LGRs (stored in this->data_[1]
        /// data_[2]). LeafView (stored in this-> data_[3]) is built with the level0-entities which weren't involded in the
        /// refinenment, together with the new born entities created in level1 and level2.
        /// Old-corners and old-faces (from coarse grid) lying on the boundary of the patches, get replaced by new-born-equivalent corners
        /// and new-born-faces.
        ///
        /// @param [in] cells_per_dim_vec      Vector of Number of (refined) cells in each direction that each
        ///                                    parent cell should be refined to.
        /// @param [in] startIJK_vec           Vector of Cartesian triplet indices where each patch starts.
        /// @param [in] endIJK_vec             Vector of Cartesian triplet indices where each patch ends.
        ///                                    Last cell part of each patch(lgr) will be
        ///                                    {endIJK_vec[<patch-number>][0]-1, ..., endIJK_vec[<patch-number>][2]-1}.
        /// @param [in] lgr_name_vec           Names (std::string) for the LGRs/levels.
        void addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec);

        // @brief TO BE DONE
        const std::map<std::string,int>& getLgrNameToLevel() const;

        // @breif Compute center of an entity/element/cell in the Eclipse way:
        //        - Average of the 4 corners of the bottom face.
        //        - Average of the 4 corners of the top face.
        //        Return average of the previous computations.
        // @param [in]   int   Index of a cell.
        // @return            'eclipse centroid'
        std::array<double,3> getEclCentroid(const int& idx) const;

        // @breif Compute center of an entity/element/cell in the Eclipse way:
        //        - Average of the 4 corners of the bottom face.
        //        - Average of the 4 corners of the top face.
        //        Return average of the previous computations.
        // @param [in]   Entity<0>   Entity
        // @return                   'eclipse centroid'
        std::array<double,3> getEclCentroid(const cpgrid::Entity<0>& elem) const;

        // @brief Return parent (coarse) intersection (face) of a refined face on the leaf grid view, whose neighboring cells
        //        are two: one coarse cell (equivalent to its origin cell from level 0), and one refined cell
        //        from certain LGR
        Dune::cpgrid::Intersection getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const;

        /// --------------- Adaptivity (begin) ---------------
        /// @brief Mark entity for refinement (or coarsening).
        ///
        /// Refinement on CpGrid is partially supported for Cartesian grids, with the keyword CARFIN.
        /// Nested refinement is not supported yet, so the so-called "host grid" is defined by default
        /// equal to the GLOBAL Grid (level zero). Therefore, we mark elements (from the GLOBAL grid)
        /// for refinement.
        /// This only works for entities of codim 0.
        ///
        /// @param [in] refCount   To mark the element for
        ///                        - refinement, refCount == 1
        ///                        - doing nothing, refCount == 0
        ///                        - coarsening, refCount == -1 (not applicable yet)
        /// @param [in] element    Entity<0>. Currently, an element from the GLOBAL grid (level zero).
        /// @return true, if marking was succesfull.
        ///         false, if marking was not possible.
        bool mark(int refCount, const cpgrid::Entity<0>& element);

        /// @brief Return refinement mark for entity.
        ///
        /// @return refinement mark (1,0,-1)  Currently, only 1 (refinement), or 0 (doing nothing).
        int getMark(const cpgrid::Entity<0>& element) const;

        /// @brief Set mightVanish flags for elements that will be refined in the next adapt() call
        ///        Need to be called after elements have been marked for refinement.
        bool preAdapt();

        /// @brief Triggers the grid refinement process
        bool adapt();

        /// @brief Triggers the grid refinement process, allowing to select diffrent refined level grids.
        ///
        /// @param [in] cells_per_dim_vec    For each set of marked elements for refinement, that will belong to a same
        ///                                  refined level grid, number of (refined) cells in each direction that each
        ///                                  parent cell should be refined to.
        /// @param [in] assignRefinedLevel   Vector with size equal to total amount of cells of the starting grid where
        ///                                  the marked elements belong. In each entry, the refined level grid where the
        ///                                  refined entities of the (parent) marked element should belong is stored.
        /// @param [in] lgr_name_vector      Each refined level grid name, e.g. {"LGR1", "LGR2"}.
        /// @param [in] isCARFIN             Default false. The keyword CARFIN implies that the selected cells to be refined form a block,
        ///                                  which can be edscribed via startIJK and endIJK Cartesian indices. This bool
        ///                                  is used to define logical_cartesian_size_ of the refined level grids according
        ///                                  to this block shape.
        /// @param [in] startIJK_vec         Default empty vector. When isCARFIN, the starting ijk Cartesian index of each
        ///                                  block of cells to be refined.
        /// @param [in] endIJK_vec           Default empty vector. When isCARFIN, the final ijk Cartesian index of each
        ///                                  block of cells to be refined.
        bool adapt(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                   const std::vector<int>& assignRefinedLevel,
                   const std::vector<std::string>& lgr_name_vec,
                   bool isCARFIN = false,
                   const std::vector<std::array<int,3>>& startIJK_vec = std::vector<std::array<int,3>>{},
                   const std::vector<std::array<int,3>>& endIJK_vec = std::vector<std::array<int,3>>{});

        /// @brief Clean up refinement markers - set every element to the mark 0 which represents 'doing nothing'
        void postAdapt();
        /// --------------- Adaptivity (end) ---------------

    private:

        /// --------------- Auxiliary methods to support Adaptivity (begin) ---------------
        
        /// @brief Refine each marked element and establish relationships between corners, faces, and cells marked for refinement,
        ///        with the refined corners, refined faces, and refined cells.
        ///
        /// --- Marked elements parameters ---
        /// @param [out] markedElem_to_itsLgr:                      Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
        ///                                                         build each refined level grid containing all the refined entities assigned for that grid.
        /// @param [out] markedElem_count:                          Total amount of marked elements to be refined. It will be used to print grid info.
        /// @param [out] cornerInMarkedElemWithEquivRefinedCorner:  For each corner from level zero, we store the marked elements where the corner appears and its equivalent
        ///                                                         refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
        ///                                                         elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
        ///                                                         cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
        ///                                                         For corners not appearing in any marked element, empty vector.
        /// @param [out] markedElemAndEquivRefinedCorner_to_corner: To correctly build the level-refined and adapted-grid topology features, we need to keep track of the
        ///                                                         corners that got replaced by equivalent refined corners, in each marked element where the corner appeared,
        ///                                                         not only in its last appearance. The last appearance will be used to avoid repetition when storing.
        ///                                                         Following the example above,
        ///                                                         markedElemAndEquivRefinedCorner_to_corner[{0, 8}] = 5;
        ///                                                         markedElemAndEquivRefinedCorner_to_corner[{1, 2}] = 5;
        /// @param [out] faceInMarkedElemAndRefinedFaces:           For each face from level zero, we store the marked elements where the face appears (maximum 2 cells)
        ///                                                         and its new-born refined faces from each auxiliary marked-element-lgr. Example: face with index 9
        ///                                                         appears in marked elements 0 and 1. Then,
        ///                                                         faceInMarkedElemAndRefinedFaces[9] = {{0, {refinedFace0_0, ..., refinedFaceN_0}},
        ///                                                                                               {1, {refinedFace0_1, ..., refinedFaceM_1}}}.
        ///                                                         For faces not appearing in any marked element, empty vector.
        /// --- Refined cells parameters ---
        /// @param [out] elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell:  Each marked element has been refined in its "own elemLgr". Refined entities should be stored in
        ///                                                                   the corresponding assigned refined level grid. To keep track of the cell index relation,
        ///                                                                   associate each
        ///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
        ///                                                                   { refined level grid assigned for the marked element, refined cell index in refined level grid }.
        /// @param [out] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
        ///                                                                   cell index relation, associate each
        ///                                                                   { refined level grid assigned for the marked element, refined cell index in refined level grid }
        ///                                                                   with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [out] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
        /// @param [in] assignRefinedLevel:                                   Each marked element can be assigned to certain refined level grid. This vector has entries 0 for
        ///                                                                   non marked elements, and the corresponding integer representing a refined level grid for marked
        ///                                                                   elements.
        /// @param [out] preAdapt_parent_to_children_cells_vec:               Parent cells and their refined children. Entry is {-1, {}} when cell has no children. Othewise,
        ///                                                                   {refined grid level where children were born, {child0, child1, ...}}
        ///                                                                   Each vector entry represents an existing level grid before calling adapt.
        /// --- Adapted cells parameters ---
        /// @param [out] elemLgrAndElemLgrCell_to_adaptedCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
        ///                                                                   the corresponding leaf grid view (or adapted grid). To keep track of the cell index relation,
        ///                                                                   associate each
        ///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
        ///                                                                   refined cell index inthe leaf grid view (or adapted grid).
        /// @param [out] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
        ///                                                                   the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
        ///                                                                   relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
        ///                                                                   { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [out] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
        /// @param [out] preAdapt_level_to_leaf_cells_vec:                    For each existing grid before calling adapt, we stablish the index relation between preAdapt cells
        ///                                                                   and cells on the leaf grid view (or adapted cells).-1 means that the cell vanished.
        /// --- Additional parameters ---
        /// @param [in] cells_per_dim_vec:                                    For each set of marked elements for refinement, that will belong to a same
        ///                                                                   refined level grid, number of (refined) cells in each direction that each
        ///                                                                   parent cell should be refined to.
        void refineAndProvideMarkedRefinedRelations(/* Marked elements parameters */
                                                    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                                    int& markedElem_count,
                                                    std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                                    std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                                    std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                                    /* Refined cells parameters */
                                                    std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell,
                                                    std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                    std::vector<int>& refined_cell_count_vec,
                                                    const std::vector<int>& assignRefinedLevel,
                                                    std::vector<std::vector<std::tuple<int,std::vector<int>>>>& preAdapt_parent_to_children_cells_vec,
                                                    /* Adapted cells parameters */
                                                    std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                                    std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                                    int& cell_count,
                                                    std::vector<std::vector<int>>& preAdapt_level_to_leaf_cells_vec,
                                                    /* Additional parameters */
                                                    const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief  Define child-parent relations from the new refined cells of the new refined level grids to its parent cells (belonging to pre-existing grid,
        ///         before adapting the grid/before updating the leaf grid view). Define the index in parent cell (-1 when cell has no parent).
        ///
        /// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
        ///                                                                  cell index relation, associate each
        ///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }
        ///                                                                  with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [in] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
        /// @param [in] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
        ///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
        ///                                                                  relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
        ///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [in] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
        ///
        /// @return refined_child_to_parent_cells_vec:   Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise,
        ///                                              {level parent cell, parent cell index}. Each vector entry represents a refined level grid.
        ///         refined_cell_to_idxInParentCell_vec: Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
        ///                                              Each vector entry represents a refined level grid.
        ///         adapted_child_to_parent_cell:        Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise,
        ///                                              {level parent cell, parent cell index}
        ///         adapted_cell_to_idxInParentCell:     Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather(). -1 when has no father.
        std::tuple< std::vector<std::vector<std::array<int,2>>>,
                    std::vector<std::vector<int>>,
                    std::vector<std::array<int,2>>,
                    std::vector<int>> defineChildToParentAndIdxInParentCell( const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                             const std::vector<int>& refined_cell_count_vec,
                                                                             const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                                                             const int& cell_count) const;

        /// @brief Define refined level grid cells indices and leaf grid view (or adapted grid) cells indices relations. Namely, level_to_leaf_cells_ for each new
        ///        refined level grid, and leaf_to_level_cells_ for the updated leaf grid view. 
        ///
        /// @param [in] elemLgrAndElemLgrCell_to_refinedLevelAdRefinedCell:  Each marked element has been refined in its "own elemLgr". Refined entities should be stored in
        ///                                                                  the corresponding assigned refined level grid. To keep track of the cell index relation, we 
        ///                                                                  associate each
        ///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
        ///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }.
        /// @param [in] refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell: Each marked element has been assigned to certain refined level grid. To keep track of the "inverse"
        ///                                                                  cell index relation, associate each
        ///                                                                  { refined level grid assigned for the marked element, refined cell index in refined level grid }
        ///                                                                  with { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [in] refined_cell_count_vec:                              Total amount of refined cells, per level (i.e. in each refined level grid).
        /// @param [in] elemLgrAndElemLgrCell_to_adaptedCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
        ///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the cell index relation,
        ///                                                                  associate each
        ///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement } with
        ///                                                                  refined cell index inthe leaf grid view (or adapted grid).
        /// @param [in] adaptedCell_to_elemLgrAndElemLgrCell:                Each marked element has been refined in its "own elemLgr". Refined entities should be also stored in
        ///                                                                  the corresponding leaf grid view (or adapted grid). To keep track of the "inverse" cell index
        ///                                                                  relation, associate the refined cell index inthe leaf grid view (or adapted grid) with
        ///                                                                  { marked element index ("elemLgr"), refined cell index in the auxiliary single-cell-refinement }.
        /// @param [in] cell_count:                                          Total amount of cells on the leaf grid view (or adapted grid).
        ///
        /// @return refined_level_to_leaf_cells_vec:                         refined_level_to_leaf_cells_vec[ levelGridIdx ] [ cell idx in that level grid ] = equivalent leaf cell idx
        ///         leaf_to_level_cells:                                     leaf_to_level_cells[ leaf cell idx ] = {level where cell was born, cell idx on that level}
        std::pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>
        defineLevelToLeafAndLeafToLevelCells(const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                             const std::vector<int>& refined_cell_count_vec,
                                             const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                             const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                             const int& cell_count) const;

        /// @brief Define various corner relations. 1. refined corners from auxiliary single marked element refinement to its corresponding refined level grid, and vice versa.
        ///                                         2. refined corners from single-cell-refinements that vanish in the "storing only once each entity process". To avoid repetition,
        ///                                            we store such corners in their "last apperance". We keep track of all the appearances since that is needed for correctly
        ///                                            define CpGridData attributes such as cell_to_point_ and face_to_point_.
        ///
        /// @param [out] elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
        ///                                                                       the corresponding assigned refined level grid. To keep track of the corner index relation, we
        ///                                                                       associate each
        ///                                                                       { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement } with
        ///                                                                       { refined level grid assigned for the marked element, refined corner index in refined level grid }.
        /// @param [out] refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
        ///                                                                       the corresponding assigned refined level grid. To keep track of the corner index relation, we
        ///                                                                       associate each
        ///                                                                       { refined level grid assigned for the marked element, refined corner index in refined level grid } with
        ///                                                                       { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }.
        /// @param [out] refined_corner_count_vec:                                Total amount of refined corners, per level (each vector entry corresponds to a refined level grid).
        /// @param [out] vanishedRefinedCorner_to_itsLastAppearance:              A refined corner might appear in several single-cell-refinements, we store it only in its last
        ///                                                                       appearance, but keep track of the vanishing. Example, a corner appears in total 3
        ///                                                                       single-cell-refinements, with indices { elemLgr1, elemLgr1Corner }, { elemLgr2, elemLgr2Corner },
        ///                                                                       and { elemLgr3, elemLgr3Corner }. Then, for X = 1, and X=2, we store
        ///                                                                       vanishedRefinedCorner_to_itsLastAppearance[{elemLgrX, elemLgrXCorner}] = {elemLgr3, elemLgr3Corner}.
        /// @param [in] markedElem_to_itsLgr
        /// @param [in] assignRefinedLevel
        /// @param [in] cornerInMarkedElemWithEquivRefinedCorner
        /// @param [in] faceInMarkedElemAndRefinedFaces
        /// @param [in] cells_per_dim_vec
        void identifyRefinedCornersPerLevel(std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                            std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                            std::vector<int>& refined_corner_count_vec,
                                            std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                            const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                            const std::vector<int>& assignRefinedLevel,
                                            const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                            const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                            const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief Define relations between single-cell-refinement faces and refined level faces.
        ///
        /// @param [out] elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace: Each marked element has been refined in its "own elemLgr". Refined faces should be stored in
        ///                                                                   the corresponding assigned refined level grid. To keep track of the face index relation, we
        ///                                                                   associate each
        ///                                                                   { marked element index ("elemLgr"), refined face index in the auxiliary single-cell-refinement } with
        ///                                                                   { refined level grid assigned for the marked element, refined face index in refined level grid }.
        /// @param [out] refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace: Each marked element has been refined in its "own elemLgr". Refined faces should be stored in
        ///                                                                   the corresponding assigned refined level grid. To keep track of the face index relation, we
        ///                                                                   associate each
        ///                                                                   { refined level grid assigned for the marked element, refined face index in refined level grid } with
        ///                                                                   { marked element index ("elemLgr"), refined face index in the auxiliary single-cell-refinement }.
        /// @param [out] refined_face_count_vec:                              Total amount of refined corners, per level (each vector entry corresponds to a refined level grid).
        /// @param [in] markedElem_to_itsLgr
        /// @param [in] assignRefinedLevel
        /// @param [in] faceInMarkedElemAndRefinedFaces
        /// @param [in] cells_per_dim_vec
        void identifyRefinedFacesPerLevel(std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                          std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                          std::vector<int>& refined_face_count_vec,
                                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                          const std::vector<int>& assignRefinedLevel,
                                          const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                          const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief Identify corners that appear on the leaf grid view.
        ///        Define various corner relations. preAdapt or refined corners from auxiliary single marked element refinement to the leaf grid view (or adapted grid), and vice versa.
        ///
        /// @param [out] elemLgrAndElemLgrCorner_to_adaptedCorner: Each marked element has been refined in its "own elemLgr". When the element has not been refined, elemLgr == -1.
        ///                                                        To keep track of the corner index relation, we associate each
        ///                                                        { marked element index ("elemLgr"),  corner index in the auxiliary single-cell-refinement }, or
        ///                                                        { -1 ("elemLgr"),  corner index in the starting grid}, with
        ///                                                        corner index in the leaf grid view .
        /// @param [out] adaptedCorner_to_elemLgrAndElemLgrCorner: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
        ///                                                        the corresponding assigned refined level grid. To keep track of the cell index relation, we
        ///                                                        associate each corner index in the leaf grid view (or adapted grid) with
        ///                                                        { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }, or
        ///                                                        { -1 ("elemLgr"),  corner index in the starting grid},
        /// @param [out] corner_count:                             Total amount of corners on the leaf grid view (or adapted grid).
        /// @param [in] markedElem_to_itsLgr
        /// @param [in] assignRefinedLevel
        /// @param [in] cornerInMarkedElemWithEquivRefinedCorner
        /// @param [in] vanishedRefinedCorner_to_itsLastAppearance
        /// @param [in] faceInMarkedElemAndRefinedFaces
        /// @param [in] cells_per_dim_vec
        void identifyLeafGridCorners(std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                     std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                     int& corner_count,
                                     const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                     const std::vector<int>& assignRefinedLevel,
                                     const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                     std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                     const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                     const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief Identify faces that appear on the leaf grid view.
        ///        Define various face relations. preAdapt or refined faces from auxiliary single marked element refinement to the leaf grid view (or adapted grid), and vice versa.
        ///
        /// @param [out] elemLgrAndElemLgrFace_to_adaptedFace: Each marked element has been refined in its "own elemLgr". When the element has not been refined, elemLgr == -1.
        ///                                                    To keep track of the face index relation, we associate each
        ///                                                    { marked element index ("elemLgr"),  face index in the auxiliary single-cell-refinement }, or
        ///                                                    { -1 ("elemLgr"),  face index in the starting grid }, with
        ///                                                    face index in the leaf grid view .
        /// @param [out] adaptedFace_to_elemLgrAndElemLgrFace: Each marked element has been refined in its "own elemLgr". Refined corners should be stored in
        ///                                                    the corresponding assigned refined level grid. To keep track of the face index relation, we associate each
        ///                                                    face index in the leaf grid view (or adapted grid) with
        ///                                                    { marked element index ("elemLgr"), refined corner index in the auxiliary single-cell-refinement }. 
        /// @param [out] face_count:                           Total amount of faces on the leaf grid view (or adapted grid).
        /// @param [in] markedElem_to_itsLgr
        /// @param [in] assignRefinedLevel
        /// @param [in] faceInMarkedElemAndRefinedFaces
        /// @param [in] cells_per_dim_vec
        void identifyLeafGridFaces(std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                   std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                   int& face_count,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief Define the corners (gemotry) for each refined level grid.
        void populateRefinedCorners(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                    const std::vector<int>& refined_corner_count_vec,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    const int& preAdaptMaxLevel,
                                    const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner) const;

        /// @brief Define the faces, face tags, face normarls, and face_to_point_, for each refined level grid.
        void populateRefinedFaces(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refine_face_normals_vec,
                                  std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                  const std::vector<int>& refined_face_count_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                  const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const int& preAdaptMaxLevel,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner) const;

        /// @brief Define the cells, cell_to_point_, global_cell_, cell_to_face_, face_to_cell_, for each refined level grid.
        void populateRefinedCells(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                  std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                  std::vector<std::vector<int>>& refined_global_cell_vec,
                                  const std::vector<int>& refined_cell_count_vec,
                                  std::vector<cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                  std::vector<cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                  const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  const int& preAdaptMaxLevel,
                                  const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::vector<std::array<int,3>>&  cells_per_dim_vec) const;

        /// @brief Set geometrical and topological attributes for each refined level grid.
        void setRefinedLevelGridsGeometries( /* Refined corner arguments */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                             const std::vector<int>& refined_corner_count_vec,
                                             /* Refined face arguments */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                             std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                             std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refine_face_normals_vec,
                                             std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                             const std::vector<int>& refined_face_count_vec,
                                             /* Refined cell argumets */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                             std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                             std::vector<std::vector<int>>& refined_global_cell_vec,
                                             std::vector<int>& refined_cell_count_vec,
                                             std::vector<cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                             std::vector<cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                             /* Auxiliary arguments */
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                             const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                             const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                             const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                                             const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                             const std::vector<int>& assignRefinedLevel,
                                             const int& preAdaptMaxLevel,
                                             const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                             const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                             const std::vector<std::array<int,3>>&  cells_per_dim_vec) const;

        /// @brief Define the corners (gemotry) for the leaf grid view (or adapted grid).
        void populateLeafGridCorners(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners,
                                     const int& corners_count,
                                     const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                     const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner) const;

        /// @brief Define the faces, face tags, face normarls, and face_to_point_, for the leaf grid view.
        void populateLeafGridFaces(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces,
                                   Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                   Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                   Opm::SparseTable<int>& adapted_face_to_point,
                                   const int& face_count,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                   const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                   const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const int& preAdaptMaxLevel) const;

        /// @brief Define the cells, cell_to_point_, cell_to_face_, face_to_cell_, for the leaf grid view (or adapted grid).
        void populateLeafGridCells(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells,
                                   std::vector<std::array<int,8>>& adapted_cell_to_point,
                                   const int& cell_count,
                                   cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                   cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                   const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                   const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                   const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                   const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const int& preAdaptMaxLevel) const;

        /// @brief Define geometrical and topological attributes for the leaf grid view (or adapted grid).
        void updateLeafGridViewGeometries( /* Leaf grid View Corners arguments */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners,
                                           const int& corner_count,
                                           /* Leaf grid View Faces arguments */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces,
                                           Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                           Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                           Opm::SparseTable<int>& adapted_face_to_point,
                                           const int& face_count,
                                           /* Leaf grid View Cells argumemts  */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells,
                                           std::vector<std::array<int,8>>& adapted_cell_to_point,
                                           const int& cell_count,
                                           cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                           cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                           /* Auxiliary arguments */
                                           const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                           const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                           const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                           const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                           const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                           const std::vector<int>& assignRefinedLevel,
                                           const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                           const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                           const int& preAdaptMaxLevel) const;

        void updateCornerHistoryLevels(const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                       const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                       const int& corner_count,
                                       const std::vector<std::array<int,2>>& preAdaptGrid_corner_history,
                                       const int& preAdaptMaxLevel,
                                       const int& newLevels);


        /// @brief For refined level grids created based on startIJK and endIJK values, compute the "local ijk/Cartesian index" within the LGR.
        ///
        /// It's confusing that this "localIJK" is stored in CpGridData member global_cell_. Potential explanation: level zero grid is also called
        /// "GLOBAL" grid.
        /// Example: a level zero grid with dimension 4x3x3, an LGR with startIJK = {1,2,2}, endIJK = {3,3,3}, and cells_per_dim = {2,2,2}.
        /// Then the dimension of the LGR is (3-1)*2 x (3-2)*2 x (3-2)* 2 = 4x2x2 = 16. Therefore the global_cell_lgr minimim value should be 0,
        /// the maximum should be 15.
        /// To invoke this method, each refined level grid must have 1. logical_cartesian_size_, 2. cell_to_idxInParentCell_, and 3. cells_per_dim_
        /// already populated.
        ///
        /// @param [in] level    Grid index where LGR is stored
        /// @param [out] global_cell_lgr
        void computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr);

        /// @brief For a leaf grid with with LGRs, we assign the global_cell_ values of either the parent cell or the equivalent cell from
        ///        level zero.
        ///        For nested refinement, we lookup the oldest ancestor, from level zero.
        void computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf);

        /// @brief Get the ijk index of a refined corner, given its corner index of a single-cell-refinement.
        ///
        /// Given a single-cell, we refine it in {nx, ny, nz} refined children cells (per direction). Then, this single-cell-refinement
        /// has in total (nx +1)(ny +1)(nz +1) refined corners. Each of this corners has an ijk value associated since they are stored
        /// (following order defined in Geometry::refine) with the index (j*(nx+1)(nz+1)) + (i(nz+1)) + k, where i=0,...,nx, j=0,...,ny,
        /// and k=0,...,nz. This method returns the ijk, given the cornerIdxInLgr = 0,...,(nx +1)(ny +1)(nz +1).
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        std::array<int,3> getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const;

        /// @brief Get the ijk index of a reined face, given its corner index of a single-cell-refinement.
        ///
        /// Given a single-cell, we refine it in {nx, ny, nz} refined children cells (per direction). Then, this single-cell-refinement
        /// has in total ((nx+1)*ny*nz) + (nx*(ny+1)*nz) + (nx*ny*(nz+1)) refined faces. Each of this faces has an ijk value associated since
        /// they are stored  with the index (following order defined in Geometry::refine):
        /// K_FACES  (k*nx*ny) + (j*nx) + i
        /// I_FACES  (nx*ny*(nz+1)) + (i*ny*nz) + (k*ny) + j
        /// J_FACES   (nx*ny*(nz+1)) + ((nx+1)*ny*nz) + (j*nx*nz) + (i*nz) + k
        ///  where i=0,...,nx-1, j=0,...,ny-1, and k=0,...,nz-1. This method returns the corresponding ijk
        /// given a faceIdxInLgr = 0,...,((nx+1)nynz) + (nx*(ny+1)nz) + (nxny*(nz+1))
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
        /// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
        std::array<int,3> getRefinedFaceIJK(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                            const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const;

        /// @brief Determine if a refined corner is located in the interior of the single-cell-refinement.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        bool isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const;

        /// @brief Determine if a refined face is located in the interior of the single-cell-refinement.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
        /// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
        bool isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                        const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const;

        /// @brief Determine if a refined corner is located on the boundary of the single-cell-refinement, and does not coincide with
        ///        a preAdapt-existing corner.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        bool isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const;

        /// @brief Determine if a new refined corner is located on an edge of the parent cell. In particular, it's on the boundary of
        ///        the single-cell-refinement, and does not coincide with  a preAdapt-existing corner.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        bool newRefinedCornerLiesOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const;

        /// @brief Determine if a refined face is located on the boundary of the single-cell-refinement.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
        /// @param [in] elemLgr_ptr:      Pointer to the single-cell-refinement grid.
        bool isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                        const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const;

        /// @brief Get the parent face where the new refined corner lays on.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        /// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
        int getParentFaceWhereNewRefinedCornerLiesOn(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr, int elemLgr) const;

        /// @brief Get the parent faces that containes the edge where the new refined corner lies on.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] cornerIdxInLgr:   Corner index in the single-cell-refinement.
        /// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
        std::array<int,2> getParentFacesAssocWithNewRefinedCornLyingOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr, int elemLgr) const;

        /// @brief A refined corner appears in two single-cell-refinements. Given the corner index in the first single-cell-refinement, compute the
        ///         corner index in the neighboring single-cell-refinement.
        ///
        /// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
        /// @param [in] cornerIdxInLgr1:       Corner index in the elemLgr1 single-cell-refinement.
        /// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
        int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1, int cornerIdxLgr1, const std::array<int,3>& cells_per_dim_lgr2) const;

        /// @brief A new refined corner lays on an edge and appears in at least two single-cell-refinements. Given the corner index in one single-cell-refinement, compute the
        ///        corner index in a neighboring single-cell-refinement.
        ///
        /// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
        /// @param [in] cornerIdxInLgr1:       Corner index in the elemLgr1 single-cell-refinement.
        /// @param [in] parentFaceLastAppearanceIdx: Parent face index where the refined corner appears for last time.
        /// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
        int replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1, int cornerIdxLgr1, int elemLgr1, int parentFaceLastAppearanceIdx,
                                                const std::array<int,3>& cells_per_dim_lgr2) const;

        /// @brief A new refined face lays on the boudndary of a single-cell-refinement appears in at most two single-cell-refinements. Given the face index in one
        ///        single-cell-refinement, compute the face index in a neighboring single-cell-refinement.
        ///
        /// @param [in] cells_per_dim_lgr1:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr1 single-cell-refinement.
        /// @param [in] faceIdxInLgr1:         Face index in the elemLgr1 single-cell-refinement.
        /// @param [in] elemLgr1_ptr:          Pointer to the elemLgr1 single-cell-refinement grid.
        /// @param [in] cells_per_dim_lgr2:    Total children cells in each direction (x-,y-, and z-direction) of the elemLgr2 single-cell-refinement.
        int replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1, int faceIdxInLgr1,
                                            const std::shared_ptr<cpgrid::CpGridData>& elemLgr1_ptr,
                                            const std::array<int,3>& cells_per_dim_lgr2) const;

        /// @brief Get the parent face index where the new refined face lays on.
        ///
        /// @param [in] cells_per_dim:    Total children cells in each direction (x-,y-, and z-direction) of the single-cell-refinement.
        /// @param [in] faceIdxInLgr:     Face index in the single-cell-refinement.
        /// @param [in] elemLgr_ptr:      Pointer to the elemLgr single-cell-refinement grid.
        /// @param [in] elemLgr:          Cell index from starting grid, that has been refined into a single-cell-refinement.
        int getParentFaceWhereNewRefinedFaceLiesOn(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                                   const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr,
                                                   int elemLgr)  const;
        
        /// --------------- Auxiliary methods to support Adaptivity (end) ---------------

        /// @brief Check if there are non neighboring connections on blocks of cells selected for refinement.
        ///
        /// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
        /// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
        ///
        /// @return True if all blocks of cells do not contain any NNCs (non-neighboring-connection).
        ///         False if there is a block with at least one cell with a NNC.
        bool nonNNCsSelectedCellsLGR( const std::vector<std::array<int,3>>& startIJK_vec,
                                      const std::vector<std::array<int,3>>& endIJK_vec) const;

        /// @brief Detect active LGRs in each process.
        ///
        /// Given blocks of cells selected for refinement on a level zero distributed grid, detect which LGRs are active
        /// in each process.
        ///
        /// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
        /// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
        /// @param [out] lgr_with_at_least_one_active_cell Determine if an LGR is not empty in a given process, we set
        ///                                                lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
        ///                                                at least one active cell, and to 0 otherwise.
        void detectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                              const std::vector<std::array<int,3>>& endIJK_vec,
                              std::vector<int>& lgr_with_at_least_one_active_cell);

        /// @brief Mark selected elements, assign them their corresponding level, and detect active LGRs.
        ///
        /// Given blocks of cells selected for refinement, mark selected elements and assign them their corresponding
        /// (refined) level (grid). When level zero grid is distributed before refinement, detect which LGRs are active
        /// in each process.
        ///
        /// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
        /// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
        /// @param [out] assignRefinedLevel   Assign level for the refinement of each marked cell. Example: refined element from
        ///                                   LGR1 have level 1, refined element rfom LGR2 have level 2, etc.
        /// @param [out] lgr_with_at_least_one_active_cell Determine if an LGR is not empty in a given process, we set
        ///                                                lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
        ///                                                at least one active cell, and to 0 otherwise.
        void markElemAssignLevelDetectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                                 std::vector<int>& assignRefinedLevel,
                                                 std::vector<int>& lgr_with_at_least_one_active_cell);

        /// @brief Auxilliary function to compute one or more properties on selected block of parent cells.
        ///
        /// @param [in] startIJK_vec    Vector of ijk values denoting the start of each block of cells selected for refinement.
        /// @param [in] endIJK_vec      Vector of ijk values denoting the end of each block of cells selected for refinement.
        /// @param [in] function        Lambda expression/function that computes the desired properties for each parent cell.
        template<class T>
        void computeOnLgrParents(const std::vector<std::array<int,3>>& startIJK_vec,
                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                 T func);

        /// @brief Predict minimum cell and point global ids per process.
        ///
        /// Predict how many new cells/points (born in refined level grids) need new globalIds, so we can assign unique
        /// new ids ( and anticipate the maximum). At this point, the grid is already refined according to the LGR specification.
        ///
        /// @param [in] assignRefinedLevel   Assign level for the refinement of each marked cell. Example: refined element from
        ///                                  LGR1 have level 1, refined element rfom LGR2 have level 2, etc.
        /// @param [in] cells_per_dim_vec    Total child cells in each direction (x-,y-, and z-direction) per block of cells.
        /// @param [in] lgr_with_at_least_one_active_cell  Determine if an LGR is not empty in a given process:
        ///                                                lgr_with_at_least_one_active_cell[level] = 1 if it contains
        ///                                                at least one active cell in the current process, and 0 otherwise.
        /// @param [out] min_globalId_cell_in_proc
        /// @param [out] min_globalId_point_in_proc
        void predictMinCellAndPointGlobalIdPerProcess(const std::vector<int>& assignRefinedLevel,
                                                      const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                      const std::vector<int>& lgr_with_at_least_one_active_cell,
                                                      int& min_globalId_cell_in_proc,
                                                      int& min_globalId_point_in_proc) const;

        /// @brief Assign cell global ids of new born cell from refined level grids. Assign 'candidate' point global ids
        ///        for points in refined level grids.
        ///
        /// @param [out] localToGlobal_cells_per_level    Relation local element.index() to assigned cell global id.
        /// @param [out] localToGlobal_points_per_level   Relation local point.index() to assigned 'candidate' global id.
        /// @param [in] min_globalId_cell_in_proc         Minimum cell global id per process.
        /// @param [in] min_globalId_point_in_proc        Minimum point global id per process.
        /// @param [in] cells_per_dim_vec                 Total child cells in each direction (x-,y-, and z-direction) per block of cells.
        void assignCellIdsAndCandidatePointIds( std::vector<std::vector<int>>& localToGlobal_cells_per_level,
                                                std::vector<std::vector<int>>& localToGlobal_points_per_level,
                                                int min_globalId_cell_in_proc,
                                                int min_globalId_point_in_proc,
                                                const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief Select and re-write point global ids.
        ///
        /// After assigning global IDs to points in refined-level grids, a single point may have
        /// "multiple unique" global IDs, one in each process to which it belongs.
        /// To reduce the unnucesary id assigments, since global IDs must be distinct across the global leaf view
        /// and consistent across each refined-level grid, we will rewrite the entries in
        /// localToGlobal_points_per_level. Using cell_to_point_ across all refined cells through
        /// communication: gathering the 8 corner points of each interior cell and scattering the
        /// 8 corner points of overlapping cells, for all child cells of a parent cell in level zero grid.
        ///
        /// @param [out] localToGlobal_points_per_level   Relation local point.index() to assigned 'candidate' global id.
        /// @param [in] parent_to_children                The communication step is based on level zero grid, via the relation parent-children-cells.
        /// @param [in] cells_per_dim_vec                 Total child cells in each direction (x-,y-, and z-direction) per block of cells.
        void selectWinnerPointIds(std::vector<std::vector<int>>&  localToGlobal_points_per_level,
                                  const std::vector<std::tuple<int,std::vector<int>>>& parent_to_children,
                                  const std::vector<std::array<int,3>>& cells_per_dim_vec) const;

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the cell_index_set_ of each refined level grid.
        void populateCellIndexSetRefinedGrid(int level);

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the cell_index_set_ of the leaf grid view.
        void populateCellIndexSetLeafGridView();

        /// @brief For a grid whose level zero has been distributed and then locally refined, populate the global_id_set_ of the leaf grid view.
        void populateLeafGlobalIdSet();

    public:

        /// @brief Compute for each level grid, a map from the global_cell_[ cell index in level grid ] to the leaf index of the equivalent cell
        ///        on the leaf grid view.
        ///        Notice that cells that vanished and do not appear on the leaf grid view will not be considered.
        ///        global_cell_[ cell index in level grid ] coincide with (local) Cartesian Index.
        std::vector<std::unordered_map<std::size_t, std::size_t>> mapLocalCartesianIndexSetsToLeafIndexSet() const;

        /// @brief Reverse map: from leaf index cell to { level, local/level Cartesian index of the cell }
        std::vector<std::array<int,2>> mapLeafIndexSetToLocalCartesianIndexSets() const;

        /// \brief Size of the overlap on the leaf level
        grid_size_t overlapSize(int) const;


        /// \brief Size of the ghost cell layer on the leaf level
        grid_size_t ghostSize(int) const;

        /// \brief Size of the overlap on a given level
        grid_size_t overlapSize(int, int) const;

        /// \brief Size of the ghost cell layer on a given level
        grid_size_t ghostSize(int, int) const;

        /// \brief returns the number of boundary segments within the macro grid
        grid_size_t numBoundarySegments() const;

        void setPartitioningParams(const std::map<std::string,std::string>& params);

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \warning May only be called once.
        bool loadBalance(int overlapLayers=1, int partitionMethod = Dune::PartitionMethod::zoltan, double imbalanceTol = 1.1)
        {
            using std::get;
            return get<0>(scatterGrid(defaultTransEdgeWgt, false, nullptr, {}, false, nullptr, true, overlapLayers, partitionMethod, imbalanceTol));
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \param edgeWeightMethod The edge-weighting method to be used on the graph partitioner.
        /// \warning May only be called once.
        bool loadBalanceSerial(int overlapLayers=1, int partitionMethod = Dune::PartitionMethod::zoltan, int edgeWeightMethod = Dune::EdgeWeightMethod::defaultTransEdgeWgt, double imbalanceTol = 1.1)
        {
            using std::get;
            return get<0>(scatterGrid(EdgeWeightMethod(edgeWeightMethod), false, nullptr, {}, true /*serial partitioning*/, nullptr, true, overlapLayers, partitionMethod, imbalanceTol));
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified as weights associated with its edges. The graph will be passed to the load balancer.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
        /// \param transmissibilities The transmissibilities used as the edge weights.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool,std::vector<std::pair<std::string,bool>>>
        loadBalance(const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections = {},
                    const double* transmissibilities = nullptr,
                    int overlapLayers=1, int partitionMethod=Dune::PartitionMethod::zoltan)
        {
            return scatterGrid(defaultTransEdgeWgt, false, wells, possibleFutureConnections, false, transmissibilities, false, overlapLayers, partitionMethod);
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified to calculate the  weights associated with its edges. The graph will be passed
        ///  to the load balancer.
        /// \param method The edge-weighting method to be used on the graph partitioner.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool,std::vector<std::pair<std::string,bool>>>
        loadBalance(EdgeWeightMethod method, const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections = {},
                    const double* transmissibilities = nullptr, bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1,
                    int partitionMethod = Dune::PartitionMethod::zoltan,
                    double imbalanceTol = 1.1)
        {
            return scatterGrid(method, ownersFirst, wells, possibleFutureConnections, false, transmissibilities, addCornerCells, overlapLayers, partitionMethod, imbalanceTol);
        }

        /// \brief Distributes this grid and data over the available nodes in a distributed machine.
        /// \param data A data handle describing how to distribute attached data.
        /// \param wells The wells of the eclipse  Default: null
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param overlapLayers The number of layers of overlap cells to be added
        ///        (default: 1)
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections = {},
                    const double* transmissibilities = nullptr,
                    int overlapLayers=1, int partitionMethod = 1)
        {
            auto ret = loadBalance(wells, possibleFutureConnections, transmissibilities, overlapLayers, partitionMethod);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified to calculate the  weights associated with its edges. The graph will be passed
        ///  to the load balancer.
        /// \param data A data handle describing how to distribute attached data.
        /// \param method The edge-weighting method to be used on the graph partitioner.
        /// \param wells The information about all possible wells. If null then
        ///            the wells will be neglected. Otherwise the wells will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This is done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
        /// \param serialPartitioning If true, the partitioning will be done on a single process.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \param imbalanceTol Set the imbalance tolerance used by the partitioner
        /// \param allowDistributedWells Allow the perforation of a well to be distributed to the
        ///        interior region of multiple processes.
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data, EdgeWeightMethod method,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                    bool serialPartitioning,
                    const double* transmissibilities = nullptr, bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1, int partitionMethod = Dune::PartitionMethod::zoltan,
                    double imbalanceTol = 1.1,
                    bool allowDistributedWells = false)
        {
            auto ret = scatterGrid(method, ownersFirst, wells, possibleFutureConnections, serialPartitioning, transmissibilities,
                                   addCornerCells, overlapLayers, partitionMethod, imbalanceTol, allowDistributedWells);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param data A data handle describing how to distribute attached data.
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data, const std::vector<int>& parts,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections = {},
                    bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1)
        {
            using std::get;
            auto ret = scatterGrid(defaultTransEdgeWgt,  ownersFirst, wells,
                                   possibleFutureConnections,
                                   /* serialPartitioning = */ false,
                                   /* transmissibilities = */ {},
                                   addCornerCells, overlapLayers, /* partitionMethod =*/ Dune::PartitionMethod::simple,
                                   /* imbalanceTol (ignored) = */ 0.0,
                                   /* allowDistributedWells = */ true, parts);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid and data over the available nodes in a distributed machine.
        /// \param data A data handle describing how to distribute attached data.
        /// \param overlapLayers The number of layers of overlap cells to be added
        ///        (default: 1)
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        template<class DataHandle>
        bool loadBalance(DataHandle& data,
                         decltype(data.fixedSize(0,0)) overlapLayers=1, int partitionMethod = Dune::PartitionMethod::zoltan)
        {
            // decltype usage needed to tell the compiler not to use this function if first
            // argument is std::vector but rather loadbalance by parts
            bool ret = loadBalance(overlapLayers, partitionMethod);
            if (ret)
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \warning May only be called once.
        bool loadBalance(const std::vector<int>& parts, bool ownersFirst=false,
                         bool addCornerCells=false, int overlapLayers=1)
        {
            using std::get;
            return get<0>(scatterGrid(defaultTransEdgeWgt,  ownersFirst, /* wells = */ {},
                                      {},
                                      /* serialPartitioning = */ false,
                                      /* trabsmissibilities = */ {},
                                      addCornerCells, overlapLayers, /* partitionMethod =*/ Dune::PartitionMethod::simple,
                                      /* imbalanceTol (ignored) = */ 0.0,
                                      /* allowDistributedWells = */ true, parts));
        }

        /// \brief Distributes this grid and data over the available nodes in a distributed machine
        /// \param data A data handle describing how to distribute attached data.
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \warning May only be called once.
        template<class DataHandle>
        bool loadBalance(DataHandle& data, const std::vector<int>& parts, bool ownersFirst=false,
                         bool addCornerCells=false, int overlapLayers=1)
        {
            bool ret = loadBalance(parts, ownersFirst, addCornerCells, overlapLayers);
            if (ret)
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Partitions the grid using Zoltan without decomposing and distributing it among processes.
        /// \param wells The wells of the eclipse.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
         /// \param transmissibilities The transmissibilities used to calculate the edge weights.
         /// \param numParts Number of parts in the partition.
         /// \return An array with the domain index for each cell.
         std::vector<int>
         zoltanPartitionWithoutScatter(const std::vector<cpgrid::OpmWellType>* wells,
                                       const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                       const double* transmissibilities,
                                       const int     numParts,
                                       const double  imbalanceTol) const;

        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level
        /// \param data The data handle describing the data. Has to adhere to the
        /// Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int /*level*/) const
        {
            communicate(data, iftype, dir);
        }

        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level.
        /// \tparam DataHandle The type of the data handle describing the data.
        /// \param data The data handle describing the data. Has to adhere to the Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const;

        /// \brief Get the collective communication object.
        const typename CpGridTraits::Communication& comm () const;
        //@}

        // ------------ End of Dune interface, start of simplified interface --------------

        /// \name The simplified grid interface.
        ///
        /// It provides additional methods not in the DUNE interface but needed by OPM.
        /// Vertices, faces, and cells are not represented as entities but identified by
        /// indices.
        //@{
        // enum { dimension = 3 }; // already defined

        typedef Dune::FieldVector<double, 3> Vector;


        const std::vector<double>& zcornData() const;


        // Topology
        /// \brief Get the number of cells.
        int numCells() const;

        /// \brief Get the number of faces.
        int numFaces() const;

        /// \brief Get The number of vertices.
        int numVertices() const;


        /// \brief Get the number of faces of a cell.
        ///
        /// Due to faults, and collapsing vertices (along pillars) this
        /// number is quite arbitrary. Its lower bound is 4, but there is
        /// no upper bound.
        /// \parame cell the index identifying the cell.
        int numCellFaces(int cell) const;

        /// \brief Get a specific face of a cell.
        /// \param cell The index identifying the cell.
        /// \param local_index The local index (in [0,numFaces(cell))) of the face in this cell.
        /// \return The index identifying the face.
        int cellFace(int cell, int local_index) const;

        /// \brief Get a list of indices identifying all faces of a cell.
        /// \param cell The index identifying the cell.
        const cpgrid::OrientedEntityTable<0,1>::row_type cellFaceRow(int cell) const;

        /// \brief Get the index identifying a cell attached to a face.
        ///
        /// Note that a face here is always oriented. If there are two
        /// neighboring cells then the orientation will be from local_index 0
        /// to local_index 1
        /// \param face The index identifying the face.
        /// \param local_index The local_index of the cell.
        /// \return The index identifying a cell or -1 if there is no such
        /// cell due the face being part of the grid boundary or the
        /// cell being stored on another process.
        int faceCell(int face, int local_index) const;
      
        /// \brief Get the sum of all faces attached to all cells.
        ///
        /// Each face identified by a unique index is counted as often
        /// as there are neigboring cells attached to it.
        /// \f$ numCellFaces()=\sum_{c} numCellFaces(c) \f$
        /// \see numCellFaces(int)const
        int numCellFaces() const;

        int numFaceVertices(int face) const;

        /// \brief Get the index identifying a vertex of a face.
        /// \param cell The index identifying the face.
        /// \param local_index The local_index (in [0,numFaceVertices(vertex) - 1]])
        ///  of the vertex.
        int faceVertex(int face, int local_index) const;

        /// \brief Get vertical position of cell center ("zcorn" average).
        /// \brief cell_index The index of the specific cell.
        double cellCenterDepth(int cell_index) const;


        const Vector faceCenterEcl(int cell_index, int face, const Dune::cpgrid::Intersection& intersection) const;

        const Vector faceAreaNormalEcl(int face) const;


        // Geometry
        /// \brief Get the Position of a vertex.
        /// \param cell The index identifying the cell.
        /// \return The coordinates of the vertex.
        const Vector& vertexPosition(int vertex) const;

        /// \brief Get the area of a face.
        /// \param cell The index identifying the face.
        double faceArea(int face) const;

        /// \brief Get the coordinates of the center of a face.
        /// \param cell The index identifying the face.
        const Vector& faceCentroid(int face) const;

        /// \brief Get the unit normal of a face.
        /// \param cell The index identifying the face.
        /// \see faceCell
        const Vector& faceNormal(int face) const;

        /// \brief Get the volume of the cell.
        /// \param cell The index identifying the cell.
        double cellVolume(int cell) const;

        /// \brief Get the coordinates of the center of a cell.
        /// \param cell The index identifying the face.
        const Vector& cellCentroid(int cell) const;

        /// \brief An iterator over the centroids of the geometry of the entities.
        /// \tparam codim The co-dimension of the entities.
        template<int codim>
        class CentroidIterator
            : public RandomAccessIteratorFacade<CentroidIterator<codim>,
                                                FieldVector<double, 3>,
                                                const FieldVector<double, 3>&, int>
        {
        public:
            /// \brief The type of the iterator over the codim geometries.
            typedef typename std::vector<cpgrid::Geometry<3-codim, 3> >::const_iterator
            GeometryIterator;
            /// \brief Constructs a new iterator from an iterator over the geometries.
            /// \param iter The iterator of the geometry objects.
            CentroidIterator(GeometryIterator iter)
            : iter_(iter)
            {}

            const FieldVector<double,3>& dereference() const
            {
                return iter_->center();
            }
            void increment()
            {
                ++iter_;
            }
            const FieldVector<double,3>& elementAt(int n)
            {
                return iter_[n]->center();
            }
            void advance(int n){
                iter_+=n;
            }
            void decrement()
            {
                --iter_;
            }
            int distanceTo(const CentroidIterator& o)
            {
                return o-iter_;
            }
            bool equals(const CentroidIterator& o) const{
                return o==iter_;
            }
        private:
            /// \brief The iterator over the underlying geometry objects.
            GeometryIterator iter_;
        };

        /// \brief Get an iterator over the cell centroids positioned at the first one.
        CentroidIterator<0> beginCellCentroids() const;

        /// \brief Get an iterator over the face centroids positioned at the first one.
        CentroidIterator<1> beginFaceCentroids() const;

        // Extra
        int boundaryId(int face) const;

        /// \brief Get the cartesian tag associated with a face tag.
        ///
        /// The tag tells us in which direction the face would point
        /// in the underlying cartesian grid.
        /// \param An iterator that points to the face and was obtained
        /// by iterating over Opm::UgGridHelpers::cell2Faces(grid).
        template<class Cell2FacesRowIterator>
        int
        faceTag(const Cell2FacesRowIterator& cell_face) const;

        //@}

        // ------------ End of simplified interface --------------

        //------------- methods not in the DUNE grid interface.

        /// \name Parallel grid extensions.
        /// Methods extending the DUNE's parallel grid interface.
        /// These are basically for scattering/gathering data to/from
        /// distributed views.
        //@{
        ///
        /// \brief Moves data from the global (all data on process) view to the distributed view.
        ///
        /// This method does not do communication but assumes that the global grid
        /// is present on every process and simply copies data to the distributed view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void scatterData(DataHandle& handle) const;

        ///
        /// \brief Moves data from the distributed view to the global (all data on process) view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void gatherData(DataHandle& handle) const;

        /// \brief The type of the map describing communication interfaces.
        using InterfaceMap = cpgrid::CpGridDataTraits::InterfaceMap;

        /// \brief Get an interface for gathering/scattering data attached to cells with communication.
        ///
        /// Scattering means sending data from the indices of the global grid on
        /// process 0 to the distributed grid on all ranks independent of the grid.
        /// Gathering is the other way around.
        /// The interface can be used with VariableSizeCommunicator and a custom
        /// index based data handle to scatter (forward direction of the communicator)
        /// and gather data (backward direction of the communicator).
        /// Here is a small example that prints the received values when scattering:
        /// \code
        /// struct Handle{
        ///   typedef int DataType;
        ///   const std::vector<int>& vals;
        ///   bool fixedsize() { return true; }
        ///   size_t size(std::size_t) { return 1; }
        ///   void gather(auto& B buf, size_t i)[ buf.write(vals[i]); }
        ///   void scatter(auto& B buf, size_t i, std::size_t) {
        ///     int val;
        ///     buf.read(val);
        ///     cout<<i<<": "<<val<<" "; }
        /// };
        ///
        /// Handle handle;
        /// handle.vals.resize(grid.size(0), -1);
        /// Dune::VariableSizeCommunicator<> comm(grid.comm(),
        ///                                       grid.cellScatterGatherInterface());
        /// comm.forward(handle);
        /// \endcode
        const InterfaceMap& cellScatterGatherInterface() const;

        /// \brief Get an interface for gathering/scattering data attached to points with communication.
        /// \see cellScatterGatherInterface
        const InterfaceMap& pointScatterGatherInterface() const;

        /// \brief Switch to the global view.
        void switchToGlobalView();

        /// \brief Switch to the distributed view.
        void switchToDistributedView();
        //@}

#if HAVE_MPI
        /// \brief The type of the parallel index set
        using ParallelIndexSet = cpgrid::CpGridDataTraits::ParallelIndexSet;
        /// \brief The type of the remote indices information
        using RemoteIndices = cpgrid::CpGridDataTraits::RemoteIndices;

        /// \brief The type of the owner-overlap-copy communication
        using CommunicationType = cpgrid::CpGridDataTraits::CommunicationType;

        /// \brief Get the owner-overlap-copy communication for cells
        ///
        /// Suitable e.g. for parallel linear algebra used by CCFV
        const CommunicationType& cellCommunication() const;

        ParallelIndexSet& getCellIndexSet();

        RemoteIndices& getCellRemoteIndices();

        const ParallelIndexSet& getCellIndexSet() const;

        const RemoteIndices& getCellRemoteIndices() const;
#endif

        /// \brief Get sorted active cell indices of numerical aquifer
        const std::vector<int>& sortedNumAquiferCells() const;

    private:
        /// \brief Scatter a global grid to all processors.
        /// \param method The edge-weighting method to be used on the graph partitioner.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param possibleFutureConnections An optional unordered_map<string, set<array<int,3>>>
        ///            containing possible future connections that might be opened during an ACTIONX.
        ///            The fist entry is the name of the well and the second entry is a set containing
        ///            the cartesian coordinates of the grid cells that get perforated of a possible
        ///            future connection. The possible future connections are handed over to the grid
        ///            partitioner to make sure these will be no the same partition when partitioning
        ///            the grid.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights in
        ///                           the graph partitioner. This is done to improve the numerical
        ///                           performance of the parallel preconditioner.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region.
        /// \param partitionMethod The method used to partition the grid, one of Dune::PartitionMethod
        /// \param imbalanceTol Set the imbalance tolerance used by the partitioner
        /// \param allowDistributedWells Allow the perforation of a well to be distributed to the
        ///        interior region of multiple processes.
        /// \param cell_part When using an external loadbalancer the partition number for each cell.
        ///                  If empty or not specified we use internal load balancing.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        scatterGrid(EdgeWeightMethod method,
                    bool ownersFirst,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                    bool serialPartitioning,
                    const double* transmissibilities,
                    bool addCornerCells,
                    int overlapLayers,
                    int partitionMethod = Dune::PartitionMethod::zoltan,
                    double imbalanceTol = 1.1,
                    bool allowDistributedWells = true,
                    const std::vector<int>& input_cell_part = {});

        /** @brief The data stored in the grid.
         *
         * All the data of all grids are stored there and
         * calls are forwarded to relevant grid.*/
        std::vector<std::shared_ptr<cpgrid::CpGridData>> data_;
        /** @brief A pointer to data of the current View. */
        cpgrid::CpGridData* current_view_data_;
        /** @brief The data stored for the distributed grid. */
        std::vector<std::shared_ptr<cpgrid::CpGridData>> distributed_data_;
        /** @brief A pointer to the current data used. */
        std::vector<std::shared_ptr<cpgrid::CpGridData>>* current_data_;
        /** @brief To get the level given the lgr-name. Default, {"GLOBAL", 0}. */
        std::map<std::string,int> lgr_names_ = {{"GLOBAL", 0}};
        /**
         * @brief Interface for scattering and gathering cell data.
         *
         * @warning Will only update owner cells.
         */
        std::shared_ptr<InterfaceMap> cell_scatter_gather_interfaces_;
        /*
         * @brief Interface for scattering and gathering point data.
         *
         * @warning Will only update owner cells
         */
        std::shared_ptr<InterfaceMap> point_scatter_gather_interfaces_;
        /**
         * @brief The global id set (also used as local one).
         */
        std::shared_ptr<cpgrid::GlobalIdSet> global_id_set_ptr_;


        /**
         * @brief Partitioning parameters
         */
        std::map<std::string,std::string> partitioningParams;

    }; // end Class CpGrid

} // end namespace Dune

#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/Iterators.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>


namespace Dune
{

    namespace Capabilities
    {
        /// \todo Please doc me !
        template <>
        struct hasEntity<CpGrid, 0>
        {
            static const bool v = true;
        };

        /// \todo Please doc me !
        template <>
        struct hasEntity<CpGrid, 3>
        {
            static const bool v = true;
        };

        template<>
        struct canCommunicate<CpGrid,0>
        {
            static const bool v = true;
        };

        template<>
        struct canCommunicate<CpGrid,3>
        {
            static const bool v = true;
        };

        /// \todo Please doc me !
        template <>
        struct hasBackupRestoreFacilities<CpGrid>
        {
            static const bool v = false;
        };

    }

    template<class DataHandle>
    void CpGrid::communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {
        current_view_data_->communicate(data, iftype, dir);
    }


    template<class DataHandle>
    void CpGrid::scatterData(DataHandle& handle) const
    {
#if HAVE_MPI
        if(distributed_data_.empty())
            OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balanced grid!");
        distributed_data_[0]->scatterData(handle, data_[0].get(), distributed_data_[0].get(), cellScatterGatherInterface(),
                                          pointScatterGatherInterface());
#else
        // Suppress warnings for unused argument.
        (void) handle;
#endif
    }

    template<class DataHandle>
    void CpGrid::gatherData(DataHandle& handle) const
    {
#if HAVE_MPI
        if(distributed_data_.empty())
            OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balance grid!");
        distributed_data_[0]->gatherData(handle, data_[0].get(), distributed_data_[0].get());
#else
        // Suppress warnings for unused argument.
        (void) handle;
#endif
    }


    template<class Cell2FacesRowIterator>
    int
    CpGrid::faceTag(const Cell2FacesRowIterator& cell_face) const
    {
        // Note that this relies on the following implementation detail:
        // The grid is always constructed such that the interior faces constructed
        // with orientation set to true are
        // oriented along the positive IJK direction. Oriented means that
        // the first cell attached to face has the lower index.
        // For faces along the boundary (only one cell, always  attached at index 0)
        // the orientation has to be determined by the orientation of the cell.
        // If it is true then in UnstructuredGrid it would be stored at index 0,
        // otherwise at index 1.
        const int cell = cell_face.getCellIndex();
        const int face = *cell_face;
        assert (0 <= cell);  assert (cell < numCells());
        assert (0 <= face);  assert (face < numFaces());

        typedef cpgrid::OrientedEntityTable<1,0>::row_type F2C;

        const cpgrid::EntityRep<1> f(face, true);
        const F2C&     f2c = current_view_data_->face_to_cell_[f];
        const face_tag tag = current_view_data_->face_tag_[f];

        assert ((f2c.size() == 1) || (f2c.size() == 2));

        int inside_cell = 0;

        if ( f2c.size() == 2 ) // Two cells => interior
        {
            if ( f2c[1].index() == cell )
            {
                inside_cell = 1;
            }
        }
        const bool normal_is_in = ! f2c[inside_cell].orientation();

        switch (tag) {
        case I_FACE:
            //                    LEFT : RIGHT
            return normal_is_in ? 0    : 1; // min(I) : max(I)
        case J_FACE:
            //                    BACK : FRONT
            return normal_is_in ? 2    : 3; // min(J) : max(J)
        case K_FACE:
            // Note: TOP at min(K) as 'z' measures *depth*.
            //                    TOP  : BOTTOM
            return normal_is_in ? 4    : 5; // min(K) : max(K)
        case NNC_FACE:
            // For nnc faces we return the otherwise unused value -1.
            return -1;
        default:
            OPM_THROW(std::logic_error, "Unhandled face tag. This should never happen!");
        }
    }

    template<int dim>
    cpgrid::Entity<dim> createEntity(const CpGrid&, int, bool);

} // namespace Dune

#include <opm/grid/cpgrid/PersistentContainer.hpp>
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#include "cpgrid/Intersection.hpp"
#include "cpgrid/Geometry.hpp"
#include "cpgrid/Indexsets.hpp"

#endif // OPM_CPGRID_HEADER
