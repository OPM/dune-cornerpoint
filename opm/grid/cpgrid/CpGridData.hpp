//===========================================================================
//
// File: CpGridData.hpp
//
// Created: Sep 17 21:11:41 2013
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//            Markus Blatt        <markus@dr-blatt.de>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// Comment: Major parts of this file originated in dune/grid/CpGrid.hpp
//          and got transfered here during refactoring for the parallelization.
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013, 2022-2023 Equinor ASA.
  Copyright 2013 Dr. Blatt - HPC-Simulation-Software & Services

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
/**
 * @brief Holds the implementation of the CpGrid as a pimple.
 * @author Markus Blatt <markus@dr-blatt.de>
 *         Atgeirr F Rasmussen <atgeirr@sintef.no>
 *         Bård Skaflestad     <bard.skaflestad@sintef.no>
 */

#ifndef OPM_CPGRIDDATA_HEADER
#define OPM_CPGRIDDATA_HEADER


#include <dune/common/parallel/mpihelper.hh>
#ifdef HAVE_DUNE_ISTL
#include <dune/istl/owneroverlapcopy.hh>
#endif

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/variablesizecommunicator.hh>
#include <dune/grid/common/gridenums.hh>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#endif

#include <opm/grid/cpgpreprocess/preprocess.h>

#include "Entity2IndexDataHandle.hpp"
#include "CpGridDataTraits.hpp"
//#include "DataHandleWrappers.hpp"
//#include "GlobalIdMapping.hpp"
#include "Geometry.hpp"

#include <array>
#include <initializer_list>
#include <set>
#include <vector>

namespace Opm
{
class EclipseState;
}
namespace Dune
{
class CpGrid;

namespace cpgrid
{

class IndexSet;
class IdSet;
class LevelGlobalIdSet;
class PartitionTypeIndicator;
template<int,int> class Geometry;
template<int> class Entity;
template<int> class EntityRep;
}
}


void markAndAdapt_check(Dune::CpGrid&,
                        const std::array<int,3>&,
                        const std::vector<int>&,
                        Dune::CpGrid&,
                        bool,
                        bool,
                        bool);

void refine_and_check(const Dune::cpgrid::Geometry<3, 3>&,
                      const std::array<int, 3>&,
                      bool);

void refinePatch_and_check(const std::array<int,3>&,
                           const std::array<int,3>&,
                           const std::array<int,3>&);

void refinePatch_and_check(Dune::CpGrid&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::string>&);

void check_global_refine(const Dune::CpGrid&,
                         const Dune::CpGrid&);

void fieldProp_check(const Dune::CpGrid& grid, Opm::EclipseGrid eclGrid, const std::string& deck_string);

void testInactiveCellsLgrs(const std::string&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::array<int,3>>&,
                           const std::vector<std::string>&);

namespace Dune
{
namespace cpgrid
{
namespace mover
{
template<class T, int i> struct Mover;
}

/**
 * @brief Struct that hods all the data needed to represent a
 * Cpgrid.
 */
class CpGridData
{
    template<class T, int i> friend struct mover::Mover;
    friend class GlobalIdSet;
    friend class HierarchicIterator;
    friend class Dune::cpgrid::IndexSet;
    friend class Dune::cpgrid::IdSet;
    friend class Dune::cpgrid::LevelGlobalIdSet;

    friend
    void ::markAndAdapt_check(Dune::CpGrid&,
                              const std::array<int,3>&,
                              const std::vector<int>&,
                              Dune::CpGrid&,
                              bool,
                              bool,
                              bool);

    friend
    void ::refine_and_check(const Dune::cpgrid::Geometry<3, 3>&,
                            const std::array<int, 3>&,
                            bool);
    friend
    void ::refinePatch_and_check(const std::array<int,3>&,
                                 const std::array<int,3>&,
                                 const std::array<int,3>&);

    friend
    void ::refinePatch_and_check(Dune::CpGrid&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::string>&);

    friend
    void ::check_global_refine(const Dune::CpGrid&,
                               const Dune::CpGrid&);
    friend
    void ::fieldProp_check(const Dune::CpGrid& grid, Opm::EclipseGrid eclGrid, const std::string& deck_string);
    friend
    void ::testInactiveCellsLgrs(const std::string&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::array<int,3>>&,
                                 const std::vector<std::string>&);

private:
    CpGridData(const CpGridData& g);

public:
    enum{
#ifndef MAX_DATA_COMMUNICATED_PER_ENTITY
        /// \brief The maximum data items allowed per cell (DUNE < 2.5.2)
        ///
        /// Due to a bug in DUNE < 2.5.2 we need to limit this when
        /// communicating. 1 is big enough for OPM as we always use
        /// one block for all unknowns, but some DUNE's grid checks
        /// actually need 2. So 2 it is.
        MAX_DATA_PER_CELL = 2
#else
        /// \brief The maximum data items allowed per cell (DUNE < 2.5.2)
        ///
        /// Due to a bug in DUNE < 2.5.2 we need to limit this when
        /// communicating. Uses the define MAX_DATA_COMMUNICATED_PER_ENTITY.
        MAX_DATA_PER_CELL = MAX_DATA_COMMUNICATED_PER_ENTITY
#endif
    };

    /// Constructor for parallel grid data.
    /// \param comm The MPI communicator
    /// Default constructor.
    explicit CpGridData(MPIHelper::MPICommunicator comm,  std::vector<std::shared_ptr<CpGridData>>& data);



    /// Constructor
    explicit CpGridData(std::vector<std::shared_ptr<CpGridData>>& data);
    /// Destructor
    ~CpGridData();




    /// number of leaf entities per codim in this process
    int size(int codim) const;

    /// number of leaf entities per geometry type in this process
    int size (GeometryType type) const
    {
        if (type.isCube()) {
            return size(3 - type.dim());
        } else {
            return 0;
        }
    }

    /// Read the Eclipse grid format ('grdecl').
    ///
    /// \param[in] filename the name of the file to read.
    ///
    /// \param[in] periodic_extension Whether or not to process the
    /// resulting grid in order to have intersections/faces along i and j
    /// boundaries match those on the other side. That is, i- faces will
    /// match i+ faces etc.
    ///
    /// \param[in] turn_normals Whether or not to turn all normals.  This is
    /// intended for handling inputs with wrong orientations.
    ///
    /// \param[in] edge_conformal Whether or not to construct an
    /// edge-conformal grid.  Typically useful in geo-mechanical
    /// applications.
    void readEclipseFormat(const std::string& filename,
                           bool periodic_extension,
                           bool turn_normals = false,
                           bool edge_conformal = false);

#if HAVE_ECL_INPUT
    /// Read the Eclipse grid format ('grdecl').
    ///
    /// \param[in] deck the parsed deck from opm-parser (which is a
    /// low-level object)
    ///
    /// \param[in] periodic_extension Whether or not to process the
    /// resulting grid in order to have intersections/faces along i and j
    /// boundaries match those on the other side. That is, i- faces will
    /// match i+ faces etc.
    ///
    /// \param[in] turn_normals Whether or not to turn all normals.  This is
    /// intended for handling inputs with wrong orientations.
    ///
    /// \param[in] clip_z Whether or not to clip the result grid in
    /// order to have planar top and bottom surfaces.
    ///
    /// \param[in] poreVolume pore volumes for use in MINPV processing, if
    /// asked for in deck
    ///
    /// \param[in] edge_conformal Whether or not to construct an
    /// edge-conformal grid.  Typically useful in geo-mechanical
    /// applications.
    void processEclipseFormat(const Opm::Deck& deck,
                              bool periodic_extension,
                              bool turn_normals = false,
                              bool clip_z = false,
                              const std::vector<double>& poreVolume = std::vector<double>{},
                              bool edge_conformal = false);

    /// Read the Eclipse grid format ('grdecl').
    ///
    /// \param[in] ecl_grid the high-level object from opm-parser which
    /// represents the simulation's grid.  In a parallel run this may be
    /// a nullptr on all ranks but rank zero.
    ///
    /// \param[in,out] ecl_state the object from opm-parser provide
    /// information regarding to pore volume, NNC, aquifer information
    /// when ecl_state is available.  NNC and aquifer connection
    /// nformation will also be updated during the function call when
    /// available and necessary.
    ///
    /// \param[in] periodic_extension Whether or not to process the
    /// resulting grid in order to have intersections/faces along i and
    /// j boundaries match those on the other side. That is, i- faces
    /// will match i+ faces etc.
    ///
    /// \param[in] turn_normals Whether or not to turn all normals.
    /// This is intended for handling inputs with wrong orientations.
    ///
    /// \param[in] clip_z Whether or not to clip the result grid in
    /// order to have planar top and bottom surfaces.
    ///
    /// \param[in] pinchActive Whether or not to force specific pinch
    /// behaviour.  If set, a face will connect two vertical cells, that
    /// are topological connected, even if there are cells with zero
    /// volume between them. If false these cells will not be connected
    /// despite their faces coinciding.
    ///
    /// \param[in] edge_conformal Whether or not to construct an
    /// edge-conformal grid.  Typically useful in geo-mechanical
    /// applications.
    ///
    /// \return Cells removed due low pore-volume and across which to create
    /// non-neighbouring connections if item 4 of the 'PINCH' keyword is set
    /// to 'ALL'.
    std::vector<std::size_t>
    processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                         Opm::EclipseState* ecl_state,
                         bool periodic_extension,
                         bool turn_normals = false,
                         bool clip_z = false,
                         bool pinchActive = true,
                         bool edge_conformal = false);
#endif

    /// Read the Eclipse grid format ('grdecl').
    ///
    /// \param input_data the data in grdecl format, declared in
    /// preprocess.h.
    ///
    /// \param[in,out] ecl_state the object from opm-parser provide
    /// information regarding to pore volume, NNC, aquifer information
    /// when ecl_state is available.  NNC and aquifer connection
    /// nformation will also be updated during the function call when
    /// available and necessary.
    ///
    /// \param[in,out] nnc Non-neighboring connections.
    ///
    /// \param[in] turn_normals Whether or not to turn all normals.
    /// This is intended for handling inputs with wrong orientations.
    ///
    /// \param[in] pinchActive Whether or not to force specific pinch
    /// behaviour.  If set, a face will connect two vertical cells, that are
    /// topological connected, even if there are cells with zero volume
    /// between them. If false these cells will not be connected despite
    /// their faces coinciding.
    ///
    /// \param[in] tolerance_unique_points Tolerance used to identify points
    /// based on their cooridinate.
    ///
    /// \param[in] edge_conformal Whether or not to construct an
    /// edge-conformal grid.  Typically useful in geo-mechanical
    /// applications.
    void processEclipseFormat(const grdecl& input_data,
#if HAVE_ECL_INPUT
                              Opm::EclipseState* ecl_state,
#endif
                              std::array<std::set<std::pair<int, int>>, 2>& nnc,
                              bool remove_ij_boundary,
                              bool turn_normals,
                              bool pinchActive,
                              double tolerance_unique_points,
                              bool edge_conformal);

    /// @brief
    ///    Extract Cartesian index triplet (i,j,k) of an active cell.
    ///
    /// @param [in] c
    ///    Active cell index.
    ///
    /// @param [out] ijk  Cartesian index triplet
    void getIJK(int c, std::array<int,3>& ijk) const
    {
        ijk = getIJK(global_cell_[c], logical_cartesian_size_);
    }

    int cellFace(int cell, int local_index) const
    {
        return cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
    }

    /// Return global_cell_ of any level grid, or the leaf grid view (in presence of refinement).
    /// global_cell_ has size number of cells present on a process and maps to the underlying Cartesian Grid.
    ///
    /// Note: CpGrid::globalCell() returns current_view_data_-> global_cell_ (current_view_data_ points at
    /// data_.back() or distributed_data_.back(), in general. If the grid has been refined, current_view_data_
    /// points at the "leaf grid view").
    const std::vector<int>& globalCell() const
    {
        return  global_cell_;
    }

    /// @brief Extract Cartesian index triplet (i,j,k) given an index between 0 and NXxNYxNZ -1
    ///    where NX, NY, and NZ is the total amoung of cells in each direction x-,y-,and z- respectively.
    ///
    /// @param [in] idx      Integer between 0 and cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]-1
    /// @param [in] cells_per_dim
    /// @return Cartesian index triplet.
    std::array<int,3> getIJK(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim) const
    {
        // idx = k*cells_per_dim_[0]*cells_per_dim_[1] + j*cells_per_dim_[0] + i
        // with 0<= i < cells_per_dim_[0], 0<= j < cells_per_dim_[1], 0<= k <cells_per_dim_[2].
        assert(cells_per_dim[0]);
        assert(cells_per_dim[1]);
        assert(cells_per_dim[2]);

        std::array<int,3> ijk = {0,0,0};
        ijk[0] = idx_in_parent_cell % cells_per_dim[0]; idx_in_parent_cell /= cells_per_dim[0];
        ijk[1] = idx_in_parent_cell % cells_per_dim[1];
        ijk[2] = idx_in_parent_cell /cells_per_dim[1];
        return ijk;
    }

    /// @brief Determine if a finite amount of patches (of cells) are disjoint, namely, they do not share any corner nor face.
    ///
    /// @param [in]  startIJK_vec  Vector of Cartesian triplet indices where each patch starts.
    /// @param [in]  endIJK_vec    Vector of Cartesian triplet indices where each patch ends.
    ///                            Last cell part of the lgr will be {endIJK_vec[<patch>][0]-1, ... ,endIJK_vec[<patch>][2]-1}.
    bool disjointPatches(const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const;

    /// @brief Compute cell indices of selected patches of cells (Cartesian grid required).
    ///
    /// @param [in]  startIJK_vec  Vector of Cartesian triplet indices where each patch starts.
    /// @param [in]  endIJK_vec    Vector of Cartesian triplet indices where each patch ends.
    ///                            Last cell part of the lgr will be {endIJK_vec[<patch>][0]-1, ... endIJK_vec[<patch>][2]-1}.
    ///
    /// @return allPatches_cells
    std::vector<int>
    getPatchesCells(const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const;

    /// @brief Check all cells selected for refinement have no NNCs (no neighbor connections).
    ///        Assumption: all grid cells are active.
    bool hasNNCs(const std::vector<int>& cellIndices) const;

    /// @brief Check startIJK and endIJK of each patch of cells to be refined are valid, i.e.
    ///        startIJK and endIJK vectors have the same size and, startIJK < endIJK coordenate by coordenate.
    ///
    /// @param [in]  startIJK_vec       Vector of Cartesian triplet indices where each patch starts.
    /// @param [in]  endIJK_vec         Vector of Cartesian triplet indices where each patch ends.
    ///                                 Last cell part of the lgr will be {endIJK_vec[patch][0]-1, ..., endIJK_vec[patch][2]-1}.
    void validStartEndIJKs(const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const;

    /// @brief Check that every cell to be refined has cuboid shape.
    void checkCuboidShape(const std::vector<int>& cellIdx_vec) const;

    /// @brief Determine if a finite amount of patches (of cells) share a face.
    ///
    /// @param [in]  startIJK_vec  Vector of Cartesian triplet indices where each patch starts.
    /// @param [in]  endIJK_vec    Vector of Cartesian triplet indices where each patch ends.
    ///                            Last cell part of the lgr will be {endIJK_vec[<patch>][0]-1, ... ,endIJK_vec[<patch>][2]-1}.
    bool patchesShareFace(const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const;

    int sharedFaceTag(const std::vector<std::array<int,3>>& startIJK_2Patches, const std::vector<std::array<int,3>>& endIJK_2Patches) const;


    /// @brief Mark entity for refinement or coarsening.
    ///
    /// Refinement on CpGrid is partially supported for Cartesian grids, with the keyword CARFIN.
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
    /// @return refinement mark (1 refinement, 0 doing nothing, -1 coarsening - not supported yet).
    int getMark(const cpgrid::Entity<0>& element) const;

    /// @brief Set mightVanish flags for elements that will be refined in the next adapt() call
    ///        Need to be called after elements have been marked for refinement.
    bool preAdapt();

    /// TO DO: Documentation. Triggers the grid refinement process - Currently, returns preAdapt()
    bool adapt();

    /// @brief Clean up refinement/coarsening markers - set every element to the mark 0 which represents 'doing nothing'
    void postAdapt();

private:
    /// @brief Check compatibility of number of subdivisions of neighboring LGRs.
    ///
    /// Check shared faces on boundaries of LGRs. Not optimal since the code below does not take into account
    /// active/inactive cells, instead, relies on "ijk-computations".
    ///
    /// @param [in]  cells_per_dim_vec    Vector of expected subdivisions per cell, per direction, in each LGR.
    /// @param [in]  startIJK_vec         Vector of Cartesian triplet indices where each patch starts.
    /// @param [in]  endIJK_vec           Vector of Cartesian triplet indices where each patch ends.
    ///                                   Last cell part of the lgr will be {endIJK_vec[patch][0]-1, ..., endIJK_vec[patch][2]-1}.
    /// @return True if all block of cells either do not share faces on their boundaries, or they may share faces with compatible
    ///         subdivisions. Example: block1 and block2 share an I_FACE, then number of subdivisions NY NZ should coincide, i.e.
    ///         if block1, block2 cells_per_dim values are {NX1, NY1, NZ1}, {NX2, NY2, NZ2}, respectively, then NY1 == NY2 and
    ///         NZ1 == Nz2.
    ///         False if at least two blocks share a face and their subdivions are not compatible. In the example above,
    ///         if NY1 != NY2 or NZ1 != NZ2.
    bool compatibleSubdivisions(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                const std::vector<std::array<int,3>>& startIJK_vec,
                                const std::vector<std::array<int,3>>& endIJK_vec) const;

    std::array<Dune::FieldVector<double,3>,8> getReferenceRefinedCorners(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim) const;

    /// @brief Compute amount of cells in each direction of a patch of cells. (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_dim Patch dimension {#cells in x-direction, #cells in y-direction, #cells in z-direction}.
    std::array<int,3> getPatchDim(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Compute corner indices of a patch of cells (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_corners
    std::vector<int> getPatchCorners(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Compute face indices of a patch of cells (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_faces
    std::vector<int> getPatchFaces(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Compute cell indices of a patch of cells (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endIJK[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_cells
    std::vector<int> getPatchCells(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Compute patch boundary corner indices (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_boundary_corners
    std::vector<int> getPatchBoundaryCorners(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Compute patch boundary face indices (Cartesian grid required).
    ///
    /// @param [in]  startIJK  Cartesian triplet index where the patch starts.
    /// @param [in]  endIJK    Cartesian triplet index where the patch ends.
    ///                        Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return patch_boundary_faces
    std::array<std::vector<int>,6> getBoundaryPatchFaces(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief For selected cell indices, computes the variation in x-,y-, and z-direction, assuming each cell has cubiod shape.
    std::array<std::vector<double>,3> getWidthsLengthsHeights(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    /// @brief Construct a 'fake cell (Geometry<3,3> object)' out of a patch of cells.(Cartesian grid required).
    ///
    /// cellifyPatch() builds a Geometry<3,3> object, 'a celliFIED patch', from a connected patch formed
    /// by the product of consecutive cells in each direction; selecting 8 corners of the patch boundary,
    /// computing center and volume.
    ///
    /// @param [in] startIJK                   Cartesian triplet index where the patch starts.
    /// @param [in] endIJK                     Cartesian triplet index where the patch ends.
    /// @param [in] patch_cells                Cell indices from the block-shaped patch.
    /// @param [out] cellifiedPatch_geometry   Required as an argument when creating a Geomtry<3,3> object.
    /// @param [out] cellifiedPatch_to_point   To store the 8 corners of the created cellifiedPatch.
    /// @param [out] allcorners_cellifiedPatch Required to build a Geometry<3,3> object.
    ///
    /// @return 'cellifiedPatchCell'         Geometry<3,3> object.
    Geometry<3,3> cellifyPatch(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK,
                               const std::vector<int>& patch_cells, DefaultGeometryPolicy& cellifiedPatch_geometry,
                               std::array<int,8>& cellifiedPatch_to_point,
                               std::array<int,8>& allcorners_cellifiedPatch) const;

    // @brief Compute the average of array<double,3>.
    //
    // @param [in] vector of array<double,3>
    // @return     array<double,3> (average of the entries of the given vector).
    std::array<double,3> getAverageArr(const std::vector<std::array<double,3>>& vec) const;

public:
    /// Add doc/or remove method and replace it with better approach
    int getGridIdx() const {
        // Not the nicest way of checking if "this" points at the leaf grid view of a mixed grid (with coarse and refined cells).
        // 1. When the grid has been refined at least onece, level_data_ptr_ ->size() >1. Therefore, there is a chance of "this" pointing at the leaf grid view.
        // 2. Unfortunately, level_ is default initialized by 0. This implies, in particular, that if someone wants to check the value of
        //    "this->level_" when "this" points at the leaf grid view of a grid that has been refined, this value is - unfortunately - equal to 0.
        // 3. Due to 2. we need an extra bool value to distinguish between the actual level 0 grid and such a leaf grid view (with incorrect level_ == 0). For this
        //    reason we check if child_to_parent_cells_.empty() [true for actual level 0 grid, false for the leaf grid view].
        // --- TO BE IMPROVED ---
        if ((level_data_ptr_ ->size() >1) && (level_ == 0) && (!child_to_parent_cells_.empty())) {
            return level_data_ptr_->size() -1;
        }
        return level_;
    }
    /// Add doc/or remove method and replace it with better approach
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> levelData() const
    {
        return *level_data_ptr_;
    }


    /// @brief Refine a single cell and return a shared pointer of CpGridData type.
    ///
    /// refineSingleCell() takes a cell and refines it in a chosen amount of cells (per direction); creating the
    /// geometries, topological relations, etc. Stored in a CpGridData object. Additionally, containers for
    /// parent-to-new-born entities are buil, as well as, new-born-to-parent. Maps(<int,bool>) to detect parent
    /// faces or cells are also provided. (Cell with 6 faces required).
    ///
    /// @param [in] cells_per_dim                 Number of (refined) cells in each direction that each parent cell should be refined to.
    /// @param [in] parent_idx                    Parent cell index, cell to be refined.
    ///
    /// @return refined_grid_ptr                  Shared pointer pointing at refined_grid.
    /// @return parent_to_refined_corners         For each corner of the parent cell, we store the index of the
    ///                                           refined corner that coincides with the old one.
    ///                                           We assume they are ordered 0,1,..7
    ///                                                              6---7
    ///                                                      2---3   |   | TOP FACE
    ///                                                      |   |   4---5
    ///                                                      0---1 BOTTOM FACE
    /// @return parent_to_children_faces/cell     For each parent face/cell, we store its child-face/cell indices.
    ///                                           {parent face/cell index in coarse level, {indices of its children in refined level}}
    /// @return child_to_parent_faces/cells       {child index, parent index}
    std::tuple< const std::shared_ptr<CpGridData>,
                const std::vector<std::array<int,2>>,                // parent_to_refined_corners(~boundary_old_to_new_corners)
                const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces (~boundary_old_to_new_faces)
                const std::tuple<int, std::vector<int>>,             // parent_to_children_cells
                const std::vector<std::array<int,2>>,                // child_to_parent_faces
                const std::vector<std::array<int,2>>>                // child_to_parent_cells
    refineSingleCell(const std::array<int,3>& cells_per_dim, const int& parent_idx) const;

    /// @brief Refine a (connected block-shaped) patch of cells. Based on the patch, a Geometry<3,3> object is created and refined.
    ///
    /// @param [in] cells_per_dim            Number of (refined) cells in each direction that each parent cell should be refined to.
    /// @param [in] startIJK                 Cartesian triplet index where the patch starts.
    /// @param [in] endIJK                   Cartesian triplet index where the patch ends.
    ///                                      Last cell part of the lgr will be {endijk[0]-1, ... endIJK[2]-1}.
    ///
    /// @return refined_grid_ptr                   Shared pointer of CpGridData type, pointing at the refined_grid
    /// @return boundary_old_to_new_corners/faces  Corner/face indices on the patch-boundary associated with new-born-entity indices.
    /// @return parent_to_children_faces/cell      For each parent face/cell, we store its child-face/cell indices.
    ///                                            {parent face/cell index in coarse level, {indices of its children in refined level}}
    /// @return child_to_parent_faces/cells        {child index, parent index}
    std::tuple< std::shared_ptr<CpGridData>,
                const std::vector<std::array<int,2>>,                // boundary_old_to_new_corners
                const std::vector<std::tuple<int,std::vector<int>>>, // boundary_old_to_new_faces
                const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces
                const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_cell
                const std::vector<std::array<int,2>>,                // child_to_parent_faces
                const std::vector<std::array<int,2>>>                // child_to_parent_cells
    refinePatch(const std::array<int,3>& cells_per_dim, const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const;

    // @breif Compute center of an entity/element/cell in the Eclipse way:
    //        - Average of the 4 corners of the bottom face.
    //        - Average of the 4 corners of the top face.
    //        Return average of the previous computations.
    // @param [in]   int   Index of a cell.
    // @return            'eclipse centroid'
    std::array<double,3> computeEclCentroid(const int idx) const;

    // @breif Compute center of an entity/element/cell in the Eclipse way:
    //        - Average of the 4 corners of the bottom face.
    //        - Average of the 4 corners of the top face.
    //        Return average of the previous computations.
    // @param [in]   Entity<0>   Entity
    // @return                   'eclipse centroid'
    std::array<double,3> computeEclCentroid(const Entity<0>& elem) const;

    // Make unique boundary ids for all intersections.
    void computeUniqueBoundaryIds();

    /// Is the grid currently using unique boundary ids?
    /// \return true if each boundary intersection has a unique id
    ///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
    bool uniqueBoundaryIds() const
    {
        return use_unique_boundary_ids_;
    }

    /// Set whether we want to have unique boundary ids.
    /// \param uids if true, each boundary intersection will have a unique boundary id.
    void setUniqueBoundaryIds(bool uids)
    {
        use_unique_boundary_ids_ = uids;
        if (use_unique_boundary_ids_ && unique_boundary_ids_.empty()) {
            computeUniqueBoundaryIds();
        }
    }

    /// Return the internalized zcorn copy from the grid processing, if
    /// no cells were adjusted during the minpvprocessing this can be
    /// and empty vector.
    const std::vector<double>& zcornData() const {
        return zcorn;
    }


    /// Get the index set. This is the lead as well as th level index set.
    /// \return The index set.
    const IndexSet& indexSet() const
    {
        return *index_set_;
    }

    /// Get the local index set.
    const cpgrid::IdSet& localIdSet() const
    {
        return *local_id_set_;
    }

    /// Get the global index set.
    const cpgrid::LevelGlobalIdSet& globalIdSet() const
    {
        return *global_id_set_;
    }

    /// The logical cartesian size of the grid.
    /// This function is not part of the Dune grid interface,
    /// and should be used with caution.
    const std::array<int, 3>& logicalCartesianSize() const
    {
        return logical_cartesian_size_;
    }

    /// \brief Redistribute a global grid.
    ///
    /// The whole grid must be available on all processors.
    void distributeGlobalGrid(CpGrid& grid,
                              const CpGridData& view_data,
                              const std::vector<int>& cell_part);

    /// \brief communicate objects for all codims on a given level
    /// \param data The data handle describing the data. Has to adhere to the
    /// Dune::DataHandleIF interface.
    /// \param iftype The interface to use for the communication.
    /// \param dir The direction of the communication along the interface (forward or backward).
    template<class DataHandle>
    void communicate(DataHandle& data, InterfaceType iftype, CommunicationDirection dir);

    void computeCellPartitionType();

    void computePointPartitionType();

    void computeCommunicationInterfaces(int noexistingPoints);

    /// \brief The type of the mpi communicator.
    using MPICommunicator = CpGridDataTraits::MPICommunicator ;
    /// \brief The type of the collective communication.
    using Communication = CpGridDataTraits::Communication;
    using CollectiveCommunication = CpGridDataTraits::CollectiveCommunication;

    /// \brief The type of the set of the attributes
    using AttributeSet = CpGridDataTraits::AttributeSet;
#if HAVE_MPI
    /// \brief The type of the  Communicator.
    using Communicator = CpGridDataTraits::Communicator;

    /// \brief The type of the map describing communication interfaces.
    using InterfaceMap = CpGridDataTraits::InterfaceMap;

    /// \brief type of OwnerOverlap communication for cells
    using CommunicationType = CpGridDataTraits::CommunicationType;

    /// \brief The type of the parallel index set
    using  ParallelIndexSet = CpGridDataTraits::ParallelIndexSet;

    /// \brief The type of the remote indices information
    using RemoteIndices = CpGridDataTraits::RemoteIndices;

    /// \brief Get the owner-overlap-copy communication for cells
    ///
    /// Suitable e.g. for parallel linear algebra used by CCFV
    CommunicationType& cellCommunication()
    {
        return cell_comm_;
    }

    /// \brief Get the owner-overlap-copy communication for cells
    ///
    /// Suitable e.g. for parallel linear algebra used by CCFV
    const CommunicationType& cellCommunication() const
    {
        return cell_comm_;
    }

    ParallelIndexSet& cellIndexSet()
    {
        return cellCommunication().indexSet();
    }

    const ParallelIndexSet& cellIndexSet() const
    {
        return cellCommunication().indexSet();
    }

    RemoteIndices& cellRemoteIndices()
    {
        return cellCommunication().remoteIndices();
    }

    const RemoteIndices& cellRemoteIndices() const
    {
            return cellCommunication().remoteIndices();
    }
#endif

    /// \brief Get sorted active cell indices of numerical aquifer
    const std::vector<int>& sortedNumAquiferCells() const
    {
        return aquifer_cells_;
    }

private:

    /// \brief Adds entries to the parallel index set of the cells during grid construction
    void populateGlobalCellIndexSet();

#if HAVE_MPI

    /// \brief Gather data on a global grid representation.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    template<class DataHandle>
    void gatherData(DataHandle& data, CpGridData* global_view,
                    CpGridData* distributed_view);


    /// \brief Gather data specific to given codimension on a global grid representation.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    /// \tparam codim The codimension
    template<int codim, class DataHandle>
    void gatherCodimData(DataHandle& data, CpGridData* global_data,
                         CpGridData* distributed_data);

    /// \brief Scatter data from a global grid representation
    /// to a distributed representation of the same grid.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    template<class DataHandle>
    void scatterData(DataHandle& data, const CpGridData* global_data,
                     const CpGridData* distributed_data, const InterfaceMap& cell_inf,
                     const InterfaceMap& point_inf);

    /// \brief Scatter data specific to given codimension from a global grid representation
    /// to a distributed representation of the same grid.
    /// \param data A data handle for getting or setting the data
    /// \param global_view The view of the global grid (to gather the data on)
    /// \param distributed_view The view of the distributed grid.
    /// \tparam DataHandle The type of the data handle used.
    /// \tparam codim The codimension.
    template<int codim, class DataHandle>
    void scatterCodimData(DataHandle& data, CpGridData* global_data,
                          CpGridData* distributed_data);

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                          const Interface& interface);

    /// \brief Communicates data of a given codimension
    /// \tparam codim The codimension
    /// \tparam DataHandle The type of the data handle describing, gathering,
    ///  and gathering the data.
    /// \param DataHandle The data handle describing, gathering,
    ///  and gathering the data.
    /// \param dir The direction of the communication.
    /// \param interface The information about the communication interface
    template<int codim, class DataHandle>
    void communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                          const InterfaceMap& interface);

#endif

    void computeGeometry(const CpGrid& grid,
                         const DefaultGeometryPolicy&  globalGeometry,
                         const std::vector<int>& globalAquiferCells,
                         const OrientedEntityTable<0, 1>& globalCell2Faces,
                         DefaultGeometryPolicy& geometry,
                         std::vector<int>& aquiferCells,
                         const OrientedEntityTable<0, 1>& cell2Faces,
                         const std::vector< std::array<int,8> >& cell2Points);

    // Representing the topology
    /** @brief Container for lookup of the faces attached to each cell. */
    cpgrid::OrientedEntityTable<0, 1> cell_to_face_;
    /**
     * @brief Container for the lookup of attaching cells for each face.
     *
     * All faces have two neighbours except for those at the domain boundary.
     * @warn  Note that along the front partition there are invalid neighbours
     * marked with index std::numeric_limits<int>::max()
     */
    cpgrid::OrientedEntityTable<1, 0> face_to_cell_;
    /** @brief Container for the lookup of the points for each face. */
    Opm::SparseTable<int>             face_to_point_;
    /** @brief Vector that contains an arrays of the points of each cell*/
    std::vector< std::array<int,8> >       cell_to_point_;
    /** @brief The size of the underlying logical cartesian grid.
     *
     * In a Eclipse a cornerpoint grid has the same number of cells
     * in each pillar. Note that of these some may have no volume
     * and this be inactive.
     */
    std::array<int, 3>                logical_cartesian_size_{};
    /** @brief vector with the gobal cell index for each cell.
     *
     * Note the size of this container is determined by the
     * the number of cells present on the process and the content
     * by the mapping to the underlying global cartesian mesh..
     */
    std::vector<int>                  global_cell_;
    /** @brief The tag of the faces. */
    cpgrid::EntityVariable<enum face_tag, 1> face_tag_;
    /** @brief The geometries representing the grid. */
    cpgrid::DefaultGeometryPolicy geometry_;
    /** @brief The type of a point in the grid. */
    typedef FieldVector<double, 3> PointType;
    /** @brief The face normals of the grid. */
    cpgrid::SignedEntityVariable<PointType, 1> face_normals_;
    /** @brief The boundary ids. */
    cpgrid::EntityVariable<int, 1> unique_boundary_ids_;
    /** @brief The index set of the grid (level). */
    std::unique_ptr<cpgrid::IndexSet> index_set_;
    /** @brief The internal local id set (not exported). */
    std::shared_ptr<const cpgrid::IdSet> local_id_set_;
    /** @brief The global id set (used also as local id set). */
    std::shared_ptr<LevelGlobalIdSet> global_id_set_;
    /** @brief The indicator of the partition type of the entities */
    std::shared_ptr<PartitionTypeIndicator> partition_type_indicator_;
    /** Mark elements to be refined **/
    std::vector<int> mark_;
    /** Level of the current CpGridData (0 when it's "GLOBAL", 1,2,.. for LGRs). */
    int level_{0};
    /** Copy of (CpGrid object).data_ associated with the CpGridData object. */
    std::vector<std::shared_ptr<CpGridData>>* level_data_ptr_;
    // SUITABLE FOR ALL LEVELS EXCEPT FOR LEAFVIEW
    /** Map between level and leafview cell indices. Only cells (from that level) that appear in leafview count. -1 when the cell vanished.*/
    std::vector<int> level_to_leaf_cells_; // In entry 'level cell index', we store 'leafview cell index'
    /** Parent cells and their children. Entry is {-1, {}} when cell has no children.*/ // {level LGR, {child0, child1, ...}}
    std::vector<std::tuple<int,std::vector<int>>> parent_to_children_cells_;
    /** Amount of children cells per parent cell in each direction. */ // {# children in x-direction, ... y-, ... z-}
    std::array<int,3> cells_per_dim_;
    // SUITABLE ONLY FOR LEAFVIEW
    /** Relation between leafview and (possible different) level(s) cell indices. */ // {level, cell index in that level}
    std::vector<std::array<int,2>> leaf_to_level_cells_;
    /** Corner history. corner_history_[ corner index ] = {level where the corner was born, its index there }, {-1,-1} otherwise. */
    std::vector<std::array<int,2>> corner_history_;
    // SUITABLE FOR ALL LEVELS INCLUDING LEAFVIEW
    /** Child cells and their parents. Entry is {-1,-1} when cell has no father. */ // {level parent cell, parent cell index}
    std::vector<std::array<int,2>> child_to_parent_cells_;
    /** Level-grid or Leaf-grid cell to parent cell and refined-cell-in-parent-cell index (number between zero and total amount
        of children per parent (cells_per_dim[0]_*cells_per_dim_[1]*cells_per_dim_[2])). Entry is -1 when cell has no father. */
    std::vector<int> cell_to_idxInParentCell_;



    /// \brief Object for collective communication operations.
    Communication ccobj_;

    // Boundary information (optional).
    bool use_unique_boundary_ids_;

    /// This vector contains zcorn values from the initialization
    /// process where a CpGrid instance has been created from
    /// cornerpoint input zcorn and coord. During the initialization
    /// the zcorn values will typically be modified, and we retain a
    /// copy here to be able to create an EclipseGrid for output.
    std::vector<double> zcorn;

    /// \brief Sorted vector of aquifer cell indices.
    std::vector<int> aquifer_cells_;

#if HAVE_MPI

    /// \brief OwnerOverlap communication for cells
    CommunicationType cell_comm_;

    /// \brief Communication interface for the cells.
    std::tuple<Interface,Interface,Interface,Interface,Interface> cell_interfaces_;
    /*
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    /// \brief Interface from interior and border to interior and border for the faces.
    std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>
    face_interfaces_;
    */
    /// \brief Interface from interior and border to interior and border for the faces.
    std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>
    point_interfaces_;

#endif

    // Return the geometry vector corresponding to the given codim.
    template <int codim>
    const EntityVariable<Geometry<3 - codim, 3>, codim>& geomVector() const
    {
        return geometry_.geomVector<codim>();
    }

    friend class Dune::CpGrid;
    template<int> friend class Entity;
    template<int> friend class EntityRep;
    friend class Intersection;
    friend class PartitionTypeIndicator;
};



#if HAVE_MPI

namespace
{
/// \brief Get a value from a tuple according to the interface type.
/// \tparam T The type of the values in the tuple.
/// \param iftype The interface type.
/// \param interfaces A tuple with the values order by interface type.
template<class T>
T& getInterface(InterfaceType iftype,
                std::tuple<T,T,T,T,T>& interfaces)
{
    switch(iftype)
    {
    case 0:
        return std::get<0>(interfaces);
    case 1:
        return std::get<1>(interfaces);
    case 2:
        return std::get<2>(interfaces);
    case 3:
        return std::get<3>(interfaces);
    case 4:
        return std::get<4>(interfaces);
    }
    OPM_THROW(std::runtime_error, "Invalid Interface type was used during communication");
}

} // end unnamed namespace

template<int codim, class DataHandle>
void CpGridData::communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data, CommunicationDirection dir,
                                  const Interface& interface)
{
    this->template communicateCodim<codim>(data, dir, interface.interfaces());
}

template<int codim, class DataHandle>
void CpGridData::communicateCodim(Entity2IndexDataHandle<DataHandle, codim>& data_wrapper, CommunicationDirection dir,
                                  const InterfaceMap& interface)
{
    Communicator comm(ccobj_, interface);

    if(dir==ForwardCommunication)
        comm.forward(data_wrapper);
    else
        comm.backward(data_wrapper);
}
#endif

template<class DataHandle>
void CpGridData::communicate(DataHandle& data, InterfaceType iftype,
                             CommunicationDirection dir)
{
#if HAVE_MPI
    if(data.contains(3,0))
    {
        Entity2IndexDataHandle<DataHandle, 0> data_wrapper(*this, data);
        communicateCodim<0>(data_wrapper, dir, getInterface(iftype, cell_interfaces_));
    }
    if(data.contains(3,3))
    {
        Entity2IndexDataHandle<DataHandle, 3> data_wrapper(*this, data);
        communicateCodim<3>(data_wrapper, dir, getInterface(iftype, point_interfaces_));
    }
#else
    // Suppress warnings for unused arguments.
    (void) data;
    (void) iftype;
    (void) dir;
#endif
}
}}

#if HAVE_MPI
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/Indexsets.hpp>

namespace Dune {
namespace cpgrid {

namespace mover
{
template<class T>
class MoveBuffer
{
    friend class Dune::cpgrid::CpGridData;
public:
    void read(T& data)
    {
        data=buffer_[index_++];
    }
    void write(const T& data)
    {
        buffer_[index_++]=data;
    }
    void reset()
    {
        index_=0;
    }
    void resize(std::size_t size)
    {
        buffer_.resize(size);
        index_=0;
    }
private:
    std::vector<T> buffer_;
    typename std::vector<T>::size_type index_;
};
template<class DataHandle,int codim>
struct Mover
{
};

template<class DataHandle>
struct BaseMover
{
    explicit BaseMover(DataHandle& data)
    : data_(data)
    {}
    template<class E>
    void moveData(const E& from, const E& to)
    {
        std::size_t size=data_.size(from);
        buffer.resize(size);
        data_.gather(buffer, from);
        buffer.reset();
        data_.scatter(buffer, to, size);
    }
    DataHandle& data_;
    MoveBuffer<typename DataHandle::DataType> buffer;
};


template<class DataHandle>
struct Mover<DataHandle,0> : public BaseMover<DataHandle>
{
    Mover(DataHandle& data, CpGridData* gatherView,
          CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        Entity<0> from_entity=Entity<0>(*gatherView_, from_cell_index, true);
        Entity<0> to_entity=Entity<0>(*scatterView_, to_cell_index, true);
        this->moveData(from_entity, to_entity);
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,1> : public BaseMover<DataHandle>
{
    Mover(DataHandle& data, CpGridData* gatherView,
          CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}

    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        typedef typename OrientedEntityTable<0,1>::row_type row_type;
        EntityRep<0> from_cell=EntityRep<0>(from_cell_index, true);
        EntityRep<0> to_cell=EntityRep<0>(to_cell_index, true);
        const OrientedEntityTable<0,1>& table = gatherView_->cell_to_face_;
        row_type from_faces=table.operator[](from_cell);
        row_type to_faces=scatterView_->cell_to_face_[to_cell];

        for(int i=0; i<from_faces.size(); ++i)
            this->moveData(from_faces[i], to_faces[i]);
    }
    CpGridData *gatherView_;
    CpGridData *scatterView_;
};

template<class DataHandle>
struct Mover<DataHandle,3> : public BaseMover<DataHandle>
{
    Mover(DataHandle& data, CpGridData* gatherView,
          CpGridData* scatterView)
    : BaseMover<DataHandle>(data), gatherView_(gatherView), scatterView_(scatterView)
    {}
    void operator()(std::size_t from_cell_index,std::size_t to_cell_index)
    {
        const std::array<int,8>& from_cell_points=
            gatherView_->cell_to_point_[from_cell_index];
        const std::array<int,8>& to_cell_points=
            scatterView_->cell_to_point_[to_cell_index];
        for(std::size_t i=0; i<8; ++i)
        {
            this->moveData(Entity<3>(*gatherView_, from_cell_points[i], true),
                           Entity<3>(*scatterView_, to_cell_points[i], true));
        }
    }
    CpGridData* gatherView_;
    CpGridData* scatterView_;
};

} // end mover namespace

template<class DataHandle>
void CpGridData::scatterData(DataHandle& data, const CpGridData* global_data,
                             const CpGridData* distributed_data, const InterfaceMap& cell_inf,
                             const InterfaceMap& point_inf)
{
#if HAVE_MPI
    if(data.contains(3,0))
    {
        Entity2IndexDataHandle<DataHandle, 0> data_wrapper(*global_data, *distributed_data, data);
        communicateCodim<0>(data_wrapper, ForwardCommunication, cell_inf);
    }
    if(data.contains(3,3))
    {
        Entity2IndexDataHandle<DataHandle, 3> data_wrapper(*global_data, *distributed_data, data);
        communicateCodim<3>(data_wrapper, ForwardCommunication, point_inf);
    }
#endif
}

template<int codim, class DataHandle>
void CpGridData::scatterCodimData(DataHandle& data, CpGridData* global_data,
                                  CpGridData* distributed_data)
{
    CpGridData *gather_view, *scatter_view;
    gather_view=global_data;
    scatter_view=distributed_data;

    mover::Mover<DataHandle,codim> mover(data, gather_view, scatter_view);


    for(auto index=distributed_data->cellIndexSet().begin(),
            end = distributed_data->cellIndexSet().end();
        index!=end; ++index)
    {
        std::size_t from=index->global();
        std::size_t to=index->local();
        mover(from,to);
    }
}

namespace
{

template<int codim, class T, class F>
void visitInterior(CpGridData& distributed_data, T begin, T endit, F& func)
{
    for(T it=begin; it!=endit; ++it)
    {
        Entity<codim> entity(distributed_data, it-begin, true);
        PartitionType pt = entity.partitionType();
        if(pt==Dune::InteriorEntity)
        {
            func(*it, entity);
        }
    }
}

template<class DataHandle>
struct GlobalIndexSizeGatherer
{
    GlobalIndexSizeGatherer(DataHandle& data_,
                            std::vector<int>& ownedGlobalIndices_,
                            std::vector<int>& ownedSizes_)
        : data(data_), ownedGlobalIndices(ownedGlobalIndices_), ownedSizes(ownedSizes_)
    {}

    template<class T, class E>
    void operator()(T& i, E& entity)
    {
            ownedGlobalIndices.push_back(i);
            ownedSizes.push_back(data.size(entity));
    }
    DataHandle& data;
    std::vector<int>& ownedGlobalIndices;
    std::vector<int>& ownedSizes;
};

template<class DataHandle>
struct DataGatherer
{
    DataGatherer(mover::MoveBuffer<typename DataHandle::DataType>& buffer_,
                 DataHandle& data_)
        : buffer(buffer_), data(data_)
    {}

    template<class T, class E>
    void operator()(T& /* it */, E& entity)
    {
        data.gather(buffer, entity);
    }
    mover::MoveBuffer<typename DataHandle::DataType>& buffer;
    DataHandle& data;
};

}

template<class DataHandle>
void CpGridData::gatherData(DataHandle& data, CpGridData* global_data,
                            CpGridData* distributed_data)
{
#if HAVE_MPI
    if(data.contains(3,0))
       gatherCodimData<0>(data, global_data, distributed_data);
    if(data.contains(3,3))
       gatherCodimData<3>(data, global_data, distributed_data);
#endif
}

template<int codim, class DataHandle>
void CpGridData::gatherCodimData(DataHandle& data, CpGridData* global_data,
                                 CpGridData* distributed_data)
{
#if HAVE_MPI
    // Get the mapping to global index from  the global id set
    const std::vector<int>& mapping =
        distributed_data->global_id_set_->getMapping<codim>();

    // Get the global indices and data size for the entities whose data is
    // to be sent, i.e. the ones that we own.
    std::vector<int>         owned_global_indices;
    std::vector<int> owned_sizes;
    owned_global_indices.reserve(mapping.size());
    owned_sizes.reserve(mapping.size());

    GlobalIndexSizeGatherer<DataHandle> gisg(data, owned_global_indices, owned_sizes);
    visitInterior<codim>(*distributed_data, mapping.begin(), mapping.end(), gisg);

    // communicate the number of indices that each processor sends
    int no_indices=owned_sizes.size();
    // We will take the address of the first elemet for MPI_Allgather below.
    // Make sure the containers have such an element.
    if ( owned_global_indices.empty() )
        owned_global_indices.resize(1);
    if ( owned_sizes.empty() )
        owned_sizes.resize(1);
    std::vector<int> no_indices_to_recv(distributed_data->ccobj_.size());
    distributed_data->ccobj_.allgather(&no_indices, 1, &(no_indices_to_recv[0]));
    // compute size of the vector capable for receiving all indices
    // and allgather the global indices and the sizes.
    // calculate displacements
    std::vector<int> displ(distributed_data->ccobj_.size()+1, 0);
    std::transform(displ.begin(), displ.end()-1, no_indices_to_recv.begin(), displ.begin()+1,
                   std::plus<int>());
    int global_size=displ[displ.size()-1];//+no_indices_to_recv[displ.size()-1];
    std::vector<int>         global_indices(global_size);
    std::vector<int> global_sizes(global_size);
    MPI_Allgatherv(&(owned_global_indices[0]), no_indices, MPITraits<int>::getType(),
                   &(global_indices[0]), &(no_indices_to_recv[0]), &(displ[0]),
                   MPITraits<int>::getType(),
                   distributed_data->ccobj_);
    MPI_Allgatherv(&(owned_sizes[0]), no_indices, MPITraits<int>::getType(),
                   &(global_sizes[0]), &(no_indices_to_recv[0]), &(displ[0]),
                   MPITraits<int>::getType(),
                   distributed_data->ccobj_);
    std::vector<int>().swap(owned_global_indices); // free data for reuse.
    // Compute the number of data items to send
    std::vector<int> no_data_send(distributed_data->ccobj_.size());
    for(typename std::vector<int>::iterator begin=no_data_send.begin(),
            i=begin, end=no_data_send.end(); i!=end; ++i)
        *i = std::accumulate(global_sizes.begin()+displ[i-begin],
                            global_sizes.begin()+displ[i-begin+1], std::size_t());
    // free at least some memory that can be reused.
    std::vector<int>().swap(owned_sizes);
    // compute the displacements for receiving with allgatherv
    displ[0]=0;
    std::transform(displ.begin(), displ.end()-1, no_data_send.begin(), displ.begin()+1,
                   std::plus<std::size_t>());
    // Compute the number of data items we will receive
    int no_data_recv = displ[displ.size()-1];//+global_sizes[displ.size()-1];

    // Collect the data to send, gather it
    mover::MoveBuffer<typename DataHandle::DataType> local_data_buffer, global_data_buffer;
    if ( no_data_send[distributed_data->ccobj_.rank()] )
    {
        local_data_buffer.resize(no_data_send[distributed_data->ccobj_.rank()]);
    }
    else
    {
        local_data_buffer.resize(1);
    }
    global_data_buffer.resize(no_data_recv);

    DataGatherer<DataHandle> gatherer(local_data_buffer, data);
    visitInterior<codim>(*distributed_data, mapping.begin(), mapping.end(), gatherer);
    MPI_Allgatherv(&(local_data_buffer.buffer_[0]), no_data_send[distributed_data->ccobj_.rank()],
                   MPITraits<typename DataHandle::DataType>::getType(),
                   &(global_data_buffer.buffer_[0]), &(no_data_send[0]), &(displ[0]),
                   MPITraits<typename DataHandle::DataType>::getType(),
                   distributed_data->ccobj_);
    Entity2IndexDataHandle<DataHandle, codim> edata(*global_data, data);
    int offset=0;
    for(int i=0; i< codim; ++i)
        offset+=global_data->size(i);

    typename std::vector<int>::const_iterator s=global_sizes.begin();
    for(typename std::vector<int>::const_iterator i=global_indices.begin(),
            end=global_indices.end();
        i!=end; ++s, ++i)
    {
        edata.scatter(global_data_buffer, *i-offset, *s);
    }
#endif
}

} // end namespace cpgrid
} // end namespace Dune

#endif

#endif
