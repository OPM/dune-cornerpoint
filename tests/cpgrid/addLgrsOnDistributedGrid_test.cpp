//===========================================================================
//
// File: addLgrsOnDistributedGrid_test.cpp
//
// Created: Monday 5 August 2024 14:44
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2024 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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
#include "config.h"

#define BOOST_TEST_MODULE LGRsOnDistributedGridTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/LookUpData.hh>
#include <opm/grid/common/CommunicationUtils.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);



#define CHECK_COORDINATES(c1, c2)                                       \
    for (int c = 0; c < 3; c++) {                                       \
        BOOST_TEST(c1[c] == c2[c], boost::test_tools::tolerance(1e-12)); \
    }

void refinePatch_and_check(Dune::CpGrid& coarse_grid,
                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                           const std::vector<std::array<int,3>>& startIJK_vec,
                           [[maybe_unused]] const std::vector<std::array<int,3>>& endIJK_vec,
                           const std::vector<std::string>& lgr_name_vec)
{
    auto& data = coarse_grid.currentData(); // what data current_view_data_ is pointing at (data_ or distributed_data_)

    BOOST_CHECK(data.size() == startIJK_vec.size() + 2);
    BOOST_CHECK( data[0]->child_to_parent_cells_.empty());
    BOOST_CHECK(coarse_grid.getLgrNameToLevel().at("GLOBAL") == 0);

    for (long unsigned int level = 1; level < startIJK_vec.size() +1; ++level) // only 1 when there is only 1 patch
    {
        BOOST_CHECK( (*data[level]).parent_to_children_cells_.empty());
        BOOST_CHECK(coarse_grid.getLgrNameToLevel().at(lgr_name_vec[level-1]) == static_cast<int>(level));

        // GLOBAL grid
        for (int cell = 0; cell <  data[0]-> size(0); ++cell)
        {
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[0], cell, true);
            BOOST_CHECK( entity.hasFather() == false);
            BOOST_CHECK_THROW(entity.father(), std::logic_error);
            BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
            BOOST_CHECK( entity.getOrigin() ==  entity);
            BOOST_CHECK( entity.getOrigin().level() == 0);
            auto it = entity.hbegin(coarse_grid.maxLevel());
            auto endIt = entity.hend(coarse_grid.maxLevel());
            const auto& [lgr, childrenList] = (*data[0]).parent_to_children_cells_[cell];
            if (entity.isLeaf()){ // In particular, cell has no children/is not a father.
                BOOST_CHECK_EQUAL(lgr, -1);
                BOOST_CHECK(childrenList.empty());
                BOOST_CHECK( entity.isLeaf() == true);
                // If it == endIt, then entity.isLeaf() true (when dristibuted_data_ is empty)
                BOOST_CHECK( it == endIt);
            }
            else{
                BOOST_CHECK(lgr != -1);
                BOOST_CHECK(childrenList.size() > 1);
                // Auxiliary int to check amount of children
                double referenceElemOneParent_volume = 0.;
                std::array<double,3> referenceElem_entity_center = {0.,0.,0.}; // Expected {.5,.5,.5}
                for (const auto& child : childrenList) {
                    BOOST_CHECK( child != -1);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][0] == 0);
                    BOOST_CHECK( data[lgr]-> child_to_parent_cells_[child][1] == cell);

                    const auto& childElem =  Dune::cpgrid::Entity<0>(*data[lgr], child, true);
                    BOOST_CHECK(childElem.hasFather() == true);
                    BOOST_CHECK(childElem.level() == lgr);
                    referenceElemOneParent_volume += childElem.geometryInFather().volume();
                    BOOST_CHECK( childElem.geometryInFather().volume() > 0 );
                    for (int c = 0; c < 3; ++c)  {
                        referenceElem_entity_center[c] += (childElem.geometryInFather().center())[c];
                    }
                }
                BOOST_CHECK_EQUAL( entity.isLeaf(), false); // parent cells do not appear in the LeafView
                // Auxiliary int to check hierarchic iterator functionality
                double referenceElemOneParent_volume_it = 0.;
                std::array<double,3> referenceElem_entity_center_it = {0.,0.,0.}; // Expected {.5,.5,.5}
                // If it != endIt, then entity.isLeaf() false (when dristibuted_data_ is empty)
                BOOST_CHECK( it != endIt );
                for (; it != endIt; ++it)
                {
                    // Do something with the son available through it->
                    BOOST_CHECK(it ->hasFather() == true);
                    BOOST_CHECK(it ->level() == lgr);
                    referenceElemOneParent_volume_it += it-> geometryInFather().volume();
                    for (int c = 0; c < 3; ++c)
                    {
                        referenceElem_entity_center_it[c] += (it-> geometryInFather().center())[c];
                    }
                }
                for (int c = 0; c < 3; ++c)
                {
                    referenceElem_entity_center[c]
                        /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
                    referenceElem_entity_center_it[c]
                        /= cells_per_dim_vec[lgr-1][0]*cells_per_dim_vec[lgr-1][1]*cells_per_dim_vec[lgr-1][2];
                }
                BOOST_CHECK_CLOSE(referenceElemOneParent_volume, 1, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[0], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[1], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center[2], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElemOneParent_volume_it, 1, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[0], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[1], .5, 1e-13);
                BOOST_CHECK_CLOSE(referenceElem_entity_center_it[2], .5, 1e-13);
            }
            BOOST_CHECK( entity.level() == 0);
        }

        if (!(data[level] -> global_cell_.empty()))
        {
            auto itMin = std::min_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
            auto itMax = std::max_element((data[level] -> global_cell_).begin(),  (data[level] -> global_cell_).end());
            BOOST_CHECK( *itMin >= 0); // An LGR can have cells distributed across different processes, so the minimum cell global id may not be zero in all processes.
            const auto& maxCartesianIdxLevel = data[level]->logical_cartesian_size_[0]*data[level]->logical_cartesian_size_[1]* data[level]->logical_cartesian_size_[2];
            BOOST_CHECK( *itMax < maxCartesianIdxLevel);
        }

        // LGRs
        for (int cell = 0; cell <  data[level]-> size(0); ++cell)
        {
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[level], cell, true);
            BOOST_CHECK( entity.hasFather() == true);
            BOOST_CHECK( entity.getOrigin() ==  entity.father());
            BOOST_CHECK( entity.getOrigin().level() == 0);
            BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                              1./(cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2]), 1e-6);
            BOOST_CHECK(entity.father().level() == 0);
            const auto& child_to_parent = (*data[level]).child_to_parent_cells_[cell];
            BOOST_CHECK_EQUAL( child_to_parent[0] == -1, false);
            BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
            BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
            BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[child_to_parent[1]]) == entity.level());
            BOOST_CHECK_EQUAL((std::find(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).begin(),
                                         std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end(),
                                         entity.index()) ==
                               std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end()) , false);
            // Check amount of children cells of the parent cell
            BOOST_CHECK_EQUAL(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).size(),
                              cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2]);
            BOOST_CHECK( entity.level() == static_cast<int>(level));
            BOOST_CHECK( entity.isLeaf() == true);
            auto it = entity.hbegin(coarse_grid.maxLevel());
            auto endIt = entity.hend(coarse_grid.maxLevel());
            // If entity.isLeaf(), then it == endIt (when dristibuted_data_ is empty)
            BOOST_CHECK( it == endIt);
        }

        // LeafView faces
        for (int face = 0; face <  data[startIJK_vec.size()+1]-> face_to_cell_.size(); ++face)
        {
            const auto& faceToPoint =  (*data[startIJK_vec.size() +1]).face_to_point_[face];
            BOOST_CHECK(faceToPoint.size() == 4);
            for (int i = 0; i < 4; ++i) {
                BOOST_CHECK((*data[startIJK_vec.size() +1]).face_to_point_[face][i] != -1);
            }

            Dune::cpgrid::EntityRep<1> faceEntity(face, true);
            BOOST_CHECK((*data[startIJK_vec.size() +1]).face_to_cell_[faceEntity].size() < 3);
        }

        auto itMin = std::min_element((data.back() -> global_cell_).begin(),  (data.back()-> global_cell_).end());
        auto itMax = std::max_element((data.back() -> global_cell_).begin(),  (data.back() -> global_cell_).end());
        BOOST_CHECK( *itMin >= 0);
        const auto& maxCartesianIdx = coarse_grid.logicalCartesianSize()[0]*coarse_grid.logicalCartesianSize()[1]*coarse_grid.logicalCartesianSize()[2];
        BOOST_CHECK( *itMax < maxCartesianIdx);

        // LeafView
        for (int cell = 0; cell <  data[startIJK_vec.size()+1]-> size(0); ++cell)
        {
            BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_point_[cell].size() == 8);
            for (int i = 0; i < 8; ++i)
            {
                BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_point_[cell][i] != -1);
            }
            Dune::cpgrid::Entity<0> entity = Dune::cpgrid::Entity<0>(*data[startIJK_vec.size()+1], cell, true);
            for (int i = 0; i < data[startIJK_vec.size()+1] -> cell_to_face_[entity].size(); ++i)
            {
                BOOST_CHECK( data[startIJK_vec.size()+1] -> cell_to_face_[entity][i].index() != -1);
            }
            const auto& child_to_parent = (*data[startIJK_vec.size()+1]).child_to_parent_cells_[cell];
            const auto& level_cellIdx = (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()];
            auto it = entity.hbegin(coarse_grid.maxLevel());
            auto endIt = entity.hend(coarse_grid.maxLevel());
            BOOST_CHECK(entity.isLeaf());
            BOOST_CHECK(it == endIt);
            if (entity.hasFather()){
                BOOST_CHECK_CLOSE(entity.geometryInFather().volume(),
                                  1./(cells_per_dim_vec[entity.level()-1][0]
                                      *cells_per_dim_vec[entity.level()-1][1]
                                      *cells_per_dim_vec[entity.level()-1][2]), 1e-6);
                BOOST_CHECK(entity.father().level() == 0);
                BOOST_CHECK(child_to_parent[0] != -1);
                BOOST_CHECK_EQUAL( child_to_parent[0] == 0, true);
                BOOST_CHECK_EQUAL( child_to_parent[1], entity.father().index());
                BOOST_CHECK( entity.father() == entity.getOrigin());
                BOOST_CHECK( entity.getOrigin().level() == 0);
                BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[child_to_parent[1]]) == entity.level());
                BOOST_CHECK_EQUAL((std::find(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).begin(),
                                             std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end(),
                                             level_cellIdx[1]) ==
                                   std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).end()) , false);
                // Check amount of children cells of the parent cell
                BOOST_CHECK_EQUAL(std::get<1>((*data[0]).parent_to_children_cells_[child_to_parent[1]]).size(),
                                  cells_per_dim_vec[entity.level()-1][0]*
                                  cells_per_dim_vec[entity.level()-1][1]*cells_per_dim_vec[entity.level()-1][2]);
                BOOST_CHECK( entity.father().isLeaf() == false);
                BOOST_CHECK( (entity.level() > 0) || (entity.level() < static_cast<int>(startIJK_vec.size()) +1));
                BOOST_CHECK( level_cellIdx[0] == entity.level());
            }
            else{
                BOOST_CHECK_THROW(entity.father(), std::logic_error);
                BOOST_CHECK_THROW(entity.geometryInFather(), std::logic_error);
                BOOST_CHECK_EQUAL(child_to_parent[0], -1);
                BOOST_CHECK_EQUAL(child_to_parent[1], -1);
                BOOST_CHECK( level_cellIdx[0] == 0);
                BOOST_CHECK( std::get<0>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]) == -1);
                BOOST_CHECK( std::get<1>((*data[0]).parent_to_children_cells_[level_cellIdx[1]]).empty());
                BOOST_CHECK( entity.level() == 0);
                // Get index of the cell in level 0
                const auto& entityOldIdx =  (*data[startIJK_vec.size()+1]).leaf_to_level_cells_[entity.index()][1];
                BOOST_CHECK( entity.getOrigin().index() == entityOldIdx);
                BOOST_CHECK( entity.getOrigin().level() == 0);
            } // end else
        }
    } // end-level-for-loop

    BOOST_CHECK( static_cast<int>(startIJK_vec.size()) == coarse_grid.maxLevel());
    BOOST_CHECK( (*data[data.size()-1]).parent_to_children_cells_.empty());

    auto it_min = std::min_element((data.back() -> global_cell_).begin(),  (data.back()-> global_cell_).end());
    auto it_max = std::max_element((data.back() -> global_cell_).begin(),  (data.back() -> global_cell_).end());
    auto it_min_level_zero = std::min_element((data.front() -> global_cell_).begin(),  (data.front() -> global_cell_).end());
    auto it_max_level_zero = std::max_element((data.front() -> global_cell_).begin(),  (data.front() -> global_cell_).end());
    BOOST_CHECK_EQUAL( *it_min, *it_min_level_zero);
    BOOST_CHECK_EQUAL( *it_max, *it_max_level_zero);

    for (long unsigned int l = 0; l < startIJK_vec.size() +1; ++l) // level 0,1,2,... , last patch
    {
        const auto& view = coarse_grid.levelGridView(l);
        for (const auto& element: elements(view)){
            BOOST_CHECK_EQUAL(element.level(), l);
        }
    }

    std::vector<int> leaf_to_parent_cell; // To store parent cell index, when leaf cell has a parent. Empty entry otherwise.
    leaf_to_parent_cell.resize(data[startIJK_vec.size()+1]-> size(0)); // Correct size.

    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(coarse_grid.leafGridView(), Dune::mcmgElementLayout());
    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> level0Mapper(coarse_grid.levelGridView(0), Dune::mcmgElementLayout());

    for (const auto& element: elements(coarse_grid.leafGridView())){
        BOOST_CHECK( ((element.level() >= 0) || (element.level() < static_cast<int>(startIJK_vec.size()) +1)));
        if (element.hasFather()) { // leaf_cell has a father!
            leaf_to_parent_cell[leafMapper.index(element)] = level0Mapper.index(element.father());
            const auto& parent_id = data[0]->localIdSet().id(element.father());
            BOOST_CHECK(element.index() == leafMapper.index(element));
            BOOST_CHECK(element.father().index() == leaf_to_parent_cell[element.index()]);
            BOOST_CHECK(element.father().index() == parent_id);
            BOOST_CHECK(element.father().index() == level0Mapper.index(element.father()));
        }
    }

    // Ids on the leaf grid view (local id and global id coincide ON THE LEAF GRID VIEW, might differ in level grids)
    std::set<int> allIds_set;
    std::vector<int> allIds_vec;
    allIds_vec.reserve(data.back()->size(0) + data.back()->size(3));
    for (const auto& element: elements(coarse_grid.leafGridView())){
        const auto& localId = data.back()->localIdSet().id(element);
        const auto& globalId = data.back()->globalIdSet().id(element);
        if (coarse_grid.comm().size()==0) {
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(localId, globalId);
        }
        allIds_set.insert(localId);
        allIds_vec.push_back(localId);
        // Check that the global_id_set_ptr_ has the correct id (id from the level where the entity was born).
        BOOST_CHECK_EQUAL( coarse_grid.globalIdSet().id(element), data[element.level()]->globalIdSet().id(element.getEquivLevelElem()));
    }
    // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding cell ids.
    BOOST_CHECK( allIds_set.size() == allIds_vec.size());

    for (const auto& point : vertices(coarse_grid.leafGridView())){
        const auto& localId = data.back()->localIdSet().id(point);
        const auto& globalId = data.back()->globalIdSet().id(point);
        if (coarse_grid.comm().size()==0) {
            // In serial run, local and global id coincide:
            BOOST_CHECK_EQUAL(localId, globalId);
        }
        allIds_set.insert(localId);
        allIds_vec.push_back(localId);
    }
    // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_) after adding point ids.
    BOOST_CHECK( allIds_set.size() == allIds_vec.size());

    // -------------------------------------------------------------

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
        // Check global id is not duplicated for interior cells in each refined level grid.
        std::vector<int> interior_cell_global_ids;
        int local_interior_cell_count = 0;
        interior_cell_global_ids.reserve(data[level]->size(0));

        for (const auto& element: elements(coarse_grid.levelGridView(level), Dune::Partitions::interior)){
            interior_cell_global_ids.push_back(data[level]->globalIdSet().id(element));
            ++local_interior_cell_count;
        }

        auto global_level_interior_cell_count = coarse_grid.comm().sum(local_interior_cell_count);
        auto [all_level_interior_cell_global_ids, displ] = Opm::allGatherv(interior_cell_global_ids, coarse_grid.comm());
        const std::set<int> all_level_interior_cell_global_ids_set(all_level_interior_cell_global_ids.begin(), all_level_interior_cell_global_ids.end());
        BOOST_CHECK( all_level_interior_cell_global_ids.size() == all_level_interior_cell_global_ids_set.size() );
        BOOST_CHECK( global_level_interior_cell_count == static_cast<int>(all_level_interior_cell_global_ids_set.size()) );
    }

    // Check global id is not duplicated for interior cells
    std::vector<int> localInteriorCellIds_vec;
    localInteriorCellIds_vec.reserve(data.back()->size(0)); // more than actually needed since only care about interior cells
    int local_interior_cells_count = 0;
    for (const auto& element: elements(coarse_grid.leafGridView())) {
        const auto& elemPartitionType = element.getEquivLevelElem().partitionType();
        if ( elemPartitionType == Dune::InteriorEntity) {
            localInteriorCellIds_vec.push_back(data.back()->globalIdSet().id(element));
            ++local_interior_cells_count;
        }
    }
    auto global_cells_count  = coarse_grid.comm().sum(local_interior_cells_count);
    auto [allGlobalIds_cells, displ] = Opm::allGatherv(localInteriorCellIds_vec, coarse_grid.comm());

    const std::set<int> allGlobalIds_cells_set(allGlobalIds_cells.begin(), allGlobalIds_cells.end());
    BOOST_CHECK( static_cast<int>(allGlobalIds_cells.size()) == global_cells_count);
    BOOST_CHECK( allGlobalIds_cells.size() == allGlobalIds_cells_set.size() );

    /** Uniqueness of point global ids cannot be checked in general since current code sets overlap layer size equal to 1,
        which in particular means that cells that share corners or edges (and not faces) with interior cells are not considered/
        seen by the process. Therefore, depending how the LGRs are distributed, there may be "multiple ids" for the same points.
        For some cases, achieving unique point global ids is possible by calling loadBalance( parts, false, true) where parts determines
        the rank of each cell, false represents ownerFisrt, true represents addCornerCells. The last defaulted argument is
        overlapLayerSize = 1.
    */

    // Local/Global id sets for level grids (level 0, 1, ..., maxLevel). For level grids, local might differ from global id.
    for (int level = 0; level < coarse_grid.maxLevel() +1; ++level)
    {
        std::set<int> levelIds_set;
        std::vector<int> levelIds_vec;
        levelIds_vec.reserve(data[level]->size(0) + data[level]->size(3));
        for (const auto& element: elements(coarse_grid.levelGridView(level))){
            const auto& localId = data[level]->localIdSet().id(element);
            // In parallel run, local and global id might not coincide.
            levelIds_set.insert(localId);
            levelIds_vec.push_back(localId);
            // Search in the leaf grid view elements for the element with the same id, if it exists.
            if (auto itIsLeaf = std::find_if( elements(coarse_grid.leafGridView()).begin(),
                                              elements(coarse_grid.leafGridView()).end(),
                                              [localId, data](const Dune::cpgrid::Entity<0>& leafElem)
                                              { return (localId == data.back()->localIdSet().id(leafElem)); });
                itIsLeaf != elements(coarse_grid.leafGridView()).end()) {
                BOOST_CHECK( itIsLeaf->getEquivLevelElem() == element);
            }
            if (element.isLeaf()) { // Check that the id of a cell not involved in any further refinement appears on the IdSet of the leaf grid view.
                BOOST_CHECK( std::find(allIds_set.begin(), allIds_set.end(), localId) != allIds_set.end());
            }
            else { // Check that the id of a cell that vanished during refinement does not appear on the IdSet of the leaf grid view.
                BOOST_CHECK( std::find(allIds_set.begin(), allIds_set.end(), localId) == allIds_set.end());
            }
            const auto& idx = data[level]->indexSet().index(element);
            BOOST_CHECK_EQUAL(idx, element.index());
        }

        for (const auto& point : vertices(coarse_grid.levelGridView(level))) {
            const auto& localId = data[level]->localIdSet().id(point);
            levelIds_set.insert(localId);
            levelIds_vec.push_back(localId);
            // Search in the leaf grid view elements for the element with the same id, if it exists.
            if (auto itIsLeaf = std::find_if( vertices(coarse_grid.leafGridView()).begin(),
                                              vertices(coarse_grid.leafGridView()).end(),
                                              [localId, data](const Dune::cpgrid::Entity<3>& leafPoint)
                                              { return (localId == data.back()->localIdSet().id(leafPoint)); });
                itIsLeaf != vertices(coarse_grid.leafGridView()).end()) {
                BOOST_CHECK( (*itIsLeaf).geometry().center() == point.geometry().center() );
            }
        }
        // Check injectivity of the map local_id_set_ (and, indirectly, global_id_set_)
        BOOST_CHECK( levelIds_set.size() == levelIds_vec.size());
    }
}

BOOST_AUTO_TEST_CASE(threeLgrs)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {10,8,8};
    grid.createCartesian(grid_dim, cell_sizes);

    // Fully interior LGRs
    std::vector<int> parts(640);
    for (int k = 0; k < 8; ++k) {
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 10; ++i)
            {
                const auto& elemIdx = (k*80) + (j*10) + i;
                if (i<5) {
                    if(j<4) {
                        parts[elemIdx] = 0;
                    }
                    else {
                        parts[elemIdx] = 1;
                    }
                }
                else {
                    if (j<4) {
                        parts[elemIdx] = 2;
                    }
                    else {
                        parts[elemIdx] = 3;
                    }
                }
            }
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts); // ownerFirst = false, addCornerCells = false, overlapLayerSize = 1
        // IT's not necessary to change the default values since the LGRs are fully interior.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,3}, {3,2,2}};
        const std::vector<std::array<int,3>> endIJK_vec = {{2,1,2}, {1,1,4}, {4,3,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
        // LGR1 element indices = 0,1,80,81 -> 32 refined cells                       LGR1 dim 4x2x4
        // LGR2 element indices = 240  -> refined into 3x3x3 = 27 cells               LGR2 dim 3x3x3
        // LGR3 element indices = 183 -> refined into 4x4x4 = 64 cells                LGR3 dim 4x4x4
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Leaf grid view total amount of cells: 10x8x8 -6 + 32 + 27 + 64 = 757
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 4x2x4 -> 5x3x5 = 75 points
        // LGR2 dim 3x3x3 -? 4x4x4 = 64 points
        // LGR3 dim 4x4x4 -> 5x5x5 = 125 points
        const std::vector<int>& expected_point_ids_per_lgr = { 75, 64, 125};
        for (int lgr = 1; lgr < 4; ++lgr) {
            std::vector<int> local_point_ids;
            local_point_ids.reserve(expected_point_ids_per_lgr[lgr-1]);
            for (const auto& element : elements(grid.levelGridView(lgr))) {
                for (int corner = 0; corner < 8; ++corner)
                {
                    const auto& point = element.subEntity<3>(corner);
                    local_point_ids.push_back(grid.currentData()[lgr]->globalIdSet().id(point));
                }
            }
            auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
            const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());
            BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == expected_point_ids_per_lgr[lgr-1]);
        }

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPoint ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());

        // Leaf grid grid total ampunt of points: 11x9x9 + (75-18) + (64-8) + (125-8) = 1121
        BOOST_CHECK( allGlobalIds_points_set.size() == 1121 );
    }
}

BOOST_AUTO_TEST_CASE(atLeastOneLgr_per_process_attempt)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    // Fully interior LGRs
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts); // ownerFirst = false, addCornerCells = false, overlapLayerSize =1
        // It's not necessary to change the default values since the LGRs are fully interior.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}, {2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,1,0}, {0,0,2}, {3,2,0}, {3,0,2}};
        const std::vector<std::array<int,3>> endIJK_vec = {{1,3,1}, {1,1,3}, {4,3,1}, {4,2,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
        // LGR1 element indices = 4,8 in rank 0. Total 16 refined cells, 45 points (45-12 = 33 with new global id).
        // LGR2 element indices = 24 in rank 1. Total 27 refined cells, 64 points (64-8 = 56 with new global id).
        // LGR3 element indices = 11 in rank 2. Total 64 refined cells, 125 points (125-8 = 117 with new global id).
        // LGR4 element indices = 27, 31 in rank 3.Total 16 refined cells, 45 points (45-12 = 33 with new global id).
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // LGR1 dim 2x4x2 (16 refined cells)
        // LGR2 dim 3x3x3 (27 refined cells)
        // LGR3 dim 4x4x4 (64 refined cells)
        // LGR4 dim 2x4x2 (16 refined cells)
        // Total global ids in leaf grid view for cells: 36-(6 marked cells) + 16 + 27 + 64 + 16 = 153
        // Total global ids in leaf grid view for points: 80 + 33 + 56 + 117 + 33 = 319
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 2x4x2 -> 3x5x3 = 45 points
        // LGR2 dim 3x3x3 -> 4x4x4 = 64 points
        // LGR3 dim 4x4x4 -> 5x5x5 = 125 points
        // LGR4 dim 2x4x2 -> 3x5x3 = 45 points
        const std::vector<int>& expected_point_ids_per_lgr = { 45, 64, 125, 45};
        for (int lgr = 1; lgr < 5; ++lgr) {
            std::vector<int> local_point_ids;
            local_point_ids.reserve(expected_point_ids_per_lgr[lgr-1]);
            for (const auto& element : elements(grid.levelGridView(lgr))) {
                for (int corner = 0; corner < 8; ++corner)
                {
                    const auto& point = element.subEntity<3>(corner);
                    local_point_ids.push_back(grid.currentData()[lgr]->globalIdSet().id(point));
                }
            }
            auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
            const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());
            BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == expected_point_ids_per_lgr[lgr-1]);
        }

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPoint ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + 33 + 56 + 117 + 33 = 319
        BOOST_CHECK( allGlobalIds_points_set.size() == 319 );
    }
}


BOOST_AUTO_TEST_CASE(not_fully_interior_lgr)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    // NOT fully interior LGRs
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}, {2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,1,0}, {0,0,2}, {3,1,0}, {3,0,2}};
        const std::vector<std::array<int,3>> endIJK_vec = {{1,3,1}, {1,1,3}, {4,2,1}, {4,2,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3", "LGR4"};
        // LGR1 element indices = 4,8 in rank 0. Total 16 refined cells, 45 points (45-12 = 33 with new global id).
        // LGR2 element indices = 24 in rank 1. Total 27 refined cells, 64 points (64-8 = 56 with new global id).
        // LGR3 element indices = 7 in rank 2. This cell is interior but it has a neighboring cell sharing its top face, cell 19 belonging to rank 3.
        // LGR4 element indices = 27, 31 in rank 3.Total 16 refined cells, 45 points (45-12 = 33 with new global id).

        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 2x4x2 -> 3x5x3 = 45 points
        // LGR2 dim 3x3x3 -> 4x4x4 = 64 points
        // LGR3 dim 4x4x4 -> 5x5x5 = 125 points
        // LGR4 dim 2x4x2 -> 3x5x3 = 45 points
        const std::vector<int>& expected_point_ids_per_lgr = { 45, 64, 125, 45};
        for (int lgr = 1; lgr < 5; ++lgr) {
            std::vector<int> local_point_ids;
            local_point_ids.reserve(expected_point_ids_per_lgr[lgr-1]);
            for (const auto& element : elements(grid.levelGridView(lgr))) {
                for (int corner = 0; corner < 8; ++corner)
                {
                    const auto& point = element.subEntity<3>(corner);
                    local_point_ids.push_back(grid.currentData()[lgr]->globalIdSet().id(point));
                }
            }
            auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
            const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());
            BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == expected_point_ids_per_lgr[lgr-1]);
        }

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPoint ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + 33 + 56 + 117 + 33 = 319
        BOOST_CHECK( allGlobalIds_points_set.size() == 319 );
    }
}

//Calling globalRefine on a distributed grid is supported now.
BOOST_AUTO_TEST_CASE(globalRefine1)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,2,1};
    grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> parts(8);
    std::vector<std::vector<int>> cells_per_rank = {{0,4},{1,5},{2,6}, {3,7}};
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    // Distribute the grid
    if(grid.comm().size()>1)
    {
        grid.loadBalance();

        grid.globalRefine(1);
        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{4,2,1}};
        const std::vector<std::string> lgr_name_vec = {"GR1"}; // GR stands for GLOBAL REFINEMENT

        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
}

BOOST_AUTO_TEST_CASE(globalRefine2)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

     std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    // Distribute the grid
    if(grid.comm().size()>1)
    {
        grid.loadBalance();

        grid.globalRefine(1);
        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{4,3,3}};
        const std::vector<std::string> lgr_name_vec = {"GR1"}; // GR stands for GLOBAL REFINEMENT
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
}

BOOST_AUTO_TEST_CASE(distributed_lgr)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,1,1}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = 1 (rank 0), 2 (rank 2). Total 16 refined cells, 45 points (45-12 = 33 with new global id).
        // LGR1 dim 4x2x2 (16 refined cells) (45 points - only 33 new points)

        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 4x2x2 -> 5x3x3 = 45 points
        std::vector<int> local_point_ids;
        local_point_ids.reserve(45); // expected_point_ids in LGR1
        for (const auto& element : elements(grid.levelGridView(1))) {
            for (int corner = 0; corner < 8; ++corner)
            {
                const auto& point = element.subEntity<3>(corner);
                local_point_ids.push_back(grid.currentData()[1]->globalIdSet().id(point));
            }
        }
        auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
        const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());
        BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == 45); // LGR1 dim 4x2x2 -> 5x3x3 = 45 points

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPointLeaf ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + (45 - 12) = 113
        BOOST_CHECK( allGlobalIds_points_set.size() == 113 );
    }
}

BOOST_AUTO_TEST_CASE(distributed_lgr_II)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);

    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{0,2,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,3,1}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = 8,9 (rank 0), 10 (rank 2). Total 24 refined cells, 63 points (63-16 = 47 with new global id).

        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 6x2x2 -> 7x3x3 = 63 points
        std::vector<int> local_point_ids;
        local_point_ids.reserve(63); // expected_point_ids in LGR1
        for (const auto& point : vertices(grid.levelGridView(1))) {
            local_point_ids.push_back(grid.currentData()[1]->globalIdSet().id(point));
        }
        auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
        const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());
        BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == 63);

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPointLeaf ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + (63 - 16) = 127
        BOOST_CHECK( allGlobalIds_points_set.size() == 127 );
    }
}

BOOST_AUTO_TEST_CASE(distributed_in_all_ranks_lgr)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.
        /** This is not enough to get unique point ids. */

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,2,2}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = {1,2,5,6,13,14,17,18} where
        // 1,5 in rank 0,
        // 13,17 in rank 1,
        // 2,6,18 in rank 2,
        // 14 in rank 3.
        // Block of cells to refine dim 2x2x2. LGR1 dim 4x4x4.
        // 64 new refined cells. 5x5x5 = 125 points (only 98 = 125 - 3x3x3 parent corners new points - new global ids).
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 4x4x4 -> 5x5x5 = 125 points
        std::vector<int> local_point_ids;
        local_point_ids.reserve(125); // expected_point_ids in LGR1

        for (const auto& element : elements(grid.levelGridView(1))) {
            for (int corner = 0; corner < 8; ++corner)
            {
                const auto& point = element.subEntity<3>(corner);
                local_point_ids.push_back(grid.currentData()[1]->globalIdSet().id(point));
            }
        }
        auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
        const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());

        // Remark: the LGR is distributed in ALL processes BUT
        // - P0, P1, P2 see the entire LGR (due to the argument addCornerCells = true).
        // However, P3 does NOT see cell 5, but sees all the other parent cells {1,2,6,13,14,17,18},
        // In P3, {1,6,11,17,21} are (overlap) added corner cells, and {2,3,7,13, 18,22,,25,29,33} are overlap
        // cells coming from overlapLayerSize = 1.
        // P3 creates +7 duplicated ids on the unseen-in-P3 cell 5
        BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == (125 + 7) ); // expected +  7 "multi-ids"

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPointLeaf ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + (125 - 27) = 178
        BOOST_CHECK( allGlobalIds_points_set.size() == (178 + 7) ); //  expected +  7 "multi-ids"
    }
}

BOOST_AUTO_TEST_CASE(distributed_in_all_ranks_lgr_II)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.
        /** This is not enough to get unique point ids. */

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,1}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,2,3}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = {13,14,17,18,25,26,29,30} where
        // 13,17,25,29 in rank 1,
        // 18 in rank 2,
        // 14,26,30 in rank 3.
        // Block of cells to refine dim 2x2x2. LGR1 dim 4x4x4.
        // 64 new refined cells. 5x5x5 = 125 points (only 98 = 125 - 3x3x3 parent corners new points - new global ids).
        grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

        // Check global id is not duplicated for points for each LGR
        // LGR1 dim 4x4x4 -> 5x5x5 = 125 points
        std::vector<int> local_point_ids;
        local_point_ids.reserve(125); // expected_point_ids in LGR1

        for (const auto& element : elements(grid.levelGridView(1))) {
            for (int corner = 0; corner < 8; ++corner)
            {
                const auto& point = element.subEntity<3>(corner);
                local_point_ids.push_back(grid.currentData()[1]->globalIdSet().id(point));
            }
        }
        auto [all_point_ids, displPoint ] = Opm::allGatherv(local_point_ids, grid.comm());
        const std::set<int> all_point_ids_set(all_point_ids.begin(), all_point_ids.end());

        // Remark: the LGR is distributed in P1-P3 processes BUT
        // - P0 does NOT see cells 25,26,30 R (even though the argument addCornerCells = true has been used).
        // - P2 does NOT see cell 25.
        // - P1 and P3 sees the entire LGR (due to the argument addCornerCells = true).
        // 4 extra unnecessary point ids are created.
        BOOST_CHECK( static_cast<int>(all_point_ids_set.size()) == (125 + 4) ); // expected +  4 "multi-ids"

        // Check global id is not duplicated for points
        std::vector<int> localPointIds_vec;
        localPointIds_vec.reserve(grid.currentData().back()->size(3));
        for (const auto& point : vertices(grid.leafGridView())) {
            // Notice that all partition type points are pushed back. Selecting only interior points does not bring us to the expected value.
            localPointIds_vec.push_back(grid.currentData().back()->globalIdSet().id(point));
        }
        auto [allGlobalIds_points, displPointLeaf ] = Opm::allGatherv(localPointIds_vec, grid.comm());
        const std::set<int> allGlobalIds_points_set(allGlobalIds_points.begin(), allGlobalIds_points.end());
        // Total global ids in leaf grid view for points: 80 + (125 - 27) = 178
        BOOST_CHECK( allGlobalIds_points_set.size() == (178 + 4) ); //  expected +  4 "multi-ids"
    }
}



BOOST_AUTO_TEST_CASE(call_adapt_on_distributed_grid)
{
    // Only for testing assignment of new global ids for refined entities (cells and point belonging to
    // refined level grids).

    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
                                                     {12,13,17,24,25,28,29,32,33},
                                                     {2,3,6,7,10,11,18,22,23},
                                                     {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts, false, true); // ownerFirst = false, addCornerCells = true, overlapLayerSize = 1
        // We set addCornersCells to achieve unique global ids for points belonging to refined level grids.


        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{3,2,2}};
        const std::vector<std::string> lgr_name_vec = {"LGR1"};
        // LGR1 element indices = {1,2,5,6,13,14,17,18} where
        // 1,5 in rank 0,
        // 13,17 in rank 1,
        // 2,6,18 in rank 2,
        // 14 in rank 3.
        // Block of cells to refine dim 2x2x2. LGR1 dim 4x4x4.
        // 64 new refined cells. 5x5x5 = 125 points (only 98 = 125 - 3x3x3 parent corners new points - new global ids).
        const std::vector<int>& marked_elemIdx = {1,2,5,6,13,14,17,18};
        std::vector<int> assignRefinedLevel(grid.currentData().front()->size(0));
        for (const auto& idx : marked_elemIdx)
            assignRefinedLevel[idx] = 1;

        grid.adapt(cells_per_dim_vec,
                   assignRefinedLevel,
                   lgr_name_vec,
                   true,
                   startIJK_vec,
                   endIJK_vec);

        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
}

// Call adapt() on a distribted grid, marking all the elements (equivalent to call globalRefine).
BOOST_AUTO_TEST_CASE(call_adapt_on_full_distributed_grid)
{
    // Only for testing assignment of new global ids for refined entities (cells and point belonging to
    // refined level grids).

    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    std::vector<int> parts(36);
    std::vector<std::vector<int>> cells_per_rank = { {0,1,4,5,8,9,16,20,21},
        {12,13,17,24,25,28,29,32,33},
        {2,3,6,7,10,11,18,22,23},
        {14,15,19,26,27,30,31,34,35} };
    for (int rank = 0; rank < 4; ++rank) {
        for (const auto& elemIdx : cells_per_rank[rank]) {
            parts[elemIdx] = rank;
        }
    }
    if(grid.comm().size()>1)
    {
        grid.loadBalance(parts);
        // Mark all elements -> 'indirect' global refinement
        for (const auto& element : elements(grid.leafGridView())){
            grid.mark(1, element);
        }
        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}};
        const std::vector<std::array<int,3>> startIJK_vec = {{1,0,0}};
        const std::vector<std::array<int,3>> endIJK_vec = {{4,3,3}};
        const std::vector<std::string> lgr_name_vec = {"GR1"};

        refinePatch_and_check(grid, cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    }
}

