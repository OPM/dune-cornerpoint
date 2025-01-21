// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#include <dune/common/version.hh>

#define BOOST_TEST_MODULE GraphRepresentationOfGridParallel
#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <opm/grid/CpGrid.hpp>

#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>

// #include <opm/grid/utility/OpmWellType.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#if HAVE_MPI
BOOST_AUTO_TEST_CASE(ExtendRootExportList)
{
    Dune::CpGrid grid;
    std::array<int, 3> dims { 10, 1, 1 };
    std::array<double, 3> size { 1., 1., 1. };
    grid.createCartesian(dims, size);
    const auto& cc = grid.comm();
    Opm::GraphOfGrid gog(grid);
    if (cc.rank() == 0) {
        gog.addWell(std::set<int> { 0, 1 }); // goes to rank 0
        gog.addWell(std::set<int> { 3, 4 }); // goes to rank 3
        gog.addWell(std::set<int> { 6, 7 }); // goes to rank 2
    }
    std::array<int, 4> ranks;
    for (int i = 0; i < 4; ++i) {
        ranks[i] = std::min(i, cc.size() - 1);
    }
    std::vector<int> gIDtoRank { 0, 0, 0, 3, 3, 3, 2, 2, 1, 1 };
    for (auto& v : gIDtoRank) {
        v = ranks[v];
    }
    // all cells from the root (including wells) are added manually before using extendRootExportList
    std::vector<int> exportedID { 5, 0, 6, 3, 8, 2, 9, 1 };
    std::vector<int> exportedRank { 3, 0, 2, 3, 1, 0, 1, 0 };
    for (auto& v : exportedRank) {
        v = ranks[v];
    }
    std::vector<std::tuple<int, int, char>> exportList, exportList2;
    exportList.reserve(grid.size(0));
    for (int i = 0; i < 8; ++i) {
        using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
        exportList.emplace_back(exportedID[i], exportedRank[i], AttributeSet::owner);
    }
    exportList2 = exportList;
    auto exportCells = Opm::Impl::extendRootExportList(gog, exportList, 0, gIDtoRank);
    // check if not using gIDtoRank leads to the same result
    auto exportCells2 = Opm::Impl::extendRootExportList(gog, exportList2, 0, std::vector<int> {});
    BOOST_REQUIRE(exportList == exportList2);
    BOOST_REQUIRE(exportCells == exportCells2);
    if (cc.rank() == 0) {
        switch (ranks[3]) {
        case 3:
            BOOST_REQUIRE((int)exportCells.size() == cc.size());
            BOOST_REQUIRE(exportCells[0].size() == 0);
            BOOST_REQUIRE(exportCells[1].size() == 0);
            BOOST_REQUIRE(exportCells[2].size() == 1);
            BOOST_CHECK(exportCells[2][0] == 7);
            BOOST_REQUIRE(exportCells[3].size() == 1);
            BOOST_CHECK(exportCells[3][0] == 4);
            break;
        case 2:
            BOOST_REQUIRE(exportCells.size() == 3);
            BOOST_REQUIRE(exportCells[0].size() == 0);
            BOOST_REQUIRE(exportCells[1].size() == 0);
            BOOST_REQUIRE(exportCells[2].size() == 2);
            std::sort(exportCells[2].begin(), exportCells[2].end());
            BOOST_CHECK(exportCells[2][0] == 4);
            BOOST_CHECK(exportCells[2][1] == 7);
            break;
        case 1:
            BOOST_REQUIRE(exportCells.size() == 2);
            BOOST_REQUIRE(exportCells[0].size() == 0);
            BOOST_REQUIRE(exportCells[1].size() == 2);
            std::sort(exportCells[1].begin(), exportCells[1].end());
            BOOST_CHECK(exportCells[1][0] == 4);
            BOOST_CHECK(exportCells[1][1] == 7);
            break;
        case 0:
            BOOST_REQUIRE(exportCells.size() == 1);
            BOOST_REQUIRE(exportCells[0].size() == 0);
            break;
        default:
            OPM_THROW(std::logic_error, "The test is designed for at most 4 ranks.");
        }
        if (ranks[3] > 0) {
            BOOST_REQUIRE(exportList.size() == 10);
            for (int i = 0; i < 10; ++i) {
                BOOST_CHECK(std::get<0>(exportList[i]) == i);
                BOOST_CHECK(std::get<1>(exportList[i]) == gIDtoRank[i]);
            }
        }
    } else {
        // non-root ranks return an empty vector, no matter what arguments they have got
        BOOST_REQUIRE(exportCells.size() == 0);
    }
}

BOOST_AUTO_TEST_CASE(CommunicateExportedCells)
{
    const auto& cc = Dune::MPIHelper::getCommunication();
    std::array<int, 4> ranks;
    for (int i = 0; i < 4; ++i) {
        ranks[i] = std::min(i, cc.size() - 1);
    }

    // root (rank 0) exports, other ranks do not
    std::vector<std::vector<int>> exportCells;
    if (cc.rank() == 0) {
        exportCells.resize(cc.size());
        for (int r = 1; r <= ranks[3]; ++r) {
            exportCells[r].push_back(r * 10);
            exportCells[r].push_back(r * 10 + 1);
            exportCells[r].push_back(r * 10 + 2);
        }
    }

    auto cells = Opm::Impl::communicateExportedCells(exportCells, cc, 0);

    if (cc.rank() == 0) {
        BOOST_REQUIRE(cells.size() == 0); // root receives nothing
    } else {
        if (cc.rank() <= ranks[3]) {
            BOOST_REQUIRE(cells.size() == 3);
            for (int i = 0; i < 3; ++i) {
                BOOST_CHECK(cells[i] == cc.rank() * 10 + i);
            }
        } else {
            BOOST_REQUIRE(cells.size() == 0); // nothing assigned to ranks 4 and up
        }
    }
}

namespace {
auto createWell(const std::string& name)
{
    using namespace Opm;
    return Dune::cpgrid::OpmWellType(name, name, 0, 0, 0, 0, 0., WellType(),
        Well::ProducerCMode(), Connection::Order(), UnitSystem(),
        0., 0., false, false, 0, Well::GasInflowEquation());
};
}

BOOST_AUTO_TEST_CASE(WellsOnThisRank)
{
    const auto& cc = Dune::MPIHelper::getCommunication();
    std::vector<Dune::cpgrid::OpmWellType> wells;
    wells.reserve(6);
    wells.push_back(createWell("OnRank3a"));
    wells.push_back(createWell("OnRank1a"));
    wells.push_back(createWell("OnRank0"));
    wells.push_back(createWell("OnRank2"));
    wells.push_back(createWell("OnRank1b"));
    wells.push_back(createWell("OnRank3b"));
    std::array<int, 4> ranks;
    for (int i = 0; i < 4; ++i) {
        ranks[i] = std::min(i, cc.size() - 1);
    }
    std::vector<int> wellRanks { 3, 1, 0, 2, 1, 3 };
    for (auto& v : wellRanks) {
        v = ranks[v];
    }

    const auto wotr = wellsOnThisRank(wells, wellRanks, cc, 0);

    BOOST_REQUIRE(wotr.size() == 6); // number of wells
    // output is sorted by name in computeParallelWells
    BOOST_CHECK(wotr[0].first == "OnRank0");
    BOOST_CHECK(wotr[1].first == "OnRank1a");
    BOOST_CHECK(wotr[2].first == "OnRank1b");
    BOOST_CHECK(wotr[3].first == "OnRank2");
    BOOST_CHECK(wotr[4].first == "OnRank3a");
    BOOST_CHECK(wotr[5].first == "OnRank3b");
    BOOST_CHECK(wotr[0].second == (cc.rank() == ranks[0]));
    BOOST_CHECK(wotr[1].second == (cc.rank() == ranks[1]));
    BOOST_CHECK(wotr[2].second == (cc.rank() == ranks[1]));
    BOOST_CHECK(wotr[3].second == (cc.rank() == ranks[2]));
    BOOST_CHECK(wotr[4].second == (cc.rank() == ranks[3]));
    BOOST_CHECK(wotr[5].second == (cc.rank() == ranks[3]));
}

// After partitioning, importList and exportList are not complete,
// other cells from wells need to be added.
BOOST_AUTO_TEST_CASE(ImportExportListExpansion)
{
    Dune::CpGrid grid;
    std::array<int, 3> dims { 3, 3, 2 };
    std::array<double, 3> size { 1., 1., 1. };
    grid.createCartesian(dims, size);
    const auto& cc = grid.comm();
    if (cc.size() == 1)
        return;

    Opm::GraphOfGrid gog(grid);
    // grid is nonempty only on the rank 0
    if (cc.rank() == 0) {
        gog.addWell(std::set<int> { 0, 1, 2 });
        gog.addWell(std::set<int> { 3, 4, 5 });
        gog.addWell(std::set<int> { 6, 7, 8 });
        gog.addWell(std::set<int> { 9, 13, 17 });
        BOOST_REQUIRE(gog.size() == 10);
    }

    // rank-specific export and import lists
    std::vector<std::tuple<int, int, char>> exportList, exportSolution;
    std::vector<std::tuple<int, int, char, int>> importList, importSolution;
    // this test works on any number or ranks although from rank 4 (including) all are empty
    // if ranks<4, the highest rank gobbles leftovers
    int maxrank = cc.size() - 1;
    std::vector<int> ranks { 0, std::min(maxrank, 1), std::min(maxrank, 2), std::min(maxrank, 3) };

    // cells[rank] holds solution, each rank 0..3 has 1 well
    std::vector<std::vector<int>> cells(4);
    cells[0] = std::vector<int> { 0, 1, 2, 10, 11 };
    cells[1] = std::vector<int> { 3, 4, 5, 12 };
    cells[2] = std::vector<int> { 6, 7, 8, 15, 16 };
    cells[3] = std::vector<int> { 9, 13, 14, 17 };
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
    char owner = AttributeSet::owner;

    for (int i = 0; i < 4; ++i) {
        for (const auto& c : cells[i]) {
            if (cc.rank() == ranks[i]) {
                importSolution.push_back(std::make_tuple(c, ranks[i], owner, -1));
            }
            if (cc.rank() == 0) {
                exportSolution.push_back(std::make_tuple(c, ranks[i], owner));
            }
        }
    }
    std::sort(importSolution.begin(), importSolution.end());
    std::sort(exportSolution.begin(), exportSolution.end());

    if (cc.rank() == ranks[0]) {
        // Zoltan does not include cells that remain on the rank into import and export list
        // but they are added manually to both lists (cell on root is in its import AND export)
        // before extendAndSortExportAndImportLists is called
        for (const auto& gID : cells[ranks[0]]) {
            importList.push_back(std::make_tuple(gID, ranks[0], owner, -1));
            exportList.push_back(std::make_tuple(gID, ranks[0], owner));
        }
        exportList.push_back(std::make_tuple(3, ranks[1], owner));
        exportList.push_back(std::make_tuple(12, ranks[1], owner));
        exportList.push_back(std::make_tuple(6, ranks[2], owner));
        exportList.push_back(std::make_tuple(15, ranks[2], owner));
        exportList.push_back(std::make_tuple(16, ranks[2], owner));
        exportList.push_back(std::make_tuple(9, ranks[3], owner));
        exportList.push_back(std::make_tuple(14, ranks[3], owner));
        std::sort(exportList.begin(), exportList.end());
    } else if (cc.rank() == ranks[1]) {
        // non-root ranks have empty export list
        // import list is not sorted
        importList.push_back(std::make_tuple(12, ranks[1], owner, -1));
        importList.push_back(std::make_tuple(3, ranks[1], owner, -1));
    }
    // no else branch, ranks[2]==1 when cc.size==2
    if (cc.rank() == ranks[2]) {
        importList.push_back(std::make_tuple(15, ranks[2], owner, -1));
        importList.push_back(std::make_tuple(6, ranks[2], owner, -1));
        importList.push_back(std::make_tuple(16, ranks[2], owner, -1));
    }
    if (cc.rank() == ranks[3]) {
        importList.push_back(std::make_tuple(9, ranks[3], owner, -1));
        importList.push_back(std::make_tuple(14, ranks[3], owner, -1));
    }

    extendAndSortExportAndImportLists(gog, cc, ranks[0], exportList, importList);

    BOOST_CHECK_MESSAGE(importList == importSolution, "On rank " + std::to_string(cc.rank()));
    BOOST_CHECK_MESSAGE(exportList == exportSolution, "On rank " + std::to_string(cc.rank()));
}

// Sequential Zoltan partitioner has a simpler structure and relies on
// extendGIDtoRank and makeExportListsFromGIDtoRank to make imp-/export lists.
BOOST_AUTO_TEST_CASE(SequentialZoltanSupport)
{
    Dune::CpGrid grid;
    std::array<int, 3> dims { 3, 3, 2 };
    std::array<double, 3> size { 1., 1., 1. };
    grid.createCartesian(dims, size);
    const auto& cc = grid.comm();
    if (cc.size() == 1)
        return;

    constexpr int root = 0;
    // this test works on any number or ranks although from rank 4 (including) all are empty
    // if ranks<4, the highest rank gobbles leftovers
    int maxrank = cc.size() - 1;
    std::vector<int> ranks { 0, std::min(maxrank, 1), std::min(maxrank, 2), std::min(maxrank, 3) };

    std::vector<int> gIDtoRank, importedCells, importSol;
    std::vector<std::vector<int>> exportedCells;
    std::vector<std::vector<int>> rankSol {{ 3, 4, 5, 12, 14 },
                                           { 0, 1, 2, 10, 11 },
                                           { 6, 7, 8, 15 },
                                           { 9, 13, 16, 17 }};
    for (int i = 0; i < 4; ++i) {
        if (cc.rank() == ranks[i])
            importSol.insert(importSol.end(), rankSol[i].begin(), rankSol[i].end());
    }
    std::sort(importSol.begin(),importSol.end());

    // grid is nonempty only on the rank 0
    if (cc.rank() == root) {
        Opm::GraphOfGrid gog(grid);
        gog.addWell(std::set<int> { 0, 1, 2 });
        gog.addWell(std::set<int> { 3, 4, 5 });
        gog.addWell(std::set<int> { 6, 7, 8 });
        gog.addWell(std::set<int> { 9, 13, 17 });
        BOOST_REQUIRE(gog.size() == 10);

        // setup vector of IDs to ranks
        gIDtoRank.resize(grid.numCells(), root);
        assert(gIDtoRank.size() == 18);
        // what partitioner sees
        gIDtoRank[0] = ranks[1];
        gIDtoRank[3] = ranks[0];
        gIDtoRank[6] = ranks[2];
        gIDtoRank[9] = ranks[3];
        gIDtoRank[10] = ranks[1];
        gIDtoRank[11] = ranks[1];
        gIDtoRank[12] = ranks[0];
        gIDtoRank[14] = ranks[0];
        gIDtoRank[15] = ranks[2];
        gIDtoRank[16] = ranks[3];
        // add well cells
        extendGIDtoRank(gog, gIDtoRank, root);
        // rearrange into the structure fitting for communication - vector per rank
        exportedCells = Opm::makeExportListsFromGIDtoRank(gIDtoRank, cc.size());
        // what stays on root
        importedCells.swap(exportedCells[root]);
    }

    auto impCells = Opm::Impl::communicateExportedCells(exportedCells, cc, root);
    if (cc.rank() != root)
        impCells.swap(importedCells);
    BOOST_REQUIRE(importedCells == importSol);
}
#endif // HAVE_MPI

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
