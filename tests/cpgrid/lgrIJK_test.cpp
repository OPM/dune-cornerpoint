//===========================================================================
//
// File: lgrIJK_test.cpp
//
// Created: Thursday 30.01.2025 08:17:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2025 Equinor ASA.

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

#define BOOST_TEST_MODULE LgrIJKTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif


#include <opm/grid/cpgrid/CpGridUtilities.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>


#include <array>
#include <stdexcept>
#include <string>
#include <vector>

struct Fixture {
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

// Create a grid and add one test LGR with dimension 6x6x3.
Dune::CpGrid
createGridAndAddTestLgr(const std::string& deck_string)
{
    Dune::CpGrid grid;
    // Create the starting grid (before adding LGRs)
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);


    // Add LGR1 and update grid view
    const std::vector<std::array<int, 3>> cells_per_dim_vec
        = {{3, 3, 3}}; // 3x3x3 child cells in x-,y-, and z-direction per ACTIVE parent cell
    const std::vector<std::array<int, 3>> startIJK_vec
        = {{1, 1, 0}}; // starts at (1,1,0) in coarse grid - equivalent to (I1-1, J1-1, K1-1) from its CARFIN block
    const std::vector<std::array<int, 3>> endIJK_vec
        = {{3, 3, 1}}; // ends at (3,3,1) in coarse grid - equivalent to (I2, J2, K2) from its CARFIN block
    const std::vector<std::string> lgr_name_vec = {"LGR1"};

    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);
    return grid;
}

void
checkExpectedSize(const int& expected_elements,
                  const std::size_t& numLgrCells,
                  const std::size_t& cellIdxToLgrCartesianIdxSize,
                  const std::size_t& lgrCartesianIdxToCellIdxSize,
                  const std::size_t& lgr1IJKSize)
{
    BOOST_CHECK_EQUAL(numLgrCells, expected_elements);
    BOOST_CHECK_EQUAL(cellIdxToLgrCartesianIdxSize, expected_elements);
    BOOST_CHECK_EQUAL(lgrCartesianIdxToCellIdxSize, expected_elements);
    BOOST_CHECK_EQUAL(lgr1IJKSize, expected_elements);
}


BOOST_AUTO_TEST_CASE(fullActiveParentCellsBlock)
{
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  3 3 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /
TOPS
	9*8325 /
 ACTNUM
        1 1 1
        1 1 1
        1 1 1
        /
PORO
  9*0.15 /
PERMX
  9*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";

    const auto grid = createGridAndAddTestLgr(deck_string);

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    // Get LGR level
    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");
    const int numLgrCells = grid.levelGridView(lgr1_level).size(0);

    const auto& cellIdxToLgrCartesianIdx = grid.currentData()[lgr1_level]->globalCell();

    // Verify the size matches expected elements
    const int expected_elements = 108; // 4 parent cells into 3x3x3 children each -> 108
    checkExpectedSize(expected_elements,
                      numLgrCells,
                      cellIdxToLgrCartesianIdx.size(),
                      lgrCartesianIdxToCellIdx.size(),
                      lgr1IJK.size());

    // Validate all ijk's are non-negative
    for (const auto& ijk : lgr1IJK) {
        BOOST_TEST(ijk[0] >= 0);
        BOOST_TEST(ijk[1] >= 0);
        BOOST_TEST(ijk[2] >= 0);
    }

    // Invalid LGR should throw an exception
    BOOST_CHECK_THROW(Opm::lgrIJK(grid, "LGR2DOESNOTEXIST"), std::runtime_error);

    // LGR1 dimension 6x6x3
    //  Visual representation per k-layer:
    //
    //       level Cartesian cell indices                      internal ordering, i.e. element.index() in
    //                                                         levelGridView(level);
    //
    // k = 2  | 102 103 104  | 105 106 107                     |  78  79  80  | 105 106 107
    //        |  96  97  98  |  99 100 101                     |  75  76  77  | 102 103 104
    //        |  90  91  92  |  93  94  95                     |  72  73  74  |  99 100 101
    //        ----------------------------                     ----------------------------
    //        |  84  85  86  |  87  88  89                     |  24  25  26  |  51  52  53
    //        |  78  79  80  |  81  82  83                     |  21  22  23  |  48  49  50
    //        |  72  73  74  |  75  76  77                     |  18  19  20  |  45  46  47
    //------------------------------------                     ----------------------------
    // k = 1  |  66  67  68  |  69  70  71                     |  69  70  71  |  96  97  98
    //        |  60  61  62  |  63  64  65                     |  66  67  68  |  93  94  95
    //        |  54  55  56  |  57  58  59                     |  63  64  65  |  90  91  92
    //        ----------------------------                     ----------------------------
    //        |  48  49  50  |  51  52  53                     |  15  16  17  |  42  43  44
    //        |  42  43  44  |  45  46  47                     |  12  13  14  |  39  40  41
    //        |  36  37  38  |  39  40  41                     |   9  10  11  |  36  37  38
    //------------------------------------                     ----------------------------
    // k = 0  |  30  31  32  |  33  34  35                     |  60  61  62  |  87  88  89
    //        |  24  25  26  |  27  28  29                     |  57  58  59  |  84  85  86
    //        |  18  19  20  |  21  22  23                     |  54  55  56  |  81  82  83
    //        ----------------------------                     ----------------------------
    //        |  12  13  14  |  15  16  17                     |   6   7   8  |  33  34  35
    //        |   6   7   8  |   9  10  11                     |   3   4   5  |  30  31  32
    //        |   0   1   2  |   3   4   5                     |   0   1   2  |  27  28  92


    // Verify that cell with level (local) Cartesian index 25 in LGR1 corresponds to cell index 58
    // and has local ijk = {1,4,0}
    BOOST_TEST_REQUIRE(cellIdxToLgrCartesianIdx[58] == 25);
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.count(25));
    BOOST_TEST(lgrCartesianIdxToCellIdx.at(25) == 58);

    BOOST_TEST(lgr1IJK[58][0] == 1);
    BOOST_TEST(lgr1IJK[58][1] == 4);
    BOOST_TEST(lgr1IJK[58][2] == 0);


    // Verify that cell with level (local) Cartesian index 64 in LGR1 corresponds to cell index 94
    // and has local ijk = {4,4,1}
    BOOST_TEST_REQUIRE(cellIdxToLgrCartesianIdx[94] == 64);
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.count(64));
    BOOST_TEST(lgrCartesianIdxToCellIdx.at(64) == 94);

    BOOST_TEST(lgr1IJK[94][0] == 4);
    BOOST_TEST(lgr1IJK[94][1] == 4);
    BOOST_TEST(lgr1IJK[94][2] == 1);

    // Verify that cell with level (local) Cartesian index 98 in LGR1 corresponds to cell index 77
    // and has local ijk = {2,4,2}
    BOOST_TEST_REQUIRE(cellIdxToLgrCartesianIdx[77] == 98);
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.count(98));
    BOOST_TEST(lgrCartesianIdxToCellIdx.at(98) == 77);

    BOOST_TEST(lgr1IJK[77][0] == 2);
    BOOST_TEST(lgr1IJK[77][1] == 4);
    BOOST_TEST(lgr1IJK[77][2] == 2);
}


/* Tests lgrIJK() for a grid containing inactive cells within the LGR block.

   The two tests below pass if initLgr() is not invoked in opm-common/opm/input/eclipse/EclipseState/EclipseState.cpp.
   Currently, inactive cells in LGRs are not supported in "opm-common."
*/

BOOST_AUTO_TEST_CASE(lgrWithActiveAndInactiveCells, *boost::unit_test::disabled())
{
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  3 3 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /
TOPS
	9*8325 /
 ACTNUM
        1 1 1
        1 1 1
        1 1 0
        /
PORO
  9*0.15 /
PERMX
  9*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";

    const auto grid = createGridAndAddTestLgr(deck_string);

    // Note: k = 0, indicating a single-layer grid with dimensions 3x3x1.
    // ACTNUM represents the active cell indicator for the parent grid of the LGR (Local Grid Refinement).
    // The grid structure is as follows, where '1' denotes an active cell and '0' denotes an inactive cell:
    //
    //   1  1  1
    //   1  1  1
    //   1  1  0
    //
    // The corresponding LGR block parent cell representation appears as:
    //   1  1
    //   1  0


    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");
    // Total active refined cells on the level grid
    const int numLgrCells = grid.levelGridView(lgr1_level).size(0);

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    const auto& cellIdxToLgrCartesianIdx = grid.currentData()[lgr1_level]->globalCell();

    // Verify the size matches expected elements
    const int expected_elements = 81; // 3 ACTIVE parent cellS into 3x3x3 children each -> 81
    checkExpectedSize(expected_elements,
                      numLgrCells,
                      cellIdxToLgrCartesianIdx.size(),
                      lgrCartesianIdxToCellIdx.size(),
                      lgr1IJK.size());

    // LGR1 dimension 6x6x3
    //  Visual representation per k-layer:
    //
    //       level Cartesian cell indices                      internal ordering, i.e. element.index() in
    //                                                         levelGridView(level);
    //
    // k = 2  | 102 103 104  |                                 |  78  79  80  |
    //        |  96  97  98  |  INACTIVE                       |  75  76  77  | INACTIVE
    //        |  90  91  92  |                                 |  72  73  74  |
    //        ----------------------------                     ----------------------------
    //        |  84  85  86  |  87  88  89                     |  24  25  26  |  51  52  53
    //        |  78  79  80  |  81  82  83                     |  21  22  23  |  48  49  50
    //        |  72  73  74  |  75  76  77                     |  18  19  20  |  45  46  47
    //------------------------------------                     ----------------------------
    // k = 1  |  66  67  68  |                                 |  69  70  71  |
    //        |  60  61  62  |  INACTIVE                       |  66  67  68  | INACTIVE
    //        |  54  55  56  |                                 |  63  64  65  |
    //        ----------------------------                     ----------------------------
    //        |  48  49  50  |  51  52  53                     |  15  16  17  |  42  43  44
    //        |  42  43  44  |  45  46  47                     |  12  13  14  |  39  40  41
    //        |  36  37  38  |  39  40  41                     |   9  10  11  |  36  37  38
    //------------------------------------                     ----------------------------
    // k = 0  |  30  31  32  |                                 |  60  61  62  |
    //        |  24  25  26  |  INACTIVE                       |  57  58  59  | INACTIVE
    //        |  18  19  20  |                                 |  54  55  56  |
    //        ----------------------------                     ----------------------------
    //        |  12  13  14  |  15  16  17                     |   6   7   8  |  33  34  35
    //        |   6   7   8  |   9  10  11                     |   3   4   5  |  30  31  32
    //        |   0   1   2  |   3   4   5                     |   0   1   2  |  27  28  29

    // Verify that cell with level (local) Cartesian index 97 in LGR1 corresponds to cell index 76
    // and has local ijk = {1,4,2}
    BOOST_TEST_REQUIRE(cellIdxToLgrCartesianIdx[76] == 97);
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.count(97));
    BOOST_TEST(lgrCartesianIdxToCellIdx.at(97) == 76);

    BOOST_TEST(lgr1IJK[76][0] == 1);
    BOOST_TEST(lgr1IJK[76][1] == 4);
    BOOST_TEST(lgr1IJK[76][2] == 2);


    // Verify that cell with level (local) Cartesian index 88 in LGR1 corresponds to cell index 52
    // has local ijk = {4,2,2}
    BOOST_TEST_REQUIRE(cellIdxToLgrCartesianIdx[52] == 88);
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.count(97));
    BOOST_TEST_REQUIRE(lgrCartesianIdxToCellIdx.at(88) == 52);

    BOOST_TEST(lgr1IJK[52][0] == 4);
    BOOST_TEST(lgr1IJK[52][1] == 2);
    BOOST_TEST(lgr1IJK[52][2] == 2);

    // Accessing inactive (non-existing) cells
    BOOST_CHECK_THROW(lgrCartesianIdxToCellIdx.at(70), std::out_of_range);
    BOOST_CHECK_THROW(lgrCartesianIdxToCellIdx.at(170), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(fullInactiveParentCellsBlock_lgrCOORD_throws, *boost::unit_test::disabled())
{
    const std::string deck_string = R"(
RUNSPEC
DIMENS
  3 3 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  2  3  2  3  1  1  6  6  3/
ENDFIN
DX
  9*1000 /
DY
	9*1000 /
DZ
	9*20 /
TOPS
	9*8325 /
 ACTNUM
        1 1 1
        1 0 0
        1 0 0
        /
PORO
  9*0.15 /
PERMX
  9*1 /
COPY
  PERMX PERMZ /
  PERMX PERMY /
/
EDIT
OIL
GAS
TITLE
The title
START
16 JUN 1988 /
PROPS
REGIONS
SOLUTION
SCHEDULE
)";
    const auto grid = createGridAndAddTestLgr(deck_string);

    // Note: k = 0, indicating a single-layer grid with dimensions 3x3x1.
    // ACTNUM represents the active cell indicator for the parent grid of the LGR (Local Grid Refinement).
    // The grid structure is as follows, where '1' denotes an active cell and '0' denotes an inactive cell:
    //
    //   1 1 1
    //   1 0 0
    //   1 0 0
    //
    // The corresponding LGR block parent cell representation appears as:
    //   0  0
    //   0  0

    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");
    // Total active refined cells on the level grid
    const int numLgrCells = grid.levelGridView(lgr1_level).size(0); // here, 0

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    const auto& cellIdxToLgrCartesianIdx = grid.currentData()[lgr1_level]->globalCell();

    // Verify the size matches expected elements
    const int expected_elements = 0;
    checkExpectedSize(expected_elements,
                      numLgrCells,
                      cellIdxToLgrCartesianIdx.size(),
                      lgrCartesianIdxToCellIdx.size(),
                      lgr1IJK.size());

    // Accessing inactive (non-existing) cell
    BOOST_CHECK_THROW(lgrCartesianIdxToCellIdx.at(0), std::out_of_range);
}
