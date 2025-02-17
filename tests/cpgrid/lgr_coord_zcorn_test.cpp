//===========================================================================
//
// File: lgr_coord_zcorn_test.cpp
//
// Created: Thursday 13.02.2025 17:17:00
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

#define BOOST_TEST_MODULE LgrCOORDandZCORNTests
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

    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");

    const auto [lgrCOORD, lgrZCORN] = Opm::lgrCOORDandZCORN(grid, lgr1_level, lgrCartesianIdxToCellIdx, lgr1IJK);

    std::cout<< "COORD for LGR with all active cells" << std::endl;
    for (const auto& coord : lgrCOORD)
    {
        std::cout << coord << std::endl;
    }
    std::cout<< std::endl;
}


BOOST_AUTO_TEST_CASE(singleCellGrid_easyToTestlgrCOORDandZCORN)
{
    // Single-cell-grid (with dimension 1x1x1)
    // {DX,DY,DZ} = {8,4,2} which makes easy to check the values
    // of LGR COORD if the LGR dimensions are 8x4x2.

    const std::string deck_string = R"(
RUNSPEC
DIMENS
  1 1 1 /
GRID
CARFIN
-- NAME I1-I2 J1-J2 K1-K2 NX NY NZ
'LGR1'  1  1  1  1  1  1  8  4  2/
ENDFIN
DX
  1*8 /
DY
	1*4 /
DZ
	1*2 /
TOPS
	1*0 /
 ACTNUM
        1
        /
PORO
  1*0.15 /
PERMX
  1*1 /
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

    Dune::CpGrid grid;
    // Create the starting grid (before adding LGRs)
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);

    // Add LGR1 and update grid view
    const std::vector<std::array<int, 3>> cells_per_dim_vec
        = {{8, 4, 2}}; // 8x4x2 child cells in x-,y-, and z-direction per ACTIVE parent cell
    const std::vector<std::array<int, 3>> startIJK_vec
        = {{0, 0, 0}}; // starts at (0,0,0) in coarse grid - equivalent to (I1-1, J1-1, K1-1) from its CARFIN block
    const std::vector<std::array<int, 3>> endIJK_vec
        = {{1, 1, 1}}; // ends at (1,1,1) in coarse grid - equivalent to (I2, J2, K2) from its CARFIN block
    const std::vector<std::string> lgr_name_vec = {"LGR1"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);


    const int lgr1_level = grid.getLgrNameToLevel().at("LGR1");

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    const auto [lgrCOORD, lgrZCORN] = Opm::lgrCOORDandZCORN(grid, lgr1_level, lgrCartesianIdxToCellIdx, lgr1IJK);
    const int nx = 8;
    const int ny = 4;
    const int nz = 2;
    // Pillars are ordered j*(nx+1) + i, i faster than j, i=0,...,nx, j=0, ..., ny.
    for (int j = 0;  j < ny+1; ++j) {
        for (int i = 0; i < nx+1; ++i) {
            int pillar_idx = (j*6*(nx+1)) + (6*i);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +0], i);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +1], j);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +2], nz);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +3], i);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +4], j);
            BOOST_CHECK_EQUAL( lgrCOORD[pillar_idx +5], 0 );
        }
    }

    std::cout<< "COORD for LGR with all active cells" << std::endl;
    for (const auto& coord : lgrCOORD)
    {

        std::cout << coord << std::endl;
    }
    std::cout<< std::endl;

    std::cout<< "ZCORN for LGR with all active cells" << std::endl;
    for (const auto& zcorn : lgrZCORN)
    {
        std::cout << zcorn << std::endl;
    }
    std::cout<< std::endl;

    // For a grid with nz layers, ZCORN values are ordered:
    //
    //      top layer nz-1
    //   bottom layer nz-1
    //      top layer nz-2
    //   bottom layer nz-2
    // ...
    //      top layer 1
    //   bottom layer 1
    //      top layer 0
    //   bottom layer 0
    for (int k = 0; k < nz; ++k) {
        for (int j = 0;  j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int zcorn_idx =  ((nz-1-k)*8*nx*ny) + (j*4*nx) + (2*i);
                BOOST_CHECK_EQUAL( lgrZCORN[zcorn_idx], k+1); // top layer zcorn values
                BOOST_CHECK_EQUAL( lgrZCORN[zcorn_idx + (4*nx*ny)], k); // bottom layer zcorn values
            }
        }
    }
}


/* Tests for a grid containing inactive cells within the LGR block.

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

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    // If a pillar within the LGR block is "inactive," its COORD values are set to
    // std::numeric_limits<double>::max() to indicate the inactive status
    const auto [lgrCOORD, lgrZCORN] = Opm::lgrCOORDandZCORN(grid, lgr1_level, lgrCartesianIdxToCellIdx, lgr1IJK);
    std::cout<< "COORD for LGR with active and inactive cells" << std::endl;
    for (const auto& coord : lgrCOORD)
    {
        std::cout << coord << std::endl;
    }
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(fullInactiveParentCellsBlock_lgrCOORDandZCORN_throws, *boost::unit_test::disabled())
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

    const auto [lgrCartesianIdxToCellIdx, lgr1IJK] = Opm::lgrIJK(grid, "LGR1");

    // All inactive cells, therefore lgrCOORD(...) throws an exception
    BOOST_CHECK_THROW(Opm::lgrCOORDandZCORN(grid, lgr1_level, lgrCartesianIdxToCellIdx, lgr1IJK), std::logic_error);
}
