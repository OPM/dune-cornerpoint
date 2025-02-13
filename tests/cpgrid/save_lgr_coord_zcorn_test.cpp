//===========================================================================
//
// File: save_lgr_coord_zcorn_test.cpp
//
// Created: Thursday 13.02.2025 15:17:00
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

#define BOOST_TEST_MODULE SaveLgrCoordZcorn
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
#include <opm/io/eclipse/EGrid.hpp> // check if it's needed

#include <array>
#include <filesystem>
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

    Dune::CpGrid grid;
    // Create the starting grid (before adding LGRs)
    Opm::Parser parser;
    const auto deck = parser.parseString(deck_string);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);


    // Loop over all lgrs
    for (const auto& [lgr_name, lgr_level] : grid.getLgrNameToLevel())
    {
        const auto [lgrCartesianIdxToCellIdx, lgrIJK] = Opm::lgrIJK(grid, lgr_name);

        /** To be improved: lgrCOORD will be already a vector<double> - soon! */
        // Assuming lgrCOORD is a vector of std::array<double, 6> with the content of COORD keyword for the LGR block
        const auto [lgrCOORD, lgrZCORN] = Opm::lgrCOORDandZCORN(grid, lgr_level, lgrCartesianIdxToCellIdx, lgrIJK);


        /** Flatten lgrCOORD will not be needed - soon! */
        std::vector<double> flatLgrCOORD(6*lgrCOORD.size());
        for (const auto& pillar_top_bottom : lgrCOORD) {
            flatLgrCOORD.insert(flatLgrCOORD.end(), pillar_top_bottom.begin(), pillar_top_bottom.end());
        }

        
        eclipse_grid.set_lgr_refinement(lgr_name, flatLgrCOORD, lgrZCORN);
     
        Opm::UnitSystem units(1);
        std::vector<Opm::NNCdata> vecNNC;
        eclipse_grid.save("TEST"+lgr_name+".EGRID",false, vecNNC, units);
    }
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


    // Loop over all lgrs
    for (const auto& [lgr_name, lgr_level] : grid.getLgrNameToLevel())
    {
        const auto [lgrCartesianIdxToCellIdx, lgrIJK] = Opm::lgrIJK(grid, lgr_name);

        /** To be improved: lgrCOORD will be already a vector<double> - soon! */
        // Assuming lgrCOORD is a vector of std::array<double, 6> with the content of COORD keyword for the LGR block
        const auto [lgrCOORD, lgrZCORN] = Opm::lgrCOORDandZCORN(grid, lgr_level, lgrCartesianIdxToCellIdx, lgrIJK);


        /** Flatten lgrCOORD will not be needed - soon! */
        std::vector<double> flatLgrCOORD(6*lgrCOORD.size());
        for (const auto& pillar_top_bottom : lgrCOORD) {
            flatLgrCOORD.insert(flatLgrCOORD.end(), pillar_top_bottom.begin(), pillar_top_bottom.end());
        }


        eclipse_grid.set_lgr_refinement(lgr_name, flatLgrCOORD, lgrZCORN);

        Opm::UnitSystem units(1);
        std::vector<Opm::NNCdata> vecNNC;
        eclipse_grid.save("TEST"+lgr_name+".EGRID",false, vecNNC, units);
    }

}

