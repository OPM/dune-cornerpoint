/*
  Copyright 2014 Statoil ASA

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

#include <config.h>

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE TEST_UG
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

/* --- our own headers --- */
#include <algorithm>
#include <vector>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cornerpoint_grid.h>  /* compute_geometry */
#include <opm/grid/GridManager.hpp>  /* compute_geometry */
#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>

using namespace std;



BOOST_AUTO_TEST_CASE(Equal) {
    const std::string filename1 = "CORNERPOINT_ACTNUM.DATA";
    const char *deck2Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 10 /\n"
        "GRID\n"
        "DXV\n"
        "10*0.25 /\n"
        "DYV\n"
        "10*0.25 /\n"
        "DZV\n"
        "10*0.25 /\n"
        "TOPS\n"
        "100*0.25 /\n"
        "PORO\n"
        "   1000*0.15 /\n"
        "EDIT\n"
        "\n";

    Opm::Parser parser;

    Opm::Deck deck1 = parser.parseFile( filename1);
    Opm::EclipseState es1(deck1);

    Opm::Deck deck2 = parser.parseString( deck2Data);
    Opm::EclipseState es2(deck2);

    BOOST_CHECK( deck1.hasKeyword("ZCORN") );
    BOOST_CHECK( deck1.hasKeyword("COORD") );

    const bool edge_conformal = false;
    Opm::GridManager grid1(es1.getInputGrid(), edge_conformal);
    Opm::GridManager grid2(es2.getInputGrid(), edge_conformal);

    const UnstructuredGrid* cgrid1 = grid1.c_grid();
    const UnstructuredGrid* cgrid2 = grid2.c_grid();



    BOOST_CHECK( grid_equal( cgrid1 , cgrid1 ));
    BOOST_CHECK( grid_equal( cgrid2 , cgrid2 ));
    BOOST_CHECK( !grid_equal( cgrid1 , cgrid2 ));
}



// TODO This method might be obsolete after using EclipseState to generate grid
BOOST_AUTO_TEST_CASE(EqualEclipseGrid) {
    const std::string filename = "CORNERPOINT_ACTNUM.DATA";
    Opm::Parser parser;
    Opm::Deck deck = parser.parseFile( filename);
    Opm::EclipseState es(deck);

    Opm::GridManager gridM(es.getInputGrid());
    const UnstructuredGrid* cgrid1 = gridM.c_grid();
    struct UnstructuredGrid * cgrid2;
    {
        struct grdecl g;
        const auto& dimens = deck["DIMENS"].back();
        const auto& coord  = deck["COORD"].back();
        const auto& zcorn  = deck["ZCORN"].back();
        const auto& actnum = deck["ACTNUM"].back();

        g.dims[0] = dimens.getRecord(0).getItem("NX").get< int >(0);
        g.dims[1] = dimens.getRecord(0).getItem("NY").get< int >(0);
        g.dims[2] = dimens.getRecord(0).getItem("NZ").get< int >(0);

        g.coord  = coord.getSIDoubleData().data();
        g.zcorn  = zcorn.getSIDoubleData().data();
        g.actnum = actnum.getIntData().data();


        cgrid2 = create_grid_cornerpoint(&g, 0.0, /* edge_conformal = */ false);
        if (!cgrid2)
            throw std::runtime_error("Failed to construct grid.");
    }


    BOOST_CHECK( grid_equal( cgrid1 , cgrid2 ));

    auto actnum = Opm::UgGridHelpers::createACTNUM(*cgrid1);
    BOOST_CHECK_EQUAL( actnum.size(), 500 );
    for (std::size_t i=0; i < 100; i++) {
        BOOST_CHECK_EQUAL(actnum[i + 200], 1);
        for (std::size_t j=0; j < 2; j++) {
            BOOST_CHECK_EQUAL(actnum[i + j * 100], 0);
            BOOST_CHECK_EQUAL(actnum[i + j * 100 + 300], 0);
        }
    }
    destroy_grid( cgrid2 );
}


BOOST_AUTO_TEST_CASE(TOPS_Fully_Specified) {
    const char *deck1Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 3 /\n"
        "GRID\n"
        "DX\n"
        "300*1000 /\n"
        "DY\n"
        "300*1000 /\n"
        "DZ\n"
        "100*20 100*30  100*50 /\n"
        "TOPS\n"
        "100*8325 /\n"
        "PORO\n"
        "   300*0.15 /\n"
        "EDIT\n"
        "\n";


    const char *deck2Data =
        "RUNSPEC\n"
        "\n"
        "DIMENS\n"
        " 10 10 3 /\n"
        "GRID\n"
        "DX\n"
        "300*1000 /\n"
        "DY\n"
        "300*1000 /\n"
        "DZ\n"
        "100*20 100*30  100*50 /\n"
        "TOPS\n"
        "100*8325 100*8345  100*8375/\n"
        "PORO\n"
        "   300*0.15 /\n"
        "EDIT\n"
        "\n";

    Opm::Parser parser;
    const Opm::Deck& deck1 = parser.parseString(deck1Data);
    const Opm::Deck& deck2 = parser.parseString(deck2Data);

    Opm::EclipseState es1(deck1);
    Opm::EclipseState es2(deck2);

    Opm::GridManager gridM1(es1.getInputGrid());
    Opm::GridManager gridM2(es2.getInputGrid());

    const UnstructuredGrid* cgrid1 = gridM1.c_grid();
    const UnstructuredGrid* cgrid2 = gridM2.c_grid();

    BOOST_CHECK(grid_equal(cgrid1, cgrid2));

    [[maybe_unused]] Opm::EclipseGrid grid = Opm::UgGridHelpers::createEclipseGrid( *cgrid1 , es1.getInputGrid( ) );
    auto actnum = Opm::UgGridHelpers::createACTNUM(*cgrid1);
    for (std::size_t g = 0; g < 300; g++)
        BOOST_CHECK_EQUAL(actnum[g], 1);
}
