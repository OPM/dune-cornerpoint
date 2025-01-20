//===========================================================================
//
// File: test_lookupdata_polyhedral.cpp
//
// Created: Monday 24.07.2023 14:30:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2023 Equinor ASA.

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
#include <type_traits>

#define BOOST_TEST_MODULE LookUpDataPolyhedralGridTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/LookUpData.hh>

#include <dune/common/unused.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/grid/cpgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/levelcartesianindexmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/DeckSection.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckKeyword.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>

#include <sstream>
#include <stdexcept>
#include <iostream>

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

void lookup_check(const Dune::PolyhedralGrid<3,3>& grid)
{
    std::vector<int> fake_feature(grid.size(0), 0);
    std::iota(fake_feature.begin(), fake_feature.end(), 3);

    std::vector<double> fake_feature_double(grid.size(0), 0.);
    std::iota(fake_feature_double.begin(), fake_feature_double.end(), .5);
    
    const auto leaf_view = grid.leafGridView();
    using GridView = std::remove_cv_t< typename std::remove_reference<decltype(grid.leafGridView())>::type>;
    // LookUpData
    const Opm::LookUpData<Dune::PolyhedralGrid<3,3>, GridView> lookUpData(leaf_view, false);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::PolyhedralGrid<3,3>> cartMapper(grid);
    const Opm::LevelCartesianIndexMapper< Dune::PolyhedralGrid<3,3> > levelCartMapp(cartMapper);
    const Opm::LookUpCartesianData<Dune::PolyhedralGrid<3,3>, GridView> lookUpCartesianData(leaf_view, cartMapper, false);
    // Mapper
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    for (const auto& elem : elements(leaf_view)) {
        // Search via Entity/Element
        const auto featureInElem = lookUpData(elem, fake_feature);
        const auto featureInElemDouble = lookUpData(elem, fake_feature_double);
        const auto featureInElemCartesian = lookUpCartesianData(elem, fake_feature);
        const auto featureInElemDoubleCartesian = lookUpCartesianData(elem,fake_feature_double);
        BOOST_CHECK(featureInElem == lookUpData.getFieldPropIdx(elem) +3);
        BOOST_CHECK(featureInElemDouble == lookUpData.getFieldPropIdx(elem) +.5);
        BOOST_CHECK(featureInElemCartesian == lookUpData.getFieldPropIdx(elem) +3);
        BOOST_CHECK(featureInElemDoubleCartesian == lookUpData.getFieldPropIdx(elem) +.5);
        // Search via INDEX
        const auto idx = mapper.index(elem);
        using IndexType = decltype(idx);
        using ElementType = std::remove_reference_t<decltype(elem)>;
        const auto featureInElemIDX = lookUpData(idx, fake_feature);
        const auto featureInElemDoubleIDX = lookUpData(idx, fake_feature_double);
        const auto featureInElemCartesianIDX = lookUpCartesianData(idx, fake_feature);
        const auto featureInElemDoubleCartesianIDX = lookUpCartesianData(idx, fake_feature_double);
        BOOST_CHECK(featureInElemIDX == (lookUpData.getFieldPropIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx))+3);
        BOOST_CHECK(featureInElemDoubleIDX == (lookUpData.getFieldPropIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx)) +.5);
        BOOST_CHECK(featureInElemCartesianIDX == (lookUpData.getFieldPropIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx)) +3);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == (lookUpData.getFieldPropIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx)) +.5);
        BOOST_CHECK(idx == (lookUpData.getFieldPropIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx)));
        BOOST_CHECK(featureInElemIDX == featureInElem);
        BOOST_CHECK(featureInElemDoubleIDX == featureInElemDouble);
        BOOST_CHECK(featureInElemCartesianIDX == featureInElemCartesian);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == featureInElemDoubleCartesian);
        // Extra checks related to Cartesian Index
        const auto cartIdx = cartMapper.cartesianIndex(idx);
        BOOST_CHECK(cartIdx == (lookUpCartesianData.getFieldPropCartesianIdx<ElementType, Dune::PolyhedralGrid<3,3>>(elem)));
        BOOST_CHECK(cartIdx == (lookUpCartesianData.getFieldPropCartesianIdx<IndexType, Dune::PolyhedralGrid<3,3>>(idx)));
        // Extra checks related to Cartesian Coordinate
        std::array<int,3> ijk;
        cartMapper.cartesianCoordinate(idx, ijk);
        std::array<int,3> ijkLevel;
        levelCartMapp.cartesianCoordinate(idx, ijkLevel, 0);
        BOOST_CHECK(ijk == ijkLevel);
        // Throw for level > 0 (Local grid refinement not supported for Polyhedral Grid)
        std::array<int,3> ijkThrow;
        BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(idx, ijkThrow, 4), std::invalid_argument);
        BOOST_CHECK_THROW(levelCartMapp.cartesianCoordinate(idx, ijkThrow, -3), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(PolyGridFromEcl)
{
#if HAVE_ECL_INPUT
    const char *deckString =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "4 4 4 /\n"
        "GRID\n"
        "DXV\n"
        "4*1 /\n"
        "DYV\n"
        "4*1 /\n"
        "DZ\n"
        "16*1 /\n"
        "TOPS\n"
        "16*100.0 /\n";

    Opm::Parser parser;
    auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclgrid( deck);
    std::vector<double> porv;

    Dune::PolyhedralGrid<3,3> grid(eclgrid, porv);
    lookup_check(grid);
#endif
}

void fieldProp_check(const Dune::PolyhedralGrid<3,3>& grid, Opm::EclipseGrid eclGrid, std::string deck_string)
{

    Opm::Deck deck = Opm::Parser{}.parseString(deck_string);
    Opm::FieldPropsManager fpm(deck, Opm::Phases{true, true, true}, eclGrid, Opm::TableManager());
    const auto& poro = fpm.get_double("PORO");
    const auto& eqlnum = fpm.get_int("EQLNUM");

    // LookUpData
    auto leaf_view = grid.leafGridView();
    using GridView = std::remove_cv_t< typename std::remove_reference<decltype(leaf_view)>>::type;
    const Opm::LookUpData<Dune::PolyhedralGrid<3,3>, GridView> lookUpData(leaf_view);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::PolyhedralGrid<3,3>> cartMapper(grid);
    const Opm::LookUpCartesianData<Dune::PolyhedralGrid<3,3>, GridView> lookUpCartesianData(leaf_view, cartMapper);
    // Element mapper
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(leaf_view, Dune::mcmgElementLayout());

    const auto& poroOnLeaf = lookUpData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");
    const auto& poroOnLeafCart = lookUpCartesianData.assignFieldPropsDoubleOnLeaf(fpm, "PORO");

    const auto& eqlnumOnLeaf = lookUpData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM", true);
    const auto& eqlnumOnLeafCart = lookUpCartesianData.assignFieldPropsIntOnLeaf<int>(fpm, "EQLNUM", true);

    for (const auto& elem : elements(leaf_view))
    {
        const auto elemIdx = mapper.index(elem);
        // PORO
        BOOST_CHECK_EQUAL(poro[elemIdx], lookUpData.fieldPropDouble(fpm, "PORO", elem));
        BOOST_CHECK_EQUAL(poro[elemIdx], lookUpData.fieldPropDouble<int>(fpm, "PORO", elemIdx));
        BOOST_CHECK_EQUAL(poro[elemIdx], lookUpCartesianData.fieldPropDouble(fpm, "PORO", elem));
        BOOST_CHECK_EQUAL(poro[elemIdx], lookUpCartesianData.fieldPropDouble(fpm, "PORO", elemIdx));
        BOOST_CHECK_EQUAL(poro[elemIdx], poroOnLeaf[elemIdx]);
        BOOST_CHECK_EQUAL(poro[elemIdx], poroOnLeafCart[elemIdx]);
        // EQLNUM
        BOOST_CHECK_EQUAL(eqlnum[elemIdx], lookUpData.fieldPropInt(fpm, "EQLNUM", elem));
        BOOST_CHECK_EQUAL(eqlnum[elemIdx], lookUpData.fieldPropInt<int>(fpm, "EQLNUM", elemIdx));
        BOOST_CHECK_EQUAL(eqlnum[elemIdx], lookUpCartesianData.fieldPropInt(fpm, "EQLNUM", elem));
        BOOST_CHECK_EQUAL(eqlnum[elemIdx], lookUpCartesianData.fieldPropInt(fpm, "EQLNUM", elemIdx));
        BOOST_CHECK_EQUAL(eqlnum[elemIdx]-true, eqlnumOnLeaf[elemIdx]);
        BOOST_CHECK_EQUAL(eqlnum[elemIdx]-true, eqlnumOnLeafCart[elemIdx]);
    }
}


BOOST_AUTO_TEST_CASE(fieldProp) {
    std::string deck_string = R"(
GRID

PORO
1.0 2.0 3.0 4.0 5.0 6.0
/
REGIONS

EQLNUM
1 2 3 4 5 6
/)";

    std::vector<int> actnum1 = {1,1,1,1,1,1};
    Opm::EclipseGrid eclGrid(3,2,1);
    eclGrid.resetACTNUM(actnum1);
    std::vector<double> porv;

    Dune::PolyhedralGrid<3,3> grid(eclGrid, porv);
    fieldProp_check(grid, eclGrid, deck_string);
}

