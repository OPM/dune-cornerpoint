/*
  Copyright 2018 SINTEF ICT, Applied Mathematics.
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

#include <config.h>

#define BOOST_TEST_MODULE mappingTest

#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/compressedToCartesian.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <unordered_map>
#include <iterator>

BOOST_AUTO_TEST_CASE(mapping)
{
    const std::vector<int> global_cell{0, 1, 2, 3, 5, 6, 8, 9};
    const int num_cells = global_cell.size(); // 8

    const std::unordered_multimap<int, int> cartesian_to_compressed =
        Opm::cartesianToCompressed(num_cells, global_cell.data());


    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(0), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(1), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(2), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(3), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(5), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(6), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(8), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(9), 1);

    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(0), 0);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(1), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(2), 2);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(3), 3);

    auto it5 = cartesian_to_compressed.equal_range(5);
    BOOST_CHECK_EQUAL(it5.first -> second, 4);
    auto it6 = cartesian_to_compressed.equal_range(6);
    BOOST_CHECK_EQUAL(it6.first -> second, 5);
    auto it8 = cartesian_to_compressed.equal_range(8);
    BOOST_CHECK_EQUAL(it8.first -> second, 6);
    auto it9 = cartesian_to_compressed.equal_range(9);
    BOOST_CHECK_EQUAL(it9.first -> second, 7);

    const int non_existing_index = 1829;
    auto it = cartesian_to_compressed.equal_range(non_existing_index);
    BOOST_CHECK(it.first == it.second);

    const std::vector<int> compressed_to_cartesian = Opm::compressedToCartesian(num_cells,  global_cell.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(compressed_to_cartesian.begin(), compressed_to_cartesian.end(),
                                  global_cell.begin(),             global_cell.end());
}

BOOST_AUTO_TEST_CASE(multiMapping)
{
    const std::vector<int> global_cell{0, 1, 2, 3,3, 3,3, 3,3, 5, 6,6,6,6, 6,6,6,6, 6,6,6,6, 8, 9}; // 24
    const int num_cells = global_cell.size();

    const std::unordered_multimap<int, int> cartesian_to_compressed =
        Opm::cartesianToCompressed(num_cells, global_cell.data());

    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(0), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(1), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(2), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(3), 6);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(5), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(6), 12);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(8), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.count(9), 1);

    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(0), 0);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(1), 1);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(2), 2);
    BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(3), 3);

    auto it3 = cartesian_to_compressed.equal_range(3);
    for (auto it = it3.first; it != it3.second; ++it)
    {
        BOOST_CHECK((it->second > 2) && (it->second < 9)); // it's unordered, so we cannot expect a specific position.
    }

    auto it5 = cartesian_to_compressed.equal_range(5);
    BOOST_CHECK_EQUAL(it5.first -> second, 9);

    auto it6 = cartesian_to_compressed.equal_range(6);
    for (auto it = it6.first; it != it6.second; ++it)
    {
        BOOST_CHECK((it->second > 9) && (it->second < 22)); // it's unordered, so we cannot expect a specific position.
    }

    auto it8 = cartesian_to_compressed.equal_range(8);
    BOOST_CHECK_EQUAL(it8.first -> second, 22);
    auto it9 = cartesian_to_compressed.equal_range(9);
    BOOST_CHECK_EQUAL(it9.first -> second, 23);

    const int non_existing_index = 1829;
    auto it = cartesian_to_compressed.equal_range(non_existing_index);
    BOOST_CHECK(it.first == it.second);

    const std::vector<int> compressed_to_cartesian = Opm::compressedToCartesian(num_cells,  global_cell.data());

    BOOST_CHECK_EQUAL_COLLECTIONS(compressed_to_cartesian.begin(), compressed_to_cartesian.end(),
                                  global_cell.begin(),             global_cell.end());
}

BOOST_AUTO_TEST_CASE(nullmapping)
{
    const int num_cells = 30;

    const std::unordered_multimap<int, int> cartesian_to_compressed =
        Opm::cartesianToCompressed(num_cells, nullptr);

    for (int i = 0; i < num_cells; ++i) {
        BOOST_CHECK_EQUAL(cartesian_to_compressed.bucket(i), i);
    }

    const std::vector<int> compressed_to_cartesian = Opm::compressedToCartesian(num_cells,  nullptr);

    for (int i = 0; i < num_cells; ++i) {
        BOOST_CHECK_EQUAL(compressed_to_cartesian[i], i);
    }
}
