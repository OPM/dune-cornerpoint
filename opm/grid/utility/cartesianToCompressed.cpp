/*
  Copyright 2018 SINTEF Digital, Mathematics & Cybernetics.
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


#include <opm/grid/utility/cartesianToCompressed.hpp>

namespace Opm
{

// Construct explicit mapping from local cartesian to active/compressed
// indices, either based on global_cell or as { 0, 1, 2, ....} if null.
// \param[in] num_cartesian_cells    The number of cartesian cells.
// \param[in] global_cell  Either null, or an array of size num_cells.
// \return                 A unordered_map containing the mapping from local cartesian
//                         to active/compressed (Attention! When LGRs, one key takes multiple values),
//                         or the map { {0, 0}, {1, 1}, ... , {num_cells - 1, num_cells - 1} }
//                         if global_cell was null.
std::unordered_multimap<int, int> cartesianToCompressed(const int num_cells,
                                                        const int* global_cell)
{
    std::unordered_multimap<int, int> retval;
    retval.reserve(num_cells);
    if (global_cell) {
        for (int i = 0; i < num_cells; ++i) {
            retval.emplace(global_cell[i], i);
        }
    } else {
        for (int i = 0; i < num_cells; ++i) {
            retval.emplace(i, i);
        }
    }
    return retval;
}


} // namespace Opm
