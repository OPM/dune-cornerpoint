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

#ifndef OPM_CPGRIDUTILITIES_HEADER_INCLUDED
#define OPM_CPGRIDUTILITIES_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>

namespace Opm
{

/// @brief Retrieves Cartesian indices for a specified Local Grid Refinement (LGR) level in a Dune::CpGrid.
///
/// This function extracts the mapping between active cell indices and Cartesian indices for a given LGR name.
/// If the specified name is not found, an exception is thrown.
///
/// @param [in] grid The Dune::CpGrid
/// @param [in] lgr_name The name of the LGR whose indices are to be retrieved.
/// @return A pair containing:
///   - A std::unordered_map<int, int> mapping Cartesian indices back to cell indices (handles inactive parent cells).
///   - A std::vector<std::array<int, 3>> storing the (i, j, k) Cartesian coordinates for active cells.
std::pair<std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
lgrIJK(const Dune::CpGrid& grid, const std::string& lgr_name);

/// @brief Extracts the COORD and ZCORN values for the LGR (Local Grid Refinement) block.
///
/// COORD: This retrieves a vector of std::array<double, 6>, where each element represents
/// the coordinate values of a pillar in the LGR block. The values correspond to the
/// COORD keyword, which defines the corner-point geometry of the grid.
/// ZCORN: It is initialized to inactive values, i.e. set to std::numeric_limits<double>::max(),
/// and then calculates the bottom and top ZCORN values for each cell in the LGR based on its
/// geometry. It determines the active cell layers and assigns ZCORN values accordingly.
///
/// Special Cases:
/// - If a pillar within the LGR block is "inactive", its COORD values are set to
/// std::numeric_limits<double>::max() to indicate the inactive status.
/// - If all pillars are "inactive", an exception is thrown.
///
/// @param [in] grid
/// @param [in] level The refinement level of the LGR block.
/// @param [in] lgrCartesianIdxToCellIdx Mapping from LGR Cartesian indices to level cell indices.
/// @param [in] lgrIJK The IJK indices corresponding to the LGR block active cells.
/// @return A pair containing:
///    - A std::vector<double> storing the coordinate values for each pillar.
///    - A std::vector<double> storing the depth of each corner point of the LGR pillars.
std::pair<std::vector<double>, std::vector<double>>
lgrCOORDandZCORN(const Dune::CpGrid& grid,
                 int level,
                 const std::unordered_map<int, int>& lgrCartesianIdxToCellIdx,
                 const std::vector<std::array<int, 3>>& lgrIJK);

/// @brief Sets the coordinates for a pillar.
///
/// This function calculates the pillar index based on the given (i, j) position
/// and assigns the top and bottom coordinates from the corresponding element corners.
///
/// @param [in] i Column index of the pillar.
/// @param [in] j Row index of the pillar.
/// @param [in] nx Number of elements in the x-direction.
/// @param [in] topCorner Index of the top element corner associated with the pillar.
/// @param [in] bottomCorner Index of the bottom element corner associated with the pillar.
/// @param [in] positionIdx To determine the correct pillar index (0-3).
/// @param [in] topElem Reference to the top element.
/// @param [in] bottomElem Reference to the bottom element.
/// @param [out] lgrCOORD Reference to the array storing the pillar coordinates.
void setPillarCoordinates(int i, int j, int nx,
                          int topCorner, int bottomCorner,
                          int positionIdx,
                          const Dune::cpgrid::Entity<0>& topElem,
                          const Dune::cpgrid::Entity<0>& bottomElem,
                          std::vector<double>& lgrCOORD);


/// @brief Processes and sets the coordinates for all four pillars of a given "(i,j) column of cells".
///
/// This function calls setPillarCoordinates for each of the four pillars (i,j), (i+1,j), (i,j+1),
/// and (i+1,j+1), assigning the correct top and bottom coordinates from the given elements.
///
/// @param [in] i Column index of the current element.
/// @param [in] j Row index of the current element.
/// @param [in] nx Number of elements in the x-direction.
/// @param [in] topElem Reference to the top element.
/// @param [in] bottomElem Reference to the bottom element.
/// @param [out] lgrCOORD Reference to the array storing the pillar coordinates.
void processPillars(int i, int j, int nx,
                    const Dune::cpgrid::Entity<0>& topElem,
                    const Dune::cpgrid::Entity<0>& bottomElem,
                    std::vector<double>& lgrCOORD);
} // namespace Opm

#endif // OPM_CPGRIDUTILITIES_HEADER_INCLUDED
