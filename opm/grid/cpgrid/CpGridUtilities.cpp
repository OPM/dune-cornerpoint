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

#include <opm/grid/cpgrid/CpGridUtilities.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{

std::pair<std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
lgrIJK(const Dune::CpGrid& grid, const std::string& lgr_name)
{
    // Check if lgr_name exists in lgr_names_
    const auto& lgr_names = grid.getLgrNameToLevel();
    auto it = lgr_names.find(lgr_name);
    if (it == lgr_names.end()) {
        OPM_THROW(std::runtime_error, "LGR name not found: " + lgr_name);
    }

    const auto level = it->second;
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapper(grid);
    const auto& levelView = grid.levelGridView(level);
    const auto numCells = levelView.size(0);

    std::vector<std::array<int, 3>> lgrIJK(numCells);
    std::unordered_map<int, int> lgrCartesianIdxToCellIdx;
    lgrCartesianIdxToCellIdx.reserve(numCells);

    // Iterate over (active) elements in the grid and populate the structures
    for (const auto& element : Dune::elements(grid.levelGridView(level))) {
        std::array<int, 3> ijk;
        levelCartMapper.cartesianCoordinate(element.index(), ijk, level);

        const int cellIndex = element.index();
        const int cartesianIdx = element.getLevelCartesianIdx();

        lgrIJK[cellIndex] = ijk;
        lgrCartesianIdxToCellIdx[cartesianIdx] = cellIndex;
    }

    return std::make_pair(lgrCartesianIdxToCellIdx, lgrIJK);
}

std::pair<std::vector<double>, std::vector<double>>
lgrCOORDandZCORN(const Dune::CpGrid& grid,
                 int level,
                 const std::unordered_map<int, int>&  lgrCartesianIdxToCellIdx,
                 const std::vector<std::array<int, 3>>& lgrIJK)
{
    const auto& levelGrid = *(grid.currentData()[level]);

    // Check not all cells are inactive
    const auto numCells = levelGrid.size(0);
    if (numCells == 0) {
        OPM_THROW(std::logic_error, "LGR in level " + std::to_string(level) + " has no active cells.\n");
    }

    // LGR dimensions
    const auto& lgr_dim = grid.currentData()[level]->logicalCartesianSize();
    const int nx = lgr_dim[0];
    const int ny = lgr_dim[1];
    const int nz = lgr_dim[2];

    // Initialize all pillars as inactive (setting COORD values to std::numeric_limits<double>::max()).
    std::vector<double> lgrCOORD(6*(nx+1)*(ny+1), std::numeric_limits<double>::max());

    // Initialize all ZCORN as inactive (setting values to std::numeric_limits<double>::max()).
    std::vector<double> lgrZCORN(8*nx*ny*nz, std::numeric_limits<double>::max());

    // Map to determine min and max k per cell column (i, j) (min/max_k = 0, ..., nz-1).
    // Initialized as {nz, -1} to detect inactive cell columns.
    std::vector<std::array<int,2>> minMaxPerCellPillar(nx*ny, {nz, -1});

   
    for (const auto& ijk : lgrIJK) {

         // Compute the bottom and top k per cell pillar (i, j).
        int cell_pillar_idx = ijk[1] * nx + ijk[0];
        auto& minMax = minMaxPerCellPillar[cell_pillar_idx];

        minMax[0] = std::min(ijk[2], minMax[0]);
        minMax[1] = std::max(ijk[2], minMax[1]);
    }

    for (const auto& elem : elements(grid.levelGridView(level))) {
        const auto& elemIJK = lgrIJK[elem.index()];

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

        int zcorn_top_00_idx = ((nz-1-elemIJK[2])*8*nx*ny) + (elemIJK[1]*4*nx) + (2*elemIJK[0]); // assoc. w. elem corner 4

        // Bottom indices
        int zcorn_top_10_idx = zcorn_top_00_idx + 1;  // assoc. w. elem corner 5
        int zcorn_top_01_idx = zcorn_top_00_idx + (2*nx);  // assoc. w. elem corner 6
        int zcorn_top_11_idx = zcorn_top_01_idx + 1; // assoc. w. elem corner 7

        // Top indices
        int zcorn_bottom_00_idx = zcorn_top_00_idx + (4*nx*ny); // assoc. w. elem corner 0
        int zcorn_bottom_10_idx = zcorn_bottom_00_idx + 1;  // assoc. w. elem corner 1
        int zcorn_bottom_01_idx = zcorn_bottom_00_idx + (2*nx); // assoc. w. elem corner 2
        int zcorn_bottom_11_idx = zcorn_bottom_01_idx + 1;  // assoc. w. elem corner

        // Note: zcorn_idx + 1 moves to the next position along the x-axis (i+1, j, k)
        //       zcorn_idx + (2*nx) moves to the next position along the y-axis (i, j+1, k)
        //       zcorn_idx + (4*nx*ny) moves to the next position along the z-axis (i,j, k+1)

        // Assign ZCORN values
        lgrZCORN[zcorn_top_00_idx] = elem.subEntity<3>(4).geometry().center()[2];
        lgrZCORN[zcorn_top_10_idx] = elem.subEntity<3>(5).geometry().center()[2];
        lgrZCORN[zcorn_top_01_idx] = elem.subEntity<3>(6).geometry().center()[2];
        lgrZCORN[zcorn_top_11_idx] = elem.subEntity<3>(7).geometry().center()[2];

        lgrZCORN[zcorn_bottom_00_idx] = elem.subEntity<3>(0).geometry().center()[2];
        lgrZCORN[zcorn_bottom_10_idx] = elem.subEntity<3>(1).geometry().center()[2];
        lgrZCORN[zcorn_bottom_01_idx] = elem.subEntity<3>(2).geometry().center()[2];
        lgrZCORN[zcorn_bottom_11_idx] = elem.subEntity<3>(3).geometry().center()[2];
    }

    // Rewrite values for active pillars
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int cell_pillar_idx = (j*nx) + i;

            // Get min/max k for pillar at (i,j)
            const auto& [bottom_k, top_k] = minMaxPerCellPillar[cell_pillar_idx];

            if ( bottom_k == nz ) {
                 continue; // no active pillar at (i,j)
            }

            const auto bottom_lgr_cartesian_idx = (bottom_k*nx*ny) + cell_pillar_idx;
            const auto top_lgr_cartesian_idx = (top_k*nx*ny) + cell_pillar_idx;

            const auto& bottomElemIdx = lgrCartesianIdxToCellIdx.at(bottom_lgr_cartesian_idx);
            const auto& topElemIdx = lgrCartesianIdxToCellIdx.at(top_lgr_cartesian_idx);

            const auto& bottomElem = Dune::cpgrid::Entity<0>(levelGrid, bottomElemIdx, true);
            const auto& topElem = Dune::cpgrid::Entity<0>(levelGrid, topElemIdx, true);

            Opm::processPillars(i,j, nx, topElem, bottomElem, lgrCOORD);
        }
    }
    return std::make_pair(lgrCOORD, lgrZCORN);
}

void setPillarCoordinates(int i, int j, int nx,
                          int topCorner, int bottomCorner, int positionIdx,
                          const Dune::cpgrid::Entity<0>& topElem,
                          const Dune::cpgrid::Entity<0>& bottomElem,
                          std::vector<double>& lgrCOORD)
{
    // positionIdx (0, 1, 2, or 3) is used to distinguish the 4 corner pillars in a cell column.
    //
    // Corner pillar position mapping:
    //
    //   positionIdx   corresponding pillar   positionIdx / 2   positionIdx % 2
    //   ---------------------------------------------------------------------
    //       0          (i, j)                     0                  0
    //       1          (i+1, j)                   0                  1
    //       2          (i, j+1)                   1                  0
    //       3          (i+1, j+1)                 1                  1
    //
    // - positionIdx / 2 determines the position at the y-axis (0 for j, 1 for j+1).
    // - positionIdx % 2 determines the position at the x-axis (0 for i, 1 for i+1).

    const int pillar = ((j + positionIdx / 2) *6* (nx + 1)) + 6*(i + positionIdx % 2);

    // Top pillar's COORD values
    const auto& top_point = topElem.subEntity<3>(topCorner).geometry().center();
    lgrCOORD[pillar] = top_point[0];
    lgrCOORD[pillar+1] = top_point[1];
    lgrCOORD[pillar+2] = top_point[2];

    // Bottom pillar's COORD values
    const auto& bottom_point = bottomElem.subEntity<3>(bottomCorner).geometry().center();
    lgrCOORD[pillar+3] = bottom_point[0];
    lgrCOORD[pillar+4] = bottom_point[1];
    lgrCOORD[pillar+5] = bottom_point[2];
}

void processPillars(int i, int j, int nx,
                    const Dune::cpgrid::Entity<0>& topElem,
                    const Dune::cpgrid::Entity<0>& bottomElem,
                    std::vector<double>& lgrCOORD)
{
    // Recall that a cell has 8 corners:
    //        6 --- 7
    //       /     /   TOP FACE
    //      4 --- 5
    //        2 --- 3
    //       /     /   BOTTOM FACE
    //      0 --- 1

    // To take into account inactive cells, consider for each (i,j) column of cells, 4 pillars:
    // (i,j)     pillar associated with bottom element corner 0 and top element corner 4
    // (i+1,j)   pillar associated with bottom element corner 1 and top element corner 5
    // (i,j+1)   pillar associated with bottom element corner 2 and top element corner 6
    // (i+1,j+1) pillar associated with bottom element corner 3 and top element corner 7
    setPillarCoordinates(i, j, nx, 4 /*topCorner*/, 0 /*bottomCorner*/,  0 /*to select pillar (i,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 5 /*topCorner*/, 1 /*bottomCorner*/,  1 /*to select pillar (i+1,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 6 /*topCorner*/, 2 /*bottomCorner*/,  2 /*to select pillar (i,j+1)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 7 /*topCorner*/, 3 /*bottomCorner*/,  3 /*to select pillar (i+1,j+1)*/, topElem, bottomElem, lgrCOORD);
}

} // namespace Opm
