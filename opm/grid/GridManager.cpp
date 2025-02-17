/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/GridManager.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/grid/cart_grid.h>
#include <opm/grid/cornerpoint_grid.h>
#include <opm/grid/MinpvProcessor.hpp>
#include <opm/common/ErrorMacros.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#endif

#include <array>
#include <algorithm>
#include <numeric>

namespace Opm
{

#if HAVE_ECL_INPUT
    /// Construct a 3d corner-point grid from a deck.
    GridManager::GridManager(const Opm::EclipseGrid& inputGrid,
                             const bool edge_conformal)
        : ug_(0)
    {
        initFromEclipseGrid(inputGrid, std::vector<double>(), edge_conformal);
    }


    GridManager::GridManager(const Opm::EclipseGrid& inputGrid,
                             const std::vector<double>& poreVolumes,
                             const bool edge_conformal)
        : ug_(0)
    {
        initFromEclipseGrid(inputGrid, poreVolumes, edge_conformal);
    }
#endif


    /// Construct a 2d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny)
    {
        ug_ = create_grid_cart2d(nx, ny, 1.0, 1.0);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }

    GridManager::GridManager(int nx, int ny,double dx, double dy)
    {
        ug_ = create_grid_cart2d(nx, ny, dx, dy);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }


    /// Construct a 3d cartesian grid with cells of unit size.
    GridManager::GridManager(int nx, int ny, int nz)
    {
        ug_ = create_grid_cart3d(nx, ny, nz);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }




    /// Construct a 3d cartesian grid with cells of size [dx, dy, dz].
    GridManager::GridManager(int nx, int ny, int nz,
                             double dx, double dy, double dz)
    {
        ug_ = create_grid_hexa3d(nx, ny, nz, dx, dy, dz);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }
    }




    /// Construct a grid from an input file.
    /// The file format used is currently undocumented,
    /// and is therefore only suited for internal use.
    GridManager::GridManager(const std::string& input_filename)
    {
        ug_ = read_grid(input_filename.c_str());
        if (!ug_) {
            OPM_THROW(std::runtime_error,
                      "Failed to read grid from file " + input_filename);
        }
    }

    /// Destructor.
    GridManager::~GridManager()
    {
        destroy_grid(ug_);
    }




    /// Access the managed UnstructuredGrid.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const UnstructuredGrid* GridManager::c_grid() const
    {
        return ug_;
    }



#if HAVE_ECL_INPUT
    // Construct corner-point grid from EclipseGrid.
    void GridManager::initFromEclipseGrid(const Opm::EclipseGrid& inputGrid,
                                          const std::vector<double>& poreVolumes,
                                          const bool edge_conformal)
    {
        grdecl g{};

        g.dims[0] = inputGrid.getNX();
        g.dims[1] = inputGrid.getNY();
        g.dims[2] = inputGrid.getNZ();

        std::vector<double> coord = inputGrid.getCOORD( );
        std::vector<double> zcorn = inputGrid.getZCORN( );
        std::vector<int> actnum = inputGrid.getACTNUM(  );

        g.coord = coord.data();
        g.zcorn = zcorn.data();
        g.actnum = actnum.data();

        const double z_tolerance = inputGrid.isPinchActive() ? inputGrid.getPinchThresholdThickness() : 0.0;

        if (!poreVolumes.empty() && (inputGrid.getMinpvMode() != MinpvMode::Inactive)) {
            MinpvProcessor mp(g.dims[0], g.dims[1], g.dims[2]);
            const std::vector<double>& minpvv  = inputGrid.getMinpvVector();
            const size_t cartGridSize = g.dims[0] * g.dims[1] * g.dims[2];
            std::vector<double> thickness(cartGridSize);
            for (size_t i = 0; i < cartGridSize; ++i) {
                thickness[i] = inputGrid.getCellThickness(i);
            }

            // The legacy code only supports the opmfil option
            bool opmfil = true; //inputGrid.getMinpvMode() == MinpvMode::OpmFIL;
            mp.process(thickness, z_tolerance, inputGrid.getPinchMaxEmptyGap(), poreVolumes, minpvv,
                       actnum, opmfil, zcorn.data());
        }

        ug_ = create_grid_cornerpoint(&g, z_tolerance, edge_conformal);
        if (!ug_) {
            OPM_THROW(std::runtime_error, "Failed to construct grid.");
        }

        attach_zcorn_copy( ug_ , zcorn.data() );

    }

#endif






} // namespace Opm
