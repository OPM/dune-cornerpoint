//===========================================================================
//
// File: BuildCpGrid.cpp
//
// Created: Wed Aug 26 12:30:54 2009
//
// Author(s): B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Porous Media project  (OPM).

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

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <dune/grid/CpGrid.hpp>
#include <dune/common/mpihelper.hh>
#include "BuildCpGrid.hpp"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    Opm::parameter::ParameterGroup param(argc, argv);
    int nx = param.getDefault<int>("nx", 100);
    int ny = param.getDefault<int>("ny",  60);
    int nz = param.getDefault<int>("nz",  15);

    double hx = param.getDefault<double>("hx", 1.0);
    double hy = param.getDefault<double>("hy", 1.0);
    double hz = param.getDefault<double>("hz", 0.5);

    double drop = param.getDefault<double>("drop", 0.0);
    double z_tol = param.getDefault<double>("z_tol", 0.0);

    // Dune::SimpleFault model(nx, ny, nz, hx, hy, hz, drop);
    Dune::SlopingFault model(nx, ny, nz, hx, hy, hz, drop);
    Dune::CpGrid grid;

    model.build(grid, z_tol);

    return 0;
}
