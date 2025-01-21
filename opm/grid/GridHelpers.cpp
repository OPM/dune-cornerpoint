/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS
  Copyright 2015 NTNU

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

#include <opm/common/utility/ActiveGridCells.hpp>

#include <opm/grid/common/Volumes.hpp>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#endif

namespace Opm
{
namespace UgGridHelpers
{
int numCells(const UnstructuredGrid& grid)
{
    return grid.number_of_cells;
}

int numFaces(const UnstructuredGrid& grid)
{
    return grid.number_of_faces;
}
int dimensions(const UnstructuredGrid& grid)
{
    return grid.dimensions;
}
int numCellFaces(const UnstructuredGrid& grid)
{
    return grid.cell_facepos[grid.number_of_cells];
}

const int* globalCell(const UnstructuredGrid& grid)
{
    return grid.global_cell;
}

const int* cartDims(const UnstructuredGrid& grid)
{
    return grid.cartdims;
}

#if HAVE_ECL_INPUT
std::vector<int> createACTNUM(const UnstructuredGrid& grid) {
    const int* dims = cartDims(grid);
    return ActiveGridCells(dims[0], dims[1], dims[2], globalCell(grid), numCells(grid)).actNum();
}
#endif

const double* beginCellCentroids(const UnstructuredGrid& grid)
{
    return grid.cell_centroids;
}

double cellCenterDepth(const UnstructuredGrid& grid, int cell_index)
{
    // This method is an alternative to the method cellCentroidCoordinate(...) below.
    // The cell center depth is computed as a raw average of cell corner depths.
    // For cornerpoint grids, this is likely to give slightly different depths that seem
    // to agree with eclipse.
    assert(grid.dimensions == 3);
    const int nd = 3; // Assuming 3-dimensional grid ...
    const int nv = 8; // Assuming 2*4 vertices ...
    double zz = 0.0;
    // Traverse the bottom and top cell-face
    for (grid_size_t i = grid.cell_facepos[cell_index+1] - 2; i < grid.cell_facepos[cell_index+1]; ++i) {
        // Traverse the vertices associated with each face
        assert(grid.face_nodepos[grid.cell_faces[i]+1] - grid.face_nodepos[grid.cell_faces[i]] == nv/2);
        for (grid_size_t j = grid.face_nodepos[grid.cell_faces[i]]; j < grid.face_nodepos[grid.cell_faces[i]+1]; ++j) {
            zz += (grid.node_coordinates+nd*(grid.face_nodes[j]))[nd-1];
        }
    }
    return zz/nv;
}

Dune::FieldVector<double,3> faceCenterEcl(const UnstructuredGrid& grid, int cell_index, int face_tag)
{
    // This method is an alternative to the method faceCentroid(...) below.
    // The face center is computed as a raw average of cell corners.
    // For faulted cells this gives different results then average of face nodes
    // that seems to agree more with eclipse.
    // This assumes that the top and bottom face nodes are ordered
    // 0--1
    // |  |
    // 3--2

    assert(grid.dimensions == 3);
    const int nd = 3; // Assuming 3-dimensional grid ...
    const int nv = 4; // Assuming 4 vertices ...
    Dune::FieldVector<double,3> center(0.0);
    //Vector center(0.0);
    // Traverse the bottom and top cell-face
    for (grid_size_t i = grid.cell_facepos[cell_index+1] - 2; i < grid.cell_facepos[cell_index+1]; ++i) {
        // Traverse the vertices associated with each face
        assert(grid.face_nodepos[grid.cell_faces[i]+1] - grid.face_nodepos[grid.cell_faces[i]] == nv);

        int start = grid.face_nodepos[grid.cell_faces[i]];

        // pick the right nodes. See order assumption above
        switch(face_tag) {
        case 0: {
            for (int indx = 0; indx < nd; ++indx) {
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+3]))[indx];
            }
        }
            break;
        case 1: {
            for (int indx = 0; indx < nd; ++indx) {
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[indx];
            }
        }
            break;
        case 2: {
            for (int indx = 0; indx < nd; ++indx) {


                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[indx];

            }
        }
            break;
        case 3: {
            for (int indx = 0; indx < nd; ++indx) {

                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+3]))[indx];

            }
        }
            break;
        case 4: {
            if (i == grid.cell_facepos[cell_index+1] - 2) {
                for (int indx = 0; indx < nd; ++indx) {
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+3]))[indx];
                }
            }
        }
            break;
        case 5: {
            if (i == grid.cell_facepos[cell_index+1] - 1) {
                for (int indx = 0; indx < nd; ++indx) {
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[indx];
                center[indx] += (grid.node_coordinates+nd*(grid.face_nodes[start+3]))[indx];
                }
            }
        }
            break;
        }
    }
    for (int indx = 0; indx < nd; ++indx) {
        center[indx] /= nv;
    }

    return center;
}


Dune::FieldVector<double,3> faceAreaNormalEcl(const UnstructuredGrid& grid, int face_index)
{
    // This method is an alternative to the method faceNormal(...) below.
    // The face Normal area is computed based on the face corners without introducing
    // a center point.
    // For cornerpoint grids, this is likely to give slightly different depths that seem
    // to agree with eclipse.
    assert(grid.dimensions == 3);
    const int nd = 3; // Assuming 3-dimensional grid ...
    const int nv = grid.face_nodepos[face_index+1] - grid.face_nodepos[face_index];
    const int start = grid.face_nodepos[face_index];

    typedef Dune::FieldVector<double,3> Vector;

    switch (nv)
    {
    case 0:
    case 1:
    case 2:
        {
            return Vector( 0 );
        }
        break;
    case 3:
        {
        Vector a, b;
        for (int i = 0; i < 3; ++i) {
            a[i] = (grid.node_coordinates+nd*(grid.face_nodes[start] ))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[i];
            b[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[i];
        }
        Vector areaNormal  = cross( a, b);
        areaNormal *= 0.5;
        return areaNormal;
        }
        break;
    case 4:
        {
        Vector a, b;
        for (int i = 0; i < 3; ++i) {
            a[i] = (grid.node_coordinates+nd*(grid.face_nodes[start] ))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+2]))[i];
            b[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+1]))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+3]))[i];
        }
        Vector areaNormal  = Dune::cross( a, b);
        areaNormal *= 0.5;
        return areaNormal;
        }
        break;
    default:
        {
            int h = (nv - 1)/2;
            int k = (nv % 2) ? 0 : nv - 1;

            Vector areaNormal ( 0 );
            Vector a, b;
            // First quads
            for (int j = 1; j < h; ++j)
            {
                for (int i = 0; i < 3; ++i) {
                    a[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+2*j] ))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start]))[i];
                    b[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+2*j + 1]))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+ 2*j-1]))[i];
                }
                areaNormal += cross( a , b ) ;
            }

            // Last triangle or quad
            for (int i = 0; i < 3; ++i) {
                a[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+2*h] ))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start]))[i];
                b[i] = (grid.node_coordinates+nd*(grid.face_nodes[start+k]))[i] - (grid.node_coordinates+nd*(grid.face_nodes[start+ 2*h-1]))[i];
            }
            areaNormal += cross( a , b ) ;
            areaNormal *= 0.5;
            return areaNormal;
        }
    }

}

double cellCentroidCoordinate(const UnstructuredGrid& grid, int cell_index,
                                 int coordinate)
{
    return grid.cell_centroids[grid.dimensions*cell_index+coordinate];
}

const double*
cellCentroid(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_centroids+(cell_index*grid.dimensions);
}

const double* beginCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes;
}
const double* endCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes+numCells(grid);
}

const double* beginFaceCentroids(const UnstructuredGrid& grid)
{
    return grid.face_centroids;
}

const double* faceCentroid(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_centroids+face_index*grid.dimensions;
}

const double* faceNormal(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_normals+face_index*grid.dimensions;
}

double faceArea(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_areas[face_index];
}

SparseTableView cell2Faces(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.cell_faces, grid.cell_facepos, numCells(grid));
}

SparseTableView face2Vertices(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.face_nodes, grid.face_nodepos, numFaces(grid));
}

const double* vertexCoordinates(const UnstructuredGrid& grid, int index)
{
    return grid.node_coordinates+dimensions(grid)*index;
}

double cellVolume(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_volumes[cell_index];
}

FaceCellTraits<UnstructuredGrid>::Type faceCells(const UnstructuredGrid& grid)
{
    return FaceCellsProxy(grid);
}


#if HAVE_ECL_INPUT
Opm::EclipseGrid createEclipseGrid(const UnstructuredGrid& grid, const Opm::EclipseGrid& inputGrid ) {
    const int * dims = UgGridHelpers::cartDims( grid );

    if ((inputGrid.getNX( ) == static_cast<size_t>(dims[0])) &&
        (inputGrid.getNY( ) == static_cast<size_t>(dims[1])) &&
        (inputGrid.getNZ( ) == static_cast<size_t>(dims[2]))) {
        std::vector<int> updatedACTNUM;
        const int* global_cell = UgGridHelpers::globalCell( grid );

        if (global_cell) {
            updatedACTNUM.assign( inputGrid.getCartesianSize( ) , 0 );
            for (int c = 0; c < numCells( grid ); c++) {
                updatedACTNUM[global_cell[c]] = 1;
            }
        }

        return Opm::EclipseGrid( inputGrid, grid.zcorn, updatedACTNUM );
    } else {
        throw std::invalid_argument("Size mismatch - dimensions of inputGrid argument and current UnstructuredGrid instance disagree");
    }
}
#endif

}
}
