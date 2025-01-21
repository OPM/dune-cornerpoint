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
#include <config.h>
#include <opm/grid/CellQuadrature.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/common/ErrorMacros.hpp>
#include <algorithm>
#include <cmath>

namespace {

/// Calculates the determinant of a 3 x 3 matrix, represented as
/// three three-dimensional arrays.
inline double determinantOf(const double* a0,
                            const double* a1,
                            const double* a2)
{
    return
        a0[0] * (a1[1] * a2[2] - a2[1] * a1[2]) -
        a0[1] * (a1[0] * a2[2] - a2[0] * a1[2]) +
        a0[2] * (a1[0] * a2[1] - a2[0] * a1[1]);
}

/// Computes the volume of a tetrahedron consisting of 4 vertices
/// with 3-dimensional coordinates
inline double tetVolume(const double* p0,
                        const double* p1,
                        const double* p2,
                        const double* p3)
{
    double a[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
    double b[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
    double c[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
    return std::fabs(determinantOf(a, b, c) / 6.0);
}

/// Calculates the area of a triangle consisting of 3 vertices
/// with 2-dimensional coordinates
inline double triangleArea2d(const double* p0,
                             const double* p1,
                             const double* p2)
{
    double a[2] = { p1[0] - p0[0], p1[1] - p0[1] };
    double b[2] = { p2[0] - p0[0], p2[1] - p0[1] };
    double a_cross_b = a[0]*b[1] - a[1]*b[0];
    return 0.5*std::fabs(a_cross_b);
}

} // anonymous namespace

namespace Opm {

CellQuadrature::CellQuadrature(const UnstructuredGrid& grid,
                               const int cell,
                               const int degree)
    : grid_(grid), cell_(cell), degree_(degree)
{
    if (grid.dimensions > 3) {
        OPM_THROW(std::runtime_error, "CellQuadrature only implemented for up to 3 dimensions.");
    }
    if (degree > 2) {
        OPM_THROW(std::runtime_error, "CellQuadrature exact for polynomial degrees > 1 not implemented.");
    }
}

int CellQuadrature::numQuadPts() const
{
    if (degree_ < 2 || grid_.dimensions == 1) {
        return 1;
    }
    // Degree 2 case.
    if (grid_.dimensions == 2) {
        return 3*(grid_.cell_facepos[cell_ + 1] - grid_.cell_facepos[cell_]);
    }
    assert(grid_.dimensions == 3);
    int sumnodes = 0;
    for (grid_size_t hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
        const int face = grid_.cell_faces[hf];
        sumnodes += grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
    }
    return 4*sumnodes;
}

void CellQuadrature::quadPtCoord(const int index, double* coord) const
{
    const int dim = grid_.dimensions;
    const double* cc = grid_.cell_centroids + dim*cell_;
    if (degree_ < 2) {
        std::copy(cc, cc + dim, coord);
        return;
    }
    // Degree 2 case.
    if (dim == 2) {
        if (index % 3 == 0) {
            // Boundary midpoint. This is the face centroid.
            const int hface = grid_.cell_facepos[cell_] + index/3;
            const int face = grid_.cell_faces[hface];
            const double* fc = grid_.face_centroids + dim*face;
            std::copy(fc, fc + dim, coord);
        } else {
            // Interiour midpoint. This is the average of the
            // cell centroid and a face node (they should
            // always have two nodes in 2d).
            const int hface = grid_.cell_facepos[cell_] + index/3;
            const int face = grid_.cell_faces[hface];
            const int nodeoff = (index % 3) - 1; // == 0 or 1
            const int node = grid_.face_nodes[grid_.face_nodepos[face] + nodeoff];
            const double* nc = grid_.node_coordinates + dim*node;
            for (int dd = 0; dd < dim; ++dd) {
                coord[dd] = 0.5*(nc[dd] + cc[dd]);
            }
        }
        return;
    }
    assert(dim == 3);
    int tetindex = index / 4;
    const int subindex = index % 4;
    const double* nc = grid_.node_coordinates;
    for (grid_size_t hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
        const int face = grid_.cell_faces[hf];
        const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
        if (nfn <= tetindex) {
            // Our tet is not associated with this face.
            tetindex -= nfn;
            continue;
        }
        const double* fc = grid_.face_centroids + dim*face;
        const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
        const int node0 = fnodes[tetindex];
        const int node1 = fnodes[(tetindex + 1) % nfn];
        const double* n0c = nc + dim*node0;
        const double* n1c = nc + dim*node1;
        const double a = 0.138196601125010515179541316563436;
        // Barycentric coordinates of our point in the tet.
        double baryc[4] = { a, a, a, a };
        baryc[subindex] = 1.0 - 3.0*a;
        for (int dd = 0; dd < dim; ++dd) {
            coord[dd] = baryc[0]*cc[dd] + baryc[1]*fc[dd] + baryc[2]*n0c[dd] + baryc[3]*n1c[dd];
        }
        return;
    }
    OPM_THROW(std::runtime_error, "Should never reach this point.");
}

double CellQuadrature::quadPtWeight(const int index) const
{
    if (degree_ < 2) {
        return grid_.cell_volumes[cell_];
    }
    // Degree 2 case.
    const int dim = grid_.dimensions;
    const double* cc = grid_.cell_centroids + dim*cell_;
    if (dim == 2) {
        const int hface = grid_.cell_facepos[cell_] + index/3;
        const int face = grid_.cell_faces[hface];
        const int* nptr = grid_.face_nodes + grid_.face_nodepos[face];
        const double* nc0 = grid_.node_coordinates + dim*nptr[0];
        const double* nc1 = grid_.node_coordinates + dim*nptr[1];
        return triangleArea2d(nc0, nc1, cc)/3.0;
    }
    assert(dim == 3);
    int tetindex = index / 4;
    const double* nc = grid_.node_coordinates;
    for (grid_size_t hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
        const int face = grid_.cell_faces[hf];
        const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
        if (nfn <= tetindex) {
            // Our tet is not associated with this face.
            tetindex -= nfn;
            continue;
        }
        const double* fc = grid_.face_centroids + dim*face;
        const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
        const int node0 = fnodes[tetindex];
        const int node1 = fnodes[(tetindex + 1) % nfn];
        const double* n0c = nc + dim*node0;
        const double* n1c = nc + dim*node1;
        return 0.25*tetVolume(cc, fc, n0c, n1c);
    }
    OPM_THROW(std::runtime_error, "Should never reach this point.");
}

} // namespace Opm
