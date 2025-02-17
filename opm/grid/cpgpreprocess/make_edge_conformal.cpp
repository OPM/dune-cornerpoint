/*
  Copyright 2025 Equinor AS

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

#include <opm/grid/cpgpreprocess/make_edge_conformal.hpp>

#include <opm/grid/cpgpreprocess/preprocess.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void my_assert(const bool fine, std::string_view message = "")
{
    if (! fine) {
        std::cerr << message << std::endl;
        throw std::runtime_error { message.data() };
    }
}

struct less_than_key {
    inline bool operator()(const std::array<int,2>& edge1,
                           const std::array<int,2>& edge2) const
    {
        std::array<int,2> sedge1 = edge1;
        std::array<int,2> sedge2 = edge2;
        std::sort(sedge1.begin(), sedge1.end());
        std::sort(sedge2.begin(), sedge2.end());

        if (sedge1[0] < sedge2[0]) {
            return true;
        }
        else if (sedge1[0] == sedge2[0]) {
            return sedge1[1] < sedge2[1];
        }
        else {
            return false;
        }

        my_assert(false, "compare");
    }
};

struct opposite {
    inline bool operator()(const std::array<int,2>& edge1,
                           const std::array<int,2>& edge2) const
    {
        return (edge1[0] == edge2[1])
            && (edge1[1] == edge2[0]);
    }
};

std::vector<int>
sorted_outer_boundary(const processed_grid&                  grid,
                      const std::array<std::vector<int>, 6>& dir_faces,
                      const std::array<std::vector<int>, 6>& dir_hfaces,
                      const int                              dir,
                      const int                              cell)
{
    std::vector<std::array<int,2>> edges{};

    for (auto locind = 0*dir_faces[dir].size();
         locind < dir_faces[dir].size(); ++locind)
    {
        const int hface = dir_hfaces[dir][locind];
        const int face = grid.cell_faces[hface];

        std::vector<std::array<int,2>> face_edges{};

        // Add edges
        for (auto fnodePos = grid.face_ptr[face + 0];
             fnodePos < grid.face_ptr[face + 1] - 1;
             ++fnodePos)
        {
            face_edges.emplace_back() = std::array {
                grid.face_nodes[fnodePos + 0],
                grid.face_nodes[fnodePos + 1],
            };
        }

        face_edges.emplace_back() = std::array {
            grid.face_nodes[grid.face_ptr[face + 1] - 1],
            grid.face_nodes[grid.face_ptr[face + 0] + 0],
        };

        if (cell == grid.face_neighbors[2*face + 1]) {
            // Reverse edges
            for (auto& edge : face_edges) {
                std::swap(edge[0], edge[1]);
            }
        }

        edges.insert(edges.end(), face_edges.begin(), face_edges.end());
    }

    // Sort so intenal edges is after each other
    std::sort(edges.begin(), edges.end(), less_than_key{});

    std::vector<std::array<int,2>> new_edges;

    // Remove internal edges.
    auto iter = edges.begin();
    auto iternext = iter;
    ++iternext;

    const opposite is_opposite{};
    for (; iternext != edges.end();) {
        if (is_opposite(*iter, *iternext)) {
            iter = iternext;
            ++iter;
            if (iter != edges.end()) {
                iternext = iter;
                ++iternext;
            }
            else {
                iternext = iter;
            }
        }
        else {
            new_edges.push_back(*iter);
            ++iter;
            ++iternext;
        }
    }

    if (iter != edges.end()) {
        new_edges.push_back(*iter);
    }

    // sort so intenal edges is after each other
    edges = new_edges;

    // Now edges should only contain oriented outer edges.  Now put them
    // after each other.
    std::vector<std::array<int,2>> sedges{};

    sedges.push_back(edges.front());
    edges.erase(edges.begin());
    while (!edges.empty()) {
        auto cedges = sedges.back();
        auto nextelem = std::find_if(edges.begin(), edges.end(),
                                     [cedges](std::array<int,2> &s)
                                     { return s[0] == cedges[1]; });

        sedges.push_back(*nextelem);
        edges.erase(nextelem);
    }

    std::vector<int> sedge;
    for (const auto& edge: sedges) {
        sedge.push_back(edge[0]);
    }

    return sedge;
}

std::vector<int>
new_tb(const std::vector<int>&  bfnodes_in,
       const std::vector<int>&  sedge,
       const std::array<int,2>& bedge,
       const int                fsigntb)
{
    std::vector<int> newedge{};
    std::vector<int> bfnodes = bfnodes_in;

    const auto ind2 = find_if(sedge.begin(), sedge.end(), [&bedge](const int s) { return s == bedge[1]; });
    const auto ind1 = find_if(sedge.begin(), sedge.end(), [&bedge](const int s) { return s == bedge[0]; });

    my_assert(ind1 != sedge.end(), "ind1 not found");
    my_assert(ind2 != sedge.end(), "ind2 not found");

    int addnode = 0;
    if ((ind1 == std::prev(sedge.end())) && (ind2 == sedge.begin())) {
        addnode = 0;
    }
    else if (ind2 - ind1 < 0) {
        newedge.insert(newedge.end(), ind1, sedge.end());
        newedge.insert(newedge.end(), sedge.begin(), std::next(ind2));
        addnode = newedge.size() -2;
    }
    else if (ind2 - ind1 > 0) {
        newedge.insert(newedge.end(), ind1, std::next(ind2));
        addnode = newedge.size() -2;
    }
    else {
        my_assert(false,"Do not exist");
    }

    if (fsigntb == 1) {
        std::reverse(newedge.begin(), newedge.end());
    }

    // ett possibly modified nodes
    if (addnode > 0) {
        auto iterstart = newedge.begin();
        ++iterstart;

        auto iterend = newedge.end();
        --iterend;

        auto node2 = std::find_if(bfnodes.begin(), bfnodes.end(),
                                  [newedge](int s) { return newedge.back() == s; });

        auto node1 = std::find_if(bfnodes.begin(), bfnodes.end(),
                                  [newedge](int s) {return *(newedge.begin()) == s; });

        if ((node1 == std::prev(bfnodes.end())) &&
            (node2 == bfnodes.begin())) {
            bfnodes.insert(bfnodes.end(), iterstart, iterend);
        }
        else if (node2 - node1 == 1) {
            bfnodes.insert(std::next(node1), iterstart, iterend);
        }
        else {
            // nodes should al readdy be added
            auto node_end = node2;
            auto it1 = node1;
            auto it2 = newedge.begin();
            while ((it1 != node_end) && (it1 != bfnodes.end())) {
                my_assert(*it1 == *it2, "Existing face wrong");
                ++it1;
                ++it2;
            }

            if ((it1 == bfnodes.end()) && (node2 != std::prev(bfnodes.end()))) {
                it1 = bfnodes.begin();
                while ((it1 != node_end) && (it1 != bfnodes.end())) {
                    my_assert(*it1 == *it2, "Existing face wrong cyclic case");
                    ++it1;
                    ++it2;
                }
            }

            addnode = 0;
        }
    }

    return bfnodes;
}

void fix_edges_at_top(const struct processed_grid& grid,
                      std::vector<int>& nodes,
                      std::vector<int>& nodePos)
{
    // are going to be the new face nodes
    std::vector<std::vector<int>> face_nodes{};
    face_nodes.reserve(grid.number_of_faces);
    for (auto i = 0*grid.number_of_faces; i < grid.number_of_faces; ++i) {
        face_nodes.emplace_back(grid.face_nodes + grid.face_ptr[i + 0],
                                grid.face_nodes + grid.face_ptr[i + 1]);
    }

    const auto nhf = grid.cell_facePos[grid.number_of_cells];

    for (auto cell = 0*grid.number_of_cells; cell < grid.number_of_cells; ++cell) {
        // Process a cells and top and bottom face number.
        std::array<std::vector<int>, 6> dir_faces{};
        std::array<std::vector<int>, 6> dir_hfaces{};

        for (auto hface = grid.cell_facePos[cell + 0];
             hface < grid.cell_facePos[cell + 1]; ++hface)
        {
            const auto hface_tag = grid.cell_faces[1*nhf + hface];

            dir_faces [hface_tag].push_back(grid.cell_faces[0*nhf + hface]);
            dir_hfaces[hface_tag].push_back(hface);
        }

        my_assert(dir_faces[4].size() == 1, "face size wrong bottom");
        my_assert(dir_faces[5].size() == 1, "face size wrong top");

        for (int dir = 0; dir < 4; ++dir) {
            if (dir_faces[dir].size() <= 1) {
                continue;
            }

            // Cannot be intersection with only one face process each
            // vertical direction.  Find all oriented edges in each
            // direction
            const auto sedge = sorted_outer_boundary
                (grid, dir_faces, dir_hfaces, dir, cell);

            // Now sedge should be orient orderd list of edges find
            // top/bottom edge to be considered
            std::array<int,2> bedge{};

            // loop over top and bootm
            for (int tb = 4; tb < 6; ++tb) {
                my_assert(dir_faces[tb].size() == 1, "face size wrong tb/bottom");

                const int bface = dir_faces[tb][0];
                const auto org_bfnodes = std::vector<int> {
                    grid.face_nodes + grid.face_ptr[bface + 0],
                    grid.face_nodes + grid.face_ptr[bface + 1]
                };

                const auto bfnodes = face_nodes[bface]; //bd

                {
                    const std::array<int,4> odir {0, 2, 1, 3};
                    if (odir[dir] == 0) {
                        bedge[0] = org_bfnodes[3];
                        bedge[1] = org_bfnodes[0];
                    }
                    else {
                        bedge[0] = org_bfnodes[odir[dir] - 1];
                        bedge[1] = org_bfnodes[odir[dir] + 0];
                    }
                }

                const int fsigntb = (grid.face_neighbors[2*bface + 1] == cell)
                    ? -1 : 1;

                if (fsigntb == 1) {
                    std::reverse(bedge.begin(), bedge.end());
                }

                face_nodes[bface] = new_tb(bfnodes, sedge, bedge, fsigntb);
            } // end tb
        } // end dir
    } // end cell

    // Add the others
    for (auto face = 0*face_nodes.size(); face < face_nodes.size(); ++face) {
        auto& fnodes = face_nodes[face];

        if (! fnodes.empty()) {
            continue;
        }

        fnodes.assign(grid.face_nodes + grid.face_ptr[face + 0],
                      grid.face_nodes + grid.face_ptr[face + 1]);
    }

    nodePos.resize(face_nodes.size() + 1);
    nodePos[0] = 0;
    int count = 0;
    for (const auto& face : face_nodes) {
        nodes.insert(nodes.end(), face.begin(), face.end());
        nodePos[count + 1] = nodePos[count] + face.size();
        ++count;
    }
}

} // Anonymous namespace

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

int make_edge_conformal(struct processed_grid* grid)
{
    std::cerr << "Fixing edge grid to be edge conformal\n";

    try {
        std::vector<int> nodes{};
        std::vector<int> nodePos{};

        fix_edges_at_top(*grid, nodes, nodePos);

        std::copy_n(nodePos.begin(), grid->number_of_faces + 1, grid->face_ptr);
        std::copy_n(nodes.begin(), grid->face_ptr[grid->number_of_faces], grid->face_nodes);
    }
    catch (const std::exception& e) {
        std::cerr << "Failed to form edge-conformal grid structure\n";
        return 0;
    }

    return 1;
}

#ifdef __cplusplus
}
#endif // __cplusplus
