//===========================================================================
//
// File: GridPartitioning.cpp
//
// Created: Mon Sep  7 10:18:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2013, 2022 Equinor ASA.
  Copyright 2013, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015       NTNU
  Copyright 2020, 2021 OPM-OP AS
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

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/istl/owneroverlapcopy.hh>
#include "GridPartitioning.hpp"
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridDataTraits.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <stack>

#ifdef HAVE_MPI
#include "mpi.h"
#endif
namespace Dune
{


    typedef std::array<int, 3> coord_t;

    namespace
    {

        struct IndexToIJK
        {
            IndexToIJK(const coord_t& lc_size)
                : num_i(lc_size[0]),
                  num_ij(lc_size[0]*lc_size[1])
            {
            }
            coord_t operator()(int index)
            {
                coord_t retval = {{ index % num_i,
                                    (index % num_ij) / num_i,
                                    index / num_ij }};
                return retval;
            }
        private:
            int num_i;
            int num_ij;
        };


        int initialPartition(const coord_t& c, const coord_t& lc_size, const coord_t& initial_split)
        {
            coord_t p_coord;
            for (int i = 0; i < 3; ++i) {
                int n = lc_size[i]/initial_split[i];
                int extra = lc_size[i] % initial_split[i];
                if (c[i] < (n+1)*extra) {
                    p_coord[i] = c[i]/(n+1);
                } else {
                    p_coord[i] = (c[i] - (n+1)*extra)/n + extra;
                }
            }
            return p_coord[0] + initial_split[0]*(p_coord[1] + initial_split[1]*p_coord[2]);
        }

        template<class Entity>
        void colourMyComponentRecursive(const CpGrid& grid,
                                        const Entity& c,
                                        const int colour,
                                        const std::vector<int>& cell_part,
                                        std::vector<int>& cell_colour)
        {
            const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
            int my_index = ix.index(c);
            cell_colour[my_index] = colour;
            // For each neighbour...
            for (CpGrid::LeafIntersectionIterator it = c.ileafbegin(); it != c.ileafend(); ++it) {
                if (it->neighbor()) {
                    int nb_index = ix.index(it->outside());
                    if (cell_part[my_index] == cell_part[nb_index] && cell_colour[nb_index] == -1) {
                        colourMyComponentRecursive(grid, it->outside(), colour, cell_part, cell_colour);
                    }
                }
            }
        }

        template<class Entity>
        void colourMyComponent(const CpGrid& grid,
                               const Entity& c,
                               const int colour,
                               const std::vector<int>& cell_part,
                               std::vector<int>& cell_colour)
        {
            typedef CpGrid::LeafIntersectionIterator NbIter;
            typedef std::pair<int, std::pair<NbIter, NbIter> > VertexInfo;
            std::stack<VertexInfo> v_stack;
            const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
            int index = ix.index(*c);
            cell_colour[index] = colour;
            NbIter cur = c->ileafbegin();
            NbIter end = c->ileafend();
            v_stack.push(std::make_pair(index, std::make_pair(cur, end)));
            while (!v_stack.empty()) {
                index = v_stack.top().first;
                cur = v_stack.top().second.first;
                end = v_stack.top().second.second;
                v_stack.pop();
                while (cur != end) {
                    bool visit_nb = false;
                    if (cur->neighbor()) {
                        int nb_index = ix.index(cur->outside());
                        if (cell_part[index] == cell_part[nb_index] && cell_colour[nb_index] == -1) {
                            visit_nb = true;
                        }
                    }
                    if (visit_nb) {
                        NbIter cur_cp = cur;
                        v_stack.push(std::make_pair(index, std::make_pair(++cur, end)));
                        index = ix.index(cur_cp->outside());
                        cur = cur_cp->outside().ileafbegin();
                        end = cur_cp->outside().ileafend();
                        cell_colour[index] = colour;
                    } else {
                        ++cur;
                    }
                }
            }
        }


        void ensureConnectedPartitions(const CpGrid& grid,
                                       int& num_part,
                                       std::vector<int>& cell_part,
                                       bool recursive = false)
        {
            std::vector<int> cell_colour(cell_part.size(), -1);
            std::vector<int> partition_used(num_part, 0);
            int max_part = num_part;
            const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
            for (CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>(); it != grid.leafend<0>(); ++it) {
                int index = ix.index(*it);
                if (cell_colour[index] == -1) {
                    int part = cell_part[index];
                    int current_colour = part;
                    if (partition_used[part]) {
                        current_colour = max_part++;
                    } else {
                        partition_used[part] = true;
                    }
                    if (recursive) {
                        colourMyComponentRecursive(grid, it, current_colour, cell_part, cell_colour);
                    } else {
                        colourMyComponent(grid, it, current_colour, cell_part, cell_colour);
                    }
                }
            }
            if (max_part != num_part) {
                num_part = max_part;
                cell_part.swap(cell_colour);
            }
        }

    } // anon namespace


    void partition(const CpGrid& grid,
                   const coord_t& initial_split,
                   int& num_part,
                   std::vector<int>& cell_part,
                   bool recursive,
                   bool ensureConnectivity)
    {
        // Checking that the initial split makes sense (that there may be at least one cell
        // in each expected partition).
        const coord_t& lc_size = grid.logicalCartesianSize();
        for (int i = 0; i < 3; ++i) {
            if (initial_split[i] > lc_size[i]) {
                OPM_THROW(std::runtime_error,
                          "In direction "  + std::to_string(i) +
                          " requested splitting " + std::to_string(initial_split[i]) +
                          " size " + std::to_string(lc_size[i]));
            }
        }

        // Initial partitioning depending on (ijk) coordinates.
        std::vector<int>::size_type  num_initial =
            initial_split[0]*initial_split[1]*initial_split[2];
        const std::vector<int>& lc_ind = grid.globalCell();
        std::vector<int> num_in_part(num_initial, 0); // no cells of partitions
        std::vector<int> my_part(grid.size(0), -1); // contains partition number of cell
        IndexToIJK ijk_coord(lc_size);
        for (int i = 0; i < grid.size(0); ++i) {
            coord_t ijk = ijk_coord(lc_ind[i]);
            int part = initialPartition(ijk, lc_size, initial_split);
            my_part[i] = part;
            ++num_in_part[part];
        }

        // Renumber partitions.
        std::vector<int> num_to_subtract(num_initial); // if partitions are empty they do not get a number.
        num_to_subtract[0] = 0;
        for (std::vector<int>::size_type i = 1; i < num_initial; ++i) {
            num_to_subtract[i] = num_to_subtract[i-1];
            if (num_in_part[i-1] == 0) {
                ++num_to_subtract[i];
            }
        }
        for (int i = 0; i < grid.size(0); ++i) {
            my_part[i] -= num_to_subtract[my_part[i]];
        }

        num_part = num_initial - num_to_subtract.back();
        cell_part.swap(my_part);

        // Check the connectivity, split.
        if ( ensureConnectivity )
        {
            ensureConnectedPartitions(grid, num_part, cell_part, recursive);
        }
    }

/// \brief Adds cells to the overlap that just share a point with an owner cell.
void addOverlapCornerCell(const CpGrid& grid, int owner,
                          const CpGrid::Codim<0>::Entity& from,
                          const CpGrid::Codim<0>::Entity& neighbor,
                          const std::vector<int>& cell_part,
                          std::vector<std::set<int> >& cell_overlap)
{
    const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
    int my_index = ix.index(from);
    int nb_index = ix.index(neighbor);
    const int num_from_subs = from.subEntities(CpGrid::dimension);
    for ( int i = 0; i < num_from_subs ; i++ )
    {
        int mypoint = ix.index(from.subEntity<CpGrid::dimension>(i));
        const int num_nb_subs = neighbor.subEntities(CpGrid::dimension);
        for ( int j = 0; j < num_nb_subs; j++)
        {
            int otherpoint = ix.index(neighbor.subEntity<CpGrid::dimension>(j));
            if ( mypoint == otherpoint )
            {
                cell_overlap[nb_index].insert(owner);
                cell_overlap[my_index].insert(cell_part[nb_index]);
                return;
            }
        }
    }
}

/// \brief Adds cells to the overlap that just share a point with an owner cell.
void addOverlapCornerCell(const CpGrid& grid, int owner,
                          const CpGrid::Codim<0>::Entity& from,
                          const CpGrid::Codim<0>::Entity& neighbor,
                          const std::vector<int>& cell_part,
                          std::vector<std::tuple<int,int,char>>& exportList)
{
    // Add corner cells to the overlap layer. Example of a subdomain of a 4x4 grid
    // with and without corner cells in the overlap is given below. Note that the
    // corner cell is not needed for cell centered finite volume schemes. 
    // I = interior cells, O = overlap cells and E = exterior cells.
    //
    //  With corner     Without corner
    //  I I O E         I I O E 
    //  I I O E         I I O E 
    //  O O O E         O O E E
    //  E E E E         E E E E
    using AttributeSet = Dune::cpgrid::CpGridDataTraits::AttributeSet;
    const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
    int my_index = ix.index(from);
    int nb_index = ix.index(neighbor);
    const int num_from_subs = from.subEntities(CpGrid::dimension);
    for ( int i = 0; i < num_from_subs ; i++ )
    {
        int mypoint = ix.index(from.subEntity<CpGrid::dimension>(i));
        const int num_nb_subs = neighbor.subEntities(CpGrid::dimension);
        for ( int j = 0; j < num_nb_subs; j++)
        {
            int otherpoint = ix.index(neighbor.subEntity<CpGrid::dimension>(j));
            if ( mypoint == otherpoint )
            {
                // Note: multiple adds for same process are possible
                exportList.emplace_back(nb_index, owner, AttributeSet::copy);
                exportList.emplace_back(my_index, cell_part[nb_index],  AttributeSet::copy);
                return;
            }
        }
    }
}

void addOverlapLayer(const CpGrid& grid, int index, const CpGrid::Codim<0>::Entity& e,
                     const int owner, const std::vector<int>& cell_part,
                     std::vector<std::set<int> >& cell_overlap, int recursion_deps)
    {
        const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
        for (CpGrid::LeafIntersectionIterator iit = e.ileafbegin(); iit != e.ileafend(); ++iit) {
            if ( iit->neighbor() ) {
                int nb_index = ix.index(iit->outside());
                if ( cell_part[nb_index]!=owner )
                {
                    cell_overlap[nb_index].insert(owner);
                    cell_overlap[index].insert(cell_part[nb_index]);
                    if ( recursion_deps>0 )
                    {
                        // Add another layer
                        addOverlapLayer(grid, nb_index, iit->outside(), owner,
                                        cell_part, cell_overlap, recursion_deps-1);
                    }
                    else
                    {
                        // Add cells to the overlap that just share a corner with e.
                        for (CpGrid::LeafIntersectionIterator iit2 = iit->outside().ileafbegin();
                             iit2 != iit->outside().ileafend(); ++iit2)
                        {
                           if ( iit2->neighbor() )
                           {
                               int nb_index2 = ix.index(iit2->outside( ));
                               if( cell_part[nb_index2]==owner ) continue;
                               addOverlapCornerCell(grid, owner, e, iit2->outside(),
                                                    cell_part, cell_overlap);
                           }
                       }
                    }
                }
            }
        }
    }

    void addOverlapLayer(const CpGrid& grid, const std::vector<int>& cell_part,
                         std::vector<std::set<int> >& cell_overlap, int mypart,
                         int layers, bool all)
    {
        cell_overlap.resize(cell_part.size());
        const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
        for (CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>();
             it != grid.leafend<0>(); ++it) {
            int index = ix.index(*it);
            int owner = -1;
            if(cell_part[index]==mypart)
                owner = mypart;
            else
            {
                if(all)
                    owner=cell_part[index];
                else
                    continue;
            }
            addOverlapLayer(grid, index, *it, owner, cell_part, cell_overlap, layers-1);
        }
    }

    void addOverlapLayer(const CpGrid& grid, int index, const CpGrid::Codim<0>::Entity& e,
                         const int owner, const std::vector<int>& cell_part,
                         std::vector<std::tuple<int,int,char>>& exportList,
                         bool addCornerCells, int recursion_deps)
    {
        using AttributeSet = Dune::cpgrid::CpGridDataTraits::AttributeSet;
        const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
        for (CpGrid::LeafIntersectionIterator iit = e.ileafbegin(); iit != e.ileafend(); ++iit) {
            if ( iit->neighbor() ) {
                int nb_index = ix.index(iit->outside());
                if ( cell_part[nb_index]!=owner )
                {
                    // Note: multiple adds for same process are possible
                    exportList.emplace_back(nb_index, owner, AttributeSet::copy);
                    exportList.emplace_back(index, cell_part[nb_index],  AttributeSet::copy);
                    if ( recursion_deps>0 )
                    {
                        // Add another layer
                        addOverlapLayer(grid, nb_index, iit->outside(), owner,
                                        cell_part, exportList, addCornerCells, recursion_deps-1);
                    }
                    else if (addCornerCells)
                    {
                        // Add cells to the overlap that just share a corner with e.
                        for (CpGrid::LeafIntersectionIterator iit2 = iit->outside().ileafbegin();
                             iit2 != iit->outside().ileafend(); ++iit2)
                        {
                            if ( iit2->neighbor() )
                            {
                                int nb_index2 = ix.index(iit2->outside());
                                if( cell_part[nb_index2]!=owner ) {
                                    addOverlapCornerCell(grid, owner, e, iit2->outside(),
                                                         cell_part, exportList);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void addOverlapLayerNoZeroTrans(const CpGrid& grid, int index, const CpGrid::Codim<0>::Entity& e,
                                    const int owner, const std::vector<int>& cell_part,
                                    std::vector<std::tuple<int,int,char>>& exportList,
                                    bool addCornerCells, int recursion_deps, const double* trans)
    {
        using AttributeSet = Dune::OwnerOverlapCopyAttributeSet::AttributeSet;
        const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
        for (CpGrid::LeafIntersectionIterator iit = e.ileafbegin(); iit != e.ileafend(); ++iit) {
            if ( iit->neighbor() ) {
                int faceId = iit->id();
                int nb_index = ix.index(iit->outside());
                if ( cell_part[nb_index]!=owner )
                {
                    // Note: multiple adds for same process are possible
                    exportList.emplace_back(nb_index, owner, AttributeSet::copy);
                    exportList.emplace_back(index, cell_part[nb_index],  AttributeSet::copy);

                    if ( recursion_deps>0 )
                    {
                        // Add another layer
                        addOverlapLayerNoZeroTrans(grid, nb_index, e, owner, cell_part,
                                                   exportList, addCornerCells, recursion_deps-1, trans);
                    }
                    else if (addCornerCells)
                    {
                        // Add cells to the overlap that just share a corner with e.
                        for (CpGrid::LeafIntersectionIterator iit2 = iit->outside().ileafbegin();
                             iit2 != iit->outside().ileafend(); ++iit2)
                        {
                            if ( iit2->neighbor() )
                            {
                                int nb_index2 = ix.index(iit2->outside());
                                if( cell_part[nb_index2]!=owner ) {
                                    addOverlapCornerCell(grid, owner, e, iit2->outside(),
                                                         cell_part, exportList);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int addOverlapLayer(const CpGrid& grid, const std::vector<int>& cell_part,
                        std::vector<std::tuple<int,int,char>>& exportList,
                        std::vector<std::tuple<int,int,char,int>>& importList,
                        const Communication<Dune::MPIHelper::MPICommunicator>& cc,
                        [[maybe_unused]] bool addCornerCells,
                        [[maybe_unused]] const double* trans, int layers)
    {
#ifdef HAVE_MPI
        using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
        auto ownerSize = exportList.size();
        const CpGrid::LeafIndexSet& ix = grid.leafIndexSet();
        std::map<int,int> exportProcs, importProcs;

        for (CpGrid::Codim<0>::LeafIterator it = grid.leafbegin<0>();
             it != grid.leafend<0>(); ++it) {
            int index = ix.index(*it);
            auto owner = cell_part[index];
            exportProcs.insert(std::make_pair(owner, 0));
            if ( trans ) {
                addOverlapLayerNoZeroTrans(grid, index, *it, owner, cell_part, exportList, addCornerCells, layers-1, trans);
            }
            else {
                addOverlapLayer(grid, index, *it, owner, cell_part, exportList, addCornerCells, layers-1);
            }
        }
        // remove multiple entries
        auto compare = [](const std::tuple<int,int,char>& t1, const std::tuple<int,int,char>& t2)
                       {
                           return (std::get<0>(t1) < std::get<0>(t2)) ||
                                                    ( ! (std::get<0>(t2) < std::get<0>(t1))
                                                      && (std::get<1>(t1) < std::get<1>(t2)) );
                       };

        auto ownerEnd = exportList.begin() + ownerSize;
        std::sort(ownerEnd, exportList.end(), compare);

        auto newEnd = std::unique(ownerEnd, exportList.end());
        exportList.resize(newEnd - exportList.begin());

        for(const auto& entry: importList)
            importProcs.insert(std::make_pair(std::get<1>(entry), 0));
        //count entries to send
        std::for_each(ownerEnd, exportList.end(),
                      [&exportProcs](const std::tuple<int,int,char>& t){ ++exportProcs[std::get<1>(t)];});

        // communicate number of entries
        std::vector<MPI_Request> requests(importProcs.size());
        auto req = requests.begin();
        int tag = 2387;
        for(auto&& proc : importProcs)
        {
            MPI_Irecv(&(proc.second), 1, MPI_INT, proc.first, tag, cc, &(*req));
            ++req;
        }

        for(const auto& proc: exportProcs)
        {
            MPI_Send(&(proc.second), 1, MPI_INT, proc.first, tag, cc);
        }
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), requests.data(), statuses.data());

        // Communicate overlap entries
        ++tag;
        std::vector<std::vector<int> > receiveBuffers(importProcs.size());
        auto buffer = receiveBuffers.begin();
        req = requests.begin();

        for(auto&& proc: importProcs)
        {
            buffer->resize(proc.second);
            MPI_Irecv(buffer->data(), proc.second, MPI_INT, proc.first, tag, cc, &(*req));
            ++req; ++buffer;
        }

        for(const auto& proc: exportProcs)
        {
            std::vector<int> sendBuffer;
            sendBuffer.reserve(proc.second);
            std::for_each(ownerEnd, exportList.end(),
                          [&sendBuffer, &proc](const std::tuple<int,int,char>& t)
                          {
                              if ( std::get<1>(t) == proc.first )
                                  sendBuffer.push_back(std::get<0>(t));
                          });
            MPI_Send(sendBuffer.data(), proc.second, MPI_INT, proc.first, tag, cc);
        }

        std::inplace_merge(exportList.begin(), ownerEnd, exportList.end());

        MPI_Waitall(requests.size(), requests.data(), statuses.data());
        buffer = receiveBuffers.begin();
        auto importOwnerSize = importList.size();

        for(const auto& proc: importProcs)
        {
            for(const auto& index: *buffer)
                importList.emplace_back(index, proc.first, AttributeSet::copy, -1);
            ++buffer;
        }
        std::sort(importList.begin() + importOwnerSize, importList.end(),
                  [](const std::tuple<int,int,char,int>& t1, const std::tuple<int,int,char,int>& t2)
                  { return std::get<0>(t1) < std::get<0>(t2);});
        return importOwnerSize;
#else
        (void) grid;
        (void) cell_part;
        (void) exportList;
        (void) importList;
        (void) cc;
        (void) layers;
        DUNE_THROW(InvalidStateException, "MPI is missing from the system");

        return 0;
#endif
    }

namespace cpgrid
{
#if HAVE_MPI

   std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
              std::vector<std::tuple<int,int,char> >,
              std::vector<std::tuple<int,int,char,int> >,
              WellConnections>
    createListsFromParts(const CpGrid& grid, const std::vector<cpgrid::OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                               const double* transmissibilities, const std::vector<int>& parts,
                               bool allowDistributedWells,
                               std::shared_ptr<cpgrid::CombinedGridWellGraph> gridAndWells)
    {
        std::vector<int> exportGlobalIds;
        std::vector<int> exportLocalIds;
        std::vector<int> exportToPart;
        std::vector<int> importGlobalIds;
        std::size_t numExport = 0;
        int root = 0;

        if (grid.comm().rank() == 0)
        {
            // Create export lists as from Zoltan output, do not include part 0!
            auto numCells = grid.size(0);
            exportGlobalIds.reserve(numCells);
            exportLocalIds.reserve(numCells);
            exportToPart.reserve(numCells);
            for (auto cell = grid.leafbegin<0>(), cellEnd = grid.leafend<0>();
                 cell != cellEnd; ++cell)
            {
                const auto& gid = grid.globalIdSet().id(*cell);
                const auto& lid = grid.localIdSet().id(*cell);
                const auto& index = grid.leafIndexSet().index(cell);
                const auto& part = parts[index];
                if (part != 0 )
                {
                    exportGlobalIds.push_back(gid);
                    exportLocalIds.push_back(lid);
                    exportToPart.push_back(part);
                    ++numExport;
                }
            }
        }

        int numImport = 0;
        std::tie(numImport, importGlobalIds) =
            scatterExportInformation(numExport, exportGlobalIds.data(),
                                     exportToPart.data(), 0,
                                     grid.comm());
        if (wells and !gridAndWells)
        {
            bool partitionIsEmpty = (grid.size(0) == 0);
            EdgeWeightMethod method{}; // We don't care which method is used, we only need the graph.
            gridAndWells.reset(new cpgrid::CombinedGridWellGraph(grid,
                                                                 wells,
                                                                 possibleFutureConnections,
                                                                 transmissibilities,
                                                                 partitionIsEmpty,
                                                                 method));
        }
        return makeImportAndExportLists(grid, grid.comm(),
                                        wells,
                                        possibleFutureConnections,
                                        gridAndWells.get(),
                                        root,
                                        numExport,
                                        numImport,
                                        exportLocalIds.data(),
                                        exportGlobalIds.data(),
                                        exportToPart.data(),
                                        importGlobalIds.data(),
                                        allowDistributedWells);
    }

    std::tuple<std::vector<int>, std::vector<std::pair<std::string,bool>>,
               std::vector<std::tuple<int,int,char> >,
               std::vector<std::tuple<int,int,char,int> >,
               WellConnections>
    vanillaPartitionGridOnRoot(const CpGrid& grid, const std::vector<cpgrid::OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                               const double* transmissibilities, bool allowDistributedWells = false)
    {
        int root = 0;
        const auto& cc = grid.comm();
        std::vector<int> parts;

        if (cc.rank() == root)
        {
            std::cout<<"WARNING: Using poor man's load balancer"<<std::endl;
            parts.resize(grid.size(0));
            int  numParts=-1;
            std::array<int, 3> initialSplit;
            initialSplit[1]=initialSplit[2]=std::pow(cc.size(), 1.0/3.0);
            initialSplit[0]=cc.size()/(initialSplit[1]*initialSplit[2]);
            partition(grid, initialSplit, numParts, parts, false, false);
        }
        return createListsFromParts(grid, wells, possibleFutureConnections, transmissibilities, parts, allowDistributedWells);
    }
#endif
} // namespace cpgrid
} // namespace Dune

