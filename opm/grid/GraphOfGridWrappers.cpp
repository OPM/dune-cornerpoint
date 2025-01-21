// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include "GraphOfGridWrappers.hpp"
#include <opm/grid/common/CommunicationUtils.hpp>
#include <opm/grid/common/ZoltanPartition.hpp> // function scatterExportInformation

namespace Opm {

#if HAVE_MPI
int getGraphOfGridNumVertices(void* pGraph, int *err)
{
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int size = gog.size();
    *err = ZOLTAN_OK;
    return size;
}

void getGraphOfGridVerticesList(void* pGraph,
               [[maybe_unused]] int dimGlobalID,
               [[maybe_unused]] int dimLocalID,
                                ZOLTAN_ID_PTR gIDs,
               [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                                int weightDim,
                                float *objWeights,
                                int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    assert(weightDim==1); // vertex weight is a single float
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int i=0;
    for (const auto& v : gog)
    {
        gIDs[i] = v.first;
        // lIDs are left unused
        objWeights[i] = v.second.weight;
        ++i;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridNumEdges(void *pGraph,
           [[maybe_unused]] int dimGlobalID,
           [[maybe_unused]] int dimLocalID,
                            int numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            int *numEdges,
                            int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    for (int i=0; i<numCells; ++i)
    {
        int nE = gog.numEdges(gIDs[i]);
        if (nE== -1)
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridNumEdges error: Vertex with ID " << gIDs[i] << " is not in graph.";
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        numEdges[i] = nE;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridEdgeList(void *pGraph,
           [[maybe_unused]] int dimGlobalID,
           [[maybe_unused]] int dimLocalID,
                            int numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            int *numEdges,
                            ZOLTAN_ID_PTR nborGIDs,
                            int *nborProc,
                            int weightDim,
                            float *edgeWeights,
                            int *err)
{
    assert(dimGlobalID==1); // ID is a single int
    assert(weightDim==1); // edge weight is a single float
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    int id=0;
    for (int i=0; i<numCells; ++i)
    {
        const auto& eList = gog.edgeList(gIDs[i]);
        if ((int)eList.size()!=numEdges[i])
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridEdgeList error: Edge number disagreement"
                 << " between Zoltan (" << numEdges[i] << ") and Graph ("
                 << eList.size() << ") for vertex with ID " << gIDs[i] << std::endl;
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        for (const auto& e : eList)
        {
            nborGIDs[id]= e.first;
            nborProc[id]= gog.getVertex(e.first).nproc;
            edgeWeights[id]= e.second;
            ++id;
        }
    }
    *err = ZOLTAN_OK;
}

template<typename Zoltan_Struct>
void setGraphOfGridZoltanGraphFunctions(Zoltan_Struct *zz,
                                        GraphOfGrid<Dune::CpGrid>& gog,
                                        bool pretendNull)
{
    GraphOfGrid<Dune::CpGrid>* pGraph = &gog;
    if (pretendNull)
    {
        Zoltan_Set_Num_Obj_Fn(zz, Dune::cpgrid::getNullNumCells, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, Dune::cpgrid::getNullVertexList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, Dune::cpgrid::getNullNumEdgesList, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, Dune::cpgrid::getNullEdgeList, pGraph);
    }
    else
    {
        Zoltan_Set_Num_Obj_Fn(zz, getGraphOfGridNumVertices, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, getGraphOfGridVerticesList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getGraphOfGridNumEdges, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getGraphOfGridEdgeList, pGraph);
    }
}
#endif // HAVE_MPI

void addFutureConnectionWells(GraphOfGrid<Dune::CpGrid>& gog,
                              const std::unordered_map<std::string, std::set<int>>& wells,
                              bool checkWellIntersections)
{
    // create compressed lookup from cartesian.
    const auto& grid = gog.getGrid();
    const auto& cpgdim = grid.logicalCartesianSize();
    std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);
    for( int i=0; i < grid.numCells(); ++i )
    {
        cartesian_to_compressed[grid.globalCell()[i]] = i;
    }

    for (const auto& w: wells)
    {
        std::set<int> wellsgID;
        for (const int& cell : w.second)
        {
            int gID = cartesian_to_compressed[cell];
            assert(gID!=-1); // well should be an active cell
            wellsgID.insert(gID);
        }
        gog.addWell(wellsgID, checkWellIntersections);
    }
}

void addWellConnections(GraphOfGrid<Dune::CpGrid>& gog,
                        const Dune::cpgrid::WellConnections& wells,
                        bool checkWellIntersections)
{
    for (const auto& w : wells)
    {
        gog.addWell(w, checkWellIntersections);
    }
}

void extendGIDtoRank(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<int>& gIDtoRank,
                     const int& root)
{
    for (const auto& w : gog.getWells())
    {
        auto wellRank = gIDtoRank[*w.begin()];
        if (wellRank!=root)
        {
            for (const auto& gID : w)
            {
                gIDtoRank[gID] = wellRank;
            }
        }
    }
}

#if HAVE_MPI
namespace Impl{

std::vector<std::vector<int>>
extendRootExportList(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<std::tuple<int,int,char>>& exportList,
                     int root,
                     const std::vector<int>& gIDtoRank)
{
    const auto& cc = gog.getGrid().comm();
    // non-root ranks have empty export lists.
    std::vector<std::vector<int>> exportedCells;
    if (cc.rank()!=root)
    {
        return exportedCells;
    }

    using iter = std::set<int>::const_iterator;
    // store which wells are exported. Contains {begin, end, destination rank}
    std::vector<std::tuple<iter,iter,int>> wellsToExport;
    // track how many cells per rank are exported - to reserve the vector size
    std::vector<int> sizesOfExport(cc.size(), 0);
    const auto& gogWells = gog.getWells();
    wellsToExport.reserve(gogWells.size());

    // with gIDtoRank we can directly tell which wells belong where
    if (gIDtoRank.size()>0)
    {
        for (const auto& well : gogWells)
        {
            auto wellID = *well.begin();
            if (gIDtoRank[wellID]!=root)
            {
                wellsToExport.emplace_back(well.begin(), well.end(), gIDtoRank[wellID]);
                sizesOfExport[gIDtoRank[wellID]] += well.size()-1;  // one of well's cells is already in the list
            }
        }
    }
    else
    {
        // make a list of wells for easy identification during search. Contains {begin, end, well ID}
        std::map<int, std::tuple<iter,iter,int>> wellMap;
        for (const auto& well : gogWells)
        {
            wellMap[*well.begin()] = std::make_tuple(well.begin(), well.end(), well.size());
        }
        // iterate once through the original exportList and identify exported wells
        for (const auto& cellProperties : exportList)
        {
            auto pWell = wellMap.find(std::get<0>(cellProperties));
            if (pWell!=wellMap.end())
            {
                int rankToExport = std::get<1>(cellProperties);
                if (rankToExport!=root) // wells on root are not exported
                {
                    const auto& [begin, end, wSize] = pWell->second;
                    wellsToExport.emplace_back(begin, end, rankToExport);
                    sizesOfExport[rankToExport] += wSize-1; // one cell is already in the list
                }
                wellMap.erase(pWell);
                if (wellMap.empty())
                {
                    break;
                }
            }
        }
    }

    // create the output: the cells that are missing from importList on non-root ranks
    // also add new cells to the exportList and sort it
    exportedCells.resize(cc.size());
    int addedExportsSize = 0;
    for (int i=0; i<cc.size(); ++i)
    {
        exportedCells[i].reserve(sizesOfExport[i]);
        addedExportsSize += sizesOfExport[i];
    }
    exportList.reserve(exportList.size()+addedExportsSize);
    for (auto& pWell : wellsToExport)
    {
        auto& [begin, end, rank] = pWell;
        // remember to skip the well's first cell which already is in the importList
        for (auto& pgID = ++begin; pgID!=end; ++pgID)
        {
            using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
            exportList.emplace_back(*pgID, rank, AttributeSet::owner);

            exportedCells[rank].push_back(*pgID);
        }
    }
    std::sort(exportList.begin(), exportList.end());
    return exportedCells;
}

std::vector<int> communicateExportedCells(const std::vector<std::vector<int>>& exportedCells,
                                          const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                          int root)
{
    // send data from root
    std::vector<int> result;
    if (cc.rank() == root) {
        assert((int)exportedCells.size() == cc.size());
        assert((int)exportedCells[root].size() == 0);
        // check send success for each well
        std::vector<MPI_Request> requestSize(2 * (cc.size() - 1));
        std::vector<int> numWells(cc.size() - 1); // store well sizes (can not send temporary)

        for (int i = 0; i < cc.size() - 1; ++i) {
            int ii = i + (int)(i >= root); // ii takes values {0,...,mpisize-1} but skips root
            int tag = 37; // a random number
            numWells[i] = exportedCells[ii].size();
            MPI_Isend(&numWells[i], 1, MPI_INT, ii, tag, cc, &requestSize[2 * i]);
            MPI_Isend(exportedCells[ii].data(), numWells[i], MPI_INT, ii, tag + 1, cc, &requestSize[2 * i + 1]);
        }
        // wait till all messages are completed
        MPI_Waitall(requestSize.size(), requestSize.data(), MPI_STATUS_IGNORE);
    } else {
        // receive data from root
        int tag = 37; // a random number
        int numWells;
        MPI_Recv(&numWells, 1, MPI_INT, root, tag, cc, MPI_STATUS_IGNORE);
        result.resize(numWells);
        MPI_Recv(result.data(), numWells, MPI_INT, root, tag + 1, cc, MPI_STATUS_IGNORE);
    }
    return result;
}

void extendAndSortImportList(std::vector<std::tuple<int,int,char,int>>& importList,
                             const std::vector<int>& extraCells)
{
    if (importList.empty()) {
        // An empty importList implies that extraCells is empty
        assert(extraCells.empty());
        return;
    }
    int rank = std::get<1>(importList[0]);

    importList.reserve(importList.size() + extraCells.size());
    for (const auto& cell : extraCells) {
        // Beware: AttributeSet and localID are hard-coded
        using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
        importList.emplace_back(cell, rank, AttributeSet::owner, -1);
    }

    std::sort(importList.begin(), importList.end());
}

} // end namespace Impl

void extendAndSortExportAndImportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                                       const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                       int root,
                                       std::vector<std::tuple<int,int,char>>& exportList,
                                       std::vector<std::tuple<int,int,char,int>>& importList,
                                       const std::vector<int>& gIDtoRank)
{
    // extend root's export list and get sets of well cells for other ranks
    auto expListToComm = Impl::extendRootExportList(gog, exportList, root, gIDtoRank);
    // obtain wells on this rank from root
    auto extraWells = Impl::communicateExportedCells(expListToComm, cc, root);
    if (cc.rank()!=root)
    {
        Impl::extendAndSortImportList(importList, extraWells);
    }
}
#endif // HAVE_MPI

std::vector<int> getWellRanks(const std::vector<int>& gIDtoRank,
                              const Dune::cpgrid::WellConnections& wellConnections)
{
    std::vector<int> wellIndices(wellConnections.size());
    for (std::size_t wellIndex = 0; wellIndex < wellConnections.size(); ++wellIndex)
    {
        int wellID = *(wellConnections[wellIndex].begin());
        wellIndices[wellIndex] = gIDtoRank[wellID];
    }
    return wellIndices;
}

#if HAVE_MPI
std::vector<std::pair<std::string, bool>>
wellsOnThisRank(const std::vector<Dune::cpgrid::OpmWellType>& wells,
                const std::vector<int>& wellRanks,
                const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                int root)
{
    auto numProcs = cc.size();
    std::vector<std::vector<int>> wells_on_proc(numProcs);
    for (std::size_t i=0; i<wellRanks.size(); ++i)
    {
        wells_on_proc[wellRanks[i]].push_back(i);
    }
    return Dune::cpgrid::computeParallelWells(wells_on_proc, wells, cc, root);
}

template<class Id>
std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> > >
makeImportAndExportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                         const Dune::Communication<MPI_Comm>& cc,
                         const std::vector<Dune::cpgrid::OpmWellType> * wells,
                         const Dune::cpgrid::WellConnections& wellConnections,
                         int root,
                         int numExport,
                         int numImport,
        [[maybe_unused]] const Id* exportLocalGids,
                         const Id* exportGlobalGids,
                         const int* exportToPart,
                         const Id* importGlobalGids)
{
    const auto& cpgrid = gog.getGrid();
    int size = cpgrid.numCells();
    int rank  = cc.rank();
    std::vector<int> gIDtoRank(size, rank);
    std::vector<std::vector<int> > wellsOnProc;

    // List entry: process to export to, (global) index, process rank, attribute there (not needed?)
    std::vector<std::tuple<int,int,char>> myExportList;
    // List entry: process to import from, global index, process rank, attribute here, local index (determined later)
    std::vector<std::tuple<int,int,char,int>> myImportList;
    float buffer = 1.05; // to allocate extra space for wells in myExportList and myImportList
    assert(rank==root || numExport==0);
    assert(rank!=root || numImport==0);
    // all cells on root are added to its export and its import list
    std::size_t reserveEx = rank!=root ? 0 : cpgrid.size(0);
    std::size_t reserveIm = rank!=root ? buffer*numImport : cpgrid.size(0)*buffer/cc.size();
    myExportList.reserve(reserveEx);
    myImportList.reserve(reserveIm);
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

    for ( int i=0; i < numImport; ++i )
    {
        myImportList.emplace_back(importGlobalGids[i], root, static_cast<char>(AttributeSet::owner), -1);
    }

    if (rank==root)
    {
        for ( int i=0; i < numExport; ++i )
        {
            gIDtoRank[exportGlobalGids[i]] = exportToPart[i];
            myExportList.emplace_back(exportGlobalGids[i], exportToPart[i], static_cast<char>(AttributeSet::owner));
        }
        // partitioner sees only one cell per well, modify remaining
        extendGIDtoRank(gog, gIDtoRank, rank);

        // Add cells that stay here to the lists. Somehow I could not persuade Zoltan to do this.
        // This also adds all well cells that were missing in the importGlobalIDs.
        for ( std::size_t i = 0; i < gIDtoRank.size(); ++i)
        {
            if ( gIDtoRank[i] == rank )
            {
                myExportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner) );
                myImportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner), -1 );
            }
        }
    }


    std::vector<std::pair<std::string, bool>> parallel_wells;
    if( wells )
    {
        // complete root's export and other's import list by adding remaining well cells
        extendAndSortExportAndImportLists(gog, cc, root, myExportList, myImportList, gIDtoRank);

        auto wellRanks = getWellRanks(gIDtoRank, wellConnections);
        parallel_wells = wellsOnThisRank(*wells, wellRanks, cc, root);
    }
    return std::make_tuple( std::move(gIDtoRank),
                            std::move(parallel_wells),
                            std::move(myExportList),
                            std::move(myImportList) );
}

namespace {
void setDefaultZoltanParameters(Zoltan_Struct* zz)
{
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
#ifndef NDEBUG
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
#else
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "0");
#endif
}

} // anon namespace

std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           Dune::cpgrid::WellConnections>
zoltanPartitioningWithGraphOfGrid(const Dune::CpGrid& grid,
                                  const std::vector<Dune::cpgrid::OpmWellType> * wells,
                                  const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                  const double* transmissibilities,
                                  const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                  Dune::EdgeWeightMethod edgeWeightMethod,
                                  int root,
                                  const double zoltanImbalanceTol,
                                  const std::map<std::string, std::string>& params)
{
    int rc = ZOLTAN_OK - 1;
    float ver = 0;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int argc=0;
    char** argv = 0 ;
    rc = Zoltan_Initialize(argc, argv, &ver);
    zz = Zoltan_Create(cc);
    if ( rc != ZOLTAN_OK )
    {
        OPM_THROW(std::runtime_error, "Could not initialize Zoltan!");
    }
    setDefaultZoltanParameters(zz);
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(zoltanImbalanceTol).c_str());
    int layers = 0; // extra layers of cells attached to wells to distance them from boundary
    for (const auto& [key, value] : params)
    {
        if (key=="EnvelopeWellLayers")
            layers = std::stoi(value);
        else
            Zoltan_Set_Param(zz, key.c_str(), value.c_str());
    }

    // root process has the whole grid, other ranks nothing
    bool partitionIsEmpty = cc.rank()!=root;

    // prepare graph and contract well cells
    // non-root processes have empty grid and no wells
    GraphOfGrid gog(grid, transmissibilities, edgeWeightMethod);
    assert(gog.size()==0 || !partitionIsEmpty);
    auto wellConnections=partitionIsEmpty ? Dune::cpgrid::WellConnections()
                                          : Dune::cpgrid::WellConnections(*wells, possibleFutureConnections, grid);
    addWellConnections(gog, wellConnections);
    gog.addNeighboringCellsToWells(layers);

    // call partitioner
    setGraphOfGridZoltanGraphFunctions(zz, gog, partitionIsEmpty);
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */

    // arrange output into tuples and add well cells
    auto importExportLists = makeImportAndExportLists(gog,
                                                      cc,
                                                      wells,
                                                      wellConnections,
                                                      root,
                                                      numExport,
                                                      numImport,
                                                      exportLocalGids,
                                                      exportGlobalGids,
                                                      exportProcs,
                                                      importGlobalGids);

    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    // add wellConnections to the importExportLists and return it
    auto result = std::tuple(std::move(std::get<0>(importExportLists)),
                             std::move(std::get<1>(importExportLists)),
                             std::move(std::get<2>(importExportLists)),
                             std::move(std::get<3>(importExportLists)),
                             std::move(wellConnections));
    return result;
}

std::vector<std::vector<int> >
makeExportListsFromGIDtoRank(const std::vector<int>& gIDtoRank, int ccsize)
{
    std::vector<std::vector<int>> exportedCells(ccsize);
    // reserve vector size, we expect a balanced spread
    for (int i = 0; i < ccsize; ++i) {
        exportedCells[i].reserve(floor(gIDtoRank.size() * 1.1 / ccsize));
    }

    for (int i = 0; i < (int)gIDtoRank.size(); ++i) {
        exportedCells[gIDtoRank[i]].push_back(i);
    }
    return exportedCells;
}

namespace {
std::tuple<int, std::vector<int>>
applySerialZoltan (const Dune::CpGrid& grid,
                   const Dune::cpgrid::WellConnections& wellConnections,
                   const double* transmissibilities,
                   int numParts,
                   Dune::EdgeWeightMethod edgeWeightMethod,
                   int root,
                   const double zoltanImbalanceTol,
                   const std::map<std::string, std::string>& params)
{
    int rc = ZOLTAN_OK;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int numExport = 0, numImport = 0;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    struct Zoltan_Struct* zz;
    int changes, numGidEntries, numLidEntries;

    int argc = 0;
    char** argv = 0;
    float ver = 0;
    std::vector<int> gIDtoRank;

    rc = Zoltan_Initialize(argc, argv, &ver);
    if (rc != ZOLTAN_OK)
        return std::make_tuple(ZOLTAN_OK + 1, gIDtoRank);
    zz = Zoltan_Create(MPI_COMM_SELF);
    if (!zz)
        return std::make_tuple(ZOLTAN_OK + 2, gIDtoRank);
    setDefaultZoltanParameters(zz);
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(zoltanImbalanceTol).c_str());
    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(numParts).c_str());
    int layers = 0; // extra layers of cells attached to wells to distance them from boundary
    for (const auto& [key, value] : params) {
        if (key == "EnvelopeWellLayers")
            layers = std::stoi(value);
        else
            Zoltan_Set_Param(zz, key.c_str(), value.c_str());
    }

    // prepare graph and contract well cells
    GraphOfGrid gog(grid, transmissibilities, edgeWeightMethod);
    addWellConnections(gog, wellConnections);
    gog.addNeighboringCellsToWells(layers);

    // call partitioner
    setGraphOfGridZoltanGraphFunctions(zz, gog, false);
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes, /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries, /* Number of integers used for a global ID */
                             &numLidEntries, /* Number of integers used for a local ID */
                             &numImport, /* Number of vertices to be sent to me */
                             &importGlobalGids, /* Global IDs of vertices to be sent to me */
                             &importLocalGids, /* Local IDs of vertices to be sent to me */
                             &importProcs, /* Process rank for source of each incoming vertex */
                             &importToPart, /* New partition for each incoming vertex */
                             &numExport, /* Number of vertices I must send to other processes*/
                             &exportGlobalGids, /* Global IDs of the vertices I must send */
                             &exportLocalGids, /* Local IDs of the vertices I must send */
                             &exportProcs, /* Process to which I send each of the vertices */
                             &exportToPart); /* Partition to which each vertex will belong */
    numImport = 0;
    if (rc == ZOLTAN_OK) {
        gIDtoRank.resize(grid.numCells(), root);
        for (int i = 0; i < numExport; ++i) {
            gIDtoRank[exportGlobalGids[i]] = exportToPart[i];
        }
        // partitioner sees only one cell per well, modify remaining
        extendGIDtoRank(gog, gIDtoRank, root);
    } else {
        rc = ZOLTAN_OK + 3; // distinguish Zoltan failures
    }
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);
    return std::make_tuple(rc, gIDtoRank);
}
} // end anonymous namespace

std::tuple<std::vector<int>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<int,int,char> >,
           std::vector<std::tuple<int,int,char,int> >,
           Dune::cpgrid::WellConnections>
zoltanSerialPartitioningWithGraphOfGrid(const Dune::CpGrid& grid,
                                        const std::vector<Dune::cpgrid::OpmWellType> * wells,
                                        const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                        const double* transmissibilities,
                                        const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                        Dune::EdgeWeightMethod edgeWeightMethod,
                                        int root,
                                        const double zoltanImbalanceTol,
                                        const std::map<std::string, std::string>& params)
{
    // root process has the whole grid, other ranks nothing
    bool partitionIsEmpty = cc.rank() != root;
    int rc = ZOLTAN_OK;
    std::vector<int> gIDtoRank;
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;
    std::vector<std::tuple<int, int, char>> myExportList;
    std::vector<std::tuple<int, int, char, int>> myImportList;
    std::vector<std::vector<int>> exportedCells;
    auto wellConnections = partitionIsEmpty ? Dune::cpgrid::WellConnections()
                                            : Dune::cpgrid::WellConnections(*wells, possibleFutureConnections, grid);

    if (cc.rank() == root) {
        std::tie(rc, gIDtoRank) = applySerialZoltan(grid,
                                                    wellConnections,
                                                    transmissibilities,
                                                    cc.size(),
                                                    edgeWeightMethod,
                                                    root,
                                                    zoltanImbalanceTol,
                                                    params);
    }

    cc.broadcast(&rc, 1, root);
    if (rc != ZOLTAN_OK) {
        switch (rc) {
        case ZOLTAN_OK+1:
            OPM_THROW(std::runtime_error, "Could not initialize Zoltan!");
        case ZOLTAN_OK+2:
            OPM_THROW(std::runtime_error, "Could not create Zoltan!");
        case ZOLTAN_OK+3:
            OPM_THROW(std::runtime_error, "Partitioning with Zoltan failed!");
        default:
            OPM_THROW(std::runtime_error, "Unknown error reported by Zoltan!");
        }
    }

    if (cc.rank() == root) {
        // prepare exportedCells for communication
        exportedCells = makeExportListsFromGIDtoRank(gIDtoRank, cc.size());
        myImportList.reserve(exportedCells[root].size());
        for (const auto& cell : exportedCells[root]) {
            myImportList.emplace_back(cell, root, static_cast<char>(AttributeSet::owner), -1);
        }
        // exclude root's own cells from communication
        exportedCells[root].resize(0);
    }
    // communicate and create import+export lists
    auto importedCells = Opm::Impl::communicateExportedCells(exportedCells, cc, root);
    if (cc.rank() == root) {
        myExportList.reserve(grid.numCells());
        for (int i = 0; i < grid.numCells(); ++i) {
            myExportList.emplace_back(i, gIDtoRank[i], static_cast<char>(AttributeSet::owner));
        }
    } else {
        myImportList.reserve(importedCells.size());
        for (const auto& cell : importedCells) {
            myImportList.emplace_back(cell, root, static_cast<char>(AttributeSet::owner), -1);
        }
    }

    // get the distribution of wells
    std::vector<std::pair<std::string, bool>> parallel_wells;
    if (wells) {
        auto wellRanks = getWellRanks(gIDtoRank, wellConnections);
        parallel_wells = wellsOnThisRank(*wells, wellRanks, cc, root);
    }

    return std::make_tuple(std::move(gIDtoRank),
                           std::move(parallel_wells),
                           std::move(myExportList),
                           std::move(myImportList),
                           std::move(wellConnections));
}
#endif // HAVE_MPI

// explicit template instantiations

} // end namespace Opm
