/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <limits>


#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
#include <opm/grid/utility/OpmWellType.hpp>

#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <dune/common/parallel/indexset.hh>

namespace Dune
{
namespace cpgrid
{

void getNullVertexList(void* cpGridPointer, int numGlobalIdEntries,
                       int numLocalIdEntries, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR lids, int wgtDim,
                       float *objWgts, int *err)
{
    (void) cpGridPointer; (void) numGlobalIdEntries;
    (void) numLocalIdEntries; (void) gids; (void) lids; (void) objWgts;
    (void) wgtDim;
    // We do nothing as we pretend to not have any grid cells.
    *err = ZOLTAN_OK;
}

void getCpGridVertexList(void* cpGridPointer, int numGlobalIdEntries,
                         int numLocalIdEntries, ZOLTAN_ID_PTR gids,
                         ZOLTAN_ID_PTR lids, int wgtDim,
                         float *objWgts, int *err)
{
    (void) wgtDim; (void) objWgts;
    const Dune::CpGrid&  grid = *static_cast<const Dune::CpGrid*>(cpGridPointer);
    auto& globalIdSet         =  grid.globalIdSet();
    auto& localIdSet          =  grid.localIdSet();

    if ( numGlobalIdEntries != numLocalIdEntries || numGlobalIdEntries != 1 )
    {
        std::cerr<<"numGlobalIdEntries="<<numGlobalIdEntries<<" numLocalIdEntries="<<numLocalIdEntries<<" grid cells="
                 <<grid.numCells()<<std::endl;
        *err = ZOLTAN_FATAL;
        return;
    }
    int idx = 0;
    for (auto cell = grid.leafbegin<0>(), cellEnd = grid.leafend<0>();
         cell != cellEnd; ++cell)
    {
        gids[idx]   = globalIdSet.id(*cell);
        lids[idx++] = localIdSet.id(*cell);
    }
    *err = ZOLTAN_OK;
}

void getNullNumEdgesList(void *cpGridPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err)
{
    (void) sizeGID; (void) sizeLID; (void) numCells; (void) globalID;
    (void) localID; (void) numEdges; (void) cpGridPointer;
    // Pretend that there are no edges
    if (numCells!=0)
        *err = ZOLTAN_FATAL;
    else
        *err = ZOLTAN_OK;
}


int getNumberOfEdgesForSpecificCell(const CombinedGridWellGraph& graph, int localCellId) {
    const Dune::CpGrid&  grid = graph.getGrid();
    // Initial set of faces is the ones of the well completions
    auto edges = graph.getWellsGraph()[localCellId];
    // For the graph there is an edge only if the face has two neighbors.
    // Therefore we need to check each face
    for ( int local_face = 0; local_face < grid.numCellFaces(localCellId); ++local_face )
    {
        const int face  = grid.cellFace(localCellId, local_face); //Get the index of the face of the cell with the local index local_face
        const int face0 = grid.faceCell(face, 0);
        const int face1 = grid.faceCell(face, 1);

        if ( face0 != -1 && face1 != -1 )
        {
            if ( face0 != localCellId )
                edges.insert(face0);
            else
                edges.insert(face1);
        }
    }
    return edges.size();
}

void getCpGridWellsNumEdgesList(void *graphPointer, int sizeGID, int sizeLID,
                           int numCells,
                           ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int *numEdges, int *err)
{
    (void) globalID;
    const CombinedGridWellGraph& graph =
        *static_cast<CombinedGridWellGraph*>(graphPointer);
    const Dune::CpGrid&  grid = graph.getGrid();
    if ( sizeGID != 1 || sizeLID != 1 || numCells != grid.numCells() )
    {
        *err = ZOLTAN_FATAL;
        return;
    }
    for( int i = 0; i < numCells;  i++ )
    {
        numEdges[i] = getNumberOfEdgesForSpecificCell(graph, localID[i]);
    }
    *err = ZOLTAN_OK;
}

void getNullEdgeList(void *cpGridPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *numEdges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgtDim, float *ewgts, int *err)
{
    (void) cpGridPointer; (void) sizeGID; (void) sizeLID; (void) numCells;
    (void) globalID; (void) localID; (void) numEdges; (void) nborGID;
    (void) nborProc; (void) wgtDim; (void) ewgts;
    *err = ZOLTAN_OK;
}

template<typename ID, typename weightType>
void fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounter(const CombinedGridWellGraph& graph, const int localCellId, ID globalID, int& neighborCounter, ID& nborGID, weightType *ewgts) {
    const Dune::CpGrid&  grid = graph.getGrid();
    // First the strong edges of the well completions.
    auto wellEdges = graph.getWellsGraph()[localCellId];
    for( auto edge : wellEdges)
    {
        nborGID[neighborCounter] = edge;
        ewgts[neighborCounter++] = std::numeric_limits<weightType>::max();
    }

    // Now the ones of the grid that are not handled by the well completions
    for ( int local_face = 0 ; local_face < grid.numCellFaces(localCellId); ++local_face )
    {
        const int face  = grid.cellFace(localCellId, local_face);
        int otherCell   = grid.faceCell(face, 0);
        if ( otherCell == localCellId || otherCell == -1 )
        {
            otherCell = grid.faceCell(face, 1);
            if ( otherCell == localCellId || otherCell == -1 )
            {
                // no real face or already handled by well
                continue;
            }
            else
            {
                if ( wellEdges.find(otherCell) == wellEdges.end() )
                {
                    nborGID[neighborCounter] = globalID[otherCell];
                    ewgts[neighborCounter++] = graph.edgeWeight(face);
                }
                continue;
            }
        }
        if ( wellEdges.find(otherCell) == wellEdges.end() )
        {
            nborGID[neighborCounter] = globalID[otherCell];
            ewgts[neighborCounter++] = graph.edgeWeight(face);
        }
    }
}

void getCpGridWellsEdgeList(void *graphPointer, int sizeGID, int sizeLID,
                       int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int *numEdges,
                       ZOLTAN_ID_PTR nborGID, int *nborProc,
                       int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;
    assert(wgtDim==1);
    const CombinedGridWellGraph& graph =
        *static_cast<const CombinedGridWellGraph*>(graphPointer);
    const Dune::CpGrid&  grid = graph.getGrid();

    if ( sizeGID != 1 || sizeLID != 1 || numCells != grid.numCells() )
    {
        *err = ZOLTAN_FATAL;
        return;
    }
#ifndef NDEBUG
    int oldNeighborCounter = 0;
#endif
    int neighborCounter = 0;

    for( int cell = 0; cell < numCells;  cell++ )
    {
        fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounter(graph, localID[cell], globalID, neighborCounter, nborGID, ewgts);
#ifndef NDEBUG
        assert(neighborCounter-oldNeighborCounter==numEdges[cell]);
        oldNeighborCounter = neighborCounter;
#endif
    }

    const int myrank = grid.comm().rank();

    for ( int i = 0; i < neighborCounter; ++i )
    {
        nborProc[i] = myrank;
    }
#if defined(DEBUG) && false // The index set will not be initialized here!
    // The above relies heavily on the grid not being distributed already.
    // Therefore we check here that all cells are owned by us.
    GlobalLookupIndexSet<Dune::CpGrid::ParallelIndexSet>
        globalIdxSet(grid.getCellIndexSet(),
                     grid.numCells());
    for ( int cell = 0; cell < numCells;  cell++ )
    {
        if ( globalIdxSet.pair(cell)->local().attribute() !=
             Dune::CpGrid::ParallelIndexSet::LocalIndex::Attribute::owner )
        {
            *err = ZOLTAN_FATAL;
        }
    }
#endif
}

CombinedGridWellGraph::CombinedGridWellGraph(const CpGrid& grid,
                                             const std::vector<OpmWellType> * wells,
                                             const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                             const double* transmissibilities,
                                             bool pretendEmptyGrid,
                                             EdgeWeightMethod edgeWeightsMethod)
    : grid_(grid), transmissibilities_(transmissibilities), edgeWeightsMethod_(edgeWeightsMethod)
{
    if ( pretendEmptyGrid )
    {
        // wellsGraph not needed
        return;
    }
    wellsGraph_.resize(grid.numCells());
    if (wells) {
        const auto& cpgdim = grid.logicalCartesianSize();
        // create compressed lookup from cartesian.
        std::vector<int> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);

        for( int i=0; i < grid.numCells(); ++i )
        {
            cartesian_to_compressed[grid.globalCell()[i]] = i;
        }
        // If we're not taking the wells into account here we don't init the connections
        well_indices_.init(*wells, possibleFutureConnections, cpgdim, cartesian_to_compressed);
        std::vector<int>().swap(cartesian_to_compressed); // free memory.
        addCompletionSetToGraph();
    }
    if (edgeWeightsMethod == logTransEdgeWgt)
        findMaxMinTrans();
}


void setCpGridZoltanGraphFunctions(Zoltan_Struct *zz,
                                   const CombinedGridWellGraph& graph,
                                   bool pretendNull)
{
    Dune::CpGrid *gridPointer = const_cast<Dune::CpGrid*>(&graph.getGrid());
    if ( pretendNull )
    {
        Zoltan_Set_Num_Obj_Fn(zz, getNullNumCells, gridPointer);
        Zoltan_Set_Obj_List_Fn(zz, getNullVertexList, gridPointer);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getNullNumEdgesList, gridPointer);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getNullEdgeList, gridPointer);
    }
    else
    {
        CombinedGridWellGraph* graphPointer = const_cast<CombinedGridWellGraph*>(&graph);
        Zoltan_Set_Num_Obj_Fn(zz, getCpGridNumCells, gridPointer);
        Zoltan_Set_Obj_List_Fn(zz, getCpGridVertexList, gridPointer);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getCpGridWellsNumEdgesList, graphPointer);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getCpGridWellsEdgeList, graphPointer);
    }
}
// Explicit template instantiation for METIS
#if HAVE_METIS
template
void fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounter(const CombinedGridWellGraph&, const int, int*, int&, int*&, int*);

template
void fillNBORGIDAndWeightsForSpecificCellAndIncrementNeighborCounter(Dune::cpgrid::CombinedGridWellGraph const&, int, long*, int&, long*&, long*);

#endif


} // end namespace cpgrid
} // end namespace Dune
#endif // defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
