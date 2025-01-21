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

#ifndef GRAPH_OF_GRID_WRAPPERS_HEADER
#define GRAPH_OF_GRID_WRAPPERS_HEADER

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/grid/GraphOfGrid.hpp>
#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp> // defines Zoltan and null-callback-functions

namespace Opm {
/*
  This file contains wrappers for GraphOfGrid that satisfy interface
  requirements of graph partitioners like Zoltan and (TODO!) Metis.

  Additionally, parsing wells is done here.
*/

#if HAVE_MPI
/// \brief callback function for ZOLTAN_NUM_OBJ_FN
///
/// returns the number of vertices in the graph
int getGraphOfGridNumVertices(void* pGraph, int *err);

/// \brief callback function for ZOLTAN_OBJ_LIST_FN
///
/// fills the vector gIDs with vertex global IDs
///  and the vector objWeights with their weights
void getGraphOfGridVerticesList(void* pGraph,
               [[maybe_unused]] int dimGlobalID,
               [[maybe_unused]] int dimLocalID,
                                ZOLTAN_ID_PTR gIDs,
               [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                                int weightDim,
                                float *objWeights,
                                int *err);

/// \brief callback function for ZOLTAN_NUM_EDGES_MULTI_FN
///
/// takes the list of global IDs (gIDs) and fills (consecutively)
/// vector numEdges with the number of their edges
void getGraphOfGridNumEdges(void *pGraph,
           [[maybe_unused]] int dimGlobalID,
           [[maybe_unused]] int dimLocalID,
                            int numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            int *numEdges,
                            int *err);

/// \brief callback function for ZOLTAN_EDGE_LIST_MULTI_FN
///
/// takes the list of global IDs (gIDs) and fills (consecutively):
/// vector nborGIDs with the list of neighbors (all into 1 vector),
/// vector nborProc with neighbors' process numbers,
/// vector edgeWeights with edge weights.
/// The vector numEdges provides the number of edges for each gID
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
                            int *err);

/// \brief Register callback functions to Zoltan
template<typename Zoltan_Struct>
void setGraphOfGridZoltanGraphFunctions(Zoltan_Struct *zz,
                                        GraphOfGrid<Dune::CpGrid>& gog,
                                        bool pretendNull);
#endif

/// \brief Adds well to the GraphOfGrid
///
/// Translates wells' cartesian ID to global ID used in the graph.
/// Adding the well contracts vertices of the well into one vertex.
///
/// checkWellIntersections==true makes the algorithm check if wells
/// intersect and if their cell IDs are present in the graph.
/// Setting it to false makes the algorithm faster but leaves user
/// responsible for keeping wells disjoint.
void addFutureConnectionWells(GraphOfGrid<Dune::CpGrid>& gog,
                              const std::unordered_map<std::string, std::set<int>>& wells,
                              bool checkWellIntersections=true);

/// \brief Add WellConnections to the GraphOfGrid
///
/// checkWellIntersections==true makes the algorithm check if wells
/// intersect and if their cell IDs are present in the graph.
/// Setting it to false makes the algorithm faster but leaves user
/// responsible for keeping wells disjoint.
void addWellConnections(GraphOfGrid<Dune::CpGrid>& gog,
                        const Dune::cpgrid::WellConnections& wells,
                        bool checkWellIntersections=true);

/// \brief Correct gIDtoRank's data about well cells
///
/// gIDtoRank's entries come from Zoltan partitioner's export list
/// that does not contain all well cells. Default value is root's rank.
/// parameter root allows skipping wells that are correct.
void extendGIDtoRank(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<int>& gIDtoRank,
                     const int& root = -1);

#if HAVE_MPI
namespace Impl{
/// \brief Add cells to the import list and sort it
///
/// Helper function for extendAndSortExportAndImportLists.
void extendAndSortImportList(std::vector<std::tuple<int,int,char,int>>& importList,
                             const std::vector<int>& extraCells);

/// \brief Add well cells' global IDs to the root's export list and output cells missing in other rank's import lists
///
/// Helper function for extendAndSortExportAndImportLists.
/// On non-root ranks, it does nothing and returns an empty vector.
/// On root, exportList is extended by well cells that are hidden from the partitioner.
/// These cells are also collected and returned so they can be communicated to other ranks.
/// \return vector[rank][cell] Each entry contains vector of cells exported to that rank.
std::vector<std::vector<int>>
extendRootExportList(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<std::tuple<int,int,char>>& exportList,
                     int root,
                     const std::vector<int>& gIDtoRank);
/// \brief Communicate cells exported from root, needed for extending other rank's import lists
///
/// Helper function for extendAndSortExportAndImportLists.
/// Only the root rank has access to the grid and wells. Partitioning
/// with GraphOfGrid hides some well cells from the partitioner and
/// the importLists are incomplete. This function serves to communicate
/// the missing (well) cells to non-root ranks.
/// Input (the prepared container with cell IDs) is thus relevant only
/// on the root rank, whereas output (cells received by communication)
/// only on the non-root ranks.
/// Sending is nonblocking to quicken the root,
/// receiving is blocking - non-root ranks have nothing else to do anyway.
/// \param exportedCells Contains for each rank the cells that are exported there,
///                      empty on non-root ranks
/// \param cc Communication object
/// \param root The root's rank
/// \return Vector of cells necessary to extend this rank's import lists,
///         empty on the root rank
std::vector<int> communicateExportedCells(const std::vector<std::vector<int>>& exportedCells,
                                          const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                          int root);
} // end namespace Impl

/// \brief Add well cells' global IDs to the root's export and others' import list
///
/// Output of the partitioning is missing vertices that were contracted.
/// This function fills in omitted gIDs and gives them the properties
/// (like process number and ownership) of their representative cell (well ID).
/// Root is the only rank with information about wells, and communicates
/// the necessary information to other ranks.
/// On root ImportList has been already extended with all cells on the current rank.
void extendAndSortExportAndImportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                                       const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                       int root,
                                       std::vector<std::tuple<int, int, char>>& exportList,
                                       std::vector<std::tuple<int, int, char, int>>& importList,
                                       const std::vector<int>& gIDtoRank = {});
#endif // HAVE_MPI

/// \brief Find to which ranks wells were assigned
///
/// returns the vector of ranks, ordering is given by wellConnections
/// \param gIDtoRank Takes global ID and returns rank
/// \param wellConnections Has global IDs of well cells
std::vector<int> getWellRanks(const std::vector<int>& gIDtoRank,
                              const Dune::cpgrid::WellConnections& wellConnections);

#if HAVE_MPI
/// \brief Get rank-specific information about which wells are present
///
/// \param wells Vector of wells containing names (and other...)
/// \param wellRanks Tells on which (single) rank is the well placed,
///                  ordering in the vector is given by wells
/// \param cc The communication object
/// \param root Rank holding the information about the grid
/// @return vector of pairs string-bool that hold the name of the well
///         and whether it is situated on this process rank
///
/// This function only gets the information from wellRanks into proper
/// format to call computeParallelWells.
std::vector<std::pair<std::string, bool>>
wellsOnThisRank(const std::vector<Dune::cpgrid::OpmWellType>& wells,
                const std::vector<int>& wellRanks,
                const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                int root);

/// \brief Transform Zoltan output into tuples
///
/// \param gog GraphOfGrid, has ref. to CpGrid and knows how well-cells were contracted
/// \param cc Communication object
/// \param wells Used to extract well names
/// \param wellConnections Contains wells' global IDs, ordered as \param wells.
/// \param root Rank of the process executing the partitioning (usually 0)
/// \param numExport Number of cells in the export list
/// \param numImport Number of cells in the import list
/// \param exportLocalGids Unused. Partitioning is performed on root
///        process that has access to all cells.
/// \param exportGlobalGids Zoltan output: Global IDs of exported cells
/// \param exportToPart     Zoltan output: ranks to which cells are exported
/// \param importGlobalGids Zoltan output: Global IDs of cells imported to this rank
/// \return gIDtoRank A vector indexed by global ID storing the rank of cell
///         parallel_wells A vector of pairs wells.name and bool of "Is wells.name on this rank?"
///         myExportList vector of cells to be moved from this rank
///         myImportList vector of cells to be moved to this rank
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
                         const Id* importGlobalGids);

/// \brief Call Zoltan partitioner on GraphOfGrid
///
/// GraphOfGrid represents a well by one vertex, so wells can not be
/// spread over several processes.
std::tuple<std::vector<int>, std::vector<std::pair<std::string, bool>>,
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
                                  const std::map<std::string,std::string>& params);

/// \brief Make complete export lists from a vector holding destination rank for each global ID
///
/// Intended to use on the root, as other ranks can not construct gIDtoRank.
/// Export lists include root's cells that stay on root.
/// \param gIDtoRank a vector indexed by global ID holding rank to which it is exported
/// \param ccsize the number of ranks, cc.size()
/// \return vector of vectors, vector[i] holds a vector of global IDs exported to rank i
std::vector<std::vector<int> >
makeExportListsFromGIDtoRank(const std::vector<int>& gIDtoRank, int ccsize);

/// \brief Call serial Zoltan partitioner on GraphOfGrid
///
/// GraphOfGrid represents a well by one vertex, so wells can not be
/// spread over several processes.
std::tuple<std::vector<int>, std::vector<std::pair<std::string, bool>>,
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
                                        const std::map<std::string,std::string>& params);
#endif // HAVE_MPI

} // end namespace Opm

#endif // GRAPH_OF_GRID_WRAPPERS_HEADER
