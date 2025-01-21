//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//            Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#endif

#include "../CpGrid.hpp"
#include "ParentToChildrenCellGlobalIdHandle.hpp"
#include "ParentToChildCellToPointGlobalIdHandle.hpp"
#include <opm/grid/common/MetisPartition.hpp>
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/GraphOfGridWrappers.hpp>
//#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <opm/grid/common/GridPartitioning.hpp>
//#include <opm/grid/common/WellConnections.hpp>
#include <opm/grid/common/CommunicationUtils.hpp>

//#include <fstream>
//#include <iostream>
#include <iomanip>
#include <tuple>

namespace
{

#if HAVE_MPI

using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

template<typename Tuple, bool first>
void reserveInterface(const std::vector<Tuple>& list, Dune::CpGrid::InterfaceMap& interface,
                      const std::integral_constant<bool, first>&)
{
    std::map<int, std::size_t> proc_to_no_cells;
    for(const auto& entry: list)
    {
        ++proc_to_no_cells[std::get<1>(entry)];
    }
    for(const auto& proc: proc_to_no_cells)
    {
        auto& entry = interface[proc.first];
        if ( first )
            entry.first.reserve(proc.second);
        else
            entry.second.reserve(proc.second);
    }
}

void setupSendInterface(const std::vector<std::tuple<int, int, char> >& list, Dune::CpGrid::InterfaceMap& interface)
{
    reserveInterface(list, interface, std::integral_constant<bool, true>());
    int cellIndex=-1;
    int oldIndex = std::numeric_limits<int>::max();
    for(const auto& entry: list)
    {
        auto index = std::get<0>(entry);
        assert(oldIndex == std::numeric_limits<int>::max() || index >= oldIndex);

        if (index != oldIndex )
        {
            oldIndex = index;
            ++cellIndex;
        }
        interface[std::get<1>(entry)].first.add(cellIndex);
    }
}

void setupRecvInterface(const std::vector<std::tuple<int, int, char, int> >& list, Dune::CpGrid::InterfaceMap& interface)
{
    reserveInterface(list, interface, std::integral_constant<bool, false>());
    for(const auto& entry: list)
    {
        auto index = std::get<3>(entry);
        interface[std::get<1>(entry)].second.add(index);
    }
}
#endif // HAVE_MPI

/// Release memory resources from CpGrid::InterfaceMap.  Used as custom
/// deleter for std::shared_ptr<InterfaceMap>.
struct FreeInterfaces
{
#if !HAVE_MPI

    /// Release memory resources for InterfaceMap object
    ///
    /// \param[in] interfaces Object for which to release memory resources.
    void operator()([[maybe_unused]] Dune::CpGrid::InterfaceMap* interfaces) const
    {
        // Nothing to do in the sequential case as the CpGrid::InterfaceMap
        // handles interface deletion in its destructor in this case.
    }

#else // HAVE_MPI

    /// Release memory resources for InterfaceMap object
    ///
    /// \param[in] interfaces Object for which to release memory resources.
    void operator()(Dune::CpGrid::InterfaceMap* interfaces) const
    {
        if (interfaces == nullptr) {
            return;
        }

        for (auto& interface : *interfaces) {
            auto& [scatter, gather] = interface.second;
            scatter.free();
            gather.free();
        }
    }

#endif // HAVE_MPI
};
}

namespace Dune
{

CpGrid::CpGrid()
    : current_view_data_(),
      distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(data_));
    current_view_data_ = data_[0].get();
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
    
}

CpGrid::CpGrid(MPIHelper::MPICommunicator comm)
    : current_view_data_(),
      distributed_data_(),
      cell_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      point_scatter_gather_interfaces_(new InterfaceMap, FreeInterfaces{}),
      global_id_set_ptr_()
{
    data_.push_back(std::make_shared<cpgrid::CpGridData>(comm, data_));
    current_view_data_ = data_[0].get();
    current_data_ = &data_;
    global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
    
}

std::vector<int>
CpGrid::zoltanPartitionWithoutScatter([[maybe_unused]] const std::vector<cpgrid::OpmWellType>* wells,
                                      [[maybe_unused]] const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                                      [[maybe_unused]] const double* transmissibilities,
                                      [[maybe_unused]] const int numParts,
                                      [[maybe_unused]] const double zoltanImbalanceTol) const
{
#if HAVE_MPI && HAVE_ZOLTAN
    const auto met = EdgeWeightMethod(1);

    return cpgrid::zoltanGraphPartitionGridForJac(*this, wells, possibleFutureConnections, transmissibilities,
                                                  this->data_[0]->ccobj_, met,
                                                  0, numParts, zoltanImbalanceTol);
#else
    return std::vector<int>(this->numCells(), 0);
#endif
}


std::pair<bool, std::vector<std::pair<std::string,bool> > >
CpGrid::scatterGrid(EdgeWeightMethod method,
                    [[maybe_unused]] bool ownersFirst,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    [[maybe_unused]] const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                    [[maybe_unused]] bool serialPartitioning,
                    const double* transmissibilities,
                    [[maybe_unused]] bool addCornerCells,
                    int overlapLayers,
                    [[maybe_unused]] int partitionMethod,
                    double imbalanceTol,
                    [[maybe_unused]] bool allowDistributedWells,
                    [[maybe_unused]] const std::vector<int>& input_cell_part)
{
    // Silence any unused argument warnings that could occur with various configurations.
    static_cast<void>(wells);
    static_cast<void>(transmissibilities);
    static_cast<void>(overlapLayers);
    static_cast<void>(method);
    static_cast<void>(imbalanceTol);

    if(!distributed_data_.empty())
    {
        std::cerr<<"There is already a distributed version of the grid."
                 << " Maybe scatterGrid was called before?"<<std::endl;
        return std::make_pair(false, std::vector<std::pair<std::string,bool> >());
    }

    if (data_.size() > 1)
    {
        if (comm().rank() == 0)
        {
            OPM_THROW(std::logic_error, "Loadbalancing a grid with local grid refinement is not supported, yet.");
        }
        else
        {
            OPM_THROW_NOLOG(std::logic_error, "Loadbalancing a grid with local grid refinement is not supported, yet.");
        }
    }

#if HAVE_MPI
    auto& cc = data_[0]->ccobj_;

    if (cc.size() > 1)
    {
        std::vector<int> computedCellPart;
        std::vector<std::pair<std::string,bool>> wells_on_proc;
        std::vector<std::tuple<int,int,char>> exportList;
        std::vector<std::tuple<int,int,char,int>> importList;
        cpgrid::WellConnections wellConnections;

        auto inputNumParts = input_cell_part.size();
        inputNumParts = this->comm().max(inputNumParts);

        if ( inputNumParts > 0 )
        {
            std::vector<int> errors;
            std::vector<std::string> errorMessages =
                { "More parts than MPI Communicator can handle",
                  "Indices of parts need to zero starting",
                  "Indices of parts need to be consecutive",
                  "Only rank 0 should provide partitioning information for each cell"};

            std::set<int> existingParts;

            if (comm().rank() == 0)
            {
                for(const auto& part: input_cell_part)
                {
                    existingParts.insert(part);
                }
                if (*input_cell_part.rbegin() >= comm().size())
                {
                    errors.push_back(0);
                }

                int i = 0;
                if (*existingParts.begin() != i)
                {
                    errors.push_back(1);
                }
                for (const auto& part: existingParts)
                {
                    if (part != i++)
                    {
                        errors.push_back(2);
                        break;
                    }
                }
                if (std::size_t(size(0)) != input_cell_part.size())
                {
                    errors.push_back(3);
                }
            }
            auto size = errors.size();
            comm().broadcast(&size, 1, 0);
            errors.resize(size);

            if (!errors.empty())
            {
                comm().broadcast(errors.data(), size, 0);
                std::string message("Loadbalance: ");
                for ( const auto& e: errors)
                {
                    message.append(errorMessages[e]).append(". ");
                }
                if (comm().rank() == 0)
                {
                    OPM_THROW(std::logic_error, message);
                }
                else
                {
                    OPM_THROW_NOLOG(std::logic_error, message);
                }
            }


            // Partitioning given externally
            std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                cpgrid::createListsFromParts(*this, wells, possibleFutureConnections, nullptr, input_cell_part,
                                                   true);
        }
        else
        {
            if (partitionMethod == Dune::PartitionMethod::zoltan)
            {
#ifdef HAVE_ZOLTAN
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections)
                    = serialPartitioning
                    ? cpgrid::zoltanSerialGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams)
                    : cpgrid::zoltanGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN if useZoltan is true. Please install!");
#endif // HAVE_ZOLTAN
            }
            else if (partitionMethod == Dune::PartitionMethod::metis)
            {
#ifdef HAVE_METIS
                if (!serialPartitioning)
                    OPM_MESSAGE("Warning: Serial partitioning is set to false and METIS was selected to partition the grid, but METIS is a serial partitioner. Continuing with serial partitioning...");
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) = cpgrid::metisSerialGraphPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, allowDistributedWells, partitioningParams);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on METIS if useMetis is true. Please install!");
#endif // HAVE_METIS

            }
            else if (partitionMethod == Dune::PartitionMethod::zoltanGoG)
            {
#ifdef HAVE_ZOLTAN
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections)
                    = serialPartitioning
                    ? Opm::zoltanSerialPartitioningWithGraphOfGrid(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, partitioningParams)
                    : Opm::zoltanPartitioningWithGraphOfGrid(*this, wells, possibleFutureConnections, transmissibilities, cc, method, 0, imbalanceTol, partitioningParams);
#else
                OPM_THROW(std::runtime_error, "Parallel runs depend on ZOLTAN if useZoltan is true. Please install!");
#endif // HAVE_ZOLTAN
            }
            else
            {
                std::tie(computedCellPart, wells_on_proc, exportList, importList, wellConnections) =
                    cpgrid::vanillaPartitionGridOnRoot(*this, wells, possibleFutureConnections, transmissibilities, allowDistributedWells);
            }
        }
        comm().barrier();

        // first create the overlap
        auto noImportedOwner = addOverlapLayer(*this, computedCellPart, exportList, importList, cc, addCornerCells,
                                               transmissibilities);
        // importList contains all the indices that will be here.
        auto compareImport = [](const std::tuple<int,int,char,int>& t1,
                                const std::tuple<int,int,char,int>&t2)
        {
            return std::get<0>(t1) < std::get<0>(t2);
        };

        if ( ! ownersFirst )
        {
            // merge owner and overlap sorted by global index
            std::inplace_merge(importList.begin(), importList.begin()+noImportedOwner,
                               importList.end(), compareImport);
        }
        // assign local indices
        int localIndex = 0;
        for(auto&& entry: importList)
            std::get<3>(entry) = localIndex++;

        if ( ownersFirst )
        {
            // merge owner and overlap sorted by global index
            std::inplace_merge(importList.begin(), importList.begin()+noImportedOwner,
                               importList.end(), compareImport);
        }

        int procsWithZeroCells{};

        if (cc.rank()==0)
        {
            // Print some statistics without communication
            std::vector<int> ownedCells(cc.size(), 0);
            std::vector<int> overlapCells(cc.size(), 0);
            for (const auto& entry: exportList)
            {
                if(std::get<2>(entry) == AttributeSet::owner)
                {
                    ++ownedCells[std::get<1>(entry)];
                }
                else
                {
                    ++overlapCells[std::get<1>(entry)];
                }
            }

            for(const auto& cellsOnProc: ownedCells)
            {
                procsWithZeroCells += (cellsOnProc == 0);
            }
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributes " << data_[0]->size(0)
                 << " active cells on " << cc.size() << " processes as follows:\n";
            ostr << "  rank   owned cells   overlap cells   total cells\n";
            ostr << "--------------------------------------------------\n";
            for (int i = 0; i < cc.size(); ++i) {
                ostr << std::setw(6) << i
                     << std::setw(14) << ownedCells[i]
                     << std::setw(16) << overlapCells[i]
                     << std::setw(14) << ownedCells[i] + overlapCells[i] << "\n";
            }
            ostr << "--------------------------------------------------\n";
            ostr << "   sum";
            auto sumOwned = std::accumulate(ownedCells.begin(), ownedCells.end(), 0);
            ostr << std::setw(14) << sumOwned;
            auto sumOverlap = std::accumulate(overlapCells.begin(), overlapCells.end(), 0);
            ostr << std::setw(16) << sumOverlap;
            ostr << std::setw(14) << (sumOwned + sumOverlap) << "\n";
            Opm::OpmLog::info(ostr.str());
        }

        // Print well distribution
        std::vector<std::pair<int,int> > procWellPairs;

        // range filters would be nice here. so C++20.
        procWellPairs.reserve(std::count_if(std::begin(wells_on_proc),
                                            std::end(wells_on_proc),
                                            [](const std::pair<std::string, bool>& p){ return p.second; }));
        int wellIndex = 0;
        for ( const auto& well: wells_on_proc)
        {
            if ( well.second )
            {
                procWellPairs.emplace_back(cc.rank(), wellIndex);
            }
            ++wellIndex;
        }

        std::tie(procWellPairs, std::ignore) = Opm::gatherv(procWellPairs, cc, 0);

        if (cc.rank() == 0)
        {
            std::sort(std::begin(procWellPairs), std::end(procWellPairs),
                      [](const std::pair<int,int>& p1, const std::pair<int,int>& p2)
                      { return p1.second < p2.second;});
            std::ostringstream ostr;
            ostr << "\nLoad balancing distributed the wells as follows:\n"
                 << "  well name            ranks with perforated cells\n"
                 << "---------------------------------------------------\n";
            auto procWellPair = std::begin(procWellPairs);
            auto endProcWellPair = std::end(procWellPairs);
            int wellIdx = 0;
            for ( const auto& well: wells_on_proc)
            {
                ostr << std::setw(16) << well.first;
                while (procWellPair != endProcWellPair && procWellPair->second < wellIdx)
                {
                    ++procWellPair;
                }
                for ( ; procWellPair != endProcWellPair && procWellPair->second == wellIdx;
                      ++procWellPair)
                {
                    ostr << " "<< std::setw(7) << procWellPair->first;
                }
                ostr << "\n";
                ++wellIdx;
            }
            Opm::OpmLog::info(ostr.str());
        }

        procsWithZeroCells = cc.sum(procsWithZeroCells);

        if (procsWithZeroCells) {
            std::string msg = "At least one process has zero cells. Aborting. \n"
                " Try decreasing the imbalance tolerance with the argument \n"
                " --imbalance-tolerance. The current value is "
                + std::to_string(imbalanceTol);
            if (cc.rank()==0)
            {
                OPM_THROW(std::runtime_error, msg );
            }
            else
            {
                OPM_THROW_NOLOG(std::runtime_error, msg);
            }
        }


        // distributed_data should be empty at this point.
        distributed_data_.push_back(std::make_shared<cpgrid::CpGridData>(cc, distributed_data_));
        distributed_data_[0]->setUniqueBoundaryIds(data_[0]->uniqueBoundaryIds());

        // Just to be sure we assume that only master knows
        cc.broadcast(&distributed_data_[0]->use_unique_boundary_ids_, 1, 0);


        // Create indexset
        distributed_data_[0]->cellIndexSet().beginResize();
        for(const auto& entry: importList)
        {
            distributed_data_[0]->cellIndexSet()
                .add(std::get<0>(entry),ParallelIndexSet::LocalIndex(std::get<3>(entry),AttributeSet(std::get<2>(entry)), true));
        }
        distributed_data_[0]->cellIndexSet().endResize();
        // add an interface for gathering/scattering data with communication
        // forward direction will be scatter and backward gather
        // Interface will communicate from owner to all
        setupSendInterface(exportList, *cell_scatter_gather_interfaces_);
        setupRecvInterface(importList, *cell_scatter_gather_interfaces_);

        distributed_data_[0]->distributeGlobalGrid(*this,*this->current_view_data_, computedCellPart);
        (*global_id_set_ptr_).insertIdSet(*distributed_data_[0]);
        distributed_data_[0]-> index_set_.reset(new cpgrid::IndexSet(distributed_data_[0]->cell_to_face_.size(),
                                                                     distributed_data_[0]-> geomVector<3>().size()));



        current_view_data_ = distributed_data_[0].get();
        current_data_ = &distributed_data_;
        return std::make_pair(true, wells_on_proc);
    }
    else
    {
        std::cerr << "CpGrid::scatterGrid() only makes sense in a parallel run. "
                  << "This run only uses one process.\n";
        return std::make_pair(false, std::vector<std::pair<std::string,bool>>());
    }
#else // #if HAVE_MPI
    std::cerr << "CpGrid::scatterGrid() is non-trivial only with "
              << "MPI support and if the target Dune platform is "
              << "sufficiently recent.\n";
    return std::make_pair(false, std::vector<std::pair<std::string,bool>>());
#endif
}


void CpGrid::createCartesian(const std::array<int, 3>& dims,
                             const std::array<double, 3>& cellsize,
                             const std::array<int, 3>& shift)
{
    if ( current_view_data_->ccobj_.rank() != 0 )
    {
        // global grid only on rank 0
        current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                             current_view_data_->logical_cartesian_size_.size(),
                                             0);
        return;
    }

    // Make the grdecl format arrays.
    // Pillar coords.
    std::vector<double> coord;
    coord.reserve(6*(dims[0] + 1)*(dims[1] + 1));
    double bot = 0.0+shift[2]*cellsize[2];
    double top = (dims[2]+shift[2])*cellsize[2];
    // i runs fastest for the pillars.
    for (int j = 0; j < dims[1] + 1; ++j) {
        double y = (j+shift[1])*cellsize[1];
        for (int i = 0; i < dims[0] + 1; ++i) {
            double x = (i+shift[0])*cellsize[0];
            double pillar[6] = { x, y, bot, x, y, top };
            coord.insert(coord.end(), pillar, pillar + 6);
        }
    }
    std::vector<double> zcorn(8*dims[0]*dims[1]*dims[2]);
    const int num_per_layer = 4*dims[0]*dims[1];
    double* offset = &zcorn[0];
    for (int k = 0; k < dims[2]; ++k) {
        double zlow = (k+shift[2])*cellsize[2];
        std::fill(offset, offset + num_per_layer, zlow);
        offset += num_per_layer;
        double zhigh = (k+1+shift[2])*cellsize[2];
        std::fill(offset, offset + num_per_layer, zhigh);
        offset += num_per_layer;
    }
    std::vector<int> actnum(dims[0]*dims[1]*dims[2], 1);

    // Process them.
    grdecl g;
    g.dims[0] = dims[0];
    g.dims[1] = dims[1];
    g.dims[2] = dims[2];
    g.coord = &coord[0];
    g.zcorn = &zcorn[0];
    g.actnum = &actnum[0];
    using NNCMap = std::set<std::pair<int, int>>;
    using NNCMaps = std::array<NNCMap, 2>;
    NNCMaps nnc;
    current_view_data_->processEclipseFormat(g,
#if HAVE_ECL_INPUT
                                             nullptr,
#endif
                                             nnc, false, false, false, 0.0);
    // global grid only on rank 0
    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
                                         0);
}

const std::array<int, 3>& CpGrid::logicalCartesianSize() const
{
    // Temporary. For a grid with LGRs, we set the logical cartesian size of the LeafGridView as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return current_view_data_ -> logical_cartesian_size_;
}

const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData() const
{
    return *current_data_;
}

std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& CpGrid::currentData()
{
    return *current_data_;
}

const std::vector<int>& CpGrid::globalCell() const
{
    // Temporary. For a grid with LGRs, we set the globalCell() of the as the one for level 0.
    //            Goal: CartesianIndexMapper well-defined for CpGrid LeafView with LGRs.
    return currentData().back() -> global_cell_;
}

void CpGrid::computeGlobalCellLgr(const int& level, const std::array<int,3>& startIJK, std::vector<int>& global_cell_lgr)
{
    assert(level);
    for (const auto& element : elements(levelGridView(level))) {
        // Element belogns to an LGR, therefore has a father. Get IJK of the father in the level grid the father was born.
        // For CARFIN, parent cells belong to level 0.
        std::array<int,3> parentIJK = {0,0,0};
        currentData()[element.father().level()]->getIJK(element.father().index(), parentIJK);
        // Each parent cell has been refined in cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2] child cells.
        // element has certain 'position' inside its parent cell that can be described with 'IJK' indices, let's denote them by ijk,
        // where 0<= i < cells_per_dim[0], 0<= j < cells_per_dim[1], 0<= k < cells_per_dim[2].
        const auto& cells_per_dim = currentData()[level]->cells_per_dim_;
        //
        // Refined cell (here 'element') has "index in parent cell": k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
        // and it's stored in  cell_to_idxInParentCell_.
        auto idx_in_parent_cell =  currentData()[level]-> cell_to_idxInParentCell_[element.index()];
        // Find ijk.
        std::array<int,3> childIJK = currentData()[level]-> getIJK(idx_in_parent_cell, cells_per_dim);
        // The corresponding lgrIJK can be computed as follows:
        const std::array<int,3>& lgrIJK = { ( (parentIJK[0] - startIJK[0])*cells_per_dim[0] ) + childIJK[0],  // Shift parent index according to the startIJK of the LGR.
                                            ( (parentIJK[1] - startIJK[1])*cells_per_dim[1] ) + childIJK[1],
                                            ( (parentIJK[2] - startIJK[2])*cells_per_dim[2] ) + childIJK[2] };
        // Dimensions of the "patch of cells" formed when providing startIJK and endIJK for an LGR
        const auto& lgr_logical_cartesian_size = currentData()[level]->logical_cartesian_size_;
        global_cell_lgr[element.index()] = (lgrIJK[2]*lgr_logical_cartesian_size[0]*lgr_logical_cartesian_size[1]) + (lgrIJK[1]*lgr_logical_cartesian_size[0]) + lgrIJK[0];
    }
}

void CpGrid::computeGlobalCellLeafGridViewWithLgrs(std::vector<int>& global_cell_leaf)
{
    for (const auto& element: elements(leafGridView()))
    {
        // When refine via CpGrid::addLgrsUpdateGridView(/*...*/), level-grid to lookup global_cell_ is equal to level-zero-grid
        // In the context of allowed nested refinement, we lookup for the oldest ancestor, also belonging to level-zero-grid.
        auto ancestor = element.getOrigin();
        int origin_in_level_zero = ancestor.index();
        while (ancestor.level()>0){
            ancestor = ancestor.getOrigin();
            origin_in_level_zero = ancestor.getOrigin().index();
        }
        assert(ancestor.level()==0);
        assert(origin_in_level_zero < currentData().front()->size(0));
        global_cell_leaf[element.index()] = currentData().front()-> global_cell_[origin_in_level_zero];
    }
}

std::vector<std::unordered_map<std::size_t, std::size_t>> CpGrid::mapLocalCartesianIndexSetsToLeafIndexSet() const
{
    std::vector<std::unordered_map<std::size_t, std::size_t>> localCartesianIdxSets_to_leafIdx(maxLevel()+1); // Plus level 0
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getEquivLevelElem().index()];
        localCartesianIdxSets_to_leafIdx[element.level()][global_cell_level] = element.index();
    }
    return localCartesianIdxSets_to_leafIdx;
}

std::vector<std::array<int,2>> CpGrid::mapLeafIndexSetToLocalCartesianIndexSets() const
{
    std::vector<std::array<int,2>> leafIdx_to_localCartesianIdxSets(currentData().back()->size(0));
    for (const auto& element : elements(leafGridView())) {
        const auto& global_cell_level = currentData()[element.level()]->globalCell()[element.getEquivLevelElem().index()];
        leafIdx_to_localCartesianIdxSets[element.index()] = {element.level(), global_cell_level};
    }
    return leafIdx_to_localCartesianIdxSets;
}

void CpGrid::getIJK(const int c, std::array<int,3>& ijk) const
{
    current_view_data_->getIJK(c, ijk);
}

bool CpGrid::uniqueBoundaryIds() const
{
    return current_view_data_->uniqueBoundaryIds();
}

void CpGrid::setUniqueBoundaryIds(bool uids)
{
    current_view_data_->setUniqueBoundaryIds(uids);
}

std::string CpGrid::name() const
{
    return "CpGrid";
}

int CpGrid::maxLevel() const
{
    if (currentData().size() == 1){
        return 0; // "GLOBAL" grid is the unique one
    }
    else {  // There are multiple LGRs
        return this -> currentData().size() - 2; // last entry is leafView, and it starts in level 0 = GLOBAL grid.
    }
}

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lbegin (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, All_Partition>( *(*current_data_)[level], 0, true);
}
template typename CpGridTraits::template Codim<0>::LevelIterator CpGrid::lbegin<0>(int) const;
template typename CpGridTraits::template Codim<1>::LevelIterator CpGrid::lbegin<1>(int) const;
template typename CpGridTraits::template Codim<3>::LevelIterator CpGrid::lbegin<3>(int) const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LevelIterator CpGrid::lend (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, All_Partition>( *(*current_data_)[level], size(level, codim), true);
}
template typename CpGridTraits::template Codim<0>::LevelIterator CpGrid::lend<0>(int) const;
template typename CpGridTraits::template Codim<1>::LevelIterator CpGrid::lend<1>(int) const;
template typename CpGridTraits::template Codim<3>::LevelIterator CpGrid::lend<3>(int) const;

template<int codim>
typename CpGridTraits::template Codim<codim>::LeafIterator CpGrid::leafbegin() const
{
    return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
}
template typename CpGridTraits::template Codim<0>::LeafIterator CpGrid::leafbegin<0>() const;
template typename CpGridTraits::template Codim<1>::LeafIterator CpGrid::leafbegin<1>() const;
template typename CpGridTraits::template Codim<3>::LeafIterator CpGrid::leafbegin<3>() const;


template<int codim>
typename CpGridTraits::template Codim<codim>::LeafIterator CpGrid::leafend() const
{
    return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, size(codim), true);
}
template typename CpGridTraits::template Codim<0>::LeafIterator CpGrid::leafend<0>() const;
template typename CpGridTraits::template Codim<1>::LeafIterator CpGrid::leafend<1>() const;
template typename CpGridTraits::template Codim<3>::LeafIterator CpGrid::leafend<3>() const;

template<int codim, PartitionIteratorType PiType>
typename CpGridTraits::template Codim<codim>::template Partition<PiType>::LevelIterator CpGrid::lbegin (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, PiType>( *(*current_data_)[level], 0, true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<0,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<1,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lbegin<3,Dune::Ghost_Partition>(int) const;

template<int codim, PartitionIteratorType PiType>
typename CpGridTraits::template Codim<codim>::template Partition<PiType>::LevelIterator CpGrid::lend (int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return  cpgrid::Iterator<codim, PiType>( *(*current_data_)[level], size(level, codim), true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<0,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<0,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<0,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<0,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<0,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<0,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<1,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<1,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<1,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<1,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<1,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<1,Dune::Ghost_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LevelIterator
CpGrid::lend<3,Dune::Interior_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LevelIterator
CpGrid::lend<3,Dune::InteriorBorder_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LevelIterator
CpGrid::lend<3,Dune::Overlap_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LevelIterator
CpGrid::lend<3,Dune::OverlapFront_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LevelIterator
CpGrid::lend<3,Dune::All_Partition>(int) const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LevelIterator
CpGrid::lend<3,Dune::Ghost_Partition>(int) const;

template<int codim, PartitionIteratorType PiType>
typename CpGridFamily::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator CpGrid::leafbegin() const
{
    return cpgrid::Iterator<codim, PiType>(*current_view_data_, 0, true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<0,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<1,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafbegin<3,Dune::Ghost_Partition>() const;

template<int codim, PartitionIteratorType PiType>
typename CpGridFamily::Traits::template Codim<codim>::template Partition<PiType>::LeafIterator CpGrid::leafend() const
{
    return cpgrid::Iterator<codim, PiType>(*current_view_data_, size(codim), true);
}
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<0,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<0,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<0,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<0>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<0,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<1,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<1,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<1,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<1>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<1,Dune::Ghost_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Interior_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Interior_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::InteriorBorder_Partition>::LeafIterator
CpGrid::leafend<3,Dune::InteriorBorder_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Overlap_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Overlap_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::OverlapFront_Partition>::LeafIterator
CpGrid::leafend<3,Dune::OverlapFront_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::All_Partition>::LeafIterator
CpGrid::leafend<3,Dune::All_Partition>() const;
template typename CpGridTraits::template Codim<3>::template Partition<Dune::Ghost_Partition>::LeafIterator
CpGrid::leafend<3,Dune::Ghost_Partition>() const;

int CpGrid::size (int level, int codim) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return currentData()[level]->size(codim);
}

int CpGrid::size (int codim) const
{
    return current_view_data_->size(codim);
}

int CpGrid::size (int level, GeometryType type) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return currentData()[level]->size(type);
}

int CpGrid::size (GeometryType type) const
{
    return current_view_data_->size(type);
}

const CpGridFamily::Traits::GlobalIdSet& CpGrid::globalIdSet() const
{
    return  *global_id_set_ptr_;
}

const CpGridFamily::Traits::LocalIdSet& CpGrid::localIdSet() const
{
    return *global_id_set_ptr_;
}

const CpGridFamily::Traits::LevelIndexSet& CpGrid::levelIndexSet(int level) const
{
    if (level<0 || level>maxLevel())
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
    return *currentData()[level]->index_set_;
}

const CpGridFamily::Traits::LeafIndexSet& CpGrid::leafIndexSet() const
{
    return *current_view_data_->index_set_;
}

void CpGrid::globalRefine (int refCount)
{
    if (refCount < 0) {
        OPM_THROW(std::logic_error, "Invalid argument. Provide a nonnegative integer for global refinement.");
    }
    // Throw if the grid has been already partially refined, i.e., there exist coarse cells with more than 6 faces.
    // This is the case when a coarse cell has not been marked for refinement, but at least one of its neighboring cells
    // got refined. Therefore, the coarse face that they share got replaced by refined-faces. In this case, we do not
    // support yet global refinement.
    if(currentData().size() >1) {
        bool isOnlyGlobalRefined = true;
        for (std::size_t level = 0; level < currentData().size(); ++level) {
            // When the grid has been refined only via global refinement, i.e., each cell has been refined into 2x2x2 children cells,
            // then the quotient between the total amount of two consecutive refined level grids is equal to 8 = 2x2x2.
            isOnlyGlobalRefined = isOnlyGlobalRefined && ( (currentData()[level+1]->size(0)) / (currentData()[level]->size(0)) == 8 );
        }
        if (!isOnlyGlobalRefined) {
            OPM_THROW(std::logic_error, "Global refinement of a mixed grid with coarse and refined cells is not supported yet.");
        }
    }
    if (refCount>0) {
        const std::vector<std::array<int,3>>& cells_per_dim_vec = {{2,2,2}}; // Arbitrary chosen values.
        for (int refinedLevel = 0; refinedLevel < refCount; ++refinedLevel) {

            std::vector<int> assignRefinedLevel(current_view_data_-> size(0));
            const auto& preAdaptMaxLevel = this ->maxLevel();
            std::vector<std::string> lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
            bool isCARFIN = true;
            const std::array<int,3>& endIJK = currentData().back()->logicalCartesianSize();

            for(const auto& element: elements(this-> leafGridView())) {
                // Mark all the elements of the current leaf grid view for refinement
                mark(1, element);
                assignRefinedLevel[element.index()] = preAdaptMaxLevel +1;
            }

            preAdapt();
            adapt(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec, isCARFIN, {{0,0,0}}, {endIJK});
            postAdapt();
        }
    }
}

const std::vector< Dune :: GeometryType >& CpGrid::geomTypes( const int codim ) const
{
    return leafIndexSet().geomTypes( codim );
}

template <int codim>
cpgrid::Entity<codim> CpGrid::entity( const cpgrid::Entity< codim >& seed ) const
{
    return cpgrid::Entity<codim>( *(this->current_view_data_), seed );
}

template cpgrid::Entity<0> CpGrid::entity<0>( const cpgrid::Entity<0>&) const;
template cpgrid::Entity<3> CpGrid::entity<3>( const cpgrid::Entity<3>&) const;


/// \brief Size of the overlap on the leaf level
unsigned int CpGrid::overlapSize(int) const {
    return 1;
}


/// \brief Size of the ghost cell layer on the leaf level
unsigned int CpGrid::ghostSize(int) const {
    return 0;
}


/// \brief Size of the overlap on a given level
unsigned int CpGrid::overlapSize(int, int) const {
    return 1;
}


/// \brief Size of the ghost cell layer on a given level
unsigned int CpGrid::ghostSize(int, int) const {
    return 0;
}

unsigned int CpGrid::numBoundarySegments() const
{
    if( uniqueBoundaryIds() )
    {
        return current_view_data_->unique_boundary_ids_.size();
    }
    else
    {
        unsigned int numBndSegs = 0;
        const int num_faces = numFaces();
        for (int i = 0; i < num_faces; ++i) {
            cpgrid::EntityRep<1> face(i, true);
            if (current_view_data_->face_to_cell_[face].size() == 1) {
                ++numBndSegs;
            }
        }
        return numBndSegs;
    }
}

void CpGrid::setPartitioningParams(const std::map<std::string,std::string>& params)
{
    partitioningParams = params;
}

const typename CpGridTraits::Communication& Dune::CpGrid::comm () const
{
    return current_view_data_->ccobj_;
}

//

const std::vector<double>& CpGrid::zcornData() const {
    return current_view_data_->zcornData();
}

int CpGrid::numCells() const
{
    return current_view_data_->cell_to_face_.size();
}
/// \brief Get the number of faces.
int CpGrid::numFaces() const
{
    return current_view_data_->face_to_cell_.size();
}
/// \brief Get The number of vertices.
int CpGrid::numVertices() const
{
    return current_view_data_->geomVector<3>().size();
}

int CpGrid::numCellFaces(int cell) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
}

int CpGrid::cellFace(int cell, int local_index) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
}

const cpgrid::OrientedEntityTable<0,1>::row_type CpGrid::cellFaceRow(int cell) const
{
    return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
}

int CpGrid::faceCell(int face, int local_index) const
{
    // In the parallel case we store non-existent cells for faces along
    // the front region. Theses marked with index std::numeric_limits<int>::max(),
    // orientation might be arbitrary, though.
    cpgrid::OrientedEntityTable<1,0>::row_type r
        = current_view_data_->face_to_cell_[cpgrid::EntityRep<1>(face, true)];
    bool a = (local_index == 0);
    bool b = r[0].orientation();
    bool use_first = a ? b : !b;
    // The number of valid cells.
    int r_size = r.size();
    // In the case of only one valid cell, this is the index of it.
    int index = 0;
    if(r[0].index()==std::numeric_limits<int>::max()){
        assert(r_size==2);
        --r_size;
        index=1;
    }
    if(r.size()>1 && r[1].index()==std::numeric_limits<int>::max())
    {
        assert(r_size==2);
        --r_size;
    }
    if (r_size == 2) {
        return use_first ? r[0].index() : r[1].index();
    } else {
        return use_first ? r[index].index() : -1;
    }
}

int CpGrid::numCellFaces() const
{
    return current_view_data_->cell_to_face_.dataSize();
}

int CpGrid::numFaceVertices(int face) const
{
    return current_view_data_->face_to_point_[face].size();
}

int CpGrid::faceVertex(int face, int local_index) const
{
    return current_view_data_->face_to_point_[face][local_index];
}

Dune::cpgrid::Intersection CpGrid::getParentIntersectionFromLgrBoundaryFace(const Dune::cpgrid::Intersection& intersection) const
{
    if ( intersection.neighbor()) {
        if ((intersection.inside().level() != intersection.outside().level())) {
            // one coarse and one refined neighboring cell
            /** Now, it could also be two refined cells. In that case, any of them will fit to search for the parent face */
            const auto& cellIn = intersection.inside();
            const auto& cellOut = intersection.outside();

            // Identify the coarse and the refined neighboring cell
            const auto coarseCell =  (cellIn.level() == 0) ? cellIn : cellOut;
            const auto refinedCell =  (coarseCell == cellIn) ? cellOut : cellIn;
            assert(coarseCell.level() != refinedCell.level());

            // Get parent cell (on level zero) of the refined cell
            const auto& parentCell = refinedCell.father();
            assert(refinedCell.father().level() == 0);

            // Get the index inside and orientation from the leaf grid (refined) face
            const auto& intersectionIdxInInside = intersection.indexInInside();

            for(const auto& parentIntersection : intersections(this->levelGridView(0), parentCell)){
                // Get the inInsideIdx and orientation from the parent intersection
                const auto& parentIdxInInside = parentIntersection.indexInInside();
                if (parentIdxInInside == intersectionIdxInInside) {
                    return parentIntersection;
                }
            }
        }
        OPM_THROW(std::invalid_argument, "Parent intersection not found for face with index: " + std::to_string(intersection.id()) +
                  " and index in inside: " + std::to_string(intersection.indexInInside()));
    }
    OPM_THROW(std::invalid_argument, "Face is on the boundary of the grid");
}

bool CpGrid::nonNNCsSelectedCellsLGR( const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const
{
    // Non neighboring connections: Currently, adding LGRs whose cells have NNCs is not supported yet.
    // Find out which (ACTIVE) elements belong to the block cells defined by startIJK and endIJK values
    // and check if they have NNCs.
    // Note: at this point, the (level zero) grid might be already distributed, before adding the LGRs. Therefore,
    // we get the active index of each element and its corresponding ijk from the global grid, to be able to determine
    // if the cell bolengs to certain block of cells selected for refinement (comparing element's ijk with start/endIJK values).
    // It is not correct to make the comparasion with element.index() and minimum/maximum index of each block of cell, since
    // element.index() is local and the block of cells are defined with global values.
    for(const auto& element: elements(this->leafGridView())) {
        std::array<int,3> ijk;
        getIJK(element.index(), ijk);
        for (std::size_t level = 0; level < startIJK_vec.size(); ++level) {
            bool belongsToLevel = ( (ijk[0] >= startIJK_vec[level][0]) && (ijk[0] < endIJK_vec[level][0]) )
                && ( (ijk[1] >= startIJK_vec[level][1]) && (ijk[1] < endIJK_vec[level][1]) )
                && ( (ijk[2] >= startIJK_vec[level][2]) && (ijk[2] < endIJK_vec[level][2]) );
            if(belongsToLevel) {
                // Check that the cell to be marked for  refinement has no NNC (no neighbouring connections).
                if (this->currentData().back()->hasNNCs({element.index()})){
                    return false;
                }
            }
        }
    }
    return true;
}

template<class T>
void CpGrid::computeOnLgrParents(const std::vector<std::array<int,3>>& startIJK_vec,
                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                 T func)
{
    // Find out which (ACTIVE) elements belong to the block cells defined by startIJK and endIJK values.
    for(const auto& element: elements(this->leafGridView())) {
        std::array<int,3> ijk;
        getIJK(element.index(), ijk);
        for (std::size_t level = 0; level < startIJK_vec.size(); ++level) {
            bool belongsToLevel = true;
            for (int c = 0; c < 3; ++c) {
                belongsToLevel = belongsToLevel && ( (ijk[c] >= startIJK_vec[level][c]) && (ijk[c] < endIJK_vec[level][c]) );
                if (!belongsToLevel)
                    break;
            }
            if(belongsToLevel) {
                func(element, level);
            }
        }
    }
}

void CpGrid::detectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                              const std::vector<std::array<int,3>>& endIJK_vec,
                              std::vector<int>& lgr_with_at_least_one_active_cell)
{
    auto markLgr = [&lgr_with_at_least_one_active_cell]([[maybe_unused]] const cpgrid::Entity<0>& element, int level)
    {
        // shifted since starting grid is level 0, and refined grids levels are >= 1.
        lgr_with_at_least_one_active_cell[level] = 1;
    };
    computeOnLgrParents(startIJK_vec, endIJK_vec, markLgr);
}

void CpGrid::markElemAssignLevelDetectActiveLgrs(const std::vector<std::array<int,3>>& startIJK_vec,
                                                 const std::vector<std::array<int,3>>& endIJK_vec,
                                                 std::vector<int>& assignRefinedLevel,
                                                 std::vector<int>& lgr_with_at_least_one_active_cell)
{
    auto assignAndDetect = [this, &assignRefinedLevel, &lgr_with_at_least_one_active_cell](const cpgrid::Entity<0>& element, int level)
    {
        mark(1, element);
        assignRefinedLevel[element.index()] = level+1;
        // shifted since starting grid is level 0, and refined grids levels are >= 1.
        lgr_with_at_least_one_active_cell[level] = 1;
    };
    computeOnLgrParents(startIJK_vec, endIJK_vec, assignAndDetect);
}

void CpGrid::predictMinCellAndPointGlobalIdPerProcess([[maybe_unused]] const std::vector<int>& assignRefinedLevel,
                                                      [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                                      [[maybe_unused]] const std::vector<int>& lgr_with_at_least_one_active_cell,
                                                      [[maybe_unused]] int& min_globalId_cell_in_proc,
                                                      [[maybe_unused]] int& min_globalId_point_in_proc) const
{
#if HAVE_MPI
    // Maximum global id from level zero. (Then, new entities get global id values greater than max_globalId_levelZero).
    // Recall that only cells and points are taken into account; faces are ignored (do not have any global id).
    auto max_globalId_levelZero = comm().max(current_data_->front()->global_id_set_->getMaxGlobalId());

    // Predict how many new cell ids per process are needed.
    std::vector<std::size_t> cell_ids_needed_by_proc(comm().size());
    std::size_t local_cell_ids_needed = 0;
    for ( const auto& element : elements( levelGridView(0), Dune::Partitions::interior) ) {
        // Get old mark (from level zero). After calling adapt, all marks are set to zero.
        bool hasBeenMarked = currentData().front()->getMark(element) == 1;
        if ( hasBeenMarked ) {
            const auto& level = assignRefinedLevel[element.index()];
            // Shift level (to level -1) since cells_per_dim_vec stores number of subdivisions in each direction (xyz)
            // per parent cell, per level, starting from level 1, ..., maxLevel.
            local_cell_ids_needed += cells_per_dim_vec[level-1][0]*cells_per_dim_vec[level-1][1]*cells_per_dim_vec[level-1][2];
        }
    }
    comm().allgather(&local_cell_ids_needed, 1, cell_ids_needed_by_proc.data());

    // Overestimate ('predict') how many new point ids per process are needed.
    // Assign for all partition type points a 'candidate of global id' (unique in each process).
    std::vector<std::size_t> point_ids_needed_by_proc(comm().size());
    std::size_t local_point_ids_needed = 0;
    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level){
        if(lgr_with_at_least_one_active_cell[level-1]>0) {
            // Amount of local_point_ids_needed might be overestimated.
            for (const auto& point : vertices(levelGridView(level))){
                // If point coincides with an existing corner from level zero, then it does not need a new global id.
                if ( !(*current_data_)[level]->corner_history_.empty() ) {
                    const auto& bornLevel_bornIdx =  (*current_data_)[level]->corner_history_[point.index()];
                    if (bornLevel_bornIdx[0] == -1)  { // Corner is new-> it needs a new(candidate) global id
                        local_point_ids_needed += 1;
                    }
                }
            }
        }
    }
    comm().allgather(&local_point_ids_needed, 1, point_ids_needed_by_proc.data());

    auto expected_max_globalId_cell = std::accumulate(cell_ids_needed_by_proc.begin(),
                                                      cell_ids_needed_by_proc.end(),
                                                      max_globalId_levelZero + 1);
    min_globalId_cell_in_proc = std::accumulate(cell_ids_needed_by_proc.begin(),
                                                cell_ids_needed_by_proc.begin()+comm().rank(),
                                                max_globalId_levelZero + 1);
    min_globalId_point_in_proc = std::accumulate(point_ids_needed_by_proc.begin(),
                                                 point_ids_needed_by_proc.begin()+ comm().rank(),
                                                 expected_max_globalId_cell+1);
#endif
}

void CpGrid::assignCellIdsAndCandidatePointIds( std::vector<std::vector<int>>& localToGlobal_cells_per_level,
                                                 std::vector<std::vector<int>>& localToGlobal_points_per_level,
                                                 int min_globalId_cell_in_proc,
                                                 int min_globalId_point_in_proc,
                                                 const std::vector<std::array<int,3>>& cells_per_dim_vec ) const
{
    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
        localToGlobal_cells_per_level[level-1].resize((*current_data_)[level]-> size(0));
        localToGlobal_points_per_level[level-1].resize((*current_data_)[level]-> size(3));
        // Notice that in general, (*current_data_)[level]-> size(0) != local owned cells/points.

        // Global ids for cells (for owned cells)
        for (const auto& element : elements(levelGridView(level))) {
            // At his point, partition_type_indicator_ of refined level grids is not set. However, for refined cells,
            // element.partitionType() returns the partition type of the parent cell. Therefore, all child cells of
            // an interior/overlap parent cell are also interior/overlap.
            if (element.partitionType() == InteriorEntity) {
                localToGlobal_cells_per_level[level - 1][element.index()] = min_globalId_cell_in_proc;
                ++min_globalId_cell_in_proc;
            }
        }
        for (const auto& point : vertices(levelGridView(level))) {
            // Check if it coincides with a corner from level zero. In that case, no global id is needed.
            const auto& bornLevel_bornIdx =  (*current_data_)[level]->corner_history_[point.index()];
            if (bornLevel_bornIdx[0] != -1)  { // Corner in the refined grid coincides with a corner from level 0.
                // Therefore, search and assign the global id of the previous existing equivalent corner.
                const auto& equivPoint = cpgrid::Entity<3>(*( (*current_data_)[bornLevel_bornIdx[0]]), bornLevel_bornIdx[1], true);
                localToGlobal_points_per_level[level-1][point.index()] = current_data_->front()->global_id_set_->id( equivPoint );
            }
            else {
                // Assign CANDIDATE global id to (all partition type) points that do not coincide with
                // any corners from level zero.
                // TO DO after invoking this method: make a final decision on a unique id for points of refined level grids,
                // via a communication step.
                localToGlobal_points_per_level[level-1][point.index()] = min_globalId_point_in_proc;
                ++min_globalId_point_in_proc;
            }
        }
    }
}

void CpGrid::selectWinnerPointIds([[maybe_unused]] std::vector<std::vector<int>>&  localToGlobal_points_per_level,
                                  [[maybe_unused]] const std::vector<std::tuple<int,std::vector<int>>>& parent_to_children,
                                  [[maybe_unused]] const std::vector<std::array<int,3>>& cells_per_dim_vec) const
{
#if HAVE_MPI
    // To store cell_to_point_ information of all refined level grids.
    std::vector<std::vector<std::array<int,8>>> level_cell_to_point(cells_per_dim_vec.size());
    // To decide which "candidate" point global id wins, the rank is stored. The smallest ranks wins,
    // i.e., the other non-selected candidates get rewritten with the values from the smallest (winner) rank.
    std::vector<std::vector<int>> level_winning_ranks(cells_per_dim_vec.size());

    for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {

        level_cell_to_point[level -1] = currentData()[level]->cell_to_point_;
        // Set std::numeric_limits<int>::max() to make sure that, during communication, the rank of the interior cell
        // wins (int between 0 and comm().size()).
        level_winning_ranks[level-1].resize(currentData()[level]->size(3), std::numeric_limits<int>::max());

        for (const auto& element : elements(levelGridView(level))) {
            // For interior cells, rewrite the rank value - later used in "point global id competition".
            if (element.partitionType() == InteriorEntity) {
                for (const auto& corner : currentData()[level]->cell_to_point_[element.index()]){
                    int rank = comm().rank();
                    level_winning_ranks[level -1][corner] = rank;
                }
            }
        }
    }
    ParentToChildCellToPointGlobalIdHandle parentToChildCellToPointGlobalId_handle(comm(),
                                                                                   parent_to_children,
                                                                                   level_cell_to_point,
                                                                                   level_winning_ranks,
                                                                                   localToGlobal_points_per_level);
    currentData().front()->communicate(parentToChildCellToPointGlobalId_handle,
                                       Dune::InteriorBorder_All_Interface,
                                       Dune::ForwardCommunication );
#endif
}

void CpGrid::populateCellIndexSetRefinedGrid([[maybe_unused]] int level)
{
#if HAVE_MPI
    const auto& level_global_id_set =  (*current_data_)[level]->global_id_set_;
    auto& level_index_set =  currentData()[level]->cellIndexSet();

    level_index_set.beginResize();
    
    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    //
    // For refined level grids, we check if the parent cell is fully interior. Then, its children won't be seen
    // by any other process. Therefore, the boolean is set to false.

    for(const auto& element : elements(levelGridView(level))) {
        if ( element.partitionType() == InteriorEntity) {

            // Check if it has an overlap neighbor
            bool parentFullyInterior = true;
            const auto& parent_cell = element.father();

            for (const auto& intersection : intersections(levelGridView(parent_cell.level()), parent_cell)) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        parentFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> true
                        level_index_set.add( level_global_id_set->id(element),
                                             ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(parentFullyInterior) { // Fully interior cells do not appear in any other process -> false
                level_index_set.add( level_global_id_set->id(element),
                                     ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(element.partitionType() == OverlapEntity);
            level_index_set.add( level_global_id_set->id(element),
                                 ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }

    level_index_set.endResize();

    currentData()[level]->cellRemoteIndices().template rebuild<false>(); // Do not ignore bool in ParallelIndexSet!
#endif
}

void CpGrid::populateCellIndexSetLeafGridView()
{
#if HAVE_MPI
    auto& leaf_index_set =  (*current_data_).back()->cellIndexSet();

    leaf_index_set.beginResize();

    // ParallelIndexSet::LocalIndex( elemIdx, attribute /* owner or copy */, true/false)
    // The bool indicates whether the global index might exist on other processes with other attributes.
    // RemoteIndices::rebuild has a Boolean as a template parameter telling the method whether to ignore this
    // Boolean on the indices or not when building.
    for(const auto& element : elements(leafGridView())) {
        const auto& elemPartitionType = element.getEquivLevelElem().partitionType();
        if ( elemPartitionType == InteriorEntity) {
            // Check if it has an overlap neighbor
            bool isFullyInterior = true;
            for (const auto& intersection : intersections(leafGridView(), element)) {
                if ( intersection.neighbor() ) {
                    const auto& neighborPartitionType = intersection.outside().getEquivLevelElem().partitionType();
                    // To help detection of fully interior cells, i.e., without overlap neighbors
                    if (neighborPartitionType == OverlapEntity )  {
                        isFullyInterior = false;
                        // Interior cells with overlap neighbor may appear in other processess -> false
                        leaf_index_set.add(globalIdSet().id(element),
                                           ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), true));
                        // Store it only once
                        break;
                    }
                }
            }
            if(isFullyInterior) { // Fully interior cells do not appear in any other process -> false
                leaf_index_set.add(globalIdSet().id(element),
                                   ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::owner), false));
            }
        }
        else { // overlap cell
            assert(elemPartitionType == OverlapEntity);
            leaf_index_set.add(globalIdSet().id(element),
                               ParallelIndexSet::LocalIndex(element.index(), AttributeSet(AttributeSet::copy), true));
        }
    }
    leaf_index_set.endResize();

    (*current_data_).back()->cellRemoteIndices().template rebuild<false>();

#endif
}

void CpGrid::populateLeafGlobalIdSet()
{
    // Global id for the cells in leaf grid view
    std::vector<int> leafCellIds(current_data_->back()->size(0));
    for(const auto& element: elements(leafGridView())){
        // Notice that for level zero cells the global_id_set_ is given, for refined level grids was defined
        // under the assumption of each lgr being fully contained in the interior of a process.
        // Therefore, it is not needed here to distingish between owned and overlap cells.
        auto equivElem = element.getEquivLevelElem();
        leafCellIds[element.index()] = (*current_data_)[element.level()]->global_id_set_->id(equivElem);
    }

    // Global id for the faces in leaf grid view. Empty vector (Entity<1> not supported for CpGrid).
    std::vector<int> leafFaceIds{};

    // Global id for the points in leaf grid view
    std::vector<int> leafPointIds(current_data_->back()->size(3));
    for(const auto& point : vertices(leafGridView())){
        const auto& level_pointLevelIdx = current_data_->back()->corner_history_[point.index()];
        assert(level_pointLevelIdx[0] != -1);
        assert(level_pointLevelIdx[1] != -1);
        const auto& pointLevelEntity =  cpgrid::Entity<3>(*( (*current_data_)[level_pointLevelIdx[0]]), level_pointLevelIdx[1], true);
        leafPointIds[point.index()] = (*current_data_)[level_pointLevelIdx[0]]->global_id_set_->id(pointLevelEntity);
    }

    current_data_->back()->global_id_set_->swap(leafCellIds, leafFaceIds, leafPointIds);
}

double CpGrid::cellCenterDepth(int cell_index) const
{
    // Here cell center depth is computed as a raw average of cell corner depths.
    // This generally gives slightly different results than using the cell centroid.
    double zz = 0.0;
    const int nv = current_view_data_->cell_to_point_[cell_index].size();
    const int nd = 3;
    for (int i=0; i<nv; ++i) {
        zz += vertexPosition(current_view_data_->cell_to_point_[cell_index][i])[nd-1];
    }
    return zz/nv;
}

const Dune::FieldVector<double,3> CpGrid::faceCenterEcl(int cell_index, int face, const Dune::cpgrid::Intersection& intersection) const
{
    // This method is an alternative to the method faceCentroid(...).
    // The face center is computed as a raw average of cell corners.
    // For faulted cells this gives different results then average of face nodes
    // that seems to agree more with eclipse.
    // This assumes the cell nodes are ordered
    // 6---7
    // | T |
    // 4---5
    //   2---3
    //   | B |
    //   0---1

    // this follows the DUNE reference cube
    static const int faceVxMap[ 6 ][ 4 ] = { {0, 2, 4, 6}, // face 0 - I_FACE false
                                             {1, 3, 5, 7}, // face 1 - I_FACE true
                                             {0, 1, 4, 5}, // face 2 - J_FACE false
                                             {2, 3, 6, 7}, // face 3 - J_FACE true
                                             {0, 1, 2, 3}, // face 4 - K_FACE false
                                             {4, 5, 6, 7}  // face 5 - K_FACE true
    };


    assert (current_view_data_->cell_to_point_[cell_index].size() == 8);
    Dune::FieldVector<double,3> center(0.0);

    bool isCoarseCellInside = (intersection.inside().level() == 0);
    bool isCoarseCellOutside = false;
    if (intersection.neighbor()){
        isCoarseCellOutside = (intersection.outside().level() == 0);
    }
    bool twoCoarseNeighboringCells = isCoarseCellInside && isCoarseCellOutside;
    bool isOnGridBoundary_coarseNeighboringCell = intersection.boundary() && isCoarseCellInside && (!intersection.neighbor());

    // For CpGrid with LGRs, a refined face with a coarse neighboring cell and a refined neighboring cell
    // (that is when the face belongs to the boundary of an LGR and is located in the interior of the grid),
    // unfortunately leads us to a different order of the faces, in cell_to_face_, depending on if the
    // neighboring cell, here with cell_index index, is the coarse one or the refined one. Preceisely,
    // cell_to_face_[cell_index - coarse neighboring cell] = { left, right, front, back, bottom, top} = {0,1,2,3,4,5} with
    // the notation above, and
    // cell_to_face_[cell_index - refined neighboring cell] = {bottom, front, left, right, back, top} = {2,3,1,4,0,5} with
    // the notation used in faceVxMap.

    for( int i=0; i<4; ++i ) {
        if ((maxLevel() == 0) || twoCoarseNeighboringCells || isOnGridBoundary_coarseNeighboringCell) {
            center += vertexPosition(current_view_data_->cell_to_point_[cell_index][ faceVxMap[ face ][ i ] ]);
        }
        else { //  (refined) intersection with one coarse neighboring cell and one refined neighboring cell
            center += vertexPosition(current_view_data_->face_to_point_[intersection.id()][i]);
        }
    }

    for (int i=0; i<3; ++i) {
        center[i] /= 4;
    }
    return center;

}

const Dune::FieldVector<double,3> CpGrid::faceAreaNormalEcl(int face) const
{
    // same implementation as ResInsight
    const int nd = Dune::FieldVector<double,3>::dimension;
    const int nv =  numFaceVertices(face);
    switch (nv)
    {
    case 0:
    case 1:
    case 2:
        {
            return Dune::FieldVector<double,3>(0.0);
        }
        break;
    case 3:
        {
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][0])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][1])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> areaNormal = cross(a,b);
            for (int i=0; i<nd; ++i) {
                areaNormal[i] /= 2;
            }
            return areaNormal;
        }
        break;
    case 4:
        {
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][0])
                - vertexPosition(current_view_data_->face_to_point_[face][2]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][1])
                - vertexPosition(current_view_data_->face_to_point_[face][3]);
            Dune::FieldVector<double,3> areaNormal = cross(a,b);
            areaNormal *= 0.5;
            return areaNormal;
        }
        break;
    default:
        {
            int h = (nv - 1)/2;
            int k = (nv % 2) ? 0 : nv - 1;

            Dune::FieldVector<double,3> areaNormal(0.0);
            // First quads
            for (int i = 1; i < h; ++i)
            {
                Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][2*i])
                    - vertexPosition(current_view_data_->face_to_point_[face][0]);
                Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][2*i+1])
                    - vertexPosition(current_view_data_->face_to_point_[face][2*i-1]);
                areaNormal += cross(a,b);
            }

            // Last triangle or quad
            Dune::FieldVector<double,3> a = vertexPosition(current_view_data_->face_to_point_[face][2*h])
                - vertexPosition(current_view_data_->face_to_point_[face][0]);
            Dune::FieldVector<double,3> b = vertexPosition(current_view_data_->face_to_point_[face][k])
                - vertexPosition(current_view_data_->face_to_point_[face][2*h-1]);
            areaNormal += cross(a,b);

            areaNormal *= 0.5;

            return areaNormal;
        }

    }
}

const Dune::FieldVector<double,3>& CpGrid::vertexPosition(int vertex) const
{
    return current_view_data_->geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
}

double CpGrid::faceArea(int face) const
{
    return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::faceCentroid(int face) const
{
    return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
}

const Dune::FieldVector<double,3>& CpGrid::faceNormal(int face) const
{
    return current_view_data_->face_normals_.get(face);
}

double CpGrid::cellVolume(int cell) const
{
    return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
}

const Dune::FieldVector<double,3>& CpGrid::cellCentroid(int cell) const
{
    return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
}

CpGrid::CentroidIterator<0> CpGrid::beginCellCentroids() const
{
    return CentroidIterator<0>(current_view_data_->geomVector<0>().begin());
}

CpGrid::CentroidIterator<1> CpGrid::beginFaceCentroids() const
{
    return CentroidIterator<1>(current_view_data_->geomVector<1>().begin());
}

const std::vector<int>& CpGrid::sortedNumAquiferCells() const{
           return current_view_data_->sortedNumAquiferCells();
}

int CpGrid::boundaryId(int face) const
{
    // Note that this relies on the following implementation detail:
    // The grid is always construct such that the faces where
    // orientation() returns true are oriented along the positive IJK
    // direction. Oriented means that the first cell attached to face
    // has the lower index.
    int ret = 0;
    cpgrid::EntityRep<1> f(face, true);
    if (current_view_data_->face_to_cell_[f].size() == 1) {
        if (current_view_data_->uniqueBoundaryIds()) {
            // Use the unique boundary ids.
            ret = current_view_data_->unique_boundary_ids_[f];
        } else {
            // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
            const bool normal_is_in =
                !(current_view_data_->face_to_cell_[f][0].orientation());
            enum face_tag tag = current_view_data_->face_tag_[f];
            switch (tag) {
            case I_FACE:
                //                   LEFT : RIGHT
                ret = normal_is_in ? 1    : 2; // min(I) : max(I)
                break;
            case J_FACE:
                //                   BACK : FRONT
                ret = normal_is_in ? 3    : 4; // min(J) : max(J)
                break;
            case K_FACE:
                // Note: TOP at min(K) as 'z' measures *depth*.
                //                   TOP  : BOTTOM
                ret = normal_is_in ? 5    : 6; // min(K) : max(K)
                break;
            case NNC_FACE:
                // This should not be possible, as NNC "faces" always
                // have two cell neighbours.
                OPM_THROW(std::logic_error, "NNC face at boundary. This should never happen!");
            }
        }
    }
    return ret;
}

const CpGrid::InterfaceMap& CpGrid::cellScatterGatherInterface() const
{
    return *cell_scatter_gather_interfaces_;
}

const CpGrid::InterfaceMap& CpGrid::pointScatterGatherInterface() const
{
    return *point_scatter_gather_interfaces_;
}

void CpGrid::switchToGlobalView()
{
    current_view_data_ = data_.back().get();
    current_data_ = &data_;
}

void CpGrid::switchToDistributedView()
{
    if (distributed_data_.empty())
        OPM_THROW(std::logic_error, "No distributed view available in grid");
    current_view_data_ = distributed_data_.back().get();
    current_data_ = &distributed_data_;
}

#if HAVE_MPI

const cpgrid::CpGridDataTraits::CommunicationType& CpGrid::cellCommunication() const
{
    return current_view_data_->cellCommunication();
}

cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet()
{
    return current_view_data_->cellIndexSet();
}

cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices()
{
    return current_view_data_->cellRemoteIndices();
}

const cpgrid::CpGridDataTraits::ParallelIndexSet& CpGrid::getCellIndexSet() const
{
    return current_view_data_->cellIndexSet();
}

const cpgrid::CpGridDataTraits::RemoteIndices& CpGrid::getCellRemoteIndices() const
{
    return current_view_data_->cellRemoteIndices();
}

#endif

#if HAVE_ECL_INPUT
std::vector<std::size_t>
CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                             Opm::EclipseState* ecl_state,
                             bool periodic_extension,
                             bool turn_normals, bool clip_z,
                             bool pinchActive)
{
    auto removed_cells = current_view_data_->processEclipseFormat(ecl_grid, ecl_state, periodic_extension,
                                                                  turn_normals, clip_z, pinchActive);
    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
                                         0);
    return removed_cells;
}

std::vector<std::size_t>
CpGrid::processEclipseFormat(const Opm::EclipseGrid* ecl_grid_ptr,
                             Opm::EclipseState* ecl_state,
                             bool periodic_extension, bool turn_normals, bool clip_z)
{
    return processEclipseFormat(ecl_grid_ptr, ecl_state, periodic_extension, turn_normals, clip_z,
                                !ecl_grid_ptr || ecl_grid_ptr->isPinchActive());
}

#endif

void CpGrid::processEclipseFormat(const grdecl& input_data,
                                  bool remove_ij_boundary, bool turn_normals)
{
    using NNCMap = std::set<std::pair<int, int>>;
    using NNCMaps = std::array<NNCMap, 2>;
    NNCMaps nnc;
    current_view_data_->processEclipseFormat(input_data,
#if HAVE_ECL_INPUT
                                             nullptr,
#endif
                                             nnc,
                                             remove_ij_boundary, turn_normals, false, 0.0);
    current_view_data_->ccobj_.broadcast(current_view_data_->logical_cartesian_size_.data(),
                                         current_view_data_->logical_cartesian_size_.size(),
                                         0);
}

template<int dim>
cpgrid::Entity<dim> createEntity(const CpGrid& grid,int index,bool orientation)
{
    return cpgrid::Entity<dim>(*grid.current_view_data_, index, orientation);
}
template cpgrid::Entity<0> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<3> createEntity(const CpGrid&, int, bool);
template cpgrid::Entity<1> createEntity(const CpGrid&, int, bool); // needed in distribution_test.cpp

bool CpGrid::mark(int refCount, const cpgrid::Entity<0>& element)
{
    // For serial run, mark elements also in the level they were born.
    if(currentData().size()>1) {
        // Mark element in its level
        currentData()[element.level()] -> mark(refCount, element.getEquivLevelElem());
    }
    // Mark element (also in the serial run case) in current_view_data_. Note that if scatterGrid has been invoked, then
    // current_view_data_ == distributed_data_[0].
    return current_view_data_-> mark(refCount, element);
}

int CpGrid::getMark(const cpgrid::Entity<0>& element) const
{
    return current_view_data_->getMark(element);
}

bool CpGrid::preAdapt()
{
    // Set the flags mighVanish for elements that have been marked for refinement/coarsening.

    // Check if elements in pre-adapt existing grids have been marked for refinment.
    // Serial run: currentData() = data_. Parallel run: currentData() = distributed_data_.
    bool isPreAdapted = false;
    for (const auto& preAdaptGrid : currentData()) {
        isPreAdapted = isPreAdapted || (preAdaptGrid -> preAdapt());
    }
    return isPreAdapted;
}

bool CpGrid::adapt()
{
    const std::vector<std::array<int,3>>& cells_per_dim_vec = {{2,2,2}}; // Arbitrary chosen values.
    std::vector<int> assignRefinedLevel(current_view_data_-> size(0));
    const auto& preAdaptMaxLevel = this ->maxLevel();

    int local_marked_elem_count = 0;
    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        const auto& element = cpgrid::Entity<0>(*current_view_data_, elemIdx, true);
        assignRefinedLevel[elemIdx] = (this->getMark(element) == 1) ? (preAdaptMaxLevel +1) : 0;
        if (this->getMark(element) == 1) {
            ++local_marked_elem_count;
        }
    }

    // Check if its a global refinement
    bool is_global_refine = false;
    std::vector<std::string> lgr_name_vec = { "LGR" + std::to_string(preAdaptMaxLevel +1) };
    // Rewrite if global refinement
#if HAVE_MPI
    auto global_marked_elem_count = comm().sum(local_marked_elem_count);
    auto global_cell_count_before_adapt = comm().sum(current_view_data_-> size(0)); // Recall overlap cells are also marked
    if (global_marked_elem_count == global_cell_count_before_adapt) {
        // GR stands for GLOBAL REFINEMET
        lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
        is_global_refine = true; // parallel
    }
#endif
    if ( (comm().size() == 0) && (local_marked_elem_count == current_view_data_-> size(0)) ) {
        // GR stands for GLOBAL REFINEMET
        lgr_name_vec = { "GR" + std::to_string(preAdaptMaxLevel +1) };
        is_global_refine = true; // sequential
    }
    if (is_global_refine) { // parallel or sequential
        const std::array<int,3>& endIJK = currentData().back()->logicalCartesianSize();
        return this->adapt(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec, is_global_refine, {{0,0,0}}, {endIJK});
    }
    return this-> adapt(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec);
}

bool CpGrid::adapt(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                   const std::vector<int>& assignRefinedLevel,
                   const std::vector<std::string>& lgr_name_vec,
                   bool isCARFIN,
                   const std::vector<std::array<int,3>>& startIJK_vec,
                   const std::vector<std::array<int,3>>& endIJK_vec)
{
    // To do: support coarsening.
    assert( static_cast<int>(assignRefinedLevel.size()) == current_view_data_->size(0));
    assert(cells_per_dim_vec.size() == lgr_name_vec.size());

    // Each marked element has its assigned level where its refined entities belong.
    const int& levels = cells_per_dim_vec.size();
    // Notice that "levels" represents also the total amount of new (after calling adapt) refined level grids.
    const int& preAdaptMaxLevel = this->maxLevel();
    // Copy corner history - needed to compute later ids, empty vector if the grid to be adapted is level 0 grid, or the grid has been distributed.
    const auto& preAdaptGrid_corner_history = (preAdaptMaxLevel>0) ? current_view_data_->corner_history_ : std::vector<std::array<int,2>>();

    auto& data = currentData(); // data pointed by current_view_data_ (data_ or distributed_data_[if loadBalance() has been invoked before adapt()]).

    // To determine if an LGR is not empty in a given process, we set
    // lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
    // at least one active cell, and to 0 otherwise.
    std::vector<int> lgr_with_at_least_one_active_cell(levels);
    detectActiveLgrs(startIJK_vec, endIJK_vec, lgr_with_at_least_one_active_cell);

    // To store/build refined level grids.
    std::vector<std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>> refined_data_vec(levels, data);
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> refined_grid_ptr_vec(levels);

    std::vector<Dune::cpgrid::DefaultGeometryPolicy> refined_geometries_vec(levels);
    std::vector<std::vector<std::array<int,8>>> refined_cell_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<0,1>> refined_cell_to_face_vec(levels);
    std::vector<Opm::SparseTable<int>> refined_face_to_point_vec(levels);
    std::vector<cpgrid::OrientedEntityTable<1,0>> refined_face_to_cell_vec(levels);
    std::vector<cpgrid::EntityVariable<enum face_tag,1>> refined_face_tags_vec(levels);
    std::vector<cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>> refined_face_normals_vec(levels);

    // Mutable containers for refined corners, faces, cells, face tags, and face normals.
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>> refined_corners_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>> refined_faces_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>> refined_cells_vec(levels);
    std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>> mutable_refined_face_tags_vec(levels);
    typedef Dune::FieldVector<double,3> PointType;
    std::vector<Dune::cpgrid::EntityVariableBase<PointType>> mutable_refined_face_normals_vec(levels);

    std::vector<std::vector<int>> refined_global_cell_vec(levels);


    // To store adapted grid
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& adaptedData = data;
#if HAVE_MPI
    auto adaptedGrid_ptr =
        std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, adaptedData);
#else
    // DUNE 2.7 is missing convertion to NO_COMM
    auto adaptedGrid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(adaptedData);
#endif
    auto& adaptedGrid = *adaptedGrid_ptr;
    Dune::cpgrid::DefaultGeometryPolicy&                         adapted_geometries = adaptedGrid.geometry_;
    std::vector<std::array<int,8>>&                              adapted_cell_to_point = adaptedGrid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>&                            adapted_cell_to_face = adaptedGrid.cell_to_face_;
    Opm::SparseTable<int>&                                       adapted_face_to_point = adaptedGrid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>&                            adapted_face_to_cell = adaptedGrid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>&                     adapted_face_tags = adaptedGrid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& adapted_face_normals = adaptedGrid.face_normals_;
    // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners =
        *(adapted_geometries.geomVector(std::integral_constant<int,3>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces =
        *(adapted_geometries.geomVector(std::integral_constant<int,1>()));
    Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells =
        *(adapted_geometries.geomVector(std::integral_constant<int,0>()));
    Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags = adapted_face_tags;
    typedef Dune::FieldVector<double,3> PointType;
    Dune::cpgrid::EntityVariableBase<PointType>& mutable_face_normals = adapted_face_normals;



    // Refine marked cells and provide marked-corner/face/cell - refined-corner/faces/cells relations.
    //
    // ------------------------ Marked elements parameters
    // -- markedElem_to_itsLgr :
    // Each marked element gets refined and we store this "auxiliary markedElementLGR", to later
    // build a unique level containing all the refined entities from all the marked elements.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData> > markedElem_to_itsLgr;
    markedElem_to_itsLgr.resize(current_view_data_->size(0));
    // -- markedElem_count: Total amount of marked elements to be refined. It will be used to print grid info.
    int markedElem_count = 0;
    // -- cornerInMarkedElemWithEquivRefinedCorner :
    // For each corner from level zero, we store the marked elements where the corner appears and its equivalent
    // refined corner in  each auxiliary marked-element-lgr. Example: corner with index 5 appears in marked
    // elements 0 and 1, with refined equivalent corner indices 8 and 2 respectively. Then,
    // cornerInMarkedElemWithEquivRefinedCorner[5] = {{0, 8}, {1, 2}}.
    // For corners not appearing in any marked element, empty vector.
    std::vector<std::vector<std::array<int,2>>> cornerInMarkedElemWithEquivRefinedCorner;
    cornerInMarkedElemWithEquivRefinedCorner.resize(current_view_data_->size(3) );
    // -- markedElemAndEquivRefinedCorner_to_corner :
    // To correctly build the level-refined and adapted-grid topology features, we need to keep track of the
    // corners that got replaced by equivalent refined corners, in each marked element where the corner appeared,
    // not only in its last appearance. The last appearance will be used to avoid repetition when storing.
    // Following the example above,
    // markedElemAndEquivRefinedCorner_to_corner[{0, 8}] = 5;
    // markedElemAndEquivRefinedCorner_to_corner[{1, 2}] = 5;
    std::map< std::array<int,2>, int > markedElemAndEquivRefinedCorn_to_corner;
    // -- faceInMarkedElemAndRefinedFaces :
    // For each face from level zero, we store the marked elements where the face appears (maximum 2 cells)
    // and its new-born refined faces from each auxiliary marked-element-lgr. Example: face with index 9
    // appears in marked elements 0 and 1. Then,
    // faceInMarkedElemAndRefinedFaces[9] = {{0, {refinedFace0_0, ..., refinedFaceN_0}},
    //                                       {1, {refinedFace0_1, ..., refinedFaceM_1}}}.
    // For faces not appearing in any marked element, empty vector.
    std::vector<std::vector<std::pair<int, std::vector<int>>>> faceInMarkedElemAndRefinedFaces;
    faceInMarkedElemAndRefinedFaces.resize(current_view_data_->face_to_cell_.size());
    // ------------------------ Refined cells parameters
    // --- Refined cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell;
    // Integer to count only REFINED cells (new-born refined cells from ANY marked element).
    std::vector<int> refined_cell_count_vec(levels, 0);
    // -- Parent-child relations --
    // Relation between the grids before adapt() and new refined cells (on the refined grid level - not in each individual lgr).
    std::vector<std::vector<std::tuple<int,std::vector<int>>>> preAdapt_parent_to_children_cells_vec(preAdaptMaxLevel +1);
    // ------------------------ Adapted cells parameters
    // --- Adapted cells and PreAdapt cells relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCell_to_adaptedCell;
    std::unordered_map<int,std::array<int,2>> adaptedCell_to_elemLgrAndElemLgrCell;
    // Integer to count adapted cells (mixed between cells from level0 (not involved in LGRs), and (new-born) refined cells).
    int cell_count = 0;
    // -- Some extra indices relations between preAdapt-grid and adapted-grid --
    // Relation between the grids before adapt() and leafview cell indices.
    std::vector<std::vector<int>> preAdapt_level_to_leaf_cells_vec(preAdaptMaxLevel +1);
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        // Resize with the corresponding amount of cells of the preAdapt level. Deafualt {-1, empty vector} when the cell has no children.
        if ( (*data[preAdaptLevel]).parent_to_children_cells_.empty()){
            preAdapt_parent_to_children_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), std::make_pair(-1, std::vector<int>{}));
        }
        else {
            preAdapt_parent_to_children_cells_vec[preAdaptLevel] =  (*data[preAdaptLevel]).parent_to_children_cells_;
        }
        // Resize with the corresponding amount of cell of the preAdapt level. Dafualt -1 when the cell vanished and does not appear on the leaf grid view.
        // In entry 'level cell index', we store 'leafview cell index', or -1 when the cell vanished.
        preAdapt_level_to_leaf_cells_vec[preAdaptLevel].resize(data[preAdaptLevel]->size(0), -1);
    }
    //
    refineAndProvideMarkedRefinedRelations( /* Marked elements parameters */
                                            markedElem_to_itsLgr,
                                            markedElem_count,
                                            cornerInMarkedElemWithEquivRefinedCorner,
                                            markedElemAndEquivRefinedCorn_to_corner,
                                            faceInMarkedElemAndRefinedFaces,
                                            /* Refined cells parameters */
                                            elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                            refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                            refined_cell_count_vec,
                                            assignRefinedLevel,
                                            preAdapt_parent_to_children_cells_vec,
                                            /* Adapted cells parameters */
                                            elemLgrAndElemLgrCell_to_adaptedCell,
                                            adaptedCell_to_elemLgrAndElemLgrCell,
                                            cell_count,
                                            preAdapt_level_to_leaf_cells_vec,
                                            /* Additional parameters */
                                            cells_per_dim_vec);

    // Update/define parent_to_children_cells_ and level_to_leaf_cells_ for all the existing level grids (level 0, 1, ..., preAdaptMaxLevel), before this call of adapt.
    for (int preAdaptLevel = 0; preAdaptLevel < preAdaptMaxLevel +1; ++preAdaptLevel) {
        (*data[preAdaptLevel]).parent_to_children_cells_ = preAdapt_parent_to_children_cells_vec[preAdaptLevel];
        (*data[preAdaptLevel]).level_to_leaf_cells_ =  preAdapt_level_to_leaf_cells_vec[preAdaptLevel];
    }

    // -- Child-parent relations --
    // refined_child_to_parent_cells_vec:   Refined child cells and their parents. Entry is {-1,-1} when cell has no father.
    //                                      Otherwise, {level parent cell, parent cell index}. Each entry represents a refined level.
    // refined_cell_to_idxInParentCell_vec: Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      Each entry represents a refined level.
    // adapted_child_to_parent_cells:       Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    // adapted_cell_to_idxInParentCell:     Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    //                                      When the cell has not been refined, -1.
    const auto& [refined_child_to_parent_cells_vec,
                 refined_cell_to_idxInParentCell_vec,
                 adapted_child_to_parent_cells,
                 adapted_cell_to_idxInParentCell] = defineChildToParentAndIdxInParentCell(refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                                          refined_cell_count_vec,adaptedCell_to_elemLgrAndElemLgrCell,
                                                                                          cell_count);

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    // refined_level_to_leaf_cells_vec:  Relation between the refined grid and leafview cell indices.
    // leaf_to_level_cells:              Relation between an adapted cell and its equivalent cell coming either from current_view_data_ or from the refined grid (level)
    const auto& [refined_level_to_leaf_cells_vec,
                 leaf_to_level_cells] = defineLevelToLeafAndLeafToLevelCells(elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                                             refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                                             refined_cell_count_vec,
                                                                             elemLgrAndElemLgrCell_to_adaptedCell,
                                                                             adaptedCell_to_elemLgrAndElemLgrCell,
                                                                             cell_count);

    // CORNERS
    // Stablish relationships between PreAdapt corners and refined or adapted ones ---
    //
    // --- Refined corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner;
    std::map<std::array<int,2>,std::array<int,2>> vanishedRefinedCorner_to_itsLastAppearance;
    // Integer to count only refined corners.
    std::vector<int> refined_corner_count_vec(levels, 0);
    identifyRefinedCornersPerLevel(/* Refined grid parameters */
                                   elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                   refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                   refined_corner_count_vec,
                                   vanishedRefinedCorner_to_itsLastAppearance,
                                   /* Additional parameters */
                                   markedElem_to_itsLgr,
                                   assignRefinedLevel,
                                   cornerInMarkedElemWithEquivRefinedCorner,
                                   faceInMarkedElemAndRefinedFaces,
                                   cells_per_dim_vec);

    // --- Adapted corners and PreAdapt corners relations ---
    std::map<std::array<int,2>,int>           elemLgrAndElemLgrCorner_to_adaptedCorner;
    std::unordered_map<int,std::array<int,2>> adaptedCorner_to_elemLgrAndElemLgrCorner;
    // Integer to count adapted corners (mixed between corners from current_view_data_ (not involved in LGRs), and (new-born) refined corners).
    int corner_count = 0;
    identifyLeafGridCorners(/* Adapted grid parameters */
                            elemLgrAndElemLgrCorner_to_adaptedCorner,
                            adaptedCorner_to_elemLgrAndElemLgrCorner,
                            corner_count,
                            /* Additional parameters */
                            markedElem_to_itsLgr,
                            assignRefinedLevel,
                            cornerInMarkedElemWithEquivRefinedCorner,
                            vanishedRefinedCorner_to_itsLastAppearance,
                            faceInMarkedElemAndRefinedFaces,
                            cells_per_dim_vec);

    // FACES
    // Stablish relationships between PreAdapt faces and refined or adapted ones ---
    // --- Refined faces and PreAdapt faces relations ---
    std::map<std::array<int,2>,std::array<int,2>> elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace;
    std::map<std::array<int,2>,std::array<int,2>> refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from level0 (not involved in LGRs), and (new-born) refined faces).
    std::vector<int> refined_face_count_vec(levels, 0);
    identifyRefinedFacesPerLevel( elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  refined_face_count_vec,
                                  markedElem_to_itsLgr,
                                  assignRefinedLevel,
                                  faceInMarkedElemAndRefinedFaces,
                                  cells_per_dim_vec);

    // --- Adapted faces and PreAdapt faces relations ---
    std::map< std::array<int,2>, int >           elemLgrAndElemLgrFace_to_adaptedFace;
    std::unordered_map< int, std::array<int,2> > adaptedFace_to_elemLgrAndElemLgrFace;
    // Integer to count adapted faces (mixed between faces from current_view_data_ (not involved in LGRs), and (new-born) refined faces).
    int face_count = 0;
    identifyLeafGridFaces(elemLgrAndElemLgrFace_to_adaptedFace,
                          adaptedFace_to_elemLgrAndElemLgrFace,
                          face_count,
                          markedElem_to_itsLgr,
                          assignRefinedLevel,
                          faceInMarkedElemAndRefinedFaces,
                          cells_per_dim_vec);

    setRefinedLevelGridsGeometries( /* Refined corner arguments */
                                    refined_corners_vec,
                                    refined_corner_count_vec,
                                    /* Refined face arguments */
                                    refined_faces_vec,
                                    mutable_refined_face_tags_vec,
                                    mutable_refined_face_normals_vec,
                                    refined_face_to_point_vec,
                                    refined_face_count_vec,
                                    /* Refined cell argumets */
                                    refined_cells_vec,
                                    refined_cell_to_point_vec,
                                    refined_global_cell_vec,
                                    refined_cell_count_vec,
                                    refined_cell_to_face_vec,
                                    refined_face_to_cell_vec,
                                    /* Auxiliary arguments */
                                    refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                    refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                    refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                    elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                    elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                    faceInMarkedElemAndRefinedFaces,
                                    refined_geometries_vec,
                                    vanishedRefinedCorner_to_itsLastAppearance,
                                    markedElem_to_itsLgr,
                                    assignRefinedLevel,
                                    preAdaptMaxLevel,
                                    markedElemAndEquivRefinedCorn_to_corner,
                                    cornerInMarkedElemWithEquivRefinedCorner,
                                    cells_per_dim_vec);

    updateLeafGridViewGeometries( /* Leaf grid View Corners arguments */
                                  adapted_corners,
                                  corner_count,
                                  /* Leaf grid View Faces arguments */
                                  adapted_faces,
                                  mutable_face_tags,
                                  mutable_face_normals,
                                  adapted_face_to_point,
                                  face_count,
                                  /* Leaf grid View Cells argumemts  */
                                  adapted_cells,
                                  adapted_cell_to_point,
                                  cell_count,
                                  adapted_cell_to_face,
                                  adapted_face_to_cell,
                                  /* Auxiliary arguments */
                                  adaptedCorner_to_elemLgrAndElemLgrCorner,
                                  adaptedFace_to_elemLgrAndElemLgrFace,
                                  adaptedCell_to_elemLgrAndElemLgrCell,
                                  elemLgrAndElemLgrFace_to_adaptedFace,
                                  faceInMarkedElemAndRefinedFaces,
                                  adapted_geometries,
                                  elemLgrAndElemLgrCorner_to_adaptedCorner,
                                  vanishedRefinedCorner_to_itsLastAppearance,
                                  markedElem_to_itsLgr,
                                  assignRefinedLevel,
                                  markedElemAndEquivRefinedCorn_to_corner,
                                  cornerInMarkedElemWithEquivRefinedCorner,
                                  cells_per_dim_vec,
                                  preAdaptMaxLevel);


    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
#if HAVE_MPI
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>((*(data[0])).ccobj_, refined_data_vec[level]);
#else
        // DUNE 2.7 is missing convertion to NO_COMM
        refined_grid_ptr_vec[level] = std::make_shared<Dune::cpgrid::CpGridData>(refined_data_vec[level]);
#endif
        // Store refined grid
        if ((level == 0) && (preAdaptMaxLevel>0)) { // Overwrite the leaf-grid-view with the first new-refined-level-grid
            data[preAdaptMaxLevel+1] = refined_grid_ptr_vec[level];
        }
        else {
            data.push_back(refined_grid_ptr_vec[level]);
        }

        Dune::cpgrid::DefaultGeometryPolicy&  refinedLevel_geometries = (*data[refinedLevelGridIdx]).geometry_;
        // Mutable containers for adapted corners, faces, cells, face tags, and face normals.
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& level_corners =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,3>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& level_faces =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,1>()));
        Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& level_cells =
            *(refinedLevel_geometries.geomVector(std::integral_constant<int,0>()));

        level_corners.swap(refined_corners_vec[level]);
        level_faces.swap(refined_faces_vec[level]);
        level_cells.swap(refined_cells_vec[level]);

        (*data[refinedLevelGridIdx]).cell_to_point_.swap(refined_cell_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).cell_to_face_.swap(refined_cell_to_face_vec[level]);

        (*data[refinedLevelGridIdx]).face_to_point_.swap(refined_face_to_point_vec[level]);
        (*data[refinedLevelGridIdx]).face_to_cell_.swap(refined_face_to_cell_vec[level]);

        cpgrid::EntityVariable<enum face_tag,1>& level_face_tags =   (*data[refinedLevelGridIdx]).face_tag_;
        Dune::cpgrid::EntityVariableBase<enum face_tag>& level_mutable_face_tags = level_face_tags;
        level_mutable_face_tags.swap(mutable_refined_face_tags_vec[level]);

        cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>&  level_face_normals =   (*data[refinedLevelGridIdx]).face_normals_;
        Dune::cpgrid::EntityVariableBase<PointType>& level_mutable_face_normals = level_face_normals;
        level_mutable_face_normals.swap(mutable_refined_face_normals_vec[level]);

        // Further Refined grid Attributes
        //
        // Populate some attributes of the level LGR
        (*data[refinedLevelGridIdx]).level_data_ptr_ = &(this -> currentData());
        (*data[refinedLevelGridIdx]).level_ = refinedLevelGridIdx;
        this -> lgr_names_[lgr_name_vec[level]] = refinedLevelGridIdx; // {"name_lgr", level}
        (*data[refinedLevelGridIdx]).child_to_parent_cells_ = refined_child_to_parent_cells_vec[level];
        (*data[refinedLevelGridIdx]).cell_to_idxInParentCell_ = refined_cell_to_idxInParentCell_vec[level];
        (*data[refinedLevelGridIdx]).level_to_leaf_cells_ =  refined_level_to_leaf_cells_vec[level];
        (*data[refinedLevelGridIdx]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[refinedLevelGridIdx]->size(0),
                                                                                     data[refinedLevelGridIdx]->size(3));
        // Determine the amount of cells per direction, per parent cell, of the corresponding LGR.
        (*data[refinedLevelGridIdx]).cells_per_dim_ = cells_per_dim_vec[level];
        // TO DO: This new code for refinement do not assume Cartesian Shape. How does logical_cartesian_size_ should be defined then?
        // When the refined level grid has been originated from a block of cells, then its logical Cartesian size
        // corresponds to the inner product between cells_per_dim_vec[level] and the dimension of the block (amount of cells in each direction).
        // In the case of a block of cells, e.g., when CARFIN keyword is used, we need the following:
        if (isCARFIN) {
            const auto& blockDim = (*data[0]).getPatchDim(startIJK_vec[level], endIJK_vec[level]);
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = { cells_per_dim_vec[level][0]*blockDim[0],
                                                                     cells_per_dim_vec[level][1]*blockDim[1],
                                                                     cells_per_dim_vec[level][2]*blockDim[2] };
        }
        else {
            (*data[refinedLevelGridIdx]).logical_cartesian_size_ = (*data[0]).logical_cartesian_size_;
            (*data[refinedLevelGridIdx]).global_cell_.swap(refined_global_cell_vec[level]);
        }
        // One alternative definition for logical_cartesian_size_ in the case where the marked elements for refinement do not form a block of cells,
        // therefore, are not associated with the keyword CARFIN, is to imagine that we put all the marked elements one next to the other, along
        // the x-axis. Then, the "imaginary" logical Cartesian size of the refined level grid would be
        // { (# marked elemnts)x cells_per_dim_vec[level][0], cells_per_dim_vec[level][1], cells_per_dim_vec[level][2]}.
        /** To do: how the definition of refined level grids logical_cartesian_size_ affects LookUpData class (and LookUpCartesianData)*/
    }

    // Store adapted grid
    data.push_back(adaptedGrid_ptr);

    // Further Adapted  grid Attributes
    (*data[levels + preAdaptMaxLevel +1]).child_to_parent_cells_ = adapted_child_to_parent_cells;
    (*data[levels + preAdaptMaxLevel +1]).cell_to_idxInParentCell_ = adapted_cell_to_idxInParentCell;
    (*data[levels + preAdaptMaxLevel +1]).leaf_to_level_cells_ =  leaf_to_level_cells;
    (*data[levels + preAdaptMaxLevel +1]).index_set_ = std::make_unique<cpgrid::IndexSet>(data[levels + preAdaptMaxLevel +1]->size(0),
                                                                                          data[levels + preAdaptMaxLevel +1]->size(3));
    (*data[levels + preAdaptMaxLevel +1]).logical_cartesian_size_ =  (*data[0]).logical_cartesian_size_;

    // Update the leaf grid view
    current_view_data_ = data.back().get();

    // When the refinement is determined by startIJK and endIJK values, the LGR has a (local) Cartesian size.
    // Therefore, each refined cell belonging to the LGR can be associated with a (local) IJK and its (local) Cartesian index.
    // If the LGR has NXxNYxNZ dimension, then the Cartesian indices take values
    // k*NN*NY + j*NX + i, where i<NX, j<Ny, k<NZ.
    // This index is stored in <refined-level-grid>.global_cell_[ refined cell index (~element.index()) ] =  k*NN*NY + j*NX + i.
    if (isCARFIN) {
        for (int level = 0; level < levels; ++level) {
            const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
            std::vector<int> global_cell_lgr(data[refinedLevelGridIdx]->size(0));
            computeGlobalCellLgr(refinedLevelGridIdx, startIJK_vec[level], global_cell_lgr);
            (*data[refinedLevelGridIdx]).global_cell_.swap(global_cell_lgr);
        }
    }

    std::vector<int> global_cell_leaf( data[levels + preAdaptMaxLevel +1]->size(0));
    computeGlobalCellLeafGridViewWithLgrs(global_cell_leaf);
    (*data[levels + preAdaptMaxLevel +1]).global_cell_.swap(global_cell_leaf);

    updateCornerHistoryLevels(cornerInMarkedElemWithEquivRefinedCorner,
                              elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                              adaptedCorner_to_elemLgrAndElemLgrCorner,
                              corner_count,
                              preAdaptGrid_corner_history,
                              preAdaptMaxLevel,
                              levels);

    this->global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*current_view_data_);
    for (int level = 0; level < levels; ++level) {
        const int refinedLevelGridIdx = level + preAdaptMaxLevel +1;
        this->global_id_set_ptr_->insertIdSet(*data[refinedLevelGridIdx]);
    }

     // Only for parallel runs
    // - Define global ids for refined level grids (level 1, 2, ..., maxLevel)
    // - Define GlobalIdMapping (cellMapping, faceMapping, pointMapping required per level)
    // - Define ParallelIndex for overlap cells and their neighbors
    if(comm().size()>1) {
#if HAVE_MPI
        // Prediction min cell and point global ids per process
        //
        // Predict how many new cells/points (born in refined level grids) need new globalIds, so we can assign unique
        // new ids ( and anticipate the maximum).
        // The grid is already refined according to the LGR specification
        // At this point, neither cell_index_set_ nor partition_type_indicator_ are populated.
        // Refined level grid cells:
        //    1. Inherit their partition type from their parent cell (i.e., element.father().partitionType()).
        //    2. Assign global ids only for interior cells.
        //    3. Communicate the already assigned cell global ids from interior to overlap refined cells.
        // Refined level grid points/vertices:
        // There are 4 partition types: interior, border, front, overlap. This classification requires that both
        // cell and face partition types are already defined, not available yet for refined level grids.
        //    1. Assign for all partition type points a 'candidate of global id' (unique in each process).
        //       Except the points that coincide with a point from level zero.
        //    2. Re-write the values for points that are corners of overlap refined cells, via communication.
        // Under the assumption of LGRs fully-interior, no communication is needed. In the general case, communication will be used
        // to populate overlap cell/point global ids on the refined level grids.
        /** Warning: due to the overlap layer size (equal to 1) cells that share corners or edges (not faces) with interior cells
            are not included/seen by the process. This, in some cases, ends up in duplicated point ids. */
        //
        int min_globalId_cell_in_proc = 0;
        int min_globalId_point_in_proc = 0;
        predictMinCellAndPointGlobalIdPerProcess(assignRefinedLevel,
                                                 cells_per_dim_vec,
                                                 lgr_with_at_least_one_active_cell,
                                                 min_globalId_cell_in_proc,
                                                 min_globalId_point_in_proc);

        // Only for level 1,2,.., maxLevel grids.
        // For each level, define the local-to-global maps for cells and points (for faces: empty).
        // 1) Assignment of new global ids is done only for owned cells and non-overlap points.
        // 2) For overlap cells and points: communicate. Not needed under the assumption of fully interior LGRs.
        std::vector<std::vector<int>> localToGlobal_cells_per_level(cells_per_dim_vec.size());
        std::vector<std::vector<int>> localToGlobal_points_per_level(cells_per_dim_vec.size());
        // Ignore faces - empty vectors.
        std::vector<std::vector<int>> localToGlobal_faces_per_level(cells_per_dim_vec.size());

        assignCellIdsAndCandidatePointIds(localToGlobal_cells_per_level,
                                          localToGlobal_points_per_level,
                                          min_globalId_cell_in_proc,
                                          min_globalId_point_in_proc,
                                          cells_per_dim_vec);


        const auto& parent_to_children = current_data_->front()->parent_to_children_cells_;
        ParentToChildrenCellGlobalIdHandle parentToChildrenGlobalId_handle(parent_to_children, localToGlobal_cells_per_level);
        currentData().front()->communicate(parentToChildrenGlobalId_handle,
                                           Dune::InteriorBorder_All_Interface,
                                           Dune::ForwardCommunication );

        // After assigning global IDs to points in refined-level grids, a single point may have
        // a "unique" global ID in each local leaf grid view for every process to which it belongs.
        // To ensure true uniqueness, since global IDs must be distinct across the global leaf view
        // and consistent across each refined-level grid, we will rewrite the entries in
        // localToGlobal_points_per_level.
        //
        // This correction is done using cell_to_point_ across all refined cells through
        // communication: gathering the 8 corner points of each interior cell and scattering the
        // 8 corner points of overlapping cells, for all child cells of a parent cell in level zero grid.
        //
        // Design decision: Why we communicate via level zero grid instead of in each refined level grid.
        // The reason is that how children are stored (the ordering) in parent_to_children_cells_
        // is always the same, accross all processes.
        // Even though the ordering of the corners in cell_to_point_ is the same accross all processes,
        // this may not be enough to correctly overwrite the "winner" point global ids for refined cells.
        //
        /** Current approach avoids duplicated point ids when
         // 1. the LGR is distributed in P_{i_0}, ..., P_{i_n}, with n+1 < comm().size(),
         // AND
         // 2. there is no coarse cell seen by a process P with P != P_{i_j}, j = 0, ..., n.
         // Otherwise, there will be duplicated point ids.
         //
         // Reason: neighboring cells that only share corners (not faces) are NOT considered in the
         // overlap layer of the process.*/
        selectWinnerPointIds(localToGlobal_points_per_level,
                             parent_to_children,
                             cells_per_dim_vec);

        for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
            // Global id set for each (refined) level grid.
            if(lgr_with_at_least_one_active_cell[level-1]>0) { // Check if LGR is active in currect process.
                (*current_data_)[level]->global_id_set_->swap(localToGlobal_cells_per_level[level-1],
                                                              localToGlobal_faces_per_level[level-1],
                                                              localToGlobal_points_per_level[level-1]);
            }
        }

        for (std::size_t level = 1; level < cells_per_dim_vec.size()+1; ++level) {
            
            populateCellIndexSetRefinedGrid(level);
            // Compute the partition type for cell
            currentData()[level]->computeCellPartitionType();
            // Compute the partition type for point
            currentData()[level]->computePointPartitionType();
            // Now we can compute the communication interface.
            currentData()[level]->computeCommunicationInterfaces(currentData()[level]->size(3));
        }

        populateLeafGlobalIdSet();

        this->global_id_set_ptr_ = std::make_shared<cpgrid::GlobalIdSet>(*(current_data_->back()));
        for (std::size_t level = 0; level < cells_per_dim_vec.size()+1; ++level) {
            this->global_id_set_ptr_->insertIdSet(*(*current_data_)[level]);
        }

        populateCellIndexSetLeafGridView();
        
        // Compute the partition type for cell
        (*current_data_).back()->computeCellPartitionType();

        // Compute the partition type for point
        (*current_data_).back()->computePointPartitionType();

        // Now we can compute the communication interface.
        current_data_->back()->computeCommunicationInterfaces(current_data_->back()->size(3));
        assert(static_cast<std::size_t>(current_data_->back()->cellIndexSet().size()) == static_cast<std::size_t>(current_data_->back()->size(0)) );
#endif
    } // end-if-comm().size()>1


    // Print total amount of cells on the adapted grid
    Opm::OpmLog::info(std::to_string(markedElem_count) + " elements have been marked (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(levels)  + " (new) refined level grid(s) (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(cell_count)  + " total cells on the leaf grid view (in " + std::to_string(comm().rank()) + " rank).\n");

    return preAdapt();
}

void CpGrid::postAdapt()
{
    // - Resize with the new amount of cells on the leaf grid view
    // - Set marks equal to zero (representing 'doing nothing')
    current_view_data_ -> postAdapt();
}

void CpGrid::addLgrsUpdateLeafView(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec,
                                   const std::vector<std::string>& lgr_name_vec)
{
    // For parallel run, level zero grid is stored in distributed_data_[0]. If CpGrid::scatterGrid has been invoked,
    // then current_view_data_ == distributed_data_[0].
    // For serial run, level zero grid is stored in data_[0]. In this case, current_view_data_ == data_[0].
    // Note: currentData() returns data_ (if grid is not distributed) or distributed_data_ otherwise.

    // Check startIJK_vec and endIJK_vec have same size, and "startIJK[patch][coordinate] < endIJK[patch][coordinate]"
    current_view_data_->validStartEndIJKs(startIJK_vec, endIJK_vec);

    // Sizes of provided vectors (number of subivisions per cells and lgrs name) should coincide.
    bool matchingSizeHasFailed = false;
    if ( (cells_per_dim_vec.size() != startIJK_vec.size())  || (lgr_name_vec.size() != startIJK_vec.size())) {
        matchingSizeHasFailed = true;
    }
    matchingSizeHasFailed = comm().max(matchingSizeHasFailed);
    if (matchingSizeHasFailed) {
        OPM_THROW(std::invalid_argument, "Sizes of provided vectors with subdivisions per cell and LGR names need to match.");
    }

    // Compatibility of number of subdivisions of neighboring LGRs: Check shared faces on boundaries of LGRs.
    //                                                              Not optimal since the code below does not take into account
    //                                                              active/inactive cells, instead, relies on "ijk-computations".
    //                                                              TO DO: improve/remove.
    // To check "Compatibility of numbers of subdivisions of neighboring LGRs".
    // The method compatibleSubdivision returns a bool. We convert it into an int since MPI within DUNE does not support bool directly.
    int compatibleSubdivisions = current_view_data_->compatibleSubdivisions(cells_per_dim_vec, startIJK_vec, endIJK_vec);
    compatibleSubdivisions = comm().min(compatibleSubdivisions); // 0 when at least one process returns false (un-compatible subdivisions).
    if(!compatibleSubdivisions) {
        if (comm().rank()==0){
            OPM_THROW(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
        }
        else{
            OPM_THROW_NOLOG(std::logic_error, "Subdivisions of neighboring LGRs sharing at least one face do not coincide. Not suppported yet.");
        }
    }

    // Non neighboring connections: Currently, adding LGRs whose cells have NNCs is not supported yet.
    // The method nonNNCsSelectedCellsLGR returns a bool. We convert it into an int since MPI within DUNE does not support bool directly.
    int nonNNCs = this->nonNNCsSelectedCellsLGR(startIJK_vec, endIJK_vec);
    // To check "Non-NNCs (non non neighboring connections)" for all processes.
    nonNNCs = comm().min(nonNNCs); // 0 when at least one process returns false (there are NNCs on a selected cell for refinement).
    if(!nonNNCs) {
        OPM_THROW(std::logic_error, "NNC face on a cell containing LGR is not supported yet.");
    }

    // Determine the assigned level for the refinement of each marked cell
    std::vector<int> assignRefinedLevel(current_view_data_->size(0));
    // To determine if an LGR is not empty in a given process, we set
    // lgr_with_at_least_one_active_cell[in that level] to 1 if it contains
    // at least one active cell, and to 0 otherwise.
    std::vector<int> lgr_with_at_least_one_active_cell(startIJK_vec.size());
    markElemAssignLevelDetectActiveLgrs(startIJK_vec,
                                        endIJK_vec,
                                        assignRefinedLevel,
                                        lgr_with_at_least_one_active_cell);

    int non_empty_lgrs = 0;
    for (std::size_t level = 0; level < startIJK_vec.size(); ++level) {
        // Do not throw if all cells of an LGR are inactive in a parallel run (The process might not 'see' those cells.)
        if (lgr_with_at_least_one_active_cell[level] == 0) {
            Opm::OpmLog::warning("LGR" + std::to_string(level+1) + " contains only inactive cells (in " + std::to_string(comm().rank()) + " rank).\n");
        }
        else {
            ++non_empty_lgrs;
        }
    }

    // Notice that in a parallel run, non_empty_lgrs represents the local active lgrs, i.e. the lgrs containing active cells which also belong
    // to the current process.
    auto globalActiveLgrs = comm().sum(non_empty_lgrs);
    if(globalActiveLgrs == 0) {
        Opm::OpmLog::warning("All the LGRs contain only inactive cells.\n");
    }

    preAdapt();
    adapt(cells_per_dim_vec, assignRefinedLevel, lgr_name_vec, true, startIJK_vec, endIJK_vec);
    postAdapt();

    // Print total refined level grids and total cells on the leaf grid view
    Opm::OpmLog::info(std::to_string(non_empty_lgrs) + " (new) refined level grid(s) (in " + std::to_string(comm().rank()) + " rank).\n");
    Opm::OpmLog::info(std::to_string(current_view_data_->size(0)) + " total cells on the leaf grid view (in " + std::to_string(comm().rank()) + " rank).\n");
}


const std::map<std::string,int>& CpGrid::getLgrNameToLevel() const{
    return lgr_names_;
}

std::array<double,3> CpGrid::getEclCentroid(const int& elemIdx) const
{
    return this-> current_view_data_ -> computeEclCentroid(elemIdx);
}

std::array<double,3> CpGrid::getEclCentroid(const cpgrid::Entity<0>& elem) const
{
    return this-> getEclCentroid(elem.index());
}

void CpGrid::refineAndProvideMarkedRefinedRelations( /* Marked elements parameters */
                                                     std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                                     int& markedElem_count,
                                                     std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                                     std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                                     std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                                     /* Refined cells parameters */
                                                     std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                                     std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                                     std::vector<int>& refined_cell_count_vec,
                                                     const std::vector<int>& assignRefinedLevel,
                                                     std::vector<std::vector<std::tuple<int,std::vector<int>>>>& preAdapt_parent_to_children_cells_vec,
                                                     /* Adapted cells parameters */
                                                     std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                                     std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                                     int& cell_count,
                                                     std::vector<std::vector<int>>& preAdapt_level_to_leaf_cells_vec,
                                                     /* Additional parameters */
                                                     const std::vector<std::array<int,3>>& cells_per_dim_vec) const 
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // Each marked element for refinement (mark equal to 1), will be refined individuality, creating its own Lgr. The element index will
    // be also used to identify its lgr. Even though, in the end, all the refined entities will belong to a unique level grid.
    // For this reason, we associate "-1" with those elements that are not involved in any refinement and will appear
    // as "coarse" cells in the leaf-grid-view (adapted-grid).

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();

    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        const auto& element = Dune::cpgrid::Entity<0>(*current_view_data_, elemIdx, true);
        // When the element is marked with 0 ("doing nothing"), it will appear in the adapted grid with same geometrical features (center, volume).
        if (getMark(element) ==  0) {
            elemLgrAndElemLgrCell_to_adaptedCell[{-1, elemIdx}] = cell_count;
            adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {-1, elemIdx};
            cell_count +=1;
            preAdapt_level_to_leaf_cells_vec[element.level()][element.getEquivLevelElem().index()] = cell_count;
        }
        
        // When the element is marked for refinement, we also mark its corners and faces
        // since they will get replaced by refined ones.
        if (getMark(element) ==  1) {
            markedElem_count +=1;
            const auto& markedElemLevel = assignRefinedLevel[elemIdx];
            assert(markedElemLevel > preAdaptMaxLevel);
            // Shift the markedElemRefinedLevel to access data containers
            const auto& shiftedLevel = markedElemLevel - preAdaptMaxLevel-1;
            // Build auxiliary LGR for the refinement of this element
            const auto& [elemLgr_ptr,
                         parentCorners_to_equivalentRefinedCorners,
                         parentFace_to_itsRefinedFaces,
                         parentCell_to_itsRefinedCells,
                         refinedFace_to_itsParentFace,
                         refinedCell_to_itsParentCell]
                = current_view_data_->refineSingleCell(cells_per_dim_vec[shiftedLevel], elemIdx);
            markedElem_to_itsLgr[ elemIdx ] = elemLgr_ptr;

            const auto& childrenCount = cells_per_dim_vec[shiftedLevel][0]*cells_per_dim_vec[shiftedLevel][1]*cells_per_dim_vec[shiftedLevel][2];
            std::vector<int> refinedChildrenList(childrenCount);

            for (int refinedCell = 0; refinedCell < childrenCount; ++refinedCell) {

                elemLgrAndElemLgrCell_to_adaptedCell[{elemIdx, refinedCell}] = cell_count;
                adaptedCell_to_elemLgrAndElemLgrCell[cell_count] = {elemIdx, refinedCell};
                cell_count +=1;

                elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell[{elemIdx, refinedCell}] = { markedElemLevel, refined_cell_count_vec[shiftedLevel]};
                refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell[{markedElemLevel, refined_cell_count_vec[shiftedLevel]}] = {elemIdx, refinedCell};
                refinedChildrenList[refinedCell] = refined_cell_count_vec[shiftedLevel];
                refined_cell_count_vec[shiftedLevel] +=1;

            }
            
            preAdapt_parent_to_children_cells_vec[element.level()][element.getEquivLevelElem().index()] = std::make_pair( markedElemLevel, refinedChildrenList);
            for (const auto& [markedCorner, lgrEquivCorner] : parentCorners_to_equivalentRefinedCorners) {
                cornerInMarkedElemWithEquivRefinedCorner[markedCorner].push_back({elemIdx, lgrEquivCorner});
                markedElemAndEquivRefinedCorn_to_corner[ {elemIdx, lgrEquivCorner}] = markedCorner;
            }
            for (const auto& [markedFace, itsRefinedFaces] : parentFace_to_itsRefinedFaces) {
                faceInMarkedElemAndRefinedFaces[markedFace].push_back({elemIdx, itsRefinedFaces});
            }
        } // end-if-elemMark==1
    } // end-elem-for-loop
}

std::tuple<std::vector<std::vector<std::array<int,2>>>, std::vector<std::vector<int>>, std::vector<std::array<int,2>>, std::vector<int>>
CpGrid::defineChildToParentAndIdxInParentCell(const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                              const std::vector<int>& refined_cell_count_vec,
                                              const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                              const int& cell_count) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // ------------------------ Refined grid parameters
    // Refined child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    // Each entry represents a refined level.
    std::vector<std::vector<std::array<int,2>>> refined_child_to_parent_cells_vec(refined_cell_count_vec.size());
    // Each entry represents a refined level. Each refined child cell has a unique index in its parent cell, to be used to build geometryInFather().
    std::vector<std::vector<int>> refined_cell_to_idxInParentCell_vec(refined_cell_count_vec.size());
    // ------------------------ Adapted grid parameters
    // Adapted child cells and their parents. Entry is {-1,-1} when cell has no father. Otherwise, {level parent cell, parent cell index}
    std::vector<std::array<int,2>> adapted_child_to_parent_cells;
    std::vector<int> adapted_cell_to_idxInParentCell;

    adapted_child_to_parent_cells.resize(cell_count, std::array<int,2>{-1,-1}); //  Entry is {-1,-1} when cell has no father, otherwise, {level parent cell, parent cell index}.
    adapted_cell_to_idxInParentCell.resize(cell_count, -1);

    // Rewrite only the entries of adapted cells that have a parent cell
    for (int cell = 0; cell < cell_count; ++cell) {
        // For the element with index "cell" in the adapted grid (or updated-leaf-grid-view),
        // - if "cell" is a new born refined cell (i.e. born in the current adapt-call), we get the parent cell from the preAdapt-leaf-grid-view ("current_view_data_").
        //   In this case, "preAdapt_parent_or_elem" represents "preAdapt_parent".
        // - if "cell" is either a coarse cell or a refined cell that was born in a preAdapt-refined-level-grid, we get the equivalent cell in the
        //   preAdapt-leaf-grid-view ("current_view_data_"). In this case, "preAdapt_parent_or_elem" represents "preAdapt_elem".
        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        // Get the element of either a parent cell of a new born refined cell with index "cell" or an equivalent cell, in the preAdapt-leaf-grid-view ("current_view_data_").
        const auto& preAdapt_parent_or_elem = Dune::cpgrid::Entity<0>(*current_view_data_, ((elemLgr != -1) ? elemLgr : elemLgrCell), true);
        if (elemLgr != -1) { // "cell" is a new born refined cell
            adapted_child_to_parent_cells[cell] = {preAdapt_parent_or_elem.level(), preAdapt_parent_or_elem.getEquivLevelElem().index()};
            adapted_cell_to_idxInParentCell[cell] = elemLgrCell;
        }
        else {// "cell" is either a coarse cell or a refined cell that was born in a preAdapt-refined-level-grid
            // Only populate the entries of refined cells that were born in preAdapt-refined-level-grids.
            if (preAdapt_parent_or_elem.hasFather()) {
                adapted_child_to_parent_cells[cell] =  {preAdapt_parent_or_elem.father().level(), preAdapt_parent_or_elem.father().index() };
                adapted_cell_to_idxInParentCell[cell] = currentData()[preAdapt_parent_or_elem.level()]->
                    cell_to_idxInParentCell_[preAdapt_parent_or_elem.getEquivLevelElem().index()];
            }
        }
    }

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {
        refined_child_to_parent_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_cell_to_idxInParentCell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        const int& level =  shiftedLevel + preAdaptMaxLevel +1;
        // Every refined cell has a parent cell
        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            const auto& [elemLgr, elemLgrCell] = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({level, cell});
            // (elemLgr == -1) means that the adapted cell is equivalent to a cell from the starting grid (current_view_data_)
            assert(elemLgr != -1);
            // Search for the level where the parent cell was born, and its index in that level grid.
            // Notice that elemLgr is a cell index in current_view_data_ (the starting grid where elements got marked for (further) refinement).
            const auto& element = Dune::cpgrid::Entity<0>(*current_view_data_, elemLgr, true);  // elemLgr == parent cell index in starting grid.
            refined_child_to_parent_cells_vec[shiftedLevel][cell] = {element.level(), element.getEquivLevelElem().index()};
            refined_cell_to_idxInParentCell_vec[shiftedLevel][cell] = elemLgrCell;
        }
    }

    return std::make_tuple<std::vector<std::vector<std::array<int,2>>>, std::vector<std::vector<int>>,
                           std::vector<std::array<int,2>>, std::vector<int>>(std::move(refined_child_to_parent_cells_vec),
                                                                             std::move(refined_cell_to_idxInParentCell_vec),
                                                                             std::move(adapted_child_to_parent_cells),
                                                                             std::move(adapted_cell_to_idxInParentCell));

}

std::pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>
CpGrid::defineLevelToLeafAndLeafToLevelCells(const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                             const std::vector<int>& refined_cell_count_vec,
                                             const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCell_to_adaptedCell,
                                             const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                             const int& cell_count) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // -- Refined to Adapted cells and Adapted-cells to {level where the cell was born, cell index on that level} --
    // Relation between the refined grid and leafview cell indices.
    std::vector<std::vector<int>> refined_level_to_leaf_cells_vec(refined_cell_count_vec.size());
    // Relation between an adapted cell and its equivalent cell coming either from current_view_data_ or from the refined grid (level)
    std::vector<std::array<int,2>> leaf_to_level_cells;
    leaf_to_level_cells.resize(cell_count);

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();

    // -- Adapted to {level, cell index in that level}  --
    for (int cell = 0; cell < cell_count; ++cell) {
        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        // elemLgr == -1 means that this adapted cell is equivalent to a cell from the starting grid. So we need to find out the level where that equivalent
        // cell was born, as well as its cell index in that level.
        if (elemLgr == -1) {
            const auto& element = Dune::cpgrid::Entity<0>(*current_view_data_, elemLgrCell, true);
            leaf_to_level_cells[cell] = { element.level(), element.getEquivLevelElem().index()};
        }
        else {
            leaf_to_level_cells[cell] = elemLgrAndElemLgrCell_to_refinedLevelAndRefinedCell.at({elemLgr, elemLgrCell});
        }
    }
    // -- Refined to adapted cells --
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel){
        refined_level_to_leaf_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {
            refined_level_to_leaf_cells_vec[shiftedLevel][cell] =
                elemLgrAndElemLgrCell_to_adaptedCell.at(refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1, cell}));
        }
    }
    return std::make_pair<std::vector<std::vector<int>>, std::vector<std::array<int,2>>>(std::move(refined_level_to_leaf_cells_vec), std::move(leaf_to_level_cells));
}

void CpGrid::identifyRefinedCornersPerLevel(std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                            std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                            std::vector<int>& refined_corner_count_vec,
                                            std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                            const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                            const std::vector<int>& assignRefinedLevel,
                                            const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                            const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                            const std::vector<std::array<int,3>>& cells_per_dim_vec) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();

    // Step 1. Replace the corners from the preAdapt grid involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the preAdapt corner
    //         in the LGRs.
    for (int corner = 0; corner < current_view_data_->size(3); ++corner) {
        if (!cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) {
            // corner involved in refinement, so we search for it in one LGR (the last one where it appears)
            // Get the lgr corner that replaces the marked corner from level zero.
            // Note: Recall that lgr coincides with the marked element index from the preAdapt grid that got refined.
            //       Since the container is a map, the lgr and the lgr corner index correspond to the last
            //       appearance of the marked corner (from the starting grid - where elements got marked).
            const auto& [lastAppearanceLgr, lastAppearanceLgrCorner] = cornerInMarkedElemWithEquivRefinedCorner[corner].back();

            const auto& lastAppearanceLgrLevel = assignRefinedLevel[lastAppearanceLgr];
            assert(lastAppearanceLgrLevel>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = lastAppearanceLgrLevel - preAdaptMaxLevel -1;

            elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{lastAppearanceLgr, lastAppearanceLgrCorner}] =
                {lastAppearanceLgrLevel, refined_corner_count_vec[shiftedLevel]};
            refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{lastAppearanceLgrLevel, refined_corner_count_vec[shiftedLevel]}] =
                {lastAppearanceLgr, lastAppearanceLgrCorner};
            refined_corner_count_vec[shiftedLevel] +=1;

            if (cornerInMarkedElemWithEquivRefinedCorner[corner].size()>1) {
                for (const auto& [elemLgr, elemLgrCorner] : cornerInMarkedElemWithEquivRefinedCorner[corner]) {
                    const auto& elemLgrLevel = assignRefinedLevel[elemLgr];
                    if (elemLgrLevel != lastAppearanceLgrLevel) {
                        const auto& shiftedElemLgrLevel = elemLgrLevel - preAdaptMaxLevel -1;
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemLgr, elemLgrCorner}] = {elemLgrLevel, refined_corner_count_vec[shiftedElemLgrLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{elemLgrLevel, refined_corner_count_vec[shiftedElemLgrLevel]}] = {elemLgr, elemLgrCorner};
                        refined_corner_count_vec[shiftedElemLgrLevel] +=1;
                    }
                }
            }
        }
    } // end corner-forloop

    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        if (markedElem_to_itsLgr.at(elemIdx)!= nullptr) {
            const auto& level = assignRefinedLevel[elemIdx];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int corner = 0; corner < markedElem_to_itsLgr.at(elemIdx) ->size(3); ++corner) {
                // Discard marked corners. Store (new born) refined corners

                // INTERIOR
                if (isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], corner)) { // It's a refined interior corner, so we store it.
                    // In this case, the corner is a new born refined corner that does not
                    // coincide with any corner from the GLOBAL grid (level 0). Therefore,
                    // it has to be stored.
                    elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level, refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                    refined_corner_count_vec[shiftedLevel] +=1;
                }

                // LYING ON EDGES
                //
                // Refined corners lying on edges - Refined edge has a 'coarse' parent edge (line between 2 corners of the parent cell)
                // To avoid repetition, we distinguish the case where the refined corner lies on an edge of its parent cell.
                // We detect the two coarse faces involved (Notice that the extremes of the parent cell have been stored previously).
                // When the marked faces appears only once, we store the corner now. Otherwise, we store the refined corner on its
                // last appearence associated with one of these parent faces, taking also into account the elemLgr. For example, when
                // the refined corners lie on an edge connecting I_FACE false and K_FACE true of the parent cell, let's say iFaceIdx,
                // kFaceIdx, with each of those faces appearing twice (maximum) :
                // iFaceIdx appearing in current "elem" and elemLgr1
                // kFaceIdx appearing in current "elem" and elemLgr2
                // Then, we take the max(elemLgr1, elemLgr2) and store the refined corner only if this maximum equals elem.
                if (newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    const auto& markedFacesTouchingEdge = getParentFacesAssocWithNewRefinedCornLyingOnEdge(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    const auto& [markedFace1, markedFace2] = markedFacesTouchingEdge;

                    int lastAppearanceMarkedFace1 = faceInMarkedElemAndRefinedFaces[markedFace1].back().first; // elemLgr1
                    int lastAppearanceMarkedFace2 = faceInMarkedElemAndRefinedFaces[markedFace2].back().first; // elemLgr2

                    int maxLastAppearance = std::max(lastAppearanceMarkedFace1, lastAppearanceMarkedFace2);
                    int faceAtMaxLastAppearance = (maxLastAppearance == lastAppearanceMarkedFace1) ? markedFace1 : markedFace2;

                    // Save the relationship between the vanished refined corner and its last appearance
                    const auto& maxLastAppearanceLevel = assignRefinedLevel[maxLastAppearance];
                    const auto& maxLastAppearanceLevelShifted = assignRefinedLevel[maxLastAppearance] - preAdaptMaxLevel -1;

                    bool atLeastOneFaceAppearsTwice = (faceInMarkedElemAndRefinedFaces[markedFace1].size()>1) ||
                        (faceInMarkedElemAndRefinedFaces[markedFace2].size()>1);
                    if (atLeastOneFaceAppearsTwice && (maxLastAppearance != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel],
                                                                                                  corner, elemIdx, faceAtMaxLastAppearance,
                                                                                                  cells_per_dim_vec[maxLastAppearanceLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {maxLastAppearance, neighboringLgrCornerIdx};
                        // Notice that, when we use these container to locate vanished corners, we might need a while-loop,
                        // since {elem, corner} leads to {lastMaxAppearance, neighboringLgrCornerIdx}, which can also vanish.
                        // So we need something like:
                        // if (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({elem, corner}) == 0)
                        //    int updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][0];
                        //    int updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][1];
                        //     while (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 0)
                        //        int tempElemLgr =  updateElemLgr;
                        //        int tempElemLgrCorner =  updateElemLgrCorner;
                        //        updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                        //        updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        // Then, use the lastest update to search for the corner in teh refined/adapted grid (which would be the one that
                        // gives elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 1).
                    }
                    if ((maxLastAppearance == elemIdx) || (level!= maxLastAppearanceLevel)) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level,refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }

                // LYING ON BOUNDARY LGR - NOT ON AN EDGE - NOT COINCIDING WITH A MARKED CORNER
                //
                // If the refined corner lies on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this corner now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                if ( isRefinedNewBornCornerOnLgrBoundary(cells_per_dim_vec[shiftedLevel], corner) &&
                     !newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    // Get the index of the marked face where the refined corner was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedCornerLiesOn(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    // check how many times marked face appearn
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;

                    const auto& lastLgrLevel = assignRefinedLevel[lastLgrWhereMarkedFaceAppeared ];
                    const auto& lastLgrLevelShifted = lastLgrLevel - preAdaptMaxLevel -1;
                    // Save the relationship between the vanished refined corner and its last appearance
                    if ((faceInMarkedElemAndRefinedFaces[markedFace].size()>1) && (lastLgrWhereMarkedFaceAppeared != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner,
                                                                                                  cells_per_dim_vec[lastLgrLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {lastLgrWhereMarkedFaceAppeared, neighboringLgrCornerIdx};
                    }

                    if ((lastLgrWhereMarkedFaceAppeared == elemIdx) || (lastLgrLevel != level)) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner[{elemIdx, corner}] = {level, refined_corner_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner[{level, refined_corner_count_vec[shiftedLevel]}] = {elemIdx, corner};
                        refined_corner_count_vec[shiftedLevel] +=1;
                    }
                }
            } // end-corner-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void CpGrid::identifyRefinedFacesPerLevel(std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                          std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                          std::vector<int>& refined_face_count_vec,
                                          const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                          const std::vector<int>& assignRefinedLevel,
                                          const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                          const std::vector<std::array<int,3>>& cells_per_dim_vec) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();

    // Step 1. Add the LGR faces, for each LGR
    for (int elem = 0; elem < current_view_data_->size(0); ++elem) {
        if (markedElem_to_itsLgr[elem]!=nullptr)  {
            const auto& level = assignRefinedLevel[elem];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int face = 0; face < markedElem_to_itsLgr[elem] ->face_to_cell_.size(); ++face) {
                // Discard marked faces. Store (new born) refined faces
                bool isNewRefinedFaceOnLgrBoundary = isRefinedFaceOnLgrBoundary(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem]);
                if (!isNewRefinedFaceOnLgrBoundary) { // It's a refined interior face, so we store it
                    // In this case, the face is a new born refined face that does not
                    // have any "parent face" from the GLOBAL grid (level 0).
                    elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elem, face}] = {level, refined_face_count_vec[shiftedLevel]};
                    refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level, refined_face_count_vec[shiftedLevel]}] = {elem, face};
                    refined_face_count_vec[shiftedLevel] +=1;
                }
                // If the refined face lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this face now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                else {
                    // Get the index of the marked face where the refined corner was born.
                    int markedFace = getParentFaceWhereNewRefinedFaceLiesOn(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem], elem);
                    assert(!faceInMarkedElemAndRefinedFaces[markedFace].empty());
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    if (lastLgrWhereMarkedFaceAppeared == elem) {
                        // Store the refined face in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{elem, face}] = {level, refined_face_count_vec[shiftedLevel]};
                        refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{level, refined_face_count_vec[shiftedLevel]}] = {elem, face};
                        refined_face_count_vec[shiftedLevel] +=1;
                    }
                    if(faceInMarkedElemAndRefinedFaces[markedFace].size()>1) { // maximum size is 2
                        const auto& firstMarkedElem = faceInMarkedElemAndRefinedFaces[markedFace][0].first;
                        const auto& firstMarkedElemLevel = assignRefinedLevel[firstMarkedElem];
                        if (firstMarkedElemLevel != level) {
                            const auto& shiftedFirstMarkedElemLevel = firstMarkedElemLevel - preAdaptMaxLevel -1;
                            elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace[{firstMarkedElem, face}] = {firstMarkedElemLevel, refined_face_count_vec[shiftedFirstMarkedElemLevel]};
                            refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace[{firstMarkedElemLevel, refined_face_count_vec[shiftedFirstMarkedElemLevel]}] = {firstMarkedElem, face};
                            refined_face_count_vec[shiftedFirstMarkedElemLevel] +=1;
                        }
                    }
                }
            } // end-face-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void CpGrid::identifyLeafGridCorners(std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                     std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                     int& corner_count,
                                     const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                     const std::vector<int>& assignRefinedLevel,
                                     const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                     std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                     const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                     const std::vector<std::array<int,3>>& cells_per_dim_vec) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // Step 1. Select/store the corners from the starting grid not involved in any (new) LGR.
    //         Replace the corners from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int corner = 0; corner < current_view_data_->size(3); ++corner) {
        if (cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the grid before adapt is called with the value -1.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{-1, corner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {-1, corner};
            corner_count +=1;
        }
        else { // corner involved in refinement, so we search it in one LGR (the last one where it appears)
            assert(!cornerInMarkedElemWithEquivRefinedCorner[corner].empty());
            // Get the lgr corner that replaces the marked corner from level zero.
            // Note: Recall that lgr coincides with the marked element index from level 0 that got refined.
            //       Since the container is a map, the lgr and the lgr corner index correspond to the last
            //       appearance of the marked corner (from level 0).
            const auto& [lastAppearanceLgr, lastAppearanceLgrCorner] = cornerInMarkedElemWithEquivRefinedCorner[corner].back();
            // Build the relationships between adapted corner and level corner, for future search due topology aspects.
            elemLgrAndElemLgrCorner_to_adaptedCorner[{lastAppearanceLgr, lastAppearanceLgrCorner}] = corner_count;
            adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {lastAppearanceLgr, lastAppearanceLgrCorner};
            corner_count +=1;
        }
    } // end corner-forloop

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel();

    for (int elemIdx = 0; elemIdx < current_view_data_->size(0); ++elemIdx) {
        if (markedElem_to_itsLgr[elemIdx]!= nullptr) {
            const auto& level = assignRefinedLevel[elemIdx];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int corner = 0; corner < markedElem_to_itsLgr[elemIdx] ->size(3); ++corner) {
                // Discard marked corners. Store (new born) refined corners

                // INTERIOR
                if (isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], corner)) { // It's a refined interior corner, so we store it.
                    // In this case, the corner is a new born refined corner that does not
                    // coincide with any corner from the GLOBAL grid (level 0). Therefore,
                    // it has to be stored.
                    elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                    adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                    corner_count += 1;
                }

                // LYING ON EDGES
                //
                // Refined corners lying on edges - Refined edge has a 'coarse' parent edge (line between 2 corners of the parent cell)
                // To avoid repetition, we distinguish the case where the refined corner lies on an edge of its parent cell.
                // We detect the two coarse faces involved (Notice that the extremes of the parent cell have been stored previously).
                // When the marked faces appeares only once, we store the corner now. Otherwise, we store the refined corner on its
                // last apparance assocaited to one of these parent faces, taking also into account the elemLgr. For example, when
                // the refined corners lies on an edge connecting I_FACE false and K_FACE true of the parent cell, let's say iFaceIdx,
                // kFaceIdx, with each of those faces appearing twice (maximum) :
                // iFaceIdx appearing in current "elem" and elemLgr1
                // kFaceIdx appearing in current "elem" and elemLgr2
                // Then, we take the max(elemLgr1, elemLgr2) and store the refined corner only if this maximum equals elem.
                if (newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    const auto& markedFacesTouchingEdge = getParentFacesAssocWithNewRefinedCornLyingOnEdge(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    const auto& [markedFace1, markedFace2] = markedFacesTouchingEdge;

                    int lastAppearanceMarkedFace1 = faceInMarkedElemAndRefinedFaces[markedFace1].back().first; // elemLgr1
                    int lastAppearanceMarkedFace2 = faceInMarkedElemAndRefinedFaces[markedFace2].back().first; // elemLgr2

                    int maxLastAppearance = std::max(lastAppearanceMarkedFace1, lastAppearanceMarkedFace2); // max(elemLgr1, elemLgr2)
                    int faceAtMaxLastAppearance = (maxLastAppearance == lastAppearanceMarkedFace1) ? markedFace1 : markedFace2;

                    int maxLastAppearanceLevel = assignRefinedLevel[maxLastAppearance];
                    int maxLastAppearanceLevelShifted = assignRefinedLevel[maxLastAppearance] - preAdaptMaxLevel -1;


                    // Save the relationship between the vanished refined corner and its last appearance
                    bool atLeastOneFaceAppearsTwice = (faceInMarkedElemAndRefinedFaces[markedFace1].size()>1) ||
                        (faceInMarkedElemAndRefinedFaces[markedFace2].size()>1);
                    if ((atLeastOneFaceAppearsTwice && (maxLastAppearance != elemIdx)) || (maxLastAppearanceLevel != level)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel],
                                                                                                  corner, elemIdx, faceAtMaxLastAppearance,
                                                                                                  cells_per_dim_vec[maxLastAppearanceLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {maxLastAppearance, neighboringLgrCornerIdx};

                        // Notice that, when we use this container to locate vanished corners, we might need a while-loop,
                        // since {elem, corner} leads to {lastMaxAppearance, neighboringLgrCornerIdx}, which can also vanish.
                        // So we need something like:
                        // if (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({elem, corner}) == 0)
                        //    int updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][0];
                        //    int updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{elem, corner}][1];
                        //     while (elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 0)
                        //        int tempElemLgr =  updateElemLgr;
                        //        int tempElemLgrCorner =  updateElemLgrCorner;
                        //        updateElemLgr =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][0];
                        //        updateElemLgrCorner =  vanishedRefinedCorner_to_itsLastAppearance[{ tempElemLgr ,  tempElemLgrCorner}][1];
                        // Then, use the lastest update to search for the corner in teh refined/adapted grid (which would be the one that
                        // gives elemLgrAndElemLgrCorner_to_adapted/refinedCorner.count({updateElemLgr, updateElemLgCorner}) == 1).
                    }
                    if (maxLastAppearance == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;
                    }
                }

                // LYING ON BOUNDARY LGR - NOT ON AN EDGE - NOT COINCIDING WITH A MARKED CORNER
                //
                // If the refined corner lies on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this corner now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                if ( isRefinedNewBornCornerOnLgrBoundary(cells_per_dim_vec[shiftedLevel], corner) &&
                     !newRefinedCornerLiesOnEdge(cells_per_dim_vec[shiftedLevel], corner)) {
                    // Get the index of the marked face where the refined corner was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedCornerLiesOn(cells_per_dim_vec[shiftedLevel], corner, elemIdx);
                    // check how many times marked face appearn
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;

                    const auto& lastLgrLevel = assignRefinedLevel[lastLgrWhereMarkedFaceAppeared ];
                    const auto& lastLgrLevelShifted = lastLgrLevel - preAdaptMaxLevel -1;

                    // Save the relationship between the vanished refined corner and its last appearance
                    if ((faceInMarkedElemAndRefinedFaces[markedFace].size()>1) && (lastLgrWhereMarkedFaceAppeared != elemIdx)) {
                        const auto& neighboringLgrCornerIdx = replaceLgr1CornerIdxByLgr2CornerIdx(cells_per_dim_vec[shiftedLevel], corner,
                                                                                                  cells_per_dim_vec[lastLgrLevelShifted]);
                        vanishedRefinedCorner_to_itsLastAppearance[{elemIdx, corner}] = {lastLgrWhereMarkedFaceAppeared, neighboringLgrCornerIdx};
                    }

                    if (lastLgrWhereMarkedFaceAppeared == elemIdx) {
                        // Store the refined corner in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrCorner_to_adaptedCorner[{elemIdx, corner}] = corner_count;
                        adaptedCorner_to_elemLgrAndElemLgrCorner[corner_count] = {elemIdx, corner};
                        corner_count += 1;
                    }
                }
            } // end-corner-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
}

void CpGrid::identifyLeafGridFaces(std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                   std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                   int& face_count,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.

    // Max level before calling adapt.
    const int& preAdaptMaxLevel = this->maxLevel(); 
    
    // Step 1. Add the LGR faces, for each LGR
    for (int elem = 0; elem < current_view_data_->size(0); ++elem) {
        if (markedElem_to_itsLgr[elem]!=nullptr)  {
            const auto& level = assignRefinedLevel[elem];
            assert(level>0);
            // To access containers with refined level grid information
            const auto& shiftedLevel = level - preAdaptMaxLevel -1;
            for (int face = 0; face < markedElem_to_itsLgr[elem] ->face_to_cell_.size(); ++face) {
                // Discard marked faces. Store (new born) refined faces
                bool isNewRefinedFaceOnLgrBoundary = isRefinedFaceOnLgrBoundary(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem]);
                if (!isNewRefinedFaceOnLgrBoundary) { // It's a refined interior face, so we store it
                    //  if (isRefinedFaceInInteriorLgr(cells_per_dim, face, markedElem_to_itsLgr[elem])) {
                    // In this case, the face is a new born refined face that does not
                    // have any "parent face" from the GLOBAL grid (level 0).
                    elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                    adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                    face_count += 1;
                }
                // If the refined face lays on the boundary of the LGR, e.i., it was born on one of the faces
                // of the marked element that got refined, then, we have two cases:
                // - the marked face appears only in one marked element -> then, we store this face now.
                // - the marked face appears twice (maximum times) in two marked elements -> we store it later.
                else {
                    // Get the index of the marked face where the refined corner was born.
                    int markedFace = getParentFaceWhereNewRefinedFaceLiesOn(cells_per_dim_vec[shiftedLevel], face, markedElem_to_itsLgr[elem], elem);
                    assert(!faceInMarkedElemAndRefinedFaces[markedFace].empty());
                    // Get the last LGR (marked element) where the marked face appeared.
                    int lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    if (lastLgrWhereMarkedFaceAppeared == elem) {
                        // Store the refined face in its last appearence - to avoid repetition.
                        elemLgrAndElemLgrFace_to_adaptedFace[{elem, face}] = face_count;
                        adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {elem, face};
                        face_count += 1;
                    }
                }
            } // end-face-for-loop
        } // end-if-nullptr
    } // end-elem-for-loop
    // Step 2. Select/store the faces from the starting grid (where cells got marked) not involved in any LGR.
    //         Replace the faces from level zero involved in LGR by the equivalent ones, born in LGRs.
    //         In this case, we avoid repetition considering the last appearance of the level zero corner
    //         in the LGRs.
    for (int face = 0; face < current_view_data_->face_to_cell_.size(); ++face) {
        if (faceInMarkedElemAndRefinedFaces[face].empty()) { // save it
            // Note: Since we are associating each LGR with its parent cell index, and this index can take
            //       the value 0, we will represent the current_view_data_ with the value -1
            elemLgrAndElemLgrFace_to_adaptedFace[{-1, face}] = face_count;
            adaptedFace_to_elemLgrAndElemLgrFace[face_count] = {-1, face};
            face_count +=1;
        }
    } // end face-forloop
}


void CpGrid::populateLeafGridCorners(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners,
                                     const int& corner_count,
                                     const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                     const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    adapted_corners.resize(corner_count);
    for (int corner = 0; corner < corner_count; ++corner) {
        const auto& [elemLgr, elemLgrCorner] = adaptedCorner_to_elemLgrAndElemLgrCorner.at(corner);
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, we represent the current_view_data_ with the value -1
        adapted_corners[corner] = ((elemLgr == -1) ? *current_view_data_ : *markedElem_to_itsLgr[elemLgr]).geometry_.geomVector(std::integral_constant<int,3>())-> get(elemLgrCorner);
    }
}

void CpGrid::populateRefinedCorners(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                    const std::vector<int>& refined_corner_count_vec,
                                    const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                    const int& preAdaptMaxLevel,
                                    const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner) const
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_corner_count_vec.size(); ++shiftedLevel) {
        refined_corners_vec[shiftedLevel].resize(refined_corner_count_vec[shiftedLevel]);
        for (int corner = 0; corner < refined_corner_count_vec[shiftedLevel]; ++corner) {
            const auto& [elemLgr, elemLgrCorner] = refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1,corner});
            refined_corners_vec[shiftedLevel][corner] =  markedElem_to_itsLgr[elemLgr] -> geometry_.geomVector(std::integral_constant<int,3>()) -> get(elemLgrCorner);
        }
    }
}

void CpGrid::populateLeafGridFaces(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces,
                                   Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                   Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                   Opm::SparseTable<int>& adapted_face_to_point,
                                   const int& face_count,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                   const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                   const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const int& preAdaptMaxLevel) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    adapted_faces.resize(face_count);
    mutable_face_tags.resize(face_count);
    mutable_face_normals.resize(face_count);

    // Auxiliary integer to count all the points in leaf_face_to_point.
    int num_points = 0;
    // Auxiliary vector to store face_to_point with non consecutive indices.
    std::vector<std::vector<int>> aux_face_to_point;
    aux_face_to_point.resize(face_count);
    for (int face = 0; face < face_count; ++face) {

        const auto& [elemLgr, elemLgrFace] = adaptedFace_to_elemLgrAndElemLgrFace.at(face);
        // Note: Since we are associating each LGR with its parent cell index, and this index can take
        //       the value 0, with the value -1 we refer to current_view_data_ (grid where the elements have been marked).
        const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);
        const auto& grid_or_elemLgr_data = (elemLgr == -1) ? *current_view_data_ : *markedElem_to_itsLgr[elemLgr];

        // Get the face geometry.
        adapted_faces[face] = (*(grid_or_elemLgr_data.geometry_.geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
        // Get the face tag.
        mutable_face_tags[face] = grid_or_elemLgr_data.face_tag_[elemLgrFaceEntity];
        // Get the face normal.
        mutable_face_normals[face] = grid_or_elemLgr_data.face_normals_[elemLgrFaceEntity];
        // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
        const auto& preAdapt_face_to_point = grid_or_elemLgr_data.face_to_point_[elemLgrFace];
        // Add the amount of points to the count num_points.
        num_points += preAdapt_face_to_point.size();

        for (std::size_t corn = 0; corn < preAdapt_face_to_point.size(); ++corn) {
            std::size_t adaptedCorn = 0; // It'll get rewritten.
            const auto& elemLgrCorn = preAdapt_face_to_point[corn];
            // Corner is stored in adapted_corners
            if (auto candidate = elemLgrAndElemLgrCorner_to_adaptedCorner.find({elemLgr, elemLgrCorn}); candidate !=  elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                adaptedCorn = candidate->second;
            }
            else{
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].empty())) {
                    lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[elemLgrCorn].back();
                }
                if (elemLgr>-1) { // It represents the refinement of the element with index "elemIdx == elemLgr"
                    if (auto corner_candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, elemLgrCorn}); corner_candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner =  cornerInMarkedElemWithEquivRefinedCorner[corner_candidate->second].back();
                    }
                    else {
                        const auto& shiftedLevel = assignRefinedLevel[elemLgr] - preAdaptMaxLevel -1; // Assigned level > preAdapt maxLevel
                        bool isNewRefinedCornInInteriorLgr = isRefinedCornerInInteriorLgr(cells_per_dim_vec[shiftedLevel], elemLgrCorn);
                        assert(!isNewRefinedCornInInteriorLgr);
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count(lastAppearanceLgr_lgrEquivCorner) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index, using find instead of count.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, elemLgrCorn});
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.find( lastAppearanceLgr_lgrEquivCorner ) == elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                            auto& tempLgr_lgrCorner = lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                } // end-if(elemLgr>-1)
                adaptedCorn =   elemLgrAndElemLgrCorner_to_adaptedCorner.at(lastAppearanceLgr_lgrEquivCorner);
            }
            aux_face_to_point[face].push_back(adaptedCorn);
        }
    } // end-adapted-faces
    // Adapted/Leaf-Grid-View face_to_point.
    adapted_face_to_point.reserve(face_count, num_points);
    for (int face = 0; face < face_count; ++face) {
        adapted_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
    }
}

void CpGrid::populateRefinedFaces(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                  std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refined_face_normals_vec,
                                  std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                  const std::vector<int>& refined_face_count_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                  const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const int& preAdaptMaxLevel,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner) const
{
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_face_count_vec.size(); ++shiftedLevel) {

        // Store the refined faces
        refined_faces_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_tags_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);
        mutable_refined_face_normals_vec[shiftedLevel].resize(refined_face_count_vec[shiftedLevel]);

        // Auxiliary integer to count all the points in refined_face_to_point.
        int refined_num_points = 0;
        // Auxiliary vector to store refined_face_to_point with non consecutive indices.
        std::vector<std::vector<int>> aux_refined_face_to_point;
        aux_refined_face_to_point.resize(refined_face_count_vec[shiftedLevel]);
        for (int face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {

            const auto& [elemLgr, elemLgrFace] = refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1,face});
            const auto& elemLgrFaceEntity =  Dune::cpgrid::EntityRep<1>(elemLgrFace, true);

            // Get the face geometry.
            refined_faces_vec[shiftedLevel][face] = (*(markedElem_to_itsLgr.at(elemLgr)->geometry_.geomVector(std::integral_constant<int,1>())))[elemLgrFaceEntity];
            // Get the face tag.
            mutable_refined_face_tags_vec[shiftedLevel][face] = markedElem_to_itsLgr.at(elemLgr)->face_tag_[elemLgrFaceEntity];
            // Get the face normal.
            mutable_refined_face_normals_vec[shiftedLevel][face] = markedElem_to_itsLgr.at(elemLgr)->face_normals_[elemLgrFaceEntity];
            // Get face_to_point_ before adapting - we need to replace the level corners by the adapted ones.
            Opm::SparseTable<int>::mutable_row_type preAdapt_face_to_point = markedElem_to_itsLgr.at(elemLgr)->face_to_point_[elemLgrFace];
            // Add the amount of points to the count num_points.
            refined_num_points += preAdapt_face_to_point.size();

            // Face_to_point
            for (std::size_t corn = 0; corn < preAdapt_face_to_point.size(); ++corn) {
                const auto& elemLgrCorn = preAdapt_face_to_point[corn];
                std::size_t refinedCorn = 0; // It'll be rewritten.
                // Corner is stored in adapted_corners
                if (auto candidate = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find({elemLgr, elemLgrCorn});
                    candidate != elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end())  {
                    refinedCorn = candidate->second[1];
                }
                else{
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                    if(auto corner_candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, elemLgrCorn});
                       corner_candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[corner_candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_refinedCorner.count(lastAppearanceLgr_lgrCorner) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index, using find instead of count.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, elemLgrCorn});
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find(lastAppearanceLgr_lgrEquivCorner) ==
                               elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                            auto& tempLgr_lgrCorner  = lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                    refinedCorn =  elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(lastAppearanceLgr_lgrEquivCorner)[1];
                }
                aux_refined_face_to_point[face].push_back(refinedCorn);
            }
        }
        // Refined face_to_point.
        refined_face_to_point_vec[shiftedLevel].reserve(refined_face_count_vec[shiftedLevel], refined_num_points);
        for (int face = 0; face < refined_face_count_vec[shiftedLevel]; ++face) {
            refined_face_to_point_vec[shiftedLevel].appendRow(aux_refined_face_to_point[face].begin(), aux_refined_face_to_point[face].end());
        }
    }
}

void CpGrid::populateLeafGridCells(Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells,
                                   std::vector<std::array<int,8>>& adapted_cell_to_point,
                                   const int& cell_count,
                                   cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                   cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                   const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                   const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                   const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                                   const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                   const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                   const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                   const std::vector<int>& assignRefinedLevel,
                                   const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                   const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                   const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                   const int& preAdaptMaxLevel) const
{
    // If the (level zero) grid has been distributed, then the preAdaptGrid is data_[0]. Otherwise, preApaptGrid is current_view_data_.
    
    // --- Adapted cells ---
    // Store the adapted cells. Main difficulty: to lookup correctly the indices of the corners and faces of each cell.
    adapted_cells.resize(cell_count);
    adapted_cell_to_point.resize(cell_count);
    for (int cell = 0; cell < cell_count; ++cell) {

        const auto& [elemLgr, elemLgrCell] = adaptedCell_to_elemLgrAndElemLgrCell.at(cell);
        const auto& elemLgrCellEntity =  Dune::cpgrid::EntityRep<0>(elemLgrCell, true);

        // Auxiliary cell_to_face
        std::vector<cpgrid::EntityRep<1>> aux_cell_to_face;

        const auto& allCorners = adapted_geometries.geomVector(std::integral_constant<int,3>());
        const auto& grid_or_elemLgr_data = (elemLgr == -1) ? *current_view_data_ : *markedElem_to_itsLgr.at(elemLgr);

        // Get the cell geometry.
        const auto& cellGeom = (*(grid_or_elemLgr_data.geometry_.geomVector(std::integral_constant<int,0>()) ) )[elemLgrCellEntity];
        // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
        const auto& preAdapt_cell_to_point  =  grid_or_elemLgr_data.cell_to_point_[elemLgrCell];
        // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
        const auto& preAdapt_cell_to_face =  grid_or_elemLgr_data.cell_to_face_[elemLgrCellEntity];

        // Cell to point.
        for (int corn = 0; corn < 8; ++corn) {
            int adaptedCorn = 0; // It'll get rewritten.
            const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
            auto adapted_candidate =  elemLgrAndElemLgrCorner_to_adaptedCorner.find({elemLgr, preAdaptCorn});
            if ( adapted_candidate == elemLgrAndElemLgrCorner_to_adaptedCorner.end() ) {
                // Corner might have vanished - Search its equivalent lgr-corner in that case -
                // last lgr where the corner appears -
                std::array<int,2> lastAppearanceLgr_lgrEquivCorner = {0, 0}; // It'll get rewritten.
                if ((elemLgr ==-1) && (!cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].empty())) {
                    lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[preAdaptCorn].back();
                }
                if (elemLgr > -1) {
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    if(auto candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, preAdaptCorn}); candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrEquivCorner = cornerInMarkedElemWithEquivRefinedCorner[candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count( lastAppearanceLgr_lgrEquivCorner ) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgr_lgrEquivCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, preAdaptCorn});
                        while (elemLgrAndElemLgrCorner_to_adaptedCorner.find(lastAppearanceLgr_lgrEquivCorner) == elemLgrAndElemLgrCorner_to_adaptedCorner.end()) {
                            auto& tempLgr_lgrCorner =  lastAppearanceLgr_lgrEquivCorner;
                            lastAppearanceLgr_lgrEquivCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                }
                adaptedCorn =  elemLgrAndElemLgrCorner_to_adaptedCorner.at(lastAppearanceLgr_lgrEquivCorner);
            }
            // Corner is stored in adapted_corners
            else {
                assert(adapted_candidate != elemLgrAndElemLgrCorner_to_adaptedCorner.end());
                adaptedCorn = adapted_candidate->second;
            }
            adapted_cell_to_point[cell][corn] = adaptedCorn;
        } // end-cell_to_point

        // Cell to face.
        for (const auto& face : preAdapt_cell_to_face) {
            const auto& preAdaptFace = face.index();
            int adaptedFace = 0; // It will be rewritten.
            // Face is stored in adapted_faces
            if (auto candidate = elemLgrAndElemLgrFace_to_adaptedFace.find({elemLgr, preAdaptFace}); candidate != elemLgrAndElemLgrFace_to_adaptedFace.end()) {
                adaptedFace = candidate->second;
                aux_cell_to_face.push_back({adaptedFace, face.orientation()});
            }
            else{
                // Face might have vanished - Search its refined lgr-children faces in that case -
                // last lgr where the face appears
                if (elemLgr ==-1) { // Coarse face got replaced by its children - from the last appearance of the marked face.
                    assert(!faceInMarkedElemAndRefinedFaces[preAdaptFace].empty());
                    const auto& [lastAppearanceLgr, lastAppearanceLgrFaces] = faceInMarkedElemAndRefinedFaces[preAdaptFace].back();
                    for (const auto& refinedFace : lastAppearanceLgrFaces) {
                        adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace.at({lastAppearanceLgr, refinedFace});
                        aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                    }
                }
                if (elemLgr>-1) { // Refined face vanished and its equivalent refined face from a neighboring lgr got stored.
                    // Get shifted level
                    const auto& shiftedLevel = assignRefinedLevel[elemLgr] - preAdaptMaxLevel -1; // Assigned level > preAdapt maxLevel
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLiesOn(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    const auto& lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel], preAdaptFace, markedElem_to_itsLgr[elemLgr],
                                                                                             cells_per_dim_vec[assignRefinedLevel[lastLgrWhereMarkedFaceAppeared] - preAdaptMaxLevel -1]);
                    adaptedFace = elemLgrAndElemLgrFace_to_adaptedFace.at({lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace});
                    aux_cell_to_face.push_back({adaptedFace, face.orientation()});
                }
            }
        } // end-cell_to_face
        // Adapted/Leaf-grid-view cell to face.
        adapted_cell_to_face.appendRow(aux_cell_to_face.begin(), aux_cell_to_face.end());

        // Create a pointer to the first element of "adapted_cell_to_point" (required as the fourth argement to construct a Geometry<3,3> type object).
        int* indices_storage_ptr = adapted_cell_to_point[cell].data();
        adapted_cells[cell] = cpgrid::Geometry<3,3>(cellGeom.center(), cellGeom.volume(), allCorners, indices_storage_ptr);
    } // adapted_cells

    // Adapted/Leaf-grid-view face to cell.
    adapted_cell_to_face.makeInverseRelation(adapted_face_to_cell);
}


void CpGrid::populateRefinedCells(std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                  std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                  std::vector<std::vector<int>>& refined_global_cell_vec,
                                  const std::vector<int>& refined_cell_count_vec,
                                  std::vector<cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                  std::vector<cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                  const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                  const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                                  const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                  const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                  const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                  const std::vector<int>& assignRefinedLevel,
                                  const int& preAdaptMaxLevel,
                                  const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                  const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                  const std::vector<std::array<int,3>>&  cells_per_dim_vec) const
{   
    // --- Refined cells ---
    for (std::size_t shiftedLevel = 0; shiftedLevel < refined_cell_count_vec.size(); ++shiftedLevel) {

        refined_cells_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_cell_to_point_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);
        refined_global_cell_vec[shiftedLevel].resize(refined_cell_count_vec[shiftedLevel]);

        const auto& allLevelCorners = refined_geometries_vec[shiftedLevel].geomVector(std::integral_constant<int,3>());

        for (int cell = 0; cell < refined_cell_count_vec[shiftedLevel]; ++cell) {

            const auto& [elemLgr, elemLgrCell] = refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell.at({static_cast<int>(shiftedLevel) + preAdaptMaxLevel +1, cell});
            assert(elemLgr >-1);
            const auto& elemLgrCellEntity = Dune::cpgrid::EntityRep<0>(elemLgrCell, true);
            // Auxiliary cell_to_face
            std::vector<cpgrid::EntityRep<1>> aux_refined_cell_to_face;

            // To supoort the simulation of a mixed grid (coarse and refined cells), global_cell_ values of refined-level-grids
            // is a vector of consecutive indices from 0 to total cells on the level-grid. In case we want to modify this, take into
            // account that LookUpData and LookUpCartesianData will be affected, therefore, code in opm-simulators related to
            // searching for field features when the grid has LGRs will be affected as well.
            refined_global_cell_vec[shiftedLevel][cell] = (preAdaptMaxLevel>0) ? current_view_data_->global_cell_[elemLgr] : cell;
            // Get pre-adapt corners of the cell that will be replaced with leaf view ones.
            const auto& preAdapt_cell_to_point = markedElem_to_itsLgr.at(elemLgr)->cell_to_point_[elemLgrCell];
            // Get pre-adapt faces of the cell that will be replaced with leaf view ones.
            const auto& preAdapt_cell_to_face = markedElem_to_itsLgr.at(elemLgr)->cell_to_face_[elemLgrCellEntity];

            // Cell to point.
            for (int corn = 0; corn < 8; ++corn) {
                int refinedCorn;
                const auto& preAdaptCorn = preAdapt_cell_to_point[corn];
                auto refined_candidate = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find({elemLgr, preAdaptCorn});
                if (refined_candidate == elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                    // Corner might have vanished - Search its equivalent lgr-corner in that case -
                    // last lgr where the corner appears -
                    std::array<int,2> lastAppearanceLgr_lgrCorner = {0, 0}; // It'll be rewritten
                    if(auto candidate = markedElemAndEquivRefinedCorn_to_corner.find({elemLgr, preAdaptCorn}); candidate != markedElemAndEquivRefinedCorn_to_corner.end()) {
                        lastAppearanceLgr_lgrCorner = cornerInMarkedElemWithEquivRefinedCorner[candidate->second].back();
                    }
                    else {
                        // To locate vanished corners, we need a while-loop, since {elemLgr, elemLgrcorner} leads to
                        // {neighboringElemLgr, neighboringElemLgrCornerIdx}, which might have also vanished.
                        // Then, use the lastest appearance of the current corner, meaning, the first (and unique one - by construction) that
                        // gives elemLgrAndElemLgrCorner_to_adaptedCorner.count( lastAppearanceLgr_lgrCorner ) == 1).
                        // This corner lies on the area occupied by a coarse face that got refined and belonged to two marked elements.
                        // Get the index of this corner with respect to the greatest marked element index.
                        lastAppearanceLgr_lgrCorner = vanishedRefinedCorner_to_itsLastAppearance.at({elemLgr, preAdaptCorn});
                        while (elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.find(lastAppearanceLgr_lgrCorner) == elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.end()) {
                            auto& tempLgr_lgrCorner =   lastAppearanceLgr_lgrCorner;
                            lastAppearanceLgr_lgrCorner =  vanishedRefinedCorner_to_itsLastAppearance.at(tempLgr_lgrCorner);
                        }
                    }
                    refinedCorn = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(lastAppearanceLgr_lgrCorner)[1];
                }
                // Corner is stored in adapted_corners
                else {
                    refinedCorn = refined_candidate->second[1];
                }
                refined_cell_to_point_vec[shiftedLevel][cell][corn] = refinedCorn;
            } // end-cell_to_point

            // Cell to face.
            for (const auto& face : preAdapt_cell_to_face) {
                const auto& preAdaptFace = face.index();
                int refinedFace = 0; // It'll be rewritten.
                // Face might have vanished - Search its refined lgr-children faces in that case -
                // last lgr where the face appears
                auto face_candidate = elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.find({elemLgr, preAdaptFace});
                if (face_candidate == elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.end()) {
                    // Get the index of the marked face where the refined face was born.
                    const auto& markedFace = getParentFaceWhereNewRefinedFaceLiesOn(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                    markedElem_to_itsLgr[elemLgr], elemLgr);
                    // Get the last LGR (marked element) where the marked face appeared.
                    const int& lastLgrWhereMarkedFaceAppeared = faceInMarkedElemAndRefinedFaces[markedFace].back().first;
                    const auto& lastAppearanceLgrEquivFace = replaceLgr1FaceIdxByLgr2FaceIdx(cells_per_dim_vec[shiftedLevel], preAdaptFace,
                                                                                             markedElem_to_itsLgr[elemLgr],
                                                                                             cells_per_dim_vec[assignRefinedLevel[lastLgrWhereMarkedFaceAppeared] - preAdaptMaxLevel -1]);
                    refinedFace = elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace.at({lastLgrWhereMarkedFaceAppeared, lastAppearanceLgrEquivFace})[1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
                // Face is stored in adapted_faces
                else {
                    refinedFace = face_candidate->second[1];
                    aux_refined_cell_to_face.push_back({refinedFace, face.orientation()});
                }
            } // end-cell_to_face
            // Refined cell to face.
            refined_cell_to_face_vec[shiftedLevel].appendRow(aux_refined_cell_to_face.begin(), aux_refined_cell_to_face.end());
            // Get the cell geometry.
            const auto& elemLgrGeom =  (*( markedElem_to_itsLgr.at(elemLgr)->geometry_.geomVector(std::integral_constant<int,0>())))[elemLgrCellEntity];

            // Create a pointer to the first element of "refined_cell_to_point" (required as the fourth argement to construct a Geometry<3,3> type object).
            int* indices_storage_ptr = refined_cell_to_point_vec[shiftedLevel][cell].data();
            refined_cells_vec[shiftedLevel][cell] = cpgrid::Geometry<3,3>(elemLgrGeom.center(), elemLgrGeom.volume(), allLevelCorners, indices_storage_ptr);
        } // refined_cells
        // Refined face to cell.
        refined_cell_to_face_vec[shiftedLevel].makeInverseRelation(refined_face_to_cell_vec[shiftedLevel]);
    } // end-shiftedLevel-for-loop
}

void CpGrid::setRefinedLevelGridsGeometries( /* Refined corner arguments */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>>& refined_corners_vec,
                                             const std::vector<int>& refined_corner_count_vec,
                                             /* Refined face arguments */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>>& refined_faces_vec,
                                             std::vector<Dune::cpgrid::EntityVariableBase<enum face_tag>>& mutable_refined_face_tags_vec,
                                             std::vector<Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>>& mutable_refined_face_normals_vec,
                                             std::vector<Opm::SparseTable<int>>& refined_face_to_point_vec,
                                             const std::vector<int>& refined_face_count_vec,
                                             /* Refined cell argumets */
                                             std::vector<Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>>& refined_cells_vec,
                                             std::vector<std::vector<std::array<int,8>>>& refined_cell_to_point_vec,
                                             std::vector<std::vector<int>>& refined_global_cell_vec,
                                             std::vector<int>& refined_cell_count_vec,
                                             std::vector<cpgrid::OrientedEntityTable<0,1>>& refined_cell_to_face_vec,
                                             std::vector<cpgrid::OrientedEntityTable<1,0>>& refined_face_to_cell_vec,
                                             /* Auxiliary arguments */
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                                             const std::map<std::array<int,2>,std::array<int,2>>& refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner,
                                             const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                             const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                                             const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                             const std::vector<Dune::cpgrid::DefaultGeometryPolicy>& refined_geometries_vec,
                                             const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                             const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                             const std::vector<int>& assignRefinedLevel,
                                             const int& preAdaptMaxLevel,
                                             const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                             const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                             const std::vector<std::array<int,3>>&  cells_per_dim_vec) const
{
    // --- Refined corners  ---
    populateRefinedCorners(refined_corners_vec,
                           refined_corner_count_vec,
                           markedElem_to_itsLgr,
                           preAdaptMaxLevel,
                           refinedLevelAndRefinedCorner_to_elemLgrAndElemLgrCorner);
    // --- Refined faces  ---
    populateRefinedFaces(refined_faces_vec,
                         mutable_refined_face_tags_vec,
                         mutable_refined_face_normals_vec,
                         refined_face_to_point_vec,
                         refined_face_count_vec,
                         refinedLevelAndRefinedFace_to_elemLgrAndElemLgrFace,
                         elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                         vanishedRefinedCorner_to_itsLastAppearance,
                         markedElem_to_itsLgr,
                         preAdaptMaxLevel,
                         cornerInMarkedElemWithEquivRefinedCorner,
                         markedElemAndEquivRefinedCorn_to_corner);
    // --- Refined cells  ---
    populateRefinedCells(refined_cells_vec,
                         refined_cell_to_point_vec,
                         refined_global_cell_vec,
                         refined_cell_count_vec,
                         refined_cell_to_face_vec,
                         refined_face_to_cell_vec,
                         refinedLevelAndRefinedCell_to_elemLgrAndElemLgrCell,
                         elemLgrAndElemLgrFace_to_refinedLevelAndRefinedFace,
                         faceInMarkedElemAndRefinedFaces,
                         refined_geometries_vec,
                         elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                         vanishedRefinedCorner_to_itsLastAppearance,
                         markedElem_to_itsLgr,
                         assignRefinedLevel,
                         preAdaptMaxLevel,
                         markedElemAndEquivRefinedCorn_to_corner,
                         cornerInMarkedElemWithEquivRefinedCorner,
                         cells_per_dim_vec);

}

void CpGrid::updateLeafGridViewGeometries( /* Leaf grid View Corners arguments */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& adapted_corners,
                                           const int& corner_count,
                                           /* Leaf grid View Faces arguments */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& adapted_faces,
                                           Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags,
                                           Dune::cpgrid::EntityVariableBase<Dune::FieldVector<double,3>>& mutable_face_normals,
                                           Opm::SparseTable<int>& adapted_face_to_point,
                                           const int& face_count,
                                           /* Leaf grid View Cells argumemts  */
                                           Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& adapted_cells,
                                           std::vector<std::array<int,8>>& adapted_cell_to_point,
                                           const int& cell_count,
                                           cpgrid::OrientedEntityTable<0,1>& adapted_cell_to_face,
                                           cpgrid::OrientedEntityTable<1,0>& adapted_face_to_cell,
                                           /* Auxiliary arguments */
                                           const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                           const std::unordered_map<int,std::array<int,2>>& adaptedFace_to_elemLgrAndElemLgrFace,
                                           const std::unordered_map<int,std::array<int,2>>& adaptedCell_to_elemLgrAndElemLgrCell,
                                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrFace_to_adaptedFace,
                                           const std::vector<std::vector<std::pair<int, std::vector<int>>>>& faceInMarkedElemAndRefinedFaces,
                                           const Dune::cpgrid::DefaultGeometryPolicy& adapted_geometries,
                                           const std::map<std::array<int,2>,int>& elemLgrAndElemLgrCorner_to_adaptedCorner,
                                           const std::map<std::array<int,2>, std::array<int,2>>& vanishedRefinedCorner_to_itsLastAppearance,
                                           const std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& markedElem_to_itsLgr,
                                           const std::vector<int>& assignRefinedLevel,
                                           const std::map<std::array<int,2>,int>& markedElemAndEquivRefinedCorn_to_corner,
                                           const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                           const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                           const int& preAdaptMaxLevel) const
{
    // --- Adapted corners ---
    populateLeafGridCorners(adapted_corners,
                            corner_count,
                            markedElem_to_itsLgr,
                            adaptedCorner_to_elemLgrAndElemLgrCorner);
    // --- Adapted faces ---
    populateLeafGridFaces(adapted_faces,
                          mutable_face_tags,
                          mutable_face_normals,
                          adapted_face_to_point,
                          face_count,
                          adaptedFace_to_elemLgrAndElemLgrFace,
                          elemLgrAndElemLgrCorner_to_adaptedCorner,
                          vanishedRefinedCorner_to_itsLastAppearance,
                          markedElem_to_itsLgr,
                          assignRefinedLevel,
                          markedElemAndEquivRefinedCorn_to_corner,
                          cornerInMarkedElemWithEquivRefinedCorner,
                          cells_per_dim_vec,
                          preAdaptMaxLevel);
    // --- Adapted cells ---
    populateLeafGridCells(adapted_cells,
                          adapted_cell_to_point,
                          cell_count,
                          adapted_cell_to_face,
                          adapted_face_to_cell,
                          adaptedCell_to_elemLgrAndElemLgrCell,
                          elemLgrAndElemLgrFace_to_adaptedFace,
                          faceInMarkedElemAndRefinedFaces,
                          adapted_geometries,
                          elemLgrAndElemLgrCorner_to_adaptedCorner,
                          vanishedRefinedCorner_to_itsLastAppearance,
                          markedElem_to_itsLgr,
                          assignRefinedLevel,
                          markedElemAndEquivRefinedCorn_to_corner,
                          cornerInMarkedElemWithEquivRefinedCorner,
                          cells_per_dim_vec,
                          preAdaptMaxLevel);
}

void CpGrid::updateCornerHistoryLevels(const std::vector<std::vector<std::array<int,2>>>& cornerInMarkedElemWithEquivRefinedCorner,
                                       const std::map<std::array<int,2>,std::array<int,2>>& elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner,
                                       const std::unordered_map<int,std::array<int,2>>& adaptedCorner_to_elemLgrAndElemLgrCorner,
                                       const int& corner_count,
                                       const std::vector<std::array<int,2>>& preAdaptGrid_corner_history,
                                       const int& preAdaptMaxLevel,
                                       const int& newLevels)
{
    for (int level = preAdaptMaxLevel+1; level < preAdaptMaxLevel + newLevels+1; ++level) {
        currentData()[level]->corner_history_.resize( currentData()[level] ->size(3), std::array<int,2>({-1,-1}));
    }
    // corner_history_ for levels 0, level 1, ..., preAdapt-maxLevel (maximum level before calling (again) adapt) should be already populated
    // corner_history_[ new corner ] = {-1,-1}
    // corner_history_[ corner equivalent to a corner in a previous level ] = { level where the corner was born, its index in that level grid}.
    for (std::size_t corner = 0; corner < cornerInMarkedElemWithEquivRefinedCorner.size(); ++corner) {
        if (!cornerInMarkedElemWithEquivRefinedCorner[corner].empty()) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at(cornerInMarkedElemWithEquivRefinedCorner[corner].back());
            currentData()[refinedLevel]->corner_history_[refinedCorner] = preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, static_cast<int>(corner)}} :  preAdaptGrid_corner_history[corner];
        }
    }
   
    // corner_history_ leaf grid view
    for ( int leafCorner = 0; leafCorner < corner_count; ++leafCorner){
        currentData().back()->corner_history_.resize(corner_count);
        const auto& [elemLgr, elemLgrCorner] = adaptedCorner_to_elemLgrAndElemLgrCorner.at(leafCorner);
        if (elemLgr != -1) {
            const auto& [refinedLevel, refinedCorner] = elemLgrAndElemLgrCorner_to_refinedLevelAndRefinedCorner.at({elemLgr, elemLgrCorner});
            currentData().back()->corner_history_[leafCorner] = { refinedLevel, refinedCorner};
        }
        else {
            currentData().back()->corner_history_[leafCorner] =  preAdaptGrid_corner_history.empty() ? std::array<int,2>{{0, elemLgrCorner}} : preAdaptGrid_corner_history[elemLgrCorner];
        }
    }
}


std::array<int,3>  CpGrid::getRefinedCornerIJK(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    std::array<int,3> ijk;
    ijk[2] = cornerIdxInLgr % (cells_per_dim[2] +1);
    cornerIdxInLgr -= ijk[2];
    cornerIdxInLgr /= (cells_per_dim[2] +1);
    ijk[0] = cornerIdxInLgr % (cells_per_dim[0]+1);
    cornerIdxInLgr -=ijk[0];
    ijk[1] = cornerIdxInLgr / (cells_per_dim[0]+1);
    return ijk;
}

std::array<int,3>  CpGrid::getRefinedFaceIJK(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                             const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{
    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(faceIdxInLgr, true);
    const auto& faceTag = elemLgr_ptr ->face_tag_[faceEntity];
    std::array<int,3> ijk;
    switch (faceTag) {
    case I_FACE:
        faceIdxInLgr -= (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1));
        // faceIdxInLgr =  (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -= ijk[1]; // (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1])
        faceIdxInLgr /= cells_per_dim[1]; // (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -=ijk[2]; // i*cells_per_dim[2]
        ijk[0] = faceIdxInLgr / cells_per_dim[2];
        break;
    case J_FACE:
        faceIdxInLgr -=  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
            + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2]);
        // faceIdxInLgr =  (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
        ijk[2] = faceIdxInLgr % cells_per_dim[2];
        faceIdxInLgr -= ijk[2]; // (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2])
        faceIdxInLgr /= cells_per_dim[2]; // (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -=ijk[0]; // j*cells_per_dim[0]
        ijk[1] = faceIdxInLgr / cells_per_dim[0];
        break;
    case K_FACE:
        //  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
        ijk[0] = faceIdxInLgr % cells_per_dim[0];
        faceIdxInLgr -= ijk[0]; // (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0])
        faceIdxInLgr /= cells_per_dim[0]; // (k*cells_per_dim[1]) + j
        ijk[1] = faceIdxInLgr % cells_per_dim[1];
        faceIdxInLgr -=ijk[1]; // k*cells_per_dim[1]
        ijk[2] = faceIdxInLgr / cells_per_dim[1];
        break;
    default:
        OPM_THROW(std::logic_error, "FaceTag is not I, J, or K!");
    }
    return ijk;
}

bool CpGrid::isRefinedCornerInInteriorLgr(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    assert(cells_per_dim[0]>0);
    assert(cells_per_dim[1]>0);
    assert(cells_per_dim[2]>0);

    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    return ((ijk[0]%cells_per_dim[0] > 0) &&  (ijk[1]%cells_per_dim[1]>0) && (ijk[2]%cells_per_dim[2]>0));
}

bool CpGrid::isRefinedFaceInInteriorLgr(const std::array<int,3>& cells_per_dim, int faceIdxInLgr, const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    return ((ijk[0]%cells_per_dim[0] > 0 && isIface) ||  (ijk[1]%cells_per_dim[1]>0 && isJface) || (ijk[2]%cells_per_dim[2]>0 && isKface));
}

bool CpGrid::isRefinedNewBornCornerOnLgrBoundary(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_I_FACE = isOnParentCell_I_FACEfalse_and_newBornCorn || isOnParentCell_I_FACEtrue_and_newBornCorn;
    bool isOnParentCell_J_FACE = isOnParentCell_J_FACEfalse_and_newBornCorn || isOnParentCell_J_FACEtrue_and_newBornCorn;
    bool isOnParentCell_K_FACE = isOnParentCell_K_FACEfalse_and_newBornCorn || isOnParentCell_K_FACEtrue_and_newBornCorn;
    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

bool CpGrid::newRefinedCornerLiesOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr) const
{
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    bool isOnEdge = isNewBornOnEdge01 || isNewBornOnEdge23 || isNewBornOnEdge02 || isNewBornOnEdge13 ||
        isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge15 || isNewBornOnEdge37 ||
        isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;

    return isOnEdge;
}


bool CpGrid::isRefinedFaceOnLgrBoundary(const std::array<int,3>& cells_per_dim, int faceIdxInLgr,
                                        const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr) const
{
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);

    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];

    bool isKface = (faceIdxInLgr < refined_k_faces);
    bool isIface = (faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces);
    bool isJface = (faceIdxInLgr >= refined_k_faces + refined_i_faces);

    bool isOnParentCell_I_FACE = isIface && (ijk[0] % cells_per_dim[0] == 0) && (ijk[1]<cells_per_dim[1]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_J_FACE = isJface && (ijk[1] % cells_per_dim[1] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[2]<cells_per_dim[2]);
    bool isOnParentCell_K_FACE = isKface && (ijk[2] % cells_per_dim[2] == 0) && (ijk[0]<cells_per_dim[0]) && (ijk[1]<cells_per_dim[1]);

    return (isOnParentCell_I_FACE || isOnParentCell_J_FACE || isOnParentCell_K_FACE);
}

std::array<int,2> CpGrid::getParentFacesAssocWithNewRefinedCornLyingOnEdge(const std::array<int,3>& cells_per_dim, int cornerIdxInLgr, int elemLgr) const
{
    assert(newRefinedCornerLiesOnEdge(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(elemLgr, true)];
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    // Corners Order defined in Geometry::refine  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);
    // Edges laying on bottom face
    bool isNewBornOnEdge01 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == 0);
    bool isNewBornOnEdge23 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == 0);
    bool isNewBornOnEdge02 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);
    bool isNewBornOnEdge13 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == 0);

    // Edges connecting bottom and top faces
    bool isNewBornOnEdge04 = (ijk[0] == 0) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge26 = (ijk[0] == 0) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge15 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == 0) && (ijk[2] % cells_per_dim[2] != 0);
    bool isNewBornOnEdge37 = (ijk[0] == cells_per_dim[0]) && (ijk[1] == cells_per_dim[1]) && (ijk[2] % cells_per_dim[2] != 0);

    // Edges laying on top face
    bool isNewBornOnEdge45 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge67 = (ijk[0] % cells_per_dim[0] != 0) && (ijk[1] == cells_per_dim[1]) && ( ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge46 = (ijk[0] == 0) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);
    bool isNewBornOnEdge57 = (ijk[0] == cells_per_dim[0]) && (ijk[1] % cells_per_dim[1] != 0) && (ijk[2] == cells_per_dim[2]);

    std::vector<int> auxFaces;
    auxFaces.reserve(2);

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        // Add I_FACE false
        bool addIfalse = isNewBornOnEdge02 || isNewBornOnEdge04 || isNewBornOnEdge26 || isNewBornOnEdge46;
        if( addIfalse && (faceTag == 0)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE false
        bool addJfalse = isNewBornOnEdge01 || isNewBornOnEdge04 || isNewBornOnEdge15 || isNewBornOnEdge45;
        if( addJfalse && (faceTag == 1)  && (!face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE false
        bool addKfalse = isNewBornOnEdge01 ||  isNewBornOnEdge13 || isNewBornOnEdge23 || isNewBornOnEdge02;
        if( addKfalse && (faceTag == 2) && (!face.orientation())) {
            auxFaces.push_back(face.index());
        }
        // Add I_FACE true
        bool addItrue = isNewBornOnEdge13 || isNewBornOnEdge15 || isNewBornOnEdge37 || isNewBornOnEdge57;
        if( addItrue && (faceTag == 0)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add J_FACE true
        bool addJtrue = isNewBornOnEdge23|| isNewBornOnEdge26 || isNewBornOnEdge37 || isNewBornOnEdge67;
        if( addJtrue && (faceTag == 1)  && (face.orientation()))  {
            auxFaces.push_back(face.index());
        }
        // Add K_FACE true
        bool addKtrue = isNewBornOnEdge45 || isNewBornOnEdge67 || isNewBornOnEdge46 || isNewBornOnEdge57;
        if(addKtrue && (faceTag == 2) && (face.orientation())) {
            auxFaces.push_back(face.index());
        }
    }
    return {auxFaces[0], auxFaces[1]};
}

int CpGrid::getParentFaceWhereNewRefinedCornerLiesOn(const std::array<int,3>& cells_per_dim,
                                                     int cornerIdxInLgr, int elemLgr) const
{
    assert(isRefinedNewBornCornerOnLgrBoundary(cells_per_dim, cornerIdxInLgr));

    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(elemLgr, true)];
    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    const auto& ijk = getRefinedCornerIJK(cells_per_dim, cornerIdxInLgr);

    bool isOnParentCell_I_FACEfalse_and_newBornCorn = ( (ijk[0] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_I_FACEtrue_and_newBornCorn = ( (ijk[0] == cells_per_dim[0]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEfalse_and_newBornCorn = ( (ijk[1] == 0) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_J_FACEtrue_and_newBornCorn = ( (ijk[1] == cells_per_dim[1]) && ((ijk[0] % cells_per_dim[0] != 0) || (ijk[2] % cells_per_dim[2] !=0) ));
    bool isOnParentCell_K_FACEfalse_and_newBornCorn = ( (ijk[2] == 0) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));
    bool isOnParentCell_K_FACEtrue_and_newBornCorn = ( (ijk[2] == cells_per_dim[2]) && ((ijk[1] % cells_per_dim[1] != 0) || (ijk[0] % cells_per_dim[0] !=0) ));

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (isOnParentCell_I_FACEfalse_and_newBornCorn && (faceTag == 0) && !face.orientation()) { // I_FACE false
            return face.index();
        }
        if (isOnParentCell_I_FACEtrue_and_newBornCorn && (faceTag == 0) && face.orientation()) { // I_FACE true
            return face.index();
        }
        if (isOnParentCell_J_FACEfalse_and_newBornCorn && (faceTag == 1) && !face.orientation()) { // J_FACE false
            return face.index();
        }
        if (isOnParentCell_J_FACEtrue_and_newBornCorn && (faceTag == 1) && face.orientation()) { // J_FACE true
            return face.index();
        }
        if (isOnParentCell_K_FACEfalse_and_newBornCorn && (faceTag == 2) && !face.orientation()) { // K_FACE false
            return face.index();
        }
        if (isOnParentCell_K_FACEtrue_and_newBornCorn && (faceTag == 2) && face.orientation()) { // K_FACE true
            return face.index();
        }
    }
    OPM_THROW(std::logic_error, "Cannot find parent face index where new refined corner lays on.");
}

int CpGrid::getParentFaceWhereNewRefinedFaceLiesOn(const std::array<int,3>& cells_per_dim,
                                                   int faceIdxInLgr,
                                                   const std::shared_ptr<cpgrid::CpGridData>& elemLgr_ptr,
                                                   int elemLgr) const
{
    assert(isRefinedFaceOnLgrBoundary(cells_per_dim, faceIdxInLgr, elemLgr_ptr));
    const auto& ijk = getRefinedFaceIJK(cells_per_dim, faceIdxInLgr, elemLgr_ptr);
    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(elemLgr, true)];

    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }

    // Order defined in Geometry::refine
    // K_FACES  (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    int refined_k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    int refined_i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    int refined_j_faces = cells_per_dim[0]*(cells_per_dim[1]+1)*cells_per_dim[2];
    assert( faceIdxInLgr < refined_k_faces + refined_i_faces + refined_j_faces);

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (faceIdxInLgr <  refined_k_faces ) { // It's a K_FACE
            if ((ijk[2] == 0) && (faceTag == 2) && !face.orientation()) { // {K_FACE, false}
                return face.index();
            }
            if ((ijk[2] == cells_per_dim[2]) && (faceTag == 2) && face.orientation()) { // {K_FACE, true}
                return face.index();
            }
        }
        if ((faceIdxInLgr >= refined_k_faces) && (faceIdxInLgr < refined_k_faces + refined_i_faces)) { // It's I_FACE
            if ((ijk[0] == 0) && (faceTag == 0) && !face.orientation()) { // {I_FACE, false}
                return face.index();
            }
            if ((ijk[0] == cells_per_dim[0]) && (faceTag == 0) && face.orientation()) { // {I_FACE, true}
                return face.index();
            }
        }
        if (faceIdxInLgr >= refined_k_faces + refined_i_faces) {// It's J_FACE
            if ((ijk[1] == 0) && (faceTag == 1) && !face.orientation()) { // {J_FACE, false}
                return face.index();
            }
            if ((ijk[1] == cells_per_dim[1]) && (faceTag == 1) && face.orientation()) { // {J_FACE, true}
                return face.index();
            }
        }
    }
    OPM_THROW(std::logic_error, "Cannot find parent face index where the new refined face lays on.");
}

int CpGrid::replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1, int cornerIdxLgr1, const std::array<int,3>& cells_per_dim_lgr2) const
{
    const auto& ijkLgr1 = getRefinedCornerIJK(cells_per_dim_lgr1, cornerIdxLgr1);
    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k
    if (ijkLgr1[0] == cells_per_dim_lgr1[0]) {
        return   (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[1] == cells_per_dim_lgr1[1]) {
        return  (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == cells_per_dim_lgr1[2]) {
        return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1));
    }
    else {
        OPM_THROW(std::logic_error, "Cannot convert corner index from one LGR to its neighboring LGR.");
    }
}

int CpGrid::replaceLgr1CornerIdxByLgr2CornerIdx(const std::array<int,3>& cells_per_dim_lgr1, int cornerIdxLgr1, int elemLgr1, int parentFaceLastAppearanceIdx,
                                                const std::array<int,3>& cells_per_dim_lgr2) const
{
    assert(newRefinedCornerLiesOnEdge(cells_per_dim_lgr1, cornerIdxLgr1));
    const auto& faces = getParentFacesAssocWithNewRefinedCornLyingOnEdge(cells_per_dim_lgr1, cornerIdxLgr1, elemLgr1);
    assert( (faces[0] == parentFaceLastAppearanceIdx) || (faces[1] == parentFaceLastAppearanceIdx)); 
 
    const auto& ijkLgr1 = getRefinedCornerIJK(cells_per_dim_lgr1, cornerIdxLgr1);
    const auto& parentCell_to_face = current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(elemLgr1, true)];

    if(parentCell_to_face.size()>6){
        OPM_THROW(std::logic_error, "The associted parent cell has more than six faces. Refinment/Adaptivity not supported yet.");
    }
    // Since lgr1 is associated with an element index smaller than "a last appearance lgr", then the only possibilities are
    // I_FACE true, J_FACE true, K_FACE true

    // Order defined in Geometry::refine
    //  (j*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (i*(cells_per_dim[2]+1)) + k

    for (const auto& face : parentCell_to_face) {
        const auto& faceEntity =  Dune::cpgrid::EntityRep<1>(face.index(), true);
        const auto& faceTag = current_view_data_->face_tag_[faceEntity];
        if (parentFaceLastAppearanceIdx == face.index() && face.orientation()) {
            if (faceTag == 0) { // I_FACE true
                // The same new born refined corner will have equal values of j and k, but i == 0 instead of cells_per_dim[0]
                return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1))  + ijkLgr1[2];
            }
            if (faceTag == 1) { // J_FACE true
                // The same new born refined corner will have equal values of i and k, but j == 0 instead of cells_per_dim[1]
                return  (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1)) + ijkLgr1[2];
            }
            if (faceTag == 2) { // K_FACE true
                // The same new born refined corner will have equal values of i and j, but k == 0 instead of cells_per_dim[2]
                return  (ijkLgr1[1]*(cells_per_dim_lgr2[0]+1)*(cells_per_dim_lgr2[2]+1)) + (ijkLgr1[0]*(cells_per_dim_lgr2[2]+1));
            }
        }
    }
    OPM_THROW(std::logic_error, "Cannot convert corner index from one LGR to its neighboring LGR.");
}

int  CpGrid::replaceLgr1FaceIdxByLgr2FaceIdx(const std::array<int,3>& cells_per_dim_lgr1, int faceIdxInLgr1,
                                             const std::shared_ptr<Dune::cpgrid::CpGridData>& elemLgr1_ptr,
                                             const std::array<int,3>& cells_per_dim_lgr2) const
{
    const auto& ijkLgr1 = getRefinedFaceIJK(cells_per_dim_lgr1, faceIdxInLgr1, elemLgr1_ptr);
    // lgr1 represents an element index < lgr2 (neighboring cells sharing a face with lgr1-element)
    // Order defined in Geometry::refine
    // K_FACES (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
    // I_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))
    //           + (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j
    // J_FACES  (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2] +1))
    //                    + ((cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2])
    //                    + (j*cells_per_dim[0]*cells_per_dim[2]) + (i*cells_per_dim[2]) + k
    const int& kFacesLgr2 = cells_per_dim_lgr2[0]*cells_per_dim_lgr2[1]*(cells_per_dim_lgr2[2]+1);
    const int& iFacesLgr2 = ((cells_per_dim_lgr2[0]+1)*cells_per_dim_lgr2[1]*cells_per_dim_lgr2[2]);

    if (ijkLgr1[0] == cells_per_dim_lgr1[0]) {
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        return  kFacesLgr2 + (ijkLgr1[2]*cells_per_dim_lgr2[1]) + ijkLgr1[1];
    }
    if (ijkLgr1[1] == cells_per_dim_lgr1[1]) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[2] == cells_per_dim_lgr2[2]);
        return kFacesLgr2 + iFacesLgr2 + (ijkLgr1[0]*cells_per_dim_lgr2[2]) + ijkLgr1[2];
    }
    if (ijkLgr1[2] == cells_per_dim_lgr1[2]) {
        assert( cells_per_dim_lgr1[0] == cells_per_dim_lgr2[0]);
        assert( cells_per_dim_lgr1[1] == cells_per_dim_lgr2[1]);
        return  (ijkLgr1[1]*cells_per_dim_lgr2[0]) + ijkLgr1[0];
    }
    else {
        OPM_THROW(std::logic_error, "Cannot convert face index from one LGR to its neighboring LGR.");
    }
}




} // namespace Dune


