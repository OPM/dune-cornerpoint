/*
  Copyright 2025 SINTEF Digital

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

#include <opm/grid/utility/createThreadIterators.hpp>
#include <opm/grid/utility/StopWatch.hpp>
#include <opm/grid/CpGrid.hpp>

#include <omp.h>

#include <iostream>
#include <vector>


template <class GridView>
class GridThreadedIteration
{
public:
    GridThreadedIteration(const GridView& gv, const std::size_t num_threads)
        : gv_(gv), nt_(num_threads)
    {
        [[maybe_unused]] std::size_t chunk_size;
        grid_chunk_iterators_ = Opm::createThreadIterators(gv_, nt_, 1000, chunk_size);
        const std::size_t num_chunks = grid_chunk_iterators_.size() - 1;
        chunks_.reserve(num_chunks);
        for (std::size_t ii = 0; ii < num_chunks; ++ii) {
            chunks_.emplace_back(grid_chunk_iterators_[ii], grid_chunk_iterators_[ii + 1]);
        }
    }
    using Iterator = typename GridView::template Codim<0>::Iterator;

    struct Chunk
    {
        Chunk() {}
        Chunk(const Iterator& i1, const Iterator& i2)
            : pi_(i1, i2)
        {
        }
        Iterator begin() const { return pi_.first; }
        Iterator end() const { return pi_.second; }
        std::pair<Iterator, Iterator> pi_;
    };

    auto begin() const
    {
        return chunks_.begin();
    }
    auto end() const
    {
        return chunks_.end();
    }
    auto size() const
    {
        return chunks_.size();
    }

private:
    const GridView& gv_;
    const std::size_t nt_;
    std::vector<Iterator> grid_chunk_iterators_;
    std::vector<Chunk> chunks_;
};


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    std::array<int, 3> dims = { 200, 200, 100 };
    std::cout << "Creating grid with " << dims[0]*dims[1]*dims[2]/1000000 << "M cells." << std::endl;
    std::array<double, 3> cellsz = { 1.0, 1.0, 1.0 };
    Dune::CpGrid grid;
    grid.createCartesian(dims, cellsz);
    const auto& gv = grid.leafGridView();

    // Plain iteration
    {
        std::cout << "Running plain loop test." << std::endl;
        std::vector<double> vols(grid.size(0));
        Opm::time::StopWatch clock;
        clock.start();
        for (const auto& elem : elements(gv)) {
            vols[elem.index()] = elem.geometry().volume();
        }
        clock.stop();
        std::cout << "Time: " << clock.secsSinceLast() << std::endl;
    }

    // Iteration with chunks
    {
        const int num_threads = 5;
        std::cout << "Running chunked loop test with " << num_threads << " threads." << std::endl;
        std::vector<double> vols(grid.size(0));
        omp_set_num_threads(num_threads);
        Opm::time::StopWatch clock;
        clock.start();
        [[maybe_unused]] std::size_t chunk_size;
        const auto grid_chunk_iterators = Opm::createThreadIterators(gv, num_threads, 1000, chunk_size);
        const std::size_t num_chunks = grid_chunk_iterators.size() - 1;
#pragma omp parallel for
        for (std::size_t chunk = 0; chunk < num_chunks; ++chunk) {
            for (auto it = grid_chunk_iterators[chunk]; it != grid_chunk_iterators[chunk+1]; ++it) {
                const auto& elem = *it;
                vols[elem.index()] = elem.geometry().volume();
            }
        }
        clock.stop();
        std::cout << "Time: " << clock.secsSinceLast() << std::endl;
    }

    // Iteration with chunks using helper class
    {
        const int num_threads = 5;
        std::cout << "Running chunked loop test with " << num_threads << " threads and helper class." << std::endl;
        std::vector<double> vols(grid.size(0));
        omp_set_num_threads(num_threads);
        Opm::time::StopWatch clock;
        clock.start();
        GridThreadedIteration gti(gv, num_threads);
#pragma omp parallel for
        for (const auto& chunk : gti) {
            for (const auto& elem : chunk) {
                vols[elem.index()] = elem.geometry().volume();
            }
        }
        clock.stop();
        std::cout << "Time: " << clock.secsSinceLast() << std::endl;
    }

}
