/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_CREATETHREADITERATORS_HEADER_INCLUDED
#define OPM_CREATETHREADITERATORS_HEADER_INCLUDED

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace Opm
{

    /// Create a vector containing a spread of iterators into the
    /// elements of the range, to facilitate for example OpenMP
    /// parallelization of iterations over the elements.
    /// While the original range may have only forward iterators to
    /// the individual elements, the returned vector will provide
    /// random access to the chunks.
    /// \tparam     Range          Range of elements that supports (multipass) forward iteration.
    /// \param[in]  r              Range to be iterated over.
    /// \param[in]  num_elem       The number of elements in r.
    /// \param[in]  num_threads    The number of threads to target.
    /// \param[in]  max_chunk_size The maximum allowed chunk size.
    /// \param[out] chunk_size     The chunk size found by the algorithm.
    /// \return                    A vector of iterators, with the range's begin and end iterators
    ///                            being the first and last ones, and filling in iterators in between
    ///                            such that the distance from one to the next is chunk_size, except
    ///                            for possibly the last interval. If num_threads is 1, there will only
    ///                            be the begin and end iterators and no chunks in between.
    template <class Range>
    auto createThreadIterators(const Range& r,
                               const std::size_t num_elem,
                               const std::size_t num_threads,
                               const std::size_t max_chunk_size,
                               std::size_t& chunk_size)
    {
        if (num_threads < 1) {
            throw std::logic_error("createThreadIterators() called with num_threads arguments = "
                                   + std::to_string(num_threads));
        }
        std::vector<decltype(std::begin(r))> chunk_iterators;
        auto it = std::begin(r);
        const auto end = std::end(r);
        if (num_threads == 1) {
            chunk_iterators.push_back(it);
            chunk_iterators.push_back(end);
        } else {
            chunk_size = std::clamp(num_elem / num_threads, 1ul, max_chunk_size);
            chunk_iterators.reserve(num_elem / chunk_size + 2);
            for (int count = 0; it != end; ++it, ++count) {
                if (count % chunk_size == 0) {
                    chunk_iterators.push_back(it);
                }
            }
            chunk_iterators.push_back(end);
        }
        return chunk_iterators;
    }


    /// Convenience function for using createThreadIterators() with
    /// the cells in a grid view.
    template <class GridView>
    auto createThreadIterators(const GridView& gv,
                               const std::size_t num_threads,
                               const std::size_t max_chunk_size,
                               std::size_t& chunk_size)
    {
        return createThreadIterators(elements(gv),
                                     gv.size(0),
                                     num_threads,
                                     max_chunk_size,
                                     chunk_size);
    }

} // namespace Opm

#endif // OPM_CREATETHREADITERATORS_HEADER_INCLUDED
