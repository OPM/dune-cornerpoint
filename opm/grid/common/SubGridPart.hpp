// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=2 sts=2:

// Note: this file is based on defaultgridview.hh from Dune, and is therefore
// licensed under the Dune license (GPLv2 + runtime exception),
// see https://dune-project.org/about/license/
// rather than the OPM license (GPLv3+)
// Copyright 2021 Dune contributors.
// Copyright 2021 SINTEF Digital, Mathematics and Cybernetics.

#ifndef OPM_SUBGRIDPART_HEADER
#define OPM_SUBGRIDPART_HEADER

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Dune
{

template <class GridImp>
class SubGridPart;

template <class GridImp>
struct SubGridPartTraits {
    using GridPartImp = SubGridPart<GridImp>;

    /** \brief type of the grid */
    using Grid = typename std::remove_const<GridImp>::type;

    /** \brief type of the index set */
    using IndexSet = typename Grid ::Traits ::LeafIndexSet;

    /** \brief type of the intersection */
    using Intersection = typename Grid ::Traits ::LeafIntersection;

    /** \brief type of the intersection iterator */
    using IntersectionIterator = typename Grid ::Traits ::LeafIntersectionIterator;

    /** \brief type of the collective communication */
    using CollectiveCommunication = typename Grid ::Traits ::Communication;


    template <class BaseEntityType>
    class SubEntity : public BaseEntityType
    {
    public:
        SubEntity()
            : BaseEntityType()
            , is_owned_(false)
        {
        }
        SubEntity(const BaseEntityType& base, const bool owned)
            : BaseEntityType(base)
            , is_owned_(owned)
        {
        }
        auto partitionType() const
        {
            if (is_owned_) {
                return Dune::InteriorEntity;
            } else {
                return Dune::OverlapEntity;
            }
        }
    private:
        bool is_owned_;
    };

    template <int cd>
    struct Codim {
        using BaseEntity = typename Grid ::Traits ::template Codim<cd>::Entity;
        using Entity = SubEntity<BaseEntity>;
        using EntitySeed = typename Grid ::Traits ::template Codim<cd>::EntitySeed;

        using Geometry = typename Grid ::template Codim<cd>::Geometry;
        using LocalGeometry = typename Grid ::template Codim<cd>::LocalGeometry;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template <PartitionIteratorType pit>
        struct Partition {
            /** \brief iterator over a given codim and partition type */
            using BaseIterator = typename Grid ::template Codim<cd>::template Partition<pit>::LeafIterator;
        };
    };

    enum { conforming = Capabilities ::isLeafwiseConforming<Grid>::v };
};


/// \brief A class to represent a part of a grid, similar to a GridView.
///
/// The differences from a GridView are:
///   - The SubGridPart consists of a set of elements (codim 0 entities),
///     considered to be Interior, and their neighbours, considered to be Overlap.
///   - When iterating over intersections on the elements, and accessing the outside()
///     elements, this can give you access to grid entities that are not in the
///     SubGridPart itself. This can only happen for intersections in the Overlap part.
///     For intersections of elements in the Interior part, the outside() element will
///     be either Interior or Overlap.
template <class GridImp>
class SubGridPart
{
    using ThisType = SubGridPart<GridImp>;

public:
    using Traits = SubGridPartTraits<GridImp>;

    /** \brief type of the grid */
    using Grid = typename Traits::Grid;

    /** \brief type of the index set */
    using IndexSet = typename Traits ::IndexSet;

    /** \brief type of the intersection */
    using Intersection = typename Traits ::Intersection;

    /** \brief type of the intersection iterator */
    using IntersectionIterator = typename Traits ::IntersectionIterator;

    /** \brief type of the collective communication */
    using CollectiveCommunication = typename Traits ::CollectiveCommunication;


    /** \brief Codim Structure */
    template <int cd>
    struct Codim : public Traits ::template Codim<cd> {
        using Entity = typename Traits::template Codim<cd>::Entity;
        class SubIterator
        {
        public:
            SubIterator(const SubGridPart& view, std::size_t index)
                : view_(&view)
                , index_(index)
            {
            }
            const Entity& operator*() const
            {
                entity_ = this->view_->get(index_);
                return entity_;
            }
            const Entity* operator->() const
            {
                entity_ = this->view_->get(index_);
                return &entity_;
            }
            SubIterator operator++()
            {
                ++index_;
                return *this;
            }
            SubIterator operator++(int)
            {
                Iterator copy(*this);
                ++index_;
                return copy;
            }
            bool operator==(const SubIterator& other) const
            {
                assert(view_ == other.view_);
                return index_ == other.index_;
            }
            bool operator!=(const SubIterator& other) const
            {
                assert(view_ == other.view_);
                return index_ != other.index_;
            }
        private:
            const SubGridPart* view_;
            std::size_t index_;
            mutable Entity entity_; // This may be low-performing for grids with large Entity objects.
        };
        using Iterator = SubIterator;

        /** \brief Define types needed to iterate over entities of a given partition type */
        template <PartitionIteratorType pit>
        struct Partition {
            /** \brief iterator over a given codim and partition type */
            using Iterator = SubIterator;
        };

    };

    enum { conforming = Traits::conforming };
    enum { dimension = GridImp::dimension };

public:

    /// Construct a view of the codim 0 entities that can be constructed from the seeds input.
    ///
    /// The seeds input is moved from and will be in a valid but indeterminate state after the call.
    SubGridPart(const Grid& grid,
                std::vector<typename Codim<0>::Entity::EntitySeed>&& seeds,
                const bool overlap = true)
        : grid_(&grid)
        , subset_(std::move(seeds))
        , num_owned_(subset_.size())
    {
        // Nothing more to do if we do not want to have overlap entities.
        if (!overlap) {
            return;
        }

        // Add neighbouring not-owned entities to subset_
        using Seed = typename Codim<0>::Entity::EntitySeed;
        std::unordered_set<int> owned;
        std::unordered_map<int, Seed> neighbors;
        const auto& iset = grid_->leafIndexSet();
        const auto& leaf_view = grid_->leafGridView();
        for (const auto& seed : subset_) {
            // Add this entity to the set of owned indices.
            const auto& entity = grid_->entity(seed);
            owned.insert(iset.index(entity));
            // Iterating over all intersections, ...
            const auto end = leaf_view.iend(entity);
            for (auto it = leaf_view.ibegin(entity); it != end; ++it) {
                if (it->boundary()) {
                    continue;
                }
                if (it->neighbor()) {
                    const auto outside_entity = it->outside();
                    // ...for all neighbour entities, add to neighbors.
                    neighbors.try_emplace(iset.index(outside_entity), outside_entity.seed());
                }
            }
        }
        // Now that owned is complete, we can eliminate any owned entries.
        std::map<int, Seed> unowned_neighbors;
        for (const auto& nb : neighbors) {
            if (owned.count(nb.first) == 0) {
                unowned_neighbors.insert(nb);
            }
        }
        subset_.resize(subset_.size() + unowned_neighbors.size());
        std::size_t count = num_owned_;
        for (const auto& neighbor : unowned_neighbors) {
            subset_[count] = neighbor.second;
            ++count;
        }
        assert(count == subset_.size());
    }

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid& grid() const
    {
        assert(grid_);
        return *grid_;
    }

    /** \brief obtain the index set */
    // const IndexSet& indexSet() const // Not implemented

    /** \brief obtain number of entities in a given codimension */
    int size(int codim) const
    {
        if (codim == 0) {
            return subset_.size();
        } else {
            return 0;
        }
    }

    /** \brief obtain number of entities with a given geometry type */
    // int size(const GeometryType& type) const // Not implemented

    /** \brief obtain begin iterator for this view */
    template <int cd>
    typename Codim<cd>::Iterator begin() const
    {
        static_assert(cd == 0, "Only codimension 0 iterators for SubGridPart.");
        using Iterator = typename Codim<cd>::Iterator;
        return Iterator(*this, 0);
    }

    /** \brief obtain end iterator for this view */
    template <int cd>
    typename Codim<cd>::Iterator end() const
    {
        static_assert(cd == 0, "Only codimension 0 iterators for SubGridPart.");
        using Iterator = typename Codim<cd>::Iterator;
        return Iterator(*this, subset_.size());
    }

    // We support iterating over Interior_Partition, Overlap_Partition and All_Partition

    /** \brief obtain begin iterator for this view */
    template <int cd, PartitionIteratorType pit>
    typename Codim<cd>::template Partition<pit>::Iterator begin() const
    {
        static_assert(cd == 0, "Only codimension 0 iterators for SubGridPart.");
        static_assert(pit == Interior_Partition || pit == Overlap_Partition || pit == All_Partition);
        if constexpr (pit == Interior_Partition || pit == All_Partition) {
            return begin<0>();
        } else {
            // Overlap partition starts at index num_owned_.
            // Note that it may be empty, i.e. begin() == end().
            return typename Codim<cd>::Iterator(*this, num_owned_);
        }
    }

    /** \brief obtain end iterator for this view */
    template <int cd, PartitionIteratorType pit>
    typename Codim<cd>::template Partition<pit>::Iterator end() const
    {
        static_assert(cd == 0, "Only codimension 0 iterators for SubGridPart.");
        static_assert(pit == Interior_Partition || pit == Overlap_Partition || pit == All_Partition);
        if constexpr (pit == Overlap_Partition || pit == All_Partition) {
            return end<0>();
        } else {
            // Interior partition ends before index num_owned_.
            return typename Codim<cd>::Iterator(*this, num_owned_);
        }
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
    {
        return entity.impl().ileafbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
    {
        return entity.impl().ileafend();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication& comm() const
    {
        return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
        if (codim == 0) {
            return subset_.size() - num_owned_;
        } else {
            return 0;
        }
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize([[maybe_unused]] int codim) const
    {
        return 0;
    }

    /** communicate data on this view */
    // template <class DataHandleImp, class DataType>
    // void
    // communicate(CommDataHandleIF<DataHandleImp, DataType>& data, InterfaceType iftype, CommunicationDirection dir) const
    // Not implemented.

private:
    using Entity0 = typename Codim<0>::Entity;
    Entity0 get(std::size_t ii) const
    {
        const bool owned = ii < num_owned_;
        return Entity0(grid_->entity(subset_[ii]), owned);
    }
    const Grid* grid_;
    std::vector<typename Entity0::EntitySeed> subset_;
    const std::size_t num_owned_;
};

} // namespace Dune

#endif // OPM_SUBGRIDPART_HEADER
