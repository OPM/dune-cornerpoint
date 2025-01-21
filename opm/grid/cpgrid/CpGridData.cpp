#include"config.h"
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include"CpGridData.hpp"
#include"DataHandleWrappers.hpp"
#include"Intersection.hpp"
#include"Entity.hpp"
#include"OrientedEntityTable.hpp"
#include"Indexsets.hpp"
#include"PartitionTypeIndicator.hpp"

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/common/GridPartitioning.hpp>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/enumset.hh>
#include <opm/grid/utility/SparseTable.hpp>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
namespace cpgrid
{

CpGridData::CpGridData(const CpGridData& g)
    : index_set_(new IndexSet(g.cell_to_face_.size(), g.geomVector<3>().size())),
      local_id_set_(new IdSet(*this)),
      global_id_set_(new LevelGlobalIdSet(local_id_set_, this)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      ccobj_(g.ccobj_), use_unique_boundary_ids_(g.use_unique_boundary_ids_)
#if HAVE_MPI
    , cell_comm_(g.ccobj_)
#endif
{
#if HAVE_MPI
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
}

CpGridData::CpGridData(std::vector<std::shared_ptr<CpGridData>>& data)
    : index_set_(new IndexSet()), local_id_set_(new IdSet(*this)),
      global_id_set_(new LevelGlobalIdSet(local_id_set_, this)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      level_data_ptr_(),
      ccobj_(Dune::MPIHelper::getCommunicator()), use_unique_boundary_ids_(false)
#if HAVE_MPI
    , cell_comm_(Dune::MPIHelper::getCommunicator())
#endif
{
#if HAVE_MPI
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
    level_data_ptr_ = &data;
}

CpGridData::CpGridData(MPIHelper::MPICommunicator comm,  std::vector<std::shared_ptr<CpGridData>>& data)
    : index_set_(new IndexSet()), local_id_set_(new IdSet(*this)),
      global_id_set_(new LevelGlobalIdSet(local_id_set_, this)), partition_type_indicator_(new PartitionTypeIndicator(*this)),
      level_data_ptr_(),
      ccobj_(comm), use_unique_boundary_ids_(false)
#if HAVE_MPI
    , cell_comm_(comm)
#endif
{
#if HAVE_MPI
    cell_interfaces_=std::make_tuple(Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_),Interface(ccobj_));
#endif
    level_data_ptr_ = &data;
}

#if HAVE_MPI
template<class InterfaceMap>
void freeInterfaces(InterfaceMap& map)
{
    typedef typename InterfaceMap::iterator Iter;
    for(Iter i=map.begin(), end=map.end(); i!=end; ++i)
    {
        i->second.first.free();
        i->second.second.free();
    }
}

template<class InterfaceMap>
void freeInterfaces(std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                    interfaces)
{
    freeInterfaces(std::get<0>(interfaces));
    freeInterfaces(std::get<1>(interfaces));
    freeInterfaces(std::get<2>(interfaces));
    freeInterfaces(std::get<3>(interfaces));
    freeInterfaces(std::get<4>(interfaces));
}
#endif

CpGridData::~CpGridData()
{
#if HAVE_MPI
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    //freeInterfaces(face_interfaces_);
    freeInterfaces(point_interfaces_);
#endif
}

void CpGridData::populateGlobalCellIndexSet()
{
#if HAVE_MPI
    auto& cell_indexset = cellIndexSet();
    cell_indexset.beginResize();
    for (int index = 0, end = size(0); index != end ; ++index){
        cell_indexset.add(global_id_set_->id(Entity<0>(*this, EntityRep<0>(index, true))),
                          ParallelIndexSet::LocalIndex(index, AttributeSet::owner, true));
    }
    cell_indexset.endResize();
#endif
}

void CpGridData::computeUniqueBoundaryIds()
{
    // Perhaps we should make available a more comprehensive interface
    // for EntityVariable, so that we don't have to build a separate
    // vector and assign() to unique_boundary_ids_ at the end.
    int num_faces = face_to_cell_.size();
    std::vector<int> ids(num_faces, 0);
    int count = 0;
    for (int i = 0; i < num_faces; ++i) {
        cpgrid::EntityRep<1> face(i, true);
        if (face_to_cell_[face].size() == 1) {
            // It's on the boundary.
            // Important! Since boundary ids run from 1 to n,
            // we use preincrement instead of postincrement below.
            ids[i] = ++count;
        }
    }
    unique_boundary_ids_.assign(ids.begin(), ids.end());
#ifdef VERBOSE
    std::cout << "computeUniqueBoundaryIds() gave all boundary intersections\n"
              << "unique boundaryId()s ranging from 1 to " << count << std::endl;
#endif
}

int CpGridData::size(int codim) const
{
    switch (codim) {
    case 0: return cell_to_face_.size();
    case 1: return 0;
    case 2: return 0;
    case 3: return geomVector<3>().size();
    default: return 0;
    }
}

#if HAVE_MPI

// A functor that counts existent entries and renumbers them.
struct CountExistent
{
    CountExistent() : count() {}

    void operator()(int& i)
    {
        if(i < std::numeric_limits<int>::max())
            count++;
    }
    int count;
};

struct AssignAndIncrement
{
    AssignAndIncrement() : i_(){}
    void operator()(int& val){ if(val<std::numeric_limits<int>::max()) val=i_++; }
    int i_;
} assigner;

/**
 * @brief Counts the number of global ids and sets them up.
 * @param indicator A vector indicating whether an entity exists.
 * @param ids A vector to the the global ids in.
 * @param idSet The idSet of the global grid.
 * @return the number of entities that exist.
 */
template<int codim>
int setupAndCountGlobalIds(const std::vector<int>& indicator, std::vector<int>& ids,
                           const IdSet& idSet)
{
    int count = std::count_if(indicator.begin(),
                              indicator.end(),
                              [](int x) { return x < std::numeric_limits<int>::max(); });
    ids.resize(count);
    typedef typename std::vector<int>::const_iterator VIter;
    for(VIter ibegin=indicator.begin(), i=ibegin, iend= indicator.end();
        i!=iend; ++i)
    {
        if(*i<std::numeric_limits<int>::max())
            ids[*i]=idSet.id(EntityRep<codim>(i-ibegin,true));
    }
    return count;
}

template<class T>
struct GetRowType
{};

template<class T>
struct GetRowType<Opm::SparseTable<T> >
{
    typedef typename Opm::SparseTable<T>::row_type type;
};
template<class E, class A>
struct GetRowType<std::vector<E,A> >
{
    typedef typename std::vector<E,A>::value_type type;
};

PartitionType getPartitionType(const PartitionTypeIndicator& p, const EntityRep<1>& f,
                               const CpGridData&)
{
    return p.getPartitionType(f);
}

PartitionType getPartitionType(const PartitionTypeIndicator& p, int i,
                               const CpGridData& grid)
{
    return p.getPartitionType(Entity<3>(grid, i, true));
}

int getIndex(const int* i)
{
    return *i;
}

template<class T>
int getIndex(T i)
{
    return i->index();
}

template<class C>
struct DefaultContainerHandle
{
    using DataType = typename C::value_type;

    DefaultContainerHandle(const C& gatherCont, C& scatterCont)
        : gatherCont_(gatherCont), scatterCont_(scatterCont)
    {}
    bool fixedSize(std::size_t, std::size_t)
    {
        return true;
    }
    bool contains(std::size_t, std::size_t codim)
    {
        return codim == 0;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        buffer.write(gatherCont_[t.index()]);
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        buffer.read(scatterCont_[t.index()]);
    }
private:
    const C& gatherCont_;
    C& scatterCont_;
};

/// \brief Handle for face tag, normal and boundary id
struct FaceTagNormalBIdHandle
{
    using DataType = std::tuple<face_tag,FieldVector<double, 3>,int>;
    using TagContainer = EntityVariable<enum face_tag, 1>;
    using BIdContainer = EntityVariable<int, 1>;
    using NormalContainer = SignedEntityVariable<FieldVector<double, 3>, 1>;

    FaceTagNormalBIdHandle(const TagContainer& gatherTags, const NormalContainer& gatherNormals, const BIdContainer& gatherBIds,
                           TagContainer& scatterTags,  NormalContainer& scatterNormals, BIdContainer& scatterBIds)
        : gatherTags_(gatherTags), gatherNormals_(gatherNormals), gatherBIds_(gatherBIds),
          scatterTags_(scatterTags), scatterNormals_(scatterNormals), scatterBIds_(scatterBIds)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        if ( t.orientation() )
        {
            buffer.write(DataType(gatherTags_[t], gatherNormals_[t],
                                  gatherBIds_[t]));
        }
        else
        {
            auto normal = -gatherNormals_[t];
            buffer.write(DataType(gatherTags_[t], normal,
                                  gatherBIds_[t]));
        }
    }
    template<class B, class T>
    void scatter(B& buffer, T& t, std::size_t )
    {
        DataType tmp;
        buffer.read(tmp);
        scatterTags_[t]=std::get<0>(tmp);
        scatterNormals_.get(t.index()) = std::get<1>(tmp);
        scatterBIds_[t] = std::get<2>(tmp);
    }
private:
    const TagContainer& gatherTags_;
    const NormalContainer& gatherNormals_;
    const BIdContainer& gatherBIds_;
    TagContainer& scatterTags_;
    NormalContainer& scatterNormals_;
    BIdContainer& scatterBIds_;
};
/// \brief Handle for face tag, normal and boundary id
struct FaceTagNormalHandle
{
    using DataType = std::tuple<face_tag,FieldVector<double, 3> >;
    using TagContainer = EntityVariable<enum face_tag, 1>;
    using NormalContainer = SignedEntityVariable<FieldVector<double, 3>, 1>;

    FaceTagNormalHandle(const TagContainer& gatherTags, const NormalContainer& gatherNormals,
                        TagContainer& scatterTags,  NormalContainer& scatterNormals)
        : gatherTags_(gatherTags), gatherNormals_(gatherNormals),
          scatterTags_(scatterTags), scatterNormals_(scatterNormals)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 4; // 3 coordinates + 1 volume
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        if ( t.orientation() )
        {
            buffer.write(std::make_tuple(gatherTags_[t], gatherNormals_[t]));
        }
        else
        {
            auto normal = - gatherNormals_[t];
            buffer.write(std::make_tuple(gatherTags_[t], normal));
        }
    }
    template<class B, class T>
    void scatter(B& buffer, T& t, std::size_t )
    {
        DataType tmp;
        buffer.read(tmp);
        scatterTags_[t]=std::get<0>(tmp);
        scatterNormals_.get(t.index()) = std::get<1>(tmp);
    }
private:
    const TagContainer& gatherTags_;
    const NormalContainer& gatherNormals_;
    TagContainer& scatterTags_;
    NormalContainer& scatterNormals_;
};
struct PointGeometryHandle
{
    using DataType = double;
    using Container = EntityVariable<cpgrid::Geometry<0, 3>, 3>;

    PointGeometryHandle(const Container& gatherCont, Container& scatterCont)
        : gatherPoints_(gatherCont), scatterPoints_(scatterCont)
    {}
    bool fixedSize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 3;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 3;
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 3, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<3>& t)
    {
        auto& geom = gatherPoints_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 3, void>::type
    scatter(B&, const T&, std::size_t )
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void scatter(B& buffer, const EntityRep<3>& t, std::size_t )
    {
        using Vector = typename Geometry<0, 3>::GlobalCoordinate;
        Vector pos;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);
        scatterPoints_[t] = Geometry<0, 3>(pos);
    }
private:
    const Container& gatherPoints_;
    Container& scatterPoints_;
};

struct FaceGeometryHandle
{
    using DataType = double;
    using Geom = Geometry<2, 3>;
    using Container = EntityVariable<Geom, 1>;

    FaceGeometryHandle(const Container& gatherCont, Container& scatterCont)
        : gatherPoints_(gatherCont), scatterPoints_(scatterCont)
    {}
    bool fixedsize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 1;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 4; // 3 coordinates + 1 volume
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        auto& geom = gatherPoints_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
        buffer.write(geom.volume());
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        using Vector = typename Geom::GlobalCoordinate;
        Vector pos;
        double vol;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);

        buffer.read(vol);
        scatterPoints_[t] = Geom(pos, vol);
    }
private:
    const Container& gatherPoints_;
    Container& scatterPoints_;
};


struct CellGeometryHandle
{
    using DataType = double;
    using Geom = Geometry<3, 3>;
    using Container = EntityVariable<Geom, 0>;

    CellGeometryHandle(const Container& gatherCont, Container& scatterCont,
                       const std::vector<int>& gatherAquiferCells,
                       std::vector<int>& scatterAquiferCells,
                       std::shared_ptr<const EntityVariable<cpgrid::Geometry<0, 3>, 3>> pointGeom,
                       const std::vector< std::array<int,8> >& cell2Points)
        : gatherCont_(gatherCont), scatterCont_(scatterCont),
          gatherAquiferCells_(gatherAquiferCells),scatterAquiferCells_(scatterAquiferCells),
          pointGeom_(std::move(pointGeom)), cell2Points_(cell2Points)
    {}

    ~CellGeometryHandle()
    {
        std::sort(scatterAquiferCells_.begin(), scatterAquiferCells_.end());
    }

    bool fixedSize(int, int)
    {
        return true;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 5; // 3 coordinates + 1 volume + 1for indicating aquifer cells (1/0)
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<0>& t)
    {
        auto& geom = gatherCont_[t];
        for (int i = 0; i < 3; i++)
            buffer.write(geom.center()[i]);
        buffer.write(geom.volume());
        auto aquiferCell = std::lower_bound(gatherAquiferCells_.begin(),
                                            gatherAquiferCells_.end(), t.index());
        double isAquifer = (aquiferCell != gatherAquiferCells_.end() && *aquiferCell == t.index());
        buffer.write(isAquifer);
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void scatter(B& buffer, const EntityRep<0>& t, std::size_t )
    {
        using Vector = typename Geom::GlobalCoordinate;
        Vector pos;
        double vol;

        for (int i = 0; i < 3; i++)
            buffer.read(pos[i]);

        buffer.read(vol);
        scatterCont_[t] = Geom(pos, vol, pointGeom_, cell2Points_[t.index()].data());
        double isAquifer;
        buffer.read(isAquifer);
        if (isAquifer == 1.0)
        {
            scatterAquiferCells_.push_back(t.index());
        }
    }
private:
    const Container& gatherCont_;
    Container& scatterCont_;
    const std::vector<int>& gatherAquiferCells_;
    std::vector<int>& scatterAquiferCells_;
    std::shared_ptr<const EntityVariable<cpgrid::Geometry<0, 3>, 3>> pointGeom_;
    const std::vector< std::array<int,8> >& cell2Points_;
};

struct Cell2PointsDataHandle
{
    using DataType = int;
    using Vector = std::vector<std::array<int,8> >;
    Cell2PointsDataHandle(const Vector& globalCell2Points,
                          const LevelGlobalIdSet& globalIds,
                          const std::vector<std::set<int> >& globalAdditionalPointIds,
                          Vector& localCell2Points,
                          std::vector<int>& flatGlobalPoints,
                          std::vector<std::set<int> >& additionalPointIds)
        : globalCell2Points_(globalCell2Points), globalIds_(globalIds), globalAdditionalPointIds_(globalAdditionalPointIds),
          localCell2Points_(localCell2Points), flatGlobalPoints_(flatGlobalPoints), additionalPointIds_(additionalPointIds)
    {}
    bool fixedSize(int, int)
    {
        return false;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == 0;
    }
    template<class T>
    std::size_t size(const T& t)
    {
        return 8+globalAdditionalPointIds_[t.index()].size();
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        std::size_t i = t.index();
        assert(i < globalCell2Points_.size());
        const auto& points = globalCell2Points_[i];
        std::for_each(points.begin(), points.end(),
                      [&buffer, this](const int& point){
                          buffer.write(globalIds_.id(EntityRep<3>(point, true)));});
        for (const auto& point: globalAdditionalPointIds_[i])
        {
            buffer.write(point);
        }
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t s)
    {
        auto i = t.index();
        auto& points = localCell2Points_[i];
        std::for_each(points.begin(), points.end(),
                      [&buffer, this](int& point){
                          buffer.read(point);
                          this->flatGlobalPoints_.push_back(point);
                      });
        for (std::size_t p = 8; p < s; ++p)
        {
            int pi{};
            buffer.read(pi);
            this->flatGlobalPoints_.push_back(pi);
            additionalPointIds_[i].insert(pi);
        }
    }
private:
    const Vector& globalCell2Points_;
    const LevelGlobalIdSet& globalIds_;
    const std::vector<std::set<int> >& globalAdditionalPointIds_;
    Vector& localCell2Points_;
    std::vector<int>& flatGlobalPoints_;
    std::vector<std::set<int> >& additionalPointIds_;
};

template<class Table, int from>
struct RowSizeDataHandle
{
    using DataType = int;
    RowSizeDataHandle(const Table& global,
                      std::vector<int>& noEntries)
        : global_(global), noEntries_(noEntries)
    {}
    bool fixedSize(int, int)
    {
        return true;
    }
    template<class T>
    std::size_t size(const T&)
    {
        return 1;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != from, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw. Only available for other codims.");
    }
    template<class B>
    void gather(B& buffer, const EntityRep<from>& t)
    {
        buffer.write(global_.rowSize(t));
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t)
    {
        buffer.read(noEntries_[t.index()]);
    }
private:
    const Table& global_;
    std::vector<int>& noEntries_;
};

template<int from>
struct SparseTableEntity
{
    SparseTableEntity(const Opm::SparseTable<int>& table)
        : table_(table)
    {}
    int rowSize(const EntityRep<from>& index) const
    {
        return table_.rowSize(index.index());
    }
private:
    const Opm::SparseTable<int>& table_;
};

struct SparseTableDataHandle
{
    using Table = Opm::SparseTable<int>;
    using DataType = int;
    static constexpr int from = 1;
    SparseTableDataHandle(const Table& global,
                          const LevelGlobalIdSet& globalIds,
                          Table& local,
                          const std::map<int,int>& global2Local)
        : global_(global), globalIds_(globalIds), local_(local), global2Local_(global2Local)
    {}
    bool fixedSize(int, int)
    {
        return false;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class T>
    std::size_t size(const T& t)
    {
        return global_.rowSize(t.index());
    }
    template<class B, class T>
    void gather(B& buffer, const T& t)
    {
        const auto& entries = global_[t.index()];
        std::for_each(entries.begin(), entries.end(), [&buffer, this](const DataType& i){buffer.write(globalIds_.id(EntityRep<3>(i, true)));});
    }
    template<class B, class T>
    void scatter(B& buffer, const T& t, std::size_t )
    {
        const auto& entries = local_[t.index()];
        for (auto&& point : entries)
        {
            int i{};
            buffer.read(i);
            if ( point != std::numeric_limits<int>::max() )
            {
                // face already processed
                continue;
            }
            auto candidate = global2Local_.find(i);
            assert(candidate != global2Local_.end());
            point = candidate->second;
        }
    }
private:
    const Table& global_;
    const LevelGlobalIdSet& globalIds_;
    Table& local_;
    const std::map<int,int>& global2Local_;
};

template<class IdSet, int from, int to>
struct OrientedEntityTableDataHandle
{
    using DataType = int;
    using Table = OrientedEntityTable<from, to>;
    using ToEntity = typename Table::ToType;
    using FromEntity = typename Table::FromType;
    OrientedEntityTableDataHandle(const Table& global, Table& local,
                                  const IdSet* globalIds = nullptr)
        : global_(global), local_(local), globalIds_(globalIds)
    {}
    bool fixedSize(int, int)
    {
        return false;
    }
    bool contains(std::size_t dim, std::size_t codim)
    {
        return dim==3 && codim == from;
    }
    template<class T>
    typename std::enable_if<T::codimension != 0, std::size_t>::type
    size(const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
        return 0;
    }
    std::size_t size(const FromEntity& t)
    {
        return global_.rowSize(t);
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    gather(B&, const T&)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
    template<class B>
    void gather(B& buffer, const FromEntity& t)
    {
        const auto& entries = global_[t];
        if (globalIds_)
        {
            std::for_each(entries.begin(), entries.end(),
                          [&buffer, this](const ToEntity& i){
                              int id = globalIds_->id(i);
                              if (!i.orientation())
                                  id = ~id;
                              buffer.write(id);});
        }
        else
        {
            std::for_each(entries.begin(), entries.end(), [&buffer](const ToEntity& i){buffer.write(i.signedIndex());});
        }
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
protected:
    const Table& global_;
    Table& local_;
    const IdSet* globalIds_;
};

class C2FDataHandle
    : public OrientedEntityTableDataHandle<LevelGlobalIdSet,0,1>
{
public:
    C2FDataHandle(const Table& global, const LevelGlobalIdSet& globalIds, Table& local,
                  std::vector<int>& grid_size_tGlobalFaceIds)
        : OrientedEntityTableDataHandle<LevelGlobalIdSet,0,1>(global, local, &globalIds),
          grid_size_tGlobalFaceIds_(grid_size_tGlobalFaceIds)
    {}
    template<class B>
    void scatter(B& buffer, const FromEntity& t, std::size_t)
    {
        auto entries = local_.row(t);
        int i{};
        for (auto&& entry : entries)
        {
            buffer.read(i);
            if ( i < 0 )
            {
                entry = ToEntity(~i, false);
                grid_size_tGlobalFaceIds_.push_back(~i);
            }
            else
            {
                entry = ToEntity(i, true);
                grid_size_tGlobalFaceIds_.push_back(i);
            }
        }
    }
    template<class B, class T>
    typename std::enable_if<T::codimension != 0, void>::type
    scatter(B&, const T&, std::size_t)
    {
        OPM_THROW(std::logic_error, "This should never throw!");
    }
private:
    std::vector<int>& grid_size_tGlobalFaceIds_;
};

template<class IndexSet>
struct IndexSet2IdSet
{
    IndexSet2IdSet(const IndexSet& indexSet)
    {
        map_.resize(indexSet.size());
        for (const auto& entry: indexSet)
            map_[entry.local()] = entry.global();
    }
    template<class T>
    int id(const T& t) const
    {
        return map_[t.index()];
    }

    std::vector<int> map_;
};

template<class IndexSet>
class F2CDataHandle
    : public OrientedEntityTableDataHandle<IndexSet2IdSet<IndexSet>,1,0>
{
public:
    using Table = typename OrientedEntityTableDataHandle<IndexSet2IdSet<IndexSet>,1,0>::Table;
    using FromEntity = typename OrientedEntityTableDataHandle<IndexSet2IdSet<IndexSet>,1,0>::FromEntity;
    using ToEntity = typename OrientedEntityTableDataHandle<IndexSet2IdSet<IndexSet>,1,0>::ToEntity;

    F2CDataHandle(const Table& global, Table& local,
                  const IndexSet2IdSet<IndexSet>& gatherLocal2Global,
                  const IndexSet& global2Local)
        : OrientedEntityTableDataHandle<IndexSet2IdSet<IndexSet>,1,0>(global, local, &gatherLocal2Global),
          global2Local_(global2Local)
    {}
    template<class B>
    void scatter(B& buffer, const FromEntity& t, std::size_t )
    {
        auto entries = this->local_.row(t);
        int i{};
        for (auto&& entry : entries)
        {
            buffer.read(i);
            if (entry.index() != std::numeric_limits<int>::max())
            {
                // face already processed, continue to save map lookup
                continue;
            }
            bool orientation = true;
            if ( i < 0 )
            {
                i = ~i;
                orientation = false;
            }
            using IndexPair = typename IndexSet::IndexPair;
            auto candidate = std::lower_bound(global2Local_.begin(), global2Local_.end(),
                                              IndexPair(i),
                                              [](const IndexPair& p1,
                                                 const IndexPair& p2){
                                                  return p1.global() < p2.global();});
            if (candidate == global2Local_.end() || i != candidate->global())
            {
                // mark cell as being stored elsewhere
                i = std::numeric_limits<int>::max();
            }
            else
            {
                i = candidate->local();
            }
            entry = ToEntity(i, orientation);
        }
    }
private:
    const IndexSet& global2Local_;
};

template<class T>
struct AttributeDataHandle
{
    typedef std::pair<int,char> DataType;

    AttributeDataHandle(int rank, const PartitionTypeIndicator& indicator,
                        std::vector<std::map<int, char> >& vals,
                        const T& cell_to_entity,
                        const CpGridData& grid)
        : rank_(rank), indicator_(indicator), vals_(vals),
          c2e_(cell_to_entity), grid_(grid)
    {}

    bool fixedSize()
    {
        return true;
    }
    std::size_t size(std::size_t i)
    {
        return c2e_[i].size();
    }
    template<class B>
    void gather(B& buffer, std::size_t i)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f)
        {
            char t=getPartitionType(indicator_, *f, grid_);
            buffer.write(std::make_pair(rank_, t));
        }
    }

    template<class B>
    void scatter(B& buffer, std::size_t i, std::size_t s)
    {
        typedef typename GetRowType<T>::type::const_iterator RowIter;
        for(RowIter f=c2e_[i].begin(), fend=c2e_[i].end();
            f!=fend; ++f, --s)
        {
            std::pair<int,char> rank_attr;
            buffer.read(rank_attr);
            vals_[getIndex(f)].insert(rank_attr);
        }
    }
    int rank_;
    const PartitionTypeIndicator& indicator_;
    std::vector<std::map<int, char> >& vals_;
    const T& c2e_;
    const CpGridData& grid_;
};


template<class T, class Functor, class FromSet, class ToSet>
struct InterfaceFunctor
{
    InterfaceFunctor(std::map<int,std::pair<T,T> >& m)
        : map_(m)
    {}
    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        if(from.contains(mine) && to.contains(other))
            func(map_[rank].first, index);
        if(from.contains(other) && to.contains(mine))
            func(map_[rank].second, index);
    }
    std::map<int,std::pair<T,T> >& map_;
    FromSet from;
    ToSet   to;
    Functor func;
};

struct InterfaceIncrementor
{
    template<class T>
    void operator()(T& t, std::size_t)
    {
        ++t;
    }
};

struct InterfaceAdder
{
    void operator()(InterfaceInformation& info, std::size_t index)
    {
        info.add(index);
    }
};

template<class Tuple>
struct InterfaceTupleFunctor
{
    InterfaceTupleFunctor(Tuple& t)
        : t_(t)
    {}

    void operator()(int rank, std::size_t index, PartitionType mine, PartitionType other)
    {
        std::get<0>(t_)(rank, index, mine, other);
        std::get<1>(t_)(rank, index, mine, other);
        std::get<2>(t_)(rank, index, mine, other);
        std::get<3>(t_)(rank, index, mine, other);
        std::get<4>(t_)(rank, index, mine, other);
    }
    Tuple& t_;
};

struct Converter
{
    typedef EnumItem<PartitionType, InteriorEntity> Interior;
    typedef EnumItem<PartitionType, BorderEntity> Border;
    typedef EnumItem<PartitionType, OverlapEntity> Overlap;
    typedef EnumItem<PartitionType, FrontEntity> Front;
    typedef Combine<Interior,Border> InteriorBorder;
    typedef Combine<Overlap,Front> OverlapFront;

    typedef std::tuple<InteriorBorder, InteriorBorder,        Overlap     , Overlap,
                       AllSet<PartitionType> > SourceTuple;
    typedef std::tuple<InteriorBorder, AllSet<PartitionType>, OverlapFront, AllSet<PartitionType>,
                       AllSet<PartitionType> > DestinationTuple;
};

/**
 * \brief Reserves space for an interface.
 * \tparam i The index of the interface.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<std::size_t i, class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    typedef typename std::map<int,std::pair<std::size_t,std::size_t> >::const_iterator Iter;
    const std::map<int,std::pair<std::size_t,std::size_t> >& sizeMap=sizes[i];
    InterfaceMap& interfaceMap=std::get<i>(interfaces);

    for(Iter iter=sizeMap.begin(), end =sizeMap.end(); iter!=end; ++iter)
    {
        std::pair<InterfaceInformation,InterfaceInformation>& interface=interfaceMap[iter->first];
        interface.first.reserve(iter->second.first);
        interface.second.reserve(iter->second.second);
    }
}

/**
 * \brief Reserves space for the interfaces.
 * \param sizes A vector with the sizes of the interface for each neighboring rank.
 * \param interface The communication interfaces.
 */
template<class InterfaceMap>
void reserve(const std::vector<std::map<int,std::pair<std::size_t,std::size_t> > >& sizes,
             std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
             interfaces)
{
    reserve<0>(sizes, interfaces);
    reserve<1>(sizes, interfaces);
    reserve<2>(sizes, interfaces);
    reserve<3>(sizes, interfaces);
    reserve<4>(sizes, interfaces);
}

/**
 * \brief A functor that calculates the size of the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct SizeFunctor :
        public InterfaceFunctor<std::size_t, InterfaceIncrementor,
                                typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                                typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<std::size_t, InterfaceIncrementor,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
    Base;
    SizeFunctor(std::map<int,std::pair<std::size_t,std::size_t> >& m)
        :Base(m)
    {}
};

/**
 * \brief A functor that adds indices to the interface.
 * \tparam i The indentifier of the interface.
 */
template<std::size_t i>
struct AddFunctor :
        public InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                                typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                                typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
{
    typedef InterfaceFunctor<InterfaceInformation, InterfaceAdder,
                             typename std::tuple_element<i,typename Converter::SourceTuple>::type,
                             typename std::tuple_element<i,typename Converter::DestinationTuple>::type>
    Base;
    AddFunctor(std::map<int,std::pair<InterfaceInformation,InterfaceInformation> >& m)
        : Base(m)
    {}
};

class FacePartitionTypeIterator
{
public:
    FacePartitionTypeIterator(const PartitionTypeIndicator* part)
        : indicator_(part), index_()
    {}
    void operator++()
    {
        ++index_;
    }
    PartitionType operator*()
    {
        return indicator_->getFacePartitionType(index_);
    }
private:
    const PartitionTypeIndicator* indicator_;
    int index_;
};


/**
 * \brief Applies a functor the each pair of the interface.
 * \tparam Functor The type of the functor to apply.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param func The functor.
 */
template<class Functor, class T>
void iterate_over_attributes(std::vector<std::map<int,char> >& attributes,
                             T my_attribute_iter, Functor& func)
{
    typedef typename std::vector<std::map<int,char> >::const_iterator Iter;
    for(Iter begin=attributes.begin(), i=begin, end=attributes.end(); i!=end;
        ++i, ++my_attribute_iter)
    {
        typedef typename std::map<int,char>::const_iterator MIter;
        for(MIter m=i->begin(), mend=i->end(); m!=mend; ++m)
        {
            func(m->first,i-begin,PartitionType(*my_attribute_iter), PartitionType(m->second));
        }
    }
}


/**
 * \brief Creates the communication interface for either faces or points.
 * \param attributes[in] A vector that contains for each index a map from other
 * process ranks to the attribute there.
 * \param my_attributes[in] A vector with the attributes of each index on this process.
 * \param[out] interfaces The tuple with the interface maps for communication.
 */
template<class InterfaceMap,class T>
void createInterfaces(std::vector<std::map<int,char> >& attributes,
                      T partition_type_iterator,
                      std::tuple<InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap,InterfaceMap>&
                      interfaces)
{
    // calculate sizes
    std::vector<std::map<int,std::pair<std::size_t,std::size_t> > > sizes(5);
    typedef std::tuple<SizeFunctor<0>,SizeFunctor<1>,SizeFunctor<2>,SizeFunctor<3>,
                       SizeFunctor<4> > SizeTuple;

    SizeTuple
        size_functor_tuple = std::make_tuple(SizeFunctor<0>(sizes[0]),
                                             SizeFunctor<1>(sizes[1]),
                                             SizeFunctor<2>(sizes[2]),
                                             SizeFunctor<3>(sizes[3]),
                                             SizeFunctor<4>(sizes[4]));
    InterfaceTupleFunctor<SizeTuple> size_functor(size_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, size_functor);
    // reserve space
    reserve(sizes, interfaces);
    // add indices to the interface
    typedef std::tuple<AddFunctor<0>,AddFunctor<1>,AddFunctor<2>,AddFunctor<3>,
                       AddFunctor<4> > AddTuple;
    AddTuple
        add_functor_tuple(AddFunctor<0>(std::get<0>(interfaces)),
                          AddFunctor<1>(std::get<1>(interfaces)),
                          AddFunctor<2>(std::get<2>(interfaces)),
                          AddFunctor<3>(std::get<3>(interfaces)),
                          AddFunctor<4>(std::get<4>(interfaces)));
    InterfaceTupleFunctor<AddTuple> add_functor(add_functor_tuple);
    iterate_over_attributes(attributes, partition_type_iterator, add_functor);

}

void CpGridData::computeGeometry(CpGrid& grid,
                                 const DefaultGeometryPolicy&  globalGeometry,
                                 const std::vector<int>& globalAquiferCells,
                                 const OrientedEntityTable<0, 1>& globalCell2Faces,
                                 DefaultGeometryPolicy& geometry,
                                 std::vector<int>& aquiferCells,
                                 const OrientedEntityTable<0, 1>& cell2Faces,
                                 const std::vector< std::array<int,8> >& cell2Points)
{
    FaceGeometryHandle faceGeomHandle(*globalGeometry.geomVector(std::integral_constant<int,1>()),
                                      *geometry.geomVector(std::integral_constant<int,1>()));
    FaceViaCellHandleWrapper<FaceGeometryHandle>
        wrappedFaceGeomHandle(faceGeomHandle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedFaceGeomHandle);

    PointGeometryHandle pointGeomHandle(*globalGeometry.geomVector(std::integral_constant<int,3>()),
                                        *geometry.geomVector(std::integral_constant<int,3>()));
    grid.scatterData(pointGeomHandle);

    CellGeometryHandle cellGeomHandle(*globalGeometry.geomVector(std::integral_constant<int,0>()),
                                      *geometry.geomVector(std::integral_constant<int,0>()),
                                      globalAquiferCells, aquiferCells,
                                      geometry.geomVector(std::integral_constant<int,3>()),
                                      cell2Points);
    grid.scatterData(cellGeomHandle);
}

void computeFace2Point(CpGrid& grid,
                       const OrientedEntityTable<0, 1>& globalCell2Faces,
                       const LevelGlobalIdSet& globalIds,
                       const OrientedEntityTable<0, 1>& cell2Faces,
                       const Opm::SparseTable<int>& globalFace2Points,
                       Opm::SparseTable<int>& face2Points,
                       const std::map<int,int>& global2local,
                       std::size_t noFaces)
{
    std::vector<int> rowSizes(noFaces);
    using EntityTable = SparseTableEntity<1>;
    using RowSizeDataHandle = RowSizeDataHandle<EntityTable, 1>;
    EntityTable wrappedGlobal(globalFace2Points);
    RowSizeDataHandle rowSizeHandle(wrappedGlobal, rowSizes);
    FaceViaCellHandleWrapper<RowSizeDataHandle>
        wrappedSizeHandle(rowSizeHandle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedSizeHandle);
    face2Points.allocate(rowSizes.begin(), rowSizes.end());
    // Use entity with index INT_MAX to mark unprocessed row entries
    for (int row = 0, size = face2Points.size(); row < size; ++row)
    {
        for (auto&& point : face2Points[row])
        {
            point = std::numeric_limits<int>::max();
        }
    }
    SparseTableDataHandle handle(globalFace2Points, globalIds, face2Points, global2local);
    FaceViaCellHandleWrapper<SparseTableDataHandle>
        wrappedHandle(handle, globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedHandle);
}

template<class IndexSet>
void computeFace2Cell(CpGrid& grid,
                      const OrientedEntityTable<0, 1>& globalCell2Faces,
                      const OrientedEntityTable<0, 1>& cell2Faces,
                      const OrientedEntityTable<1, 0>& globalFace2Cells,
                      OrientedEntityTable<1, 0>& face2Cells,
                      const IndexSet& global2local,
                      const IndexSet& globalIndexSet,
                      std::size_t noFaces)
{
    std::vector<int> rowSizes(noFaces);
    using Table = OrientedEntityTable<1,0>;
    RowSizeDataHandle<Table,1> rowSizeHandle(globalFace2Cells, rowSizes);
    FaceViaCellHandleWrapper<RowSizeDataHandle<Table,1> > wrappedSizeHandle(rowSizeHandle,
                                                                            globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedSizeHandle);
    face2Cells.allocate(rowSizes.begin(), rowSizes.end());
    // Use entity with index INT_MAX to mark unprocessed row entries
    for (int row = 0, size = face2Cells.size(); row < size; ++row)
    {
        for (auto&& face : face2Cells.row(EntityRep<1>(row, true)))
        {
            face = EntityRep<0>(std::numeric_limits<int>::max(), true);
        }
    }
    IndexSet2IdSet<IndexSet> local2Global(globalIndexSet);
    F2CDataHandle<IndexSet> entryHandle(globalFace2Cells, face2Cells, local2Global, global2local);
    FaceViaCellHandleWrapper<F2CDataHandle<IndexSet> > wrappedEntryHandle(entryHandle,
                                                                          globalCell2Faces, cell2Faces);
    grid.scatterData(wrappedEntryHandle);
#ifndef NDEBUG
    for (int row = 0, size = face2Cells.size(); row < size; ++row)
    {
        bool oneValid = false;
        for (auto&& face : face2Cells.row(EntityRep<1>(row, true)))
        {
            oneValid = oneValid || face.index() != std::numeric_limits<int>::max();
        }
        assert(oneValid);
    }
#endif
}


std::map<int,int> computeCell2Face(CpGrid& grid,
                                   const OrientedEntityTable<0, 1>& globalCell2Faces,
                                   const LevelGlobalIdSet& globalIds,
                                   OrientedEntityTable<0, 1>& cell2Faces,
                                   std::vector<int>& map2Global,
                                   std::size_t noCells)
{
    std::vector<int> rowSizes(noCells);
    using Table = OrientedEntityTable<0,1>;
    RowSizeDataHandle<Table,0> rowSizeHandle(globalCell2Faces, rowSizes);
    grid.scatterData(rowSizeHandle);
    cell2Faces.allocate(rowSizes.begin(), rowSizes.end());
    map2Global.reserve((noCells*6)*1.1);
    C2FDataHandle handle(globalCell2Faces, globalIds, cell2Faces,
                         map2Global);
    grid.scatterData(handle);
    // make map2Global a map from local index to global id
    std::sort(map2Global.begin(),map2Global.end());
    auto newEnd = std::unique(map2Global.begin(),map2Global.end());
    map2Global.resize(newEnd - map2Global.begin());
    // Convert face ids to local ones
    std::map<int, int> map2Local;
    auto current_global = map2Global.begin();
    int localId = 0;
    // \todo improve since we are inserting values by sorted keys.
    std::generate_n(std::inserter(map2Local, map2Local.begin()),
                    map2Global.size(),
                    [&localId, &current_global](){ return std::make_pair(*(current_global++), localId++); });
    // translate global to local ids
    for (int row = 0, size = cell2Faces.size(); row < size; ++row)
    {
        for (auto&& face : cell2Faces.row(EntityRep<0>(row, true)))
        {
            auto index = face.signedIndex();
            if (index < 0 )
            {
                face = EntityRep<1>(map2Local[~index], false);
            }
            else
            {
                face = EntityRep<1>(map2Local[index], true);
            }
        }
    }
    return map2Local;
}

std::vector<std::set<int> > computeAdditionalFacePoints(const std::vector<std::array<int,8> >& globalCell2Points,
                                                        const OrientedEntityTable<0, 1>& globalCell2Faces,
                                                        const Opm::SparseTable<int>& globalFace2Points,
                                                        const LevelGlobalIdSet& globalIds)
{
    std::vector<std::set<int> > additionalFacePoints(globalCell2Points.size());

    for ( std::size_t c = 0; c < globalCell2Points.size(); ++c)
    {
        const auto& points = globalCell2Points[c];
        for(const auto& face: globalCell2Faces[EntityRep<0>(c, true)])
            for(const auto& point: globalFace2Points[face.index()])
            {
                auto candidate = std::find(points.begin(), points.end(), point);
                if(candidate == points.end())
                    // point is not a corner of the cell
                    additionalFacePoints[c].insert(globalIds.id(EntityRep<3>(point,true)));
            }
    }
    return additionalFacePoints;
}

template<bool send, class Map2Global, class Map2Local>
void createInterfaceList(const typename CpGridData::InterfaceMap::value_type& procCellLists,
                         const std::vector<std::array<int,8> >& cell2Points,
                         const std::vector<std::set<int> >& additionalPoints,
                         const Map2Global& local2Global,
                         Map2Local& map2Local,
                         typename CpGridData::InterfaceMap::mapped_type& pointLists
                         )
{
    const auto& cellList = send? procCellLists.second.first : procCellLists.second.second;

    // Create list of global ids from cell lists
    std::vector<int> tmpPoints;
    std::size_t noAdditional{};
    for (auto const& addPoints: additionalPoints)
        noAdditional += addPoints.size();

    tmpPoints.reserve(cell2Points.size()*8+noAdditional);

    for (std::size_t c = 0; c < cellList.size(); ++c)
    {
        for (const auto& point: cell2Points[cellList[c]])
            tmpPoints.push_back(local2Global(point));

        for (const auto& point: additionalPoints[cellList[c]])
            tmpPoints.push_back(point); // is already a global id
    }

    // sort list and make it unique
    std::sort(tmpPoints.begin(), tmpPoints.end());
    auto newEnd = std::unique(tmpPoints.begin(), tmpPoints.end());
    tmpPoints.resize(newEnd - tmpPoints.begin());

    // map entries to local indices
    for ( auto&& point: tmpPoints)
        point = map2Local[point];

    auto& pointList = send ? pointLists.first : pointLists.second;
    pointList.reserve(tmpPoints.size());

    for ( auto && point: tmpPoints)
        pointList.add(point);
}

std::map<int,int> computeCell2Point(CpGrid& grid,
                                    const std::vector<std::array<int,8> >& globalCell2Points,
                                    const LevelGlobalIdSet& globalIds,
                                    const OrientedEntityTable<0, 1>& globalCell2Faces,
                                    const Opm::SparseTable<int>& globalFace2Points,
                                    std::vector<std::array<int,8> >& cell2Points,
                                    std::vector<int>& map2Global,
                                    std::size_t noCells,
                                    const typename CpGridData::InterfaceMap& cellInterfaces,
                                    typename CpGridData::InterfaceMap& pointInterfaces
                                    )
{
    cell2Points.resize(noCells);
    map2Global.reserve(noCells*8*1.1);
    auto globalAdditionalPoints = computeAdditionalFacePoints(globalCell2Points, globalCell2Faces,
                                                              globalFace2Points,
                                                              globalIds);
    std::vector<std::set<int> > additionalPoints(noCells);
    Cell2PointsDataHandle handle(globalCell2Points, globalIds,
                                 globalAdditionalPoints,
                                 cell2Points,
                                 map2Global,
                                 additionalPoints);
    grid.scatterData(handle);
    // make map2Global a map from local index to global id
    std::sort(map2Global.begin(),map2Global.end());
    auto newEnd = std::unique(map2Global.begin(),map2Global.end());
    map2Global.resize(newEnd - map2Global.begin());
    // Convert point ids to local ones
    std::map<int, int> map2Local;
    auto current_global = map2Global.begin();
    int localId = 0;
    // \todo improve since we are inserting values by sorted keys.
    std::generate_n(std::inserter(map2Local, map2Local.begin()),
                    map2Global.size(),
                    [&localId, &current_global](){ return std::make_pair(*(current_global++), localId++); });
    for (auto&& points : cell2Points)
    {
        for (auto&& point : points)
        {
            point = map2Local[point];
        }
    }

    // Create interfaces for point communication

    for ( const auto& procCellLists: cellInterfaces)
    {
        // The send list
        ReversePointGlobalIdSet globalMap2Local(globalIds);
        createInterfaceList<true>(procCellLists, globalCell2Points,
                                  globalAdditionalPoints,
                                  [&globalIds](int i){
                                      return globalIds.id(EntityRep<3>(i, true));
                                  },
                                  globalMap2Local,
                                  pointInterfaces[procCellLists.first]);
        globalMap2Local.release();
        // The receive list
        createInterfaceList<false>(procCellLists, cell2Points,
                                   additionalPoints,
                                   [&map2Global](int i)
                                   {
                                       return map2Global[i];
                                   },
                                   map2Local,
                                   pointInterfaces[procCellLists.first]);
    }
    return map2Local;
}

#endif // #if HAVE_MPI

void CpGridData::distributeGlobalGrid(CpGrid& grid,
                                      const CpGridData& view_data,
                                      const std::vector<int>& /* cell_part */)
{
#if HAVE_MPI
    auto& cell_indexset = cellIndexSet();
    auto& cell_remote_indices = cellRemoteIndices();
    // setup the remote indices.
    cell_remote_indices.template rebuild<false>(); // We could probably also compute this on our own, like before?

    // We can identify existing cells with the help of the index set.
    // Now we need to compute the existing faces and points. Either exist
    // if they are reachable from an existing cell.
    // We use std::numeric_limits<int>::max() to indicate non-existent entities.
    std::vector<int> map2GlobalFaceId;
    std::vector<int> map2GlobalPointId;
    std::map<int,int> point_indicator =
        computeCell2Point(grid, view_data.cell_to_point_, *view_data.global_id_set_, view_data.cell_to_face_,
                          view_data.face_to_point_, cell_to_point_,
                          map2GlobalPointId, cell_indexset.size(),
                          *grid.cell_scatter_gather_interfaces_,
                          *grid.point_scatter_gather_interfaces_);

    // create global ids array for cells. The parallel index set uses the global id
    // as the global index.
    std::vector<int> map2GlobalCellId(cell_indexset.size());
    for(const auto& i: cell_indexset)
    {
        map2GlobalCellId[i.local()]=i.global();
    }

    std::map<int,int> face_indicator =
        computeCell2Face(grid, view_data.cell_to_face_, *view_data.global_id_set_, cell_to_face_,
                         map2GlobalFaceId, cell_indexset.size());

    auto noExistingPoints = map2GlobalPointId.size();
    auto noExistingFaces = map2GlobalFaceId.size();

    global_id_set_->swap(map2GlobalCellId, map2GlobalFaceId, map2GlobalPointId);

    computeFace2Cell(grid, view_data.cell_to_face_, cell_to_face_,
                     view_data.face_to_cell_, face_to_cell_, cell_indexset, view_data.cellIndexSet(), noExistingFaces);
    computeFace2Point(grid,  view_data.cell_to_face_, *view_data.global_id_set_, cell_to_face_,
                      view_data.face_to_point_, face_to_point_, point_indicator,
                      noExistingFaces);


    logical_cartesian_size_=view_data.logical_cartesian_size_;

    // Set up the new topology arrays
    geometry_.geomVector(std::integral_constant<int,1>()) -> resize(noExistingFaces);
    geometry_.geomVector(std::integral_constant<int,0>()) -> resize(cell_to_face_.size());
    geometry_.geomVector(std::integral_constant<int,3>()) -> resize(noExistingPoints);

    computeGeometry(grid, view_data.geometry_, view_data.aquifer_cells_, view_data.cell_to_face_,
                    geometry_, aquifer_cells_, cell_to_face_, cell_to_point_);

    global_cell_.resize(cell_indexset.size());

    // communicate global cell
    DefaultContainerHandle<std::vector<int> > indexHandle(view_data.global_cell_, global_cell_);
    grid.scatterData(indexHandle);

    // Scatter face tags, normals, and boundary ids.
    auto noBids = view_data.unique_boundary_ids_.size();
    bool hasBids = ccobj_.max(noBids);
    face_tag_.resize(noExistingFaces);
    face_normals_.resize(noExistingFaces);

    if (hasBids)
    {
        unique_boundary_ids_.resize(noExistingFaces);
        FaceTagNormalBIdHandle faceHandle(view_data.face_tag_, view_data.face_normals_, view_data.unique_boundary_ids_,
                                          face_tag_, face_normals_, unique_boundary_ids_);
        FaceViaCellHandleWrapper<FaceTagNormalBIdHandle>
            wrappedFaceHandle(faceHandle, view_data.cell_to_face_, cell_to_face_);
        grid.scatterData(wrappedFaceHandle);
    }
    else
    {
        FaceTagNormalHandle faceHandle(view_data.face_tag_, view_data.face_normals_,
                                       face_tag_, face_normals_);
        FaceViaCellHandleWrapper<FaceTagNormalHandle>
            wrappedFaceHandle(faceHandle, view_data.cell_to_face_, cell_to_face_);
        grid.scatterData(wrappedFaceHandle);
    }

    // Compute the partition type for cell
    computeCellPartitionType();

    // Compute partition type for points
    computePointPartitionType();
   
    computeCommunicationInterfaces(noExistingPoints);   
#else // #if HAVE_MPI
    static_cast<void>(grid);
    static_cast<void>(view_data);
#endif
}

void CpGridData::computeCellPartitionType()
{
#if HAVE_MPI
    // Compute the partition type for cell
    auto& cell_indexset = cellIndexSet();
    partition_type_indicator_->cell_indicator_.resize(cell_indexset.size());
    for(const auto& i: cell_indexset)
    {
        partition_type_indicator_->cell_indicator_[i.local()]=
            i.local().attribute()==AttributeSet::owner?
            InteriorEntity:OverlapEntity;
    }
#endif
}

void CpGridData::computePointPartitionType()
{
#if HAVE_MPI
    // We initialize all points with interior. Then we loop over the faces. If a face is of
    // type border, then the type of the point is overwritten with border. In the other cases
    // we set the type of the point to the one of the face as long as the type of the point is
    // not border.
    partition_type_indicator_->point_indicator_.resize(geometry_.geomVector<3>().size(),
                                                       InteriorEntity);
    for(int i=0; i<face_to_point_.size(); ++i)
    {
        for(auto p=face_to_point_[i].begin(),
                pend=face_to_point_[i].end(); p!=pend; ++p)
        {
            PartitionType new_type=partition_type_indicator_->getFacePartitionType(i);
            PartitionType old_type=PartitionType(partition_type_indicator_->point_indicator_[*p]);
            if(old_type==InteriorEntity)
            {
                if(new_type!=OverlapEntity)
                    partition_type_indicator_->point_indicator_[*p]=new_type;
            }
            if(old_type==OverlapEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
            if(old_type==FrontEntity && new_type==BorderEntity)
                partition_type_indicator_->point_indicator_[*p]=new_type;
        }
    }
#endif
}

void CpGridData::computeCommunicationInterfaces([[maybe_unused]] int noExistingPoints)
{
#if HAVE_MPI
    // Compute the interface information for cells
    std::get<InteriorBorder_All_Interface>(cell_interfaces_)
        .build(cellRemoteIndices(), EnumItem<AttributeSet, AttributeSet::owner>(),
               AllSet<AttributeSet>());
    std::get<Overlap_OverlapFront_Interface>(cell_interfaces_)
        .build(cellRemoteIndices(), EnumItem<AttributeSet, AttributeSet::copy>(),
               EnumItem<AttributeSet, AttributeSet::copy>());
    std::get<Overlap_All_Interface>(cell_interfaces_)
        .build(cellRemoteIndices(), EnumItem<AttributeSet, AttributeSet::copy>(),
               AllSet<AttributeSet>());
    std::get<All_All_Interface>(cell_interfaces_)
        .build(cellRemoteIndices(), AllSet<AttributeSet>(), AllSet<AttributeSet>());

    // Now we use the all_all communication of the cells to compute which faces and points
    // are also present on other processes and with what attribute.
    const auto& all_all_cell_interface = std::get<All_All_Interface>(cell_interfaces_);

    Communicator comm(all_all_cell_interface.communicator(),
                      all_all_cell_interface.interfaces());

    /*
    // code deactivated, because users cannot access face indices and therefore
    // communication on faces makes no sense!
    std::vector<std::map<int,char> > face_attributes(noExistingFaces);
    AttributeDataHandle<Opm::SparseTable<EntityRep<1> > >
    face_handle(ccobj_.rank(), *partition_type_indicator_,
    face_attributes, static_cast<Opm::SparseTable<EntityRep<1> >&>(cell_to_face_),
    *this);
    if( std::get<All_All_Interface>(cell_interfaces_).interfaces().size() )
    {
    comm.forward(face_handle);
    }
    createInterfaces(face_attributes, FacePartitionTypeIterator(partition_type_indicator_),
    face_interfaces_);
    std::vector<std::map<int,char> >().swap(face_attributes);
    */
    std::vector<std::map<int,char> > point_attributes(noExistingPoints);
    AttributeDataHandle<std::vector<std::array<int,8> > >
        point_handle(ccobj_.rank(), *partition_type_indicator_,
                     point_attributes, cell_to_point_, *this);
    if( static_cast<const Dune::Interface&>(std::get<All_All_Interface>(cell_interfaces_))
        .interfaces().size() )
    {
        comm.forward(point_handle);
    }
    createInterfaces(point_attributes, partition_type_indicator_->point_indicator_.begin(),
                     point_interfaces_);
#endif
}

std::array<Dune::FieldVector<double,3>,8> CpGridData::getReferenceRefinedCorners(int idx_in_parent_cell, const std::array<int,3>& cells_per_dim) const
{
    // Refined cells in parent cell: k*cells_per_dim[0]*cells_per_dim[1] + j*cells_per_dim[0] + i
    std::array<int,3> ijk = getIJK(idx_in_parent_cell, cells_per_dim);

    std::array<Dune::FieldVector<double,3>,8> corners_in_parent_reference_elem = { // corner '0'
        {{ double(ijk[0])/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '1'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '2'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '3'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2])/cells_per_dim[2] },
         // corner '4'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '5'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1])/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '6'
         { double(ijk[0])/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] },
         // corner '7'
         { double(ijk[0]+1)/cells_per_dim[0], double(ijk[1]+1)/cells_per_dim[1], double(ijk[2]+1)/cells_per_dim[2] }
        }
    };
    return corners_in_parent_reference_elem;
}

std::array<int,3> CpGridData::getPatchDim(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    return {endIJK[0]-startIJK[0], endIJK[1]-startIJK[1], endIJK[2]-startIJK[2]};
}

std::vector<int> CpGridData::getPatchCorners(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
    const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
    // Get grid dimension (total cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    /// PATCH CORNERS
    std::vector<int> patch_corners;
    patch_corners.reserve((patch_dim[0]+1)*(patch_dim[1]+1)*(patch_dim[2]+1));
    for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                patch_corners.push_back((j*(grid_dim[0]+1)*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+k);
            } // end i-for-loop
        } // end j-for-loop
    } // end k-for-loop
    return patch_corners;
}

std::vector<int> CpGridData::getPatchFaces(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
    const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
    // Get grid dimension (total cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    /// PATCH FACES
    std::vector<int> patch_faces;
    patch_faces.reserve(((patch_dim[0]+1)*patch_dim[1]*patch_dim[2])     // i_patch_faces
                        + (patch_dim[0]*(patch_dim[1]+1)*patch_dim[2])   // j_patch_faces
                        + (patch_dim[0]*patch_dim[1]*(patch_dim[2]+1))); // k_patch_faces
    int face_idx;
    // I_FACES
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                face_idx = (j*(grid_dim[0]+1)*grid_dim[2]) +(i*grid_dim[2]) + k;
                patch_faces.push_back(face_idx);
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // J_FACES
    for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                face_idx = ((grid_dim[0]+1)*grid_dim[1]*grid_dim[2]) // i_grid_faces
                    + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2]) + k;
                patch_faces.push_back(face_idx);
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // K_FACES
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                face_idx = (grid_dim[0]*(grid_dim[1]+1)*grid_dim[2]) //j_grid_faces
                    + ((grid_dim[0]+1)*grid_dim[1]*grid_dim[2])          // i_grid_faces
                    + (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ k;
                patch_faces.push_back(face_idx);
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    return patch_faces;
}

std::vector<int> CpGridData::getPatchCells(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
    const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
    // Get grid dimension (total cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    std::vector<int> patch_cells;
    patch_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
    /// PATCH CELLS
    for (int k = startIJK[2]; k < endIJK[2]; ++k) {
        for (int j = startIJK[1]; j < endIJK[1]; ++j) {
            for (int i = startIJK[0]; i < endIJK[0]; ++i) {
                patch_cells.push_back((k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) +i);
            } // end i-for-loop
        } // end j-for-loop
    } // end k-for-loop
    return patch_cells;
}

void CpGridData::checkCuboidShape(const std::vector<int>& cellIdx_vec) const
{
    bool cuboidShape = true;
    for (const auto cellIdx : cellIdx_vec)
    {
        const auto cellToPoint = cell_to_point_[cellIdx]; // bottom face corners {0,1,2,3}, top face corners {4,5,6,7}
        // Compute 'cuboid' volume with corners: |corn[1]-corn[0]|x|corn[3]-corn[1]|x|corn[5]-corn[1]|
        std::vector<cpgrid::Geometry<0,3>::GlobalCoordinate> aFewCorners;
        aFewCorners.resize(4); // {'0', '1', '3', '5'}
        aFewCorners[0] = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[0]).center();
        aFewCorners[1] = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[1]).center();
        aFewCorners[2] = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[3]).center();
        aFewCorners[3] = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[5]).center();
        //  l = length. b = breadth. h = height.
        double  length, breadth, height;
        length = std::sqrt( ((aFewCorners[1][0] -aFewCorners[0][0])*(aFewCorners[1][0] -aFewCorners[0][0])) +
                            ((aFewCorners[1][1] -aFewCorners[0][1])*(aFewCorners[1][1] -aFewCorners[0][1])) +
                            ((aFewCorners[1][2] -aFewCorners[0][2])*(aFewCorners[1][2] -aFewCorners[0][2])));
        breadth = std::sqrt( ((aFewCorners[1][0] -aFewCorners[2][0])*(aFewCorners[1][0] -aFewCorners[2][0])) +
                             ((aFewCorners[1][1] -aFewCorners[2][1])*(aFewCorners[1][1] -aFewCorners[2][1])) +
                             ((aFewCorners[1][2] -aFewCorners[2][2])*(aFewCorners[1][2] -aFewCorners[2][2])));
        height = std::sqrt( ((aFewCorners[1][0] -aFewCorners[3][0])*(aFewCorners[1][0] -aFewCorners[3][0])) +
                            ((aFewCorners[1][1] -aFewCorners[3][1])*(aFewCorners[1][1] -aFewCorners[3][1])) +
                            ((aFewCorners[1][2] -aFewCorners[3][2])*(aFewCorners[1][2] -aFewCorners[3][2])));
        const double cuboidVolume = length*breadth*height;
        const auto cellVolume =  (*(this -> geometry_.geomVector(std::integral_constant<int,0>())))[EntityRep<0>(cellIdx, true)].volume();
        cuboidShape = cuboidShape && (std::abs(cuboidVolume - cellVolume) <  1e-6);
        if (!cuboidShape){
            OPM_THROW(std::logic_error, "At least one cell has no cuboid shape. Its refinement is not supported yet.\n");
        }
    }
}

std::array<std::vector<double>,3> CpGridData::getWidthsLengthsHeights(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    std::vector<double> widthsX;
    widthsX.reserve(endIJK[0] - startIJK[0]);
    std::vector<double> lengthsY;
    lengthsY.reserve(endIJK[1] - startIJK[1]);
    std::vector<double> heightsZ;
    heightsZ.reserve(endIJK[2] - startIJK[2]);

    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();

    for (int i = startIJK[0]; i < endIJK[0]; ++i) {
        int cellIdx = (startIJK[2]*grid_dim[0]*grid_dim[1]) + (startIJK[1]*grid_dim[0]) + i;
        const auto cellToPoint = cell_to_point_[cellIdx]; // bottom face corners {0,1,2,3}, top face corners {4,5,6,7}
        // x = |corn[1]-corn[0]|
        // Compute difference and dot using DUNE functionality
        auto difference = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[0]).center();
        difference -= (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[1]).center();
        auto x = difference.two_norm();
        widthsX.push_back(x);
    }
    for (int j = startIJK[1]; j < endIJK[1]; ++j)
    {
        int cellIdx = (startIJK[2]*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) + startIJK[0];
        const auto cellToPoint = cell_to_point_[cellIdx]; // bottom face corners {0,1,2,3}, top face corners {4,5,6,7}
        // y = |corn[3]-corn[1]|
        // Compute difference and dot using DUNE functionality
        auto difference = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[3]).center();
        difference -= (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[1]).center();
        auto y = difference.two_norm();
        lengthsY.push_back(y);
    }
    for (int k = startIJK[2]; k < endIJK[2]; ++k)
    {
        int cellIdx = (k*grid_dim[0]*grid_dim[1]) + (startIJK[1]*grid_dim[0]) + startIJK[0];
        const auto cellToPoint = cell_to_point_[cellIdx]; // bottom face corners {0,1,2,3}, top face corners {4,5,6,7}
        // z = |corn[4]-corn[0]|
        // Compute difference and dot using DUNE functionality
        auto difference = (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[4]).center();
        difference -= (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellToPoint[0]).center();
        auto z = difference.two_norm();
        heightsZ.push_back(z);
    }
    return {widthsX, lengthsY, heightsZ};
}

std::vector<int> CpGridData::getPatchBoundaryCorners(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
    const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
    // Get grid dimension (total cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    /// PATCH BOUNDARY CORNERS
    std::vector<int> patch_boundary_corners;
    patch_boundary_corners.reserve(((patch_dim[0]+1)*(patch_dim[2]+1)*2) + ((patch_dim[1]-1)*(patch_dim[2]+1)*2)
                                   + ((patch_dim[0]-1)*(patch_dim[1]-1)*2));
    for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                if ( (j == startIJK[1]) || (j == endIJK[1])
                     ||  (i == startIJK[0]) || (i == endIJK[0])
                     ||  (k == startIJK[2]) || (k == endIJK[2])) {
                    patch_boundary_corners.push_back((j*(grid_dim[0]+1)*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+k);
                }
            } // end i-for-loop
        } // end j-for-loop
    } // end k-for-loop
    return patch_boundary_corners;
}

std::array<std::vector<int>,6> CpGridData::getBoundaryPatchFaces(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK) const
{
    // Get the patch dimension (total cells in each direction). Used to 'reserve vectors'.
    const std::array<int,3>& patch_dim = getPatchDim(startIJK, endIJK);
    // Get grid dimension (total cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    // Auxiliary integers to simplify notation.
    const int& i_grid_faces =  (grid_dim[0]+1)*grid_dim[1]*grid_dim[2];
    const int& j_grid_faces =  grid_dim[0]*(grid_dim[1]+1)*grid_dim[2];

    std::array<std::vector<int>,6> boundary_patch_faces;
    // { I_FACE false vector, I_FACE true vector, J_FACE false vector, J_FACE true vector, K_FACE false vector, K_FACE true vector}
    boundary_patch_faces[0].reserve(2*patch_dim[1]*patch_dim[2]); // I_FACE false vector
    boundary_patch_faces[1].reserve(2*patch_dim[1]*patch_dim[2]); // I_FACE true vector
    boundary_patch_faces[2].reserve(2*patch_dim[0]*patch_dim[2]); // J_FACE false vector (front)
    boundary_patch_faces[3].reserve(2*patch_dim[0]*patch_dim[2]); // J_FACE true vector  (back)
    boundary_patch_faces[4].reserve(2*patch_dim[0]*patch_dim[1]); // K_FACE false vector (bottom)
    boundary_patch_faces[5].reserve(2*patch_dim[0]*patch_dim[1]); // K_FACE true vector  (top)
    // Boundary I_FACE faces
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int k = startIJK[2]; k < endIJK[2]; ++k) {
            boundary_patch_faces[0].push_back( (j*(grid_dim[0]+1)*grid_dim[2]) + (startIJK[0]*grid_dim[2])+ k); // I_FACE false
            boundary_patch_faces[1].push_back( (j*(grid_dim[0]+1)*grid_dim[2]) + (endIJK[0]*grid_dim[2])+ k); // I_FACE true
        }
    }
    // Boundary J_FACE faces
    for (int i = startIJK[0]; i < endIJK[0]; ++i) {
        for (int k = startIJK[2]; k < endIJK[2]; ++k) {
            boundary_patch_faces[2].push_back(i_grid_faces + (startIJK[1]*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k); // J_FACE false
            boundary_patch_faces[3].push_back(i_grid_faces + (endIJK[1]*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k); // J_FACE true
        }
    }
    // Boundary K_FACE faces
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]; ++i) {
            boundary_patch_faces[4].push_back( i_grid_faces + j_grid_faces +
                                               (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ startIJK[2] ); // K_FACE false
            boundary_patch_faces[5].push_back( i_grid_faces + j_grid_faces +
                                               (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ endIJK[2]); // K_FACE true
        }
    }
    return boundary_patch_faces;
}

bool CpGridData::disjointPatches(const std::vector<std::array<int,3>>& startIJK_vec,
                                 const std::vector<std::array<int,3>>& endIJK_vec) const
{
    assert(!startIJK_vec.empty());
    assert(!endIJK_vec.empty());
    if ((startIJK_vec.size() == 1) && (endIJK_vec.size() == 1)){
        return true;
    }
    if (startIJK_vec.size() != endIJK_vec.size() ){
        OPM_THROW(std::logic_error, "Sizes of the arguments differ. Not enough information provided.");
    }
    for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch){
        bool valid_patch = true;
        for (int c = 0; c < 3; ++c){
            valid_patch = valid_patch && (startIJK_vec[patch][c] < endIJK_vec[patch][c]);
        }
        if (!valid_patch){
            OPM_THROW(std::logic_error, "There is at least one invalid block of cells.");
        }
    }
    bool are_disjoint = true;
    for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch) {
        bool patch_disjoint_with_otherPatches = true;
        for (grid_size_t other_patch = patch+1; other_patch < startIJK_vec.size(); ++other_patch) {
            bool otherPatch_on_rightOrLeft_of_patch = (startIJK_vec[other_patch][0] > endIJK_vec[patch][0]) ||
                (endIJK_vec[other_patch][0] < startIJK_vec[patch][0]);
            if (!otherPatch_on_rightOrLeft_of_patch) {
                bool otherPatch_on_frontOrBack_of_patch = (startIJK_vec[other_patch][1] > endIJK_vec[patch][1]) ||
                    (endIJK_vec[other_patch][1] < startIJK_vec[patch][1]);
                if (!otherPatch_on_frontOrBack_of_patch) {
                    bool otherPatch_on_topOrBottom_of_patch = (startIJK_vec[other_patch][2] > endIJK_vec[patch][2]) ||
                        (endIJK_vec[other_patch][2] < startIJK_vec[patch][2]);
                    patch_disjoint_with_otherPatches = patch_disjoint_with_otherPatches && otherPatch_on_topOrBottom_of_patch;
                    // true for disjoint patches, false for overlapping ones.
                }
            }
            else{
                patch_disjoint_with_otherPatches = patch_disjoint_with_otherPatches && otherPatch_on_rightOrLeft_of_patch;
                // true for disjoint patches
            }
            if (!patch_disjoint_with_otherPatches){
                return patch_disjoint_with_otherPatches; // should be false
            }
        }
        are_disjoint = are_disjoint && patch_disjoint_with_otherPatches;
    }
    return are_disjoint; // should be true
}

bool CpGridData::patchesShareFace(const std::vector<std::array<int,3>>& startIJK_vec,
                                  const std::vector<std::array<int,3>>& endIJK_vec) const
{
    assert(!startIJK_vec.empty());
    assert(!endIJK_vec.empty());
    if ((startIJK_vec.size() == 1) && (endIJK_vec.size() == 1)){
        return false;
    }
    if (startIJK_vec.size() != endIJK_vec.size() ){
        OPM_THROW(std::logic_error, "Sizes of the arguments differ. Not enough information provided.");
    }
    for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch){
        bool valid_patch = true;
        for (int c = 0; c < 3; ++c){
            valid_patch = valid_patch && (startIJK_vec[patch][c] < endIJK_vec[patch][c]);
        }
        if (!valid_patch){
            OPM_THROW(std::logic_error, "There is at least one invalid block of cells.");
        }
    }

    const auto& detectSharing = [](std::vector<int> faceIdxs, std::vector<int> otherFaceIdxs){
        bool faceIsShared = false;
        for (const auto& face : faceIdxs) {
            for (const auto& otherFace : otherFaceIdxs) {
                faceIsShared = faceIsShared || (face == otherFace);
                if (faceIsShared) {
                    return faceIsShared; // should be true here
                }
            }
        }
        return faceIsShared; // should be false here
    };

    for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch) {
        const auto& [iFalse, iTrue, jFalse, jTrue, kFalse, kTrue] = this->getBoundaryPatchFaces(startIJK_vec[patch], endIJK_vec[patch]);
        for (grid_size_t other_patch = patch+1; other_patch < startIJK_vec.size(); ++other_patch) {
            const auto& [iFalseOther, iTrueOther, jFalseOther, jTrueOther, kFalseOther, kTrueOther] =
                getBoundaryPatchFaces(startIJK_vec[other_patch], endIJK_vec[other_patch]);
            bool isShared = false;
            if (startIJK_vec[other_patch][0] == endIJK_vec[patch][0]) {
                isShared = isShared || detectSharing(iTrue, iFalseOther);
            }
            if (endIJK_vec[other_patch][0] == startIJK_vec[patch][0]) {
                isShared = isShared || detectSharing(iFalse, iTrueOther);
            }
            if (startIJK_vec[other_patch][1] == endIJK_vec[patch][1]) {
                isShared = isShared || detectSharing(jTrue, jFalseOther);
            }
            if (endIJK_vec[other_patch][1] == startIJK_vec[patch][1]) {
                isShared = isShared || detectSharing(jFalse, jTrueOther);
            }
            if (startIJK_vec[other_patch][2] == endIJK_vec[patch][2]) {
                isShared = isShared || detectSharing(kTrue, kFalseOther);
            }
            if (endIJK_vec[other_patch][2] == startIJK_vec[patch][2]) {
                isShared = isShared || detectSharing(kFalse, kTrueOther);
            }
            if (isShared) {
                return isShared;
            }
        } // other patch for-loop
    } // patch for-loop
    return false;
}

int CpGridData::sharedFaceTag(const std::vector<std::array<int,3>>& startIJK_2Patches, const std::vector<std::array<int,3>>& endIJK_2Patches) const
{
    assert(startIJK_2Patches.size() == 2);
    assert(endIJK_2Patches.size() == 2);

    int faceTag = -1; // 0 represents I_FACE, 1 J_FACE, and 2 K_FACE. Use -1 for no sharing face case.
     
    if (patchesShareFace(startIJK_2Patches, endIJK_2Patches)) {
        
        const auto& detectSharing = [](const std::vector<int>& faceIdxs, const std::vector<int>& otherFaceIdxs){
            bool faceIsShared = false;
            for (const auto& face : faceIdxs) {
                for (const auto& otherFace : otherFaceIdxs) {
                    faceIsShared = faceIsShared || (face == otherFace);
                    if (faceIsShared) {
                        return faceIsShared; // should be true here
                    }
                }
            }
            return faceIsShared; // should be false here
        };
     
        const auto& [iFalse, iTrue, jFalse, jTrue, kFalse, kTrue] = this->getBoundaryPatchFaces(startIJK_2Patches[0], endIJK_2Patches[0]);
        const auto& [iFalseOther, iTrueOther, jFalseOther, jTrueOther, kFalseOther, kTrueOther] =
            this->getBoundaryPatchFaces(startIJK_2Patches[1], endIJK_2Patches[1]);


        bool isShared = false;

        // Check if patch1 lays on the left of patch2, so they might share an I_FACE that
        // for patch1 is false-oriented (contained in iFalse) and for patch2 is true-oriented (contained in iTrueOther).
        // patch2 | patch1
        if (startIJK_2Patches[0][0] == endIJK_2Patches[1][0]) {
            isShared = isShared || detectSharing(iFalse, iTrueOther);
            if (isShared) {
                faceTag = 0;
            }
        }
        // Check if patch1 lays on the right of patch2, so they might share an I_FACE that
        // for patch1 is true-oriented (contained in iTrue) and for patch2 is false-oriented (contained in iFalseOther).
        // patch1 | patch2
        if (endIJK_2Patches[0][0] == startIJK_2Patches[1][0]) {
            isShared = isShared || detectSharing(iTrue, iFalseOther);
            if (isShared) {
                faceTag = 0;
            }
        }
        // Check if patch1 lays in front of patch2, so they might share an J_FACE that
        // for patch1 is true-oriented (contained in jTrue) and for patch2 is false-oriented (contained in jFalseOther).
        //      patch2
        //   -----
        // patch1
        if (endIJK_2Patches[0][1] == startIJK_2Patches[1][1]) {
            isShared = isShared || detectSharing(jTrue, jFalseOther);
            if (isShared) {
                faceTag = 1;
            }
        }
        // Check if patch1 lays in back of patch2, so they might share an J_FACE that
        // for patch1 is false-oriented (contained in jFalse) and for patch2 is true-oriented (contained in jTrueOther).
        //      patch1
        //   -----
        // patch2
        if (startIJK_2Patches[0][1] == endIJK_2Patches[1][1]) {
            isShared = isShared || detectSharing(jFalse, jTrueOther);
            if (isShared) {
                faceTag = 1;
            }
        }
        // Check if patch1 lays on the bottom of patch2, so they might share an K_FACE that
        // for patch1 is true-oriented (contained in kTrue) and for patch2 is false-oriented (contained in kFalseOther).
        // patch2
        // -----
        // patch1
        if (endIJK_2Patches[0][2] == startIJK_2Patches[1][2]) {
            isShared = isShared || detectSharing(kTrue, kFalseOther);
            if (isShared) {
                faceTag = 2;
            }
        }
        // Check if patch1 lays on the top of patch2, so they might share an K_FACE that
        // for patch1 is false-oriented (contained in kFalse) and for patch2 is true-oriented (contained in kTrueOther).
        // patch1
        // -----
        // patch2
        if (startIJK_2Patches[0][2] == endIJK_2Patches[1][2]) {
            isShared = isShared || detectSharing(kFalse, kTrueOther);
            if (isShared) {
                faceTag = 2;
            }
        }
    }
    return faceTag; // -1 when no face is shared, otherwise: 0 (shared I_FACE), 1 (shared J_FACE), 2 (shared K_FACE)
}


std::vector<int>
CpGridData::getPatchesCells(const std::vector<std::array<int,3>>& startIJK_vec, const std::vector<std::array<int,3>>& endIJK_vec) const
{
    std::vector<int> all_cells;
    for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch){
        /// PATCH CELLS
        const auto& patch_cells = CpGridData::getPatchCells(startIJK_vec[patch], endIJK_vec[patch]);
        all_cells.insert(all_cells.end(), patch_cells.begin(), patch_cells.end());
    }
    return all_cells;
}

bool CpGridData::hasNNCs(const std::vector<int>& cellIndices) const
{
    bool hasNNC = false;
    for (const auto cellIdx : cellIndices)
    {
        bool cellHasNNC = false;
        Dune::cpgrid::Entity<0> entity(*this, cellIdx, true);
        auto cellFaces = this->cell_to_face_[entity];
        for (const auto face : cellFaces) {
            Dune::cpgrid::Entity<1> f(face.index(), true);
            enum face_tag tag = this->face_tag_[f];
            if (tag == face_tag::NNC_FACE) {
                cellHasNNC = true;
                break;
            }
        }
        hasNNC = hasNNC || cellHasNNC;
        if (hasNNC){
            break;
        }
    }
    return hasNNC;
}

void CpGridData::validStartEndIJKs(const std::vector<std::array<int,3>>& startIJK_vec,
                                   const std::vector<std::array<int,3>>& endIJK_vec) const
{
    if (startIJK_vec.size() == endIJK_vec.size()) {
        for (grid_size_t patch = 0; patch < startIJK_vec.size(); ++patch) {
            bool validPatch = true;
            for (int c = 0; c < 3; ++c) {
                // valid startIJK and endIJK for each patch
                validPatch = validPatch && (startIJK_vec[patch][c] < endIJK_vec[patch][c]);
                if (!validPatch) {
                    OPM_THROW(std::invalid_argument, "Invalid IJK-indices in LGR"+std::to_string(patch +1)+", end I/J/K need to be larger than start I/J/K.\n");
                }
            }
        }
    }
    else {
        OPM_THROW(std::invalid_argument, "Sizes of provided vectors with start and end IJK need to match.\n");
    }
}

bool CpGridData::compatibleSubdivisions(const std::vector<std::array<int,3>>& cells_per_dim_vec,
                                        const std::vector<std::array<int,3>>& startIJK_vec,
                                        const std::vector<std::array<int,3>>& endIJK_vec) const
{
    bool compatibleSubdivisions = true;
    if (startIJK_vec.size() > 1) {
        bool notAllowedYet = false;
        for (std::size_t level = 0; level < startIJK_vec.size(); ++level) {
            for (std::size_t otherLevel = level+1; otherLevel < startIJK_vec.size(); ++otherLevel) {
                const auto& sharedFaceTag = this-> sharedFaceTag({startIJK_vec[level], startIJK_vec[otherLevel]}, {endIJK_vec[level],endIJK_vec[otherLevel]});
                if(sharedFaceTag == -1){
                    break; // Go to the next "other patch"
                }
                if (sharedFaceTag == 0 ) {
                    notAllowedYet = notAllowedYet ||
                        ((cells_per_dim_vec[level][1] != cells_per_dim_vec[otherLevel][1]) || (cells_per_dim_vec[level][2] != cells_per_dim_vec[otherLevel][2]));
                }
                if (sharedFaceTag == 1) {
                    notAllowedYet = notAllowedYet ||
                        ((cells_per_dim_vec[level][0] != cells_per_dim_vec[otherLevel][0]) || (cells_per_dim_vec[level][2] != cells_per_dim_vec[otherLevel][2]));
                }
                if (sharedFaceTag == 2) {
                    notAllowedYet = notAllowedYet ||
                        ((cells_per_dim_vec[level][0] != cells_per_dim_vec[otherLevel][0]) || (cells_per_dim_vec[level][1] != cells_per_dim_vec[otherLevel][1]));
                }
                if (notAllowedYet){
                    compatibleSubdivisions = false;
                    break;
                }
            } // end-otherLevel-for-loop
        } // end-level-for-loop
    }// end-if-patchesShareFace
    return compatibleSubdivisions;
}

Geometry<3,3> CpGridData::cellifyPatch(const std::array<int,3>& startIJK, const std::array<int,3>& endIJK,
                                       const std::vector<int>& patch_cells,
                                       DefaultGeometryPolicy& cellifiedPatch_geometry,
                                       std::array<int,8>& cellifiedPatch_to_point,
                                       std::array<int,8>& allcorners_cellifiedPatch) const
{
    if (patch_cells.empty()){
        OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
    }
    if (patch_cells.size() == 1){
        return (*(this -> geometry_.geomVector(std::integral_constant<int,0>())))[EntityRep<0>(patch_cells[0], true)];
    }
    else{
        checkCuboidShape(patch_cells);
        // Get grid dimension.
        const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
        // Select 8 corners of the patch boundary to be the 8 corners of the 'cellified patch'.
        cellifiedPatch_to_point = { // Corner-index: (J*(grid_dim[0]+1)*(grid_dim[2]+1)) + (I*(grid_dim[2]+1)) +K
            // Index of corner '0' {startI, startJ, startK}
            (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + startIJK[2],
            // Index of corner '1' '{endI, startJ, startK}'
            (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + startIJK[2],
            // Index of corner '2' '{startI, endJ, startK}'
            (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + startIJK[2],
            // Index of corner '3' '{endI, endJ, startK}'
            (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + startIJK[2],
            // Index of corner '4' '{startI, startJ, endK}'
            (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1))+ endIJK[2],
            // Index of corner '5' '{endI, startJ, endK}'
            (startIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + endIJK[2],
            // Index of corner '6' '{startI, endJ, endK}'
            (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (startIJK[0]*(grid_dim[2]+1)) + endIJK[2],
            // Index of corner '7' {endI, endJ, endK}
            (endIJK[1]*(grid_dim[0]+1)*(grid_dim[2]+1)) + (endIJK[0]*(grid_dim[2]+1)) + endIJK[2]};
        EntityVariableBase<cpgrid::Geometry<0,3>>& cellifiedPatch_corners =
            *(cellifiedPatch_geometry.geomVector(std::integral_constant<int,3>()));
        cellifiedPatch_corners.resize(8);
        // Compute the center of the 'cellified patch' and its corners.
        Geometry<0,3>::GlobalCoordinate cellifiedPatch_center = {0., 0.,0.};
        for (int corn = 0; corn < 8; ++corn) {
            // FieldVector in DUNE 2.6 is missing operator/ using a loop
            for(int i=0; i < 3; ++i){
                cellifiedPatch_center[i] +=
                    (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellifiedPatch_to_point[corn]).center()[i]/8.;
            }
            cellifiedPatch_corners[corn] =
                (*(this -> geometry_.geomVector(std::integral_constant<int,3>()))).get(cellifiedPatch_to_point[corn]);
        }
        // Compute the volume of the 'cellified patch'.
        double cellifiedPatch_volume = 0.;
        for (const auto& idx : patch_cells) {
            cellifiedPatch_volume +=
                (*(this -> geometry_.geomVector(std::integral_constant<int,0>())))[EntityRep<0>(idx, true)].volume();
        }
        // Indices of 'all the corners', in this case, 0-7 (required to construct a Geometry<3,3> object).
        allcorners_cellifiedPatch = {0,1,2,3,4,5,6,7};
        // Create a pointer to the first element of "cellfiedPatch_to_point" (required to construct a Geometry<3,3> object).
        const int* cellifiedPatch_indices_storage_ptr = &allcorners_cellifiedPatch[0];
        // Construct (and return) the Geometry<3,3> of the 'cellified patch'.
        return Geometry<3,3>(cellifiedPatch_center, cellifiedPatch_volume,
                             cellifiedPatch_geometry.geomVector(std::integral_constant<int,3>()),
                             cellifiedPatch_indices_storage_ptr);
    }
}

std::tuple< const std::shared_ptr<CpGridData>,
            const std::vector<std::array<int,2>>,                // parent_to_refined_corners(~boundary_old_to_new_corners)
            const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces (~boundary_old_to_new_faces)
            const std::tuple<int, std::vector<int>>,             // parent_to_children_cells
            const std::vector<std::array<int,2>>,                // child_to_parent_faces
            const std::vector<std::array<int,2>>>                // child_to_parent_cells
CpGridData::refineSingleCell(const std::array<int,3>& cells_per_dim, const int& parent_idx) const
{
    // To store the LGR/refined-grid.
    std::vector<std::shared_ptr<CpGridData>> refined_data;
    std::shared_ptr<CpGridData> refined_grid_ptr = std::make_shared<CpGridData>(refined_data); // ccobj_
    auto& refined_grid = *refined_grid_ptr;
    DefaultGeometryPolicy& refined_geometries = refined_grid.geometry_;
    std::vector<std::array<int,8>>& refined_cell_to_point = refined_grid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face = refined_grid.cell_to_face_;
    Opm::SparseTable<int>& refined_face_to_point = refined_grid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell = refined_grid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>& refined_face_tags = refined_grid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;
    // Get parent cell
    const cpgrid::Geometry<3,3>& parent_cell = (*(geometry_.geomVector(std::integral_constant<int,0>())))[EntityRep<0>(parent_idx, true)];
    // Get parent cell corners.
    const std::array<int,8>& parent_to_point = this->cell_to_point_[parent_idx];
    const std::set<int> nonRepeated_parentCorners(parent_to_point.begin(), parent_to_point.end());
    if (nonRepeated_parentCorners.size() != 8){
        OPM_THROW(std::logic_error, "Cell is not a hexahedron. Cannot be refined (yet).");
    }
    // Refine parent cell
    parent_cell.refineCellifiedPatch(cells_per_dim, refined_geometries, refined_cell_to_point, refined_cell_to_face,
                                     refined_face_to_point, refined_face_to_cell, refined_face_tags, refined_face_normals,
                                     {1,1,1}, /*widthX, lengthY, heightZ*/ {1.}, {1.}, {1.});
    const std::vector<std::array<int,2>>& parent_to_refined_corners{
        // corIdx (J*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (I*(cells_per_dim[2]+1)) +K
        // replacing parent-cell corner '0' {0,0,0}
        {parent_to_point[0], 0},
        // replacing parent-cell corner '1' {cells_per_dim[0], 0, 0}
        {parent_to_point[1], cells_per_dim[0]*(cells_per_dim[2]+1)},
        // replacing parent-cell corner '2' {0, cells_perd_dim[1], 0}
        {parent_to_point[2], cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)},
        // replacing parent-cell corner '3' {cells_per_dim[0], cells_per_dim[1], 0}
        {parent_to_point[3], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))},
        // replacing parent-cell corner '4' {0, 0, cells_per_dim[2]}
        {parent_to_point[4], cells_per_dim[2]},
        // replacing parent-cell corner '5' {cells_per_dim[0], 0, cells_per_dim[2]}
        {parent_to_point[5], (cells_per_dim[0]*(cells_per_dim[2]+1)) + cells_per_dim[2]},
        // replacing parent-cell corner '6' {0, cells_per_dim[1], cells_per_dim[2]}
        {parent_to_point[6], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + cells_per_dim[2]},
        // replacing parent-cell corner '7' {cells_per_dim[0], cells_per_dim[1], cells_per_dim[2]}
        {parent_to_point[7], (cells_per_dim[1]*(cells_per_dim[0]+1)*(cells_per_dim[2]+1)) + (cells_per_dim[0]*(cells_per_dim[2]+1))
         + cells_per_dim[2]}};
    // Get parent_cell_to_face = { {face, orientation}, {another face, its orientation}, ...}.
    const auto& parent_cell_to_face = (this-> cell_to_face_[EntityRep<0>(parent_idx, true)]);
    // To store relation old-face to new-born-faces (children faces).
    std::vector<std::tuple<int,std::vector<int>>>  parent_to_children_faces;
    parent_to_children_faces.reserve(6);
    // To store child-to-parent-face relation. Child-faces ordered with the criteria introduced in refine()(Geometry.hpp)K,I,Jfaces.
    std::vector<std::array<int,2>> child_to_parent_faces;
    child_to_parent_faces.reserve(refined_face_to_cell.size());
    // Auxiliary integers to simplify new-born-face-index notation.
    const int& k_faces = cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1);
    const int& i_faces = (cells_per_dim[0]+1)*cells_per_dim[1]*cells_per_dim[2];
    // Populate parent_to_children_faces and child_to_parent_faces.
    for (const auto& face : parent_cell_to_face) {
        // Check face tag to identify the type of face (bottom, top, left, right, front, or back).
        auto& parent_face_tag = (this-> face_tag_[Dune::cpgrid::EntityRep<1>(face.index(), true)]);
        // To store the new born faces for each face.
        std::vector<int> children_faces; // Cannot reserve/resize "now", it depends of the type of face.
        // K_FACES
        if (parent_face_tag == face_tag::K_FACE) {
            children_faces.reserve(cells_per_dim[0]*cells_per_dim[1]);
            for (int j = 0; j < cells_per_dim[1]; ++j) {
                for (int i = 0; i < cells_per_dim[0]; ++i) {
                    int child_face;
                    if (!face.orientation()) // false -> BOTTOM FACE -> k=0
                        child_face = (j*cells_per_dim[0]) + i;
                    else // true -> TOP FACE -> k=cells_per_dim[2]
                        child_face = (cells_per_dim[2]*cells_per_dim[0]*cells_per_dim[1]) +(j*cells_per_dim[0]) + i;
                    children_faces.push_back(child_face);
                    child_to_parent_faces.push_back({child_face, face.index()});
                } // i-for-lopp
            } //j-for-loop
        } // if-K_FACE
        // I_FACES
        if (parent_face_tag == face_tag::I_FACE) {
            children_faces.reserve(cells_per_dim[1]*cells_per_dim[2]);
            for (int k = 0; k < cells_per_dim[2]; ++k) {
                for (int j = 0; j < cells_per_dim[1]; ++j) {
                    int child_face;
                    if (!face.orientation()) // false -> LEFT FACE -> i=0
                        child_face = k_faces + (k*cells_per_dim[1]) + j;
                    else // true -> RIGHT FACE -> i=cells_per_dim[0]
                        child_face = k_faces + (cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[1]) + j;
                    children_faces.push_back(child_face);
                    child_to_parent_faces.push_back({child_face, face.index()});
                } // j-for-loop
            } // k-for-loop
        } // if-I_FACE
        // J_FACES
        if (parent_face_tag == face_tag::J_FACE) {
            children_faces.reserve(cells_per_dim[0]*cells_per_dim[2]);
            for (int i = 0; i < cells_per_dim[0]; ++i) {
                for (int k = 0; k < cells_per_dim[2]; ++k) {
                    int child_face;
                    if (!face.orientation()) // false -> FRONT FACE -> j=0
                        child_face = k_faces + i_faces + (i*cells_per_dim[2]) + k;
                    else  // true -> BACK FACE -> j=cells_per_dim[1]
                        child_face = k_faces + i_faces  + (cells_per_dim[1]*cells_per_dim[0]*cells_per_dim[2])
                            + (i*cells_per_dim[2]) + k;
                    children_faces.push_back(child_face);
                    child_to_parent_faces.push_back({child_face, face.index()});
                } // k-for-loop
            } // i-for-loop
        } // if-J_FACE
        parent_to_children_faces.push_back(std::make_tuple(face.index(), children_faces));
    }
    std::tuple<int, std::vector<int>> parent_to_children_cells; // {parent cell index (in level0), {child0,...,childN (in level1)}}
    auto& [ parent_index, children_cells ] = parent_to_children_cells;
    children_cells.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
    // To store the child to parent cell relation.
    std::vector<std::array<int,2>> child_to_parent_cell; // {child index (in level1), parent cell index (in level0)}
    child_to_parent_cell.reserve(cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]);
    // Populate children_cells and child_to_parent_cell.
    for (int cell = 0; cell < cells_per_dim[0]*cells_per_dim[1]*cells_per_dim[2]; ++cell) {
        children_cells.push_back(cell);
        child_to_parent_cell.push_back({cell, parent_idx});
    }
    return {refined_grid_ptr, parent_to_refined_corners, parent_to_children_faces, parent_to_children_cells,
            child_to_parent_faces, child_to_parent_cell};
}

std::tuple< std::shared_ptr<CpGridData>,
            const std::vector<std::array<int,2>>,                // boundary_old_to_new_corners
            const std::vector<std::tuple<int,std::vector<int>>>, // boundary_old_to_new_faces
            const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_faces
            const std::vector<std::tuple<int,std::vector<int>>>, // parent_to_children_cell
            const std::vector<std::array<int,2>>,                // child_to_parent_faces
            const std::vector<std::array<int,2>>>                // child_to_parent_cells
CpGridData::refinePatch(const std::array<int,3>& cells_per_dim, const std::array<int,3>& startIJK,
                        const std::array<int,3>& endIJK) const
{
    // Coarse grid dimension (amount of cells in each direction).
    const std::array<int,3>& grid_dim = this -> logicalCartesianSize();
    // Check that the grid is a Cartesian one.
    grid_size_t gXYZ = grid_dim[0]*grid_dim[1]*grid_dim[2];
    if (global_cell_.size() != gXYZ){
        OPM_THROW(std::logic_error, "Grid is not Cartesian. Patch cannot be refined.");
    }
    // To store LGR/refined-grid.
    std::vector<std::shared_ptr<CpGridData>> refined_data;
    std::shared_ptr<CpGridData> refined_grid_ptr = std::make_shared<CpGridData>(refined_data); // ccobj_
    auto& refined_grid = *refined_grid_ptr;
    DefaultGeometryPolicy& refined_geometries = refined_grid.geometry_;
    std::vector<std::array<int,8>>& refined_cell_to_point = refined_grid.cell_to_point_;
    cpgrid::OrientedEntityTable<0,1>& refined_cell_to_face = refined_grid.cell_to_face_;
    Opm::SparseTable<int>& refined_face_to_point = refined_grid.face_to_point_;
    cpgrid::OrientedEntityTable<1,0>& refined_face_to_cell = refined_grid.face_to_cell_;
    cpgrid::EntityVariable<enum face_tag,1>& refined_face_tags = refined_grid.face_tag_;
    cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& refined_face_normals = refined_grid.face_normals_;
    // Patch dimension (amount of cells in each direction).
    const auto& patch_dim = getPatchDim(startIJK, endIJK);
    const auto& patch_faces = getPatchFaces(startIJK, endIJK);
    const auto& patch_cells = getPatchCells(startIJK, endIJK);
    checkCuboidShape(patch_cells);
    // Get dx,dy,dz for each cell of the patch to be refined
    const auto& [widthsX, lengthsY, heightsZ] = getWidthsLengthsHeights(startIJK, endIJK);
    // Construct the Geometry of the cellified patch.
    DefaultGeometryPolicy cellified_patch_geometry;
    std::array<int,8> cellifiedPatch_to_point;
    std::array<int,8> allcorners_cellifiedPatch;
    cpgrid::Geometry<3,3> cellified_patch = this -> cellifyPatch(startIJK, endIJK, patch_cells, cellified_patch_geometry,
                                                                 cellifiedPatch_to_point, allcorners_cellifiedPatch);

    // Refine the "cellified_patch".
    cellified_patch.refineCellifiedPatch(cells_per_dim,
                                         refined_geometries,
                                         refined_cell_to_point,
                                         refined_cell_to_face,
                                         refined_face_to_point,
                                         refined_face_to_cell,
                                         refined_face_tags,
                                         refined_face_normals,
                                         patch_dim,
                                         widthsX, lengthsY, heightsZ);

    // Some integers to reduce notation later.
    const int& xfactor = cells_per_dim[0]*patch_dim[0];
    const int& yfactor = cells_per_dim[1]*patch_dim[1];
    const int& zfactor = cells_per_dim[2]*patch_dim[2];
    // To store the relation between old-corner-indices and the equivalent new-born ones (laying on the patch boundary).
    std::vector<std::array<int,2>> boundary_old_to_new_corners;
    boundary_old_to_new_corners.reserve((2*(cells_per_dim[0]+1)*(cells_per_dim[2]+1))
                                        + (2*(cells_per_dim[1]-1)*(cells_per_dim[2]+1))
                                        + (2*(cells_per_dim[0]-1)*(cells_per_dim[1]-1)));
    for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                if ( (j == startIJK[1]) || (j == endIJK[1]) ){ // Corners in the front/back of the patch.
                    boundary_old_to_new_corners.push_back({
                            // Old corner index
                            (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                            // New-born corner index (equivalent corner).
                            (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                            + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                            + (cells_per_dim[2]*(k-startIJK[2])) });
                }
                if ( (i == startIJK[0]) || (i == endIJK[0]) ) { // Corners in the left/right of the patch.
                    boundary_old_to_new_corners.push_back({
                            // Old corner index.
                            (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                            // New-born corner index (equivalent corner).
                            (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                            + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                            + (cells_per_dim[2]*(k-startIJK[2]))});
                }
                if ( (k == startIJK[2]) || (k == endIJK[2]) ) { // Corners in the bottom/top of the patch.
                    boundary_old_to_new_corners.push_back({
                            // Old corner index.
                            (j*(grid_dim[2]+1)*(grid_dim[0]+1)) + (i*(grid_dim[2]+1)) + k,
                            // New-born corner index (equivalent corner)
                            (cells_per_dim[1]*(j-startIJK[1])*(xfactor +1)*(zfactor +1))
                            + (cells_per_dim[0]*(i-startIJK[0])*(zfactor +1))
                            + (cells_per_dim[2]*(k-startIJK[2]))});
                }
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // To store face-indices of faces on the boundary of the patch.
    std::vector<int> boundary_patch_faces;
    // Auxiliary integers to simplify notation.
    const int& bound_patch_faces = (2*patch_dim[1]*patch_dim[2]) + (patch_dim[0]*2*patch_dim[2]) + (patch_dim[0]*patch_dim[1]*2);
    boundary_patch_faces.reserve(bound_patch_faces);
    // To store relation between old-face-index and its new-born-face indices.
    std::vector<std::tuple<int, std::vector<int>>> boundary_old_to_new_faces; // {face index, its children-indices}
    boundary_old_to_new_faces.reserve(bound_patch_faces);
    // Auxiliary integers to simplify notation.
    const int& i_grid_faces =  (grid_dim[0]+1)*grid_dim[1]*grid_dim[2];
    const int& j_grid_faces =  grid_dim[0]*(grid_dim[1]+1)*grid_dim[2];
    // To store relation bewteen parent face and its children (all faces of the patch, not only the ones on the boundary).
    std::vector<std::tuple<int,std::vector<int>>> parent_to_children_faces;
    parent_to_children_faces.reserve(patch_faces.size());
    // To store relation child-face-index and its parent-face-index.
    std::vector<std::array<int,2>> child_to_parent_faces; // {child index (in 'level 1'), parent index (in 'level 0')}
    child_to_parent_faces.reserve(refined_face_to_cell.size());
    // Populate child_to_parent_faces, parent_to_children_faces, boundary_old_to_new_faces, boundary_faces.
    // I_FACES
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]+1; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                int face_idx = (j*(grid_dim[0]+1)*grid_dim[2]) + (i*grid_dim[2])+ k;
                int l =  (i-startIJK[0])*cells_per_dim[0]; // l playing the role of the corresponding "i index" in the LGR
                // To store new born faces, per face. CHILDREN-FACES ARE ORDERED AS IN refine(), Geometry.hpp
                std::vector<int> children_list;  // I_FACE ikj (xzy-direction)
                // l,m,n play the role of 'x,y,z-direction', lnm = fake ikj (how I_FACES are 'ordered' in refine())
                for (int n = (k-startIJK[2])*cells_per_dim[2];n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                    for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                        children_list.push_back((xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor) + (n*yfactor) + m);
                        child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) +(l*yfactor*zfactor)
                                + (n*yfactor) + m, face_idx});
                    } // end m-for-loop
                } // end n-for-loop
                // Add parent information of each face to "parent_to_children_faces".
                parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                if ((i == startIJK[0]) || (i == endIJK[0])) { // Detecting if the face is on the patch boundary.
                    boundary_patch_faces.push_back(face_idx);
                    // Associate each old face on the boundary of the patch with the new born ones.
                    boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                }
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // J_FACES
    for (int j = startIJK[1]; j < endIJK[1]+1; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]; ++k) {
                int face_idx = i_grid_faces + (j*grid_dim[0]*grid_dim[2]) + (i*grid_dim[2])+ k;
                int m =  (j-startIJK[1])*cells_per_dim[1]; // m playing the role of the corresponding "j index" in the LGR
                // To store new born faces, per face. CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                std::vector<int> children_list;  // J_FACE jik (yxz-direction)
                // l,m,n play the role of 'x,y,z-direction', mln = fake jik (how J_FACES are 'ordered' in refine())
                for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                    for (int n = (k-startIJK[2])*cells_per_dim[2]; n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                        children_list.push_back((xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                                + (m*xfactor*zfactor) + (l*zfactor)+n);
                        child_to_parent_faces.push_back({(xfactor*yfactor*(zfactor+1)) + ((xfactor+1)*yfactor*zfactor)
                                + (m*xfactor*zfactor) + (l*zfactor)+n, face_idx});
                    } // end n-for-loop
                } // end l-for-loop
                // Add parent information of each face to "parent_to_children_faces".
                parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                if ((j == startIJK[1]) || (j == endIJK[1])) { // Detecting if face is on the patch boundary.
                    boundary_patch_faces.push_back(face_idx);
                    // Associate each old face on the boundary of the patch with the new born ones.
                    boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                }
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // K_FACES
    for (int j = startIJK[1]; j < endIJK[1]; ++j) {
        for (int i = startIJK[0]; i < endIJK[0]; ++i) {
            for (int k = startIJK[2]; k < endIJK[2]+1; ++k) {
                int face_idx = i_grid_faces + j_grid_faces + (j*grid_dim[0]*(grid_dim[2]+1)) + (i*(grid_dim[2]+1))+ k;
                int n =  (k-startIJK[2])*cells_per_dim[2]; // n playing the role of the corresponding "k index" in the LGR
                // To store new born faces, per face. CHILDREN FACES ARE ORDERED AS IN refine(), Geometry.hpp
                std::vector<int> children_list;  // K_FACE kji (zyx-direction)
                // l,m,n play the role of 'x,y,z-direction', nml = fake kji (how K_FACES are 'ordered' in refine())
                for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                    for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                        children_list.push_back((n*xfactor*yfactor) + (m*xfactor)+ l);
                        child_to_parent_faces.push_back({(n*xfactor*yfactor) + (m*xfactor)+ l, face_idx});
                    } // end l-for-loop
                } // end m-for-loop
                // Add parent information of each face to "parent_to_children_faces".
                parent_to_children_faces.push_back(std::make_tuple(face_idx, children_list));
                if ((k == startIJK[2]) || (k == endIJK[2])) { // Detecting if the face is on the patch boundary.
                    boundary_patch_faces.push_back(face_idx);
                    // Associate each old face on the boundary of the patch with the new born ones.
                    boundary_old_to_new_faces.push_back(std::make_tuple(face_idx, children_list));
                }
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop
    // To store the relation between parent cell and its new-born-cells.
    // {parent index (coarse grid), {child 0 index, child 1 index, ... (refined grid)}}
    std::vector<std::tuple<int,std::vector<int>>> parent_to_children_cells;
    parent_to_children_cells.reserve(patch_dim[0]*patch_dim[1]*patch_dim[2]);
    // To store the relation between a new-born-cell and its parent cell.
    std::vector<std::array<int,2>> child_to_parent_cells; // {child index (refined grid), parent cell index (coarse grid)}
    child_to_parent_cells.reserve(xfactor*yfactor*zfactor);
    for (int k = 0; k < grid_dim[2]; ++k) {
        for (int j = 0; j < grid_dim[1]; ++j) {
            for (int i = 0; i < grid_dim[0]; ++i) {
                int cell_idx = (k*grid_dim[0]*grid_dim[1]) + (j*grid_dim[0]) +i;
                std::vector<int> children_list;
                if ( (i > startIJK[0]-1) && (i < endIJK[0]) && (j > startIJK[1]-1) && (j < endIJK[1])
                     && (k > startIJK[2]-1) && (k < endIJK[2])) {
                    for (int n = (k-startIJK[2])*cells_per_dim[2]; n < (k-startIJK[2]+1)*cells_per_dim[2]; ++n) {
                        for (int m = (j-startIJK[1])*cells_per_dim[1]; m < (j-startIJK[1]+1)*cells_per_dim[1]; ++m) {
                            for (int l = (i-startIJK[0])*cells_per_dim[0]; l < (i-startIJK[0]+1)*cells_per_dim[0]; ++l) {
                                children_list.push_back((n*xfactor*yfactor) + (m*xfactor) + l);
                                child_to_parent_cells.push_back({(n*xfactor*yfactor) + (m*xfactor) + l, cell_idx});
                            }// end l-for-loop
                        } // end m-for-loop
                    } // end n-for-loop
                    parent_to_children_cells.push_back(std::make_tuple(cell_idx, children_list));
                }// end if 'patch cells'
            } // end i-for-loop
        } // end j-for-loop
    } // end k-for-loop

    return {refined_grid_ptr, boundary_old_to_new_corners, boundary_old_to_new_faces, parent_to_children_faces,
            parent_to_children_cells, child_to_parent_faces, child_to_parent_cells};
}

bool CpGridData::mark(int refCount, const cpgrid::Entity<0>& element)
{
    if (refCount == -1) {
        OPM_THROW(std::logic_error, "Coarsening is not supported yet.");
    }
    // Check the cell to be marked for refinement has no NNC (no neighbouring connections). Throw otherwise.
    if (hasNNCs({element.index()}) && (refCount == 1)) {
        OPM_THROW(std::logic_error, "Refinement of cells with face representing an NNC is not supported, yet.");
    }
    assert((refCount == 0) || (refCount == 1)); // Do nothing (0), Refine (1), Coarsen (-1) not supported yet.
    if (mark_.empty()) {
        mark_.resize(this->size(0));
    }
    mark_[element.index()] = refCount;
    return (mark_[element.index()] == refCount);
}

int CpGridData::getMark(const cpgrid::Entity<0>& element) const
{
    return mark_.empty() ? 0 : mark_[element.index()];
}

bool CpGridData::preAdapt()
{
    // [Indirectly] Set mightVanish flags for elements that have been marked for refinement
    if(mark_.empty()) {
        return false;
    }
    else {
        for (int elemIdx = 0; elemIdx <  this-> size(0); ++elemIdx) {
            const auto& element = Dune::cpgrid::Entity<0>(*this, elemIdx, true);
            if (getMark(element) != 0)  // 1 (to be refined), 0 (do nothing), -1 (to be coarsened - not supported yet)
                return true;
        }
    }
    return false;
}

bool CpGridData::adapt()
{
    return preAdapt();
}

void CpGridData::postAdapt()
{
    mark_.resize(this->size(0), 0);
}

std::array<double,3> CpGridData::computeEclCentroid(const int idx) const
{
    // The following computation is the same as the one used in Eclipse Grid.
    const auto& cell_to_point_indices = this -> cell_to_point_[idx];
    std::array<double,8> X;
    std::array<double,8> Y;
    std::array<double,8> Z;
    for (int cornIdx = 0; cornIdx < 8; ++cornIdx) {
        X[cornIdx] = (this-> geometry_.geomVector(std::integral_constant<int,3>())
                      -> get(cell_to_point_indices[cornIdx])).center()[0];
        Y[cornIdx] = (this-> geometry_.geomVector(std::integral_constant<int,3>())
                      -> get(cell_to_point_indices[cornIdx])).center()[1];
        Z[cornIdx] = (this-> geometry_.geomVector(std::integral_constant<int,3>())
                      -> get(cell_to_point_indices[cornIdx])).center()[2];

    }
    return std::array<double,3> { { std::accumulate(X.begin(), X.end(), 0.0) / 8.0,
                                        std::accumulate(Y.begin(), Y.end(), 0.0) / 8.0,
                                        std::accumulate(Z.begin(), Z.end(), 0.0) / 8.0 } };
}

std::array<double,3> CpGridData::computeEclCentroid(const Entity<0>& elem) const
{
    return this->computeEclCentroid(elem.index());
}

} // end namespace cpgrid
} // end namespace Dune
