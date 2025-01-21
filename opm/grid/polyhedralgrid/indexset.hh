// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_INDEXSET_HH
#define DUNE_POLYHEDRALGRID_INDEXSET_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

#include <opm/grid/grid_size.h>
#include <opm/grid/polyhedralgrid/declaration.hh>

namespace Dune
{

  // PolyhedralGridIndexSet
  // --------------

  template< int dim, int dimworld, typename coord_t >
  class PolyhedralGridIndexSet
      : public IndexSet< PolyhedralGrid< dim, dimworld, coord_t >, PolyhedralGridIndexSet< dim, dimworld, coord_t >, int >
  {
    typedef PolyhedralGrid<dim, dimworld, coord_t > GridType;

  protected:
    typedef PolyhedralGridIndexSet< dim, dimworld, coord_t > This;
      typedef IndexSet< GridType, This, int > Base;

    typedef typename std::remove_const< GridType >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;

    typedef typename Base::IndexType IndexType;

    PolyhedralGridIndexSet ( const GridType& grid )
        : grid_(&grid)
    {
    }

    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity::codimension >( entity );
    }

    template< int cd >
    IndexType index ( const typename Traits::template Codim< cd >::Entity &entity ) const
    {
      return entity.impl().index();
    }

    template< int cd >
    IndexType subIndex ( const typename Traits::template Codim< cd >::Entity &entity, int i, grid_size_t codim ) const
    {
      return subIndex( entity, i, codim );
    }

    template< class Entity >
    IndexType subIndex ( const Entity &entity, int i, grid_size_t codim ) const
    {
      if( codim == 0 )
        return index( entity );
      else if ( codim == 1 )
        return index( entity.impl().template subEntity< 1 > ( i ) );
      else if ( codim == dimension )
      {
        return index( entity.impl().template subEntity< dimension > ( i ) );
      }
      else
      {
        DUNE_THROW(NotImplemented,"codimension not available");
        return IndexType( -1 );
      }
    }

    IndexType size ( GeometryType type ) const
    {
      return grid().size( type );
    }

    int size ( int codim ) const
    {
      return grid().size( codim );
    }

    template< class Entity >
    bool contains ( const Entity &entity ) const
    {
        return index(entity) >= 0 && index(entity) < size(Entity::codimension);
    }

    const std::vector< GeometryType > &geomTypes ( int codim ) const
    {
        return grid().geomTypes(codim);
    }

    const std::vector< GeometryType >& types(int codim) const
    {
        return grid().geomTypes(codim);
    }

    const GridType& grid() const { assert( grid_ ); return *grid_; }

  protected:
    const GridType *grid_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_INDEXSET_HH
