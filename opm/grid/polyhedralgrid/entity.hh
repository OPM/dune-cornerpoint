// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_ENTITY_HH
#define DUNE_POLYHEDRALGRID_ENTITY_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/entity.hh>

#include <opm/grid/grid_size.h>

namespace Dune
{

  // PolyhedralGridEntityBasic
  // -----------------

  /** \copydoc PolyhedralGridEntity
   *
   *  \nosubgrouping
   */
  template< int codim, int dim, class Grid >
  class PolyhedralGridEntityBasic
  {
  protected:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

  public:
    /** \name Attributes
     *  \{ */

    //! codimensioon of the entity
    static const int codimension = codim;
    //! dimension of the grid
    static const int dimension = Traits::dimension;
    //! dimension of the entity
    static const int mydimension = dimension - codimension;
    //! dimension of the world
    static const int dimensionworld = Traits::dimensionworld;

    /** \} */

    /** \name Types Required by DUNE
     *  \{ */

    //! coordinate type of the grid
    typedef typename Traits::ctype ctype;

    typedef typename Traits::Index Index;

    //! type of entity
    typedef typename Grid::template Codim< codimension >::Entity Entity;

    //! type of corresponding entity seed
    typedef typename Grid::template Codim< codimension >::EntitySeed EntitySeed;
    //! type of corresponding geometry
    typedef typename Traits::template Codim< codimension >::Geometry     Geometry;
    typedef typename Traits::template Codim< codimension >::GeometryImpl GeometryImpl;

    /** \} */

  protected:
    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Traits::ExtraData ExtraData;

  public:
    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    PolyhedralGridEntityBasic ()
    : data_( nullptr ),
      seed_()
    {}

    /** \brief construct a null entity with data pointer */
    explicit PolyhedralGridEntityBasic ( ExtraData data )
    : data_( data ),
      seed_()
    {}

    /** \brief construct an initialized entity
     */
    PolyhedralGridEntityBasic ( ExtraData data, const EntitySeed& seed )
    : data_( data )
    , seed_( seed )
    {}

    /** \} */

    /** \name Methods Shared by Entities of All Codimensions
     *  \{ */

    /** \brief obtain the name of the corresponding reference element
     *
     *  This type can be used to access the DUNE reference element.
     */
    GeometryType type () const
    {
      return data()->geometryType(seed_);
    }

    /** \brief obtain the level of this entity */
    int level () const
    {
      return 0;
    }

    /** \brief obtain the partition type of this entity */
    PartitionType partitionType () const
    {
      return InteriorEntity; // data()->partitionType( *this );
    }

    /** obtain the geometry of this entity */
    Geometry geometry () const
    {
      return Geometry( GeometryImpl( data(), seed_ ) );
    }

    /** \brief return EntitySeed */
    EntitySeed seed () const { return seed_; }

    /** \} */


    /** \name Methods Supporting the Grid Implementation
     *  \{ */
    ExtraData data() const { return data_; }

    Index index () const { return seed_.index(); }

    /** \} */

    bool equals(const PolyhedralGridEntityBasic& other) const
    {
      return seed_.equals(other.seed_);
    }

  protected:
    ExtraData   data_;
    EntitySeed  seed_;
  };



  // PolyhedralGridEntity
  // ------------

  template< int codim, int dim, class Grid >
  class PolyhedralGridEntity : public PolyhedralGridEntityBasic< codim, dim, Grid >
  {
    typedef PolyhedralGridEntityBasic< codim, dim, Grid > Base ;
  protected:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

  protected:
    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Traits::ExtraData ExtraData;

    using Base::seed_;
  public:
    typedef typename Base :: EntitySeed EntitySeed;
    using Base :: codimension ;
    using Base :: data ;

    PolyhedralGridEntity ()
    : Base()
    {}

    explicit PolyhedralGridEntity ( ExtraData data_param )
    : Base( data_param )
    {}

    PolyhedralGridEntity ( ExtraData data_param, const EntitySeed& seed )
    : Base( data_param, seed )
    {}

    grid_size_t subEntities( const grid_size_t cd ) const
    {
      if( cd == Base :: codimension )
        return 1;
      else
        return data()->subEntities( seed_, cd );
    }

    template< int cd >
    typename Grid::template Codim< cd >::EntityPointer
    subEntity ( int i ) const
    {
      typedef typename Traits::template Codim< cd >::EntityPointerImpl EntityPointerImpl;
      typedef typename Traits::template Codim< cd >::EntityImpl        EntityImpl;
      return EntityPointerImpl( EntityImpl( data(), data()->template subEntitySeed< cd >( seed_, i ) ) );
    }
  };


  // PolyhedralGridEntity for codimension 0
  // ----------------------------------

  /** \copydoc PolyhedralGridEntity
   *
   *  \nosubgrouping
   */
  template< int dim, class Grid >
  class PolyhedralGridEntity< 0, dim, Grid > : public PolyhedralGridEntityBasic< 0, dim, Grid >
  {
    typedef PolyhedralGridEntityBasic< 0, dim, Grid > Base ;
    typedef PolyhedralGridEntity< 0, dim, Grid > This;

  protected:
    typedef typename Base::Traits Traits;

    // type of extra data, e.g. a pointer to grid (here empty)
    typedef typename Base::ExtraData ExtraData;

    using Base::seed_;
  public:
    using Base::codimension ;
    using Base::data ;

  protected:
    typedef typename Traits :: LeafIntersectionIteratorImpl  LeafIntersectionIteratorImpl;
    typedef typename Traits::template Codim< codimension >::LocalGeometryImpl LocalGeometryImpl;

  public:
    /** \name Types Required by DUNE
     *  \{ */

    typedef typename Base :: EntitySeed EntitySeed;

    //! type of corresponding local geometry
    typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
    //! type of corresponding entity
    typedef typename Traits::template Codim< codimension >::Entity Entity;
    //! type of corresponding entity
    typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

    //! type of hierarchic iterator
    typedef typename Traits::HierarchicIterator        HierarchicIterator;
    //! type of leaf intersection iterator
    typedef typename Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    //! type of level intersection iterator
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    PolyhedralGridEntity ()
    : Base()
    {}

    /** \brief construct a null entity with data pointer */
    explicit PolyhedralGridEntity ( ExtraData data_param )
    : Base( data_param )
    {}

    /** \brief construct an initialized entity */
    PolyhedralGridEntity ( ExtraData data_param, const EntitySeed& seed )
    : Base( data_param, seed )
    {}

    /** \} */

    const This& dereference() const
    { return *this; }

    This& dereference()
    { return *this; }

    grid_size_t subEntities( const grid_size_t codim ) const
    {
      if( codim == 0 )
        return 1;
      else
        return data()->subEntities( seed_, codim );
    }

    template< int codim >
    int count () const
    {
      return subEntities( codim );
    }

    template< int codim >
    typename Grid::template Codim< codim >::EntityPointer
    subEntity ( int i ) const
    {
      typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
      typedef typename Traits::template Codim< codim >::EntityImpl        EntityImpl;
      return EntityPointerImpl( EntityImpl( data(), data()->template subEntitySeed< codim >( seed_, i ) ) );
    }

    bool hasBoundaryIntersections () const
    {
      return data()->hasBoundaryIntersections( this->seed() );
    }

    LeafIntersectionIterator ibegin () const
    {
      return LeafIntersectionIterator( LeafIntersectionIteratorImpl( data(), seed_, true ) );
    }

    LeafIntersectionIterator iend () const
    {
      return LeafIntersectionIterator( LeafIntersectionIteratorImpl( data(), seed_, false ) );
    }

    LeafIntersectionIterator  ileafbegin  () const { return ibegin(); }
    LevelIntersectionIterator ilevelbegin () const { return ibegin(); }

    LeafIntersectionIterator  ileafend  () const { return iend(); }
    LevelIntersectionIterator ilevelend () const { return iend(); }

    bool isLeaf () const
    {
      return true;
    }

    EntityPointer father () const
    {
      DUNE_THROW(InvalidStateException,"no father available");
      typedef typename Traits::template Codim< 0 >::EntityImpl EntityImpl;
      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
      return EntityPointer( EntityPointerImpl( EntityImpl( data() ) ) );
    }

    bool hasFather () const
    {
      return false;
    }

    LocalGeometry geometryInFather () const
    {
      DUNE_THROW(InvalidStateException,"no father available");
      return LocalGeometry( LocalGeometryImpl( data() ) );
    }

    HierarchicIterator hbegin ( int maxLevel ) const
    {
      return hend( maxLevel );
    }

    HierarchicIterator hend ( int ) const
    {
      typedef typename Traits :: HierarchicIteratorImpl HierarchicIteratorImpl ;
      return HierarchicIterator( HierarchicIteratorImpl( data(), false ) );
    }

    bool isRegular () const
    {
      return true;
    }

    bool isNew () const
    {
      return false;
    }

    bool mightVanish () const
    {
      return false;
    }

    /** \} */

  };


} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITY_HH
