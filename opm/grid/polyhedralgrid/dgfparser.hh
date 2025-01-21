// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
#define DUNE_POLYHEDRALGRID_DGFPARSER_HH

#include <algorithm>
#include <numeric>
#include <memory>
#include <utility>

#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/grid/io/file/dgfparser/blocks/polyhedron.hh>

#include <opm/grid/polyhedralgrid/gridfactory.hh>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#endif

namespace Dune
{
  // DGFGridFactory for PolyhedralGrid
  // ---------------------------------

  template< int dim, int dimworld, class coord_t >
  struct DGFGridFactory< PolyhedralGrid< dim, dimworld, coord_t > >
  {
    typedef PolyhedralGrid< dim, dimworld, coord_t > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input, MPICommunicator = MPIHelper::getCommunicator() )
      : gridPtr_(),
        grid_( nullptr )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream" );
      generate( input );
    }

    explicit DGFGridFactory ( const std::string &filename, MPICommunicator /* comm */ = MPIHelper::getCommunicator() )
      : gridPtr_(),
        grid_( nullptr )
    {
      std::ifstream input( filename );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found" );

#if HAVE_ECL_INPUT
      if( !DuneGridFormatParser::isDuneGridFormat( input ) )
      {
        Opm::Parser parser;
        const auto deck = parser.parseString( filename );
        std::vector<double> porv;

        gridPtr_.reset( new Grid( deck, porv ) );
        return ;
      }
      else
#endif
      {
        generate( input );
      }
    }

    Grid *grid () const
    {
      if( ! grid_ )
      {
        // set pointer to grid and avoid grid being deleted
        grid_ = gridPtr_.release();
      }
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection& /*intersection*/ ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection& ) const
    {
      return false;
    }

    bool haveBoundaryParameters () const { return false; }

    template< int codim >
    int numParameters () const
    {
      //return (codim == dimension ? numVtxParams_ : 0);;
      return 0;
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection& ) const
    {
      return DGFBoundaryParameter::defaultValue();;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity& )
    {
      static std::vector< double > dummy;
      return dummy;
    }

  private:
    int readVertices ( std::istream &input, std::vector< std::vector< double > > &vertices )
    {
      int dimWorld = Grid::dimensionworld ;
      dgf::VertexBlock vtxBlock( input, dimWorld );
      if( !vtxBlock.isactive() )
        DUNE_THROW( DGFException, "Vertex block not found" );

      vtxBlock.get( vertices, vtxParams_, numVtxParams_ );
      return vtxBlock.offset();
    }

    std::vector< std::vector< int > > readPolygons ( std::istream &input, int numVtx, int vtxOfs )
    {
      dgf::PolygonBlock polygonBlock( input, numVtx, vtxOfs );
      if( !polygonBlock.isactive() )
        DUNE_THROW( DGFException, "Polygon block not found" );

      std::vector< std::vector< int > > polygons;
      polygonBlock.get( polygons );
      return polygons;
    }

    std::vector< std::vector< int > > readPolyhedra ( std::istream &input, int numPolygons )
    {
      dgf::PolyhedronBlock polyhedronBlock( input, numPolygons );
      std::vector< std::vector< int > > polyhedra;
      if( polyhedronBlock.isactive() )
      {
        polyhedronBlock.get( polyhedra );
      }
      return polyhedra;
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, double *dest )
    {
      for( ; begin != end; ++begin )
        dest = std::copy( begin->begin(), begin->end(), dest );
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, int *dest, int *offset )
    {
      int size = 0;
      for( ; begin != end; ++begin )
      {
        *(offset++) = size;
        size += begin->size();
        dest = std::copy( begin->begin(), begin->end(), dest );
      }
      *offset = size;
    }

    void generate ( std::istream &input )
    {
      // check whether an interval block is active, otherwise read polyhedrons
      dgf::IntervalBlock intervalBlock( input );
      if( intervalBlock.isactive() )
      {
        if( intervalBlock.numIntervals() != 1 )
          DUNE_THROW( DGFException, "Currently, CpGrid can only handle 1 interval block." );

        if( intervalBlock.dimw() != dimworld )
          DUNE_THROW( DGFException, "CpGrid cannot handle an interval of dimension " << intervalBlock.dimw() << "." );
        const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

        std::vector< double > spacing( dimworld );
        for( int i=0; i<dimworld; ++i )
          spacing[ i ] = (interval.p[ 1 ][ i ] - interval.p[ 0 ][ i ]) / interval.n[ i ];

        gridPtr_.reset( new Grid( interval.n, spacing ) );
        return ;
      }
      else // polyhedral input
      {
        typedef std::vector< std::vector< double > > CoordinateVectorType;
        CoordinateVectorType nodes;

        typedef std::vector< std::vector< int > > IndexVectorType;
        IndexVectorType faces;
        IndexVectorType cells;

        const int vtxOfs = readVertices( input, nodes );

        faces = readPolygons ( input, nodes.size(), vtxOfs );
        cells = readPolyhedra( input, faces.size() );

        if( cells.empty() )
        {
          DUNE_THROW( DGFException, "Polyhedron block not found!" );
        }

        typedef GridFactory< Grid > GridFactoryType;
        typedef typename GridFactoryType :: Coordinate Coordinate ;

        GridFactoryType gridFactory;

        const int nNodes = nodes.size();
        Coordinate node( 0 );
        for( int i=0; i<nNodes; ++i )
        {
          for( int d=0; d<Coordinate::dimension; ++d )
            node[ d ] = nodes[ i ][ d ];

          gridFactory.insertVertex( node );
        }
        //nodes.swap( CoordinateVectorType() );

        // insert faces with type none/dim-1
        GeometryType type;
        type = Dune::GeometryTypes::none(Grid::dimension-1);
        std::vector< grid_size_t > numbers;

        const int nFaces = faces.size();
        for(int i = 0; i < nFaces; ++ i )
        {
          // copy values into appropriate data type
          std::vector<int>& face = faces[ i ];
          numbers.resize( face.size() );
          std::copy( face.begin(), face.end(), numbers.begin() );
          gridFactory.insertElement( type, numbers );
        }

        // free memory
        //faces.swap( IndexVectorType() );

        // insert cells with type none/dim
        type = Dune::GeometryTypes::none(Grid::dimension);

        const int nCells = cells.size();
        for(int i = 0; i < nCells; ++ i )
        {
          // copy values into appropriate data type
          std::vector<int>& cell = cells[ i ];
          numbers.resize( cell.size() );
          std::copy( cell.begin(), cell.end(), numbers.begin() );
          gridFactory.insertElement( type, numbers );
        }
        //cells.swap( IndexVectorType() );

        gridPtr_ = gridFactory.createGrid();
      } // end else branch

     // alternative conversion to polyhedral format that does not work yet.
#if 0
        else // convert dgf input to polyhedral
        {
          nodes.resize( dgf.nofvtx );
          // copy vertices
          std::copy( dgf.vtx.begin(), dgf.vtx.end(), nodes.begin() );

          for( const auto& node : nodes )
          {
            for( size_t i=0; i<node.size(); ++i )
              std::cout << node[i] << " ";
            std::cout << std::endl;
          }

          const grid_size_t nVx = dgf.elements[ 0 ].size();

          typedef std::vector< int > face_t;
          std::map< face_t, int > tmpFaces;

          const int nFaces = (nVx == dim+1) ? dim+1 : 2*dim;

          Dune::GeometryType type( (nVx == dim+1) ?
              Impl :: SimplexTopology< dim > :: type :: id :
              Impl :: CubeTopology   < dim > :: type :: id,
              dim );

          const auto& refElem = Dune::referenceElement< typename Grid::ctype, dim >( type );

          cells.resize( dgf.nofelements );

          face_t face;
          int faceNo = 0;
          for( int n = 0; n < dgf.nofelements; ++n )
          {
            const auto& elem = dgf.elements[ n ];
            auto& cell = cells[ n ];
            assert( elem.size() == nVx );
            cell.resize( nFaces );
            for(int f=0; f<nFaces; ++f )
            {
              const int nFaceVx = refElem.size(f, 1, dim);
              face.resize( nFaceVx );
              for( int j=0; j<nFaceVx; ++j )
              {
                face[ j ] = elem[ refElem.subEntity(f, 1, j , dim) ];
              }
              std::sort( face.begin(), face.end() );
              auto it = tmpFaces.find( face );
              int myFaceNo = -1;
              if( it == tmpFaces.end() )
              {
                myFaceNo = faceNo++;
                tmpFaces.insert( std::make_pair( face, myFaceNo ) );
              }
              else
                myFaceNo = it->second;

              cell[ f ] = myFaceNo;
            }

            /*
            for( const auto& c : cell )
            {
              std::cout << c << " ";
            }
            std::cout << std::endl;
            */
          }
#endif
    }

    mutable std::unique_ptr< Grid > gridPtr_;
    mutable Grid* grid_;
    int numVtxParams_;
    std::vector< std::vector< double > > vtxParams_;
  };



  // DGFGridInfo for PolyhedralGrid
  // ------------------------------

  template< int dim, int dimworld >
  struct DGFGridInfo< PolyhedralGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
      return 0;
    }

    static double refineWeight ()
    {
      return 0;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
