// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_CAPABILITIES_HH
#define DUNE_POLYHEDRALGRID_CAPABILITIES_HH

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>

//- dune-metagrid includes
#include <opm/grid/polyhedralgrid/declaration.hh>

//#include <opm/grid/grid_size.h>

namespace Dune
{

  namespace Capabilities
  {

    // Capabilities from dune-grid
    // ---------------------------

    template< int dim, int dimworld, class coord_t >
    struct hasSingleGeometryType< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
      static const size_t topologyId = ~0u;
    };


    template< int dim, int dimworld, class coord_t >
    struct isCartesian< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };


    template< int dim, int dimworld, class coord_t, int codim >
    struct hasEntity< PolyhedralGrid< dim, dimworld, coord_t >, codim >
    {
      static const bool v = (codim == 0 || codim == 1 || codim == dim);
    };


    template< int dim, int dimworld, class coord_t, int codim >
    struct hasEntityIterator< PolyhedralGrid< dim, dimworld, coord_t >, codim >
    {
      static const bool v = (codim == 0 || codim == 1 || codim == dim);
    };

    template< int dim, int dimworld, class coord_t, int codim >
    struct canCommunicate< PolyhedralGrid< dim, dimworld, coord_t >, codim >
    {
        static const bool v = false;
    };


    template< int dim, int dimworld, class coord_t >
    struct hasBackupRestoreFacilities< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class coord_t >
    struct isLevelwiseConforming< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class coord_t >
    struct isLeafwiseConforming< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class coord_t >
    struct threadSafe< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };

    template< int dim, int dimworld, class coord_t >
    struct viewThreadSafe< PolyhedralGrid< dim, dimworld, coord_t > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_CAPABILITIES_HH
