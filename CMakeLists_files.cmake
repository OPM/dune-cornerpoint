# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up six lists:
# MAIN_SOURCE_FILES     List of compilation units which will be included in
#                       the library. If it isn't on this list, it won't be
#                       part of the library. Please try to keep it sorted to
#                       maintain sanity.
#
# TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
# TEST_DATA_FILES       Files from the source three that should be made
#                       available in the corresponding location in the build
#                       tree in order to run tests there.
#
# EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#                       build, but which is not part of the library nor is
#                       run as tests.
#
# PUBLIC_HEADER_FILES   List of public header files that should be
#                       distributed together with the library. The source
#                       files can of course include other files than these;
#                       you should only add to this list if the *user* of
#                       the library needs it.
#
# ATTIC_FILES           Unmaintained files. This for the projects developers
#                       only. Don't expect these files to build.

# originally generated with the command:
# find dune -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
  opm/grid/cpgrid/Intersection.cpp
  opm/grid/cpgrid/CpGridData.cpp
  opm/grid/cpgrid/CpGrid.cpp
  opm/grid/cpgrid/DataHandleWrappers.cpp
  opm/grid/cpgrid/GridHelpers.cpp
  opm/grid/cpgrid/Iterators.cpp
  opm/grid/cpgrid/Indexsets.cpp
  opm/grid/cpgrid/PartitionTypeIndicator.cpp
  opm/grid/cpgrid/CpGridUtilities.cpp
  opm/grid/cpgrid/processEclipseFormat.cpp
  opm/grid/common/GeometryHelpers.cpp
  opm/grid/common/GridPartitioning.cpp
  opm/grid/common/MetisPartition.cpp
  opm/grid/common/WellConnections.cpp
  opm/grid/common/ZoltanGraphFunctions.cpp
  opm/grid/common/ZoltanPartition.cpp
  opm/grid/CellQuadrature.cpp
  opm/grid/ColumnExtract.cpp
  opm/grid/FaceQuadrature.cpp
  opm/grid/GraphOfGrid.cpp
  opm/grid/GraphOfGridWrappers.cpp
  opm/grid/GridHelpers.cpp
  opm/grid/GridManager.cpp
  opm/grid/GridUtilities.cpp
  opm/grid/MinpvProcessor.cpp
  opm/grid/cart_grid.c
  opm/grid/cornerpoint_grid.c
  opm/grid/cpgpreprocess/facetopology.c
  opm/grid/cpgpreprocess/geometry.c
  opm/grid/cpgpreprocess/preprocess.c
  opm/grid/cpgpreprocess/uniquepoints.c
  opm/grid/UnstructuredGrid.c
  opm/grid/grid_equal.cpp
  opm/grid/utility/compressedToCartesian.cpp
  opm/grid/utility/cartesianToCompressed.cpp
  opm/grid/utility/StopWatch.cpp
  opm/grid/utility/WachspressCoord.cpp
  )

if (opm-common_FOUND)
  list(APPEND MAIN_SOURCE_FILES
		opm/grid/utility/VelocityInterpolation.cpp
		opm/grid/transmissibility/trans_tpfa.c)
endif()

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
  tests/p2pcommunicator_test.cc
  tests/test_cartgrid.cpp
  tests/test_cpgrid.cpp
  tests/test_communication_utils.cpp
  tests/test_column_extract.cpp
  tests/cpgrid/addLgrsOnDistributedGrid_test.cpp
  tests/cpgrid/distribution_test.cpp
  tests/cpgrid/entityrep_test.cpp
  tests/cpgrid/entity_test.cpp
  tests/cpgrid/facetag_test.cpp
  tests/cpgrid/orientedentitytable_test.cpp
  tests/cpgrid/partition_iterator_test.cpp
  tests/cpgrid/zoltan_test.cpp
  tests/test_cellCentroid_polyhedralGrid.cpp
  tests/test_compressed_cartesian_mapping.cpp
  tests/test_elementchunks.cpp
  tests/test_geom2d.cpp
  tests/test_gridutilities.cpp
  tests/test_lookupdata_polyhedral.cpp
  tests/test_minpvprocessor.cpp
  tests/test_polyhedralgrid.cpp
  tests/test_quadratures.cpp
  tests/test_repairzcorn.cpp
  tests/test_sparsetable.cpp
  tests/test_subgridpart.cpp
	)

if(Boost_VERSION_STRING VERSION_GREATER 1.53)
	list(APPEND TEST_SOURCE_FILES
	  tests/cpgrid/adapt_cpgrid_test.cpp
	  tests/cpgrid/avoidNNCinLGRs_test.cpp
	  tests/cpgrid/avoidNNCinLGRsCpGrid_test.cpp
	  tests/cpgrid/cuboidShape_test.cpp
	  tests/cpgrid/disjointPatches_test.cpp
	  tests/cpgrid/eclCentroid_test.cpp
	  tests/cpgrid/geometry_test.cpp
	  tests/cpgrid/grid_lgr_test.cpp
	  tests/cpgrid/inactiveCell_lgr_test.cpp
	  tests/cpgrid/lgr_cartesian_idx_test.cpp
	  tests/cpgrid/lgrIJK_test.cpp
		tests/cpgrid/lgr_coord_zcorn_test.cpp
	  tests/cpgrid/lookUpCellCentroid_cpgrid_test.cpp
	  tests/cpgrid/lookupdataCpGrid_test.cpp
	  tests/cpgrid/replace_lgr1_corner_idx_by_lgr2_corner_idx_test.cpp
	  tests/cpgrid/replace_lgr1_face_idx_by_lgr2_face_idx_test.cpp
	  tests/cpgrid/shifted_cart_test.cpp
	  tests/test_graphofgrid.cpp
	  tests/test_graphofgrid_parallel.cpp
  )
endif()

if(HAVE_ECL_INPUT)
  list(APPEND TEST_SOURCE_FILES
		tests/test_regionmapping.cpp
		tests/test_ug.cpp
		tests/cpgrid/grid_nnc.cpp
		tests/cpgrid/grid_pinch.cpp
	)
endif()

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
     tests/CORNERPOINT_ACTNUM.DATA
     tests/compressed_gridproperty.data
     tests/FIVE.DATA
     tests/FIVE_ACTNUM.DATA
     tests/FIVE_PINCH.DATA
     tests/FIVE_PINCH_NOGAP.DATA
     tests/FIVE_PINCH_NOGAP2.DATA
     tests/FIVE_PINCH_NOGAP3.DATA
  )

# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
  examples/finitevolume/finitevolume.cc
  examples/mirror_grid.cpp
  examples/griditer.cpp
  )

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
  examples/mirror_grid.cpp
  )
if(HAVE_ECL_INPUT)
  list(APPEND EXAMPLE_SOURCE_FILES examples/grdecl2vtu.cpp)
  list(APPEND PROGRAM_SOURCE_FILES examples/grdecl2vtu.cpp)
endif()

# originally generated with the command:
# find dune -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
  opm/grid/common/CommunicationUtils.hpp
  opm/grid/common/GeometryHelpers.hpp
  opm/grid/common/GridAdapter.hpp
  opm/grid/common/GridPartitioning.hpp
  opm/grid/common/Volumes.hpp
  opm/grid/common/p2pcommunicator.hh
  opm/grid/common/p2pcommunicator_impl.hh
  opm/grid/common/WellConnections.hpp
  opm/grid/cpgrid/CartesianIndexMapper.hpp
  opm/grid/cpgrid/CpGridData.hpp
  opm/grid/cpgrid/CpGridDataTraits.hpp
  opm/grid/cpgrid/CpGridUtilities.hpp
  opm/grid/cpgrid/DataHandleWrappers.hpp
  opm/grid/cpgrid/DefaultGeometryPolicy.hpp
  opm/grid/cpgrid/dgfparser.hh
  opm/grid/cpgrid/Entity2IndexDataHandle.hpp
  opm/grid/cpgrid/Entity.hpp
  opm/grid/cpgrid/EntityRep.hpp
  opm/grid/cpgrid/Geometry.hpp
  opm/grid/cpgrid/GlobalIdMapping.hpp
  opm/grid/cpgrid/GridHelpers.hpp
  opm/grid/cpgrid/LevelCartesianIndexMapper.hpp
  opm/grid/CpGrid.hpp
  opm/grid/cpgrid/Indexsets.hpp
  opm/grid/cpgrid/Intersection.hpp
  opm/grid/cpgrid/Iterators.hpp
  opm/grid/LookUpCellCentroid.hh
  opm/grid/LookUpData.hh
  opm/grid/cpgrid/OrientedEntityTable.hpp
  opm/grid/cpgrid/ParentToChildrenCellGlobalIdHandle.hpp
	opm/grid/cpgrid/ParentToChildCellToPointGlobalIdHandle.hpp
  opm/grid/cpgrid/PartitionIteratorRule.hpp
  opm/grid/cpgrid/PartitionTypeIndicator.hpp
  opm/grid/cpgrid/PersistentContainer.hpp
  opm/grid/common/CartesianIndexMapper.hpp
  opm/grid/common/GridEnums.hpp
  opm/grid/common/LevelCartesianIndexMapper.hpp
  opm/grid/common/MetisPartition.hpp
  opm/grid/common/SubGridPart.hpp
  opm/grid/common/ZoltanGraphFunctions.hpp
  opm/grid/common/ZoltanPartition.hpp
  opm/grid/polyhedralgrid/capabilities.hh
  opm/grid/polyhedralgrid/cartesianindexmapper.hh
  opm/grid/polyhedralgrid/declaration.hh
  opm/grid/polyhedralgrid/dgfparser.hh
  opm/grid/polyhedralgrid/entity.hh
  opm/grid/polyhedralgrid/entitypointer.hh
  opm/grid/polyhedralgrid/entityseed.hh
  opm/grid/polyhedralgrid/geometry.hh
  opm/grid/polyhedralgrid/gridhelpers.hh
  opm/grid/polyhedralgrid/grid.hh
  opm/grid/polyhedralgrid/gridview.hh
  opm/grid/polyhedralgrid/levelcartesianindexmapper.hh
  opm/grid/polyhedralgrid.hh
  opm/grid/polyhedralgrid/idset.hh
  opm/grid/polyhedralgrid/indexset.hh
  opm/grid/polyhedralgrid/intersection.hh
  opm/grid/polyhedralgrid/intersectioniterator.hh
  opm/grid/polyhedralgrid/iterator.hh
  opm/grid/polyhedralgrid/persistentcontainer.hh
  opm/grid/UnstructuredGrid.h
  opm/grid/CellQuadrature.hpp
  opm/grid/ColumnExtract.hpp
  opm/grid/FaceQuadrature.hpp
  opm/grid/GraphOfGrid.hpp
  opm/grid/GraphOfGridWrappers.hpp
  opm/grid/GridHelpers.hpp
  opm/grid/GridManager.hpp
  opm/grid/GridUtilities.hpp
  opm/grid/MinpvProcessor.hpp
  opm/grid/RepairZCORN.hpp
  opm/grid/cart_grid.h
  opm/grid/cornerpoint_grid.h
  opm/grid/cpgpreprocess/facetopology.h
  opm/grid/cpgpreprocess/geometry.h
  opm/grid/cpgpreprocess/preprocess.h
  opm/grid/cpgpreprocess/uniquepoints.h
  opm/grid/transmissibility/trans_tpfa.h
  opm/grid/transmissibility/TransTpfa.hpp
  opm/grid/transmissibility/TransTpfa_impl.hpp
  opm/grid/utility/compressedToCartesian.hpp
  opm/grid/utility/cartesianToCompressed.hpp
  opm/grid/utility/createThreadIterators.hpp
  opm/grid/utility/ElementChunks.hpp
  opm/grid/utility/IteratorRange.hpp
  opm/grid/utility/OpmWellType.hpp
  opm/grid/utility/RegionMapping.hpp
  opm/grid/utility/SparseTable.hpp
  opm/grid/utility/StopWatch.hpp
  opm/grid/utility/VariableSizeCommunicator.hpp
  opm/grid/utility/VelocityInterpolation.hpp
  opm/grid/utility/WachspressCoord.hpp
  opm/grid/utility/platform_dependent/disable_warnings.h
  opm/grid/utility/platform_dependent/reenable_warnings.h
  )
