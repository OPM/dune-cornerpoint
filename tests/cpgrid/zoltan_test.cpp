/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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
#include <config.h>
#define ONE_TO_ALL
#define BOOST_TEST_MODULE ZoltanTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/common/ZoltanGraphFunctions.hpp>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>


class MPIError {
public:
  /** @brief Constructor. */
  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
  /** @brief The error string. */
  std::string errorstring;
  /** @brief The mpi error code. */
  int errorcode;
};

#ifdef HAVE_MPI
void MPI_err_handler(MPI_Comm *, int *err_code, ...){
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}
#endif

struct MPIFixture
{
    MPIFixture()
    {
#if HAVE_MPI
    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    helper = &Dune::MPIHelper::instance(m_argc, m_argv);
#ifdef MPI_2
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#else
        MPI_Errhandler_create(MPI_err_handler, &handler);
        MPI_Errhandler_set(MPI_COMM_WORLD, handler);
#endif
#endif
    }
    ~MPIFixture()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif
    }
    Dune::MPIHelper* helper;
#if HAVE_MPI
    MPI_Errhandler handler;
#endif
};

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(zoltan)
{

    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    {
        [[maybe_unused]] const auto& inst = Dune::MPIHelper::instance(m_argc, m_argv);

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
        int rc;
        float ver;
        int procs=1;
        struct Zoltan_Struct *zz;
        int changes, numGidEntries, numLidEntries, numImport, numExport;
        int myRank;
        ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
        int *importProcs, *importToPart, *exportProcs, *exportToPart;
        rc = Zoltan_Initialize(m_argc, m_argv, &ver);
        BOOST_REQUIRE (rc == ZOLTAN_OK);

        //#ifdef ONE_TO_ALL
        //zz = Zoltan_Create(MPI_COMM_SELF);
        //#else
        zz = Zoltan_Create(MPI_COMM_WORLD);
        //#endif

        /* General parameters */
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
        Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
        Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
        Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "3");
#ifdef ONE_TO_ALL
        //    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string((long long)inst.size()).c_str());
#endif
        /* Graph parameters */

        Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
        Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

        MPI_Comm_size(MPI_COMM_WORLD, &procs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        Dune::CpGrid grid;
        Dune::cpgrid::CombinedGridWellGraph graph(grid, nullptr, {}, nullptr, false, Dune::EdgeWeightMethod::uniformEdgeWgt);
        std::array<int, 3> dims={{1, procs, procs}};
        std::array<double, 3> size={{ 1.0, 1.0, 1.0}};
#ifdef ONE_TO_ALL
        if (myRank==0)
#endif
            grid.createCartesian(dims, size);

        MPI_Barrier(MPI_COMM_WORLD);

        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, graph, true /*pretend empty grid, which we have since our grid is actually empty*/);

        BOOST_REQUIRE(grid.comm()==MPI_COMM_WORLD);
        if (myRank != 0)
        {
            BOOST_REQUIRE(grid.numCells()==0);
        }

        //ZOLTAN_TRACE_ENTER(zz, yo);
        rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                                 &changes,        /* 1 if partitioning was changed, 0 otherwise */
                                 &numGidEntries,  /* Number of integers used for a global ID */
                                 &numLidEntries,  /* Number of integers used for a local ID */
                                 &numImport,      /* Number of vertices to be sent to me */
                                 &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                                 &importLocalGids,   /* Local IDs of vertices to be sent to me */
                                 &importProcs,    /* Process rank for source of each incoming vertex */
                                 &importToPart,   /* New partition for each incoming vertex */
                                 &numExport,      /* Number of vertices I must send to other processes*/
                                 &exportGlobalGids,  /* Global IDs of the vertices I must send */
                                 &exportLocalGids,   /* Local IDs of the vertices I must send */
                                 &exportProcs,    /* Process to which I send each of the vertices */
                                 &exportToPart);  /* Partition to which each vertex will belong */
        //ZOLTAN_TRACE_DETAIL(zz, yo, "load balance");
        BOOST_REQUIRE (rc == ZOLTAN_OK);

        if(myRank > 0)
        {
            MPI_Status stat;
            int i=0;
            MPI_Recv(&i, 1, MPI_INT, myRank-1, 787, MPI_COMM_WORLD, &stat);
        }
        std::cout<<"Begin Rank "<<myRank<<":"<<std::endl;
        for ( int i=0; i < numExport; i++ )
        {
            std::cout<<"e"<<exportLocalGids[i]<<" ("<<exportGlobalGids[i]<<") => part="<<exportToPart[i]<<
                " proc="<<exportProcs[i]<<", ";
        }
        std::cout<<std::endl;
        for ( int i=0; i < numImport; i++ )
        {
            std::cout<<"i"<<importLocalGids[i]<<" ("<<importGlobalGids[i]<<") => part="<<importToPart[i]<<
                " proc="<<importProcs[i]<<", ";
        }
        std::cout<<"End Rank "<<myRank<<":"<<std::endl;
        std::cout<<std::endl;
        if(myRank <procs-1)
        {
            int i=0;
            MPI_Send(&i, 1, MPI_INT, myRank+1, 787, MPI_COMM_WORLD);
        }
#endif
    }
}
