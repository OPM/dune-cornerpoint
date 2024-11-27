#include "config.h"               // know what grids are present
#include <iostream>               // for input/output to shell
#include <fstream>                // for input/output to files
#include <vector>                 // STL vector class
#include <array>


// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/common/mcmgmapper.hh> // mapper class

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh> // include mpi helper class

// checks for defined gridtype and inlcudes appropriate dgfparser implementation
//#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include "vtkout.hh"
// #include"transportproblem2.hh"
#include "initialize.hh"
#include "evolve.hh"

#include "opm/grid/CpGrid.hpp"

typedef Dune::CpGrid GridType;

//===============================================================
// the time loop function working for all types of grids
//===============================================================

template<class G>
void timeloop(const G& grid, double tend)
{
    // make a mapper for codim 0 entities in the leaf grid
#if DUNE_VERSION_GTE(DUNE_GEOMETRY, 2, 8)
    Dune::MultipleCodimMultipleGeomTypeMapper<typename G::LeafGridView>  mapper(grid.leafGridView(), Dune::mcmgElementLayout());
#else
    Dune::LeafMultipleCodimMultipleGeomTypeMapper<G> mapper(grid, Dune::mcmgElementLayout());
#endif

    // allocate a vector for the concentration
    std::vector<double> c(mapper.size());

    // initialize concentration with initial values
    initialize(grid,mapper,c);                           /*@\label{fvc:init}@*/

    vtkout(grid,c,"concentration",0,0.0);

    // now do the time steps
    double t=0,dt;
    int k=0;
    const double saveInterval = 0.1;
    double saveStep = 0.1;
    int counter = 1;

    while (t<tend)                                       /*@\label{fvc:loop0}@*/
    {
        // augment time step counter
        ++k;

        // apply finite volume scheme
        evolve(grid,mapper,c,t,dt);

        // augment time
        t += dt;

        // check if data should be written
        if (t >= saveStep)
        {
            // write data
            vtkout(grid,c,"concentration",counter,t);

            // increase counter and saveStep for next interval
            saveStep += saveInterval;
            ++counter;
        }

        // print info about time, timestep size and counter
        std::cout << "s=" << grid.size(0)
        << " k=" << k << " t=" << t << " dt=" << dt << std::endl;
    }                                              /*@\label{fvc:loop1}@*/

    // output results
    vtkout(grid,c,"concentration",counter,tend);     /*@\label{fvc:file}@*/
}



void initGrid(const Dune::ParameterTree &param, GridType& grid)
{
    std::string fileformat = param.get<std::string>("fileformat");
    if (fileformat == "sintef_legacy") {
        OPM_THROW(std::runtime_error, "sintef_legacy fileformat no longer supported.");
    }
#if HAVE_OPM_PARSER
    else if (fileformat == "eclipse") {
        std::string filename = param.get<std::string>("filename");
        if (param.hasKey("z_tolerance")) {
            std::cerr << "****** Warning: z_tolerance parameter is obsolete, use PINCH in deck input instead\n";
        }
        bool periodic_extension = param.get<bool>("periodic_extension", false);
        bool turn_normals = param.get<bool>("turn_normals", false);

        Opm::ParseContext parseContext;
        Opm::Parser parser;
        auto deck = parser.parseFile(filename , parseContext);
        const int* actnum = deck.hasKeyword("ACTNUM") ? deck.getKeyword("ACTNUM").getIntData().data() : nullptr;
        Opm::EclipseGrid ecl_grid(deck , actnum);

        grid.processEclipseFormat(ecl_grid, periodic_extension, turn_normals);
    }
#endif
    else if (fileformat == "cartesian") {
        std::array<int, 3> dims = {{ param.get<int>("nx", 1),
                                param.get<int>("ny", 1),
                                param.get<int>("nz", 1) }};
        std::array<double, 3> cellsz = {{ param.get<double>("dx", 1.0),
                                     param.get<double>("dy", 1.0),
                                     param.get<double>("dz", 1.0) }};
        grid.createCartesian(dims, cellsz);
    } else {
        OPM_THROW(std::runtime_error, "Unknown file format string: " + fileformat);
    }
}


//===============================================================
// The main function creates objects and does the time loop
//===============================================================

int main(int argc , char ** argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc,argv);

    Dune::ParameterTree param;
    Dune::ParameterTreeParser::readINITree( argv[ 1 ], param );

    // start try/catch block to get error messages from dune
    try {
        using namespace Dune;

        GridType grid;
        initGrid( param , grid );
        // do time loop until end time 0.5
        timeloop(grid, 0.5);
    } catch (Dune::Exception & e) {
        std::cout << "DUNE ERROR: " << e.what() << std::endl;
        return 1;
    } catch (std::exception & e) {
        std::cout << "STL ERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cout << "Unknown ERROR" << std::endl;
        return 1;
    }
    // done
    return 0;
}
