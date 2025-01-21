//===========================================================================
//
// File: readEclipseFormat.cpp
//
// Created: Fri Jun 12 09:16:59 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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


#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "CpGridData.hpp"
#include "Geometry.hpp"

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/grid/common/GeometryHelpers.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/Indexsets.hpp>

#include <opm/grid/cpgpreprocess/preprocess.h>
#include <opm/grid/MinpvProcessor.hpp>
#include <opm/grid/RepairZCORN.hpp>
#include <opm/grid/utility/StopWatch.hpp>

#include <cstddef>
#include <fstream>
#include <iostream>
#include <initializer_list>
#include <set>
#include <string>
#include <utility>

namespace Dune
{

    using NNCMap = std::set<std::pair<int, int>>;
    using NNCMaps = std::array<NNCMap, 2>;
    enum NNCMapsIndex { PinchNNC = 0,
                        ExplicitNNC = 1 };


    // Forward declarations.
    namespace
    {
#if HAVE_ECL_INPUT
        std::vector<double>
        getSanitizedZCORN(const ::Opm::EclipseGrid& ecl_grid,
                          const ::std::vector<int>& actnum);
#endif

        typedef std::array<int, 3> coord_t;
        typedef std::array<double, 8> cellz_t;

        cellz_t getCellZvals(const coord_t& c, const coord_t& n, const double* z);

#if HAVE_ECL_INPUT
        void addOuterCellLayer(const grdecl& original,
                               std::vector<double>& new_coord,
                               std::vector<double>& new_zcorn,
                               std::vector<int>& new_actnum,
                               grdecl& output);
#endif

        void removeOuterCellLayer(processed_grid& grid);
        // void removeUnusedNodes(processed_grid& grid); // NOTE: not deleted, see comment at definition.
        void buildTopo(const processed_grid& output,
                       const NNCMaps& nnc,
                       std::vector<int>& global_cell,
                       cpgrid::OrientedEntityTable<0, 1>& c2f,
                       cpgrid::OrientedEntityTable<1, 0>& f2c,
                       Opm::SparseTable<int>& f2p,
                       std::vector<std::array<int,8> >& c2p,
                       std::vector<int>& face_to_output_face);
        void buildGeom(const processed_grid& output,
                       const cpgrid::OrientedEntityTable<0, 1>& c2f,
                       const std::vector<std::array<int,8> >& c2p,
                       const std::vector<int>& face_to_output_face,
                       const std::unordered_map<size_t, double>& aquifer_cell_volumes,
                       cpgrid::EntityVariable<cpgrid::Geometry<3, 3>, 0>& cell_geom,
                       cpgrid::EntityVariable<cpgrid::Geometry<2, 3>, 1>& face_geom,
                       std::shared_ptr<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> point_geom,
                       cpgrid::SignedEntityVariable<FieldVector<double, 3> , 1>& normals,
                       bool turn_normals);
    } // anon namespace




namespace cpgrid
{

#if HAVE_ECL_INPUT
    std::vector<std::size_t>
    CpGridData::processEclipseFormat(const Opm::EclipseGrid* ecl_grid_ptr,
                                     Opm::EclipseState* ecl_state,
                                     bool periodic_extension, bool turn_normals, bool clip_z,
                                     bool pinchActive)
    {
        std::vector<std::size_t> removed_cells;

        if (ccobj_.rank() != 0 ) {
            if (ecl_state) {
                // Handle potential exception during MINPV processing
                // Needed because later there is collective communication that will
                // otherwise deadlock
                int success = 1;
                ccobj_.broadcast(&success, 1, 0);
                if (success == 0) {
                    throw std::runtime_error("Error during MINPV processing");
                }
            }
            // Store global grid only on rank 0
            return removed_cells;
        }

        const Opm::EclipseGrid& ecl_grid = *ecl_grid_ptr;
        std::vector<double> coordData = ecl_grid.getCOORD();
        std::vector<int> actnumData = ecl_grid.getACTNUM();

        // Mutable because grdecl::zcorn is non-const.
        auto zcornData = getSanitizedZCORN(ecl_grid, actnumData);

        // Make input struct for processing code.
        grdecl g;
        g.dims[0] = ecl_grid.getNX();
        g.dims[1] = ecl_grid.getNY();
        g.dims[2] = ecl_grid.getNZ();
        g.coord = &coordData[0];
        g.zcorn = &zcornData[0];

        g.actnum = actnumData.empty() ? nullptr : &actnumData[0];
        Opm::MinpvProcessor::Result minpv_result;

        double tolerance_unique_points = 0;
        bool pinchOptionALL = false;
        NNCMaps nnc_cells;

        // Possibly process MINPV and PINCH
        // This even needs to be done if neither of them is specified.
        if (ecl_state ) {
            const size_t cartGridSize = g.dims[0] * g.dims[1] * g.dims[2];
            const auto& fp = ecl_state->fieldProps();
            const auto& permZ = [&fp, cartGridSize](){
                if(fp.has_double("PERMZ")) return fp.get_global_double("PERMZ");
                if(fp.has_double("PERMY")) return fp.get_global_double("PERMY");
                if(fp.has_double("PERMY")) return fp.get_global_double("PERMX");
                // Make this part run without PERM* for some tests
                return std::vector<double>(cartGridSize, 1);
            }();

            try {
                Opm::MinpvProcessor mp(g.dims[0], g.dims[1], g.dims[2]);
                std::vector<double> thickness(cartGridSize);
                for (size_t i = 0; i < cartGridSize; ++i) {
                    thickness[i] = ecl_grid.getCellThickness(i);
                }
                const double z_tolerance = ecl_grid.isPinchActive() ?  ecl_grid.getPinchThresholdThickness() : 0.0;
                const bool nogap = !pinchActive || ecl_grid.getPinchGapMode() ==  Opm::PinchMode::NOGAP;
                const auto& poreVolume = ecl_state->fieldProps().porv(true);
                pinchOptionALL = ecl_grid.getPinchOption() == Opm::PinchMode::ALL;
                const auto& transMult = ecl_state->getTransMult();
                auto multZ =[ &transMult] (int cartindex) {
                    return transMult.getMultiplier(cartindex, ::Opm::FaceDir::ZPlus) *
                        transMult.getMultiplier(cartindex, ::Opm::FaceDir::ZMinus);
                };
                minpv_result = mp.process(thickness, z_tolerance, ecl_grid.getPinchMaxEmptyGap(),
                                          poreVolume, ecl_grid.getMinpvVector(), actnumData, false,
                                          zcornData.data(), nogap, pinchOptionALL,
                                          permZ, multZ, tolerance_unique_points);
                if (!minpv_result.nnc.empty()) {
                    this->zcorn = zcornData;
                }
            }catch(const std::runtime_error& e){
                int success = 0;
                // comminicate failure to others.
                ccobj_.broadcast(&success, 1, 0);
                throw; // rethrow
            }
            int success = 1;
            // communicate success to others
            ccobj_.broadcast(&success, 1, 0);

            // Add PINCH NNCs.
            std::vector<Opm::NNCdata> pinchedNNCs;

            for (const auto& [cell1, cell2] : minpv_result.nnc) {
                nnc_cells[PinchNNC].insert({cell1, cell2});

                if (pinchOptionALL) {
                    auto topIJK = ecl_grid.getIJK(cell1);
                    auto bottomIJK = ecl_grid.getIJK(cell2);
                    std::vector<double> trans_between(bottomIJK[2]-topIJK[2]);
                    std::vector<std::size_t> cells_between;
                    cells_between.reserve(trans_between.size());

                    // first we need to calculate transmissibilities from permeability and geometry
                    auto trans = trans_between.begin();
                    auto bottom_cell_info = ecl_grid.getCellAndBottomCenterNormal(cell1);

                    // We calculate the transmissibilty tran for each intersection between the
                    // two active cells. For each intersection we need geometry information
                    // (distance from cell center to face center, face area and normal) of the
                    // the top (cell_top) and the bottom cell (cell_bottom).
                    // we start at the top (top cell is active) and iterate until the bottom
                    // (bottom cell is active).
                    for (std::size_t cell_top = cell1, cell_bottom = cell_top
                             + ecl_grid.getNX() * ecl_grid.getNY();
                         cell_top < static_cast<std::size_t>(cell2);
                         cell_top = cell_bottom, cell_bottom = cell_top
                             + ecl_grid.getNX() * ecl_grid.getNY(), ++trans) {

                        const auto top_cell_info = bottom_cell_info;
                        bottom_cell_info = ecl_grid.getCellAndBottomCenterNormal(cell_bottom);
                        const auto& cell_center_bottom = std::get<0>(bottom_cell_info);
                        cells_between.push_back(cell_top);

                        auto compute_half_trans =
                            [](const auto& cell_center,
                               const auto& face_center,
                               const auto& area_normal,
                               double perm)
                            {
                                auto half_trans = perm;
                                std::array<double, 3> distance;
                                std::transform(cell_center.begin(), cell_center.end(),
                                               face_center.begin(),
                                               distance.begin(), std::minus<double>());
                                half_trans *= std::abs(std::inner_product(area_normal.begin(),
                                                                          area_normal.end(),
                                                                          distance.begin(), 0.));
                                half_trans /= std::inner_product(distance.begin(),
                                                                 distance.end(),
                                                                 distance.begin(), 0.);
                                return half_trans;
                            };

                        auto half_trans_top = compute_half_trans(std::get<0>(top_cell_info),
                                                                 std::get<1>(top_cell_info),
                                                                 std::get<2>(top_cell_info),
                                                                 permZ[cell_top]);
                        auto half_trans_bottom = compute_half_trans(cell_center_bottom,
                                                                    std::get<1>(top_cell_info),
                                                                    std::get<2>(top_cell_info),
                                                                    permZ[cell_bottom]);

                        if (std::abs(half_trans_top) < 1e-30 || std::abs(half_trans_bottom) < 1e-30)
                            *trans = 0.0;
                        else
                            *trans = 1.0 / (1.0/half_trans_top + 1.0/half_trans_bottom);
                    }

                    // Possibly overwrite with specified TRANZ values.
                    if (fp.tran_active("TRANZ"))
                    {
                        fp.apply_tranz_global(cells_between, trans_between);
                    }

                    // Apply Multipliers

                    // \todo FIXME We assume here that MULTZ does not change in the SCHEDULE or such changes have
                    //       no effect here. This might be wrong and need fixing. Which basically means we need to store
                    //       all the intermediate transmissibilities for the harmonic average and later apply additional
                    //       multipliers
                    const auto& transMult = ecl_state->getTransMult();
                    trans = trans_between.begin();

                    for (std::size_t cell_top = cell1, cell_bottom = cell_top + ecl_grid.getNX() * ecl_grid.getNY();
                         cell_top < static_cast<std::size_t>(cell2);
                         cell_top = cell_bottom, cell_bottom +=  ecl_grid.getNX() * ecl_grid.getNY(), ++trans) {

                        *trans *= transMult.getMultiplier(cell_top, ::Opm::FaceDir::ZPlus) *
                            transMult.getMultiplier(cell_bottom, ::Opm::FaceDir::ZMinus);
                    }

                    //Compute harmonic average over pinched out cells.
                    double average{};
                    bool isZero = false;
                    for(const auto& trans1: trans_between)
                        if (std::abs(trans1) >= 1e-30)
                            average += 1.0 / trans1;
                        else
                            isZero = true;

                    if (isZero)
                        average = 0;
                    else
                        average = 1.0 / average;

                    // Set nnc and transmissibility, last param indicates that this from pinch
                    // It is needed to overwrite transmissibilities instead of adding to existing ones.
                    pinchedNNCs.emplace_back(cell1, cell2, average);
                }
            }

            if (!nnc_cells[PinchNNC].empty()) {
                auto suffix = std::string{(nnc_cells[PinchNNC].size() != 1)? "s" : ""};
                Opm::OpmLog::info(std::to_string(nnc_cells[PinchNNC].size()) + " pinch-out connection" + suffix + " generated");
            }

            // Add explicit NNCs.
            const auto& nncs = ecl_state->getInputNNC();
            for (const auto single_nnc : nncs.input()) {
                // Repeated NNCs will only exist in the map once (repeated
                // insertions have no effect). The code that computes the
                // transmissibilities is responsible for ensuring repeated NNC
                // transmissibilities are added.
                nnc_cells[ExplicitNNC].insert({single_nnc.cell1, single_nnc.cell2});
            }
            // Add the pinch NNCs with transmissibilties due to PINCH option 4 all
            ecl_state->setPinchNNC(std::move(pinchedNNCs));
            ecl_state->prune_global_for_schedule_run();
        }

        // this variable is only required because getCellZvals() needs
        // a coord_t instead of a plain integer pointer...
        coord_t logicalCartesianSize;
        for (int axisIdx = 0; axisIdx < 3; ++axisIdx)
            logicalCartesianSize[axisIdx] = g.dims[axisIdx];

        // Handle zcorn clipping. The g variable points to the data in
        // the clipped_zcorn variable, i.e. clipped_zcorn must remain
        // in scope.
        std::vector<double> clipped_zcorn;
        if (clip_z) {
            double minz_top = 1e100;
            double maxz_bot = -1e100;
            for (int i = 0; i < g.dims[0]; ++i) {
                for (int j = 0; j < g.dims[1]; ++j) {
                    coord_t logicalCartesianCoord;
                    logicalCartesianCoord[0] = i;
                    logicalCartesianCoord[1] = j;

                    logicalCartesianCoord[2] = 0;
                    std::array<double, 8> cellz_bot = getCellZvals(logicalCartesianCoord, logicalCartesianSize, &zcornData[0]);

                    logicalCartesianCoord[2] = g.dims[2] - 1;
                    std::array<double, 8> cellz_top = getCellZvals(logicalCartesianCoord, logicalCartesianSize, &zcornData[0]);

                    for (int dd = 0; dd < 4; ++dd) {
                        minz_top = std::min(cellz_top[dd+4], minz_top);
                        maxz_bot = std::max(cellz_bot[dd], maxz_bot);
                    }
                }
            }
            if (minz_top <= maxz_bot) {
                OPM_THROW(std::runtime_error, "Grid cannot be clipped to a shoe-box (in z): Would be empty afterwards.");
            }
            int num_zcorn = zcornData.size();
            clipped_zcorn.resize(num_zcorn);
            for (int i = 0; i < num_zcorn; ++i) {
                clipped_zcorn[i] = std::max(maxz_bot, std::min(minz_top, g.zcorn[i]));
            }
            g.zcorn = &clipped_zcorn[0];
            this->zcorn = clipped_zcorn;
        }

        if (periodic_extension) {
            // Extend grid periodically with one layer of cells in the (i, j) directions.
            std::vector<double> new_coord;
            std::vector<double> new_zcorn;
            std::vector<int> new_actnum;
            grdecl new_g;
            addOuterCellLayer(g, new_coord, new_zcorn, new_actnum, new_g);
            // Make the grid.
            processEclipseFormat(new_g, ecl_state, nnc_cells, true, turn_normals, pinchActive, tolerance_unique_points);
        } else {
            // Make the grid.
            processEclipseFormat(g, ecl_state, nnc_cells, false, turn_normals, pinchActive, tolerance_unique_points);
        }

        return minpv_result.removed_cells;
    }
#endif // #if HAVE_ECL_INPUT


    enum { NNCFace = -1 };


    /// Read the Eclipse grid format ('.grdecl').
    void CpGridData::processEclipseFormat(const grdecl& input_data,
#if HAVE_ECL_INPUT
                                          Opm::EclipseState* ecl_state,
#endif
                                          NNCMaps& nnc, bool remove_ij_boundary, bool turn_normals,
                                          bool pinchActive,
                                          double tolerance_unique_points)
    {
        if( ccobj_.rank() != 0 )
        {
            OPM_THROW(std::logic_error, "Processing  eclipse file only allowed on rank 0");
        }
        // Process.
#ifdef VERBOSE
        std::cout << "Processing eclipse data." << std::endl;
#endif

        processed_grid output;
        int process_ok;

#if HAVE_ECL_INPUT
        if (ecl_state && ecl_state->aquifer().hasNumericalAquifer()) {
            const auto aquifer_cell_volumes = ecl_state->aquifer().numericalAquifers().aquiferCellVolumes();
            const size_t global_nc = input_data.dims[0] * input_data.dims[1] * input_data.dims[2];
            std::vector<int> is_aquifer_cell(global_nc, 0);
            for ([[maybe_unused]]const auto&[global_index, volume] : aquifer_cell_volumes) {
                is_aquifer_cell[global_index] = 1;
            }
            process_ok = process_grdecl(&input_data, tolerance_unique_points, is_aquifer_cell.data(), &output, pinchActive);
        } else
#endif
        {
            process_ok = process_grdecl(&input_data, tolerance_unique_points, nullptr, &output, pinchActive);
        }

        if (process_ok == 0) {
            OPM_THROW(std::runtime_error,
                      "Failed to build unstructured "
                      "grid from COORD/ZCORN");
        }

        if (remove_ij_boundary) {
            removeOuterCellLayer(output);
            // removeUnusedNodes(output);
        }

#if HAVE_ECL_INPUT
        if (ecl_state) {
            const auto& aquifer = ecl_state->aquifer();
            if (aquifer.hasNumericalAquifer()) {
                const size_t global_nc = input_data.dims[0] * input_data.dims[1] * input_data.dims[2];
                std::vector<int> new_actnum(global_nc, 0);
                for (int i = 0; i < output.number_of_cells; ++i) {
                    new_actnum[output.local_cell_index[i]] = 1;
                }
                const auto& ecl_grid = ecl_state->getInputGrid();
                const auto& fp = ecl_state->fieldProps();
                aquifer.mutableNumericalAquifers().postProcessConnections(ecl_grid, new_actnum);
                const auto& aquifer_nnc = aquifer.numericalAquifers().aquiferConnectionNNCs(ecl_grid, fp);
                // We need to update the nnc in the ecl_state
                ecl_state->appendInputNNC(aquifer_nnc);
                for (const auto single_nnc : aquifer_nnc) {
                    nnc[ExplicitNNC].insert({single_nnc.cell1, single_nnc.cell2});
                }
            }
        }
#endif

        // Move data into the grid's structures.
#ifdef VERBOSE
        std::cout << "Building topology." << std::endl;
#endif
        std::vector<int> face_to_output_face;
        buildTopo(output, nnc, global_cell_, cell_to_face_, face_to_cell_, face_to_point_, cell_to_point_, face_to_output_face);
        std::copy(output.dimensions, output.dimensions + 3, logical_cartesian_size_.begin());

#ifdef VERBOSE
        std::cout << "Building geometry." << std::endl;
#endif
        // here we need the cell volumes based on the active index order
        std::unordered_map<size_t, double> aquifer_cell_volumes_local;
#if HAVE_ECL_INPUT
        if (ecl_state && ecl_state->aquifer().hasNumericalAquifer()) {
            const auto& aquifer_cell_volumes = ecl_state->aquifer().numericalAquifers().aquiferCellVolumes();
            aquifer_cells_.reserve(aquifer_cell_volumes.size());
            for (auto nc = this->global_cell_.size(), i = 0 * nc; i < nc; ++i) {
                auto aquCellPos = aquifer_cell_volumes.find(this->global_cell_[i]);
                if (aquCellPos != aquifer_cell_volumes.end()) {
                    aquifer_cell_volumes_local.emplace(i, aquCellPos->second);
                    aquifer_cells_.push_back(i);
                }
            }
        }
#endif
        std::sort(aquifer_cells_.begin(), aquifer_cells_.end());
        buildGeom(output, cell_to_face_, cell_to_point_, face_to_output_face, aquifer_cell_volumes_local, *(geometry_.geomVector(std::integral_constant<int,0>())),
                  *( geometry_.geomVector(std::integral_constant<int,1>())), geometry_.geomVector(std::integral_constant<int,3>()),
                  face_normals_, turn_normals);

#ifdef VERBOSE
        std::cout << "Assigning face tags." << std::endl;
#endif
        int nf = face_to_output_face.size();
        std::vector<enum face_tag> temp_tags(nf);
        for (int i = 0; i < nf; ++i) {
            const int output_face = face_to_output_face[i];
            if (output_face == -1) {
                temp_tags[i] = NNC_FACE;
            } else {
                temp_tags[i] = output.face_tag[output_face];
            }
        }
        face_tag_.assign(temp_tags.begin(), temp_tags.end());

#ifdef VERBOSE
        std::cout << "Cleaning up." << std::endl;
#endif
        // Clean up the output struct.
        free_processed_grid(&output);

        computeUniqueBoundaryIds();

        if(ccobj_.size()>1)
            populateGlobalCellIndexSet();

        index_set_ = std::make_unique<IndexSet>(cell_to_face_.size(), geomVector<3>().size());

#ifdef VERBOSE
        std::cout << "Done with grid processing." << std::endl;
#endif
    }

    } // end namespace cpgrid



    // ---- Implementation details below ----




    namespace
    {
#if HAVE_ECL_INPUT
        std::vector<double>
        getSanitizedZCORN(const ::Opm::EclipseGrid& ecl_grid,
                          const ::std::vector<int>& actnumData)
        {
            std::vector<double> zcornData = ecl_grid.getZCORN();

            auto repair = ::Opm::UgGridHelpers::RepairZCORN {
                std::move(zcornData), actnumData,
                std::vector<std::size_t>{ ecl_grid.getNX() ,
                                          ecl_grid.getNY() ,
                                          ecl_grid.getNZ() }
            };

            zcornData = repair.destructivelyGrabSanitizedValues();

            if (repair.switchedToDepth()) {
                std::cout << "ZCORN Values Switched from Elevation to "
                          << "Depth (Sign Reversal)\n";
            }

            {
                const auto& statTBB = repair.statTopBelowBottom();

                if (statTBB.cells > std::size_t{0}) {
                    std::cout << "ZCORN Changes From Top Not Below Bottom:\n"
                              << "  - Number of Cells Changed:   "
                              << statTBB.cells << '\n'
                              << "  - Number of Corners Changed: "
                              << statTBB.corners << '\n';
                }
            }

            {
                const auto& statBBLT = repair.statBottomBelowLowerTop();

                if (statBBLT.cells > std::size_t{0}) {
                    std::cout << "ZCORN Changes From Bottom Not Below Lower Top:\n"
                              << "  - Number of Cells Changed:   "
                              << statBBLT.cells << '\n'
                              << "  - Number of Corners Changed: "
                              << statBBLT.corners << '\n';
                }
            }

            return zcornData;
        }
#endif

        typedef std::array<int, 3> coord_t;
        typedef std::array<double, 8> cellz_t;

        cellz_t getCellZvals(const coord_t& c, const coord_t& n, const double* z)
        {
            // cout << c << endl;
            int delta[3] = { 1,
                             2*n[0],
                             4*n[0]*n[1] };
            int ix = 2*(c[0]*delta[0] + c[1]*delta[1] + c[2]*delta[2]);
            // cout << ix << endl;
            cellz_t cellz = {{ z[ix], z[ix + delta[0]],
                               z[ix + delta[1]], z[ix + delta[1] + delta[0]],
                               z[ix + delta[2]], z[ix + delta[2] + delta[0]],
                               z[ix + delta[2] + delta[1]], z[ix + delta[2] + delta[1] + delta[0]] }};
            return cellz;
        }



        void setCellZvals(const coord_t& c, const coord_t& n, double* z, const cellz_t& cellvals)
        {
            int delta[3] = { 1,
                             2*n[0],
                             4*n[0]*n[1] };
            int ix = 2*(c[0]*delta[0] + c[1]*delta[1] + c[2]*delta[2]);
            z[ix]                                  = cellvals[0];
            z[ix + delta[0]]                       = cellvals[1];
            z[ix + delta[1]]                       = cellvals[2];
            z[ix + delta[1] + delta[0]]            = cellvals[3];
            z[ix + delta[2]]                       = cellvals[4];
            z[ix + delta[2] + delta[0]]            = cellvals[5];
            z[ix + delta[2] + delta[1]]            = cellvals[6];
            z[ix + delta[2] + delta[1] + delta[0]] = cellvals[7];
        }

        coord_t indexToIjk(const coord_t& n, const int index)
        {
            coord_t c;
            c[2] = index/(n[0]*n[1]);
            c[1] = (index%(n[0]*n[1]))/n[0];
            c[0] = index%n[0];
            return c;
        }

        void findTopAndBottomZ(const coord_t& n, const std::vector<double>& z, double& zb, double& zt)
        {
            int numperlevel = 4*n[0]*n[1];
            zb = *std::max_element(z.begin(), z.begin() + numperlevel);
            zt = *std::min_element(z.end() - numperlevel, z.end());
        }


#if HAVE_ECL_INPUT
        /// Add an outer cell layer in the (i, j) directions,
        /// repeating the cells on the other side (for periodic
        /// boundary conditions).
        void addOuterCellLayer(const grdecl& original,
                               std::vector<double>& new_coord,
                               std::vector<double>& new_zcorn,
                               std::vector<int>& new_actnum,
                               grdecl& output)
        {
            // Based on periodic_extension.cpp from the old C++ code,
            // with a few changes:
            //  1. We want actnum of the added cells to be true.
            //  2. We do not treat other fields such as PORO, SATNUM etc.
            //     since the grid will be reduced back to its regular
            //     size before those fields are processed.

            OPM_MESSAGE("WARNING: Assuming vertical pillars in a cartesian grid.");

            // Build new-to-old cell index table.
            // First expand in x.
            coord_t n = {{ original.dims[0], original.dims[1], original.dims[2] }};
            std::vector<int> x_new2old;
            x_new2old.reserve((n[0]+2)*n[1]*n[2]);
            for (int kz = 0; kz < n[2]; ++kz) {
                for (int jy = 0; jy < n[1]; ++jy) {
                    int row_ix = kz*n[0]*n[1] + jy*n[0];
                    x_new2old.push_back(row_ix + n[0] - 1);
                    for (int ix = 1; ix < n[0] + 1; ++ix) {
                        x_new2old.push_back(row_ix + ix - 1);
                    }
                    x_new2old.push_back(row_ix);
                }
            }
            // copy(x_new2old.begin(), x_new2old.end(), ostream_iterator<int>(cout, " "));
            // cout << endl;
            // Then expand in y.
            const int num_new_cells = (n[0]+2)*(n[1]+2)*n[2];
            std::vector<int> new2old;
            new2old.reserve(num_new_cells);
            for (int kz = 0; kz < n[2]; ++kz) {
                for (int jy = 0; jy < n[1] + 2; ++jy) {
                    int offset = kz*(n[0] + 2)*n[1] + (jy - 1)*(n[0] + 2);
                    if (jy == 0) {
                        offset = kz*(n[0] + 2)*n[1] + (n[1] - 1)*(n[0] + 2);
                    } else if (jy == n[1] + 1) {
                        offset = kz*(n[0] + 2)*n[1];
                    }
                    for (int ix = 0; ix < n[0] + 2; ++ix) {
                        new2old.push_back(x_new2old[offset + ix]);
                    }
                }
            }
            assert(int(new2old.size()) == num_new_cells);
            // copy(new2old.begin(), new2old.end(), ostream_iterator<int>(cout, " "));
            // cout << endl;
            // On second thought, we should have used a multidimensional array or something...

            // Build new COORD field.
            std::vector<double> coord;
            coord.reserve(6*(n[0] + 3)*(n[1] + 3));
            const double* old_coord = original.coord;
            double dx = old_coord[6] - old_coord[0];
            double dy = old_coord[6*(n[0] + 1) + 1] - old_coord[1];
            double ox = old_coord[0] - dx;
            double oy = old_coord[1] - dy;
            for (int jy = 0; jy < n[1] + 3; ++jy) {
                double y = oy + jy*dy;
                for (int ix = 0; ix < n[0] + 3; ++ix) {
                    double x = ox + ix*dx;
                    coord.push_back(x);
                    coord.push_back(y);
                    coord.push_back(0.0);
                    coord.push_back(x);
                    coord.push_back(y);
                    coord.push_back(1.0);
                }
            }

            // Build new ZCORN field, PERMX, PORO, ACTNUM, SATNUM.
            const double* old_zcorn = original.zcorn;
            const int* old_actnum = original.actnum;
            std::vector<double> zcorn(8*num_new_cells);
            std::vector<int> actnum(num_new_cells);
            coord_t new_n = {{ n[0] + 2, n[1] + 2, n[2] }};
            for (int kz = 0; kz < new_n[2]; ++kz) {
                for (int jy = 0; jy < new_n[1]; ++jy) {
                    for (int ix = 0; ix < new_n[0]; ++ix) {
                        int new_cell_index = ix + jy*(new_n[0]) + kz*(new_n[0])*(new_n[1]);
                        int old_cell_index = new2old[new_cell_index];

                        cellz_t cellvals = getCellZvals(indexToIjk(n, old_cell_index), n, old_zcorn);
                        // cout << new_cell_index << ' ' << old_cell_index << ' ' << cellvals << endl;
                        setCellZvals(indexToIjk(new_n, new_cell_index), new_n, &zcorn[0], cellvals);
                        actnum[new_cell_index] = old_actnum?old_actnum[old_cell_index]:1;
                        if (ix == 0 || ix == new_n[0] - 1
                            || jy == 0 || jy == new_n[1] - 1) {
                            actnum[new_cell_index] = 1;  // This line is changed from the original.
                        }
                    }
                }
            }

            // Clamp z-coord to make shoe box shape
            bool clamp_z = true;
            if (clamp_z) {
                double zb;
                double zt;
                findTopAndBottomZ(new_n, zcorn, zb, zt);
                for (int i = 0; i < int(zcorn.size()); ++i) {
                    zcorn[i] =  std::min(zt, std::max(zb, zcorn[i]));
                }
            }

            // Build output.
            new_coord.swap(coord);
            new_zcorn.swap(zcorn);
            new_actnum.swap(actnum);
            output.dims[0] = new_n[0];
            output.dims[1] = new_n[1];
            output.dims[2] = new_n[2];
            output.coord = &new_coord[0];
            output.zcorn = &new_zcorn[0];
            output.actnum = &new_actnum[0];
        }
#endif




        /// Helper function used by removeOuterCellLayer().
        int newLogCartFromOld(const int idx, const int dim[3])
        {
            // Compute old (i, j, k).
            const int Nx = dim[0];
            const int Ny = dim[1];
            const int NxNy = Nx*Ny;
            int k = idx/NxNy;
            // if (k <= 0 || k >= dim[2] - 1) return -1;
            int j = (idx - NxNy*k)/Nx;
            if (j <= 0 || j >= Ny - 1) return -1;
            int i = idx - Nx*j - Nx*Ny*k;
            if (i <= 0 || i >= Nx - 1) return -1;
            // return (Nx - 2)*(Ny - 2)*(k - 1) + (Nx - 2)*(j - 1) + (i - 1);
            return (Nx - 2)*(Ny - 2)*k + (Nx - 2)*(j - 1) + (i - 1);
        }

        /// Removes all (i, j) boundary cells from a grid.
        void removeOuterCellLayer(processed_grid& grid)
        {
            // Remove outer cells as follows:
            //   1. Build a new local_cell_index (in a new variable), compute new number_of_cells.
            //   2. Build the inverse lookup: From old logical cartesian to new cell indices.
            //   3. Modify face_neighbours by replacing each entry by its new cell index (or -1).
            //   4. Modify dimensions[], number_of_cells and replace local_cell_index.
            // After this, we still have the same number of faces, it's just that some of them may
            // have only (-1, -1) as neighbours.

            // Part 1 and 2 in one pass.
            std::vector<int> new_index_to_new_lcart;
            new_index_to_new_lcart.reserve(grid.number_of_cells); // A little too large, but no problem.
            int num_old_lcart = grid.dimensions[0]*grid.dimensions[1]*grid.dimensions[2];
            std::vector<int> old_lcart_to_new_index(num_old_lcart, -1);
            for (int i = 0; i < grid.number_of_cells; ++i) {
                int old_lcart = grid.local_cell_index[i];
                int new_lcart = newLogCartFromOld(old_lcart, grid.dimensions);
                if (new_lcart != -1) {
                    old_lcart_to_new_index[old_lcart] = new_index_to_new_lcart.size();
                    new_index_to_new_lcart.push_back(new_lcart);
                } else {
                    old_lcart_to_new_index[old_lcart] = -1;
                }
            }

            // Part 3, modfying the face->cell connections.
            for (grid_size_t i = 0; i < 2*grid.number_of_faces; ++i) {
                int old_index = grid.face_neighbors[i];
                if (old_index != -1) {
                    int old_lcart = grid.local_cell_index[old_index];
                    int new_index = old_lcart_to_new_index[old_lcart];
                    grid.face_neighbors[i] = new_index; // May be -1, if cell is to be removed.
                }
            }

            // Part 4, modifying the other output data.
            grid.dimensions[0] = grid.dimensions[0] - 2;
            grid.dimensions[1] = grid.dimensions[1] - 2;
            grid.dimensions[2] = grid.dimensions[2] - 2;
            grid.number_of_cells = new_index_to_new_lcart.size();
            std::copy(new_index_to_new_lcart.begin(), new_index_to_new_lcart.end(), grid.local_cell_index);
        }





        /*

        // NOTE: This function has been commented out as it is not in use, and therefore
        //       generates warnings. It has not been deleted, as it should eventually be
        //       used, what is missing is thorough testing of the method.

        void removeUnusedNodes(processed_grid& grid)
        {
            // Nodes are considered unused if they are unreachable from a cell.

            // The following data are modified by this function:
            //    grid.face_nodes
            //    grid.number_of_nodes
            //    grid.node_coordinates[]

            // The following caveats apply:
            //   - grid.face_nodes will contain -1 where nodes have been removed
            //     (this will only happen for faces not neighbouring any cells)
            //   - grid.number_of_nodes_on_pillars is unchanged (and therefore wrong)

            // Remove unused nodes in three steps:
            //
            //   1. Build the old_to_new node index array (-1 means removed).
            //
            //      a) Initially, all are considered unused. We first signify usage
            //         by changing to a 0 the entries corresponding to reachable nodes.
            std::vector<int> old_to_new(grid.number_of_nodes, -1);
            for (int face = 0; face < grid.number_of_faces; ++face) {
                if (grid.face_neighbors[2*face] != -1 || grid.face_neighbors[2*face + 1] != -1) {
                    // Face is reachable
                    for (int ii = grid.face_ptr[face]; ii < grid.face_ptr[face + 1]; ++ii) {
                        int node = grid.face_nodes[ii];
                        old_to_new[node] = 0;
                    }
                }
            }
            //      b) Set the new indices by simple array compression.
            int nodecount = 0;
            for (int node = 0; node < grid.number_of_nodes; ++node) {
                if (old_to_new[node] != -1) {
                    assert(old_to_new[node] == 0);
                    old_to_new[node] = nodecount;
                    ++nodecount;
                }
            }

            //   2. Use old_to_new to transform grid.face_nodes and grid.node_coordinates[].
            for (int fnode = 0; fnode < grid.face_ptr[grid.number_of_faces]; ++fnode) {
                int old = grid.face_nodes[fnode];
                grid.face_nodes[fnode] = old_to_new[old];
            }
            double* nc = grid.node_coordinates;
            for (int node = 0; node < grid.number_of_nodes; ++node) {
                int newidx = old_to_new[node];
                if (newidx != -1) {
                    nc[3*newidx]     = nc[3*node];
                    nc[3*newidx + 1] = nc[3*node + 1];
                    nc[3*newidx + 2] = nc[3*node + 2];
                }
            }

            //   3. Set grid.number_of_nodes.
            grid.number_of_nodes = nodecount;
        }
        */




        std::vector<int> createGlobalToLocal(const processed_grid& output,
                                             const std::vector<int>& global_cell)
        {
            std::vector<int> global_to_local;
            int cart_size = 1;
            const int num_dims = sizeof(output.dimensions)/sizeof(*output.dimensions);
            for (int idx = 0; idx < num_dims ; ++idx) {
                cart_size *= output.dimensions[idx];
            }
            // create the inverse map of global_cell
            global_to_local.resize(cart_size, -1);
            for (size_t idx = 0; idx < global_cell.size(); ++idx) {
                global_to_local[global_cell[idx]] = idx;
            }
            return global_to_local;
        }





        NNCMap filterNNCs(const processed_grid& output,
                          const NNCMap& nnc,
                          const std::vector<int>& global_to_local)
        {
            NNCMap filtered_nnc;
            const int num_faces = output.number_of_faces;
            std::vector<std::pair<int, int>> face_cells(num_faces);
            // Sort all face->cell mappings so that lowest cell number comes first.
            for (int f = 0; f < num_faces; ++f) {
                const int c1 = output.face_neighbors[2*f];
                const int c2 = output.face_neighbors[2*f + 1];
                if (c1 < c2) {
                    face_cells[f] = { c1, c2 };
                } else {
                    face_cells[f] = { c2, c1 };
                }
            }
            // Sort face->cell mappings according to first, then second cell.
            std::sort(face_cells.begin(), face_cells.end());
            // For each nnc, add it to filtered_nnc only if not found in face->cell mappings.
            for (const auto& nncpair : nnc) {
                if (nncpair.first < 0 || nncpair.second < 0 ||
                    nncpair.first >= static_cast<int>(global_to_local.size()) ||
                    nncpair.second >= static_cast<int>(global_to_local.size())) {
                    Opm::OpmLog::warning("nnc_invalid", "NNC connection requested between invalid cells.");
                    continue;
                }
                const int c1 = global_to_local[nncpair.first];
                const int c2 = global_to_local[nncpair.second];
                if (c1 < 0 || c2 < 0) {
                    Opm::OpmLog::warning("nnc_inactive", "NNC connection requested between inactive cells.");
                    continue;
                }
                const auto beg = std::lower_bound(face_cells.begin(), face_cells.end(), std::make_pair(c1, -1));
                const auto end = std::find_if(beg, face_cells.end(), [c1](const std::pair<int, int>& p){ return p.first > c1; });
                const auto it = std::find_if(beg, end, [c2](const std::pair<int, int>& p){ return p.second == c2; });
                if (it == end) {
                    // The connection (c1, c2) was not found in the face->cell mapping.
                    filtered_nnc.insert(nncpair);
                }
            }
            return filtered_nnc;
        }





        void buildFaceToCellNNC(const processed_grid& output,
                                const NNCMap& nnc,
                                const std::vector<int>& global_to_local,
                                cpgrid::OrientedEntityTable<1, 0>& f2c,
                                std::vector<int>& face_to_output_face)
        {
            // Add nnc connections first, to preserve the
            // property that top and bottom faces come last
            // (per cell) in the cell-face mapping.
            // A complicating factor is that an NNC could be specified
            // for a connection that is also appearing as a face in
            // the geometry-based grid processing. In that case we
            // should ensure we do not add it twice, and therefore we
            // filter them out first.
            NNCMap filtered_nnc = filterNNCs(output, nnc, global_to_local);
            cpgrid::EntityRep<0> cells[2];
            for (const auto& nncpair : filtered_nnc) {
                const int c1 = global_to_local[nncpair.first];
                const int c2 = global_to_local[nncpair.second];
                cells[0].setValue(c1, true);
                cells[1].setValue(c2, false);
                std::sort(cells, cells + 2);
                f2c.appendRow(cells, cells + 2);
                face_to_output_face.push_back(cpgrid::NNCFace);
            }
        }





        void buildFaceToCell(const processed_grid& output,
                             const NNCMaps& nnc,
                             const std::vector<int>& global_cell,
                             cpgrid::OrientedEntityTable<1, 0>& f2c,
                             std::vector<int>& face_to_output_face)
        {
            std::vector<int> global_to_local;
            if (!nnc[ExplicitNNC].empty() || !nnc[PinchNNC].empty())
                global_to_local = createGlobalToLocal(output, global_cell);

            f2c.clear();
            face_to_output_face.clear();
            // Reserve to save allocation time. True required size may be smaller.
            face_to_output_face.reserve(output.number_of_faces + nnc[ExplicitNNC].size());
            if (!nnc[ExplicitNNC].empty()) {
                buildFaceToCellNNC(output, nnc[ExplicitNNC], global_to_local, f2c, face_to_output_face);
            }
            int nf = output.number_of_faces;
            cpgrid::EntityRep<0> cells[2];
            // int next_skip_zmin_face = -1; // Not currently used, see comments further down.
            for (int i = 0; i < nf; ++i) {
                const int* fnc = output.face_neighbors + 2*i;
                int cellcount = 0;
                if (fnc[0] != -1) {
                    cells[cellcount].setValue(fnc[0], true);
                    ++cellcount;
                }
                if (fnc[1] != -1) {
                    cells[cellcount].setValue(fnc[1], false);
                    ++cellcount;
                }
                if (output.face_tag[i] == K_FACE ) {
                    if (fnc[1] == -1 ) {
                        // Add the NNC created from the minpv processs
                        // at the bottom of the cell.
                        if(fnc[0] != -1)
                        {
                            auto it = nnc[PinchNNC].lower_bound({global_cell[fnc[0]], 0});
                            if (it != nnc[PinchNNC].end() && it->first == global_cell[fnc[0]]) {
                                const int other_cell = global_to_local[it->second];
                                cells[cellcount].setValue(other_cell, false);
                                ++cellcount;
                                // Now we must ensure that the face connecting
                                // the second cell to the boundary is skipped,
                                // it is considered replaced by the connection
                                // added here.
                                // Note: not done anymore, see comment below.
                                // next_skip_zmin_face = other_cell;
                            }
                        }
                    } else if (fnc[0] == -1) {
                        // This block has been disabled, as removing the original top face
                        // of a cell changed the cell corners, thereby changing TRANX and TRANY
                        // for example. With the block disabled, cells that receive a PINCH
                        // connection face from above will retain their original top face, in addition
                        // to the new one. A typical such cell will therefore have 7 faces rather than
                        // the usual 6 (away from faults).
                        /*
                        // Check if this is a cell we should skip.
                        if (fnc[1] == next_skip_zmin_face) {
                            next_skip_zmin_face = -1;
                            continue; // Do not add this face to f2c.
                        }
                        */
                    }
                }

                // Assertation below is no longer true, due to periodic_extension etc.
                // Instead, the appendRow() is put inside an if test.
                // assert(cellcount == 1 || cellcount == 2);
                if (cellcount > 0) {
                    std::sort(cells, cells + cellcount);
                    f2c.appendRow(cells, cells + cellcount);
                    face_to_output_face.push_back(i);
                }
            }
        }





        void buildTopo(const processed_grid& output,
                       const NNCMaps& nnc,
                       std::vector<int>& global_cell,
                       cpgrid::OrientedEntityTable<0, 1>& c2f,
                       cpgrid::OrientedEntityTable<1, 0>& f2c,
                       Opm::SparseTable<int>& f2p,
                       std::vector<std::array<int,8> >& c2p,
                       std::vector<int>& face_to_output_face)
        {
            // Map local to global cell index.
            global_cell.assign(output.local_cell_index,
                               output.local_cell_index + output.number_of_cells);

            // Build face to cell mapping.
            buildFaceToCell(output, nnc, global_cell, f2c, face_to_output_face);
            int num_faces = f2c.size();

            // Build cell to face.
            f2c.makeInverseRelation(c2f);
            int num_cells = c2f.size();

            // Build face to point
            const int* fn = output.face_nodes;
            const grid_size_t* fp = output.face_ptr;
            for (int face = 0; face < num_faces; ++face) {
                int output_face = face_to_output_face[face];
                if (output_face == cpgrid::NNCFace) {
                    // Add an empty row, i.e. no points are associated with the NNC face.
                    f2p.appendRow(fn, fn);
                } else {
                    f2p.appendRow(fn + fp[output_face], fn + fp[output_face+1]);
                }
            }

            // Build cell to point
            c2p.clear();
            c2p.reserve(num_cells);
            for (int i = 0; i < num_cells; ++i) {
                cpgrid::OrientedEntityTable<0, 1>::row_type cf = c2f[cpgrid::EntityRep<0>(i, true)];
                // We know that the bottom and top faces come last.
                int numf = cf.size();
                int bot_face = face_to_output_face[cf[numf - 2].index()];
                assert(output.face_tag[bot_face] == K_FACE);
                int bfbegin = output.face_ptr[bot_face];
                assert(output.face_ptr[bot_face + 1] - bfbegin == 4);
                int top_face = face_to_output_face[cf[numf - 1].index()];
                assert(output.face_tag[top_face] == K_FACE);
                int tfbegin = output.face_ptr[top_face];
                assert(output.face_ptr[top_face + 1] - tfbegin == 4);
                // We want the corners in 'x fastest, then y, then z' order,
                // so we need to take the face_nodes in noncyclic order: 0 1 3 2.
                std::array<int,8> corners = {{ output.face_nodes[bfbegin],
                                               output.face_nodes[bfbegin + 1],
                                               output.face_nodes[bfbegin + 3],
                                               output.face_nodes[bfbegin + 2],
                                               output.face_nodes[tfbegin],
                                               output.face_nodes[tfbegin + 1],
                                               output.face_nodes[tfbegin + 3],
                                               output.face_nodes[tfbegin + 2] }};
                c2p.push_back(corners);
            }
#ifndef NDEBUG
#ifdef VERBOSE
            std::cout << "Doing extra topology integrity check." << std::endl;
#endif
            cpgrid::OrientedEntityTable<1, 0> f2c_again;
            c2f.makeInverseRelation(f2c_again);
            assert(f2c == f2c_again);
#endif
        }

        /// Encapsulate a vector<T>, and a permutation array used for access.
        template <typename T>
        class IndirectArray
        {
        public:
            IndirectArray(const std::vector<T>& data, const int* beg, const int* end)
                : data_(data), beg_(beg), end_(end)
            {
            }
            const T& operator[](int index) const
            {
                assert(index >= 0 && index < size());
                return data_[beg_[index]];
            }
            int size() const
            {
                return end_ - beg_;
            }
            typedef T value_type;
        private:
            const std::vector<T>& data_;
            const int* beg_;
            const int* end_;
        };


        // Helper template for making Geometry objects.
        // The generic one is suitable for dim == 0 (vertices).
        template <int dim>
        struct MakeGeometry
        {
            cpgrid::Geometry<dim, 3> operator()(const FieldVector<double, 3>& pos)
            {
                return cpgrid::Geometry<dim, 3>(pos);
            }
        };

        // Specialization for cells.
        template <>
        struct MakeGeometry<3>
        {
            std::shared_ptr<const cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> allcorners_;
            MakeGeometry(std::shared_ptr<const cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> allcorners)
                : allcorners_(std::move(allcorners))
            {
            }
            cpgrid::Geometry<3, 3> operator()(const FieldVector<double, 3>& pos,
                                                      double vol,
                                                      const std::array<int,8>& corner_indices)
            {
                return cpgrid::Geometry<3, 3>(pos, vol, allcorners_, &corner_indices[0]);
            }
        };

        // Specialization for intersections.
        template <>
        struct MakeGeometry<2>
        {
            cpgrid::Geometry<2, 3> operator()(const FieldVector<double, 3>& pos, double vol)
            {
                return cpgrid::Geometry<2, 3>(pos, vol);
            }
        };



        void buildGeom(const processed_grid& output,
                       const cpgrid::OrientedEntityTable<0, 1>& c2f,
                       const std::vector<std::array<int,8> >& c2p,
                       const std::vector<int>& face_to_output_face,
                       const std::unordered_map<size_t, double>& aquifer_cell_volumes,
                       cpgrid::EntityVariable<cpgrid::Geometry<3, 3>, 0>& cell_geom,
                       cpgrid::EntityVariable<cpgrid::Geometry<2, 3>, 1>& face_geom,
                       std::shared_ptr<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> point_geom_ptr,
                       cpgrid::SignedEntityVariable<FieldVector<double, 3>, 1>& normals,
                       bool turn_normals)
        {
            typedef FieldVector<double, 3> point_t;
            std::vector<point_t> points;
            std::vector<point_t> face_normals;
            std::vector<point_t> face_centroids;
            std::vector<double>  face_areas;
            std::vector<point_t> cell_centroids;
            std::vector<double>  cell_volumes;
            auto& point_geom = *point_geom_ptr;
            using namespace GeometryHelpers;
#ifdef VERBOSE
            Opm::time::StopWatch clock;
            clock.start();
#endif
            // Get the points.
            int np = output.number_of_nodes;
            points.clear();
            points.reserve(np);
            for (int i = 0; i < np; ++i) {
                // \TODO add a convenience explicit constructor
                // for FieldVector taking an iterator.
                point_t pt;
                for (int dd = 0; dd < 3; ++dd) {
                    pt[dd] = output.node_coordinates[3*i + dd];
                }
                points.push_back(pt);
            }
#ifdef VERBOSE
            std::cout << "Points:             " << clock.secsSinceLast() << std::endl;
#endif

            // Get the face data.
            // \TODO Both the face and (especially) the cell section
            // is not very efficient. It could be rewritten easily
            // (focus on the polygonCellXXX methods first).
            // \TODO Use exact geometry instead of these approximations.
            int nf = face_to_output_face.size();
            const int* fn = output.face_nodes;
            const grid_size_t* fp = output.face_ptr;
            for (int face = 0; face < nf; ++face) {
                // Computations in this loop could be speeded up
                // by doing more of them simultaneously.
                int output_face = face_to_output_face[face];
                if (output_face == cpgrid::NNCFace) {
                    // NNC faces are purely topological constructs,
                    // and do not have any embedded geometry.
                    // However, since the ewoms code will multiply and
                    // divide by the face area even if not necessary
                    // for the cell-centered FV discretization (because
                    // it wants to deal with velocities rather than fluxes),
                    // we have to set the areas to 1 to avoid trouble.
                    face_normals.push_back({-1e100, -1e100, -1e100});
                    face_centroids.push_back({-1e100, -1e100, -1e100});
                    face_areas.push_back(1.0);
                } else {
                    IndirectArray<point_t> face_pts(points, fn + fp[output_face], fn + fp[output_face+1]);
                    point_t avg = average(face_pts);
                    point_t centroid = polygonCentroid(face_pts, avg);
                    point_t normal = polygonNormal(face_pts, centroid);
                    double area = polygonArea(face_pts, centroid);
                    face_normals.push_back(normal);
                    face_centroids.push_back(centroid);
                    face_areas.push_back(area);
                }
            }
#ifdef VERBOSE
            std::cout << "Faces:              " << clock.secsSinceLast() << std::endl;
#endif
            // Get the cell data.
            int nc = output.number_of_cells;
            std::vector<int> face_indices;
            for (int cell = 0; cell < nc; ++cell) {
                cpgrid::EntityRep<0> cell_ent(cell, true);
                cpgrid::OrientedEntityTable<0, 1>::row_type cf = c2f[cell_ent];
                face_indices.clear();
                for (int local_index = 0; local_index < cf.size(); ++local_index) {
                    if (face_to_output_face[cf[local_index].index()] != cpgrid::NNCFace) {
                        face_indices.push_back(cf[local_index].index());
                    }
                }
                IndirectArray<point_t> cell_pts(face_centroids, &face_indices[0], &face_indices[0] + face_indices.size());
                point_t cell_avg = average(cell_pts);
                point_t cell_centroid(0.0);
                double tot_cell_vol = 0.0;
                for (int local_index = 0; local_index < cf.size(); ++local_index) {
                    int face = cf[local_index].index();
                    int output_face = face_to_output_face[face];
                    if (output_face == cpgrid::NNCFace) {
                        // Skip NNC face, do not contribute to cell geometry.
                        continue;
                    }
                    IndirectArray<point_t> face_pts(points, fn + fp[output_face], fn + fp[output_face+1]);
                    double small_vol = polygonCellVolume(face_pts, face_centroids[face], cell_avg);
                    tot_cell_vol += small_vol;
                    point_t face_contrib = polygonCellCentroid(face_pts, face_centroids[face], cell_avg);
                    face_contrib *= small_vol;
                    cell_centroid += face_contrib;
                }
                // TODO: for the cell with zero volumes, data in cell_centroid is NaN already. We need to
                // find solution to deal with those cells to get sensible cell_centroid.
                if (tot_cell_vol > 0.) {
                    cell_centroid /= tot_cell_vol;
                }
// #define HACK_CELL_CENTROIDS     // when this is defined, you get the average of top and bottom face centroids.
#ifdef HACK_CELL_CENTROIDS
                int numf = cf.size();
                cell_centroid = face_centroids[face_indices[numf - 2]];
                cell_centroid += face_centroids[face_indices[numf - 1]];
                cell_centroid *= 0.5;
#endif
                cell_centroids.push_back(cell_centroid);
                cell_volumes.push_back(tot_cell_vol);
            }
            // update the volumes of numerical aquifer cells
            for (const auto& [index, volume] : aquifer_cell_volumes) {
                cell_volumes[index] = volume;
            }



#ifdef VERBOSE
            std::cout << "Cells:              " << clock.secsSinceLast() << std::endl;
#endif


            // \TODO Fix the code below , as it is:
            // A) wasteful
            // B) wordy,
            // C) slow
            // D) copied from readSintefLegacyFormat.cpp
            // Points
            // Points
            point_geom.reserve(points.size());
            MakeGeometry<0> mpointg;
            std::transform(points.begin(), points.end(),
                           std::back_inserter(point_geom), mpointg);
            // Cells
            cell_geom.reserve(nc);
            MakeGeometry<3> mcellg(point_geom_ptr);
//             std::transform(cell_centroids.begin(), cell_centroids.end(),
//                            cell_volumes.begin(),
//                            std::back_inserter(cell_geom, mcellg);
            for (int c = 0;  c < nc; ++c) {
                cell_geom.push_back(mcellg(cell_centroids[c], cell_volumes[c], c2p[c]));
            }
            // Faces
            face_geom.reserve(face_centroids.size());
            MakeGeometry<2> mfaceg;
            std::transform(face_centroids.begin(), face_centroids.end(),
                           face_areas.begin(),
                           std::back_inserter(face_geom), mfaceg);
#ifdef VERBOSE
            std::cout << "Transforms/copies:  " << clock.secsSinceLast() << std::endl;
#endif

            // The final, combined object (yes, a lot of copying goes on here).
            if (turn_normals) {
                int num_normals = face_normals.size();
                for (int i = 0; i < num_normals; ++i) {
                    face_normals[i] *= -1.0;
                }
            }
            normals.assign(face_normals.begin(), face_normals.end());
#ifdef VERBOSE
            std::cout << "Final construction: " << clock.secsSinceLast() << std::endl;
#endif
        }
    } // anon namespace
} // namespace Dune

