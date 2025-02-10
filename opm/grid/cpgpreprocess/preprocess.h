/*===========================================================================
//
// File: preprocess.h
//
// Created: Fri Jun 19 08:43:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//==========================================================================*/

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

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

#ifndef OPM_PREPROCESS_HEADER
#define OPM_PREPROCESS_HEADER

/**
 * \file
 * Low-level corner-point processing routines and supporting data structures.
 *
 * User code should typically employ higher-level routines such as
 * create_grid_cornerpoint() in order to construct fully formed UnstructuredGrid
 * data structures from a corner-point specification. Incidentally, the routines
 * provided by this module are used to implement function
 * create_grid_cornerpoint().
 */

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * Raw corner-point specification of a particular geological model.
     */
    struct grdecl {
        int           dims[3]; /**< Cartesian box dimensions. */
        const double *coord;   /**< Pillar end-points. */
        const double *zcorn;   /**< Corner-point depths. */
        const int    *actnum;  /**< Explicit "active" map.  May be NULL.*/
    };

    /**
     * Connection taxonomy.
     */
    enum face_tag {
        I_FACE,        /**< Connection topologically normal to J-K plane. */
        J_FACE,        /**< Connection topologically normal to I-K plane. */
        K_FACE,        /**< Connection topologically normal to I-J plane. */
        NNC_FACE       /**< Arbitrary non-neighbouring connection. */
    };

    /**
     * Result structure representing minimal derived topology and geometry of
     * a geological model in corner-point format.
     */
    struct processed_grid {
        int m; /**< Upper bound on "number_of_faces".  For internal use in
                    function process_grid()'s memory management. */
        int n; /**< Upper bound on "number_of_nodes".  For internal use in
                    function process_grid()'s memory management. */

        int    dimensions[3];     /**< Cartesian box dimensions. */

        unsigned number_of_faces; /**< Total number of unique grid faces
                                   * (i.e., connections). */
        int    *face_nodes;     /**< Node (vertex) numbers of each face,
                                 * stored sequentially. */
        unsigned int *face_ptr; /**< Start position for each face's
                                 * `face_nodes'. */
        int    *face_neighbors;   /**< Global cell numbers.  Two elements
                                   * per face, stored sequentially. */
        enum face_tag *face_tag;  /**< Classification of grid's individual
                                   * connections (faces). */

        /* Cell-to-face mapping, the transpose of the face-to-cell mapping,
         * needed for edge conformal processing */
        int    *cell_facePos;   /**< Start position for each cell's
                                 * `cell_faces'. */

        int    *cell_faces;     /**< Face numbers of each cell, stored
                                 * sequentially. Indexed by cell_facePos. */

        int    number_of_nodes;   /**< Number of unique grid vertices. */
        int    number_of_nodes_on_pillars; /**< Total number of unique cell
                                            * vertices that lie on
                                            * pillars. */
        double *node_coordinates; /**< Vertex coordinates.  Three doubles
                                   * (\f$x\f$, \f$y\f$, \f$z\f$) per vertex,
                                   * stored sequentially. */

        int    number_of_cells;   /**< Number of active grid cells. */
        int    *local_cell_index; /**< Deceptively named local-to-global cell
                                       index mapping. */
    };

    /**
     * Construct a prototypical grid representation from a corner-point
     * specification.
     *
     * Pinched cells will be removed irrespective of any explicit "active" map
     * in the geological model input specification. On input, the result
     * structure "out" must point to a valid management structure. In other
     * words, the result structure must point to a region of memory that is
     * typically backed by automatic or allocated (dynamic) storage duration.
     *
     * @param[in] pinchActive Whether cells with zero volume should be
     *                    pinched out and neighboring cells should be
     *                    connected.
     *
     * @param[in] edge_conformal Whether or not to create an edge-conformal
     *                    grid.  This is an experimental feature, aimed at
     *                    supporting geo-mechanical workflows, that should
     *                    typically not be used in production runs of
     *                    traditional reservoir simulations.  Non-zero to
     *                    enable edge-conformal processing, zero to disable
     *                    this mode.
     *
     * @param[in]     tol Absolute tolerance of node-coincidence.
     *
     * @param[in]     g   Corner-point specification.  If "actnum" is NULL, then
     *                    the specification is interpreted as if all cells are
     *                    initially active.
     *
     * @param[in,out] out Minimal grid representation featuring face-to-cell
     *                    neighbourship definition, vertex geometry, face's
     *                    constituent vertices, and local-to-global cell
     *                    mapping.
     *
     * @return One (1, true) if grid successfully generated, zero (0, false)
     * otherwise.
     */
    int process_grdecl(int                    pinchActive,
                       int                    edge_conformal,
                       double                 tol,
                       const struct grdecl   *g,
                       const int             *is_aquifer_cell,
                       struct processed_grid *out);

    /**
     * Release memory resources acquired in previous grid processing using
     * function process_grdecl().
     *
     * Note: This function releases the resources associated to the individual
     * fields of the processed_grid, but does not free() the structure itself.
     *
     * @param[in,out] g Prototypical grid representation obtained in an earlier
     *                  call to function process_grdecl().
     */
    void free_processed_grid(struct processed_grid *g);

    /**
     * Populate cell-to-face mapping.
     *
     * @param[in,out] grid On input, the neighbourship structure from
     * process_grdecl(), with the 'edge_conformal' flag set to true.  On
     * output, the same structure, but with added cell-to-face mappings.
     */
    void add_cells(struct processed_grid *grid);

#ifdef __cplusplus
}
#endif

#endif /* OPM_PREPROCESS_HEADER */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
