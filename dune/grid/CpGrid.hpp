//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri May 29 20:26:36 2009
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

#ifndef OPM_CPGRID_HEADER
#define OPM_CPGRID_HEADER

#include <string>
#include <map>
#include <array>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridenums.hh>

#include "cpgrid/CpGridData.hpp"
#include "cpgrid/Intersection.hpp"
#include "cpgrid/Entity.hpp"
#include "cpgrid/Geometry.hpp"
#include "cpgrid/Iterators.hpp"
#include "cpgrid/Indexsets.hpp"
#include "cpgrid/DefaultGeometryPolicy.hpp"
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <iostream>


namespace Opm {
    class EclipseGridParser;
    namespace parameter {
	class ParameterGroup;
    }
}

namespace Dune
{

    class CpGrid;

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridTraits
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridTraits
    {
	/// \brief The type that implements the grid.
	typedef CpGrid Grid;

	/// \brief The type of the intersection at the leafs of the grid.
	typedef cpgrid::Intersection LeafIntersection;
	/// \brief The type of the intersection at the levels of the grid.
	typedef cpgrid::Intersection LevelIntersection;
	/// \brief The type of the intersection iterator at the leafs of the grid.
	typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
	/// \brief The type of the intersection iterator at the levels of the grid.
	typedef cpgrid::IntersectionIterator LevelIntersectionIterator;

	/// \brief The type of the  hierarchic iterator.
	typedef cpgrid::HierarchicIterator HierarchicIterator;

	/// \brief Traits associated with a specific codim.
	/// \tparam cd The codimension.
	template <int cd>
	struct Codim
	{
	    /// \brief The type of the geometry associated with the entity.
	    /// IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
 	    typedef cpgrid::Geometry<3-cd, 3> Geometry;
	    //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> Geometry;
	    /// \brief The type of the local geometry associated with the entity.
 	    typedef cpgrid::Geometry<3-cd, 3> LocalGeometry;
	    //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> LocalGeometry;
	    /// \brief The type of the entity.
	    typedef cpgrid::Entity<cd> Entity;

	    /// \brief The type of the iterator over all level entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition> LevelIterator;

	    /// \brief The type of the iterator over all leaf entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition> LeafIterator;

	    /// \brief The type of the entity pointer for entities of this codim.
	    typedef cpgrid::EntityPointer<cd> EntityPointer;

	    /// \brief The type of the entity pointer for entities of this codim.
	    typedef cpgrid::EntityPointer<cd> EntitySeed;

	    /// \brief Traits associated with a specific grid partition type.
	    /// \tparam pitype The type of the grid partition.
	    template <PartitionIteratorType pitype>
	    struct Partition
	    {
		/// \brief The type of the iterator over the level entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype> LevelIterator;
		/// \brief The type of the iterator over the leaf entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype> LeafIterator;
	    };
	};

	/// \brief Traits associated with a specific grid partition type.
	/// \tparam pitype The type of the grid partition.
	template <PartitionIteratorType pitype>
	struct Partition
	{
	    /// \brief The type of the level grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LevelGridView;
	    typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid, pitype> > LevelGridView;

	    /// \brief The type of the leaf grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LeafGridView;
	    typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid, pitype> > LeafGridView;
	};

	/// \brief The type of the level index set.
	typedef cpgrid::IndexSet LevelIndexSet;
	/// \brief The type of the leaf index set.
	typedef cpgrid::IndexSet LeafIndexSet;
	/// \brief The type of the global id set.
	typedef cpgrid::GlobalIdSet GlobalIdSet;
	/// \brief The type of the local id set.
	typedef cpgrid::IdSet LocalIdSet;

	/// \brief The type of the collective communication.

    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    };



    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridFamily
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridFamily
    {
	typedef CpGridTraits Traits;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGrid
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> provides \ref Dune::Grid </em>]
    class CpGrid
	: public GridDefaultImplementation<3, 3, double, CpGridFamily >
    {
    public:

	// --- Typedefs ---


	/// Family typedef, why is this not defined by Grid<>?
	typedef CpGridFamily GridFamily;


	// --- Methods ---


	/// Default constructor
	CpGrid();

        ~CpGrid();
 
        /// Initialize the grid from parameters.
	void init(const Opm::parameter::ParameterGroup& param);


	/// Read the Sintef legacy grid format ('topogeom').
	/// \param grid_prefix the grid name, such that topology is
	/// found in <grid_prefix>-topo.dat etc.
	void readSintefLegacyFormat(const std::string& grid_prefix);
                

	/// Write the Sintef legacy grid format ('topogeom').
	/// \param grid_prefix the grid name, such that topology will be
	/// found in <grid_prefix>-topo.dat etc.
        void writeSintefLegacyFormat(const std::string& grid_prefix) const;
        

	/// Read the Eclipse grid format ('grdecl').
	/// \param filename the name of the file to read.
	/// \param z_tolerance points along a pillar that are closer together in z
	///        coordinate than this parameter, will be replaced by a single point.
	/// \param periodic_extension if true, the grid will be (possibly) refined, so that
	///        intersections/faces along i and j boundaries will match those on the other
	///        side. That is, i- faces will match i+ faces etc.
	void readEclipseFormat(const std::string& filename, double z_tolerance, bool periodic_extension, bool turn_normals = false);

        /// Read the Eclipse grid format ('grdecl').
        /// \param input_data the data contained in a parser object.
        /// \param z_tolerance points along a pillar that are closer together in z
        ///        coordinate than this parameter, will be replaced by a single point.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        void processEclipseFormat(const Opm::EclipseGridParser& input_parser, double z_tolerance, bool periodic_extension, bool turn_normals = false, bool clip_z = false);

        /// Read the Eclipse grid format ('grdecl').
        /// \param input_data the data contained in a parser object.
        /// \param z_tolerance points along a pillar that are closer together in z
        ///        coordinate than this parameter, will be replaced by a single point.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        void processEclipseFormat(Opm::DeckConstPtr deck, double z_tolerance, bool periodic_extension, bool turn_normals = false, bool clip_z = false);

        /// Read the Eclipse grid format ('grdecl').
        /// \param input_data the data in grdecl format, declared in preprocess.h.
        /// \param z_tolerance points along a pillar that are closer together in z
        ///        coordinate than this parameter, will be replaced by a single point.
        /// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
        void processEclipseFormat(const grdecl& input_data, double z_tolerance, bool remove_ij_boundary, bool turn_normals = false);

	/// Create a cartesian grid.
	/// \param dims the number of cells in each cartesian direction.
	/// \param cellsize the size of each cell in each dimension.
	void createCartesian(const array<int, 3>& dims,
			     const array<double, 3>& cellsize);

	/// The logical cartesian size of the global grid.
	/// This function is not part of the Dune grid interface,
	/// and should be used with caution.
        const std::array<int, 3>& logicalCartesianSize() const
        {
            return current_view_data_->logical_cartesian_size_;
        }

	/// Retrieve mapping from internal ("compressed") active grid
	/// cells to external ("uncompressed") cells.  Specifically,
	/// @code globalCell()[i] @endcode is the linearized Cartesian
	/// index of grid cell @code i @endcode.  This method should
	/// only be used by classes which really need it, such as
	/// those dealing with permeability fields from the input deck
	/// from whence the current CpGrid was constructed.
        const std::vector<int>& globalCell() const
        {
            return current_view_data_->global_cell_;
        }

        /// @brief
        ///    Extract Cartesian index triplet (i,j,k) of an active cell.
        ///
        /// @param [in] c
        ///    Active cell index.
        ///
        /// @param [out] ijk  Cartesian index triplet
        void getIJK(const int c, std::array<int,3>& ijk) const
        {
            if(current_view_data_==distributed_data_ && comm().size() > 1)
                OPM_THROW(std::runtime_error, "getIJK is only usable"
                          << " on a non-distributed grid!");

            // TODO: replace "data_" by "current_view_data_" once this
            // method works for parallel simulations
            data_->getIJK(c, ijk);
        }

	/// Is the grid currently using unique boundary ids?
	/// \return true if each boundary intersection has a unique id
	///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
	bool uniqueBoundaryIds() const
	{
	    return current_view_data_->uniqueBoundaryIds();
	}

	/// Set whether we want to have unique boundary ids.
	/// \param uids if true, each boundary intersection will have a unique boundary id.
	void setUniqueBoundaryIds(bool uids)
	{
            current_view_data_->setUniqueBoundaryIds(uids);
	}

	// --- Dune interface below ---


        /// Return grid name. It's the same as the class name.
	/// What did you expect, something funny?
        std::string name() const
        {
            return "CpGrid";
        }


        /// Return maximum level defined in this grid. Levels are numbered
        /// 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const
	{
            return 0;
        }


        /// Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
        }


        /// one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,All_Partition>(*current_view_data_, size(codim), true );

        }


        /// Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*current_view_data_, 0, true );
        }


        /// one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*current_view_data_, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim, PiType>(*current_view_data_, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim, PiType>(*current_view_data_, size(codim), true);
        }


        /// \brief Number of grid entities per level and codim
        int size (int level, int codim) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(codim);
        }


        /// number of leaf entities per codim in this process
        int size (int codim) const
	{
            return current_view_data_->size(codim);
        }


        /// number of entities per level and geometry type in this process
        int size (int level, GeometryType type) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(type);
        }


        /// number of leaf entities per geometry type in this process
        int size (GeometryType type) const
        {
            return current_view_data_->size(type);
        }


        /// \brief Access to the GlobalIdSet
        const Traits::GlobalIdSet& globalIdSet() const
	{
            return *current_view_data_->global_id_set_;
        }


        /// \brief Access to the LocalIdSet
        const Traits::LocalIdSet& localIdSet() const
	{
            return *current_view_data_->local_id_set_;
        }


        /// \brief Access to the LevelIndexSets
        const Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return *current_view_data_->index_set_;
        }


        /// \brief Access to the LeafIndexSet
        const Traits::LeafIndexSet& leafIndexSet() const
        {
            return *current_view_data_->index_set_;
        }


        /// global refinement
        void globalRefine (int)
        {
            std::cout << "Warning: Global refinement not implemented, yet." << std::endl;
        }

        const std::vector< Dune :: GeometryType >& geomTypes( const int codim ) const 
        {
          return leafIndexSet().geomTypes( codim );
        }

	/*  No refinement implemented. GridDefaultImplementation's methods will be used.

        /// \brief Mark entity for refinement
	///
	/// This only works for entities of codim 0.
	/// The parameter is currently ignored
	///
	/// \return <ul>
	/// <li> true, if marking was succesfull </li>
	/// <li> false, if marking was not possible </li>
	/// </ul>

        bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
        {
            return hostgrid_->mark(refCount, getHostEntity<0>(*e));
        }

        /// \brief Return refinement mark for entity
	///
	/// \return refinement mark (1,0,-1)

        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            return hostgrid_->getMark(getHostEntity<0>(*e));
        }

        /// \todo Please doc me !
        bool preAdapt() {
            return hostgrid_->preAdapt();
        }


        /// Triggers the grid refinement process
        bool adapt()
        {
            return hostgrid_->adapt();
        }

        /// \brief Clean up refinement markers
        void postAdapt() {
            return hostgrid_->postAdapt();
        }

	end of refinement section */

        /// \brief Size of the overlap on the leaf level
        unsigned int overlapSize(int) const {
            return 1;
        }


        /// \brief Size of the ghost cell layer on the leaf level
        unsigned int ghostSize(int) const {
            return 0;
        }


        /// \brief Size of the overlap on a given level
        unsigned int overlapSize(int, int) const {
            return 1;
        }


        /// \brief Size of the ghost cell layer on a given level
        unsigned int ghostSize(int, int) const {
            return 0;
        }
        
        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \warning May only be called once.
        bool loadBalance()
        {
            return scatterGrid();
        }
        
        
        template<class DataHandle>
        bool loadBalance(DataHandle& data)
        {
            bool ret = scatterGrid();
            scatterData(data);
            return ret;
        }
        
        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level
        /// \param data The data handle describing the data. Has to adhere to the 
        /// Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int /*level*/) const
        {
            communicate(data, iftype, dir);
        }

        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level.
        /// \tparam DataHandle The type of the data handle describing the data.
        /// \param data The data handle describing the data. Has to adhere to the Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {
            current_view_data_->communicate(data, iftype, dir);
        }
        
        /// dummy collective communication
        const CollectiveCommunication& comm () const
        {
            return current_view_data_->ccobj_;
        }

        // ------------ End of Dune interface, start of simplified interface --------------

        // enum { dimension = 3 }; // already defined

        typedef Dune::FieldVector<double, 3> Vector;

        // Topology
        int numCells() const
        {
            return current_view_data_->cell_to_face_.size();
        }
        int numFaces() const
        {
            return current_view_data_->face_to_cell_.size();
        }
        int numVertices() const
        {
            return current_view_data_->geomVector<3>().size();
        }

        int numCellFaces(int cell) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
        }
                
        int cellFace(int cell, int local_index) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
        }

        const cpgrid::OrientedEntityTable<0,1>::row_type cellFaceRow(int cell) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
        }

        int faceCell(int face, int local_index) const
        {
            // In the parallel case we store non-existent cells for faces along
            // the front region. Theses marked with index std::numeric_limits<int>::max(),
            // orientation might be arbitrary, though.
            cpgrid::OrientedEntityTable<1,0>::row_type r
                = current_view_data_->face_to_cell_[cpgrid::EntityRep<1>(face, true)];
            bool a = (local_index == 0);
            bool b = r[0].orientation();
            bool use_first = a ? b : !b;
            // The number of valid cells.
            int size = r.size();
            // In the case of only one valid cell, this is the index of it.
            int index = 0;
            if(r[0].index()==std::numeric_limits<int>::max()){
                assert(size==2);
                --size;
                index=1;
            }
            if(r[1].index()==std::numeric_limits<int>::max())
            {
                assert(size==2);
                --size;
            }
            if (size == 2) {
                return use_first ? r[0].index() : r[1].index();
            } else {
                return use_first ? r[index].index() : -1;
            }
        }
        int numCellFaces() const
        {
            return current_view_data_->cell_to_face_.dataSize();
        }
        int numFaceVertices(int face) const
        {
            return current_view_data_->face_to_point_[face].size();
        }
        int faceVertex(int face, int local_index) const
        {
            return current_view_data_->face_to_point_[face][local_index];
        }

        // Geometry
        const Vector& vertexPosition(int vertex) const
        {
            return current_view_data_->geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
        }
        double faceArea(int face) const
        {
            return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
        }
        const Vector& faceCentroid(int face) const
        {
            return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
        }
        const Vector& faceNormal(int face) const
        {
            return current_view_data_->face_normals_.get(face);
        }
        double cellVolume(int cell) const
        {
            return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
        }
        const Vector& cellCentroid(int cell) const
        {
            return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
        }

        /// \brief An iterator over the centroids of the geometry of the entities.
        /// \tparam codim The co-dimension of the entities.
        template<int codim>
        class CentroidIterator
            : public RandomAccessIteratorFacade<CentroidIterator<codim>,
                                                FieldVector<double, 3>,
                                                const FieldVector<double, 3>&, int>
        {
        public:
            /// \brief The type of the iterator over the codim geometries.
            typedef typename std::vector<cpgrid::Geometry<3-codim, 3> >::const_iterator
            GeometryIterator;
            /// \brief Constructs a new iterator from an iterator over the geometries.
            /// \param iter The iterator of the geometry objects.
            CentroidIterator(GeometryIterator iter)
            : iter_(iter)
            {}

            const FieldVector<double, 3>& dereference() const
            {
                return iter_->center();
            }
            void increment()
            {
                ++iter_;
            }
            const FieldVector<double, 3>& elementAt(int n)
            {
                return iter_[n]->center();
            }
            void advance(int n)
            {
                iter_+=n;
            }
            void decrement()
            {
                --iter_;
            }
            int distanceTo(const CentroidIterator& o)
            {
                return o-iter_;
            }
            bool equals(const CentroidIterator& o) const
            {
                return o==iter_;
            }
        private:
            /// \brief The iterator over the underlying geometry objects.
            GeometryIterator iter_;
        };

        CentroidIterator<0> beginCellCentroids() const
        {
            return CentroidIterator<0>(current_view_data_->geomVector<0>().begin());
        }

        CentroidIterator<1> beginFaceCentroids() const
        {
            return CentroidIterator<1>(current_view_data_->geomVector<1>().begin());
        }

        // Extra
        int boundaryId(int face) const
        {
            int ret = 0;
            cpgrid::EntityRep<1> f(face, true);
            if (current_view_data_->face_to_cell_[f].size() == 1) {
                if (current_view_data_->uniqueBoundaryIds()) {
                    // Use the unique boundary ids.
                    ret = current_view_data_->unique_boundary_ids_[f];
                } else {
                    // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
                    const bool normal_is_in = 
                        !(current_view_data_->face_to_cell_[f][0].orientation());
                    enum face_tag tag = current_view_data_->face_tag_[f];
                    switch (tag) {
                    case LEFT:
                        //                   LEFT : RIGHT
                        ret = normal_is_in ? 1    : 2; // min(I) : max(I)
                        break;
                    case BACK:
                        //                   BACK : FRONT
                        ret = normal_is_in ? 3    : 4; // min(J) : max(J)
                        break;
                    case TOP:
                        // Note: TOP at min(K) as 'z' measures *depth*.
                        //                   TOP  : BOTTOM
                        ret = normal_is_in ? 5    : 6; // min(K) : max(K)
                        break;
                    }
                }
            }
            return ret;
        }


        // ------------ End of simplified interface --------------

        //------------- methods not in the DUNE grid interface.
        
        ///
        /// \brief Moves data from the global (all data on process) view to the distributed view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void scatterData(DataHandle& handle)
        {
#if HAVE_MPI && DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
            if(!distributed_data_)
                OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balanced grid!");
            distributed_data_->scatterData(handle, data_, distributed_data_);
#else
            // Suppress warnings for unused argument.
            (void) handle;
#endif
        }

        ///
        /// \brief Moves data from the distributed view to the global (all data on process) view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void gatherData(DataHandle& handle)
        {
#if HAVE_MPI && DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
            if(!distributed_data_)
                OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balance grid!");
            distributed_data_->gatherData(handle, data_, distributed_data_);
#else
            // Suppress warnings for unused argument.
            (void) handle;
#endif
        }

        /// \brief Switch to the global view.
        void switchToGlobalView()
        {
            current_view_data_=data_;
        }
        
        /// \brief Switch to the distributed view.
        void switchToDistributedView()
        {
            current_view_data_=distributed_data_;
        }
        
    private:
        /// Scatter a global grid to all processors.
        bool scatterGrid();
        
        /** @brief The data stored in the grid. 
         * 
         * All the data of the grid is stored there and
         * calls are forwarded to it.*/
        cpgrid::CpGridData *data_;
        /** @brief A pointer to data of the current View. */
        cpgrid::CpGridData *current_view_data_;
        /** @brief The data stored for the distributed grid. */
        cpgrid::CpGridData * distributed_data_;
    }; // end Class CpGrid



    namespace Capabilities
    {
	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 0>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 3>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
	template <>
	struct isParallel<CpGrid>
	{
	    static const bool v = true;
	};

        template<>
        struct canCommunicate<CpGrid,0>
        {
            static const bool v = true;
        };

        template<>
        struct canCommunicate<CpGrid,3>
        {
            static const bool v = true;
        };

	/// \todo Please doc me !
	template <>
	struct hasBackupRestoreFacilities<CpGrid>
	{
	    static const bool v = false;
	};

    }

} // namespace Dune

#include <dune/grid/cpgrid/PersistentContainer.hpp>
#endif // OPM_CPGRID_HEADER
