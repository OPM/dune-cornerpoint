/*
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
Copyright 2009, 2010, 2022 Equinor ASA.

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
#ifdef HAVE_CONFIG_H
#include"config.h"
#endif
#include"Intersection.hpp"
#include"EntityRep.hpp"
#include"Entity.hpp"
#include "CpGridData.hpp"

namespace Dune
{
namespace cpgrid
{

Intersection::Intersection(const CpGridData& grid, const EntityRep<0>& cell, int subindex, bool update_now)

                : pgrid_(&grid),
                  index_(cell.index()),
                  subindex_(subindex),
                  faces_of_cell_(grid.cell_to_face_[cell]),
                  nbcell_(cell.index()), // Init to self, which is invalid.
                  is_on_boundary_(false)
            {
                assert(index_ >= 0);
                if (update_now) {
                    update();
                }
            }
int Intersection::boundaryId() const
            {
                if (boundary()) {
                    if (pgrid_->uniqueBoundaryIds()) {
                        // Use the unique boundary ids.
                        OrientedEntityTable<0,1>::ToType face = faces_of_cell_[subindex_];
                        return pgrid_->unique_boundary_ids_[face];
                    } else {
                        // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
                        if (pgrid_->face_tag_[faces_of_cell_[subindex_]] == NNC_FACE)
                        {
                            // This should not be possible, as NNC "faces" always
                            // have two cell neighbours and thus are not on the boundary.
                            OPM_THROW(std::logic_error, "NNC face at boundary. This should never happen!");
                        }

                        return indexInInside() + 1;
                    }
                }
                return 0;
            }

int Intersection::boundarySegmentIndex() const
            {
                // Since this is almost the same that we did for
                // 'unique boundary ids' we use those numbers, although since
                // they are 1-based and not 0-based we must be careful.
                if (!boundary()) {
                    OPM_THROW(std::runtime_error, "Cannot call boundarySegmentIndex() on non-boundaries.");
                }
                assert(!pgrid_->unique_boundary_ids_.empty());
                // Use the unique boundary ids (subtract 1).
                const EntityRep<1>& face = faces_of_cell_[subindex_];
                return pgrid_->unique_boundary_ids_[face] - 1;
            }
void Intersection::update()
            {
                const EntityRep<1>& face = faces_of_cell_[subindex_];
                OrientedEntityTable<1,0>::row_type cells_of_face = pgrid_->face_to_cell_[face];
                is_on_boundary_ = cells_of_face.size() == 1;
                // Wether there is no nother nbcell for this intersection
                // i.e. either this on the boundary or a front intersection
                bool has_no_nbcell = is_on_boundary_ ||
                    cells_of_face[0].index()==std::numeric_limits<int>::max() ||
                    cells_of_face[1].index()==std::numeric_limits<int>::max();
                if (has_no_nbcell) {
                    nbcell_ = std::numeric_limits<int>::max(); // neighbor is not within this process
                } else {
                    assert(cells_of_face.size() == 2);
                    if (cells_of_face[0].index() == index_) {
                        nbcell_ = cells_of_face[1].index();
                    } else {
                        nbcell_ = cells_of_face[0].index();
                    }
                }
            }

void Intersection::increment()
            {
                ++subindex_;
                if (subindex_ < faces_of_cell_.size()) {
                    update();
                }
            }

int Intersection::indexInInside() const
{
    // Use the face tags (I_FACE = 0, J_FACE = 1, K_FACE = 2)
    // to decide if an intersection is
    // on an x, y, or z face and use orientations to decide
    // if its (for example) an xmin or xmax face.
    const auto& face = faces_of_cell_[subindex_];
    const auto tag = pgrid_->face_tag_[face];
    if (tag == NNC_FACE) {
        return -1;
    }
    else {
        return 2 * static_cast<int>(tag) +
            static_cast<int>(face.orientation());
    }
}

FieldVector<Intersection::ctype, 3> Intersection::outerNormal(const FieldVector<ctype, 2>&) const
{
    return pgrid_->face_normals_[faces_of_cell_[subindex_]];
}

FieldVector<Intersection::ctype, 3> Intersection::integrationOuterNormal(const FieldVector<ctype, 2>& unused) const
{
    FieldVector<ctype, 3> n = pgrid_->face_normals_[faces_of_cell_[subindex_]];
    return n*=geometry().integrationElement(unused);
}

FieldVector<Intersection::ctype, 3> Intersection::unitOuterNormal(const FieldVector<ctype, 2>&) const
{
    return pgrid_->face_normals_[faces_of_cell_[subindex_]];
}

FieldVector<Intersection::ctype, 3> Intersection::centerUnitOuterNormal() const
{
    return pgrid_->face_normals_[faces_of_cell_[subindex_]];
}

Intersection::Entity Intersection::inside() const
{
    return Entity(*pgrid_, index_, true);
}

Intersection::Entity Intersection::outside() const
{
    return Entity(*pgrid_, nbcell_, true);
}

Intersection::Geometry Intersection::geometry() const
{
    return pgrid_-> geometry_.geomVector<1>()[faces_of_cell_[subindex_]];
}


} // end namespace cpgrid
} // end namespace Dune
