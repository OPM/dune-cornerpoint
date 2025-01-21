//===========================================================================
//
// File: LookUpData.hh
//
// Created: Tue May 23 14:44:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.

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
#ifndef OPM_LOOKUPDATA_HH
#define OPM_LOOKUPDATA_HH

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/grid/cpgrid/Entity.hpp>

#include <functional>
#include <string>
#include <type_traits>
#include <vector>

namespace Dune
{
class CpGrid;
}

namespace Opm
{

/// LookUpData class - To search field properties of leaf grid view elements via element/elementIndex
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template <typename Grid, typename GridView>
class LookUpData
{
public:
    /// \brief:     Constructor taking a GridView and a bool
    /// \param [in] GridView
    /// \param [in] isFieldPropInLgr   bool: default false (search field property in unrefined grid)
    ///                                      true (search field property in refined grid; LGR-id/level required)
    ///                                      Currently, isFieldPropInLgr_ == false means that all the field
    ///                                      properties are given in the unrefined grid (level 0).
    explicit LookUpData(const  GridView& gridView, bool isFieldPropInLgr = false) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout()),
        isFieldPropInLgr_(isFieldPropInLgr)
    {
    }

    /// \brief: Get field propertry for an element in the leaf grid view, from a vector and element index.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    template<typename FieldPropType>
    FieldPropType operator()(const int& elemIdx, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Get field propertry for an element in the leaf grid view, from a vector.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    template<typename EntityType, typename FieldPropType>
    typename std::enable_if_t<!std::is_same_v<EntityType, grid_size_t>, FieldPropType>
    operator()(const EntityType& elem, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Get field property of type double from field properties manager by name.
    std::vector<double> assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                     const std::string& propString) const;

    /// \brief: Get field property of type int from field properties manager by name.
    template<typename IntType>
    std::vector<IntType> assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                   const std::string& propString,
                                                   const bool& needsTranslation,
                                                   std::function<void(IntType, int)> valueCheck = [](IntType, int){}) const;

    /// \brief: Get property of type double from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const ElemOrIndex& elemOrIndex) const;

    /// \brief: Get property of type int from field properties manager by name, via element.
    template<typename ElemOrIndex>
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const ElemOrIndex& elemOrIndex) const;

    /// \brief: Return the same element index for all grids different from CpGrid.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
    getFieldPropIdx(const EntityType& elem) const;

    /// \brief: Return index to search for the field propertries, for CpGrids.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
    getFieldPropIdx(const EntityType& elem) const;

    /// \brief: Return the same element index for all grids different from CpGrid.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename GridType>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdx(const int& elemIdx) const;

    /// \brief: Return the index to search for the field properties, for CpGrids.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropIdx(const int& elemIdx) const;


protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    bool isFieldPropInLgr_;
}; // end LookUpData class

/// LookUpCartesianData - To search field properties of leaf grid view elements via CartesianIndex (cartesianMapper)
///
/// Instead of using a specialitation for Dune::CpGrid, we implement std::enable_if
/// to overload methods with different definitions: for Dune:CpGrid and for other
/// Grid types. An auxiliary defualt template parameter (GridType = Grid) is added
/// to deal with the dependent names at template instantiation.
template<typename Grid, typename GridView>
class LookUpCartesianData
{
public:
    /// \brief: Constructor taking a GridView, a CartesianIndexMapper, and a bool
    ///
    /// \param [in] gridView
    /// \param [in] mapper   Dune::CartesianIndexMapper<Grid>.
    /// \param [in] isFieldPropInLgr bool: default false (search field property in unrefined grid)
    ///                                    true (search field property in refined grid; LGR-id/level required)
    ///                                    Currently, isFieldPropInLgr_ == false means that all the field
    ///                                    properties are given in the unrefined grid (level 0).
    explicit LookUpCartesianData(const GridView& gridView,
                                 const Dune::CartesianIndexMapper<Grid>& mapper,
                                 bool isFieldPropInLgr = false) :
        gridView_(gridView),
        elemMapper_(gridView, Dune::mcmgElementLayout()),
        cartMapper_(&mapper),
        isFieldPropInLgr_(isFieldPropInLgr)
    {
    }

    /// \brief: Get field property for an element in the leaf grid view, from a vector, via Cartesian Index.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    template<typename FieldPropType>
    FieldPropType operator()(const int& elemIdx, const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Get field property for an element in the leaf grid view, from a vector, via Cartesian Index.
    ///
    ///         For general grids, the field property vector is assumed to be given for the gridView_.
    ///         For CpGrid, the field property vector is assumed to be given for level 0 when isFieldPropLgr_ == false,
    ///         and for certain LGR/level > 0 when isFieldPropLgr_ == true.
    template<typename EntityType, typename FieldPropType>
    typename std::enable_if_t<!std::is_same_v<EntityType, grid_size_t>, FieldPropType>
    operator()(const EntityType& elem,const std::vector<FieldPropType>& fieldProp) const;

    /// \brief: Get field property of type double from field properties manager by name.
    std::vector<double> assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                     const std::string& propString) const;

    /// \brief: Get field property of type int from field properties manager by name.
    template<typename IntType>
    std::vector<IntType> assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                   const std::string& propString,
                                                   const bool& needsTranslation,
                                                   std::function<void(IntType, int)> valueCheck = [](IntType, int){}) const;

    /// \brief: Get property of type double from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    double fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                           const std::string& propString,
                           const ElemOrIndex& elemOrIndex) const;

    /// \brief: Get property of type int from field properties manager by name, via element or its index.
    template<typename ElemOrIndex>
    int fieldPropInt(const FieldPropsManager& fieldPropsManager,
                     const std::string& propString,
                     const ElemOrIndex& elemOrIndex) const;

    /// \brief: Return the same element index for all grids different from CpGrid.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>  && !std::is_same_v<EntityType, grid_size_t>,int>
    getFieldPropCartesianIdx(const EntityType& elem) const;

    /// \brief: Return index to search for the field propertries, for CpGrids.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename EntityType, typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
    getFieldPropCartesianIdx(const EntityType& elem) const;

    /// \brief: Return the same element index for all grids different from CpGrid.
    ///
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename GridType = Grid>
    typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdx(const int& elemIdx) const;

    /// \brief: Return index to search for the field propertries, for CpGrids.
    ///
    ///         When isFieldPropInLgr_ == false : fieldPropIdx == Index of the origin cell (parent/equivalent cell when element
    ///                                           has nofather) in level 0.
    ///         When isFieldPropInLgr_ == true  : fieldPropIdx == Index of the equivalent cell in the LGR (level>0).
    /// \tparam     GridType    Auxiliary type to overload the method, distinguishing general grids from CpGrid, with std::enable_if.
    ///                         Default: GridType = Grid.
    template<typename GridType = Grid>
    typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
    getFieldPropCartesianIdx(const int& elemIdx) const;

protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
    bool isFieldPropInLgr_;
}; // end LookUpCartesianData class
}
// end namespace Opm



/// LookUpData

template<typename Grid, typename GridView>
template<typename FieldPropType>
FieldPropType Opm::LookUpData<Grid,GridView>::operator()(const int& elemIdx,
                                                         const std::vector<FieldPropType>& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx);
    assert(0 <= fieldPropIdx && static_cast<int>(fieldProp.size()) > fieldPropIdx);
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
template<typename EntityType, typename FieldPropType>
typename std::enable_if_t<!std::is_same_v<EntityType, grid_size_t>,FieldPropType>
Opm::LookUpData<Grid,GridView>::operator()(const EntityType& elem,
                                           const std::vector<FieldPropType>& fieldProp) const
{
    const auto& fieldPropIdx = this->getFieldPropIdx<EntityType,Grid>(elem);
    assert( (0 <= fieldPropIdx) && (static_cast<int>(fieldProp.size()) > fieldPropIdx));
    return fieldProp[fieldPropIdx];
}

template<typename Grid, typename GridView>
std::vector<double> Opm::LookUpData<Grid,GridView>::assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                 const std::string& propString) const
{
    std::vector<double> fieldPropOnLeaf;
    grid_size_t numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_double(propString);
    if ( (propString == "PORV") && (gridView_.grid().maxLevel() > 0)) {
        // PORV poreVolume. LGRs supported (so far) only for CpGrid.
        // For CpGrid with LGRs, poreVolume of a cell on the leaf grid view which has a parent cell on level 0,
        // is computed as  porv[parent] * leafCellVolume / parentCellVolume. In this way, the sum of the pore
        // volume of a parent cell coincides with the sum of the pore volume of its children.
        for (const auto& element : elements(gridView_)) {
            const auto& elemIdx = this-> elemMapper_.index(element);
            const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
            if (element.hasFather()) {
                const auto fatherVolume = element.father().geometry().volume();
                const auto& elemVolume = element.geometry().volume();
                fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx] * elemVolume / fatherVolume;
            }
            else {
                fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx];
            }
        }
    }
    else {
        for (const auto& element : elements(gridView_)) {
            const auto& elemIdx = this-> elemMapper_.index(element);
            const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
            fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx];
        }
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename IntType>
std::vector<IntType> Opm::LookUpData<Grid,GridView>::assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                               const std::string& propString,
                                                                               const bool& needsTranslation,
                                                                               std::function<void(IntType, int)> valueCheck) const
{
    std::vector<IntType> fieldPropOnLeaf;
    grid_size_t numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_int(propString);
    for (const auto& element : elements(gridView_)) {
        const auto& elemIdx = this-> elemMapper_.index(element);
        const auto& fieldPropIdx = this->getFieldPropIdx<Grid>(elemIdx); // gets parentIdx (or (lgr)levelIdx) for CpGrid with LGRs
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropIdx] - needsTranslation;
        valueCheck(fieldProp[fieldPropIdx], fieldPropIdx);
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
double Opm::LookUpData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                       const std::string& propString,
                                                       const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
int Opm::LookUpData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                 const std::string& propString,
                                                 const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    assert(elem.level() == 0); // LGRs (level>0) only supported for CpGrid.
    return this-> elemMapper_.index(elem);
}

template<typename Grid, typename GridView>
template<typename EntityType,typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    static_assert(std::is_same_v<EntityType,Dune::cpgrid::Entity<0>>);
    if (isFieldPropInLgr_ && elem.level()) { // level > 0 == true ; level == 0 == false
        // In case some LGRs do not have refined field properties, the next line need to be modified.
        return elem.getLevelElem().index();
    }
    else {
        return elem.getOrigin().index();
    }
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
    assert(gridView_.grid().maxLevel() == 0);
    return elemIdx;
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpData<Grid,GridView>::getFieldPropIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().currentData().back()), elemIdx, true);
    if (isFieldPropInLgr_ && elem.level()) { // level > 0 == true ; level == 0 == false
        // In case some LGRs do not have refined field properties, the next line need to be modified.
        return elem.getLevelElem().index();
    }
    else {
        return elem.getOrigin().index();
    }
}



/// LookUpCartesianData

template<typename Grid, typename GridView>
template<typename FieldPropType>
FieldPropType Opm::LookUpCartesianData<Grid,GridView>::operator()(const int& elemIdx,
                                                                  const std::vector<FieldPropType>& fieldProp) const
{
    assert(cartMapper_);
    const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
    assert(0 <=  fieldPropCartIdx && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx));
    return fieldProp[fieldPropCartIdx];
}

template<typename Grid, typename GridView>
template<typename EntityType, typename FieldPropType>
typename std::enable_if_t<!std::is_same_v<EntityType, grid_size_t>,FieldPropType>
Opm::LookUpCartesianData<Grid,GridView>::operator()(const EntityType& elem,
                                                    const std::vector<FieldPropType>& fieldProp) const
{
    assert(cartMapper_);
    const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<EntityType,Grid>(elem);
    assert( (0 <= fieldPropCartIdx) && (static_cast<int>(fieldProp.size()) > fieldPropCartIdx) );
    return fieldProp[fieldPropCartIdx];
}

template<typename Grid, typename GridView>
std::vector<double> Opm::LookUpCartesianData<Grid,GridView>::assignFieldPropsDoubleOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                          const std::string& propString) const
{
    std::vector<double> fieldPropOnLeaf;
    grid_size_t numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_double(propString);
    for (grid_size_t elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropCartIdx];
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename IntType>
std::vector<IntType> Opm::LookUpCartesianData<Grid,GridView>::assignFieldPropsIntOnLeaf(const FieldPropsManager& fieldPropsManager,
                                                                                        const std::string& propString,
                                                                                        const bool& needsTranslation,
                                                                                        std::function<void(IntType, int)> valueCheck) const
{
    std::vector<IntType> fieldPropOnLeaf;
    grid_size_t numElements = gridView_.size(0);
    fieldPropOnLeaf.resize(numElements);
    const auto& fieldProp = fieldPropsManager.get_int(propString);
    for (grid_size_t elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        const auto fieldPropCartIdx = this->getFieldPropCartesianIdx<Grid>(elemIdx);
        fieldPropOnLeaf[elemIdx] = fieldProp[fieldPropCartIdx] - needsTranslation;
        valueCheck(fieldProp[fieldPropCartIdx], fieldPropCartIdx);
    }
    return fieldPropOnLeaf;
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
double Opm::LookUpCartesianData<Grid,GridView>::fieldPropDouble(const FieldPropsManager& fieldPropsManager,
                                                                const std::string& propString,
                                                                const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_double(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename ElemOrIndex>
int Opm::LookUpCartesianData<Grid,GridView>::fieldPropInt(const FieldPropsManager& fieldPropsManager,
                                                          const std::string& propString,
                                                          const ElemOrIndex& elemOrIndex) const
{
    const auto& fieldPropVec = fieldPropsManager.get_int(propString);
    return this ->operator()(elemOrIndex,fieldPropVec);
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    // Check there are no LGRs. LGRs (level>0) only supported for CpGrid.
    assert(gridView_.grid().maxLevel() == 0);
    return cartMapper_-> cartesianIndex(this->elemMapper_.index(elem));
}

template<typename Grid, typename GridView>
template<typename EntityType, typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid> && !std::is_same_v<EntityType, grid_size_t>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const EntityType& elem) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    if (isFieldPropInLgr_ && elem.level()) { // level == 0 false; level > 0 true
        return elem.getLevelCartesianIdx();
    }
    else {
        return cartMapper_-> cartesianIndex(this->elemMapper_.index(elem));
    }
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<!std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    return  cartMapper_-> cartesianIndex(elemIdx);
}

template<typename Grid, typename GridView>
template<typename GridType>
typename std::enable_if_t<std::is_same_v<GridType,Dune::CpGrid>,int>
Opm::LookUpCartesianData<Grid,GridView>::getFieldPropCartesianIdx(const int& elemIdx) const
{
    static_assert(std::is_same_v<Grid,GridType>);
    const auto& elem = Dune::cpgrid::Entity<0>(*(gridView_.grid().currentData().back()), elemIdx, true);
    return this -> getFieldPropCartesianIdx<Dune::cpgrid::Entity<0>,Dune::CpGrid>(elem);
}

#endif
