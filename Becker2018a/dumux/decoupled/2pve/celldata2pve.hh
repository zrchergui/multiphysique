// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_CELLDATA2PVE_HH
#define DUMUX_CELLDATA2PVE_HH

#include <dumux/porousmediumflow/2p/sequential/celldata.hh>

/**
 * \file
 * \brief  Class including data of one grid cell
 */

namespace Dumux
{

/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations for one grid cell.
/*! The variables of two-phase flow, which are phase pressures and saturations are stored in this class. Further, resulting cell values for constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure are stored.
 * Additionally, data assigned to cell-cell interfaces, so-called flux-data are stored.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization: in case of incompressible flow bool = <tt>false</tt>
 */
template<class TypeTag>
class CellData2PVE: public CellData2P<TypeTag, false>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

private:
    Scalar minGasPlumeDist_;
    Scalar gasPlumeDist_;
    Scalar residualSegSaturation_;
    Scalar oldSaturation_;
public:

    //! Constructs a CellData2P object
    CellData2PVE()
    {
        minGasPlumeDist_ = 1e100;
        gasPlumeDist_ = 0.0;
        residualSegSaturation_ = 0.0;
        oldSaturation_ = 0.0;
    }

    /*!\brief Returns gasPlumeDist_, the height of plume distance from bottom
     */
    Scalar gasPlumeDist() const
    {
        return gasPlumeDist_;
    }

    /*!\brief Sets gasPlumeDist_, the height of plume distance from bottom
     *
     * \param Xi gasPlumeDist
     */
    void setGasPlumeDist(Scalar Xi)
    {
        gasPlumeDist_ = Xi;
    }

    /*!\brief Returns minGasPlumeDist_, the minimum height of plume distance from bottom
     */
    Scalar minGasPlumeDist() const
    {
        return minGasPlumeDist_;
    }

    /*!\brief Sets minGasPlumeDist_, the minimum height of plume distance from bottom
     *
     * \param Xi gasPlumeDist
     */
    void setMinGasPlumeDist(Scalar Xi)
    {
        minGasPlumeDist_ = Xi;
    }

    /*!\brief Returns residualSegSaturation_, the pseudo wetting phase saturation inside plume
     */
   Scalar residualSegSaturation() const
   {
       return residualSegSaturation_;
   }

   /*!\brief Sets residualSegSaturation_, the pseudo wetting phase saturation inside plume
    *
    * \param Xi gasPlumeDist
    */
   void setResidualSegSaturation(Scalar resSegSaturation)
   {
       residualSegSaturation_ = resSegSaturation;
   }

  /*!\brief Sets oldResidualSegSaturation_, saturation in plume from last time step
   *
   * \param oldSaturation wetting phase saturation
   */
  void setOldSaturation(Scalar oldSaturation)
  {
      oldSaturation_ = oldSaturation;
  }

  /*!\brief Returns oldResidualSegSaturation_, saturation in plume from last time step
   */
 Scalar oldSaturation() const
 {
     return oldSaturation_;
 }

};
}
#endif
