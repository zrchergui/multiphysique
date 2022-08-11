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

/*!
 * \file
 *
 * \brief Specification of the material parameters
 *       for the Brooks Corey constitutive relations.
 */
#ifndef DUMUX_BROOKS_COREY_PARAMS_HH
#define DUMUX_BROOKS_COREY_PARAMS_HH

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \brief Specification of the material parameters
 *       for the Brooks Corey constitutive relations.
 *
 *        \ingroup fluidmatrixinteractionsparams
 *
 *\see BrooksCorey
 */
template <class ScalarT>
class BrooksCoreyParams
{
public:
    typedef ScalarT Scalar;

    BrooksCoreyParams()
    {
        Valgrind::SetUndefined(*this);
    }

    BrooksCoreyParams(Scalar pe, Scalar lambdaKr, Scalar lambdaPc)
        : pe_(pe), lambdaKr_(lambdaKr), lambdaPc_(lambdaPc)
    {
    }

    /*!
     * \brief Returns the entry pressure in \f$\mathrm{[Pa]}\f$
     */
    Scalar pe() const
    { return pe_; }

    /*!
     * \brief Set the entry pressure in \f$\mathrm{[Pa]}\f$]
     */
    void setPe(Scalar v)
    { pe_ = v; }


    /*!
     * \brief Returns the lambda shape parameter for relative permeability \f$\mathrm{[-]}\f$
     */
    Scalar lambdaKr() const
    { return lambdaKr_; }

    /*!
     * \brief Set the lambda shape parameter for relative permeability \f$\mathrm{[-]}\f$
     */
    void setLambdaKr(Scalar v)
    { lambdaKr_ = v; }

        /*!
     * \brief Returns the lambda shape parameter for capillary pressure \f$\mathrm{[-]}\f$
     */
    Scalar lambdaPc() const
    { return lambdaPc_; }

    /*!
     * \brief Set the lambda shape parameter for capillary pressure \f$\mathrm{[-]}\f$
     */
    void setLambdaPc(Scalar v)
    { lambdaPc_ = v; }

private:
    Scalar pe_;
    Scalar lambdaKr_;
    Scalar lambdaPc_;
};
} // namespace Dumux

#endif
