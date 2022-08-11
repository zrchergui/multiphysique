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
 * \brief   Parameters for the linear capillary pressure and
 *          relative permeability <-> saturation relations
 */
#ifndef EXPONENTIAL_MATERIAL_PARAMS_HH
#define EXPONENTIAL_MATERIAL_PARAMS_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of params for the linear material
 *        law.
 *
 *        \ingroup fluidmatrixinteractionsparams
 */
template<class ScalarT>
class ExponentialMaterialParams
{
public:
    typedef ScalarT Scalar;

    ExponentialMaterialParams()
    {}

    ExponentialMaterialParams(Scalar entryPc, Scalar lambda, Scalar exponent)
    {
        setEntryPc(entryPc);
        setMaxPc(lambda);
        setExponent(exponent);
    }


    /*!
     * \brief Return the entry pressure for the linear material law in \f$\mathrm{[Pa]}\f$.
     *
     * The entry pressure is reached at \f$\mathrm{\overline{S}_w = 1}\f$
     */
    Scalar pe() const
    { return entryPc_; }

    /*!
     * \brief Set the entry pressure for the linear material law in \f$\mathrm{[Pa]}\f$.
     *
     * The entry pressure is reached at \f$\mathrm{\overline{S}_w = 1}\f$
     */
    void setPe(Scalar v)
    { entryPc_ = v; }

    /*!
     * \brief Return the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    Scalar lambda() const
    { return lambda_; }

    /*!
     * \brief Set the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    void setLambda(Scalar v)
    { lambda_ = v; }

    /*!
     * \brief Return the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    Scalar exponent() const
    { return exponent_; }

    /*!
     * \brief Set the maximum capillary pressure for the linear material law in \f$\mathrm{[Pa]}\f$..
     *
     * The maximum capillary pressure is reached at \f$\mathrm{\overline{S}_w = 0}\f$
     */
    void setExponent(Scalar exponent)
    { exponent_ = exponent; }

    /*!
     * \brief Dummy
     */
//    Scalar lambda() const
//    { return 0.0; }

private:
    Scalar entryPc_;
    Scalar lambda_;
    Scalar exponent_;
};
} // namespace Dumux

#endif
