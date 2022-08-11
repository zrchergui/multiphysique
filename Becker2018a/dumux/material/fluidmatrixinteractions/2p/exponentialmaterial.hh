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
 * \brief   Linear capillary pressure and
 *          relative permeability <-> saturation relations
 */
#ifndef EXPONENTIAL_MATERIAL_HH
#define EXPONENTIAL_MATERIAL_HH

#include <dune/common/deprecated.hh>
#include "exponentialmaterialparams.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Linear capillary pressure and
 * relative permeability <-> saturation relations
 *
 * The entry pressure is reached at \f$\mathrm{ \overline{S}_w = 1}\f$, the maximum
 * capillary pressure is observed at \f$\mathrm{ \overline{S}_w = 0}\f$.
 *
 * For general info: EffToAbsLaw
 *
 * \see LinearMaterialParams
 */
template <class ScalarT, class ParamsT = ExponentialMaterialParams<ScalarT> >
class ExponentialMaterial
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f$\mathrm{
     p_C = (1 - \overline{S}_w ) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     }\f$
     *
     * \param swe Effective saturation of the wetting phase \f$\overline{S}_w\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Capillary pressure calculated by linear constitutive relation.
     */
    static Scalar pc(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return params.pe()*pow(swe, -1.0/params.lambda());
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f$\mathrm{
         S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     }\f$
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$ in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Effective wetting phase saturation calculated as inverse of the linear constitutive relation.
     */
    static Scalar sw(const Params &params, Scalar pc)
    {
        assert(pc >= 0);

        Scalar tmp = pow(pc/params.pe(), -params.lambda());
        return std::min(std::max(tmp, Scalar(0.0)), Scalar(1.0));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * This is equivalent to
     * \f$\mathrm{
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     }\f$
     * \param swe  Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$\mathrm{[p_c]}\f$ w.r.t. effective saturation according to linear material relation.
    */
    static Scalar dpc_dswe(const Params &params, Scalar swe)
    {
        assert(0 <= swe && swe <= 1);

        return - params.pe()/params.lambda() * pow(swe, -1/params.lambda() - 1);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     *
     * \param pc Capillary pressure \f$\mathrm{[p_C]}\f$  in \f$\mathrm{[Pa]}\f$.
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return Partial derivative of effective saturation w.r.t. \f$\mathrm{[p_c]}\f$ according to linear relation.
     */
    static Scalar dswe_dpc(const Params &params, Scalar pc)
    {
        assert(pc >= 0);

        return -params.lambda()/params.pe() * pow(pc/params.pe(), - params.lambda() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \return Relative permeability of the wetting phase calculated as linear relation.
     */
    static Scalar krw(const Params &params, Scalar swe)
    {
        if (swe <= 0.0)
            return 0.0;
        else if (swe >= 1.0)
            return 1.0;

        if(params.exponent() == 1)
            return swe;
        else if (params.exponent() == 2)
            return swe*swe;
        else if (params.exponent() == 3)
            return swe*swe*swe;
        else if (params.exponent() == 4)
            return swe*swe*swe*swe;
        else if (params.exponent() == 5)
            return swe*swe*swe*swe*swe;
        else
            return std::pow(swe, params.exponent());
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param params A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \param swe Effective saturation of the wetting phase \f$\mathrm{[\overline{S}_w]}\f$ conversion from absolute saturation happened in EffToAbsLaw.
     * \return Relative permeability of the non-wetting phase calculated as linear relation.
     */
    static Scalar krn(const Params &params, Scalar swe)
    {
        if (swe >= 1.0)
            return 0.0;
        else if (swe <= 0.0)
            return 1.0;

        Scalar sne = 1 - swe;
        if(params.exponent() == 1)
            return sne;
        else if (params.exponent() == 2)
            return sne*sne;
        else if (params.exponent() == 3)
            return sne*sne*sne;
        else if (params.exponent() == 4)
            return sne*sne*sne*sne;
        else if (params.exponent() == 5)
            return sne*sne*sne*sne*sne;
        else
            return std::pow(sne, params.exponent());
    }
};
}

#endif
