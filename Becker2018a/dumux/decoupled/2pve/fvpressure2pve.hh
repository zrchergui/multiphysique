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
#ifndef DUMUX_FVPRESSURE2PVE_HH
#define DUMUX_FVPRESSURE2PVE_HH

#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressure.hh>

/**
 * \file
 * \brief  Finite Volume discretization of a two-phase flow pressure equation.
 */

namespace Dumux {
//! \ingroup FVPressure2p
/*!  \brief Finite Volume discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * This model solves equations of the form
 * \f[
 * \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) + \text{div}\, \boldsymbol{v}_{total} = q.
 * \f]
 * The definition of the total velocity \f$\boldsymbol{v}_{total}\f$ depends on the choice of the primary pressure variable.
 * Further, fluids can be assumed to be compressible or incompressible (Property: <tt>EnableCompressibility</tt>).
 * In the incompressible case a wetting \f$(w) \f$ phase pressure as primary variable leads to
 *
 * \f[
 * - \text{div}\,  \left[\lambda \boldsymbol K \left(\textbf{grad}\, p_w + f_n \textbf{grad}\,
 *   p_c + \sum f_\alpha \rho_\alpha \, g \, \textbf{grad}\, z\right)\right] = q,
 * \f]
 *
 * a non-wetting (\f$ n \f$) phase pressure yields
 * \f[
 *  - \text{div}\,  \left[\lambda \boldsymbol K  \left(\textbf{grad}\, p_n - f_w \textbf{grad}\,
 *    p_c + \sum f_\alpha \rho_\alpha \, g  \, \textbf{grad}\, z\right)\right] = q,
 *  \f]
 * and a global pressure leads to
 * \f[
 * - \text{div}\, \left[\lambda \boldsymbol K \left(\textbf{grad}\,
 *   p_{global} + \sum f_\alpha \rho_\alpha \, g \, \textbf{grad}\, z\right)\right] = q.
 * \f]
 * Here, \f$ p_\alpha \f$ is a phase pressure, \f$ p_ {global} \f$ the global pressure of a classical fractional flow formulation
 * (see e.g. P. Binning and M. A. Celia, ''Practical implementation of the fractional flow approach to multi-phase flow simulation'',
 *  Advances in water resources, vol. 22, no. 5, pp. 461-478, 1999.),
 * \f$ p_c = p_n - p_w \f$ is the capillary pressure, \f$ \boldsymbol K \f$ the absolute permeability,
 * \f$ \lambda = \lambda_w +  \lambda_n \f$ the total mobility depending on the
 * saturation (\f$ \lambda_\alpha = k_{r_\alpha} / \mu_\alpha \f$),
 * \f$ f_\alpha = \lambda_\alpha / \lambda \f$ the fractional flow function of a phase,
 * \f$ \rho_\alpha \f$ a phase density, \f$ g \f$ the gravity constant and \f$ q \f$ the source term.
 *
 * For all cases, \f$ p = p_D \f$ on \f$ \Gamma_{Dirichlet} \f$, and \f$ \boldsymbol v_{total} \cdot  \boldsymbol n  = q_N \f$
 * on \f$ \Gamma_{Neumann} \f$.
 *
 * The slightly compressible case is only implemented for phase pressures! In this case for a wetting
 * \f$(w) \f$ phase pressure as primary variable the equations are formulated as
 * \f[
 * \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,
 * \left[\lambda \boldsymbol{K} \left(\textbf{grad}\, p_w + f_n \, \textbf{grad}\, p_c + \sum f_\alpha \rho_\alpha \,
 * g \, \textbf{grad}\, z\right)\right] = q,
 * \f]
 * and for a non-wetting (\f$ n \f$) phase pressure as
 *  \f[
 *  \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,
 * \left[\lambda \boldsymbol{K}  \left(\textbf{grad}\, p_n - f_w \textbf{grad}\, p_c + \sum f_\alpha \rho_\alpha \,
 * g \, \textbf{grad}\, z\right)\right] = q,
 *  \f]
 * In this slightly compressible case the following definitions are valid:
 * \f$ \lambda = \rho_w \lambda_w + \rho_n \lambda_n \f$, \f$ f_\alpha = (\rho_\alpha \lambda_\alpha) / \lambda \f$
 * This model assumes that temporal changes in density are very small and thus terms of temporal derivatives are negligible in the pressure equation.
 * Depending on the formulation the terms including time derivatives of saturations are simplified by inserting  \f$ S_w + S_n = 1 \f$.
 *
 *  In the IMPES models the default setting is:
 *
 *  - formulation: \f$ p_w-S_w \f$ (Property: \a Formulation defined as \a DecoupledTwoPCommonIndices::pwsw)
 *
 *  - compressibility: disabled (Property: \a EnableCompressibility set to \a false)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2PVE: public FVPressure2P<TypeTag>
{
    typedef FVPressure2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, PressureSolutionVector) PressureSolutionVector;

    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        sw = Indices::saturationW,
        sn = Indices::saturationNw,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressureEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    enum VEModel
    {
        sharpInterface,
        capillaryFringe
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

protected:
    //! \cond \private
    typedef typename ParentType::EntryType EntryType;
    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };
    //! \endcond

public:
    // Function which calculates the boundary flux entry
    void getFluxOnBoundary(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first);

//    // updates and stores constitutive relations
    void updateMaterialLaws();

    /*! \brief Initializes the pressure model
     *
     * \copydetails ParentType::initialize()
     *
     * \param solveTwice indicates if more than one iteration is allowed to get an initial pressure solution
     */
    void initialize(bool solveTwice = true)
    {
        if (!compressibility_)
        {
            const auto element = *problem_.gridView().template begin<0>();
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
            fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
            fluidState.setTemperature(problem_.temperature(element));
            fluidState.setSaturation(wPhaseIdx, 1.);
            fluidState.setSaturation(nPhaseIdx, 0.);
            density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
            density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
            viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
            viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
        }

        ParentType::initialize();
    }

    /*!\brief Returns the reconstructed phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPressure(Scalar height, int phaseIdx, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar coarsePressureW = cellData.pressure(wPhaseIdx);
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = cellData.gasPlumeDist();

        Scalar reconstPressure[numPhases];
        reconstPressure[wPhaseIdx] = coarsePressureW;//reconstruct phase pressures for no ve model
        reconstPressure[nPhaseIdx] = cellData.pressure(nPhaseIdx);
        if(veModel_ == sharpInterface && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx];
        }
        else if(veModel_ == capillaryFringe && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx] + entryP_;
        }
        else if (veModel_ == sharpInterface && height > gasPlumeDist) //reconstruct non-wetting phase pressure for sharp interface ve model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist);
            reconstPressure[nPhaseIdx]= reconstPressure[wPhaseIdx];
        }
        else if (veModel_ == capillaryFringe && height > gasPlumeDist) //reconstruct non-wetting phase pressure for capillary fringe model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist) + entryP_;
        }

        return reconstPressure[phaseIdx];
    }

    /*!\brief Returns the reconstructed capillary phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstCapillaryPressure(Scalar height, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = cellData.gasPlumeDist();

        Scalar reconstCapillaryPressure = 0.0;//reconstruct capillary pressure for sharp interface ve model
        if (veModel_ == capillaryFringe && height <= gasPlumeDist) //reconstruct capillary pressure for capillary fringe model
        {
            reconstCapillaryPressure = entryP_;
        }
        else if (veModel_ == capillaryFringe && height > gasPlumeDist) //reconstruct capillary pressure for capillary fringe model
        {
            reconstCapillaryPressure = densityW * problem_.gravity().two_norm() * (height-gasPlumeDist)
                    + entryP_ - densityNw * problem_.gravity().two_norm() * (height-gasPlumeDist);
        }

        return reconstCapillaryPressure;
    }

    /*!\brief Returns the reconstructed phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstSaturation(Scalar height, int phaseIdx, const Element element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar minGasPlumeDist = cellData.minGasPlumeDist();
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(residualSegSaturation(element), cellData.saturation(wPhaseIdx));
        Scalar resSatN = problem_.spatialParams().materialLawParams(element).snr();

        Scalar reconstSaturation[numPhases];
        if (veModel_ == sharpInterface) //reconstruct phase saturation for sharp interface ve model
        {
            reconstSaturation[wPhaseIdx] = 1.0;
            reconstSaturation[nPhaseIdx] = 0.0;
            if(height > minGasPlumeDist && height < gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = 1.0-resSatN;
                reconstSaturation[nPhaseIdx] = resSatN;
            }
            else if(height >= gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = resSatW;
                reconstSaturation[nPhaseIdx] = 1.0-resSatW;
            }
        }
        else if (veModel_ == capillaryFringe) //reconstruct phase saturation for capillary fringe model
        {
            Scalar densityW = density_[wPhaseIdx];
            Scalar densityNw = density_[nPhaseIdx];

            reconstSaturation[wPhaseIdx] = 1.0;
            reconstSaturation[nPhaseIdx] = 0.0;
            if(height > minGasPlumeDist && height < gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = 1.0-resSatN;
                reconstSaturation[nPhaseIdx] = resSatN;
            }
            else if(height >= gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = std::pow(((height - gasPlumeDist) * (densityW - densityNw) * problem_.gravity().two_norm() + entryP_), (-lambda_))
                * std::pow(entryP_, lambda_) * (1.0 - resSatW - resSatN) + resSatW;
                reconstSaturation[nPhaseIdx] = 1.0 - reconstSaturation[wPhaseIdx];
            }
        }
        else //reconstruct phase saturation for no ve model
        {
            reconstSaturation[wPhaseIdx] = cellData.saturation(wPhaseIdx);
            reconstSaturation[nPhaseIdx] = cellData.saturation(nPhaseIdx);
        }
        return reconstSaturation[phaseIdx];
    }

    /*! \brief Calculates wetting saturation integral over height given two bounds
     */
    Scalar saturationIntegral(Scalar lowerBound, Scalar upperBound, const Element element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar minGasPlumeDist = cellData.minGasPlumeDist();
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(residualSegSaturation(element), cellData.saturation(wPhaseIdx));
        Scalar resSatN = problem_.spatialParams().materialLawParams(element).snr();

        //different integration areas: 1. from bottom to minGasPlumeDist, 2. from minGasPlumeDist to gasPlumeDist, 3. from gasPlumeDist to top
        Scalar integral = 0.0;

        if(lowerBound < minGasPlumeDist)
        {
            Scalar upperPartBound = std::min(minGasPlumeDist, upperBound);
            integral += upperPartBound - lowerBound;
        }
        if(upperBound > minGasPlumeDist && lowerBound < gasPlumeDist )
        {
            Scalar lowerpartBound = std::max(minGasPlumeDist, lowerBound);
            Scalar upperPartBound = std::min(gasPlumeDist, upperBound);
            integral += (upperPartBound - lowerpartBound)*(1.0-resSatN);
        }
        if(upperBound > gasPlumeDist)
        {
            Scalar lowerPartBound = std::max(gasPlumeDist, lowerBound);
            if(veModel_ == sharpInterface)
            {
            integral += (upperBound - lowerPartBound)*resSatW;
            }

            if(veModel_ == capillaryFringe)//numerical integration
            {
                int intervalNumber = 10;
                Scalar deltaZ = (upperBound - lowerPartBound)/intervalNumber;
                for(int count=0; count<intervalNumber; count++ )
                {
                integral += (reconstSaturation(lowerPartBound + count*deltaZ, wPhaseIdx, element) +
                            reconstSaturation(lowerPartBound + (count+1)*deltaZ, wPhaseIdx, element))/2.0*deltaZ;
                }
            }
        }
        return integral;
    }


    /*! \brief Calculates new residual saturation for VE-model depending on segregation time
     */
    Scalar residualSegSaturation(const Element& element)
    {
//        Scalar residualSegSaturation = problem_.residualSegSaturation();
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar residualSegSaturation = cellData.residualSegSaturation();
        return residualSegSaturation;
    }

    /*! \brief Calculates new residual saturation for VE-model depending on segregation time
     */
    void calculateResidualSegSaturation(Element element)
    {
        Scalar residualSegSat = 0.0;
        Scalar gasPlumeTip = 0.0; //dummy

        int correction = 0;
        if (ParameterTree::tree().hasKey("VE.correction"))
        {
            correction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, correction);
        }

        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        if(correction == 1 || correction == 2 || correction == 3 || correction == 4)
        {
            GlobalPosition globalPos = element.geometry().center();
            Scalar aquiferHeight = problem_.bBoxMax()[dim - 1];
            Scalar charHeight = aquiferHeight;
            Scalar porosity = problem_.spatialParams().porosity(element);
            Scalar viscosityW = viscosity_[wPhaseIdx];
            Scalar viscosityNw = viscosity_[nPhaseIdx];
            Scalar time = problem_.timeManager().time() + problem_.timeManager().timeStepSize();
            Scalar permeability = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, PermeabilityVertical);
            Scalar densityW = density_[wPhaseIdx];
            Scalar densityN = density_[nPhaseIdx];
            Scalar gravity = gravity_.two_norm();
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            Scalar resSatN = problem_.spatialParams().materialLawParams(element).snr();
            Scalar exponent = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Exponent);
            Scalar lambda = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LambdaKr);
            int relPermModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, Model);

            if(correction == 2)//correct the domain height with analytic solution
            {
                charHeight = viscosityNw/viscosityW * aquiferHeight;
            }
            else if(correction == 3)//correct the domain height with numerical solution
            {
                Scalar gasVolume = -GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, Scalar, "BoundaryConditions", Injectionrate) * aquiferHeight
                                  * time / (densityN*porosity);
                charHeight = gasVolume/gasPlumeTip;
            }
            if(isnan(charHeight) || isinf(charHeight))
                charHeight = 0.0;

            if(correction == 4)//calculate saturation explicitly
            {
                Scalar sat = cellData.saturation(wPhaseIdx);
                Scalar satOld = cellData.oldSaturation();
                residualSegSat = cellData.residualSegSaturation();
                Scalar gasPlumeDist = cellData.gasPlumeDist();
                Scalar domainHeight = problem_.bBoxMax()[dim - 1];
                Scalar corrGasPlumeHeight = 0.0;
                Scalar CTZ = problem_.capillaryTransitionZone();

                if(gasPlumeDist + CTZ < 0.0)//all of column is gas plume not including capillary fringe
                    corrGasPlumeHeight = domainHeight;
                else if((domainHeight - (gasPlumeDist + CTZ)) > 0.0)//otherwise all is capillary fringe
                {
                    corrGasPlumeHeight = domainHeight - (gasPlumeDist + CTZ);
                }

                Scalar timeStepSize = problem_.timeManager().timeStepSize();
                Scalar mobilityW = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), residualSegSat)/viscosityW;
                Scalar mobilityN = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), residualSegSat)/viscosityNw;
                Scalar drainageFlux = mobilityW*mobilityN/(mobilityW+mobilityN)*permeability*gravity*(densityW-densityN);

                Scalar residualSegSatAverage = problem_.residualSegSatAverage();

                if(residualSegSat > sat)
                {
                    residualSegSat = sat;
                    corrGasPlumeHeight = domainHeight;
                }

                if(satOld > 1.0 - 1e-6)
                {
                    residualSegSat = residualSegSatAverage;
                }

                if(problem_.timeManager().timeStepIndex() < 2)//no segregation in first time step
                {
                    if(corrGasPlumeHeight > 1e-6)
                        residualSegSat = sat;
                    else
                        residualSegSat = 0.9;
                }
                else if(corrGasPlumeHeight < 1e-6 || sat > 1.0 - 1e-6)
                {
                    residualSegSat = residualSegSatAverage;
                }
                else if(drainageFlux > corrGasPlumeHeight * residualSegSat * (1.0-residualSegSat) * porosity/timeStepSize)
                {
                    std::cout << "attention: time step too large local vertical!!" << std::endl;
                    std::cout << "drainageFlux " << drainageFlux << std::endl;
                    std::cout << "corrGasPlumeHeight " << corrGasPlumeHeight << std::endl;
                    std::cout << "residualSegSatvorher " << residualSegSat << std::endl;
                    std::cout << "residualSegSatnachher " << (corrGasPlumeHeight * porosity * residualSegSat - drainageFlux  * timeStepSize/(1.0-residualSegSat)) /
                            (corrGasPlumeHeight * porosity - drainageFlux * timeStepSize / (1.0-residualSegSat)) << std::endl;
                    std::cout << "residualSegSatinstead " << residualSegSatAverage << std::endl;
                    residualSegSat = residualSegSatAverage;
                }
                else
                {
                    residualSegSat = (corrGasPlumeHeight * porosity * residualSegSat - drainageFlux  * timeStepSize/(1.0-residualSegSat)) /
                            (corrGasPlumeHeight * porosity - drainageFlux * timeStepSize / (1.0-residualSegSat));
//                    std::cout << "residualSegSat " << residualSegSat << std::endl;
//                    std::cout << "corrGasPlumeHeight " << corrGasPlumeHeight << std::endl;
//                    std::cout << "gasPlumeDist " << gasPlumeDist << std::endl;
                }

                if(problem_.timeManager().timeStepIndex() == 500)
                    residualSegSat = residualSegSatAverage;

                //hack: for gasPlumeDist smaller than -CTZ, gasPlumeDist is set to -CTZ and a new resSatW according to that is calculated
                Scalar fullIntegralHack =
                        1.0 / (1.0 - lambda_) * std::pow(((domainHeight + CTZ) * (densityW - densityN) * gravity_.two_norm() + entryP_),(1.0 - lambda_))
                        * (1.0 - residualSegSat - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityN) * gravity_.two_norm())
                        + residualSegSat * domainHeight
                        - 1.0 / (1.0 - lambda_) * std::pow((entryP_ + CTZ*(densityW - densityN) * gravity_.two_norm()), (1.0-lambda_))
                        * (1.0 - residualSegSat - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityN) * gravity_.two_norm());
                if(fullIntegralHack > sat * domainHeight && sat < 1.0-1e-8)
                {
                       residualSegSat =
                            (sat * domainHeight
                            - 1.0/(1.0 - lambda_) * std::pow(((domainHeight + CTZ)*(densityW - densityN)*gravity_.two_norm() + entryP_), (1.0-lambda_))
                            * std::pow(entryP_, lambda_) * (1.0 - resSatN) / ((densityW - densityN) * gravity_.two_norm())
                            + 1.0/(1.0 - lambda_) * std::pow((CTZ * (densityW - densityN) * gravity_.two_norm() + entryP_), (1.0 - lambda_))
                            * (1.0 - resSatN)/((densityW - densityN) * gravity_.two_norm()) * std::pow(entryP_, lambda_))/
                            (- 1.0/(1.0 - lambda_) * std::pow(((domainHeight + CTZ) * (densityW - densityN) * gravity_.two_norm() + entryP_), (1.0 - lambda_))
                            * std::pow(entryP_, lambda_)/((densityW - densityN) * gravity_.two_norm()) + domainHeight
                            + 1.0/(1.0 - lambda_) * std::pow((entryP_ + CTZ * (densityW - densityN) * gravity_.two_norm()), (1.0-lambda_))
                            * std::pow(entryP_, lambda_)/((densityW - densityN) * gravity_.two_norm())) - 1e-6;
                }
            }
            else if(relPermModel == 0 || (relPermModel == 1 && exponent == 1))
            {
                residualSegSat = (charHeight*porosity*viscosityW/(time*permeability*(densityW-densityN)*gravity))*(1.0-resSatW-resSatN) + resSatW;
            }
            else
            {
                if(relPermModel == 2)
                {
                    exponent = 2.0/lambda + 3.0;
                }
                residualSegSat = std::pow((charHeight*porosity*viscosityW/(time*permeability*(densityW-densityN)*gravity)), (1.0/(exponent)))*
                        (1.0-resSatW-resSatN) + resSatW;
            }
        }
        else
        {
            residualSegSat = problem_.spatialParams().materialLawParams(element).swr();
            Scalar PseudoResidualSaturation = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, VE, PseudoResidualSaturation);
            residualSegSat = residualSegSat + PseudoResidualSaturation;
        }
        cellData.setResidualSegSaturation(residualSegSat);
    }



    /*! \brief Calculates gasPlumeDist, distance of gas plume from bottom
     *
     * Stores minGasPlumeDist for all grid cells
     */
    void calculateGasPlumeDist(const Element element, Scalar satW)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        Scalar minGasPlumeDist = cellData.minGasPlumeDist();
        Scalar domainHeight = problem_.bBoxMax()[dim - 1];
        if(minGasPlumeDist > domainHeight)
            minGasPlumeDist = domainHeight;
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(residualSegSaturation(element), cellData.saturation(wPhaseIdx));
        Scalar resSatN = problem_.spatialParams().materialLawParams(element).snr();

        Scalar gasPlumeDist = 0.0;

        if(satW > 1.0-1e-6)
            gasPlumeDist = domainHeight;
        else if (veModel_ == sharpInterface) //calculate gasPlumeDist for sharp interface ve model
        {
            //calculate gasPlumeDist for gasPlumeDist > minGasPlumeDist
            if (satW * domainHeight > minGasPlumeDist + (domainHeight - minGasPlumeDist) * resSatW + 1e-6)
            {
                gasPlumeDist = (domainHeight * (satW - resSatW) - minGasPlumeDist * resSatN) / (1.0 - resSatW - resSatN);
            }
            //calculate gasPlumeDist for gasPlumeDist < minGasPlumeDist and calculate new minGasPlumeDist
            else if (satW * domainHeight < minGasPlumeDist + (domainHeight - minGasPlumeDist) * resSatW - 1e-6)
            {
                gasPlumeDist = domainHeight * (satW - resSatW) / (1.0 - resSatW);
                cellData.setMinGasPlumeDist(gasPlumeDist);
            }
            //calculate gasPlumeDist for gasPlumeDist = minGasPlumeDist
            else
            {
                gasPlumeDist = minGasPlumeDist;
            }
        }

        else if (veModel_ == capillaryFringe) //calculate gasPlumeDist for capillary fringe model
        {
            Scalar densityW = density_[wPhaseIdx];
            Scalar densityNw = density_[nPhaseIdx];
            Scalar temp =
                   1.0 / (1.0 - lambda_) * std::pow(((domainHeight - minGasPlumeDist)* (densityW - densityNw) * gravity_.two_norm() + entryP_), (1.0 - lambda_))
                   * (1.0 - resSatW - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityNw) * gravity_.two_norm())
                   + resSatW * (domainHeight - minGasPlumeDist)
                   - entryP_ * (1.0 - resSatW - resSatN) / ((1.0 - lambda_) * (densityW - densityNw) * gravity_.two_norm())
                   + minGasPlumeDist;

            //calculate gasPlumeDist for gasPlumeDist > minGasPlumeDist
            Scalar residual = 1.0;
            if (satW * domainHeight > temp)
            {
                //GasPlumeDist>0 (always in this case)
                gasPlumeDist = domainHeight / 2.0; //XiStart
                //solve equation for
                int count = 0;
                for (;count < 100; count++)
                {
                    residual =
                    1.0 / (1.0 - lambda_) * std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_),(1.0 - lambda_))
                    * (1.0 - resSatW - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityNw) * gravity_.two_norm())
                    + resSatW * (domainHeight - gasPlumeDist)
                    - entryP_ * (1.0 - resSatW - resSatN) / ((1.0 - lambda_) * (densityW - densityNw) * gravity_.two_norm())
                    + minGasPlumeDist + (gasPlumeDist - minGasPlumeDist) * (1.0 - resSatN)
                    - satW * domainHeight;

                    if (fabs(residual) < 1e-10)
                    {
                        break;
                    }

                    Scalar derivation = std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_), -lambda_)
                    * (resSatN + resSatW - 1.0) * std::pow(entryP_, lambda_) - resSatW + 1.0 - resSatN;

                    gasPlumeDist = gasPlumeDist - residual / (derivation);
                }
            }
            //calculate gasPlumeDist for gasPlumeDist < minGasPlumeDist and calculate new minGasPlumeDist
            else if (satW * domainHeight < temp)
            {
                //check if GasPlumeDist>0, GasPlumeDist=0, GasPlumeDist<0
                Scalar fullIntegral =
                        1.0 / (1.0 - lambda_) * std::pow((domainHeight * (densityW - densityNw) * gravity_.two_norm() + entryP_),(1.0 - lambda_))
                        * (1.0 - resSatW - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityNw) * gravity_.two_norm())
                        + resSatW * domainHeight
                        - entryP_ * (1.0 - resSatW - resSatN) / ((1.0 - lambda_) * (densityW - densityNw) * gravity_.two_norm());
                //GasPlumeDist>0
                if (fullIntegral < satW * domainHeight)
                {
                    gasPlumeDist = domainHeight / 2.0; //XiStart
                    //solve equation for
                    int count = 0;
                    for (; count < 100; count++)
                    {
                        residual =
                        1.0 / (1.0 - lambda_) * std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_),(1.0 - lambda_))
                        * (1.0 - resSatW - resSatN) * std::pow(entryP_, lambda_) / ((densityW - densityNw) * gravity_.two_norm())
                        + resSatW * (domainHeight - gasPlumeDist)
                        - entryP_ * (1.0 - resSatW - resSatN) / ((1.0 - lambda_) * (densityW - densityNw) * gravity_.two_norm())
                        + gasPlumeDist
                        - satW * domainHeight;

                        if (fabs(residual) < 1e-12)
                        {
                            break;
                        }

                        Scalar derivation =
                               std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_), -lambda_)
                               * (- 1.0 + resSatN + resSatW) * std::pow(entryP_, lambda_) - resSatW + 1.0;

                        gasPlumeDist = gasPlumeDist - residual / (derivation);
                    }
                }
                //GasPlumeDist<0
                else if (fullIntegral > satW * domainHeight)
                {
                    gasPlumeDist = 0.0; //XiStart
                    //solve equation
                    int count = 0;
                    for (; count < 100; count++)
                    {
                        Scalar residual =
                               1.0 / (1.0 - lambda_)
                               * std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_), (1.0 - lambda_))
                               * (1.0 - resSatW - resSatN)* std::pow(entryP_, lambda_) / ((densityW - densityNw) * gravity_.two_norm())
                               + resSatW * domainHeight
                               - 1.0 / (1.0 - lambda_) * std::pow((entryP_ - gasPlumeDist * (densityW - densityNw)
                               * gravity_.two_norm()),(1.0 - lambda_)) * (1.0 - resSatW - resSatN) * std::pow(entryP_, lambda_)
                               / ((densityW - densityNw) * gravity_.two_norm())
                               - satW * domainHeight;

                        if (fabs(residual) < 1e-10)
                        {
                            break;
                        }

                        Scalar derivation =
                               std::pow(((domainHeight - gasPlumeDist) * (densityW - densityNw) * gravity_.two_norm() + entryP_), -lambda_)
                               * (- 1.0 + resSatN + resSatW) * std::pow(entryP_, lambda_)
                               + std::pow((entryP_ - gasPlumeDist * (densityW - densityNw) * gravity_.two_norm()), -lambda_)
                               * (1.0 - resSatN - resSatW) * std::pow(entryP_, lambda_);

                        gasPlumeDist = gasPlumeDist - residual / (derivation);
                    }
                }
                //GasPlumeDist=0
                else
                {
                    gasPlumeDist = 0.0;
                }
                cellData.setMinGasPlumeDist(std::max(0.0, gasPlumeDist));
            }
            //calculate gasPlumeDist for gasPlumeDist = minGasPlumeDist
            else
            {
                gasPlumeDist = minGasPlumeDist;
            }
        }
        cellData.setGasPlumeDist(gasPlumeDist);
    }

    /*! \brief Calculates relative permeability as numerical integral of fine-scale permeabilities over z
     *
     */
        Scalar calculateRelPermeabilityCoarse(Scalar lowerBound, Scalar upperBound, const Element element, int phaseIdx)
        {
            int eIdxGlobal = problem_.variables().index(element);
            CellData& cellData = problem_.variables().cellData(eIdxGlobal);
            Scalar gasPlumeDist = cellData.gasPlumeDist();

            Scalar permeabilityCoarse = 0.0;

            //Attention: This introduces a hard-coded lense
//            if(lowerBound < gasPlumeDist && phaseIdx == wPhaseIdx)
//            {
//                Scalar upperPartBound = std::min(gasPlumeDist, upperBound);
//                permeabilityCoarse += upperPartBound - lowerBound;
//            }
//            if(upperBound > gasPlumeDist)
//            {
                Scalar lowerPartBound = lowerBound;
                if(veModel_ == sharpInterface)
                {
                    //set residual saturations according to some criterium
                    Scalar resSatW = std::min(residualSegSaturation(element), cellData.saturation(wPhaseIdx));
                    switch (phaseIdx)
                    {
                    case wPhaseIdx:
                    {
                        permeabilityCoarse += (upperBound - lowerPartBound)*MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), resSatW);
                    }
                    break;
                    case nPhaseIdx:
                    {
                        permeabilityCoarse += (upperBound - lowerPartBound)*MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), resSatW);
                    }
                    break;
                    }
                }
                else if(veModel_ == capillaryFringe) //numerical integration
                {
                    int reconstruction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, Reconstruction);
                    int intervalNumber = std::pow(2, reconstruction);//TODO: find better solution, too slow
                    Scalar deltaZ = (upperBound - lowerPartBound)/intervalNumber;

                    Scalar saturationStart = reconstSaturation(lowerPartBound, wPhaseIdx, element);
                    Scalar saturationEnd = reconstSaturation(upperBound, wPhaseIdx, element);

                    GlobalPosition globalPos = element.geometry().center();
                    if (globalPos[0] > 100 && globalPos[0] < 120)//lense
                    {
                        Scalar permeabilityStart = 2e-12;
                        Scalar permeabilityEnd = 2e-15;
                    switch (phaseIdx)
                    {
                    case wPhaseIdx:
                    {
                        permeabilityCoarse += deltaZ / 2.0 * ((MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturationStart)
                        * permeabilityStart)
                        + (MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturationEnd))*permeabilityEnd);
                        for (int count = 1; count < intervalNumber; count++)
                        {
                            Scalar permeability = 2e-12;
                            if(deltaZ * count + lowerPartBound > 20)
                                permeability = 2e-15;
                            Scalar saturation = reconstSaturation(deltaZ * count + lowerPartBound, wPhaseIdx, element);
                            permeabilityCoarse += deltaZ * MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturation) * permeability;
                        }
                        permeabilityCoarse = permeabilityCoarse/1.334e-12;
                        break;
                    }
                    case nPhaseIdx:
                    {
                        permeabilityCoarse += deltaZ / 2.0 * ((MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturationStart)
                        * permeabilityStart)
                        + (MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturationEnd))*permeabilityEnd);
                        for (int count = 1; count < intervalNumber; count++)
                        {
                            Scalar permeability = 2e-12;
                            if(deltaZ * count + lowerPartBound > 20)
                                permeability = 2e-15;
                            Scalar saturation = reconstSaturation(deltaZ * count + lowerPartBound, wPhaseIdx, element);
                            permeabilityCoarse += deltaZ * MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturation) * permeability;
                        }
                        permeabilityCoarse = permeabilityCoarse/1.334e-12;
                        break;
                    }
                    }
                    }
                    else//no lense
                    {
                    switch (phaseIdx)
                    {
                    case wPhaseIdx:
                    {
                        permeabilityCoarse += deltaZ / 2.0 * (MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturationStart)
                        + MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturationEnd));
                        for (int count = 1; count < intervalNumber; count++)
                        {
                            Scalar saturation = reconstSaturation(deltaZ * count + lowerPartBound, wPhaseIdx, element);
                            permeabilityCoarse += deltaZ * MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), saturation);
                        }
                        break;
                    }
                    case nPhaseIdx:
                    {
                        permeabilityCoarse += deltaZ / 2.0 * (MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturationStart)
                        + MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturationEnd));
                        for (int count = 1; count < intervalNumber; count++)
                        {
                            Scalar saturation = reconstSaturation(deltaZ * count + lowerPartBound, wPhaseIdx, element);
                            permeabilityCoarse += deltaZ * MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), saturation);
                        }
                        break;
                    }
                    }
                    }
                }
//            }
            return (permeabilityCoarse / (upperBound - lowerBound));
        }

    //! Constructs a FVPressure2P object
    /**
     * \param problem A problem class object
     */
    FVPressure2PVE(Problem& problem) : ParentType(problem), problem_(problem), gravity_(problem.gravity())
    {
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (pressureType_ == pGlobal && compressibility_)
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported for global pressure!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, Vtk, OutputLevel);
        veModel_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
        lambda_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Lambda);
        entryP_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, EntryPressure);
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;
    int veModel_;
    Scalar lambda_;
    Scalar entryP_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    //! gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
    //! gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
};

/*! \brief Function which calculates the flux entry at a boundary
 *
 * \copydetails FVPressure::getFluxOnBoundary(EntryType&,const Intersection&,const CellData&,const bool)
 *
 * Dirichlet boundary condition for pressure depends on the formulation (\f$p_w\f$ (default), \f$p_n\f$, \f$p_{global}\f$),
 * Neumann boundary condition are the phase mass fluxes (\f$q_w\f$ and \f$q_n\f$, [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVPressure2PVE<TypeTag>::getFluxOnBoundary(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first)
{
    Element element = intersection.inside();

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = element.geometry().center();

    // center of face in global coordinates
    const GlobalPosition& globalPosJ = intersection.geometry().center();

    // get mobilities and fractional flow factors
    Scalar lambdaWI = cellData.mobility(wPhaseIdx);
    Scalar lambdaNwI = cellData.mobility(nPhaseIdx);
    Scalar fractionalWI = cellData.fracFlowFunc(wPhaseIdx);
    Scalar fractionalNwI = cellData.fracFlowFunc(nPhaseIdx);

    // get capillary pressure
    Scalar pcI = cellData.capillaryPressure();

    //get face index
    int isIndexI = intersection.indexInInside();

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal =
            intersection.centerUnitOuterNormal();

    // get face area
    Scalar faceArea = intersection.geometry().volume();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    BoundaryTypes bcType;
    problem_.boundaryTypes(bcType, intersection);
    PrimaryVariables boundValues(0.0);

    if (bcType.isDirichlet(eqIdxPress))
    {
        problem_.dirichlet(boundValues, intersection);

        //permeability vector at boundary
        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability,
                problem_.spatialParams().intrinsicPermeability(element));

        Dune::FieldVector < Scalar, dim > permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
        Scalar satW = 0;
        Scalar satNw = 0;
        if (bcType.isDirichlet(eqIdxSat))
        {
            switch (saturationType_)
            {
            case sw:
            {
                satW = boundValues[saturationIdx];
                satNw = 1.0 - boundValues[saturationIdx];
                break;
            }
            case sn:
            {
                satW = 1.0 - boundValues[saturationIdx];
                satNw = boundValues[saturationIdx];
                break;
            }
            }
        }
        else
        {
            satW = cellData.saturation(wPhaseIdx);
            satNw = cellData.saturation(nPhaseIdx);
        }
        Scalar temperature = problem_.temperature(element);

        //get dirichlet pressure boundary condition
        Scalar pressBound = boundValues[pressureIdx];

        //calculate constitutive relations depending on the kind of saturation used
        Scalar pcBound = MaterialLaw::pc(problem_.spatialParams().materialLawParams(element), satW);
        if (veModel_ == sharpInterface)
        {
            //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
            //sharp interface injection cells
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            pcBound = gravity_.two_norm() * problem_.bBoxMax()[dim - 1] * (satW - resSatW) / (1.0 - resSatW) * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            //minGasPlume = init
        }
        else if (veModel_ == capillaryFringe)
        {
            //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
            //sharp interface injection cells (with entry pressure)
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            pcBound = gravity_.two_norm() * problem_.bBoxMax()[dim - 1] * (satW - resSatW) / (1.0 - resSatW) * (density_[nPhaseIdx] - density_[wPhaseIdx]) + entryP_;
            //minGasPlume = init
        }

        //determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNw = 0;
        if (pressureType_ == pw)
        {
            pressW = pressBound;
            pressNw = pressBound + pcBound;
        }
        else if (pressureType_ == pn)
        {
            pressW = pressBound - pcBound;
            pressNw = pressBound;
        }

        Scalar densityWBound = density_[wPhaseIdx];
        Scalar densityNwBound = density_[nPhaseIdx];
        Scalar viscosityWBound = viscosity_[wPhaseIdx];
        Scalar viscosityNwBound = viscosity_[nPhaseIdx];
        Scalar rhoMeanW = 0;
        Scalar rhoMeanNw = 0;

        if (compressibility_)
        {
            FluidState fluidState;
            fluidState.setSaturation(wPhaseIdx, satW);
            fluidState.setSaturation(nPhaseIdx, satNw);
            fluidState.setTemperature(temperature);
            fluidState.setPressure(wPhaseIdx, pressW);
            fluidState.setPressure(nPhaseIdx, pressNw);

            densityWBound = FluidSystem::density(fluidState, wPhaseIdx);
            densityNwBound = FluidSystem::density(fluidState, nPhaseIdx);
            viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx)
                    / densityWBound;
            viscosityNwBound = FluidSystem::viscosity(fluidState, nPhaseIdx)
                    / densityNwBound;

            rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + densityWBound);
            rhoMeanNw = 0.5 * (cellData.density(nPhaseIdx) + densityNwBound);
        }

        Scalar lambdaWBound = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), satW) / viscosityWBound;
        Scalar lambdaNwBound = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), satW) / viscosityNwBound;
        //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
        //sharp interface injection cells
        if (veModel_ == sharpInterface || veModel_ == capillaryFringe)
        {
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            lambdaWBound = (satW - resSatW) / (1.0-resSatW) / viscosityWBound;
            lambdaNwBound = ((1.0 - satW) / (1.0-problem_.spatialParams().materialLawParams(element).swr())) / viscosityNwBound; //TODO: compressibility (like in Bo's model)
        }

        Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNwBound);
        Scalar fractionalNwBound = lambdaNwBound
                / (lambdaWBound + lambdaNwBound);

        Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
        Scalar fMeanNw = 0.5 * (fractionalNwI + fractionalNwBound);

        Scalar potentialDiffW = 0;
        Scalar potentialDiffNw = 0;

        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];

        if (!first)
        {
            potentialDiffW = cellData.fluxData().upwindPotential(wPhaseIdx, isIndexI);
            potentialDiffNw = cellData.fluxData().upwindPotential(nPhaseIdx, isIndexI);

            if (compressibility_)
            {
                densityW = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
                densityNw = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : densityNwBound;

                densityW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : densityW;
                densityNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : densityNw;
            }

            //calculate potential gradient
            switch (pressureType_)
            {
            case pw:
            {
                potentialDiffW = (cellData.pressure(wPhaseIdx) - pressBound);
                potentialDiffNw = (cellData.pressure(nPhaseIdx) - pressBound
                        - pcBound);
                break;
            }
            case pn:
            {
                potentialDiffW = (cellData.pressure(wPhaseIdx) - pressBound
                        + pcBound);
                potentialDiffNw = (cellData.pressure(nPhaseIdx) - pressBound);
                break;
            }
            case pGlobal:
            {
                potentialDiffW = (cellData.globalPressure() - pressBound
                        - fMeanNw * (pcI - pcBound));
                potentialDiffNw = (cellData.globalPressure() - pressBound
                        + fMeanW * (pcI - pcBound));
                break;
            }
            }

            potentialDiffW += densityW * (distVec * gravity_);
            potentialDiffNw += densityNw * (distVec * gravity_);
        }

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWBound;
        lambdaW =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(
                        potentialDiffW, 0.0, 1.0e-30)) ?
                        0.5 * (lambdaWI + lambdaWBound) : lambdaW;
        Scalar lambdaNw = (potentialDiffNw > 0.) ? lambdaNwI : lambdaNwBound;
        lambdaNw =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(
                        potentialDiffNw, 0.0, 1.0e-30)) ?
                        0.5 * (lambdaNwI + lambdaNwBound) : lambdaNw;

        if (compressibility_)
        {
            densityW = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
            densityW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : densityW;
            densityNw = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : densityNwBound;
            densityNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : densityNw;
        }

        Scalar scalarPerm = permeability.two_norm();
        //calculate current matrix entry
        entry[matrix] = (lambdaW + lambdaNw) * scalarPerm / dist * faceArea;
        entry[rhs] = entry[matrix] * pressBound;

        //calculate right hand side
        //calculate unit distVec
        distVec /= dist;
        Scalar areaScaling = (unitOuterNormal * distVec);
        //this treatment of g allows to account for gravity flux through faces where the face normal has no z component (e.g. parallelepiped grids)
        entry[rhs] -= (lambdaW * densityW + lambdaNw * densityNw) * scalarPerm
                * (gravity_ * distVec) * faceArea * areaScaling;

        if (pressureType_ == pw)
        {
            //add capillary pressure term to right hand side
            entry[rhs] -= 0.5 * (lambdaNwI + lambdaNwBound) * scalarPerm * (pcI - pcBound) / dist * faceArea;
        }
        else if (pressureType_ == pn)
        {
            //add capillary pressure term to right hand side
            entry[rhs] += 0.5 * (lambdaWI + lambdaWBound) * scalarPerm * (pcI - pcBound) / dist * faceArea;
        }
    }
    //set neumann boundary condition
    else if (bcType.isNeumann(eqIdxPress))
    {
        problem_.neumann(boundValues, intersection);

        if (!compressibility_)
        {
            boundValues[wPhaseIdx] /= density_[wPhaseIdx];
            boundValues[nPhaseIdx] /= density_[nPhaseIdx];
        }
        entry[rhs] = -(boundValues[wPhaseIdx] + boundValues[nPhaseIdx]) * faceArea;
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
    }
}

/*! \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 * In the compressible case additionally the densities and viscosities are stored.
 */
template<class TypeTag>
void FVPressure2PVE<TypeTag>::updateMaterialLaws()
{
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        Scalar temperature = problem_.temperature(element);

        //determine phase saturations from primary saturation variable
        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar satNw = cellData.saturation(nPhaseIdx);

        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];

        Scalar viscosityW = viscosity_[wPhaseIdx];
        Scalar viscosityNw = viscosity_[nPhaseIdx];

        //set pseudo residual segSaturation inside plume
        calculateResidualSegSaturation(element);
        //set gas plume distance
        calculateGasPlumeDist(element, satW);
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        if(isnan(gasPlumeDist))
            std::cout << "attention: gasPlumeDist isnan" << std::endl;

        Scalar pc = MaterialLaw::pc(problem_.spatialParams().materialLawParams(element), satW);
        if (veModel_ == sharpInterface)
        {
            pc = gravity_.two_norm() * gasPlumeDist * (densityNw - densityW);
        }
        else if (veModel_ == capillaryFringe)
        {
            pc = gravity_.two_norm() * gasPlumeDist * (densityNw - densityW) + entryP_; //TODO: compressibility (for the Dumux model, not necessary in Bo's slightly compressible model)
        }

        //determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNw = 0;
        if (pressureType_ == pw)
        {
            pressW = cellData.pressure(wPhaseIdx);
            pressNw = pressW + pc;
        }
        else if (pressureType_ == pn)
        {
            pressNw = cellData.pressure(nPhaseIdx);
            pressW = pressNw - pc;
        }

        if (compressibility_)
        {
            FluidState& fluidState = cellData.fluidState();
            fluidState.setTemperature(temperature);

            fluidState.setPressure(wPhaseIdx, pressW);
            fluidState.setPressure(nPhaseIdx, pressNw);

            densityW = FluidSystem::density(fluidState, wPhaseIdx);
            densityNw = FluidSystem::density(fluidState, nPhaseIdx);

            viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
            viscosityNw = FluidSystem::viscosity(fluidState, nPhaseIdx);

            //store density
            fluidState.setDensity(wPhaseIdx, densityW);
            fluidState.setDensity(nPhaseIdx, densityNw);

            //store viscosity
            fluidState.setViscosity(wPhaseIdx, viscosityW);
            fluidState.setViscosity(nPhaseIdx, viscosityNw);
        }
        else
        {
            cellData.setCapillaryPressure(pc);

            if (pressureType_ != pGlobal)
            {
                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressNw);
            }
        }

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), satW) / viscosityW;
        Scalar mobilityNw = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), satW) / viscosityNw;
        if (veModel_ == sharpInterface || veModel_ == capillaryFringe)
        {
            mobilityW = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, wPhaseIdx) / viscosityW;
            mobilityNw = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, nPhaseIdx) / viscosityNw;  //TODO: compressibility (like in Bo's model)
        }
        else if (veModel_ == capillaryFringe)
        {
            mobilityW = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, wPhaseIdx) / viscosityW;
            mobilityNw = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, nPhaseIdx) / viscosityNw;
        }

        if (compressibility_)
        {
            mobilityW *= densityW;
            mobilityNw *= densityNw;
        }

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNw);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));

        GlobalPosition globalPos = element.geometry().center();
        //pressure and potential are defined at the bottom of the domain
        if (veModel_ == sharpInterface || veModel_ == capillaryFringe)
        {
            globalPos[dim-1] = 0.0;
        }
        Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;

        Scalar potW = pressW + gravityDiff * densityW;
        Scalar potNw = pressNw + gravityDiff * densityNw;

        if (pressureType_ == pGlobal)
        {
            potW = cellData.globalPressure() - cellData.fracFlowFunc(nPhaseIdx) * pc + gravityDiff * densityW;
            potNw = cellData.globalPressure() - cellData.fracFlowFunc(wPhaseIdx) * pc + gravityDiff * densityNw;
        }

        cellData.setPotential(wPhaseIdx, potW);
        cellData.setPotential(nPhaseIdx, potNw);
    }
}

}
#endif
