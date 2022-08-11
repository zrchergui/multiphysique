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
#ifndef DUMUX_FVPRESSURE2PVEMULTIDIM_HH
#define DUMUX_FVPRESSURE2PVEMULTIDIM_HH

#include <dune/common/float_cmp.hh>

// dumux environment
#include "fvpressure2pve.hh"

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
template<class TypeTag> class FVPressure2PVEMultiDim: public FVPressure2PVE<TypeTag>
{
    typedef FVPressure2PVE<TypeTag> ParentType;

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
    typedef typename CellData::AdaptedValues AdaptedValues;
    typedef typename GET_PROP_TYPE(TypeTag, PressureSolutionVector) PressureSolutionVector;

    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

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
    typedef std::array<unsigned int, dim> CellArray;

protected:
    //! \cond \private
    typedef typename ParentType::EntryType EntryType;
    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };
    //! \endcond

public:
    /*!\brief Function which assembles the system of equations to be solved
     *
     *  This function assembles the Matrix and the right hand side (RHS) vector to solve for
     * a pressure field with a Finite-Volume (FV) discretization.
     * Implementations of this base class have to provide the methods <tt>getSource()</tt>,
     * <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
     *
     * \param first Indicates if function is called at the initialization step or during the simulation
     * (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
     */
    void assemble(bool first);

    // Function which calculates the flux entry
    void getFlux(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first);

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

    void update()
    {
        int gridSize = problem_.gridView().size(0);
        // update RHS vector, matrix
        if (problem_.gridAdapt().wasAdapted())
        {
            this->A_.setSize(gridSize, gridSize); //
            this->f_.resize(gridSize);
            this->pressure().resize(gridSize);


            for (int i = 0; i < gridSize; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);

                switch (pressureType_)
                {
                case pw:
                    this->pressure()[i] = cellData.pressure(wPhaseIdx);
                    break;
                case pn:
                    this->pressure()[i] = cellData.pressure(nPhaseIdx);
                    break;
                case pGlobal:
                    this->pressure()[i] = cellData.globalPressure();
                    break;
                }
            }

            ParentType::initializeMatrix();
        }


        ParentType::update();

//        calculateVelocity();

        return;
    }

    /*! \brief Stores the pressure solution of a cell
     *
     * Calculates secondary pressure variables and stores pressures.
     *
     * \param element Grid element
     */
    void storePressureSolution(const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        if (compressibility_)
        {
            density_[wPhaseIdx] = cellData.density(wPhaseIdx);
            density_[nPhaseIdx] = cellData.density(nPhaseIdx);
        }

        switch (pressureType_)
        {
        case pw:
        {
            Scalar pressW = this->pressure()[eIdxGlobal];
            Scalar pc = cellData.capillaryPressure();

            Scalar oldPressureW = cellData.pressure(wPhaseIdx);
            Scalar oldPressureN = cellData.pressure(nPhaseIdx);

            cellData.setOldPressure(wPhaseIdx, oldPressureW);
            cellData.setOldPressure(nPhaseIdx, oldPressureN);

            cellData.setPressure(wPhaseIdx, pressW);
            cellData.setPressure(nPhaseIdx, pressW + pc);

            GlobalPosition globalPos = element.geometry().center();

            int veModel = cellData.veModel();
            if (veModel == sharpInterface || veModel == capillaryFringe)
            {
                globalPos[dim-1] = 0.0;
            }
            Scalar gravityDiff = (problem_.bBoxMax() - globalPos) * gravity_;
            Scalar potW = pressW + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = pressW + pc + gravityDiff * density_[nPhaseIdx];

            //Attention: is this not changing the potential and thus probably the upwinding between
            //pressure and saturation step? -> unphysical behavior could occur
            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        case pn:
        {
            Scalar pressNw = this->pressure()[eIdxGlobal];
            Scalar pc = cellData.capillaryPressure();

            cellData.setPressure(nPhaseIdx, pressNw);
            cellData.setPressure(wPhaseIdx, pressNw - pc);

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potW = pressNw - pc + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = pressNw + gravityDiff * density_[nPhaseIdx];

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        case pGlobal:
        {
            Scalar press = this->pressure()[eIdxGlobal];
            cellData.setGlobalPressure(press);

            Scalar pc = cellData.capillaryPressure();
            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;

            //This is only an estimation!!! -> only used for uwind directions and time-step estimation!!!
            Scalar potW = press - cellData.fracFlowFunc(nPhaseIdx) * pc + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = press - cellData.fracFlowFunc(wPhaseIdx) * pc + gravityDiff * density_[nPhaseIdx];

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        }
        cellData.fluxData().resetVelocity();
    }

    /*!\brief Returns the reconstructed phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPressure(Scalar height, int phaseIdx, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar coarsePressureW = cellData.pressure(wPhaseIdx);
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();

        Scalar reconstPressure[numPhases];
        reconstPressure[wPhaseIdx] = coarsePressureW;//reconstruct phase pressures for no ve model
        reconstPressure[nPhaseIdx] = cellData.pressure(nPhaseIdx);
        if(veModel == sharpInterface && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx];
        }
        else if(veModel == capillaryFringe && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx] + entryP;
        }
        else if (veModel == sharpInterface && height > gasPlumeDist) //reconstruct non-wetting phase pressure for sharp interface ve model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist);
            reconstPressure[nPhaseIdx]= reconstPressure[wPhaseIdx];
        }
        else if (veModel == capillaryFringe && height > gasPlumeDist) //reconstruct non-wetting phase pressure for capillary fringe model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist) + entryP;
        }

        return reconstPressure[phaseIdx];
    }

    /*!\brief Returns the reconstructed phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPressure(Scalar coarsePressureW, Scalar height, int phaseIdx, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();

        Scalar reconstPressure[numPhases];
        if(veModel == sharpInterface && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx];
        }
        else if(veModel == capillaryFringe && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx] + entryP;
        }
        else if (veModel == sharpInterface && height > gasPlumeDist) //reconstruct non-wetting phase pressure for sharp interface ve model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist);
            reconstPressure[nPhaseIdx]= reconstPressure[wPhaseIdx];
        }
        else if (veModel == capillaryFringe && height > gasPlumeDist) //reconstruct non-wetting phase pressure for capillary fringe model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist) + entryP;
        }

        return reconstPressure[phaseIdx];
    }

    /*!\brief Returns the reconstructed phase pressure based on a container holding information about the coarse cell
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPressure(Scalar height, int phaseIdx, const AdaptedValues& valuesCoarse, const Element& elementCoarse)
    {
        int veModel = valuesCoarse.veModel;
        Scalar coarsePressureW = valuesCoarse.pressW;
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = valuesCoarse.gasPlumeDist;
        Scalar entryP = problem_.spatialParams().materialLawParams(elementCoarse).pe();

        Scalar reconstPressure[numPhases];
        if(veModel == sharpInterface && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx];
        }
        else if(veModel == capillaryFringe && height <= gasPlumeDist)
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = reconstPressure[wPhaseIdx] + entryP;
        }
        else if (veModel == sharpInterface && height > gasPlumeDist) //reconstruct non-wetting phase pressure for sharp interface ve model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist);
            reconstPressure[nPhaseIdx]= reconstPressure[wPhaseIdx];
        }
        else if (veModel == capillaryFringe && height > gasPlumeDist) //reconstruct non-wetting phase pressure for capillary fringe model
        {
            reconstPressure[wPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * height;
            reconstPressure[nPhaseIdx] = coarsePressureW - densityW * problem_.gravity().two_norm() * gasPlumeDist
                    - densityNw * problem_.gravity().two_norm() * (height - gasPlumeDist) + entryP;
        }

        return reconstPressure[phaseIdx];
    }

    /*!\brief Returns the reconstructed phase potential based on a container holding information about the coarse cell
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPotential(Scalar height, int phaseIdx, const Element& element)
    {
        Scalar density = density_[phaseIdx];
        Scalar pressure = reconstPressure(height, phaseIdx, element);
        GlobalPosition globalPos = element.geometry().center();
        globalPos[dim-1] = height;

        Scalar reconstPotential = pressure + (problem_.bBoxMax() - globalPos) * gravity_ * density;

        return reconstPotential;
    }


    /*!\brief Returns the reconstructed phase potential based on a container holding information about the coarse cell
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstPotential(Scalar height, int phaseIdx, const AdaptedValues& valuesCoarse, const Element& elementCoarse)
    {
        Scalar density = density_[phaseIdx];
        Scalar pressure = reconstPressure(height, phaseIdx, valuesCoarse, elementCoarse);
        GlobalPosition globalPos = elementCoarse.geometry().center();
        globalPos[dim-1] = height;

        Scalar reconstPotential = pressure + (problem_.bBoxMax() - globalPos) * gravity_ * density;

        return reconstPotential;
    }

    /*!\brief Returns the reconstructed capillary phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstCapillaryPressure(Scalar height, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();

        Scalar reconstCapillaryPressure = 0.0;//reconstruct capillary pressure for sharp interface ve model
        if (veModel == capillaryFringe && height <= gasPlumeDist) //reconstruct capillary pressure for capillary fringe model
        {
            reconstCapillaryPressure = entryP;
        }
        else if (veModel == capillaryFringe && height > gasPlumeDist) //reconstruct capillary pressure for capillary fringe model
        {
            reconstCapillaryPressure = densityW * problem_.gravity().two_norm() * (height-gasPlumeDist)
                    + entryP - densityNw * problem_.gravity().two_norm() * (height-gasPlumeDist);
        }

        return reconstCapillaryPressure;
    }

    /*!\brief Returns the reconstructed phase saturation
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstSaturation(Scalar height, int phaseIdx, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar minGasPlumeDist = cellData.minGasPlumeDist();
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(this->residualSegSaturation(element), cellData.saturation(wPhaseIdx));
        Scalar resSatN = problem_.spatialParams().materialLawParams(element).snr();

        Scalar reconstSaturation[numPhases];
        if (veModel == sharpInterface) //reconstruct phase saturation for sharp interface ve model
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
        else if (veModel == capillaryFringe) //reconstruct phase saturation for capillary fringe model
        {
            Scalar densityW = density_[wPhaseIdx];
            Scalar densityNw = density_[nPhaseIdx];
            Scalar lambda = problem_.spatialParams().materialLawParams(element).lambda();
            Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();

            reconstSaturation[wPhaseIdx] = 1.0;
            reconstSaturation[nPhaseIdx] = 0.0;
            if(height > minGasPlumeDist && height < gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = 1.0-resSatN;
                reconstSaturation[nPhaseIdx] = resSatN;
            }
            else if(height >= gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = std::pow(((height - gasPlumeDist) * (densityW - densityNw) * problem_.gravity().two_norm() + entryP), (-lambda))
                * std::pow(entryP, lambda) * (1.0 - resSatW - resSatN) + resSatW;
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

    /*!\brief Returns the reconstructed phase saturation given a container holding coarse cell information
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar reconstSaturation(Scalar height, int phaseIdx, const AdaptedValues& valuesCoarse, const Element& elementCoarse)
    {
        int veModel = valuesCoarse.veModel;
        Scalar gasPlumeDist = valuesCoarse.gasPlumeDist;
        Scalar minGasPlumeDist = valuesCoarse.minGasPlumeDist;
        Scalar coarseSatW = valuesCoarse.saturationW;
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(this->residualSegSaturation(elementCoarse), coarseSatW);
        Scalar resSatN = problem_.spatialParams().materialLawParams(elementCoarse).snr();

        Scalar reconstSaturation[numPhases];
        if (veModel == sharpInterface) //reconstruct phase saturation for sharp interface ve model
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
        else if (veModel == capillaryFringe) //reconstruct phase saturation for capillary fringe model
        {
            Scalar densityW = density_[wPhaseIdx];
            Scalar densityNw = density_[nPhaseIdx];
            Scalar lambda = problem_.spatialParams().materialLawParams(elementCoarse).lambda();
            Scalar entryP = problem_.spatialParams().materialLawParams(elementCoarse).pe();

            reconstSaturation[wPhaseIdx] = 1.0;
            reconstSaturation[nPhaseIdx] = 0.0;
            if(height > minGasPlumeDist && height < gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = 1.0-resSatN;
                reconstSaturation[nPhaseIdx] = resSatN;
            }
            else if(height >= gasPlumeDist)
            {
                reconstSaturation[wPhaseIdx] = std::pow(((height - gasPlumeDist) * (densityW - densityNw) * problem_.gravity().two_norm() + entryP), (-lambda))
                * std::pow(entryP, lambda) * (1.0 - resSatW - resSatN) + resSatW;
                reconstSaturation[nPhaseIdx] = 1.0 - reconstSaturation[wPhaseIdx];
            }
        }
        else //reconstruct phase saturation for no ve model
        {
            reconstSaturation[wPhaseIdx] = coarseSatW;
            reconstSaturation[nPhaseIdx] = 1.0-coarseSatW;
        }
        return reconstSaturation[phaseIdx];
    }

    /*! \brief Calculates saturation integral over height given two bounds. Attention! For saturation you have to divide by height!
     */
    Scalar saturationIntegral(Scalar lowerBound, Scalar upperBound, const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar gasPlumeDist = cellData.gasPlumeDist();
        Scalar minGasPlumeDist = cellData.minGasPlumeDist();
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(this->residualSegSaturation(element), cellData.saturation(wPhaseIdx));
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
            if(veModel == sharpInterface)
            {
                integral += (upperBound - lowerPartBound)*resSatW;
            }

            if(veModel == capillaryFringe)//numerical integration
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

    /*! \brief Calculates wetting saturation integral over height given two bounds and a container holding coarse cell information.
     * Attention! For saturation you have to divide by height!
     */
    Scalar saturationIntegral(Scalar lowerBound, Scalar upperBound, const AdaptedValues& valuesCoarse, const Element& elementCoarse)
    {
        int veModel = valuesCoarse.veModel;
        Scalar gasPlumeDist = valuesCoarse.gasPlumeDist;
        Scalar minGasPlumeDist = valuesCoarse.minGasPlumeDist;
//        std::cout << "gasPlumeDist " << gasPlumeDist << "minGasPlumeDist " << minGasPlumeDist << std::endl;
        Scalar coarseSatW = valuesCoarse.saturationW;
        //set residual saturations according to some criterium
        Scalar resSatW = std::min(this->residualSegSaturation(elementCoarse), coarseSatW);
        Scalar resSatN = problem_.spatialParams().materialLawParams(elementCoarse).snr();

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
            if(veModel == sharpInterface)
            {
                integral += (upperBound - lowerPartBound)*resSatW;
            }

            if(veModel == capillaryFringe)//numerical integration
            {
                int intervalNumber = 10;
                Scalar deltaZ = (upperBound - lowerPartBound)/intervalNumber;
                for(int count=0; count<intervalNumber; count++ )
                {
                    integral += (reconstSaturation(lowerPartBound + count*deltaZ, wPhaseIdx, valuesCoarse, elementCoarse) +
                            reconstSaturation(lowerPartBound + (count+1)*deltaZ, wPhaseIdx, valuesCoarse, elementCoarse))/2.0*deltaZ;
                }
            }
        }
//        std::cout << "coarseSatW " << coarseSatW << "integral " << integral << std::endl;
        return integral;
    }

    /*! \brief Calculates relative permeability as numerical integral of fine-scale permeabilities over z,
     *
     */
    Scalar calculateRelPermeabilityCoarse(Scalar lowerBound, Scalar upperBound, const Element& element, int phaseIdx)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();
        Scalar gasPlumeDist = cellData.gasPlumeDist();

        Scalar permeabilityCoarse = 0.0;

        //Attention: This introduces a hard-coded lense
//        if(lowerBound < gasPlumeDist && phaseIdx == wPhaseIdx)
//        {
//            Scalar upperPartBound = std::min(gasPlumeDist, upperBound);
//            permeabilityCoarse += upperPartBound - lowerBound;
//        }
//        if(upperBound > gasPlumeDist)
//        {
            Scalar lowerPartBound = lowerBound;
            if(veModel == sharpInterface)
            {
                //set residual saturations according to some criterium
                Scalar resSatW = std::min(this->residualSegSaturation(element), cellData.saturation(wPhaseIdx));
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

            if(veModel == capillaryFringe)//numerical integration
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
//        }

        return (permeabilityCoarse / (upperBound - lowerBound));
    }

    //! \brief Write data files
     /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        FVPressure2P<TypeTag>::addOutputVtkFields(writer);

        if(problem_.vtkOutputLevel()>=1)
            // add multidim stuff
        {
            int size = problem_.gridView().size(0);
            ScalarSolutionType *subdomain = writer.allocateManagedBuffer(size);
            for (int i = 0; i < size; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);
                (*subdomain)[i] = cellData.veModel();
            }
            writer.attachCellData(*subdomain, "veModel");
        }
    }

    const Scalar density(int phaseIdx) const
    {
        return density_[phaseIdx];
    }

    const Scalar gravity() const
    {
        return gravity_.two_norm();
    }

    //! Constructs a FVPressure2P object
    /**
     * \param problem A problem class object
     */
    FVPressure2PVEMultiDim(Problem& problem) : ParentType(problem), problem_(problem), velocity_(problem), gravity_(problem.gravity())
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
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant
    FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity)> velocity_;

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    //! gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
    //! gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
};

/*!\brief Function which assembles the system of equations to be solved
 *
 *  This function assembles the Matrix and the right hand side (RHS) vector to solve for
 * a pressure field with a Finite-Volume (FV) discretization.
 * Implementations of this base class have to provide the methods <tt>getSource()</tt>,
 * <tt>getStorage()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt> if the assemble() method is called!
 *
 * \param first Indicates if function is called at the initialization step or during the simulation
 *              (If <tt>first</tt> is <tt>true</tt>, no pressure field of previous iterations is required)
 */
template<class TypeTag>
void FVPressure2PVEMultiDim<TypeTag>::assemble(bool first)
{
//    if(first)
//    {
//        ParentType::assemble(true);
//        return;
//    }

    // initialization: set matrix A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        // get the global index of the cell
        int eIdxGlobalI = problem_.variables().index(element);

        // assemble interior element contributions
        if (element.partitionType() == Dune::InteriorEntity)
        {
            // get the cell data
            CellData& cellDataI = problem_.variables().cellData(eIdxGlobalI);

            EntryType entries(0.);

            /*****  source term ***********/
            this->getSource(entries, element, cellDataI, first);
            this->f_[eIdxGlobalI] += entries[rhs];

            /*****  flux term ***********/
            // iterate over all faces of the cell
            for (const auto& intersection : Dune::intersections(problem_.gridView(), element))
            {
                /************* handle interior face *****************/
                if (intersection.neighbor())
                {
                    auto elementNeighbor = intersection.outside();

                    int eIdxGlobalJ = problem_.variables().index(elementNeighbor);

                    // check for hanging nodes
                    // take a hanging node never from the element with smaller level!
                    bool haveSameLevel = (element.level() == elementNeighbor.level());
                    // calculate only from one side, but add matrix entries for both sides
                    // the last condition is needed to properly assemble in the presence
                    // of ghost elements
                    if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce)
                        && (eIdxGlobalI > eIdxGlobalJ)
                        && elementNeighbor.partitionType() == Dune::InteriorEntity)
                        continue;

                    //check for hanging nodes
                    entries = 0;
                    this->getFlux(entries, intersection, cellDataI, first);

                    //set right hand side
                    this->f_[eIdxGlobalI] -= entries[rhs];

                    // set diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];

                    // set off-diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalJ] -= entries[matrix];

                    // The second condition is needed to not spoil the ghost element entries
                    if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce)
                        && elementNeighbor.partitionType() == Dune::InteriorEntity)
                    {
                        this->f_[eIdxGlobalJ] += entries[rhs];
                        this->A_[eIdxGlobalJ][eIdxGlobalJ] += entries[matrix];
                        this->A_[eIdxGlobalJ][eIdxGlobalI] -= entries[matrix];
                    }

                } // end neighbor

                /************* boundary face ************************/
                else
                {
                    entries = 0;
                    this->getFluxOnBoundary(entries, intersection, cellDataI, first);

                    //set right hand side
                    this->f_[eIdxGlobalI] += entries[rhs];
                    // set diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
                }
            } //end interfaces loop
    //        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);

            /*****  storage term ***********/
            entries = 0;
            this->getStorage(entries, element, cellDataI, first);
            this->f_[eIdxGlobalI] += entries[rhs];
    //         set diagonal entry
            this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
        }
        // assemble overlap and ghost element contributions
        else
        {
            this->A_[eIdxGlobalI] = 0.0;
            this->A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            this->f_[eIdxGlobalI] = this->pressure()[eIdxGlobalI];
        }
    } // end grid traversal
//    printmatrix(std::cout, A_, "global stiffness matrix after assempling", "row", 11,3);
//    printvector(std::cout, f_, "right hand side", "row", 10);
}

/*! \brief Function which calculates the flux entry
 *
 * \copydetails FVPressure::getFlux(EntryType&,const Intersection&,const CellData&,const bool)
 *
 */
template<class TypeTag>
void FVPressure2PVEMultiDim<TypeTag>::getFlux(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(elementJ));

    // get global coordinates of cell centers
    GlobalPosition globalPosI = elementI.geometry().center();
    GlobalPosition globalPosJ = elementJ.geometry().center();

    // get mobilities and fractional flow factors
    Scalar lambdaWI = cellData.mobility(wPhaseIdx);
    Scalar lambdaNwI = cellData.mobility(nPhaseIdx);
    Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
    Scalar lambdaNwJ = cellDataJ.mobility(nPhaseIdx);

    int veModelI = cellData.veModel();
    int veModelJ = cellDataJ.veModel();

    //get viscosities
    Scalar viscosityW = viscosity_[wPhaseIdx];
    Scalar viscosityNw = viscosity_[nPhaseIdx];

    if(veModelI != veModelJ)
    {
        CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        Scalar refinement2D = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, 2DRefinement);
        Scalar deltaZ = problem_.bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, refinement2D));
        if(veModelI < 2)
        {
            Scalar top = globalPosJ[dim - 1] + deltaZ/2.0;
            Scalar bottom = globalPosJ[dim - 1] - deltaZ/2.0;
            lambdaWI = calculateRelPermeabilityCoarse(bottom, top, elementI, wPhaseIdx) / viscosityW;
            lambdaNwI = calculateRelPermeabilityCoarse(bottom, top, elementI, nPhaseIdx) / viscosityNw;
        }
        if(veModelJ < 2)
        {
            Scalar top = globalPosI[dim - 1] + deltaZ/2.0;
            Scalar bottom = globalPosI[dim - 1] - deltaZ/2.0;
            lambdaWJ = calculateRelPermeabilityCoarse(bottom, top, elementJ, wPhaseIdx) / viscosityW;
            lambdaNwJ = calculateRelPermeabilityCoarse(bottom, top, elementJ, nPhaseIdx) / viscosityNw;
        }
    }

    // get capillary pressure
    Scalar pcI = cellData.capillaryPressure();
    Scalar pcJ = cellDataJ.capillaryPressure();

    if(veModelI != veModelJ)
    {
        if(veModelI < 2)
        {
            pcI = problem_.pressureModel().reconstCapillaryPressure(globalPosJ[dim-1],elementI);
        }
        if(veModelJ < 2)
        {
            pcJ = problem_.pressureModel().reconstCapillaryPressure(globalPosI[dim-1],elementJ);
        }
    }

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    // get face area
    Scalar faceArea = intersection.geometry().volume();

    // distance vector between barycenters
    //TODO: only valid for horizontal domain, make more flexible
    if(veModelI != veModelJ)
    {
        if(veModelI < 2)
        {
            globalPosI[dim-1] = globalPosJ[dim-1];
        }
        if(veModelJ < 2)
        {
            globalPosJ[dim-1] = globalPosI[dim-1];
        }
    }
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    // compute vectorized permeabilities
    DimMatrix meanPermeability(0);

    problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(elementI),
            problem_.spatialParams().intrinsicPermeability(elementJ));

    Dune::FieldVector<Scalar, dim> permeability(0);
    meanPermeability.mv(unitOuterNormal, permeability);

    Scalar rhoMeanW = 0;
    Scalar rhoMeanNw = 0;
    if (compressibility_)
    {
        rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
        rhoMeanNw = 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));
    }

    //calculate potential gradients
    Scalar potentialDiffW = 0;
    Scalar potentialDiffNw = 0;

    //if we are at the very first iteration we can't calculate phase potentials
    if (!first)
    {
        Scalar potentialIW = cellData.potential(wPhaseIdx);
        Scalar potentialINw = cellData.potential(nPhaseIdx);
        Scalar potentialJW = cellDataJ.potential(wPhaseIdx);
        Scalar potentialJNw = cellDataJ.potential(nPhaseIdx);

        //this is necessary because potentials cannot be stored for ghost cells
        if(veModelI != veModelJ)
        {
            if(veModelI < 2)
            {
                potentialIW = reconstPressure(globalPosJ[dim-1],wPhaseIdx,elementI) + (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[wPhaseIdx];
                potentialINw = reconstPressure(globalPosJ[dim-1],nPhaseIdx,elementI)+ (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[nPhaseIdx];
            }
            if(veModelJ < 2)
            {
                potentialJW = reconstPressure(globalPosI[dim-1],wPhaseIdx,elementJ) + (problem_.bBoxMax() - globalPosI) * gravity_ * density_[wPhaseIdx];
                potentialJNw = reconstPressure(globalPosI[dim-1],nPhaseIdx,elementJ)+ (problem_.bBoxMax() - globalPosI) * gravity_ * density_[nPhaseIdx];
            }
        }

        potentialDiffW = potentialIW - potentialJW;
        potentialDiffNw = potentialINw - potentialJNw;

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

            density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
            density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];

            potentialDiffW = cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx);
            potentialDiffNw = cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx);

            potentialDiffW += density_[wPhaseIdx] * (distVec * gravity_);
            potentialDiffNw += density_[nPhaseIdx] * (distVec * gravity_);
        }
    }

    //do the upwinding of the mobility depending on the phase potentials
    Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWJ;
    lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
    Scalar lambdaNw = (potentialDiffNw > 0) ? lambdaNwI : lambdaNwJ;
    lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNw;

    if (compressibility_)
    {
        density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
        density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

        density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
        density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];
    }

    Scalar scalarPerm = permeability.two_norm();

    //calculate current matrix entry
    entry[matrix] = (lambdaW + lambdaNw) * scalarPerm / dist * faceArea;

    //calculate right hand side
    //calculate unit distVec
    distVec /= dist;
    Scalar areaScaling = (unitOuterNormal * distVec);
    //this treatment of g allows to account for gravity flux through faces where the face normal has no z component (e.g. parallelepiped grids)
    entry[rhs] = (lambdaW * density_[wPhaseIdx] + lambdaNw * density_[nPhaseIdx]) * scalarPerm * (gravity_ * distVec) * faceArea * areaScaling;
    //add term for VE-cell at boundary: unknown pressure has to be written as reconstructed pressure in terms of height
    if(veModelI != veModelJ)
    {
        if(veModelI < 2 && globalPosJ[dim-1] <= cellData.gasPlumeDist())
        {
            entry[rhs] -= (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * density_[wPhaseIdx] * gravity_.two_norm() * globalPosJ[dim-1];
        }
        else if(veModelI == 0 && globalPosJ[dim-1] > cellData.gasPlumeDist())
        {
            entry[rhs] -= (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * (density_[wPhaseIdx] * gravity_.two_norm() * cellData.gasPlumeDist()
                    + density_[nPhaseIdx] * gravity_.two_norm() * (globalPosJ[dim-1]-cellData.gasPlumeDist()));
        }
        else if(veModelI == 1 && globalPosJ[dim-1] > cellData.gasPlumeDist())
        {
            entry[rhs] -= (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * density_[wPhaseIdx] * gravity_.two_norm() * globalPosJ[dim-1];
        }
        if(veModelJ < 2 && globalPosI[dim-1] <= cellDataJ.gasPlumeDist())
        {
            entry[rhs] += (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * density_[wPhaseIdx] * gravity_.two_norm() * globalPosI[dim-1];
        }
        else if(veModelJ == 0 && globalPosI[dim-1] > cellDataJ.gasPlumeDist())
        {
            entry[rhs] += (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * (density_[wPhaseIdx] * gravity_.two_norm() * cellDataJ.gasPlumeDist()
                    + density_[nPhaseIdx] * gravity_.two_norm() * (globalPosI[dim-1]-cellDataJ.gasPlumeDist()));
        }
        else if(veModelJ == 1 && globalPosI[dim-1] > cellDataJ.gasPlumeDist())
        {
            entry[rhs] += (lambdaW + lambdaNw) * scalarPerm / dist * faceArea * density_[wPhaseIdx] * gravity_.two_norm() * globalPosI[dim-1];
        }
    }

    if (pressureType_ == pw)
    {
        //add capillary pressure term to right hand side
        entry[rhs] += 0.5 * (lambdaNwI + lambdaNwJ) * scalarPerm * (pcI - pcJ) / dist * faceArea;
    }
    else if (pressureType_ == pn)
    {
        //add capillary pressure term to right hand side
        entry[rhs] -= 0.5 * (lambdaWI + lambdaWJ) * scalarPerm * (pcI - pcJ) / dist * faceArea;
    }
}

/*! \brief Function which calculates the flux entry at a boundary
 *
 * \copydetails FVPressure::getFluxOnBoundary(EntryType&,const Intersection&,const CellData&,const bool)
 *
 * Dirichlet boundary condition for pressure depends on the formulation (\f$p_w\f$ (default), \f$p_n\f$, \f$p_{global}\f$),
 * Neumann boundary condition are the phase mass fluxes (\f$q_w\f$ and \f$q_n\f$, [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVPressure2PVEMultiDim<TypeTag>::getFluxOnBoundary(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first)
{
    Element element = intersection.inside();
    int veModel = cellData.veModel();

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
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

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

        problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(element));

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
        if (veModel == sharpInterface)
        {
            //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
            //sharp interface injection cells
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            pcBound = gravity_.two_norm() * problem_.bBoxMax()[dim - 1] * (satW - resSatW) / (1.0 - resSatW) * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            //minGasPlume = init
        }
        else if (veModel == capillaryFringe)
        {
            //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
            //sharp interface injection cells (with entry pressure)
            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
            Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();
            pcBound = gravity_.two_norm() * problem_.bBoxMax()[dim - 1] * (satW - resSatW) / (1.0 - resSatW) * (density_[nPhaseIdx] - density_[wPhaseIdx]) + entryP;
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
            viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx) / densityWBound;
            viscosityNwBound = FluidSystem::viscosity(fluidState, nPhaseIdx) / densityNwBound;

            rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + densityWBound);
            rhoMeanNw = 0.5 * (cellData.density(nPhaseIdx) + densityNwBound);
        }

        Scalar lambdaWBound = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), satW) / viscosityWBound;
        Scalar lambdaNwBound = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element), satW) / viscosityNwBound;
        //since there is no possibility to store values like minGasPlumeDist or gasPlumeDist for boundary ghost cells they are treated as
        //sharp interface injection cells
        if (veModel == sharpInterface || veModel == capillaryFringe)
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
void FVPressure2PVEMultiDim<TypeTag>::updateMaterialLaws()
{
    for (const auto& element : Dune::elements(problem_.gridView()))
    {
        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        int veModel = cellData.veModel();

        Scalar temperature = problem_.temperature(element);

        //determine phase saturations from primary saturation variable
        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar satNw = cellData.saturation(nPhaseIdx);

        Scalar densityW = density_[wPhaseIdx];
        Scalar densityNw = density_[nPhaseIdx];

        Scalar viscosityW = viscosity_[wPhaseIdx];
        Scalar viscosityNw = viscosity_[nPhaseIdx];

        //set pseudo residual segSaturation inside plume
        this->calculateResidualSegSaturation(element);
        //set gas plume distance
        this->calculateGasPlumeDist(element, satW);
        Scalar gasPlumeDist = cellData.gasPlumeDist();

        Scalar pc = MaterialLaw::pc(problem_.spatialParams().materialLawParams(element), satW);
        if (veModel == sharpInterface)
        {
            pc = gravity_.two_norm() * gasPlumeDist * (densityNw - densityW);
        }
        else if (veModel == capillaryFringe)
        {
            Scalar entryP = problem_.spatialParams().materialLawParams(element).pe();
            pc = gravity_.two_norm() * gasPlumeDist * (densityNw - densityW) + entryP; //TODO: compressibility (for the Dumux model, not necessary in Bo's slightly compressible model)
        }
//        std::cout << "veModel " << veModel << "satW " << satW << "pc " << pc << std::endl;

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
        if (veModel == sharpInterface || veModel == capillaryFringe)
        {
            mobilityW = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, wPhaseIdx) / viscosityW;
            mobilityNw = calculateRelPermeabilityCoarse(0.0, problem_.bBoxMax()[dim - 1], element, nPhaseIdx) / viscosityNw;  //TODO: compressibility (like in Bo's model)
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
        if (veModel == sharpInterface || veModel == capillaryFringe)
        {
            globalPos[dim-1] = 0.0;
        }
        Scalar gravityDiff = (problem_.bBoxMax() - globalPos) * gravity_;

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
