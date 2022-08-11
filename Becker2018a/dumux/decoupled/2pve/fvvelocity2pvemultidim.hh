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
#ifndef DUMUX_FVVELOCITY2PVEMULTIDIM_HH
#define DUMUX_FVVELOCITY2PVEMULTIDIM_HH

/**
 * @file
 * @brief  Velocity Field from a finite volume solution of a pressure equation.
 */

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/gridenums.hh>
//#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/velocity.hh>

namespace Dumux
{
//! \ingroup FVPressure2p
//! \brief Determines the velocity from a finite volume solution of the  pressure equation of a sequential model (IMPES).
/*! Calculates phase velocities or total velocity from a known pressure field applying a finite volume discretization.
 * The wetting or the non-wetting phase pressure, or the global pressure has to be given as piecewise constant cell values.
 * The phase velocities are calculated following  Darcy's law as
 * \f[
 * \boldsymbol v_\alpha = \lambda_\alpha \boldsymbol K \left(\textbf{grad}\, p_\alpha + \rho_\alpha g  \textbf{grad}\, z \right),
 * \f]
 * where \f$ p_\alpha \f$ denotes the pressure of phase \f$ _\alpha \f$ (wetting or non-wetting),
 * \f$ \boldsymbol K \f$ the absolute permeability, \f$ \lambda_\alpha \f$ the phase mobility,
 * \f$ \rho_\alpha \f$ the phase density and \f$ g \f$ the gravity constant.
 * The total velocity is either calculated as sum of the phase velocities
 * \f[
 * \boldsymbol v_{total} = \boldsymbol v_{wetting}+\boldsymbol v_{non-wetting},
 * \f]
 * or with a given global pressure
 * \f[
 * \boldsymbol v_{total} = \lambda_{total} \boldsymbol K \left(\textbf{grad}\,
 *  p_{global} + \sum f_\alpha \rho_\alpha g  \textbf{grad}\, z\right).
 * \f]
 *
 * \tparam TypeTag The Type Tag
 */

template<class TypeTag>
class FVVelocity2PVEMultiDim : public FVVelocity2P<TypeTag>
{
    typedef FVVelocity2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView)GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        vt = Indices::velocityTotal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef std::array<unsigned int, dim> CellArray;

public:
    /*! \brief Constructs a FVVelocity2P object
     * \param problem A Problem class object
     */
    FVVelocity2PVEMultiDim(Problem& problem) :
        ParentType(problem), problem_(problem), gravity_(problem.gravity())
    {
        if (GET_PROP_VALUE(TypeTag, EnableCompressibility) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;
    }

    //! For initialization
    void initialize()
    {
        ParentType::initialize();
        if (!compressibility_)
        {
            const auto element = *problem_.gridView().template begin<0> ();
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
    }

    // Calculates the velocity at a cell-cell interface.
    void calculateVelocity(const Intersection&, CellData&);

    /*! \brief Indicates if velocity is reconstructed in the pressure step or in the transport step
     *
     * Returns false (In the adaptive finite volume scheme the velocity has to be calculated separately
     * to make sure the hanging nodes are treated correctly.)
     */
//    bool calculateVelocityInTransport()
//    {
//            return false;
//    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    //! gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    //! gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
};

/*! \brief Calculates the velocity at a cell-cell interface.
 *
 * Calculates the velocity at a cell-cell interface from a given pressure field.
 *
 * \param intersection Intersection of two grid cells
 * \param cellData Object containing all model relevant cell data
 */
template<class TypeTag>
void FVVelocity2PVEMultiDim<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellData)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    int eIdxGlobalJ = problem_.variables().index(elementJ);

    CellData& cellDataJ = problem_.variables().cellData(eIdxGlobalJ);

    // get global coordinates of cell centers
    GlobalPosition globalPosI = (elementI).geometry().center();
    GlobalPosition globalPosJ = (elementJ).geometry().center();

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
            lambdaWI = problem_.pressureModel().calculateRelPermeabilityCoarse(bottom, top, elementI, wPhaseIdx) / viscosityW;
            lambdaNwI = problem_.pressureModel().calculateRelPermeabilityCoarse(bottom, top, elementI, nPhaseIdx) / viscosityNw;
        }
        if(veModelJ < 2)
        {
            Scalar top = globalPosI[dim - 1] + deltaZ/2.0;
            Scalar bottom = globalPosI[dim - 1] - deltaZ/2.0;
            lambdaWI = problem_.pressureModel().calculateRelPermeabilityCoarse(bottom, top, elementJ, wPhaseIdx) / viscosityW;
            lambdaNwI = problem_.pressureModel().calculateRelPermeabilityCoarse(bottom, top, elementJ, nPhaseIdx) / viscosityNw;
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

    //get face index
    int isIndexI = intersection.indexInInside();
    int isIndexJ = intersection.indexInOutside();

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    // distance vector between barycenters
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

    //calculate potential gradients
    Scalar potentialIW = cellData.potential(wPhaseIdx);
    Scalar potentialINw = cellData.potential(nPhaseIdx);
    Scalar potentialJW = cellDataJ.potential(wPhaseIdx);
    Scalar potentialJNw = cellDataJ.potential(nPhaseIdx);

    //this is necessary because potentials cannot be stored for ghost cells
    //for the velocity calculation after the pressure calculation we have to also regard for the fact that the pressure has been updated but not the potential
    if(veModelI != veModelJ)
    {
        if(veModelI < 2)
        {
//            potentialIW = problem_.pressureModel().reconstPressure(cellData.oldPressure(wPhaseIdx),globalPosJ[dim-1],wPhaseIdx,elementI)
//                        + (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[wPhaseIdx];
//            potentialINw = problem_.pressureModel().reconstPressure(cellData.oldPressure(nPhaseIdx),globalPosJ[dim-1],nPhaseIdx,elementI)
//                         + (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[nPhaseIdx];
            potentialIW = problem_.pressureModel().reconstPressure(cellData.pressure(wPhaseIdx),globalPosJ[dim-1],wPhaseIdx,elementI)
                        + (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[wPhaseIdx];
            potentialINw = problem_.pressureModel().reconstPressure(cellData.pressure(nPhaseIdx),globalPosJ[dim-1],nPhaseIdx,elementI)
                         + (problem_.bBoxMax() - globalPosJ) * gravity_ * density_[nPhaseIdx];
        }
        if(veModelJ < 2)
        {
//            potentialJW = problem_.pressureModel().reconstPressure(cellData.oldPressure(wPhaseIdx),globalPosI[dim-1],wPhaseIdx,elementI)
//                                + (problem_.bBoxMax() - globalPosI) * gravity_ * density_[wPhaseIdx];
//            potentialJNw = problem_.pressureModel().reconstPressure(cellData.oldPressure(nPhaseIdx),globalPosI[dim-1],nPhaseIdx,elementI)
//                                + (problem_.bBoxMax() - globalPosI) * gravity_ * density_[nPhaseIdx];
            potentialJW = problem_.pressureModel().reconstPressure(cellData.pressure(wPhaseIdx),globalPosI[dim-1],wPhaseIdx,elementI)
                                + (problem_.bBoxMax() - globalPosI) * gravity_ * density_[wPhaseIdx];
            potentialJNw = problem_.pressureModel().reconstPressure(cellData.pressure(nPhaseIdx),globalPosI[dim-1],nPhaseIdx,elementI)
                                + (problem_.bBoxMax() - globalPosI) * gravity_ * density_[nPhaseIdx];
        }
    }

    Scalar potentialDiffW = potentialIW - potentialJW;
    Scalar potentialDiffNw = potentialINw - potentialJNw;

    if (compressibility_)
    {
        density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
        density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

        density_[wPhaseIdx] =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                        density_[wPhaseIdx];
        density_[nPhaseIdx] =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                        density_[nPhaseIdx];

        potentialDiffW = (cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx));
        potentialDiffNw = (cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx));

        potentialDiffW += density_[wPhaseIdx] * (distVec * gravity_); //delta z/delta x in unitOuterNormal[z]
        potentialDiffNw += density_[nPhaseIdx] * (distVec * gravity_);
    }

    //store potentials for further calculations (velocity, saturation, ...)
    cellData.fluxData().setUpwindPotential(wPhaseIdx, isIndexI, potentialDiffW);
    cellData.fluxData().setUpwindPotential(nPhaseIdx, isIndexI, potentialDiffNw);

    cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, isIndexJ, -potentialDiffW);
    cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, isIndexJ, -potentialDiffNw);

    //do the upwinding of the mobility depending on the phase potentials
    Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWJ;
    lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
    Scalar lambdaNw = (potentialDiffNw > 0.) ? lambdaNwI : lambdaNwJ;
    lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNw;

    if (compressibility_)
    {
        density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
        density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

        density_[wPhaseIdx] =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                        density_[wPhaseIdx];
        density_[nPhaseIdx] =
                (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                        density_[nPhaseIdx];
    }

    Scalar scalarPerm = permeability.two_norm();

    //calculate the gravity term
    Dune::FieldVector<Scalar, dimWorld> velocityW(unitOuterNormal);
    Dune::FieldVector<Scalar, dimWorld> velocityNw(unitOuterNormal);

    //calculate unit distVec
    distVec /= dist;
    Scalar areaScaling = (unitOuterNormal * distVec);
    //this treatment of g allows to account for gravity flux through faces where the face normal
    //has no z component (e.g. parallelepiped grids)
    Scalar gravityTermW = (gravity_ * distVec) * density_[wPhaseIdx] * areaScaling;
    Scalar gravityTermNw = (gravity_ * distVec) * density_[nPhaseIdx] * areaScaling;

    //calculate velocity depending on the pressure used -> use pc = pn - pw
    switch (pressureType_)
    {
    case pw:
    {
        Scalar pressureIW = cellData.pressure(wPhaseIdx);
        Scalar pressureJW = cellDataJ.pressure(wPhaseIdx);
        if(veModelI != veModelJ)
        {
            if(veModelI < 2)
            {
                pressureIW = problem_.pressureModel().reconstPressure(globalPosJ[dim-1], wPhaseIdx, elementI);
            }
            if(veModelJ < 2)
            {
                pressureJW = problem_.pressureModel().reconstPressure(globalPosI[dim-1], wPhaseIdx, elementJ);
            }
        }
        velocityW *= lambdaW * scalarPerm
                * ((pressureIW - pressureJW) / dist + gravityTermW);
        velocityNw *= lambdaNw * scalarPerm
                * ((pressureIW - pressureJW) / dist + gravityTermNw)
                + 0.5 * (lambdaNwI + lambdaNwJ) * scalarPerm * (pcI - pcJ) / dist;
        break;
    }
    case pn:
    {
        velocityW *= lambdaW * scalarPerm
                * ((cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist + gravityTermW)
                - 0.5 * (lambdaWI + lambdaWJ) * scalarPerm * (pcI - pcJ) / dist;
        velocityNw *= lambdaNw * scalarPerm
                * ((cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist + gravityTermNw);
        break;
    }
    case pGlobal:
    {
        velocityW *= (lambdaW + lambdaNw) * scalarPerm * (cellData.globalPressure() - cellDataJ.globalPressure()) / dist
                + scalarPerm * (lambdaW * gravityTermW + lambdaNw * gravityTermNw);
        velocityNw = 0;
        break;
    }
    }

    //store velocities
    //TODO: store mean velocity for hanging nodes, right now only velocity of last visited intersection is stored (makes no difference for calculation, only for output)
    if(veModelI != veModelJ)
    {
        if(veModelI < 2)
        {
            cellDataJ.fluxData().setVelocity(wPhaseIdx, isIndexJ, velocityW);
            cellDataJ.fluxData().setVelocity(nPhaseIdx, isIndexJ, velocityNw);
            cellDataJ.fluxData().setVelocityMarker(isIndexJ);

            cellData.fluxData().setVelocity(wPhaseIdx, isIndexI, velocityW);
            cellData.fluxData().setVelocity(nPhaseIdx, isIndexI, velocityNw);
        }
        if(veModelJ < 2)
        {
            cellData.fluxData().setVelocity(wPhaseIdx, isIndexI, velocityW);
            cellData.fluxData().setVelocity(nPhaseIdx, isIndexI, velocityNw);
            cellData.fluxData().setVelocityMarker(isIndexI);

            cellDataJ.fluxData().setVelocity(wPhaseIdx, isIndexJ, velocityW);
            cellDataJ.fluxData().setVelocity(nPhaseIdx, isIndexJ, velocityNw);
        }
    }
    else
    {
    cellData.fluxData().setVelocity(wPhaseIdx, isIndexI, velocityW);
    cellData.fluxData().setVelocity(nPhaseIdx, isIndexI, velocityNw);
    cellData.fluxData().setVelocityMarker(isIndexI);

    cellDataJ.fluxData().setVelocity(wPhaseIdx, isIndexJ, velocityW);
    cellDataJ.fluxData().setVelocity(nPhaseIdx, isIndexJ, velocityNw);
    cellDataJ.fluxData().setVelocityMarker(isIndexJ);
    }

//                        printvector(std::cout, cellData.fluxData().velocity(), "velocity", "row", 4, 1, 3);
}
}
#endif
