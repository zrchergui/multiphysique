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
#ifndef DUMUX_FVSATURATION2PVE_HH
#define DUMUX_FVSATURATION2PVE_HH

//#include <dumux/decoupled/2p/transport/transportproperties2p.hh>
#include "dumux/porousmediumflow/2p/sequential/transport/cellcentered/saturation.hh"

/**
 * \file
 * \brief  Finite Volume discretization of a saturation transport equation
 */

namespace Dumux
{
//! \ingroup FVSaturation2p
//! \brief The finite volume discretization of a saturation transport equation
/*! This model solves equations of the form
 *
 *  \f[
 *  \phi \frac{\partial (\rho_\alpha S_\alpha)}{\partial t} + \text{div}\, (\rho_\alpha \boldsymbol{v_\alpha}) = q_\alpha,
 *  \f]
 *
 *  where \f$ S_\alpha \f$ is the saturation of phase \f$ \alpha \f$ (wetting \f$(w) \f$,
 *  non-wetting \f$(n) \f$) and \f$ \boldsymbol v_\alpha \f$ is the phase velocity defined by
 *  the multi-phase Darcy equation.
 *  If a phase velocity is reconstructed from the pressure solution it can be directly inserted into
 *  the previous equation. In the incompressible case the equation is further divided by the phase density
 *  \f$ \rho_\alpha \f$. If a total velocity is reconstructed the saturation equation is reformulated into:
 *
 * \f[
 *  \phi \frac{\partial S_w}{\partial t} + f_w \text{div}\, \boldsymbol{v}_{t} + f_w \lambda_n \boldsymbol{K}\left(\textbf{grad}\,
 *  p_c + (\rho_n-\rho_w) \, g \, \textbf{grad} z \right)= q_\alpha,
 * \f]
 * to get a wetting phase saturation or
 * \f[
 * \phi \frac{\partial S_n}{\partial t} + f_n \text{div}\, \boldsymbol{v}_{t} - f_n \lambda_w \boldsymbol{K}\left(\textbf{grad}\,
 * p_c + (\rho_n-\rho_w) \, g \, \textbf{grad} z \right)= q_\alpha,
 * \f]
 * if the non-wetting phase saturation is the primary transport variable.
 *
 *  The total velocity formulation is only implemented for incompressible fluids and \f$ f_\alpha \f$
 *  is the fractional flow function, \f$ \lambda_\alpha \f$ is the mobility, \f$ \boldsymbol K \f$
 *  the absolute permeability,\f$ p_c \f$ the capillary pressure, \f$ \rho \f$ the fluid density,
 *  \f$ g \f$ the gravity constant, and \f$ q \f$ the source term.
 *
 *
 *  In the IMPES models the default setting is:
 *
 * formulation: \f$ p_w \f$ - \f$ S_w \f$ (Property: \a Formulation defined as \a DecoupledTwoPCommonIndices::pwsw)
 *
 * compressibility: disabled (Property: \a EnableCompressibility set to \a false)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVSaturation2PVE: public FVSaturation2P<TypeTag>
{
    typedef FVSaturation2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TransportModel) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Velocity) Velocity;
    typedef typename GET_PROP_TYPE(TypeTag, CapillaryFlux) CapillaryFlux;
    typedef typename GET_PROP_TYPE(TypeTag, GravityFlux) GravityFlux;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        vt = Indices::velocityTotal,
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        saturationIdx = Indices::saturationIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:

    /*! \brief Globally updates the saturation solution
     *
     * \param updateVec Vector containing the global update.
     */
    void updateSaturationSolution(TransportSolutionType& updateVec)
    {
        Scalar dt = problem_.timeManager().timeStepSize();
        ElementIterator eEndIt = problem_.gridView().template end<0>();

        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt; ++eIt)
        {
            int eIdxGlobal = problem_.variables().index(*eIt);
            updateSaturationSolution(*eIt, updateVec[eIdxGlobal][0], dt);
        }
    }

    /*! \brief Globally updates the saturation solution
     *
     * \param updateVec Vector containing the global update.
     * \param dt time step for update
     */
    void updateSaturationSolution(TransportSolutionType& updateVec, Scalar dt)
    {
        ElementIterator eEndIt = problem_.gridView().template end<0>();

        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eEndIt; ++eIt)
        {
            int eIdxGlobal = problem_.variables().index(*eIt);
            updateSaturationSolution(*eIt, updateVec[eIdxGlobal][0], dt);
        }
    }

    /*! \brief Updates the saturation solution of a cell
     *
     * Calculates secondary saturation variables and stores saturations.
     * Calculation and subsequently storage of gasPlumeDist is necessary for the
     * adaptive scheme. It needs to be up to date if a VE cell is converted
     * into a full-D cell.
     *
     * \param eIdxGlobal Global cell index
     * \param update Cell saturation update
     * \param dt Current time step
     */
    void updateSaturationSolution(const Element& element, Scalar update, Scalar dt)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        switch (saturationType_)
        {
        case sw:
        {
            cellData.setOldSaturation(cellData.saturation(wPhaseIdx));
            Scalar sat = cellData.saturation(wPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx, sat);
                cellData.setSaturation(nPhaseIdx, 1 - sat);
                if(veModel_ == 0 || veModel_ == 1)
                    problem_.pressureModel().calculateGasPlumeDist(element, sat);
            break;
        }
        case sn:
        {
            Scalar sat = cellData.saturation(nPhaseIdx) + dt*update;

                cellData.setSaturation(wPhaseIdx,1 -sat);
                cellData.setSaturation(nPhaseIdx, sat);
                if(veModel_ == 0 || veModel_ == 1)
                    problem_.pressureModel().calculateGasPlumeDist(element, 1-sat);
            break;
        }
        }
    }

    /*! \brief Constructs a FVSaturation2P object
     *
     * \param problem A problem class object
     */
    FVSaturation2PVE(Problem& problem) :
            ParentType(problem), problem_(problem), threshold_(1e-6),
            switchNormals_(GET_PARAM_FROM_GROUP(TypeTag, bool, Impet, SwitchNormals))
    {
        if (compressibility_ && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        capillaryFlux_ = std::make_shared<CapillaryFlux>(problem);
        gravityFlux_ = std::make_shared<GravityFlux>(problem);
        velocity_ = std::make_shared<Velocity>(problem);

        vtkOutputLevel_ = GET_PARAM_FROM_GROUP(TypeTag, int, Vtk, OutputLevel);
        porosityThreshold_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, PorosityThreshold);
        veModel_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Problem& problem_;
    std::shared_ptr<Velocity> velocity_;
    std::shared_ptr<CapillaryFlux> capillaryFlux_;
    std::shared_ptr<GravityFlux> gravityFlux_;

    int vtkOutputLevel_;
    Scalar porosityThreshold_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);

    Scalar density_[numPhases];

    const Scalar threshold_;
    const bool switchNormals_;

    int veModel_;
};



}
#endif
