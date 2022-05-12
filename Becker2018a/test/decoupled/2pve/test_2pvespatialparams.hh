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
 * \brief spatial parameters for the sequential 2p test
 */
#ifndef TEST_2PVE_SPATIALPARAMS_HH
#define TEST_2PVE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/exponentialmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscoreyLambda.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/io/plotmateriallaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class Test2PVESpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test2PVESpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(Test2PVESpatialParams, SpatialParams, Dumux::Test2PVESpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test2PVESpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef BrooksCorey<Scalar> RawMaterialLaw;
//    typedef LinearMaterial<Scalar> RawMaterialLaw;
//    typedef ExponentialMaterial<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/*!
 *
 * \ingroup IMPETtests
 * \brief spatial parameters for the sequential 2p test
 */
template<class TypeTag>
class Test2PVESpatialParams: public FVSpatialParams<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        dimWorld=Grid::dimensionworld,
        dim = Grid::dimension,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief This is called from the problem and creates a gnuplot output
     *        of e.g the pc-Sw curve
     */
    void plotMaterialLaw()
    {
        PlotMaterialLaw<TypeTag> plotMaterialLaw;
        GnuplotInterface<Scalar> gnuplot(plotFluidMatrixInteractions_);
        gnuplot.setOpenPlotWindow(plotFluidMatrixInteractions_);
        plotMaterialLaw.addpcswcurve(gnuplot, materialLawParams_, 0.0, 1.0, "");
        plotMaterialLaw.addkrcurves(gnuplot, materialLawParams_, 0.0, 1.0, "");
    }

    const FieldMatrix& intrinsicPermeability(const Element& element) const
    {
        GlobalPosition globalPos = element.geometry().center();
        if (globalPos[0] > 100 && globalPos[0] < 120)
        {
            return permeabilityLense_;
        }
        return permeability_; //(20mD)
    }

    double porosity(const Element& element) const
    {
        return porosity_;
    }


    // return the parameter object for the Brooks-Corey material law which depends on the position
//    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    const MaterialLawParams& materialLawParams(const Element& element) const
    {
            return materialLawParams_;
    }

    //set the residual saturations for the phases
    void setResidualSaturation(Scalar residualSat[numPhases])
    {
        materialLawParams_.setSwr(residualSat[wPhaseIdx]);
        materialLawParams_.setSnr(residualSat[nPhaseIdx]);
    }


    Test2PVESpatialParams(const GridView& gridView)
    : ParentType(gridView), permeability_(0), permeabilityLense_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

        plotFluidMatrixInteractions_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions);
        Scalar lambdaPc = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Lambda);
        Scalar lambdaKr = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LambdaKr);
        Scalar entryPressure = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, EntryPressure);
        Scalar maxPc = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MaxPc);
        Scalar exponent = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Exponent);
        Scalar permeabilityH = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, PermeabilityHorizontal);
        Scalar permeabilityV = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, PermeabilityVertical);
        Scalar porosity = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Porosity);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(entryPressure);
        // Brooks-Corey shape parameter for relative permeability
        materialLawParams_.setLambdaPc(lambdaPc);
        // Brooks-Corey shape parameter for relative permeability
        materialLawParams_.setLambdaKr(lambdaKr);

        // parameters for the linear law
//        materialLawParams_.setEntryPc(entryPressure);
//        materialLawParams_.setMaxPc(maxPc);

        // parameters for the exponential law
//        materialLawParams_.setPe(entryPressure);
//        materialLawParams_.setLambda(lambdaPc);
//        materialLawParams_.setExponent(exponent);

        permeability_[0][0] = permeabilityH;
        permeability_[1][1] = permeabilityV;

        permeabilityLense_[0][0] = 1.334e-12;
        permeabilityLense_[1][1] = 1.334e-12;

        porosity_ = porosity;
    }

private:
    MaterialLawParams materialLawParams_;
    bool plotFluidMatrixInteractions_;
    FieldMatrix permeabilityLense_;
    FieldMatrix permeability_;
    Scalar porosity_;
};

} // end namespace
#endif
