/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef TEST_ENERGY_STORAGE_SPATIALPARAMS_HH
#define TEST_ENERGY_STORAGE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/io/plotmateriallaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class TestEnergyStorageSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG (TestEnergyStorageSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(TestEnergyStorageSpatialParams, SpatialParams, Dumux::TestEnergyStorageSpatialParams<TypeTag>);

// Set the material law
SET_PROP(TestEnergyStorageSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
//    typedef LinearMaterial<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};
}

/** \todo Please doc me! */

template<class TypeTag>
class TestEnergyStorageSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld, numEq = 1
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

    Scalar intrinsicPermeability (const Element& element) const
    {
        GlobalPosition globalPos = element.geometry().center();
        if (globalPos[0] > 100 && globalPos[0] < 120 && globalPos[dim-1] > 20)
	  return 2.0e-12; //2.0e-15
        return permeability_; //(20mD)
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const Element &element) const
    {
            return materialLawParams_;
    }

    TestEnergyStorageSpatialParams(const GridView& gridView) :
            ParentType(gridView), permeability_(0)
    {
        plotFluidMatrixInteractions_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions);

        // residual saturations
        materialLawParams_.setSwr(0.0);
        materialLawParams_.setSnr(0.0);

//        // parameters for the Brooks-Corey Law
//        // entry pressures
        materialLawParams_.setPe(1e5);
//        // Brooks-Corey shape parameters
        materialLawParams_.setLambda(2.0);

        // parameters for the linear
        // entry pressures function
//        materialLawParams_.setEntryPc(0);
//        materialLawParams_.setMaxPc(0);

        permeability_ = 2.0e-12;
    }

private:
    MaterialLawParams materialLawParams_;
    bool plotFluidMatrixInteractions_;
    Scalar permeability_;
};

} // end namespace
#endif
