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
 * \brief test problem for the sequential 2p model
 */
#ifndef DUMUX_TEST_TWOPVE_PROBLEM_HH
#define DUMUX_TEST_TWOPVE_PROBLEM_HH

#include <dumux/io/gnuplotinterface.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <appl/VEmod/h2o.hh>
#include <appl/VEmod/ch4.hh>
#include <dumux/material/components/simpleco2.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/decoupled/2pve/celldata2pve.hh>
#include <dumux/decoupled/2pve/fvsaturation2pve.hh>
#include <dumux/decoupled/2pve/fvpressure2pve.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

#include "test_2pvespatialparams.hh"

#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxdefault.hh>

namespace Dumux
{

template<class TypeTag>
class TWOPVETestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TWOPVETestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, Test2PVESpatialParams));

// Set the grid type
SET_TYPE_PROP(TWOPVETestProblem, Grid, Dune::UGGrid<2>);

// Set the problem property
SET_TYPE_PROP(TWOPVETestProblem, Problem, Dumux::TWOPVETestProblem<TypeTag>);

// Set the cell data
SET_TYPE_PROP(TWOPVETestProblem, CellData, CellData2PVE<TypeTag>);

//Set the transport model
SET_TYPE_PROP(TWOPVETestProblem, TransportModel, FVSaturation2PVE<TypeTag>);

//Set the pressure model
SET_TYPE_PROP(TWOPVETestProblem, PressureModel, FVPressure2PVE<TypeTag>);

////////////////////////////////////////////////////////////////////////
//Switch to a p_n-S_w formulation
//
//SET_INT_PROP(IMPESTestProblem, Formulation,
//        DecoupledTwoPCommonIndices::pnsn);
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//Switch to a p_global-S_w formulation
//
//SET_INT_PROP(IMPESTestProblem, Formulation,
//        DecoupledTwoPCommonIndices::pGlobalSw);
//
//Define the capillary pressure term in the transport equation -> only needed in case of a p_global-S_w formulation!
//SET_TYPE_PROP(IMPESTestProblem, CapillaryFlux, CapillaryDiffusion<TypeTag>);
//
//Define the gravity term in the transport equation -> only needed in case of a p_global-S_w formulation!
//SET_TYPE_PROP(IMPESTestProblem, GravityFlux, GravityPart<TypeTag>);
//
////////////////////////////////////////////////////////////////////////

// Set the wetting phase
SET_PROP(TWOPVETestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TWOPVETestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::GasPhase<Scalar, Dumux::CH4<Scalar> > type;
//    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

SET_TYPE_PROP(TWOPVETestProblem, GridView, Dune::UGGrid<2>::LevelGridView);

// Enable gravity
SET_BOOL_PROP(TWOPVETestProblem, ProblemEnableGravity, true);

//SET_TYPE_PROP(TWOPVETestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
SET_TYPE_PROP(TWOPVETestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxDefault<TypeTag>);
SET_BOOL_PROP(TWOPVETestProblem, EnableCompressibility, false);
}
/*!
 * \ingroup 2PVETest
 *
 * \brief test problem for the sequential 2p model
 *
 * Water is injected from the left side into a rectangular 2D domain also
 * filled with water. Upper and lower boundary is closed (Neumann = 0),
 * and there is free outflow on the right side.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_impes -parameterFile ./test_impes.input</tt>,
 * where the arguments define the parameter file..
 */
template<class TypeTag>
class TWOPVETestProblem: public IMPESProblem2P<TypeTag>
{
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename Grid::LeafGridView LeafGridView;
typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonWettingPhase;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    nPhaseIdx = Indices::nPhaseIdx,
    wPhaseIdx = Indices::wPhaseIdx,
    pwIdx = Indices::pwIdx,
    swIdx = Indices::swIdx,
    eqIdxPress = Indices::pressureEqIdx,
    eqIdxSat = Indices::satEqIdx
};

enum VEModel
{
    sharpInterface,
    capillaryFringe
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename LeafGridView::Traits::template Codim<0>::Entity LeafElement;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
typedef typename SpatialParams::MaterialLaw MaterialLaw;

typedef std::array<unsigned int, dim> CellArray;

public:
TWOPVETestProblem(TimeManager &timeManager, const GridView &gridView) :
ParentType(timeManager, gridView), eps_(1e-6), visualizationWriter_(GridCreator::grid().leafGridView(), "gridAfterReconstruction"),
normWriter_(GridCreator::grid().leafGridView(), "gridAfterRefinement")
{
    int outputInterval = 1e9;
    if (ParameterTree::tree().hasKey("Problem.OutputInterval"))
    {
        outputInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputInterval);
    }
    this->setOutputInterval(outputInterval);

    if (ParameterTree::tree().hasKey("Problem.OutputTimeInterval"))
    {
        Scalar outputTimeInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputTimeInterval);
        this->setOutputTimeInterval(outputTimeInterval);
    }

    const int level = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, Reconstruction);

    // store pointer to all ve-columns in a map, indicated by their global index
    // iterate over all elements
    for (const auto& element : Dune::elements(this->gridView()))
    {
        // column number equals global Index
        int eIdxGlobal = this->variables().index(element);

        mapColumns_.insert(std::make_pair(eIdxGlobal, element));
        dummy_ = element;
    }

    //calculate length of capillary transition zone (CTZ)
    Scalar swr = this->spatialParams().materialLawParams(dummy_).swr();
    Scalar snr = this->spatialParams().materialLawParams(dummy_).snr();
    Scalar satW1 = 1.0 - snr;
    Scalar satW2 = swr + 0.1*(1.0-swr-snr);
    Scalar pc1 = MaterialLaw::pc(this->spatialParams().materialLawParams(dummy_), satW1);
    Scalar pc2 = MaterialLaw::pc(this->spatialParams().materialLawParams(dummy_), satW2);
    GlobalPosition globalPos = dummy_.geometry().center();
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar tempRef = temperatureAtPos(globalPos);
    Scalar densityW = WettingPhase::density(tempRef, pRef);
    Scalar densityN = NonWettingPhase::density(tempRef, pRef);
    Scalar gravity = this->gravity().two_norm();
    CTZ_ = (pc2-pc1)/((densityW-densityN)*gravity);
    std::cout << "CTZ " << CTZ_ << std::endl;

    LeafGridView visualizationGridView = GridCreator::grid().leafGridView();
    for(int i=0;i<level;i++)
    {
        for (const auto& element : Dune::elements(visualizationGridView))
        {
            GridCreator::grid().mark(element, UG::D2::BLUE, 1);
        }
        // adapt the grid
        GridCreator::grid().preAdapt();
        GridCreator::grid().adapt();
        GridCreator::grid().postAdapt();
    }

    residualSegSaturation_ = 1.0;
    residualSegSatAverage_ = 1.0;
    corrGasPlumeHeight_ = 0.0;
    corrGasPlumeHeightOld_ = 0.0;
    gasPlumeTipOld_ = 0.0;
    averageSat_ = 0.0;
    averageSatOld_ = 0.0;
    isSegregated_ = false;

    name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);

    outputFile_.open("plumeTip.out", std::ios::trunc);
    outputFile_.close();
}

/*!
 * \brief Called by the Dumux::TimeManager in order to
 *        initialize the problem.
 */
void init()
{
    ParentType::init();

    bool plotFluidMatrixInteractions = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions);

    // plot the Pc-Sw curves, if requested
    if(plotFluidMatrixInteractions)
        this->spatialParams().plotMaterialLaw();
}

void postTimeStep()
{
    ParentType::postTimeStep();

    int outputInterval = 1e9;
    if (ParameterTree::tree().hasKey("Problem.OutputInterval"))
    {
        outputInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputInterval);
    }

    CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
    int reconstruction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, Reconstruction);
    double deltaX = this->bBoxMax()[0]/numberOfCells[0];
    double deltaZ = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, reconstruction));

    //write output for reconstructed solution
    if(this->timeManager().timeStepIndex() % outputInterval == 0 || this->timeManager().willBeFinished()
            || this->timeManager().episodeWillBeFinished() || this->timeManager().timeStepIndex() == 0)
    {
        LeafGridView visualizationGridView = GridCreator::grid().leafGridView();
        // use visualization writer
        visualizationWriter_.gridChanged();
        visualizationWriter_.beginWrite(this->timeManager().time() + this->timeManager().timeStepSize());
        //write stuff out
        typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;

        int size = visualizationGridView.size(0);
        ScalarSolutionType *pressureW = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *pressureNw = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *capillaryPressure = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *saturationW = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *saturationNw = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *mobilityW = visualizationWriter_.allocateManagedBuffer (size);
        ScalarSolutionType *mobilityNw = visualizationWriter_.allocateManagedBuffer (size);

        for (const auto& element : Dune::elements(visualizationGridView))
        {
            GlobalPosition globalPos = element.geometry().center();
            int eIdxGlobal = this->variables().index(element);
            int j = round((globalPos[0] - (deltaX/2.0))/deltaX);

            Element veElement = mapColumns_.find(j)->second;

            Scalar top = globalPos[dim - 1] + deltaZ/2.0;
            Scalar bottom = globalPos[dim - 1] - deltaZ/2.0;

            Scalar satW = this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
            Scalar satNw = 1.0-satW;
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar temp = temperatureAtPos(globalPos);

            (*pressureW)[eIdxGlobal] = this->pressureModel().reconstPressure(globalPos[dim-1], wPhaseIdx, veElement);
            (*pressureNw)[eIdxGlobal] = this->pressureModel().reconstPressure(globalPos[dim-1], nPhaseIdx, veElement);
            (*capillaryPressure)[eIdxGlobal] = this->pressureModel().reconstCapillaryPressure(globalPos[dim-1], veElement);
            (*saturationW)[eIdxGlobal] = satW;
            (*saturationNw)[eIdxGlobal] = satNw;
            (*mobilityW)[eIdxGlobal] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, wPhaseIdx)/WettingPhase::viscosity(temp, pRef);
            (*mobilityNw)[eIdxGlobal] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, nPhaseIdx)/NonWettingPhase::viscosity(temp, pRef);
        }
        visualizationWriter_.attachCellData(*pressureW, "wetting pressure");
        visualizationWriter_.attachCellData(*pressureNw, "non-wetting pressure");
        visualizationWriter_.attachCellData(*capillaryPressure, "capillary pressure");
        visualizationWriter_.attachCellData(*saturationW, "wetting saturation");
        visualizationWriter_.attachCellData(*saturationNw, "non-wetting saturation");
        visualizationWriter_.attachCellData(*mobilityW, "wetting mobility");
        visualizationWriter_.attachCellData(*mobilityNw, "non-wetting mobility");
        visualizationWriter_.endWrite();

//        //refine grid for error norm calculation
//        GridCreator::grid().preAdapt();
//        const int additionalRefinementSteps = 2;
//        for(int i=0;i<additionalRefinementSteps;i++)
//        {
//            for (const auto& element : Dune::elements(visualizationGridView))
//            {
//                GlobalPosition globalPos = element.geometry().center();
//                GridCreator::grid().mark(1, element);
//            }
//            // adapt the grid
//            GridCreator::grid().adapt();
//            GridCreator::grid().postAdapt();
//        }
//        int reconstruction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, Reconstruction);
//        const double deltaZRefined = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, (reconstruction+additionalRefinementSteps)));

//        // use error norm writer
//        normWriter_.gridChanged();
//        normWriter_.beginWrite(this->timeManager().time() + this->timeManager().timeStepSize());
//
//        size = visualizationGridView.size(0);
//        ScalarSolutionType *pressureW1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *pressureNw1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *capillaryPressure1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *saturationW1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *saturationNw1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *mobilityW1 = normWriter_.allocateManagedBuffer (size);
//        ScalarSolutionType *mobilityNw1 = normWriter_.allocateManagedBuffer (size);
//
//        //write output on refined grid
//        for (const auto& element : Dune::elements(visualizationGridView))
//        {
//            //identify column number of element and check if it is a VE column or a full-d column
//            GlobalPosition globalPos = element.geometry().center();
//            int columnNumber = floor(globalPos[0]/deltaX);//starting with 0
//
//            //identify ID of cell, same as ID of parallel run with 20 blocks
//            //determine block number starting from 0
//            int blockNumberX = floor(globalPos[0]/35);
//            int blockNumberY = floor(globalPos[dim-1]/15);
//            int blockNumber = blockNumberY*10 + blockNumberX;
//            //determine local cell ID in block, starting from 0
//            Scalar localX = globalPos[0] - (blockNumberX*35);
//            Scalar localY = globalPos[dim-1] - (blockNumberY*15);
//            int localNumberX = floor(localX/0.25);
//            int localNumberY = floor(localY/(30.0/512.0));
//            int localID = localNumberY*140 + localNumberX;
//            //determine global ID, starting from 0
//            int globalID = blockNumber*35840 + localID;
//
//            int j = round((globalPos[0] - (deltaX/2.0))/deltaX);
//            Element veElement = mapColumns_.find(j)->second;
//
//            const double deltaZRefined = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, (reconstruction+additionalRefinementSteps)));
//            Scalar top = globalPos[dim - 1] + deltaZRefined/2.0;
//            Scalar bottom = globalPos[dim - 1] - deltaZRefined/2.0;
//
//            Scalar satW = this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
//            Scalar satNw = 1.0-satW;
//            Scalar pRef = referencePressureAtPos(globalPos);
//            Scalar temp = temperatureAtPos(globalPos);
//
//            (*pressureW1)[globalID] = this->pressureModel().reconstPressure(globalPos[dim-1], wPhaseIdx, veElement);
//            (*pressureNw1)[globalID] = this->pressureModel().reconstPressure(globalPos[dim-1], nPhaseIdx, veElement);
//            (*capillaryPressure1)[globalID] = this->pressureModel().reconstCapillaryPressure(globalPos[dim-1], veElement);
//            (*saturationW1)[globalID] = satW;
//            (*saturationNw1)[globalID] = satNw;
//            (*mobilityW1)[globalID] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, wPhaseIdx)/WettingPhase::viscosity(temp, pRef);
//            (*mobilityNw1)[globalID] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, nPhaseIdx)/NonWettingPhase::viscosity(temp, pRef);
//        }
//        normWriter_.attachCellData(*pressureW1, "wetting pressure");
//        normWriter_.attachCellData(*pressureNw1, "non-wetting pressure");
//        normWriter_.attachCellData(*capillaryPressure1, "capillary pressure");
//        normWriter_.attachCellData(*saturationW1, "wetting saturation");
//        normWriter_.attachCellData(*saturationNw1, "non-wetting saturation");
//        normWriter_.attachCellData(*mobilityW1, "wetting mobility");
//        normWriter_.attachCellData(*mobilityNw1, "non-wetting mobility");
//        normWriter_.endWrite();
//
//        //coarsen grid after error norm calculation
//        for(int i=0;i<additionalRefinementSteps;i++)
//        {
//            for (const auto& element : Dune::elements(visualizationGridView))
//            {
//                GlobalPosition globalPos = element.geometry().center();
//                GridCreator::grid().mark(-1, element);
//            }
//            // adapt the grid
//            GridCreator::grid().adapt();
//            GridCreator::grid().postAdapt();
//            GridCreator::grid().preAdapt();
//        }
    }

    //check mass conservativity, calculate gas plume distance and average height for capillary fringe:
    GridView GridView = this->gridView();
    Scalar totalMassN = 0;
    Scalar gasPlumeTip = 0.0;
    Scalar corrGasPlumeHeight = 0.0;
    Scalar averageSat = 0.0;
    Scalar averageSatOldLocation = 0.0;
    Scalar volumePlumeColumns = 0.0;
    Scalar volumePlumeColumnsOldLocation = 0.0;
    for (const auto& element : Dune::elements(GridView))
    {
        int eIdxGlobal = this->variables().index(element);
        CellData& cellData = this->variables().cellData(eIdxGlobal);
        GlobalPosition globalPos = element.geometry().center();
        Scalar pRef = referencePressureAtPos(globalPos);
        Scalar temp = temperatureAtPos(globalPos);
        Scalar volume = element.geometry().volume();
        Scalar massN = cellData.saturation(nPhaseIdx) * volume * this->spatialParams().porosity(element) * NonWettingPhase::density(temp, pRef);
        totalMassN += massN;
        int j = round((globalPos[0] - (deltaX/2.0))/deltaX);
        //calculate saturation value depending on reconstruction
        Scalar top = this->bBoxMax()[dim - 1];
        Scalar bottom = top - deltaZ;
        Scalar satWTop = this->pressureModel().saturationIntegral(bottom, top, element)/(top-bottom);
        if(satWTop < (1.0-eps_))
        {
            gasPlumeTip = std::max(gasPlumeTip, (j+1)*deltaX);
            Scalar satW = cellData.saturation(wPhaseIdx);
            averageSat += satW * volume;
            volumePlumeColumns += volume;

            if(globalPos[0] < gasPlumeTipOld_)
            {
                averageSatOldLocation += satW * volume;
                volumePlumeColumnsOldLocation += volume;
            }

            Scalar gasPlumeDist = cellData.gasPlumeDist();
            if(gasPlumeDist + CTZ_ < 0.0)//all of column is gas plume not including capillary fringe
            {
                corrGasPlumeHeight += top * volume;
            }
            else if((top - (gasPlumeDist + CTZ_)) > 0.0)//otherwise all is capillary fringe
            {
                corrGasPlumeHeight += (top - (gasPlumeDist + CTZ_)) * volume;
            }
        }
    }
    averageSat_ = averageSat/volumePlumeColumns;
    averageSatOldLocation = averageSatOldLocation/volumePlumeColumnsOldLocation;
    corrGasPlumeHeight = corrGasPlumeHeight/volumePlumeColumns;
    if(isnan(averageSat_))
    {
        averageSat_ = 0.0;
//        std::cout << "attention: averageSat_ is nan " << "averageSat " << averageSat << "volumePlumeColumns " << volumePlumeColumns << std::endl;
    }
    if(isnan(averageSatOldLocation))
        averageSatOldLocation = 0.0;
    if(isnan(corrGasPlumeHeight))
        corrGasPlumeHeight = 0.0;
    calculateResidualSegSaturation(gasPlumeTip, corrGasPlumeHeight, averageSatOldLocation, averageSat_);
    gasPlumeTipOld_ = gasPlumeTip;

    corrGasPlumeHeight = 0.0;
    Scalar volumePlumeColumnsCorrected = 0.0;
    Scalar residualSegSatAverage = 0.0;
    for (const auto& element : Dune::elements(GridView))
    {
        int eIdxGlobal = this->variables().index(element);
        CellData& cellData = this->variables().cellData(eIdxGlobal);
        Scalar satW = cellData.saturation(wPhaseIdx);
        this->pressureModel().calculateGasPlumeDist(element, satW);
        GlobalPosition globalPos = element.geometry().center();
        Scalar volume = element.geometry().volume();
        //calculate saturation value depending on reconstruction
        Scalar top = this->bBoxMax()[dim - 1];
        Scalar bottom = top - deltaZ;
        Scalar satWTop = this->pressureModel().saturationIntegral(bottom, top, element)/(top-bottom);
        if(satWTop < (1.0-eps_))
        {
            Scalar gasPlumeDist = cellData.gasPlumeDist();
            if(gasPlumeDist + CTZ_ < 0.0)//all of column is gas plume not including capillary fringe
            {
                corrGasPlumeHeight += top * volume;
                residualSegSatAverage += cellData.residualSegSaturation() * volume;
                volumePlumeColumnsCorrected += volume;
            }
            else if((top - (gasPlumeDist + CTZ_)) > 0.0)//otherwise all is capillary fringe
            {
                corrGasPlumeHeight += (top - (gasPlumeDist + CTZ_)) * volume;
                residualSegSatAverage += cellData.residualSegSaturation() * volume;
                volumePlumeColumnsCorrected += volume;
            }
        }
    }
    corrGasPlumeHeightOld_ = corrGasPlumeHeight/volumePlumeColumns;
    residualSegSatAverage_ = residualSegSatAverage/volumePlumeColumnsCorrected;
    if(isnan(residualSegSatAverage_))
        residualSegSatAverage_ = 0.0;
//    std::cout << "corrGasPlumeHeightOld_: " << corrGasPlumeHeightOld_ << std::endl;

    Scalar time = this->timeManager().time()+this->timeManager().timeStepSize();
    Scalar totalMassNExpected = -GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, Injectionrate)
            * time * this->bBoxMax()[1];
    std::cout << "Error in mass: " << totalMassNExpected - totalMassN << ". ";

    if(this->timeManager().timeStepIndex() % outputInterval == 0 || this->timeManager().willBeFinished()
            || this->timeManager().episodeWillBeFinished() || this->timeManager().timeStepIndex() == 0)
    {
        //write out distance of gas plume from injection
        outputFile_.open("plumeTip.out", std::ios::app);
        outputFile_ << time << " " << gasPlumeTip << std::endl;
        outputFile_.close();
    }
}

/*! \brief Calculates new residual saturation for VE-model depending on segregation time
 */
void calculateResidualSegSaturation(Scalar gasPlumeTip, Scalar corrGasPlumeHeight, Scalar averageSatOldLocation, Scalar averageSat)
{
    int correction = 0;
    if (ParameterTree::tree().hasKey("VE.correction"))
    {
        correction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, correction);
    }

    Scalar segSatW;
    if(correction == 1 || correction == 2 || correction == 3 || correction == 4)
    {
        GlobalPosition globalPos = dummy_.geometry().center();
        Scalar pRef = referencePressureAtPos(globalPos);
        Scalar tempRef = temperatureAtPos(globalPos);
        Scalar aquiferHeight = this->bBoxMax()[dim - 1];
        Scalar charHeight = aquiferHeight;
        Scalar porosity = this->spatialParams().porosity(dummy_);
        Scalar viscosityW = WettingPhase::viscosity(tempRef, pRef);
        Scalar viscosityNw = NonWettingPhase::viscosity(tempRef, pRef);
        Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        Scalar permeability = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, PermeabilityVertical);
        Scalar densityW = WettingPhase::density(tempRef, pRef);
        Scalar densityN = NonWettingPhase::density(tempRef, pRef);
        Scalar gravity = this->gravity().two_norm();
        Scalar resSatW = this->spatialParams().materialLawParams(dummy_).swr();
        Scalar resSatN = this->spatialParams().materialLawParams(dummy_).snr();
        Scalar exponent = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, Exponent);
        Scalar lambda = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, LambdaKr);
        int relPermModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, Model);
        Scalar entryP = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, EntryPressure);
        Scalar injectionRate = -GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, Scalar, "BoundaryConditions", Injectionrate);

        if(correction == 2)//correct the domain height with analytic solution
        {
            charHeight = viscosityNw/viscosityW * aquiferHeight;
        }
        else if(correction == 3)//correct the domain height with numerical solution
        {
            Scalar gasVolume = injectionRate * aquiferHeight
                              * time / (densityN*porosity);
            charHeight = gasVolume/gasPlumeTip;
        }
        if(isnan(charHeight) || isinf(charHeight))
            charHeight = 0.0;

        if(correction == 4)//calculate saturation explicitly
        {
            Scalar timeStepSize = this->timeManager().timeStepSize();
            Scalar mobilityW = MaterialLaw::krw(this->spatialParams().materialLawParams(dummy_), residualSegSaturation_)/viscosityW;
            Scalar mobilityN = MaterialLaw::krn(this->spatialParams().materialLawParams(dummy_), residualSegSaturation_)/viscosityNw;
//            Scalar drainageFlux = mobilityW*permeability*gravity*(densityW-densityN);
            Scalar drainageFlux = mobilityW*mobilityN/(mobilityW+mobilityN)*permeability*gravity*(densityW-densityN);

            //calculate new residualSegSaturation inside plume with same corrGasPlumeHeight as before (corrGasPlumeHeightOld_) due to vertical flow
            //gas flowing into gasPlume region replaces water out of plume
//            std::cout << "residualSegSaturation vorher " << residualSegSaturation_ << std::endl;

            //neglect tip that flows out of plume (it automatically gets new residualSegSaturation), 1. bug (works perfectly), 2. correct solution (does not work)
            Scalar satDiff = averageSat - averageSatOld_;
//            Scalar satDiff = averageSatOldLocation - averageSatOld_;
//            std::cout << "satDiff " << satDiff << " averageSat " << averageSat << " averageSatOld_ " << averageSatOld_ <<  std::endl;
//            std::cout << "gasPlumeTipOld_ " << gasPlumeTipOld_ << std::endl;
//            std::cout << "gasPlumeTip " << gasPlumeTip << std::endl;
//            if(gasPlumeTipOld_ <= gasPlumeTip + 1e-6 && gasPlumeTipOld_ >= gasPlumeTip - 1e-6)
//            {
                residualSegSaturation_ = residualSegSaturation_ + satDiff*aquiferHeight/corrGasPlumeHeightOld_;
//            }
//            else
//            {
//                residualSegSaturation_ = (residualSegSaturation_ * corrGasPlumeHeightOld_ * porosity * gasPlumeTipOld_ +
//                        (gasPlumeTip - gasPlumeTipOld_) * porosity * corrGasPlumeHeightOld_ - injectionRate/densityN * timeStepSize * aquiferHeight)/
//                        (gasPlumeTip*corrGasPlumeHeightOld_*porosity);
//            }

            //entire new plume region gets new residualSegSaturation depending on corrGasPlumeHeightOld_
//            residualSegSaturation_ = 1.0 - (1.0-residualSegSaturation_) * gasPlumeTipOld_/gasPlumeTip -
//                    (1.0-averageSat)*aquiferHeight/corrGasPlumeHeightOld_ + (1.0-averageSatOld_)*aquiferHeight/corrGasPlumeHeightOld_*gasPlumeTipOld_/gasPlumeTip;
            averageSatOld_ = averageSat;
            if(isinf(residualSegSaturation_) || isnan(residualSegSaturation_))
                residualSegSaturation_ = 1.0;
//            std::cout << "residualSegSaturation nachher1 " << residualSegSaturation_ << std::endl;

            if(residualSegSaturation_ < eps_)
            {
//                std::cout << "attention: time step too large1!!" << std::endl;
                residualSegSaturation_ = 0.0;
                isSegregated_ = true;
            }

            if(this->timeManager().timeStepIndex() < 10)//no segregation in first 10 time steps
            {
                residualSegSaturation_ = 1.0;
                corrGasPlumeHeight_ = aquiferHeight;
            }
            else if(this->timeManager().timeStepIndex() == 10)//start of segregation after 10th time step
            {
                residualSegSaturation_ = averageSat;
                corrGasPlumeHeight_ = aquiferHeight;
                isSegregated_ = false;
            }
            else if(drainageFlux > corrGasPlumeHeightOld_ * residualSegSaturation_ * (1.0-residualSegSaturation_) * porosity/timeStepSize)
            {
//                std::cout << "attention: time step too large2!!" << std::endl;
                residualSegSaturation_ = 0.0;
                isSegregated_ = true;
            }
            else
            {
                //calculate and store new residualSegSaturation
                Scalar oldResidualSegSaturation = residualSegSaturation_;
                residualSegSaturation_ = (corrGasPlumeHeightOld_ * porosity * residualSegSaturation_ - drainageFlux  * timeStepSize/(1.0-residualSegSaturation_)) /
                        (corrGasPlumeHeightOld_ * porosity - drainageFlux * timeStepSize / (1.0-residualSegSaturation_));
//                std::cout << "residualSegSaturation nachher2 " << residualSegSaturation_ << std::endl;
                //calculate and store new corrGasPlumeHeight
                corrGasPlumeHeight_ = corrGasPlumeHeightOld_ - drainageFlux * timeStepSize/(porosity * (1.0-oldResidualSegSaturation));
            }

            if(isSegregated_ == true)
            {
                residualSegSaturation_ = 0.0;
            }

//            corrGasPlumeHeightOld_ = corrGasPlumeHeight;
        }
        else if(relPermModel == 0 || (relPermModel == 1 && exponent == 1))
        {
            residualSegSaturation_ = (charHeight*porosity*viscosityW/(time*permeability*(densityW-densityN)*gravity))*(1.0-resSatW-resSatN) + resSatW;
        }
        else
        {
            if(relPermModel == 2)
            {
                exponent = 2.0/lambda + 3.0;
            }
            residualSegSaturation_ = std::pow((charHeight*porosity*viscosityW/(time*permeability*(densityW-densityN)*gravity)), (1.0/(exponent)))*
                    (1.0-resSatW-resSatN) + resSatW;
        }
    }
    else
    {
        residualSegSaturation_ = this->spatialParams().materialLawParams(dummy_).swr();
        Scalar PseudoResidualSaturation = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, VE, PseudoResidualSaturation);
        residualSegSaturation_ = residualSegSaturation_ + PseudoResidualSaturation;
    }
}

/*!
 * \brief The problem name.
 *
 * This is used as a prefix for files generated by the simulation.
 */
const std::string name() const
{
    return name_;
}

bool shouldWriteRestartFile() const
{
    return false;
}

/*!
 * \brief Returns the temperature within the domain.
 *
 * This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 326.0; // -> 53°C
}

// \}

//! Returns the reference pressure for evaluation of constitutive relations
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1.0e7; // -> 100 bar
}

/*!
* \brief Returns the type of boundary condition.
*
* BC for pressure equation can be dirichlet (pressure) or neumann (flux).
*
* BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
*/
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        bcTypes.setAllDirichlet();
    }
    // all other boundaries
    else
    {
        bcTypes.setAllNeumann();
    }
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (globalPos[0] > this->bBoxMax()[0] - eps_)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar temp = temperatureAtPos(globalPos);

//            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
//            if(veModel == sharpInterface || veModel == capillaryFringe)//for ve-models the pressure is defined at the bottom of the domain
//            {
//                values[pwIdx] = (pRef + this->bBoxMax()[dim-1]
//                        * WettingPhase::density(temp, pRef)
//                        * this->gravity().two_norm());
//            }
//            else
//            {
                values[pwIdx] = (pRef + (this->bBoxMax()[dim-1] - globalPos[dim-1])
                        * WettingPhase::density(temp, pRef)
                        * this->gravity().two_norm());
//            }
        }
        else
        {
            values[pwIdx] = 1.0e7;
        }
        values[swIdx] = 1.0;
    }
    else
    {
        values[pwIdx] = 1.0e7;
        values[swIdx] = 0.0;
    }
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0.0;

    if (globalPos[0] < eps_)
    {
        Scalar injectionN = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, Injectionrate);

//        Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();

//      if(time < 2.0e6)
//      {
            values[nPhaseIdx] = injectionN;
//      }
//      else
//      {
//          const Scalar sinus_length = 2.0e6;
//          values[nPhaseIdx] = sin(2.0*M_PI*time/sinus_length) * injectionN + injectionN/3.0;
//          values[wPhaseIdx] = injectionW;
//      }
    }
}

void sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0.0;
}

//! return initial solution -> only saturation values have to be given!
void initial(PrimaryVariables &values,
        const Element& element) const
{
    values[pwIdx] = 1.0e7;
    values[swIdx] = 1.0;
}

Scalar residualSegSaturation() const
{
    return residualSegSaturation_;
}

Scalar residualSegSatAverage() const
{
    return residualSegSatAverage_;
}

Scalar averageSat() const
{
    return averageSat_;
}

Scalar corrGasPlumeHeightAverage() const
{
    return corrGasPlumeHeight_;
}



Scalar capillaryTransitionZone() const
{
    return CTZ_;
}


private:

const Scalar eps_;
std::string name_;
Dumux::VtkMultiWriter<LeafGridView> visualizationWriter_;
Dumux::VtkMultiWriter<LeafGridView> normWriter_;
std::map<int, Element> mapColumns_;
std::ofstream outputFile_;
Element dummy_;
Scalar CTZ_;
Scalar residualSegSaturation_;
Scalar residualSegSatAverage_;
Scalar corrGasPlumeHeight_;
Scalar corrGasPlumeHeightOld_;
Scalar gasPlumeTipOld_;
Scalar averageSatOld_;
Scalar averageSat_;
bool isSegregated_;
Dumux::GnuplotInterface<double> gnuplot_;
};
} //end namespace

#endif
