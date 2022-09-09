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
#ifndef DUMUX_TEST_MULTIDIM_PROBLEM_HH
#define DUMUX_TEST_MULTIDIM_PROBLEM_HH

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/simpleco2.hh>

#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressurepropertiesadaptive.hh>
#include <dumux/decoupled/2pve/celldata2pvemultidim.hh>
#include <dumux/decoupled/2pve/fvvelocity2pvemultidim.hh>
#include <dumux/decoupled/2pve/fvsaturation2pvemultidim.hh>
#include <dumux/decoupled/2pve/fvpressure2pvemultidim.hh>
#include <dumux/decoupled/2pve/gridadaptve.hh>
#include <dumux/decoupled/2pve/gridadaptionindicatorvefulld.hh>
#include <dumux/decoupled/2pve/gridadaptinitializationindicatorvefulld.hh>
#include <dumux/decoupled/2pve/variableclassadaptiveve.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

#include "test_multidimspatialparams.hh"


#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxdefault.hh>


namespace Dumux {

    template<class TypeTag>
    class MultiDimTestProblem;

    //////////
    // Specify the properties
    //////////
    namespace Properties {
        NEW_TYPE_TAG(MultiDimTestProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, TestMultiDimSpatialParams));

        // Set the grid type
        //#if HAVE_UG
        SET_TYPE_PROP(MultiDimTestProblem, Grid, Dune::UGGrid<2>);
        //#else
        //SET_TYPE_PROP(MultiDimTestProblem, Grid, Dune::YaspGrid<2>);
        //#endif

        // Set the problem property
        SET_TYPE_PROP(MultiDimTestProblem, Problem, Dumux::MultiDimTestProblem<TypeTag>);

        // Set the cell data
        SET_TYPE_PROP(MultiDimTestProblem, CellData, CellData2PVEMultiDim<TypeTag>);

        //! Set velocity reconstruction implementation
        SET_TYPE_PROP(MultiDimTestProblem, Velocity, Dumux::FVVelocity2PVEMultiDim<TypeTag>);

        //Set the transport model
        SET_TYPE_PROP(MultiDimTestProblem, TransportModel, FVSaturation2PVEMultiDim<TypeTag>);

        //Set the pressure model
        SET_TYPE_PROP(MultiDimTestProblem, PressureModel, Dumux::FVPressure2PVEMultiDim<TypeTag>);

        //Set the grid adaption model
        SET_TYPE_PROP(MultiDimTestProblem, GridAdaptModel, Dumux::GridAdaptVE<TypeTag, true>);

        //Set the variable class
        SET_TYPE_PROP(MultiDimTestProblem, Variables, Dumux::VariableClassAdaptiveVE<TypeTag>);

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
        SET_PROP(MultiDimTestProblem, WettingPhase)
        {
        private:
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        public:
            typedef FluidSystems::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
        };

        // Set the non-wetting phase
        SET_PROP(MultiDimTestProblem, NonwettingPhase)
        {
        private:
            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        public:
            typedef FluidSystems::GasPhase<Scalar, Dumux::CH4<Scalar> > type;
            //    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
        };

        //SET_TYPE_PROP(MultiDimTestProblem, GridView, Dune::UGGrid<2>::LevelGridView);
        SET_TYPE_PROP(MultiDimTestProblem, LinearSolver, Dumux::SuperLUBackend<TypeTag>);

        // Enable gravity
        SET_BOOL_PROP(MultiDimTestProblem, ProblemEnableGravity, true);

        //SET_TYPE_PROP(MultiDimTestProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);
        SET_BOOL_PROP(MultiDimTestProblem, EnableCompressibility, false);

        //Adaptivity
        SET_TYPE_PROP(MultiDimTestProblem, AdaptionIndicator, GridAdaptionIndicator2PVEFullD<TypeTag>);
        SET_TYPE_PROP(MultiDimTestProblem, AdaptionInitializationIndicator, GridAdaptInitializationIndicatorVEFullD<TypeTag>);
        SET_BOOL_PROP(MultiDimTestProblem, GridAdaptEnableInitializationIndicator, true);
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
    class MultiDimTestProblem: public IMPESProblem2P<TypeTag> {
        typedef IMPESProblem2P<TypeTag> ParentType;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename Grid::GlobalIdSet GlobalIdSet;
        typedef typename GlobalIdSet::IdType IdType;
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
        MultiDimTestProblem(TimeManager &timeManager, const GridView &gridView) :
            ParentType(timeManager, gridView), eps_(1e-6), visualizationWriter_(GridCreator::grid().leafGridView(), "gridAfterReconstruction"),
            debugWriter_(GridCreator::grid().leafGridView(), "debugGridAfterAdaptation"), normWriter_(GridCreator::grid().leafGridView(), "gridAfterRefinement")
        {
            this->setGrid(GridCreator::grid());
            //    GridCreator::grid().globalRefine(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumRefine));

            if (ParameterTree::tree().hasKey("Problem.OutputInterval"))
                {
                    outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputInterval);
                }
            this->setOutputInterval(outputInterval_);

            Scalar outputTimeInterval = 1e6;
            if (ParameterTree::tree().hasKey("Problem.OutputTimeInterval"))
                {
                    outputTimeInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputTimeInterval);
                }
            this->setOutputTimeInterval(outputTimeInterval);

            updateColumnMap();


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
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            CTZ_ = (pc2-pc1)/((densityW-densityNw)*this->gravity().two_norm());
            std::cout << "CTZ " << CTZ_ << std::endl;

            //calculate segregation time
            Scalar height = this->bBoxMax()[dim-1];
            Scalar porosity = this->spatialParams().porosity(dummy_);
            Scalar viscosityW = WettingPhase::viscosity(tempRef, pRef);
            Scalar permeability = this->spatialParams().intrinsicPermeability(dummy_);
            Scalar gravity = this->gravity().two_norm();
            segTime_ = (height*porosity*viscosityW)/(permeability*gravity*(densityW-densityNw));

            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            int numberOfColumns = numberOfCells[0];
            modelVector_.clear();
            if(modelVector_.size() != numberOfColumns)
                {
                    modelVector_.resize(numberOfColumns);
                }
            for(int i=0; i<numberOfColumns; i++)
                {
                    modelVector_[i] = 2;
                }

            name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
            averageSatInPlume_ = 1.0;

            // outputFile_.open("averageSatColumn.out", std::ios::trunc);
            // outputFile_.close();
            // outputFile_.open("errorSat.out", std::ios::trunc);
            // outputFile_.close();
            outputFile_.open("numberOfCells.out", std::ios::trunc);
            outputFile_.close();

            //// create files to store results (Zakaria 2022-08-25 and 2022-08-30)
            outputFile_.open("timeCriterion20.out", std::ios::trunc);
            outputFile_.close();      
            outputFile_.open("timeZp20.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("timeCriterion200.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("timeZp200.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("satProfiles.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("relPermProfiles.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("pressProfiles.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("capPressProfiles.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("criterion.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("zp.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("CPUTime.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("CPUNewton.out", std::ios::trunc);
            outputFile_.close();
            outputFile_.open("iterNumberNewton.out", std::ios::trunc);
            outputFile_.close();
            ////
        }

        void updateColumnMap()// TODO: more efficient way? Only change entries when column has changed
        {
            mapAllColumns_.clear();
            // column number
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            double deltaX = this->bBoxMax()[0]/numberOfCells[0];
            //assign VE model depending on cell volume
            double volumeVE = this->bBoxMax()[1] * deltaX;
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            // store pointer to all ve-elements in a map, indicated by their global index
            // iterate over all elements
            for (const auto& element : Dune::elements(this->gridView()))
                {
                    GlobalPosition globalPos = element.geometry().center();
                    int j = round((globalPos[0] - (deltaX/2.0))/deltaX);
                    mapAllColumns_.insert(std::make_pair(j, element));
                    dummy_ = element;
                }
        }

        void updateModelVector()
        {
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            int numberOfColumns = numberOfCells[0];
            modelVector_.clear();
            if(modelVector_.size() != numberOfColumns)
                {
                    modelVector_.resize(numberOfColumns);
                }
            for(int i=0; i<numberOfColumns; i++)
                {
                    typename std::map<int, Element>::iterator it = mapAllColumns_.find(i);
                    int eIdxGlobal = this->variables().index(it->second);
                    int veModel = this->variables().cellData(eIdxGlobal).veModel();
                    modelVector_[i] = veModel;
                }

        }

        void setModel()// TODO: more efficient way? Only change entries when column has changed
        {
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            double deltaX = this->bBoxMax()[0]/numberOfCells[0];
            double volumeVE = this->bBoxMax()[dim-1] * deltaX;//TODO: 3-D
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            // set Model: VE or full-D
            // iterate over all elements
            for (const auto& element : Dune::elements(this->gridView()))
                {
                    double volume = element.geometry().volume();
                    int eIdxGlobal = this->variables().index(element);
                    if(volume >= volumeVE - eps_)
                        {
                            this->variables().cellData(eIdxGlobal).setVeModel(veModel);
                        }
                    else if(volume < volumeVE - eps_)
                        {
                            this->variables().cellData(eIdxGlobal).setVeModel(2);
                        }
                }
            updateModelVector();
        }

        void init()
        {
            ParentType::init();

            // store initial CPU time (Zakaria 22-08-25)
            beginCPU_ = std::chrono::high_resolution_clock::now();
            ////

            bool plotFluidMatrixInteractions = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions);

            // plot the Pc-Sw curves, if requested
            if(plotFluidMatrixInteractions)
                this->spatialParams().plotMaterialLaw();

            // initialize the Newton first guess to calculate the gas plume height (Zakaria 22-08-25)
            Scalar domainHeight = this->bBoxMax()[dim - 1];
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            int numberOfColumns = numberOfCells[0];
            gasPlumeDistTemp_.resize(numberOfColumns); // temporary value of gasPlumeDist to be updated at each time step
            for (int i = 0; i < numberOfColumns; i++) 
                {
                    gasPlumeDistTemp_[i] = domainHeight/2.0; // we take half the domain height
                }
            ////

            // initialize column counter and cell center vector used to store solution at only one time step (Zakaria 22-08-25)
            colCounter_ = 0;
            const int maxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
            zCenter_.resize(std::pow(2,maxLevel));
            ////
        }

        void postTimeStep()
        {
            ParentType::postTimeStep();

            // store CPU times (Zakaria 2022-08-25)
            auto endCPU = std::chrono::high_resolution_clock::now(); // time of end of simulation
            auto elapsedCPU = std::chrono::duration_cast<std::chrono::nanoseconds>(endCPU - beginCPU_); // duration of simulation
            auto beginPTS = std::chrono::high_resolution_clock::now(); // start of postTimeStep
            ////

            //plot reconstructed solution
            if(this->timeManager().timeStepIndex() % outputInterval_ == 0 || this->timeManager().willBeFinished()
               || this->timeManager().episodeWillBeFinished() || this->timeManager().timeStepIndex() == 0)
                {
                    int reconstruction = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, Reconstruction);
                    const int maxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
                    LeafGridView gridView = this->gridView();

                    CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
                    const double deltaX = this->bBoxMax()[0]/numberOfCells[0];
                    const double deltaZ = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, reconstruction));
                    const double volume2D = this->bBoxMax()[dim-1] * deltaX /(std::pow(2,maxLevel));//TODO: 3-D
                    //refine VE columns
                    std::unordered_map<IdType, int> mapGlobalIdx;
                    GridCreator::grid().preAdapt();
                    for(int i=0;i<reconstruction;i++)
                        {
                            for (const auto& element : Dune::elements(gridView))
                                {
                                    if(i == 0)
                                        {
                                            IdType globalId = gridView.grid().globalIdSet().id(element);
                                            int eIdxGlobal = this->variables().index(element);
                                            mapGlobalIdx.insert(std::pair<IdType,int>(globalId, eIdxGlobal));
                                        }
                                    GlobalPosition globalPos = element.geometry().center();
                                    int columnNumber = round((globalPos[0] - (deltaX/2.0))/deltaX);
                                    if(modelVector_[columnNumber] == 0 || modelVector_[columnNumber] == 1)
                                        GridCreator::grid().mark(element, UG::D2::BLUE, 1);
                                }
                            // adapt the grid
                            GridCreator::grid().adapt();
                        }
                    GridCreator::grid().postAdapt();

                    // use visualization writer
                    visualizationWriter_.gridChanged();
                    visualizationWriter_.beginWrite(this->timeManager().time() + this->timeManager().timeStepSize());
                    //write stuff out
                    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;
                    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

                    int size = gridView.size(0);
                    ScalarSolutionType *pressureW = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *pressureNw = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *potentialW = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *potentialNw = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *capillaryPressure = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *saturationW = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *saturationNw = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *mobilityW = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *mobilityNw = visualizationWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *veModel = visualizationWriter_.allocateManagedBuffer (size);

                    int i = 0;
                    for (const auto& element : Dune::elements(gridView))
                        {
                            //identify column number of element and check if it is a VE column or a full-d column
                            GlobalPosition globalPos = element.geometry().center();
                            int columnNumber = round((globalPos[0] - (deltaX/2.0))/deltaX);//starting with 0

                            //full-d elements
                            if(modelVector_[columnNumber] == 2)
                                {
                                    IdType globalId = gridView.grid().globalIdSet().id(element);
                                    int eIdxGlobal2D = mapGlobalIdx.find(globalId)->second;
                                    CellData& cellData = this->variables().cellData(eIdxGlobal2D);
                                    (*pressureW)[i]  = cellData.pressure(wPhaseIdx);
                                    (*pressureNw)[i] = cellData.pressure(nPhaseIdx);
                                    (*potentialW)[i]  = cellData.potential(wPhaseIdx);
                                    (*potentialNw)[i] = cellData.potential(nPhaseIdx);
                                    (*capillaryPressure)[i] = cellData.capillaryPressure();
                                    (*saturationW)[i] = cellData.saturation(wPhaseIdx);
                                    (*saturationNw)[i] = cellData.saturation(nPhaseIdx);
                                    (*mobilityW)[i] = cellData.mobility(wPhaseIdx);
                                    (*mobilityNw)[i] = cellData.mobility(nPhaseIdx);
                                    (*veModel)[i] = cellData.veModel();
                                }
                            else
                                {
                                    Element veElement = mapAllColumns_.find(columnNumber)->second;
                                    int eIdxGlobalVE = this->variables().index(veElement);
                                    CellData& cellData = this->variables().cellData(eIdxGlobalVE);

                                    Scalar top = globalPos[dim - 1] + deltaZ/2.0;
                                    Scalar bottom = globalPos[dim - 1] - deltaZ/2.0;

                                    Scalar pRef = referencePressureAtPos(globalPos);
                                    Scalar temp = temperatureAtPos(globalPos);

                                    (*pressureW)[i] = this->pressureModel().reconstPressure(globalPos[dim-1], wPhaseIdx, veElement);
                                    (*pressureNw)[i] = this->pressureModel().reconstPressure(globalPos[dim-1], nPhaseIdx, veElement);
                                    (*potentialW)[i] = this->pressureModel().reconstPotential(globalPos[dim-1], wPhaseIdx, veElement);
                                    (*potentialNw)[i] = this->pressureModel().reconstPotential(globalPos[dim-1], nPhaseIdx, veElement);
                                    (*capillaryPressure)[i] = this->pressureModel().reconstCapillaryPressure(globalPos[dim-1], veElement);
                                    (*saturationW)[i] = this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
                                    (*saturationNw)[i] = 1.0 - this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
                                    (*mobilityW)[i] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, wPhaseIdx)/WettingPhase::viscosity(temp, pRef);
                                    (*mobilityNw)[i] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, nPhaseIdx)/NonWettingPhase::viscosity(temp, pRef);
                                    (*veModel)[i] = cellData.veModel();
                                }
                            i++;
                        }
                    visualizationWriter_.attachCellData(*pressureW, "wetting pressure");
                    visualizationWriter_.attachCellData(*pressureNw, "non-wetting pressure");
                    visualizationWriter_.attachCellData(*potentialW, "wetting potential");
                    visualizationWriter_.attachCellData(*potentialNw, "non-wetting potential");
                    visualizationWriter_.attachCellData(*capillaryPressure, "capillary pressure");
                    visualizationWriter_.attachCellData(*saturationW, "wetting saturation");
                    visualizationWriter_.attachCellData(*saturationNw, "non-wetting saturation");
                    visualizationWriter_.attachCellData(*mobilityW, "wetting mobility");
                    visualizationWriter_.attachCellData(*mobilityNw, "non-wetting mobility");
                    visualizationWriter_.attachCellData(*veModel, "veModel");
                    visualizationWriter_.endWrite();


                    //        //refine grid for error norm calculation
                    //        GridCreator::grid().preAdapt();
                    //        const int additionalRefinementSteps = 2;
                    //        for(int i=0;i<additionalRefinementSteps;i++)
                    //        {
                    //            for (const auto& element : Dune::elements(gridView))
                    //            {
                    //                GlobalPosition globalPos = element.geometry().center();
                    //                GridCreator::grid().mark(1, element);
                    //            }
                    //            // adapt the grid
                    //            GridCreator::grid().adapt();
                    //            GridCreator::grid().postAdapt();
                    //        }
                    //        const double deltaZRefined = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, (reconstruction+additionalRefinementSteps)));
                    //
                    //        // use error norm writer
                    //        normWriter_.gridChanged();
                    //        normWriter_.beginWrite(this->timeManager().time() + this->timeManager().timeStepSize());
                    //
                    //        size = gridView.size(0);
                    //        ScalarSolutionType *pressureW1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *pressureNw1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *potentialW1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *potentialNw1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *capillaryPressure1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *saturationW1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *saturationNw1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *mobilityW1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *mobilityNw1 = normWriter_.allocateManagedBuffer (size);
                    //        ScalarSolutionType *veModel1 = normWriter_.allocateManagedBuffer (size);
                    //
                    //        //write output on refined grid
                    //        i = 0;
                    //        for (const auto& element : Dune::elements(gridView))
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
                    //            //full-d elements
                    //            if(modelVector_[columnNumber] == 2)
                    //            {
                    //                IdType globalId = gridView.grid().globalIdSet().id(element.father().father());
                    //                int eIdxGlobal2D = mapGlobalIdx.find(globalId)->second;
                    //                CellData& cellData = this->variables().cellData(eIdxGlobal2D);
                    //                (*pressureW1)[globalID]  = cellData.pressure(wPhaseIdx);
                    //                (*pressureNw1)[globalID] = cellData.pressure(nPhaseIdx);
                    //                (*potentialW1)[globalID]  = cellData.potential(wPhaseIdx);
                    //                (*potentialNw1)[globalID] = cellData.potential(nPhaseIdx);
                    //                (*capillaryPressure1)[globalID] = cellData.capillaryPressure();
                    //                (*saturationW1)[globalID] = cellData.saturation(wPhaseIdx);
                    //                (*saturationNw1)[globalID] = cellData.saturation(nPhaseIdx);
                    //                (*mobilityW1)[globalID] = cellData.mobility(wPhaseIdx);
                    //                (*mobilityNw1)[globalID] = cellData.mobility(nPhaseIdx);
                    //                (*veModel1)[globalID] = cellData.veModel();
                    //            }
                    //            else
                    //            {
                    //                Element veElement = mapAllColumns_.find(columnNumber)->second;
                    //                int eIdxGlobalVE = this->variables().index(veElement);
                    //                CellData& cellData = this->variables().cellData(eIdxGlobalVE);
                    //
                    //                Scalar top = globalPos[dim - 1] + deltaZRefined/2.0;
                    //                Scalar bottom = globalPos[dim - 1] - deltaZRefined/2.0;
                    //
                    //                Scalar pRef = referencePressureAtPos(globalPos);
                    //                Scalar temp = temperatureAtPos(globalPos);
                    //
                    //                (*pressureW1)[globalID] = this->pressureModel().reconstPressure(globalPos[dim-1], wPhaseIdx, veElement);
                    //                (*pressureNw1)[globalID] = this->pressureModel().reconstPressure(globalPos[dim-1], nPhaseIdx, veElement);
                    //                (*potentialW1)[globalID] = this->pressureModel().reconstPotential(globalPos[dim-1], wPhaseIdx, veElement);
                    //                (*potentialNw1)[globalID] = this->pressureModel().reconstPotential(globalPos[dim-1], nPhaseIdx, veElement);
                    //                (*capillaryPressure1)[globalID] = this->pressureModel().reconstCapillaryPressure(globalPos[dim-1], veElement);
                    //                (*saturationW1)[globalID] = this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
                    //                (*saturationNw1)[globalID] = 1.0 - this->pressureModel().saturationIntegral(bottom, top, veElement)/(top-bottom);
                    //                (*mobilityW1)[globalID] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, wPhaseIdx)/WettingPhase::viscosity(temp, pRef);
                    //                (*mobilityNw1)[globalID] = this->pressureModel().calculateRelPermeabilityCoarse(bottom, top, veElement, nPhaseIdx)/NonWettingPhase::viscosity(temp, pRef);
                    //                (*veModel1)[globalID] = cellData.veModel();
                    //            }
                    //            i++;
                    //        }
                    //        normWriter_.attachCellData(*pressureW1, "wetting pressure");
                    //        normWriter_.attachCellData(*pressureNw1, "non-wetting pressure");
                    //        normWriter_.attachCellData(*potentialW1, "wetting potential");
                    //        normWriter_.attachCellData(*potentialNw1, "non-wetting potential");
                    //        normWriter_.attachCellData(*capillaryPressure1, "capillary pressure");
                    //        normWriter_.attachCellData(*saturationW1, "wetting saturation");
                    //        normWriter_.attachCellData(*saturationNw1, "non-wetting saturation");
                    //        normWriter_.attachCellData(*mobilityW1, "wetting mobility");
                    //        normWriter_.attachCellData(*mobilityNw1, "non-wetting mobility");
                    //        normWriter_.attachCellData(*veModel1, "veModel");
                    //        normWriter_.endWrite();
                    //
                    //        //coarsen grid after error norm calculation
                    //        for(int i=0;i<additionalRefinementSteps;i++)
                    //        {
                    //            for (const auto& element : Dune::elements(gridView))
                    //            {
                    //                GlobalPosition globalPos = element.geometry().center();
                    //                GridCreator::grid().mark(-1, element);
                    //            }
                    //            // adapt the grid
                    //            GridCreator::grid().adapt();
                    //            GridCreator::grid().postAdapt();
                    //            GridCreator::grid().preAdapt();
                    //        }

                    //coarsen ve-columns
                    for(int i=0;i<reconstruction;i++)
                        {
                            for (const auto& element : Dune::elements(gridView))
                                {
                                    GlobalPosition globalPos = element.geometry().center();
                                    int columnNumber = round((globalPos[0] - (deltaX/2.0))/deltaX);
                                    if(modelVector_[columnNumber] == 0 || modelVector_[columnNumber] == 1)
                                        GridCreator::grid().mark(-1, element);
                                }
                            // adapt the grid
                            GridCreator::grid().preAdapt();
                            GridCreator::grid().adapt();
                            GridCreator::grid().postAdapt();
                        }
                }

            // store time in result files (Zakaria 2022-08-25)
            outputFile_.open("timeCriterion20.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("timeZp20.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("timeCriterion200.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("timeZp200.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("criterion.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("zp.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("CPUTime.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("CPUNewton.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            outputFile_.open("iterNumberNewton.out", std::ios::app);
            outputFile_ << this->timeManager().time()/segTime_;
            outputFile_.close();
            ////
            
            //write out average column saturation
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            std::vector<Scalar> averageSatColumn(0);
            Scalar gasPlumeVolume = 0.0;// total volume of gas plume
            averageSatColumn.resize(numberOfCells[0]);
            for (int i = 0; i != averageSatColumn.size(); ++i)
                {
                    Scalar totalVolume = 0.0;// total volume of column
                    typename std::map<int, Element>::iterator it = mapAllColumns_.lower_bound(i);
                    for (; it != mapAllColumns_.upper_bound(i); ++it)
                        {
                            int globalIdxI = this->variables().index(it->second);
                            GlobalPosition globalPos = (it->second).geometry().center();
                            Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
                            Scalar volume = it->second.geometry().volume();
                            averageSatColumn[i] += satW * volume;
                            totalVolume += volume;

                            if(globalPos[0]<=TwoDArea_ && satW<1.0-eps_)//attention!
                                {
                                    averageSatInPlume_ += satW * volume;
                                    gasPlumeVolume += volume;
                                }
                        }
                    averageSatColumn[i] = averageSatColumn[i]/totalVolume;//average wetting saturation in column (equals gasPlumeDist for SI and no compressibility)
                    // outputFile_.open("averageSatColumn.out", std::ios::app);
                    // outputFile_ << " " << averageSatColumn[i];
                    // outputFile_.close();
                }
            //    averageSatInPlume_ = averageSatInPlume_/gasPlumeVolume;//average wetting saturation in plume
            //    if(gasPlumeVolume<eps_)
            //    {
            //        averageSatInPlume_ = 1.0;
            //    }
            //
            //    //calculate error for certain column
            //    Scalar errorSat(0.0);
            //    Scalar gasPlumeHeight2D(0.0);
            //
            //    int TwoDRefinement = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, 2DRefinement);
            //    double deltaZ = this->bBoxMax()[dim - 1]/(numberOfCells[dim - 1]*std::pow(2, TwoDRefinement));
            //
            //    int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            //
            //    int columnNumber = 9;
            //    typename std::map<int, Element>::iterator it = mapAllColumns_.lower_bound(columnNumber);
            //    Scalar gasPlumeDist = calculateGasPlumeDist(it->second, averageSat[columnNumber]);
            //
            //    for (; it != mapAllColumns_.upper_bound(columnNumber); ++it)
            //    {
            //        int globalIdxI = this->variables().index(it->second);
            //        GlobalPosition globalPos = (it->second).geometry().center();
            //        Scalar top = globalPos[dim - 1] + deltaZ/2.0;
            //        Scalar bottom = globalPos[dim - 1] - deltaZ/2.0;
            //
            //        Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
            //        Scalar resSatW = this->spatialParams().materialLawParams(it->second).swr();
            //
            //        if(veModel == sharpInterface)//calculate error for VE model
            //        {
            //            if (top <= gasPlumeDist)
            //            {
            //                errorSat += std::abs(deltaZ * (satW - 1.0));
            //            }
            //            else if (bottom >= gasPlumeDist)
            //            {
            //                errorSat += std::abs(deltaZ * (satW - resSatW));
            //
            //            }
            //            else
            //            {
            //                Scalar lowerDelta = gasPlumeDist - bottom;
            //                Scalar upperDelta = top - gasPlumeDist;
            //                errorSat += std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW));
            //            }
            //        }
            //        else if(veModel == capillaryFringe)//calculate error for capillary fringe model
            //        {
            //            if (top <= gasPlumeDist)
            //            {
            //                errorSat += std::abs(deltaZ * (satW - 1.0));
            //            }
            //            else if (bottom >= gasPlumeDist)
            //            {
            //                errorSat += calculateErrorSatIntegral(bottom, top, satW, gasPlumeDist);
            //            }
            //            else
            //            {
            //                Scalar lowerDelta = gasPlumeDist - bottom;
            //                Scalar upperDelta = top - gasPlumeDist;
            //                errorSat += std::abs(lowerDelta * (satW - 1.0)) + calculateErrorSatIntegral(gasPlumeDist, top, satW, gasPlumeDist);
            //            }
            //        }
            //        if(satW<0.99)
            //        {
            //            gasPlumeHeight2D += (top-bottom);
            //        }
            //    }
            //
            //    errorSat = errorSat/(this->bBoxMax()[dim - 1]-gasPlumeDist);
            //
            //    if(averageSat[columnNumber]>1.0-eps_)
            //    {
            //        errorSat = 0.0;
            //    }
            //
            //    std::cout << "errorSat " << errorSat;
            //
            //    outputFile_.open("errorSat.out", std::ios::app);
            //    outputFile_ << this->timeManager().time() << " " << errorSat << std::endl;
            //    outputFile_.close();
            //


            // calculation of criteria (Zakaria 2022-08-25 and 2022-08-30)
            int numberOfColumns =numberOfCells[0];
            criterion_.resize(numberOfColumns);
            for (int i = 0; i < numberOfColumns; i++) {
                criterion_[i] = 0.0;
            }

            const int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            const int maxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
            const double endTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeManager, TEnd);

            for (int i=0;i<numberOfColumns;i++) {
                Scalar domainHeight = this->bBoxMax()[dim - 1];

                //// set the Newton initial guess to calculate plume height (Zakaria 22-08-25)
                typename std::map<int, Element>::iterator it = mapAllColumns_.lower_bound(i);                
                Scalar gasPlumeDist = calculateGasPlumeDist(it->second,averageSatColumn[i],gasPlumeDistTemp_[i]);
                //                gasPlumeDistTemp_[i] = std::min(0.99*domainHeight,gasPlumeDist); // update new temporary height
                gasPlumeDistTemp_[i] = domainHeight/2.0;
                ////////

                //// store bottom pressure of column i to compute hydrostatic pressure later (Zakaria 22-08-25)	 
                Element veElement = mapAllColumns_.lower_bound(i)->second;
                int eIdxGlobal = this->variables().index(veElement);
                Scalar coarsePressW = this->variables().cellData(eIdxGlobal).pressure(wPhaseIdx);
                ////////
                
                if (i < 4) // never have VE in the first four columns
                    {
                        criterion_[i] = 2;
                        
                        outputFile_.open("criterion.out", std::ios::app);
                        outputFile_ << " " << criterion_[i];
                        outputFile_.close();
                        outputFile_.open("zp.out", std::ios::app);
                        outputFile_ << " " << gasPlumeDist ;
                        outputFile_.close();

                        //// plot profiles at final time (Zakaria 2022-12-25)
                        if (this->timeManager().time() + this->timeManager().timeStepSize() >= endTime && colCounter_ < numberOfColumns)
                            {
                                colCounter_ += 1;
                                if (i == 0)
                                    {
                                        // iterate over cells in column and write z-location and store centers of cells
                                        int k = 0;
                                        typename std::map<int, Element>::iterator it1 = mapAllColumns_.lower_bound(i);
                                        for (; it1 != mapAllColumns_.upper_bound(i); ++it1)
                                            {
                                                GlobalPosition globalPos = (it1->second).geometry().center();
                                                Scalar z = globalPos[dim - 1];
                                                outputFile_.open("satProfiles.out", std::ios::app);
                                                outputFile_ << z << " ";
                                                outputFile_.close();
                                                outputFile_.open("relPermProfiles.out", std::ios::app);
                                                outputFile_ << z << " ";
                                                outputFile_.close();
                                                outputFile_.open("pressProfiles.out", std::ios::app);
                                                outputFile_ << z << " ";
                                                outputFile_.close();
                                                outputFile_.open("capPressProfiles.out", std::ios::app);
                                                outputFile_ << z << " ";
                                                outputFile_.close();

                                                zCenter_[k] = z;
                                                k += 1;
                                            }
                                    }

                                // next line
                                outputFile_.open("satProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("relPermProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("pressProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("capPressProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();

                                // iterate over cells in column and store results
                                typename std::map<int, Element>::iterator it2 = mapAllColumns_.lower_bound(i);
                                for (; it2 != mapAllColumns_.upper_bound(i); ++it2) {
                                    int globalIdxI = this->variables().index(it2->second);
                                    Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
                                    Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(it2->second), satW);
                                    Scalar pressW = this->variables().cellData(globalIdxI).pressure(wPhaseIdx);
                                    Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(it2->second), satW);
                                    outputFile_.open("satProfiles.out", std::ios::app);
                                    outputFile_ << satW << " ";
                                    outputFile_.close();
                                    outputFile_.open("relPermProfiles.out", std::ios::app);
                                    outputFile_ << krw << " ";
                                    outputFile_.close();
                                    outputFile_.open("pressProfiles.out", std::ios::app);
                                    outputFile_ << pressW << " ";
                                    outputFile_.close();
                                    outputFile_.open("capPressProfiles.out", std::ios::app);
                                    outputFile_ << pc << " ";
                                    outputFile_.close();
                                }
                            }
                        ////////
                        continue;
                    }

                int columnModel = modelVector_[i];
                if (columnModel == 0 || columnModel == 1) 
                    {
                        outputFile_.open("criterion.out", std::ios::app);
                        outputFile_ << " " << criterion_[i]; // should be zero
                        outputFile_.close();
                        outputFile_.open("zp.out", std::ios::app);
                        outputFile_ << " " << gasPlumeDist ;
                        outputFile_.close();

                        //// store errors and plume height in files for 20th and 200th columns (Zakaria 2022-08-25)
                        if (i == 19)
                            {
                                outputFile_.open("timeCriterion20.out", std::ios::app);
                                outputFile_ << " " << criterion_[i] << std::endl;
                                outputFile_.close();
                                outputFile_.open("timeZp20.out", std::ios::app);
                                outputFile_ << " " << gasPlumeDist << std::endl;
                                outputFile_.close();
                            }
                        if (i == 199)
                            {
                                outputFile_.open("timeCriterion200.out", std::ios::app);
                                outputFile_ << " " << criterion_[i] << std::endl;
                                outputFile_.close();
                                outputFile_.open("timeZp200.out", std::ios::app);
                                outputFile_ << " " << gasPlumeDist << std::endl;
                                outputFile_.close();
                            }
                        ////////

                        //// plot profiles at final time (Zakaria 2022-12-25)
                        if (this->timeManager().time() + this->timeManager().timeStepSize() >= endTime && colCounter_ < numberOfColumns)
                            {
                                colCounter_ += 1;

                                outputFile_.open("satProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("relPermProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("pressProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();
                                outputFile_.open("capPressProfiles.out", std::ios::app);
                                outputFile_ << std::endl;
                                outputFile_.close();

                                // store results (only one cell here)
                                for (int k = 0; k < std::pow(2,maxLevel); ++k) 
                                    {
                                        Scalar satW = reconstSaturation(zCenter_[k],gasPlumeDist);
                                        Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(it->second), satW);
                                        Element veElement = mapAllColumns_.lower_bound(i)->second;
                                        Scalar pressW = reconstPressureW(zCenter_[k],gasPlumeDist,coarsePressW);
                                        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(it->second), satW);
                                        outputFile_.open("satProfiles.out", std::ios::app);
                                        outputFile_ << satW << " ";
                                        outputFile_.close();
                                        outputFile_.open("relPermProfiles.out", std::ios::app);
                                        outputFile_ << krw << " ";
                                        outputFile_.close();
                                        outputFile_.open("pressProfiles.out", std::ios::app);
                                        outputFile_ << pressW << " ";
                                        outputFile_.close();
                                        outputFile_.open("capPressProfiles.out", std::ios::app);
                                        outputFile_ << pc << " ";
                                        outputFile_.close();
                                    }
                            }
                        ////////
                        continue;
                    }


                //// initialize all criteria to select from (Zakaria 2022-12-25 and 2022-08-30)
                Scalar errorSat = 0;
                Scalar errorSatNorm = 0;
                Scalar errorRelPerm = 0;
                Scalar errorRelPermNorm = 0;
                Scalar errorPress = 0;
                Scalar errorPressNorm = 0;
                Scalar errorCapPress = 0;
                Scalar errorCapPressNorm = 0;
                ////////

                double deltaZ = domainHeight/(numberOfCells[dim - 1]*std::pow(2, maxLevel));
                Scalar resSatW = this->spatialParams().materialLawParams(it->second).swr();

                Scalar bottomPressW = reconstPressureW(-deltaZ/2,gasPlumeDist,coarsePressW); //slightly different from coarsePressW when columnModel > 1 since coarsePressW is in this case the pressure at the center of the bottom cell and not the pressure at bottom

                for (; it != mapAllColumns_.upper_bound(i); ++it) {
                    Element element = it->second;
                    int globalIdxI = this->variables().index(element);
                    GlobalPosition globalPos = element.geometry().center();
                    Scalar top = globalPos[dim - 1] + deltaZ / 2.0;
                    Scalar bottom = globalPos[dim - 1] - deltaZ / 2.0;

                    Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
                    Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(element), satW);
                    Scalar pressW = this->variables().cellData(globalIdxI).pressure(wPhaseIdx);
                    Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(element), satW);
                    Scalar pe = this->spatialParams().materialLawParams(element).pe();

                    if (veModel == sharpInterface) // calculate error for VE model
                        {
                            if (top <= gasPlumeDist) {
                                errorSat += std::abs(deltaZ * (satW - 1.0));
                                errorRelPerm += std::abs(deltaZ * (krw - 1.0));
                                errorSatNorm += (1/satW)* std::abs(deltaZ * (satW - 1.0));
                                errorRelPerm += (1/krw)* std::abs(deltaZ * (krw - 1.0));
                            } else if (bottom >= gasPlumeDist) {
                                errorSat += std::abs(deltaZ * (satW - resSatW));
                                errorRelPerm += std::abs(deltaZ * (krw - 0.0));
                                errorSatNorm += (1/satW)* std::abs(deltaZ * (satW - resSatW));
                                errorRelPermNorm += (1/krw)* std::abs(deltaZ * (krw - 0.0));
                            } else {
                                Scalar lowerDelta = gasPlumeDist - bottom;
                                Scalar upperDelta = top - gasPlumeDist;
                                errorSat += std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW));
                                errorRelPerm += std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0));
                                errorSatNorm += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW)));
                                errorRelPermNorm += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0)));
                            }
                        }
                    if (veModel == capillaryFringe) // calculate error for capillary fringe model
                        {
                            if (top <= gasPlumeDist) {
                                errorSat += std::abs(deltaZ * (satW - 1.0));
                                errorRelPerm += std::abs(deltaZ * (krw - 1.0));
                                errorSatNorm +=(1/satW)* ( std::abs(deltaZ * (satW - 1.0)) ) ;
                                errorRelPermNorm += (1/krw)* ( std::abs(deltaZ * (krw - 1.0)) ) ;
                                errorCapPress += std::abs(deltaZ * (pc - pe));
                                errorCapPressNorm += (1/pc)* std::abs(deltaZ * (pc - pe));
                            } else if (bottom >= gasPlumeDist) {
                                Scalar eIntSat = calculateErrorSatIntegral(bottom, top, satW, gasPlumeDist);
                                Scalar eIntPerm = calculateErrorKrwIntegral(bottom, top, satW, gasPlumeDist);
                                Scalar eIntPc = calculateErrorCapPressIntegral(bottom, top, pc, gasPlumeDist);
                                errorSat += eIntSat;
                                errorRelPerm += eIntPerm;
                                errorSatNorm += (1/satW)* eIntSat;
                                errorRelPermNorm += (1/krw)* eIntPerm;
                                errorCapPress += eIntPc;
                                errorCapPressNorm += (1/pc)* eIntPc;
                            } else {
                                Scalar lowerDelta = gasPlumeDist - bottom;
                                Scalar upperDelta = top - gasPlumeDist;

                                Scalar eIntSat = calculateErrorSatIntegral(gasPlumeDist, top, satW, gasPlumeDist);
                                Scalar eIntPerm = calculateErrorKrwIntegral(gasPlumeDist, top, satW, gasPlumeDist);
                                Scalar eIntPc = calculateErrorCapPressIntegral(gasPlumeDist, top, pc, gasPlumeDist);
                                errorSat += std::abs(lowerDelta * (satW - 1.0)) + eIntSat;
                                errorRelPerm += std::abs(lowerDelta * (krw - 1.0)) + eIntPerm;
                                errorSatNorm += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + eIntSat);
                                errorRelPerm += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + eIntPerm); 
                                errorCapPress += std::abs(lowerDelta* (pc - pe)) + eIntPc;
                                errorCapPressNorm += (1/pc)* (std::abs(lowerDelta* (pc - pe)) + eIntPc);
                            }
                        }
                    Scalar eIntPress = calculateErrorPressIntegral(bottom, top, pressW, gasPlumeDist, bottomPressW);
                    errorPress += eIntPress;
                    errorPressNorm += (1/pressW)* eIntPress;
                }


                //// select which criterion is being used (Zakaria 2022-08-25)
                //criterion_[i] = errorSat / (domainHeight - gasPlumeDist);
                criterion_[i] = errorRelPerm / (domainHeight - gasPlumeDist);
                //criterion_[i] =  errorPress / (domainHeight - gasPlumeDist);
                //criterion_[i] = errorSatNorm / (domainHeight - gasPlumeDist);
                //criterion_[i] = errorRelPermNorm / (domainHeight - gasPlumeDist);
                //criterion_[i] =  errorPressNorm / (domainHeight - gasPlumeDist);
                //criterion_[i] =  errorCapPress / (domainHeight - gasPlumeDist);
                //criterion_[i] =  errorCapPressNorm / (domainHeight - gasPlumeDist);
                ////////
                
                //// store errors and plume height in files for 20th and 200th columns (Zakaria 2022-08-25)
                if (i == 19)
                    {
                        outputFile_.open("timeCriterion20.out", std::ios::app);
                        outputFile_ << " " << criterion_[i] << std::endl;
                        outputFile_.close();
                        outputFile_.open("timeZp20.out", std::ios::app);
                        outputFile_ << " " << gasPlumeDist << std::endl;
                        outputFile_.close();
                    }
                if (i == 199)
                    {
                        outputFile_.open("timeCriterion200.out", std::ios::app);
                        outputFile_ << " " << criterion_[i] << std::endl;
                        outputFile_.close();
                        outputFile_.open("timeZp200.out", std::ios::app);
                        outputFile_ << " " << gasPlumeDist << std::endl;
                        outputFile_.close();
                    }
                ////////

                //// store errors and plume height in files for all columns (Zakaria 2022-08-25)
                outputFile_.open("criterion.out", std::ios::app);
                outputFile_ << " " << criterion_[i];
                outputFile_.close();
                outputFile_.open("zp.out", std::ios::app);
                outputFile_ << " " << gasPlumeDist ;
                outputFile_.close();
                ////////

                //// plot profiles at final time (Zakaria 2022-12-25)
                if (this->timeManager().time() + this->timeManager().timeStepSize() >= endTime && colCounter_ < numberOfColumns)
                    {
                        colCounter_ += 1;

                        outputFile_.open("satProfiles.out", std::ios::app);
                        outputFile_ << std::endl;
                        outputFile_.close();
                        outputFile_.open("relPermProfiles.out", std::ios::app);
                        outputFile_ << std::endl;
                        outputFile_.close();
                        outputFile_.open("pressProfiles.out", std::ios::app);
                        outputFile_ << std::endl;
                        outputFile_.close();
                        outputFile_.open("capPressProfiles.out", std::ios::app);
                        outputFile_ << std::endl;
                        outputFile_.close();

                        // store results (only one cell here)
                        // store results
                        typename std::map<int, Element>::iterator it2 = mapAllColumns_.lower_bound(i);
                        for (; it2 != mapAllColumns_.upper_bound(i); ++it2)
                            {
                                int globalIdxI = this->variables().index(it2->second);
                                Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
                                Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(it2->second), satW);
                                Scalar pressW = this->variables().cellData(globalIdxI).pressure(wPhaseIdx);
                                Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(it2->second), satW);
                                outputFile_.open("satProfiles.out", std::ios::app);
                                outputFile_ << satW << " ";
                                outputFile_.close();
                                outputFile_.open("relPermProfiles.out", std::ios::app);
                                outputFile_ << krw << " ";
                                outputFile_.close();
                                outputFile_.open("pressProfiles.out", std::ios::app);
                                outputFile_ << pressW << " ";
                                outputFile_.close();
                                outputFile_.open("capPressProfiles.out", std::ios::app);
                                outputFile_ << pc << " ";
                                outputFile_.close();
                            }
                    }
                ////////
            }

            //// finalize result files (Zakaria 2022-08-25)
            outputFile_.open("criterion.out", std::ios::app);
            outputFile_ << std::endl;
            outputFile_.close();
            outputFile_.open("zp.out", std::ios::app);
            outputFile_ << std::endl ;
            outputFile_.close();
            outputFile_.open("CPUNewton.out", std::ios::app);
            outputFile_ << std::endl;
            outputFile_.close();
            outputFile_.open("iterNumberNewton.out", std::ios::app);
            outputFile_ << std::endl;
            outputFile_.close();
            ////////
            //// end of calculation of criteria (Zakaria 2022-08-25)


            //check mass conservativity:
            GridView GridView = this->gridView();
            Scalar totalMassN = 0.0;
            int numberOfCellsTotal = 0;
            for (const auto& element : Dune::elements(GridView))
                {
                    int eIdxGlobal = this->variables().index(element);
                    CellData& cellData = this->variables().cellData(eIdxGlobal);
                    GlobalPosition globalPos = element.geometry().center();
                    Scalar pRef = referencePressureAtPos(globalPos);
                    Scalar temp = temperatureAtPos(globalPos);
                    Scalar massN = cellData.saturation(nPhaseIdx) * element.geometry().volume() * this->spatialParams().porosity(element) * NonWettingPhase::density(temp, pRef);
                    totalMassN += massN;
                    numberOfCellsTotal += 1;
                }
            Scalar totalMassNExpected = -GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, InjectionrateN) * (this->timeManager().time()+this->timeManager().timeStepSize()) * this->bBoxMax()[1];
            std::cout << "error " << totalMassNExpected - totalMassN;

            outputFile_.open("numberOfCells.out", std::ios::app);
            outputFile_ << this->timeManager().time() << " " << numberOfCellsTotal << std::endl;
            outputFile_.close();


            // postTimeStep times (Zakaria 2022-08-25)
            auto endPTS = std::chrono::high_resolution_clock::now(); // end time of postTimeStep
            auto elapsedPTS = std::chrono::duration_cast<std::chrono::nanoseconds>(endPTS - beginPTS); // duration of postTimeStep
            ////
            
            // store duration of simulation and postTimeStep (Zakaria 2022-08-25)
            outputFile_.open("CPUTime.out", std::ios::app);
            outputFile_ << " " << elapsedCPU.count() << " " << elapsedPTS.count() << std::endl;
            outputFile_.close();
            ////
        }

        /*!
         * \brief Capability to introduce problem-specific routines after grid adaptation
         *
         * Function is called at the end of the standard grid
         * modification routine, GridAdapt::adaptGrid() , to allow
         * for problem-specific output etc.
         */
        void debugPlot()
        {
            // write out new grid
            if(this->timeManager().timeStepIndex() % outputInterval_ == 0)
                {
                    debugWriter_.gridChanged();
                    debugWriter_.beginWrite(this->timeManager().time());
                    //write stuff out
                    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;
                    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

                    int size = this->gridView().size(0);
                    ScalarSolutionType *pressureW = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *pressureNw = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *potentialW = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *potentialNw = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *capillaryPressure = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *saturationW = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *saturationNw = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *mobilityW = debugWriter_.allocateManagedBuffer (size);
                    ScalarSolutionType *mobilityNw = debugWriter_.allocateManagedBuffer (size);

                    int i = 0;
                    for (const auto& element : Dune::elements(this->gridView()))
                        {
                            //identify column number of element and check if it is a VE column or a full-d column
                            int eIdxGlobal = this->variables().index(element);
                            CellData& cellData = this->variables().cellData(eIdxGlobal);
                            (*pressureW)[i]  = cellData.pressure(wPhaseIdx);
                            (*pressureNw)[i] = cellData.pressure(nPhaseIdx);
                            (*potentialW)[i]  = cellData.potential(wPhaseIdx);
                            (*potentialNw)[i] = cellData.potential(nPhaseIdx);
                            (*capillaryPressure)[i] = cellData.capillaryPressure();
                            (*saturationW)[i] = cellData.saturation(wPhaseIdx);
                            (*saturationNw)[i] = cellData.saturation(nPhaseIdx);
                            (*mobilityW)[i] = cellData.mobility(wPhaseIdx);
                            (*mobilityNw)[i] = cellData.mobility(nPhaseIdx);
                            i++;
                        }
                    debugWriter_.attachCellData(*pressureW, "wetting pressure");
                    debugWriter_.attachCellData(*pressureNw, "non-wetting pressure");
                    debugWriter_.attachCellData(*potentialW, "wetting potential");
                    debugWriter_.attachCellData(*potentialNw, "non-wetting potential");
                    debugWriter_.attachCellData(*capillaryPressure, "capillary pressure");
                    debugWriter_.attachCellData(*saturationW, "wetting saturation");
                    debugWriter_.attachCellData(*saturationNw, "non-wetting saturation");
                    debugWriter_.attachCellData(*mobilityW, "wetting mobility");
                    debugWriter_.attachCellData(*mobilityNw, "non-wetting mobility");
                    debugWriter_.endWrite();
                }
        }

        /*! \brief Calculates gasPlumeDist, distance of gas plume from bottom
         *
         * Stores minGasPlumeDist for all grid cells
         */
        Scalar calculateGasPlumeDist(const Element& element, Scalar satW, Scalar initGuess)
        {
            Scalar domainHeight = this->bBoxMax()[dim - 1];
            Scalar resSatW = this->spatialParams().materialLawParams(element).swr();
            Scalar resSatN = this->spatialParams().materialLawParams(element).snr();
            Scalar gravity = this->gravity().two_norm();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar gasPlumeDist = 0.0;

            // store number of iterations in Newton method (Zakaria 2022-08-25)
            int iterNumber = 0;
            ////

            // start time of Newton method (Zakaria 2022-08-25)
            auto beginNewton = std::chrono::high_resolution_clock::now();
            ////

            if (veModel == sharpInterface) //calculate gasPlumeDist for sharp interface ve model
                {
                    gasPlumeDist = domainHeight * (satW - resSatW) / (1.0 - resSatW);
                }

            else if (veModel == capillaryFringe) //calculate gasPlumeDist for capillary fringe model
                {
                    GlobalPosition globalPos = element.geometry().center();
                    Scalar pRef = referencePressureAtPos(globalPos);
                    Scalar tempRef = temperatureAtPos(globalPos);
                    Scalar densityW = WettingPhase::density(tempRef, pRef);
                    Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
                    Scalar lambda = this->spatialParams().materialLawParams(element).lambda();
                    Scalar entryP = this->spatialParams().materialLawParams(element).pe();

                    // start from previous plume height (Zakaria 2022-08-25)
                    Scalar Xi = initGuess;
                    ////

                    Scalar fullIntegral = 1.0 / (1.0 - lambda) * (1.0 - resSatW - resSatN) / ((densityW - densityNw) * gravity) * (std::pow(entryP, lambda)
                                                                                                                                   - std::pow(entryP, 2.0 - lambda) + std::pow((domainHeight * (densityW - densityNw) * gravity + entryP), (1.0 - lambda)));
                    //GasPlumeDist>0
                    if (fullIntegral < satW * domainHeight)
                        {
                            //solve equation for
                            for (int count = 0; count < 100; count++)
                                {
                                    // update number of Newton iterations (Zakaria 2022-08-25)
                                    iterNumber += 1;
                                    ////

                                    Scalar residual = 1.0 / (1.0 - lambda) * std::pow(((domainHeight - Xi) * (densityW - densityNw) * gravity + entryP),(1.0 - lambda))
                                        * (1.0 - resSatW - resSatN) * std::pow(entryP, lambda) / ((densityW - densityNw) * gravity) + resSatW * (domainHeight - Xi)
                                        - entryP * (1.0 - resSatW - resSatN) / ((1.0 - lambda) * (densityW - densityNw) * gravity) + Xi - satW * domainHeight;

                                    if (fabs(residual) < 1e-10)
                                        break;

                                    Scalar derivation = std::pow(((domainHeight - Xi) * (densityW - densityNw) * gravity + entryP), -lambda) * (resSatN + resSatW - 1.0)
                                        * std::pow(entryP, lambda) - resSatW + 1.0;

                                    Xi = Xi - residual / (derivation);
                                }
                        }
                    //GasPlumeDist<0
                    else if (fullIntegral > satW * domainHeight)
                        {
                            //solve equation
                            for (int count = 0; count < 100; count++)
                                {
                                    // update number of Newton iterations (Zakaria 2022-08-25)
                                    iterNumber += 1;
                                    ////

                                    Scalar residual = 1.0 / (1.0 - lambda) * std::pow(((domainHeight - Xi) * (densityW - densityNw) * gravity + entryP),
                                                                                      (1.0 - lambda)) * (1.0 - resSatW - resSatN)* std::pow(entryP, lambda) / ((densityW - densityNw) * gravity)
                                        + resSatW * domainHeight - 1.0 / (1.0 - lambda) * std::pow(((-Xi) * (densityW - densityNw)
                                                                                                    * gravity + entryP),(1.0 - lambda)) * (1.0 - resSatW - resSatN) * std::pow(entryP, lambda) / ((densityW - densityNw)
                                                                                                                                                                                                  * gravity)
                                        - satW * domainHeight;
                                    if (fabs(residual) < 1e-10)
                                        break;

                                    Scalar derivation = std::pow(((domainHeight - Xi) * (densityW - densityNw) * gravity + entryP), -lambda)
                                        * (resSatN + resSatW - 1.0) * std::pow(entryP, lambda) + std::pow(((-Xi) * (densityW - densityNw) * gravity + entryP),
                                                                                                          -lambda) * (1.0 - resSatN - resSatW) * std::pow(entryP, lambda);

                                    Xi = Xi - residual / (derivation);
                                };
                        }
                    //GasPlumeDist=0
                    else
                        {
                            Xi = 0.0;
                        }
                    gasPlumeDist = Xi;
                }

            // store end time for Newton and duration for Newton (Zakaria 2022-08-25)
            auto endNewton = std::chrono::high_resolution_clock::now();
            auto elapsedNewton = std::chrono::duration_cast<std::chrono::nanoseconds>(endNewton - beginNewton);

            outputFile_.open("CPUNewton.out", std::ios::app);
            outputFile_ << " " << elapsedNewton.count();
            outputFile_.close();
            outputFile_.open("iterNumberNewton.out", std::ios::app);
            outputFile_ << " " << iterNumber;
            outputFile_.close();
            ////

            return gasPlumeDist;
        }

        /*! \brief Calculates integral of difference of wetting saturation over z
         *
         */
        Scalar calculateErrorSatIntegral(Scalar lowerBound, Scalar upperBound, Scalar satW, Scalar gasPlumeDist)
        {
            int intervalNumber = 10;
            Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;

            Scalar satIntegral = 0.0;
            for(int count=0; count<intervalNumber; count++ )
                {
                    satIntegral += std::abs((reconstSaturation(lowerBound + count*deltaZ, gasPlumeDist)
                                             + reconstSaturation(lowerBound + (count+1)*deltaZ, gasPlumeDist))/2.0 - satW);
                }
            satIntegral = satIntegral * deltaZ;

            return satIntegral;
        }
        /*! \brief Calculates integral of difference of relative permeability over z
         *
         */
        Scalar calculateErrorKrwIntegral(Scalar lowerBound, Scalar upperBound, Scalar satW, Scalar gasPlumeDist)
        {
            int intervalNumber = 10;
            Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;

            Scalar krwIntegral = 0.0;
            for(int count=0; count<intervalNumber; count++ )
                {
                    Scalar sat1 = reconstSaturation(lowerBound + count*deltaZ, gasPlumeDist);
                    Scalar sat2 = reconstSaturation(lowerBound + (count+1)*deltaZ, gasPlumeDist);
                    Scalar krw1 = MaterialLaw::krw(this->spatialParams().materialLawParams(dummy_), sat1);
                    Scalar krw2 = MaterialLaw::krw(this->spatialParams().materialLawParams(dummy_), sat2);
                    Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(dummy_), satW);
                    krwIntegral += std::abs((krw1 + krw2)/2.0 - krw);
                }
            krwIntegral = krwIntegral * deltaZ;

            return krwIntegral;
        }

        // Calculates integral of difference of pressure over z (Zakaria 2022-08-25)
        /*! \brief Calculates integral of difference of relative Pressure over z
         *
         */
        Scalar calculateErrorPressIntegral(Scalar lowerBound, Scalar upperBound, Scalar pressW, Scalar gasPlumeDist, Scalar coarsePressW)
        {
            int intervalNumber = 10;
            Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;

            Scalar pressIntegral = 0.0;
            for(int count=0; count<intervalNumber; count++ )
                {
                    pressIntegral += std::abs((reconstPressureW(lowerBound + count*deltaZ,gasPlumeDist,coarsePressW) + reconstPressureW(lowerBound + (count+1)*deltaZ,gasPlumeDist,coarsePressW))/2.0 - pressW);
                }
            pressIntegral = pressIntegral * deltaZ;

            return pressIntegral;
        }
        ////

        // Calculates integral of difference of capillary pressure over z (Zakaria 2022-08-30)
        /*! \brief Calculates integral of difference of capillary Pressure over z
         *
         */
        Scalar calculateErrorCapPressIntegral(Scalar lowerBound, Scalar upperBound, Scalar pc, Scalar gasPlumeDist)
        {
            int intervalNumber = 10;
            Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;

            Scalar capPressIntegral = 0.0;
            for(int count=0; count<intervalNumber; count++ )
                {
                    capPressIntegral += std::abs((reconstCapPressure(lowerBound + count*deltaZ, gasPlumeDist) + reconstCapPressure(lowerBound + (count+1)*deltaZ, gasPlumeDist))/2.0 - pc);
                }
            capPressIntegral = capPressIntegral * deltaZ;

            return capPressIntegral;
        }
        ////

        /*! \brief Calculates gasPlumeDist, distance of gas plume from bottom
         *
         * Stores minGasPlumeDist for all grid cells
         */
        Scalar reconstSaturation(Scalar height, Scalar gasPlumeDist)
        {
            GlobalPosition globalPos = dummy_.geometry().center();
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar tempRef = temperatureAtPos(globalPos);
            Scalar resSatW = this->spatialParams().materialLawParams(dummy_).swr();
            Scalar resSatN = this->spatialParams().materialLawParams(dummy_).snr();
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar entryP = this->spatialParams().materialLawParams(dummy_).pe();
            Scalar lambda = this->spatialParams().materialLawParams(dummy_).lambda();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar reconstSaturation = 0.0;

            if (veModel == sharpInterface) //reconstruct phase saturation for sharp interface ve model
                {
                    reconstSaturation = 1.0;
                    if(height > gasPlumeDist)
                        {
                            reconstSaturation = resSatW;
                        }
                }
            else if (veModel == capillaryFringe) //reconstruct phase saturation for capillary fringe model
                {
                    reconstSaturation = 1.0;
                    if(height > gasPlumeDist)
                        {
                            reconstSaturation = std::pow(((height - gasPlumeDist) * (densityW - densityNw) * this->gravity().two_norm() + entryP), (-lambda))
                                * std::pow(entryP, lambda) * (1.0 - resSatW - resSatN) + resSatW;
                        }
                }

            return reconstSaturation;
        }

        // Calculation of the reconstructed water pressure (Zakaria 2022-08-25)
        /*! \brief Calculation of the water reconstructed pressure
         *
         *
         */
        Scalar reconstPressureW(Scalar height,Scalar gasPlumeDist,Scalar coarsePressW)
        {
            GlobalPosition globalPos = dummy_.geometry().center();
            Scalar pRef = referencePressureAtPos(globalPos);
            Scalar tempRef = temperatureAtPos(globalPos);
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar gravity = this->gravity().two_norm();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar reconstPressure = coarsePressW; 

            if(veModel == sharpInterface && height <= gasPlumeDist)
                {
                    reconstPressure -= densityW * gravity * height;
                }
            else if (veModel == sharpInterface && height > gasPlumeDist)
                {
                    reconstPressure -= densityW * gravity * gasPlumeDist + densityNw * gravity * (height - gasPlumeDist);
                }
            else if(veModel == capillaryFringe)
                {
                    reconstPressure -= densityW * gravity * height;
                }

            return reconstPressure;
        }
        ////

        // Calculation of the reconstructed capillary pressure (Zakaria 2022-08-30)
        /*! \brief Calculation of the water reconstructed capillary pressure
         *
         *
         */
        Scalar reconstCapPressure(Scalar height, Scalar gasPlumeDist)
        {
            GlobalPosition globalPos = dummy_.geometry().center();
            Scalar pRef = this->referencePressureAtPos(globalPos);
            Scalar tempRef = this->temperatureAtPos(globalPos);
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar gravity = this->gravity().two_norm();
            Scalar entryP = this->spatialParams().materialLawParams(dummy_).pe();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar reconstCapPressure = 0.0;

            if(veModel == 1 && height <= gasPlumeDist)
                {
                    reconstCapPressure = entryP;
                }
            else if (veModel == 1 && height > gasPlumeDist)
                {
                    reconstCapPressure = entryP + (densityW - densityNw)*gravity*(height - gasPlumeDist);
                }

            return reconstCapPressure;
        }
        ////

        /*!
         *
         * return
         */
        Scalar averageSatInPlume()
        {
            return averageSatInPlume_;
        }

        /*!
         * \name Problem parameters
         */
        // \{

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

                            values[pwIdx] = (pRef + (this->bBoxMax()[dim-1] - 0.0)
                                             * WettingPhase::density(temp, pRef)
                                             * this->gravity().two_norm());
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
                    Scalar injectionN = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, InjectionrateN);
                    Scalar injectionW = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, InjectionrateW);

                    Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();

                    //      if(time > 210.)
                    //      {
                    values[nPhaseIdx] = injectionN;
                    values[wPhaseIdx] = injectionW;
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

        std::multimap<int, Element>& getColumnMap()
        {
            return mapAllColumns_;
        }

        int getColumnModel(int columnIndex)
        {
            return modelVector_[columnIndex];
        }

        Scalar residualSegSaturation() const
        {}

        Scalar residualSegSatAverage() const
        {}

        Scalar capillaryTransitionZone() const
        {
            return CTZ_;
        }

    private:

        const Scalar eps_;
        std::string name_;
        Dumux::VtkMultiWriter<LeafGridView> visualizationWriter_;
        Dumux::VtkMultiWriter<LeafGridView> normWriter_;
        Dumux::VtkMultiWriter<LeafGridView> debugWriter_;
        std::multimap<int, Element> mapAllColumns_;
        std::ofstream outputFile_;
        int outputInterval_;
        Scalar TwoDArea_;
        Scalar averageSatInPlume_;
        Element dummy_;
        Scalar CTZ_;
        std::vector<int> modelVector_;

        // additional global variables (Zakaria 2022-08-25)
        Scalar segTime_;
        std::vector<Scalar> gasPlumeDistTemp_;
        std::chrono::_V2::system_clock::time_point beginCPU_;
        int colCounter_;
        std::vector<Scalar> criterion_;
        std::vector<Scalar> zCenter_;
        ////
    };
} //end namespace

#endif
