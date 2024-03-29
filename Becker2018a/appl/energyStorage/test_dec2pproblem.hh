/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUMUX_TEST_ENERGY_STORAGE_PROBLEM_HH
#define DUMUX_TEST_ENERGY_STORAGE_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/ch4.hh>


#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/2dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/omethod/2dpressureproperties.hh>
//#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2padaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>

#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>
#include "test_dec2pspatialparams.hh"

//#include <dumux/linear/amgbackend.hh>


namespace Dumux {

  template<class TypeTag>
  class TestEnergyStorageProblem;

  //////////
  // Specify the properties for the lens problem
  //////////
  namespace Properties {
    NEW_TYPE_TAG(TestEnergyStorageProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, TestEnergyStorageSpatialParams));
    //NEW_TYPE_TAG(TestEnergyStorageProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, TestEnergyStorageSpatialParams));
    //NEW_TYPE_TAG(TestEnergyStorageProblem, INHERITS_FROM(FvMpfaL2dPressureTwoP, FVTransportTwoP, IMPESTwoP, TestEnergyStorageSpatialParams));
    //NEW_TYPE_TAG(TestEnergyStorageProblem, INHERITS_FROM(FvMpfaO2dPressureTwoP, FVTransportTwoP, IMPESTwoP, TestEnergyStorageSpatialParams));


    // Set the grid type
    //SET_TYPE_PROP(TestEnergyStorageProblem, Grid, Dune::YaspGrid<2>);
    SET_TYPE_PROP(TestEnergyStorageProblem, Grid, Dune::UGGrid<2>);

    // use the AMG backend for parallel runs
    //SET_TYPE_PROP(TestEnergyStorageProblem, LinearSolver, AMGBackend<TypeTag>);

    // Set the problem property
    SET_TYPE_PROP(TestEnergyStorageProblem, Problem, Dumux::TestEnergyStorageProblem<TypeTag>);

    // Set the wetting phase
    SET_PROP(TestEnergyStorageProblem, WettingPhase)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    public:
      typedef FluidSystems::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
    };

    // Set the non-wetting phase
    SET_PROP(TestEnergyStorageProblem, NonwettingPhase)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    public:
      typedef FluidSystems::GasPhase<Scalar, Dumux::CH4<Scalar> > type;
      //    typedef FluidSystems::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
    };

    //SET_TYPE_PROP(TestEnergyStorageProblem, FluidSystem, BrineCO2FluidSystem<TypeTag>);

    // Enable gravity
    SET_BOOL_PROP(TestEnergyStorageProblem, ProblemEnableGravity, true);
    SET_TYPE_PROP(TestEnergyStorageProblem, EvalCflFluxFunction, Dumux::EvalCflFluxDefault<TypeTag>);
    SET_SCALAR_PROP(TestEnergyStorageProblem, ImpetCFLFactor, 0.9);
    SET_SCALAR_PROP(TestEnergyStorageProblem, ImpetIterationFlag, 0);
    SET_SCALAR_PROP(TestEnergyStorageProblem, ImpetIterationNumber, 2);
    SET_BOOL_PROP(TestEnergyStorageProblem, EnableCompressibility, false);
  }
  /*!
   * \ingroup DecoupledProblems
   */
  template<class TypeTag>
  class TestEnergyStorageProblem: public IMPESProblem2P<TypeTag> {
    typedef IMPESProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView)GridView;
    typedef typename Grid::LeafGridView LeafGridView;
    typedef typename Grid::GlobalIdSet GlobalIdSet;
    typedef typename GlobalIdSet::IdType IdType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonWettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    enum
      {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
      };

    enum
      {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pwIdx = Indices::pwIdx,
        swIdx = Indices::swIdx,
        pressEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx
      };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes)::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef std::array<int, dim> CellArray;

  public:
    TestEnergyStorageProblem(TimeManager& timeManager, const GridView &gridView) :
      ParentType(timeManager, gridView), normWriter_(GridCreator::grid().leafGridView(), "gridAfterRefinement")
    {
      int outputInterval = 0;
      if (ParameterTree::tree().hasKey("Problem.OutputInterval"))
        {
	  outputInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputInterval);
        }
      this->setOutputInterval(outputInterval);

      Scalar outputTimeInterval = 1e6;
      if (ParameterTree::tree().hasKey("Problem.OutputTimeInterval"))
        {
	  outputTimeInterval = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, OutputTimeInterval);
        }
      this->setOutputTimeInterval(outputTimeInterval);

      // store pointer to all elements in a multimap, elements belonging to the same column
      // have the same key, starting with key = 0 for first column
      // TODO: only works for equidistant grids

      // iterate over all elements
      int j = 0;
      for (const auto& element : Dune::elements(this->gridView()))
        {
	  // identify column number
	  GlobalPosition globalPos = element.geometry().center();
	  CellArray numberOfCellsX = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
	  double deltaX = this->bBoxMax()[0]/numberOfCellsX[0];

	  j = round((globalPos[0] - (deltaX/2.0))/deltaX);

	  mapColumns_.insert(std::make_pair(j, element));
	  dummy_ = element;
        }

      numberOfColumns_ = j;

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

      //calculate segregation time
      Scalar height = this->bBoxMax()[dim-1];
      Scalar porosity = this->spatialParams().porosity(dummy_);
      Scalar viscosityW = WettingPhase::viscosity(tempRef, pRef);
      Scalar permeability = this->spatialParams().intrinsicPermeability(dummy_);
      segTime_ = (height*porosity*viscosityW)/(permeability*gravity*(densityW-densityN));
      std::cout << "segTime " << segTime_ << std::endl;

      //        outputFile_.open("errorTimeRelPerm20.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat20.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeVel20.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9c.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9d.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9e.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9f.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9g.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat9h.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeRelPerm200.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat200.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorVel200.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99c.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99d.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99e.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99f.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99g.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorTimeSat99h.out", std::ios::trunc);
      //        outputFile_.close();

      //        outputFile_.open("timeAverageSat.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorSat.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorRelPerm.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorVel.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("errorNumSat.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("satProfiles.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("relPermProfiles.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("averageSatColumn.out", std::ios::trunc);
      //        outputFile_.close();
      //        outputFile_.open("averageSatPlume.out", std::ios::trunc);
      //        outputFile_.close();

      //// create files to store results (Zakaria 2022-08-12 and 2022-08-30)
      outputFile_.open("errorTimeSat20.out", std::ios::trunc);
      outputFile_.close();      
      outputFile_.open("errorTimeSatNorm20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeRelPerm20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeRelPermNorm20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimePress20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimePressNorm20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeCapPress20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeCapPressNorm20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("timeZp20.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeSat200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeSatNorm200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeRelPerm200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeRelPermNorm200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimePress200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimePressNorm200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeCapPress200.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorTimeCapPressNorm200.out", std::ios::trunc);
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
      outputFile_.open("errorSat.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorSatNorm.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorRelPerm.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorRelPermNorm.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorPress.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorPressNorm.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorCapPress.out", std::ios::trunc);
      outputFile_.close();
      outputFile_.open("errorCapPressNorm.out", std::ios::trunc);
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

      veModel_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

      // store final time to save final saturation, permeability and pressure profiles (Zakaria 2022-08-26)
      endTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, TimeManager, TEnd);
      ////

      //for adaptivity:
      GridCreator::grid().globalRefine(GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel));
      this->setGrid(GridCreator::grid());
    }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
      ParentType::init();

      // store initial CPU time (Zakaria 2022-08-12)
      beginCPU_ = std::chrono::high_resolution_clock::now();
      ////

      bool plotFluidMatrixInteractions = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotFluidMatrixInteractions);

      // plot the Pc-Sw curves, if requested
      if(plotFluidMatrixInteractions)
	this->spatialParams().plotMaterialLaw();

      // initialize the Newton first guess to calculate the gas plume height (Zakaria 2022-08-12)
      Scalar domainHeight = this->bBoxMax()[dim - 1];
      gasPlumeDistTemp_.resize(numberOfColumns_); // temporary value of gasPlumeDist to be updated at each time step
      for (int i = 0; i < numberOfColumns_; i++) 
	{
	  gasPlumeDistTemp_[i] = domainHeight/2.0; // we take half the domain height
	}
      ////

      // initialize column counter used to store solution at only one time step (Zakaria 2022-08-12)
      colCounter_ = 0;
      ////
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
    const char *name() const
    {
      if (ParameterTree::tree().hasKey("Problem.Name"))
        {
	  std::string fileName(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name));
	  return fileName.c_str();
        }
      else
        {
	  return "test_impes2p";
        }
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
      return 326.0; // -> 41,85°C
    }

    // \}

    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
      return 1.0e7; //
    }

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

    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
      values = 0;

      if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
	  if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
	      Scalar pRef = referencePressureAtPos(globalPos);
	      Scalar temp = temperatureAtPos(globalPos);

	      values[pwIdx] = (pRef + (this->bBoxMax()[dim-1] - globalPos[dim-1])
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
	  values[nPhaseIdx] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, Injectionrate);
        }
    }

    void source(PrimaryVariables &values, const Element& element) const
    {
      values = 0.0;
    }

    void initialAtPos(PrimaryVariables &values,
		      const GlobalPosition& globalPos) const
    {
      values[pressEqIdx] = 1.0e7;
      values[swIdx] = 1.0;
    }

    void postTimeStep()
    {
      ParentType::postTimeStep();

      // store CPU times (Zakaria 2022-08-12)
      auto endCPU = std::chrono::high_resolution_clock::now(); // time of end of simulation
      auto elapsedCPU = std::chrono::duration_cast<std::chrono::nanoseconds>(endCPU - beginCPU_); // duration of simulation
      auto beginPTS = std::chrono::high_resolution_clock::now(); // start of postTimeStep
      ////

      //         //refine grid for error norm calculation
      //         LeafGridView gridView = this->gridView();
      //         std::unordered_map<IdType, int> mapGlobalIdx;
      //         GridCreator::grid().preAdapt();
      //         const int additionalRefinementSteps = 2;
      //         for(int i=0;i<additionalRefinementSteps;i++)
      //         {
      //             for (const auto& element : Dune::elements(gridView))
      //             {
      //                 if(i == 0)
      //                 {
      //                     IdType globalId = gridView.grid().globalIdSet().id(element);
      //                     int eIdxGlobal = this->variables().index(element);
      //                     mapGlobalIdx.insert(std::pair<IdType,int>(globalId, eIdxGlobal));
      //                 }
      //                 GlobalPosition globalPos = element.geometry().center();
      //                 GridCreator::grid().mark(1, element);
      //             }
      //             // adapt the grid
      //             GridCreator::grid().adapt();
      //             GridCreator::grid().postAdapt();
      //         }
      //
      //         // use error norm writer
      //         normWriter_.gridChanged();
      //         normWriter_.beginWrite(this->timeManager().time() + this->timeManager().timeStepSize());
      //
      //         typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::ScalarSolution ScalarSolutionType;
      //         typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;
      //
      //         int size = gridView.size(0);
      //         ScalarSolutionType *pressureW1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *pressureNw1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *potentialW1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *potentialNw1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *capillaryPressure1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *saturationW1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *saturationNw1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *mobilityW1 = normWriter_.allocateManagedBuffer (size);
      //         ScalarSolutionType *mobilityNw1 = normWriter_.allocateManagedBuffer (size);
      //
      //         //write output on refined grid
      //         for (const auto& element : Dune::elements(gridView))
      //         {
      //             //identify column number of element and check if it is a VE column or a full-d column
      //             GlobalPosition globalPos = element.geometry().center();
      //
      //             //identify ID of cell, same as ID of parallel run with 20 blocks
      //             //determine block number starting from 0
      //             int blockNumberX = floor(globalPos[0]/35);
      //             int blockNumberY = floor(globalPos[dim-1]/15);
      //             int blockNumber = blockNumberY*10 + blockNumberX;
      //             //determine local cell ID in block, starting from 0
      //             Scalar localX = globalPos[0] - (blockNumberX*35);
      //             Scalar localY = globalPos[dim-1] - (blockNumberY*15);
      //             int localNumberX = floor(localX/0.25);
      //             int localNumberY = floor(localY/(30.0/512.0));
      //             int localID = localNumberY*140 + localNumberX;
      //             //determine global ID, starting from 0
      //             int globalID = blockNumber*35840 + localID;
      //
      //             //full-d elements
      //             IdType globalId = gridView.grid().globalIdSet().id(element.father().father());
      //             int eIdxGlobal2D = mapGlobalIdx.find(globalId)->second;
      //             CellData& cellData = this->variables().cellData(eIdxGlobal2D);
      //             (*pressureW1)[globalID]  = cellData.pressure(wPhaseIdx);
      //             (*pressureNw1)[globalID] = cellData.pressure(nPhaseIdx);
      //             (*potentialW1)[globalID]  = cellData.potential(wPhaseIdx);
      //             (*potentialNw1)[globalID] = cellData.potential(nPhaseIdx);
      //             (*capillaryPressure1)[globalID] = cellData.capillaryPressure();
      //             (*saturationW1)[globalID] = cellData.saturation(wPhaseIdx);
      //             (*saturationNw1)[globalID] = cellData.saturation(nPhaseIdx);
      //             (*mobilityW1)[globalID] = cellData.mobility(wPhaseIdx);
      //             (*mobilityNw1)[globalID] = cellData.mobility(nPhaseIdx);
      //         }
      //         normWriter_.attachCellData(*pressureW1, "wetting pressure");
      //         normWriter_.attachCellData(*pressureNw1, "non-wetting pressure");
      //         normWriter_.attachCellData(*potentialW1, "wetting potential");
      //         normWriter_.attachCellData(*potentialNw1, "non-wetting potential");
      //         normWriter_.attachCellData(*capillaryPressure1, "capillary pressure");
      //         normWriter_.attachCellData(*saturationW1, "wetting saturation");
      //         normWriter_.attachCellData(*saturationNw1, "non-wetting saturation");
      //         normWriter_.attachCellData(*mobilityW1, "wetting mobility");
      //         normWriter_.attachCellData(*mobilityNw1, "non-wetting mobility");
      //         normWriter_.endWrite();
      //
      //         //coarsen grid after error norm calculation
      //         for(int i=0;i<additionalRefinementSteps;i++)
      //         {
      //             for (const auto& element : Dune::elements(gridView))
      //             {
      //                 GlobalPosition globalPos = element.geometry().center();
      //                 GridCreator::grid().mark(-1, element);
      //             }
      //             // adapt the grid
      //             GridCreator::grid().adapt();
      //             GridCreator::grid().postAdapt();
      //             GridCreator::grid().preAdapt();
      //         }
      //

      // store time in result files (Zakarai 2022-08-12 and 2022-08-30)
      outputFile_.open("errorTimeSat20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeSatNorm20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeRelPerm20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeRelPermNorm20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimePress20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimePressNorm20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeCapPress20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeCapPressNorm20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("timeZp20.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeSat200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeSatNorm200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeRelPerm200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeRelPermNorm200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimePress200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimePressNorm200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeCapPress200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorTimeCapPressNorm200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("timeZp200.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorSat.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorSatNorm.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorRelPerm.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorRelPermNorm.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorPress.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorPressNorm.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorCapPress.out", std::ios::app);
      outputFile_ << this->timeManager().time()/segTime_;
      outputFile_.close();
      outputFile_.open("errorCapPressNorm.out", std::ios::app);
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

      Scalar averageSatColumn[numberOfColumns_];
      //      Scalar averageSatPlume[numberOfColumns_];
      Scalar errorSat[numberOfColumns_];
      //      Scalar errorNumSat[numberOfColumns_];
      Scalar errorRelPerm[numberOfColumns_];
      //      Scalar velocityWVerticalMax[numberOfColumns_];
      //      Scalar velocityWVerticalMaxTotal;
      
      // initialize error arrays for pressure criteria - only used for capillary fringe model - and normalized saturation and relative permeability criteria (Zakaria 2022-08-12 and 2022-08-30)
      Scalar errorSatNorm[numberOfColumns_];
      Scalar errorRelPermNorm[numberOfColumns_];
      Scalar errorPress[numberOfColumns_];
      Scalar errorPressNorm[numberOfColumns_];
      Scalar errorCapPress[numberOfColumns_];
      Scalar errorCapPressNorm[numberOfColumns_];
      ////

      //initialize/reset average column saturation
      for (int i = 0; i < numberOfColumns_; i++)
	{
	  averageSatColumn[i] = 0.0;
	  //            averageSatPlume[i] = 0.0;
	  errorSat[i] = 0.0;
	  //            errorNumSat[i] = 0.0;
	  errorRelPerm[i] = 0.0;
	  //             velocityWVerticalMax[i] = 0.0;
	  //             velocityWVerticalMaxTotal = 0.0;

	  // fill additional error arrays with zeros (Zakaria 2022-08-12 and 2022-08-30)
	  errorSatNorm[i] = 0.0;
	  errorRelPermNorm[i] = 0.0;
	  errorPress[i] = 0.0;
	  errorPressNorm[i] = 0.0;
	  errorCapPress[i] = 0.0;
	  errorCapPressNorm[i] = 0.0;
	  /////
	}
      Scalar averageSatTotal = 0.0;
      Scalar gasPlumeVolume = 0.0;// volume of gas plume
      for (int i = 0; i != numberOfColumns_; ++i)
	{
	  Scalar totalVolume = 0.0;// total volume of column
	  Scalar gasVolume = 0.0;// volume with SatN>0.0 in one column;
	  typename std::map<int, Element>::iterator it = mapColumns_.lower_bound(i);
	  for (; it != mapColumns_.upper_bound(i); ++it)
	    {
	      int globalIdxI = this->variables().index(it->second);
	      Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
	      Scalar satNw = this->variables().cellData(globalIdxI).saturation(nPhaseIdx);
	      Scalar volume = it->second.geometry().volume();
	      averageSatColumn[i] += satW * volume;
	      totalVolume += volume;
	      if(satNw>0.0)
		{
		  //		  averageSatPlume[i] += satW * volume;
		  gasVolume += volume;
		  averageSatTotal += satW * volume;
		  gasPlumeVolume += volume;
		}
	    }
	  averageSatColumn[i] = averageSatColumn[i]/totalVolume;//average wetting saturation in column (equals gasPlumeDist for SI and no compressibility)

	  // to make sure the average saturation does not exceed 1, otherwise Newton method may not converge to compute the gas plume thickness (Zakaria 2022-09-12)
	  averageSatColumn[i] = std::min(averageSatColumn[i],1.);
	  ////

	  //averageSatPlume[i] = averageSatPlume[i]/gasVolume;//average wetting saturation in plume

	  //            outputFile_.open("averageSatColumn.out", std::ios::app);
	  //            outputFile_ << " " << averageSatColumn[i];
	  //            outputFile_.close();

	  //            outputFile_.open("averageSatPlume.out", std::ios::app);
	  //            outputFile_ << " " << averageSatPlume[i];
	  //            outputFile_.close();

	  Scalar domainHeight = this->bBoxMax()[dim - 1];
	  CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
	  double deltaZ = domainHeight/numberOfCells[dim - 1];
	  Scalar resSatW = this->spatialParams().materialLawParams(it->second).swr();

	  // set the Newton initial guess to calculate plume height (Zakaria 2022-08-12)
	  Scalar gasPlumeDist = calculateGasPlumeDist(dummy_,averageSatColumn[i],gasPlumeDistTemp_[i]);
	  //gasPlumeDistTemp_[i] = std::min(0.99*domainHeight,gasPlumeDist); // update new temporary height
	  gasPlumeDistTemp_[i] = domainHeight/2.0;
	  ////

	  //calculate error to VE situation
	  it = mapColumns_.lower_bound(i);

	  // store the bottom element of column i to compute coarse pressure later (Zakaria 2022-08-12 and 2022-08-30)	  
	  Element veElement = mapColumns_.lower_bound(i)->second;
	  int eIdxGlobal = this->variables().index(veElement);
	  Scalar coarsePressW = this->variables().cellData(eIdxGlobal).pressure(wPhaseIdx);
	  Scalar bottomPressW = reconstPressureW(-deltaZ/2,gasPlumeDist,coarsePressW); //slightly different from coarsePressW when columnModel > 1 since coarsePressW is in this case the pressure at the center of the bottom cell and not the pressure at bottom
	  ////

	  for (; it != mapColumns_.upper_bound(i); ++it)
	    {
	      Element element = it->second;
	      int globalIdxI = this->variables().index(element);
	      GlobalPosition globalPos = element.geometry().center();
	      Scalar top = globalPos[dim - 1] + deltaZ/2.0;
	      Scalar bottom = globalPos[dim - 1] - deltaZ/2.0;

	      Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
	      Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(element), satW);

	      // pressures in cell (Zakaria 2022-08-12 and 2022-08-30)
	      Scalar pressW = this->variables().cellData(globalIdxI).pressure(wPhaseIdx);
	      Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParams(element), satW);
	      Scalar pe = this->spatialParams().materialLawParams(element).pe();
	      ////

	      if(veModel_ == 0)//calculate error for VE model
		{
		  if (top <= gasPlumeDist)
		    {
		      errorSat[i] += std::abs(deltaZ * (satW - 1.0));
		      //errorNumSat[i] += std::abs(deltaZ * (satW - 1.0));
		      errorRelPerm[i] += std::abs(deltaZ * (krw - 1.0));

		      // normalized errors (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* std::abs(deltaZ * (satW - 1.0));
		      errorRelPermNorm[i] += (1/krw)* std::abs(deltaZ * (krw - 1.0));
		      ////
		    }
		  else if (bottom >= gasPlumeDist)
		    {
		      errorSat[i] += std::abs(deltaZ * (satW - resSatW));
		      //errorNumSat[i] += std::abs(deltaZ * (satW - resSatW));
		      errorRelPerm[i] += std::abs(deltaZ * (krw - 0.0));

		      // normalized errors (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* std::abs(deltaZ * (satW - resSatW));
		      errorRelPermNorm[i] += (1/krw)* std::abs(deltaZ * (krw - 0.0));
		      ////
		    }
		  else
		    {
		      Scalar lowerDelta = gasPlumeDist - bottom;
		      Scalar upperDelta = top - gasPlumeDist;
		      errorSat[i] += std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW));
		      //errorNumSat[i] += std::abs(deltaZ * satW - lowerDelta);
		      errorRelPerm[i] += std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0));

		      // normalized errors (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW)));
		      errorRelPermNorm[i] += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0)));
		      ////
		    }
		}
	      if(veModel_ == 1)//calculate error for capillary fringe model///////////////////////////////////////////////////////////////////////////////////
		{
		  if (top <= gasPlumeDist)
		    {
		      errorSat[i] += std::abs(deltaZ * (satW - 1.0));
		      //errorNumSat[i] += std::abs(deltaZ * (satW - 1.0));
		      errorRelPerm[i] += std::abs(deltaZ * (krw - 1.0));

		      // normalized errors (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* ( std::abs(deltaZ * (satW - 1.0)) ) ;
		      errorRelPermNorm[i] += (1/krw)* ( std::abs(deltaZ * (krw - 1.0)) ) ;
		      ////

		      // capillary pressure errors (Zakaria 2022-08-30)
		      errorCapPress[i] += std::abs(deltaZ * (pc - pe));
		      errorCapPressNorm[i] += (1/pc)* std::abs(deltaZ * (pc - pe));
		      ////
		    }
		  else if (bottom >= gasPlumeDist)
		    {
		      Scalar eIntSat = calculateErrorSatIntegral(bottom, top, satW, gasPlumeDist);
		      errorSat[i] += eIntSat;
		      //errorNumSat[i] += std::abs(deltaZ * satW - eIntSat);

		      Scalar eIntPerm = calculateErrorKrwIntegral(bottom, top, satW, gasPlumeDist);
		      errorRelPerm[i] += eIntPerm;
		      
		      // errors in normalized variables (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* eIntSat;
		      errorRelPermNorm[i] += (1/krw)* eIntPerm;
		      ////

		      // capillary pressure errors (Zakaria 2022-08-30)
		      Scalar eIntPc = calculateErrorCapPressIntegral(bottom, top, pc, gasPlumeDist);
		      errorCapPress[i] += eIntPc;
		      errorCapPressNorm[i] += (1/pc)* eIntPc;
		      ////
		    }
		  else
		    {
		      Scalar lowerDelta = gasPlumeDist - bottom;
		      Scalar upperDelta = top - gasPlumeDist;

		      Scalar eIntSat = calculateErrorSatIntegral(gasPlumeDist, top, satW, gasPlumeDist);
		      errorSat[i] += std::abs(lowerDelta * (satW - 1.0)) + eIntSat;
		      //errorNumSat[i] += std::abs(deltaZ * satW - (eIntSat + lowerDelta));

		      Scalar eIntPerm = calculateErrorKrwIntegral(gasPlumeDist, top, satW, gasPlumeDist);
		      errorRelPerm[i] += std::abs(lowerDelta * (krw - 1.0)) + eIntPerm;

		      // error in normalized variables (Zakaria 2022-08-12)
		      errorSatNorm[i] += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + eIntSat);
		      errorRelPermNorm[i] += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + eIntPerm);
		      ////

		      // capillary pressure errors (Zakaria 2022-08-30)
		      Scalar eIntPc = calculateErrorCapPressIntegral(gasPlumeDist, top, pc, gasPlumeDist);
		      errorCapPress[i] += std::abs(lowerDelta * (pc - pe)) + eIntPc;
		      errorCapPressNorm[i] += (1/pc)* (std::abs(lowerDelta * (pc - pe)) + eIntPc);
		      ////
		    }
		}
	      // error in water pressure (Zakaria 2022-08-26)
	      Scalar eIntPress = calculateErrorPressIntegral(bottom, top, pressW, gasPlumeDist, bottomPressW);
	      errorPress[i] += eIntPress;
	      errorPressNorm[i] += (1/pressW)* eIntPress;
	      ////

	      //                 //calculate velocity criterion
	      //                 for (const auto& intersection : intersections(this->gridView(), element))
	      //                 {
	      //                     int isIndex = intersection.indexInInside();
	      //                     GlobalPosition velocityW = this->variables().cellData(globalIdxI).fluxData().velocity(wPhaseIdx, isIndex);
	      //                     GlobalPosition gravityNormalized = this->gravity();
	      //                     gravityNormalized /= this->gravity().two_norm();
	      //                     velocityW = velocityW * gravityNormalized;
	      //                     velocityWVerticalMax[i] = std::max(velocityWVerticalMax[i], std::abs(velocityW.two_norm()));
	      //                     velocityWVerticalMaxTotal = std::max(velocityWVerticalMaxTotal, velocityWVerticalMax[i]);
	      //                 }
	    }
	  //
	  //            errorSat[i] = errorSat[i]/(domainHeight*deltaZ);
	  //            errorNumSat[i] = errorNumSat[i]/domainHeight;
	  errorSat[i] = errorSat[i]/(domainHeight-gasPlumeDist);
	  //            errorNumSat[i] = errorNumSat[i]/(domainHeight-gasPlumeDist);
	  errorRelPerm[i] = errorRelPerm[i]/(domainHeight-gasPlumeDist);
	  //errorSat[i] = errorSat[i]/satW;
	  //errorRelPerm[i] = errorRelPerm[i]/krw;

	  // update of pressure and normalized errors (Zakaria 2022-08-12 and 2022-08-30)
	  errorSatNorm[i] = errorSatNorm[i]/(domainHeight-gasPlumeDist);
	  errorRelPermNorm[i] = errorRelPermNorm[i]/(domainHeight-gasPlumeDist);
	  errorPress[i] = errorPress[i]/(domainHeight-gasPlumeDist);
	  errorPressNorm[i] = errorPressNorm[i]/(domainHeight-gasPlumeDist);
	  errorCapPress[i] = errorCapPress[i]/(domainHeight-gasPlumeDist);
	  errorCapPressNorm[i] = errorCapPressNorm[i]/(domainHeight-gasPlumeDist);
	  ////

	  if(averageSatColumn[i]>1.0-eps_)
	    {
	      errorSat[i] = 0.0;
	      errorRelPerm[i] = 0.0;
	      //errorNumSat[i] = 0.0;
	      //velocityWVerticalMax[i] = 0.0;
	      
	      // set pressure and norm errors to zero if avg saturation is 1 (Zakaria 2022-08-12 and 2022-08-30)
	      errorSatNorm[i] = 0.0;
	      errorRelPermNorm[i] = 0.0;
	      errorPress[i] = 0.0;
	      errorPressNorm[i] = 0.0;
	      errorCapPress[i] = 0.0;
	      errorCapPressNorm[i] = 0.0;
	      ////
	    }

	  if(i == 19)
	    {
	      //                Scalar errorSatb = errorRelPerm[i]/((domainHeight-gasPlumeDist)*deltaZ/CTZ_);
	      //                Scalar errorSatc = errorRelPerm[i]/(domainHeight);
	      //                Scalar errorSatd = errorRelPerm[i]/(domainHeight*deltaZ/CTZ_);
	      //                Scalar errorSatf = errorSat[i]/((domainHeight-gasPlumeDist)*deltaZ/CTZ_);
	      //                Scalar errorSatg = errorSat[i]/(domainHeight);
	      //                Scalar errorSath = errorSat[i]/(domainHeight*deltaZ/CTZ_);

	      //                outputFile_.open("errorTimeRelPerm20.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << errorRelPerm[i] << std::endl;
	      //                outputFile_.close();
	      //                outputFile_.open("errorTimeSat9b.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatb << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat9c.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatc << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat9d.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatd << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSat << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat9f.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatf << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat9g.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatg << std::endl;
	      //                outputFile_.close();

	      //                outputFile_.open("errorTimeSat20.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << errorSat[i] << std::endl;
	      //                outputFile_.close();

	      //                outputFile_.open("errorTimeVel20.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << velocityWVerticalMax[i] << std::endl;
	      //                outputFile_.close();

	      // store errors and plume height in files for 20th column (Zakaria 2022-08-12 and 2022-08-30)
	      outputFile_.open("errorTimeSat20.out", std::ios::app);
	      outputFile_ << " " << errorSat[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeSatNorm20.out", std::ios::app);
	      outputFile_ << " " << errorSatNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeRelPerm20.out", std::ios::app);
	      outputFile_ << " " << errorRelPerm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeRelPermNorm20.out", std::ios::app);
	      outputFile_ << " " << errorRelPermNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimePress20.out", std::ios::app);
	      outputFile_ << " " << errorPress[i]<< std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimePressNorm20.out", std::ios::app);
	      outputFile_ << " " << errorPressNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeCapPress20.out", std::ios::app);
	      outputFile_ << " " << errorCapPress[i]<< std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeCapPressNorm20.out", std::ios::app);
	      outputFile_ << " " << errorCapPressNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("timeZp20.out", std::ios::app);
	      outputFile_ << " " << gasPlumeDist << std::endl;
	      outputFile_.close();
	      ////
	    }

	  if(i == 199)
	    {
	      //                Scalar errorSatb = errorRelPerm[i]/((domainHeight-gasPlumeDist)*deltaZ/CTZ_);
	      //                Scalar errorSatc = errorRelPerm[i]/(domainHeight);
	      //                Scalar errorSatd = errorRelPerm[i]/(domainHeight*deltaZ/CTZ_);
	      //                Scalar errorSatf = errorSat[i]/((domainHeight-gasPlumeDist)*deltaZ/CTZ_);
	      //                Scalar errorSatg = errorSat[i]/(domainHeight);
	      //                Scalar errorSath = errorSat[i]/(domainHeight*deltaZ/CTZ_);

	      //                outputFile_.open("errorTimeRelPerm200.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << errorRelPerm[i] << std::endl;
	      //                outputFile_.close();

	      //                outputFile_.open("errorTimeSat99b.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatb << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat99c.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatc << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat99d.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatd << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat99e.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSate << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat99f.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatf << std::endl;
	      //                outputFile_.close();
	      //
	      //                outputFile_.open("errorTimeSat99g.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time() << " " << errorSatg << std::endl;
	      //                outputFile_.close();

	      //                outputFile_.open("errorTimeSat200.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << errorSat[i] << std::endl;
	      //                outputFile_.close();

	      //                outputFile_.open("errorTimeVel200.out", std::ios::app);
	      //                outputFile_ << this->timeManager().time()/segTime_ << " " << velocityWVerticalMax[i] << std::endl;
	      //                outputFile_.close();

	      // store errors and plume height in files for 200th column (Zakaria 2022-08-12 and 2022-08-30)
	      outputFile_.open("errorTimeSat200.out", std::ios::app);
	      outputFile_ << " " << errorSat[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeSatNorm200.out", std::ios::app);
	      outputFile_ << " " << errorSatNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeRelPerm200.out", std::ios::app);
	      outputFile_ << " " << errorRelPerm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeRelPermNorm200.out", std::ios::app);
	      outputFile_ << " " << errorRelPermNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimePress200.out", std::ios::app);
	      outputFile_ << " " << errorPress[i]<< std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimePressNorm200.out", std::ios::app);
	      outputFile_ << " " << errorPressNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeCapPress200.out", std::ios::app);
	      outputFile_ << " " << errorCapPress[i]<< std::endl;
	      outputFile_.close();
	      outputFile_.open("errorTimeCapPressNorm200.out", std::ios::app);
	      outputFile_ << " " << errorCapPressNorm[i] << std::endl;
	      outputFile_.close();
	      outputFile_.open("timeZp200.out", std::ios::app);
	      outputFile_ << " " << gasPlumeDist << std::endl;
	      outputFile_.close();
	      ////
	    }

	  //            outputFile_.open("errorSat.out", std::ios::app);
	  //            outputFile_ << " " << errorSat[i];
	  //            outputFile_.close();
	  //            outputFile_.open("errorNumSat.out", std::ios::app);
	  //            outputFile_ << " " << errorNumSat[i];
	  //            outputFile_.close();
	  //            outputFile_.open("errorVel.out", std::ios::app);
	  //            outputFile_ << " " << velocityWVerticalMax[i];
	  //            outputFile_.close();

	  // store errors and plume height in files for all columns (Zakaria 2022-08-12 and 2022-08-30)
	  outputFile_.open("errorSat.out", std::ios::app);
	  outputFile_ << " " << errorSat[i];
	  outputFile_.close();
	  outputFile_.open("errorSatNorm.out", std::ios::app);
	  outputFile_ << " " << errorSatNorm[i];
	  outputFile_.close();
	  outputFile_.open("errorRelPerm.out", std::ios::app);
	  outputFile_ << " " << errorRelPerm[i];
	  outputFile_.close();
	  outputFile_.open("errorRelPermNorm.out", std::ios::app);
	  outputFile_ << " " << errorRelPermNorm[i];
	  outputFile_.close();
	  outputFile_.open("errorPress.out", std::ios::app);
	  outputFile_ << " " << errorPress[i];
	  outputFile_.close();
	  outputFile_.open("errorPressNorm.out", std::ios::app);
	  outputFile_ << " " << errorPressNorm[i];
	  outputFile_.close();
	  outputFile_.open("errorCapPress.out", std::ios::app);
	  outputFile_ << " " << errorCapPress[i];
	  outputFile_.close();
	  outputFile_.open("errorCapPressNorm.out", std::ios::app);
	  outputFile_ << " " << errorCapPressNorm[i];
	  outputFile_.close();
	  outputFile_.open("zp.out", std::ios::app);
	  outputFile_ << " " << gasPlumeDist ;
	  outputFile_.close();
	  ////

	  //           //plot profiles
	  //            if((this->timeManager().timeStepIndex() == 2500 && i == 2) || (this->timeManager().timeStepIndex() == 2500 && i == 10) || (this->timeManager().timeStepIndex() == 2500 && i == 25) || (this->timeManager().timeStepIndex() == 2500 && i == 60))
	  //            {
	  //                //iterate over cells in column and write z-location
	  //                typename std::map<int, Element>::iterator it = mapColumns_.lower_bound(i);
	  //                for (; it != mapColumns_.upper_bound(i); ++it)
	  //                {
	  //                    GlobalPosition globalPos = (it->second).geometry().center();
	  //                    Scalar z1 = globalPos[dim - 1];
	  //                    outputFile_.open("satProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << z1;
	  //                    outputFile_.close();
	  //                    outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << z1;
	  //                    outputFile_.close();
	  //                }
	  //
	  //                //next line
	  //                outputFile_.open("satProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //                outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //
	  //                //iterate over cells in column and write satW and krw
	  //                typename std::map<int, Element>::iterator it2 = mapColumns_.lower_bound(i);
	  //                for (; it2 != mapColumns_.upper_bound(i); ++it2)
	  //                {
	  //                    int globalIdxI = this->variables().index(it2->second);
	  //                    Scalar satW = this->variables().cellData(globalIdxI).saturation(wPhaseIdx);
	  //                    Scalar krw = MaterialLaw::krw(this->spatialParams().materialLawParams(it2->second), satW);
	  //
	  //                    outputFile_.open("satProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << satW;
	  //                    outputFile_.close();
	  //                    outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << krw;
	  //                    outputFile_.close();
	  //                }
	  //
	  //                //next line
	  //                outputFile_.open("satProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //                outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //
	  //                //iterate over height and write z-location
	  //                int steps = 100;
	  //                Scalar deltaZ = domainHeight/steps;
	  //                for (int j=0; j <= steps; ++j)
	  //                {
	  //                    Scalar z2 = j*deltaZ;
	  //
	  //                    outputFile_.open("satProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << z2;
	  //                    outputFile_.close();
	  //                    outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << z2;
	  //                    outputFile_.close();
	  //                }
	  //
	  //                //next line
	  //                outputFile_.open("satProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //                outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //
	  ////                Scalar gasPlumeDist = calculateGasPlumeDist(it->second, averageSatColumn[i]);
	  //                //iterate over height and write reconstructed saturation
	  //                for (int j=0; j <= steps; ++j)
	  //                {
	  //                    Scalar z2 = j*deltaZ;
	  //                    Scalar reconstSat = reconstSaturation(z2, gasPlumeDist);
	  //                    Scalar reconstKrw = MaterialLaw::krw(this->spatialParams().materialLawParams(it2->second), reconstSat);
	  //
	  //                    outputFile_.open("satProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << reconstSat;
	  //                    outputFile_.close();
	  //                    outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                    outputFile_ << " " << reconstKrw;
	  //                    outputFile_.close();
	  //                }
	  //
	  //                //next line
	  //                outputFile_.open("satProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //                outputFile_.open("relPermProfiles.out", std::ios::app);
	  //                outputFile_ << " " << std::endl;
	  //                outputFile_.close();
	  //            }
	  //   }

	  // plot profiles at final time (Zakaria 2022-08-12 and 2022-08-30)
	  if (this->timeManager().time() + this->timeManager().timeStepSize() >= endTime_ && colCounter_ < numberOfColumns_)
	    {
	      colCounter_+=1;
	      if (i==0)
		{ 
		  // store z-coords of cell centers for first column (for rectangular meshes only)
		  typename std::map<int, Element>::iterator it1 = mapColumns_.lower_bound(i);
		  for (; it1 != mapColumns_.upper_bound(i); ++it1)
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
		    }
		}

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

	      // store results
	      typename std::map<int, Element>::iterator it2 = mapColumns_.lower_bound(i);
	      for (; it2 != mapColumns_.upper_bound(i); ++it2)
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
	  ////

	}

      //        averageSatTotal = averageSatTotal/gasPlumeVolume;

      //        outputFile_.open("timeAverageSat.out", std::ios::app);
      //        outputFile_ << this->timeManager().time() << " " << averageSatTotal << std::endl;
      //        outputFile_.close();
      //
      //        outputFile_.open("errorSat.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      //
      //        outputFile_.open("errorNumSat.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      //
      //        outputFile_.open("errorRelPerm.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      //        
      //        outputFile_.open("errorVel.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      // 
      //        outputFile_.open("averageSatColumn.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      // 
      //        outputFile_.open("averageSatPlume.out", std::ios::app);
      //        outputFile_ << " " << std::endl;
      //        outputFile_.close();
      //
      //       //check mass conservativity:
      //       GridView GridView = this->gridView();
      //       Scalar totalMassN = 0;
      //       for (const auto& element : Dune::elements(GridView))
      //       {
      //           int eIdxGlobal = this->variables().index(element);
      //           GlobalPosition globalPos = element.geometry().center();
      //           Scalar pRef = referencePressureAtPos(globalPos);
      //           Scalar temp = temperatureAtPos(globalPos);
      //           Scalar massN = this->variables().cellData(eIdxGlobal).saturation(nPhaseIdx) * element.geometry().volume() * this->spatialParams().porosity(element)
      //                        * NonWettingPhase::density(temp, pRef);
      //           totalMassN += massN;
      //       }
      //       Scalar totalMassNExpected = -GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, double, BoundaryConditions, Injectionrate) * (this->timeManager().time()+this->timeManager().timeStepSize())
      //               * this->bBoxMax()[1];
      //       std::cout << "Error in mass: " << totalMassNExpected - totalMassN << ". ";


      // finalize result files (Zakaria 2022-08-12 and 2022-08-30)
      outputFile_.open("errorSat.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorSatNorm.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorRelPerm.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorRelPermNorm.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorPress.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorPressNorm.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorCapPress.out", std::ios::app);
      outputFile_ << std::endl;
      outputFile_.close();
      outputFile_.open("errorCapPressNorm.out", std::ios::app);
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
      ////

      // postTimeStep times (Zakaria 2022-08-12)
      auto endPTS = std::chrono::high_resolution_clock::now(); // end time of postTimeStep
      auto elapsedPTS = std::chrono::duration_cast<std::chrono::nanoseconds>(endPTS - beginPTS); // duration of postTimeStep
      ////

      // store duration of simulation and postTimeStep (Zakaria 2022-08-12)
      outputFile_.open("CPUTime.out", std::ios::app);
      outputFile_ << " " << elapsedCPU.count() << " " << elapsedPTS.count()<< std::endl;
      outputFile_.close();
      ////
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

      Scalar gasPlumeDist = 0.0;

      // store number of iterations in Newton method (Zakaria 2022-08-12)
      int iterNumber = 0;
      ////

      // start time of Newton method (Zakaria 2022-08-12)
      auto beginNewton = std::chrono::high_resolution_clock::now();
      ////

      if (veModel_ == 0) //calculate gasPlumeDist for sharp interface ve model
        {
	  gasPlumeDist = domainHeight * (satW - resSatW) / (1.0 - resSatW);
        }

      else if (veModel_ == 1) //calculate gasPlumeDist for capillary fringe model
        {
	  GlobalPosition globalPos = element.geometry().center();
	  Scalar pRef = referencePressureAtPos(globalPos);
	  Scalar tempRef = temperatureAtPos(globalPos);
	  Scalar densityW = WettingPhase::density(tempRef, pRef);
	  Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
	  Scalar lambda = this->spatialParams().materialLawParams(element).lambda();
	  Scalar entryP = this->spatialParams().materialLawParams(element).pe();

	  // start from previous plume height (Zakaria 2022-08-12)
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
		  // update number of Newton iterations (Zakaria 2022-08-12)
		  iterNumber += 1;
		  ////

		  Scalar residual = 1.0 / (1.0 - lambda) * std::pow(((domainHeight - Xi) * (densityW - densityNw) * gravity + entryP),(1.0 - lambda)) * (1.0 - resSatW - resSatN) 
		    * std::pow(entryP, lambda) / ((densityW - densityNw) * gravity) + resSatW * (domainHeight - Xi) - entryP * (1.0 - resSatW - resSatN) / ((1.0 - lambda) * (densityW - densityNw) * gravity) 
		    + Xi - satW * domainHeight;

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
		  // update number of Newton iterations (Zakaria 2022-08-12)
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
                }
            }
	  //GasPlumeDist=0
	  else
            {
	      Xi = 0.0;
            }
	  gasPlumeDist = Xi;
	}

      // store end time for Newton and duration for Newton (Zakaria 2022-08-12)
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

    //    /*! \brief Calculates integral of wetting saturation over z
    //     *
    //     */
    //    Scalar calculateSatIntegral(Scalar lowerBound, Scalar upperBound, Scalar gasPlumeDist)
    //    {
    //        int intervalNumber = 10;
    //        Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;
    //
    //        Scalar satIntegral = 0.0;
    //        for(int count=0; count<intervalNumber; count++ )
    //        {
    //            satIntegral += std::abs((reconstSaturation(lowerBound + count*deltaZ, gasPlumeDist)
    //                    + reconstSaturation(lowerBound + (count+1)*deltaZ, gasPlumeDist))/2.0);
    //        }
    //        satIntegral = satIntegral * deltaZ;
    //
    //        return satIntegral;
    //    }

    /*! \brief Calculates integral of difference of relative wetting permeability over z
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

    // Calculates integral of difference of pressure over z (Zakaria 2022-08-12)
    /*! \brief Calculates integral of difference of relative Pressure over z
     *
     */
    Scalar calculateErrorPressIntegral(Scalar lowerBound, Scalar upperBound, Scalar pressW, Scalar gasPlumeDist,Scalar coarsePressW)
    {
      int intervalNumber = 10;
      Scalar deltaZ = (upperBound - lowerBound)/intervalNumber;

      Scalar pressIntegral = 0.0;
      for(int count=0; count<intervalNumber; count++ )
        {
	  pressIntegral += std::abs((reconstPressureW(lowerBound + count*deltaZ, gasPlumeDist, coarsePressW) + reconstPressureW(lowerBound + (count+1)*deltaZ, gasPlumeDist, coarsePressW))/2.0 - pressW);

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
      Scalar domainHeight = this->bBoxMax()[dim - 1];
      GlobalPosition globalPos = dummy_.geometry().center();
      Scalar pRef = referencePressureAtPos(globalPos);
      Scalar tempRef = temperatureAtPos(globalPos);
      Scalar resSatW = this->spatialParams().materialLawParams(dummy_).swr();
      Scalar resSatN = this->spatialParams().materialLawParams(dummy_).snr();
      Scalar densityW = WettingPhase::density(tempRef, pRef);
      Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
      Scalar entryP = this->spatialParams().materialLawParams(dummy_).pe();
      Scalar lambda = this->spatialParams().materialLawParams(dummy_).lambda();

      Scalar reconstSaturation = 0.0;

      if (veModel_ == 0) //reconstruct phase saturation for sharp interface ve model
        {
	  reconstSaturation = 1.0;
	  if(height > gasPlumeDist)
            {
	      reconstSaturation = resSatW;
            }
        }
      else if (veModel_ == 1) //reconstruct phase saturation for capillary fringe model
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

    // Calculation of the reconstructed water pressure (Zakaria 2022-08-12)
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

      Scalar reconstPressure = coarsePressW;

      if(veModel_ == 0 && height <= gasPlumeDist)
        {
	  reconstPressure -= densityW * gravity * height;
        }
      else if (veModel_ == 0 && height > gasPlumeDist)
        {
	  reconstPressure -= densityW * gravity * gasPlumeDist + densityNw * gravity * (height - gasPlumeDist);
        }
      else if(veModel_ == 1)
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
      
      Scalar reconstCapPressure = 0.0;

      if(veModel_ == 1 && height <= gasPlumeDist)
	{
	  reconstCapPressure = entryP;
	}
      else if (veModel_ == 1 && height > gasPlumeDist)
	{
	  reconstCapPressure = entryP + (densityW - densityNw)*gravity*(height - gasPlumeDist);
	}

      return reconstCapPressure;
    }
    ////

  private:
    static constexpr Scalar eps_ = 1e-6;
    static constexpr Scalar depthBOR_ = 1000;

    Dumux::VtkMultiWriter<GridView> normWriter_;
    std::multimap<int, Element> mapColumns_;
    std::ofstream outputFile_;
    int veModel_;
    Element dummy_;
    Scalar CTZ_;
    Scalar segTime_;
    int numberOfColumns_;

    // additional global variables (Zakaria 2022-08-12)
    std::vector<Scalar> gasPlumeDistTemp_;
    std::chrono::_V2::system_clock::time_point beginCPU_;
    int colCounter_;
    double endTime_;
    ////
  };
}
//end namespace

#endif
