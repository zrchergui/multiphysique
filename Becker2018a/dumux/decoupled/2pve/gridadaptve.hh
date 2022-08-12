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
 * \brief Base class for h-adaptive sequential models.
 */
#ifndef DUMUX_GRIDADAPTVE_HH
#define DUMUX_GRIDADAPTVE_HH

#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/gridadaptproperties.hh>
#include <unordered_map>

namespace Dumux
{

/*!\ingroup IMPET
 * @brief Standard Module for h-adaptive simulations
 *
 * This class is created by the problem class with the template
 * parameters <TypeTag, true> and provides basic functionality
 * for adaptive methods:
 *
 * A standard implementation adaptGrid() will prepare everything
 * to calculate the next pressure field on the new grid.
 */
template<class TypeTag, bool adaptive>
class GridAdaptVE
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::Grid Grid;
    typedef typename Grid::LeafGridView LeafGridView;
    typedef typename Grid::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptionIndicator) AdaptionIndicator;
    typedef typename GET_PROP_TYPE(TypeTag, AdaptionInitializationIndicator) AdaptionInitializationIndicator;

    typedef typename Grid::Traits::LocalIdSet IdSet;
    typedef typename IdSet::IdType IdType;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef std::array<unsigned int, dim> CellArray;


public:
    /*!
     * Constructor for h-adaptive simulations (adaptive grids)
     * @param problem The problem
     */
    GridAdaptVE (Problem& problem)
        : problem_(problem), adaptionIndicator_(problem), marked_(0), coarsened_(0)
    {
        levelMin_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MinLevel);
        levelMax_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
        adaptationInterval_ = GET_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, AdaptionInterval);

        if (levelMin_ < 0)
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
                         << " :  Dune cannot coarsen to gridlevels smaller 0! "<< std::endl;
    }

    /*!
     * @brief Initalization method of the h-adaptive module
     *
     * Prepares the grid for simulation after the initialization of the
     * problem. The applied indicator is selectable via the property
     * AdaptionInitializationIndicator
     */
    void init()
    {
        adaptionIndicator_.init();

        if (!GET_PARAM_FROM_GROUP(TypeTag, bool, GridAdapt, EnableInitializationIndicator))
            return;

        AdaptionInitializationIndicator adaptionInitIndicator(problem_, adaptionIndicator_);

        int maxIter = levelMax_;
        int iter = 0;
        while (iter < maxIter)
        {
            adaptGrid(adaptionInitIndicator);

            if (!wasAdapted())
            {
                break;
            }

            int shouldInitialize = adaptionInitIndicator.initializeModel();
            if (problem_.grid().comm().max(shouldInitialize))
                problem_.model().initialize();

            iter++;
        }
    }

    /*!
     * @brief Standard method to adapt the grid
     *
     * This method is called from IMPETProblem::preTimeStep() if
     * adaptive grids are used in the simulation. It uses the standard
     * indicator (selected by the property AdaptionIndicator) and forwards to
     * with it to the ultimate method adaptGrid(indicator), which
     * uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Mark the elements
     * 3) Store primary variables in a map
     * 4) Adapt the grid, adapt variables sizes, update mappers
     * 5) Reconstruct primary variables, regain secondary variables
     */
    void adaptGrid()
    {
        adaptGrid(adaptionIndicator_) ;
    }

    /*!
     * @brief Method to adapt the grid with individual indicator vector
     *
     * @param indicator The refinement indicator that is applied
     *
     * This method is called by an user-defined preTimeStep() of
     * the applied problem and takes a given vector with indicator
     * values.
     *
     * It uses a standard procedure for adaptivity:
     * 1) Determine the refinement indicator
     * 2) Mark the elements
     * 3) Store primary variables in a map
     * 4) Adapt the grid, adapt variables sizes, update mappers
     * 5) Reconstruct primary variables, regain secondary variables
     */
    template<class Indicator>
    void adaptGrid(Indicator& indicator)
    {
        // reset internal counter for marked elements
        marked_ = coarsened_ = 0;

        // check for adaption interval: Adapt only at certain time step indices
        if (problem_.timeManager().timeStepIndex() % adaptationInterval_ != 0)
            return;

        /**** 1) determine refining parameter if standard is used ***/
        // if not, the indicatorVector and refinement Bounds have to
        // specified by the problem through setIndicator()
        indicator.calculateIndicator();

        int maxIter = levelMax_;
        int iter = 0;
        while (iter < maxIter)
        {
            /**** 2) mark elements according to indicator     *********/
            markElements(indicator);

            /****  2b) Do pre-adaption step    *****/
            problem_.grid().preAdapt();
            problem_.preAdapt();

            if(iter < 1)
            {
                /****  3) Put primary variables in a map         *********/
                problem_.variables().storePrimVars(problem_);
            }

            /****  4) Adapt Grid and size of variable vectors    *****/
            problem_.grid().adapt();

            if (!wasAdapted())
            {
                break;
            }

            problem_.updateColumnMap();
            iter++;
        }
        //        forceRefineRatio(1);

        // update mapper to new cell indices
        problem_.variables().elementMapper().update();

        // adapt size of vectors
        problem_.variables().adaptVariableSize(problem_.variables().elementMapper().size());

        /****  5) (Re-)construct primary variables to new grid **/
        problem_.variables().reconstructPrimVars(problem_);

        // delete markers in grid
        problem_.grid().postAdapt();

        // call method in problem for potential output etc.
        problem_.postAdapt();

        //set VE or fullD model according to cell size
        if (wasAdapted())
        {
            problem_.setModel();
        }

        return;
    }

    /*!
     * Mark Elements for grid refinement according to applied Indicator
     * Add buffer zone between boundary (VE columns directly next to 2D columns are refined)
     * @return Total ammount of marked cells
     */
    template<class Indicator>
    void markElements(Indicator& indicator)
    {
        CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        int numberOfColumns = numberOfCells[0];

        std::multimap<int, Element> mapAllColumns = problem_.getColumnMap();
        for(int columnNumber=0; columnNumber<numberOfColumns; ++columnNumber)
        {
            typename std::map<int, Element>::iterator it = mapAllColumns.lower_bound(columnNumber);
            // refine
            if (indicator.refine(columnNumber))
            {
                for (; it != mapAllColumns.upper_bound(columnNumber); ++it)
                {
                    Element element = it->second;

                    // only mark non-ghost elements
                    if (element.partitionType() == Dune::GhostEntity)
                        continue;

                    if (element.level() < levelMax_)
                    {
                        problem_.grid().mark(element, UG::D2::BLUE, 1);
                        ++marked_;
                    }
                }
            }
            // coarsen
            else if (indicator.coarsen(columnNumber))
            {
                //VE columns directly and two slots next to 2D columns are refined (buffer cell)
                //TODO: 3D?
                if(indicator.refine(std::min((numberOfColumns-1),columnNumber+1)) || indicator.refine(std::max(0, columnNumber-1)) ||
                        indicator.refine(std::min((numberOfColumns-2),columnNumber+1)) || indicator.refine(std::max(0, columnNumber-2)))
                {
                    for (; it != mapAllColumns.upper_bound(columnNumber); ++it)
                    {
                        Element element = it->second;

                         // only mark non-ghost elements
                         if (element.partitionType() == Dune::GhostEntity)
                             continue;

                        if (element.level() < levelMax_)
                        {
                            problem_.grid().mark(element, UG::D2::BLUE, 1);
                            ++marked_;
                        }
                    }
                }
                else
                {
                    for (; it != mapAllColumns.upper_bound(columnNumber); ++it)
                    {
                        Element element = it->second;

                         // only mark non-ghost elements
                         if (element.partitionType() == Dune::GhostEntity)
                             continue;

                         if(element.hasFather())
                         {
                             problem_.grid().mark( -1, element );
                             ++coarsened_;
                         }
                    }
                }
            }
        }
    }

    /*!
     * @brief Returns true if grid cells have been marked for adaptation
     */
    bool wasAdapted()
    {
        int sumMarked = problem_.grid().comm().sum(marked_);
        int sumCoarsened = problem_.grid().comm().sum(coarsened_);

        return (sumMarked != 0 || sumCoarsened != 0);
    }

    /*!
     * Sets minimum and maximum refinement levels
     *
     * @param levMin minimum level for coarsening
     * @param levMax maximum level for refinement
     */
    void setLevels(int levMin, int levMax)
    {
        if (levMin < 0)
            Dune::dgrave <<  __FILE__<< ":" <<__LINE__
                         << " :  Dune cannot coarsen to gridlevels smaller 0! "<< std::endl;
        levelMin_ = levMin;
        levelMax_ = levMax;
    }

    /*!
     * @brief Returns maximum refinement level
     *
     * The value is the assign maximum possible level,
     * not the actual maximum level of the grid.
     * @return levelMax_ maximum level for refinement
     */
    const int getMaxLevel() const
    {
        return levelMax_;
    }
    /*!
     * @brief Returns minimum refinement level
     *
     * The value is the assign minimum possible level,
     * not the actual minimum level of the grid.
     * @return levelMin_ minimum level for coarsening
     */
    const int getMinLevel() const
    {
        return levelMin_;
    }

    AdaptionIndicator& adaptionIndicator()
    {
        return adaptionIndicator_;
    }

    AdaptionIndicator& adaptionIndicator() const
    {
        return adaptionIndicator_;
    }

private:
    /*!
     * @brief Method ensuring the refinement ratio of 2:1
     *
     * For any given entity, a loop over the neighbors checks weather the
     * entities refinement would require that any of the neighbors has
     * to be refined, too.
     * This is done recursively over all levels of the grid.
     *
     * @param entity Element of interest that is to be refined
     * @param level level of the refined entity: it is at least 1
     * @return true if everything was successful
     */
    bool checkNeighborsRefine_(const Element &entity, int level = 1)
    {
        // this also refines the neighbor elements
        for(const auto& intersection : intersections(problem_.gridView(), entity))
        {
            if(!intersection.neighbor())
                continue;

            auto outside = intersection.outside();

            // only mark non-ghost elements
            if (outside.partitionType() == Dune::GhostEntity)
                continue;

            if ((outside.level() < levelMax_)
                && (outside.level() < entity.level()))
            {
                problem_.grid().mark(1, outside);
                ++marked_;

                if(level != levelMax_)
                    checkNeighborsRefine_(outside, ++level);
            }
        }
        return true;
    }


    /*!
     * \brief Enforces a given refine ratio after grid was adapted
     *
     * If the refine ratio is not taken into consideration during
     * marking, then this method ensures a certain ratio.
     *
     * @param maxLevelDelta The maximum level difference (refine ratio)
     *             between neighbors.
     */
    void forceRefineRatio(int maxLevelDelta = 1)
    {
        LeafGridView leafGridView = problem_.gridView();
        // delete all existing marks
        problem_.grid().postAdapt();
        bool done;
        do
        {
            // run through all cells
            done=true;
            for (const auto& element : elements(problem_.gridView()))
            {
                // only mark non-ghost elements
                if (element.partitionType() == Dune::GhostEntity)
                    continue;

                // run through all neighbor-cells (intersections)
                for (const auto& intersection : intersections(leafGridView, element))
                {
                    if(!intersection.neighbor())
                        continue;

                    if (element.level() + maxLevelDelta < intersection.outside().level())
                    {
                        problem_.grid().mark( 1, element );
                        done=false;
                    }
                }
            }
            if (done==false)
            {
                // adapt the grid
                problem_.grid().adapt();
                // delete marks
                problem_.grid().postAdapt();
            }
        }
        while (done!=true);
    }

    // private Variables
    Problem& problem_;
    AdaptionIndicator adaptionIndicator_;

    int marked_;
    int coarsened_;

    int levelMin_;
    int levelMax_;

    int adaptationInterval_;
    int numRefine_;
};

/*!
 * @brief Class for NON-adaptive simulations
 *
 * This class provides empty methods for non-adaptive simulations
 * for compilation reasons. If adaptivity is desired, create the
 * class with template arguments <TypeTag, true> instead.
 */
template<class TypeTag>
class GridAdaptVE<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem)     Problem;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

public:
    void init()
    {}
    void adaptGrid()
    {}
    bool wasAdapted()
    {
        return false;
    }
    void setLevels(int, int)
    {}
    void setTolerance(int, int)
    {}
    const void setIndicator(const ScalarSolutionType&,
                            const Scalar&, const Scalar&)
    {}
    GridAdaptVE (Problem& problem)
    {}
};

}
#endif /* DUMUX_GRIDADAPT_HH */
