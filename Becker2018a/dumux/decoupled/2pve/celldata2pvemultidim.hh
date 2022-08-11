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
#ifndef DUMUX_CELLDATA2PVEMULTIDIM_HH
#define DUMUX_CELLDATA2PVEMULTIDIM_HH

#include "celldata2pve.hh"

/**
 * \file
 * \brief  Class including data of one grid cell
 */

namespace Dumux
{

/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations for one grid cell.
/*! The variables of two-phase flow, which are phase pressures and saturations are stored in this class. Further, resulting cell values for constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure are stored.
 * Additionally, data assigned to cell-cell interfaces, so-called flux-data are stored.
 *
 * \tparam TypeTag The problem TypeTag
 * \tparam bool Used for specialization: in case of incompressible flow bool = <tt>false</tt>
 */
template<class TypeTag>
class CellData2PVEMultiDim: public CellData2PVE<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef std::array<unsigned int, dim> CellArray;

private:
    int veModel_;
    Scalar oldPressure_[numPhases];

public:
    //! Collection of variables that have to be mapped if the grid is adapted
    /**
     * For an immiscible two-phase model, the following data has to be
     * transferred to a new grid.
     */
    struct AdaptedValues
    {
        Scalar saturationW;
        Scalar saturationNw;
        Scalar pressW;
        Scalar pressNw;
        Scalar potW;
        Scalar potNw;
        Scalar volCorr;
        Scalar gasPlumeDist;
        Scalar minGasPlumeDist;
        int count;
        int veModel;
        AdaptedValues()
        {
            saturationW = 0.;
            saturationNw = 0.;
            pressW = 0.;
            pressNw = 0.;
            potW = 0.;
            potNw = 0.;
            count = 0;
            gasPlumeDist = 0;
            minGasPlumeDist = 0;
            volCorr = 0;
            veModel = 0;

        }
    };

    typedef AdaptedValues LoadBalanceData;

    //! Constructs a CellData2P object
    CellData2PVEMultiDim()
    {
        for (int i = 0; i < numPhases;i++)
        {
            oldPressure_[i] = 0.0;
        }
        veModel_ = 2;
    }

    int& veModel()
    {
        return veModel_;
    }
    //! Return subdomain information
    /** Acess function to get subdomain information
     */
    const int& veModel() const
    {
        return veModel_;
    }
    //! Specify subdomain information and fluidStateType
    /** This function is only called if
     */
    void setVeModel(int index)
    {
        veModel_ = index;
    }

    /*!\brief Returns the cell phase pressure of the last time step
     *
     * \param phaseIdx Index of a fluid phase
     */
    Scalar oldPressure(int phaseIdx)
    {
        return oldPressure_[phaseIdx];
    }

    /*!\brief Stores the old cell phase pressure
     *
     * \param phaseIdx Index of a fluid phase
     * \param sat Phase saturation which is stored
     */
    void setOldPressure(int phaseIdx, Scalar pressure)
    {
        oldPressure_[phaseIdx] = pressure;
    }

    //! Stores values to be adapted in an adaptedValues container
    /**
     * Stores values to be adapted from the current CellData objects into
     * the adaptation container in order to be mapped on a new grid.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element to be stored
     */
    void storeAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        adaptedValues.saturationW = this->saturation(wPhaseIdx);
        adaptedValues.saturationNw = this->saturation(nPhaseIdx);
        adaptedValues.pressW = this->pressure(wPhaseIdx);
        adaptedValues.pressNw = this->pressure(nPhaseIdx);
        adaptedValues.potW = this->potential(wPhaseIdx);
        adaptedValues.potNw = this->potential(nPhaseIdx);
        adaptedValues.gasPlumeDist = this->gasPlumeDist();
        adaptedValues.minGasPlumeDist = this->minGasPlumeDist();
        adaptedValues.volCorr = this->volumeCorrection();
        adaptedValues.veModel = this->veModel();
    }
    //! Stores sons entries into father element for averaging
    /**
     * Sum up the adaptedValues (sons values) into father element. We store from leaf
     * upwards, so sons are stored first, then cells on the next leaf (=fathers)
     * can be averaged. The stored values are wrong by one factor which is the number of sons
     * associated with the current father. This is corrected in setAdaptionValues by dividing
     * there through the count.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param adaptedValuesFather Values to be adapted of father cell
     * \param fatherElement The element of the father
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather,
                                    const Element& fatherElement,
                                    Problem& problem)
    {
        adaptedValuesFather.saturationW += adaptedValues.saturationW / adaptedValues.count;
        adaptedValuesFather.saturationNw += adaptedValues.saturationNw / adaptedValues.count;
        Scalar densityW = problem.pressureModel().density(wPhaseIdx);
        Scalar gravity = problem.pressureModel().gravity();
        CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        int MaxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
        double deltaZ = problem.bBoxMax()[dim-1]/(numberOfCells[dim-1]*std::pow(2, MaxLevel));
        adaptedValuesFather.pressW = std::max(adaptedValues.pressW, adaptedValuesFather.pressW);
        if(fatherElement.level() == 0 && adaptedValuesFather.count == 2)//TODO: this is hardcoded for blue refinement
        {
            adaptedValuesFather.pressW = (adaptedValuesFather.pressW + densityW*gravity*deltaZ/2.0) * adaptedValuesFather.count;
        }
        adaptedValuesFather.pressNw = adaptedValues.pressNw / adaptedValues.count;//TODO: It is not necessary to have a correct value here, right? (non-wetting pressure and potentials are set in updateMaterialLaw directly after grid adaptation)
        adaptedValuesFather.potW += adaptedValues.potW / adaptedValues.count;
        adaptedValuesFather.potNw += adaptedValues.potNw / adaptedValues.count;
        adaptedValuesFather.volCorr += adaptedValues.volCorr / adaptedValues.count;
    }
    //! Set adapted values in CellData
    /**
     * This methods stores reconstructed values into the cellData object, by
     * this setting a newly mapped solution to the storage container of the
     * sequential models.
     *
     * \param adaptedValues Container for model-specific values to be adapted
     * \param element The element where things are stored.
     */
    void setAdaptionValues(AdaptedValues& adaptedValues, const Element& element)
    {
        this->setSaturation(wPhaseIdx, adaptedValues.saturationW / adaptedValues.count);
        this->setSaturation(nPhaseIdx, adaptedValues.saturationNw / adaptedValues.count);
        this->setPressure(wPhaseIdx, adaptedValues.pressW/ adaptedValues.count);
        this->setPressure(nPhaseIdx, adaptedValues.pressNw / adaptedValues.count);
        this->setPotential(wPhaseIdx, adaptedValues.potW / adaptedValues.count);
        this->setPotential(nPhaseIdx, adaptedValues.potNw / adaptedValues.count);
        this->setUpdate(adaptedValues.volCorr / adaptedValues.count);
    }

    //! Reconstructs sons entries from data of father cell
    /**
     * Reconstructs a new solution from a father cell into a newly
     * generated son cell. New cell is stored into the global
     * adaptationMap.
     *
     * \param adaptionMap Global map storing all values to be adapted
     * \param father Entity Pointer to the father cell
     * \param son Entity Pointer to the newly created son cell
     * \param problem The problem
     */
    static void reconstructAdaptionValues(Dune::PersistentContainer<Grid, AdaptedValues>& adaptionMap,
            const Element& father, const Element& son, Problem& problem)
    {
        AdaptedValues& adaptedValues = adaptionMap[son];
        AdaptedValues& adaptedValuesFather = adaptionMap[father];
        GlobalPosition globalPos = son.geometry().center();
        CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        int MaxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
        double deltaZ = problem.bBoxMax()[dim-1]/(numberOfCells[dim-1]*std::pow(2, MaxLevel));
        Scalar top = globalPos[dim-1] + deltaZ/2.0;
        Scalar bottom = globalPos[dim-1] - deltaZ/2.0;
        adaptedValues.saturationW = problem.pressureModel().saturationIntegral(bottom, top, adaptedValuesFather, father)/(top-bottom);
//        std:: cout << "satWCelldata " << adaptedValues.saturationW << std::endl;;
        adaptedValues.saturationNw = 1.0-adaptedValues.saturationW;
        adaptedValues.pressW = problem.pressureModel().reconstPressure(globalPos[dim-1], wPhaseIdx, adaptedValuesFather, father);
        adaptedValues.pressNw = problem.pressureModel().reconstPressure(globalPos[dim-1], nPhaseIdx, adaptedValuesFather, father);
//        std::cout << "fatherPressW " << adaptedValuesFather.pressW  << "hight " << globalPos[dim-1] << "pressW " << adaptedValues.pressW << std::endl;
        adaptedValues.potW = problem.pressureModel().reconstPotential(globalPos[dim-1], wPhaseIdx, adaptedValuesFather, father);
        adaptedValues.potNw = problem.pressureModel().reconstPotential(globalPos[dim-1], nPhaseIdx, adaptedValuesFather, father);
        adaptedValues.volCorr = adaptedValuesFather.volCorr;

//        adaptedValues.saturationW = adaptedValuesFather.saturationW / adaptedValuesFather.count;
//        adaptedValues.saturationNw = adaptedValuesFather.saturationNw / adaptedValuesFather.count;
//        adaptedValues.pressW = adaptedValuesFather.pressW / adaptedValuesFather.count;
//        adaptedValues.pressNw = adaptedValuesFather.pressNw / adaptedValuesFather.count;
//        adaptedValues.potW = adaptedValuesFather.potW / adaptedValuesFather.count;
//        adaptedValues.potNw = adaptedValuesFather.potNw / adaptedValuesFather.count;
//        adaptedValues.volCorr = adaptedValuesFather.volCorr / adaptedValuesFather.count;
    }

};
}
#endif
