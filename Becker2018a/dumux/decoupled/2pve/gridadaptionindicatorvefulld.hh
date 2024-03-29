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
#ifndef DUMUX_GRIDADAPTIONINDICATOR2PVEFULLD_HH
#define DUMUX_GRIDADAPTIONINDICATOR2PVEFULLD_HH

#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/linear/vectorexchange.hh>

/**
 * @file
 * @brief  Class defining a standard, saturation dependent indicator for grid adaption
 */
namespace Dumux
{
    /*!\ingroup IMPES
     * @brief  Class defining a standard, saturation dependent indicator for grid adaption
     *
     * \tparam TypeTag The problem TypeTag
     */
    template<class TypeTag>
    class GridAdaptionIndicator2PVEFullD
    {
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GridView::Traits::template Codim<0>::Entity Element;

        typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
        typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;
        typedef typename SolutionTypes::ElementMapper ElementMapper;

        typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
        typedef typename SpatialParams::MaterialLaw MaterialLaw;

        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

        typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
        typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonWettingPhase;

        enum
            {
                sw = Indices::saturationW,
                sn = Indices::saturationNw
            };
        enum
            {
                wPhaseIdx = Indices::wPhaseIdx,
                nPhaseIdx = Indices::nPhaseIdx
            };
        enum
            {
                dim = GridView::dimension, dimWorld = GridView::dimensionworld
            };

        typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
        typedef std::array<unsigned int, dim> CellArray;
        typedef Dune::FieldVector<Scalar, dim> DimVector;


    public:
        /*! \brief Calculates the indicator used for refinement/coarsening for each column.
         *
         * This indicator is based on the saturation profile of the wetting phase compared to the VE-saturation profile.
         */
        void calculateIndicator()
        {
            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
            int numberOfColumns = numberOfCells[0];
            Scalar domainHeight = problem_.bBoxMax()[dim - 1];
            const int maxLevel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, GridAdapt, MaxLevel);
            double deltaZ = domainHeight/(numberOfCells[dim - 1]*std::pow(2, maxLevel));
            // prepare an indicator for refinement
            if(indicatorVector_.size() != numberOfColumns)
                {
                    indicatorVector_.resize(numberOfColumns);
                }
            indicatorVector_ = 0.0;

            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            std::multimap<int, Element> mapAllColumns = problem_.getColumnMap();
            for(int i=0;i<numberOfColumns;i++)
                {
                    if(i < 4)//never have VE in the first four columns
                        {
                            indicatorVector_[i] = 2;
                            continue;
                        }
                    int columnModel = problem_.getColumnModel(i);
                    if(columnModel == 0 || columnModel == 1) 
                        continue;
                    Scalar averageSatColumn = 0.0;// average wetting saturation in column
                    Scalar totalVolume = 0.0;// total volume of column
                    Scalar gasVolume = 0.0;// volume with SatN>0.0 in one column;
                    typename std::map<int, Element>::iterator it = mapAllColumns.lower_bound(i);
                    dummy_ = it->second;
                    for (; it != mapAllColumns.upper_bound(i); ++it)
                        {
                            Element element = it->second;
                            int globalIdxI = problem_.variables().index(element);
                            Scalar satW = problem_.variables().cellData(globalIdxI).saturation(wPhaseIdx);
                            Scalar volume = element.geometry().volume();
                            averageSatColumn += satW * volume;
                            totalVolume += volume;
                        }
                    averageSatColumn = averageSatColumn/totalVolume;

                    // to make sure the average saturation does not exceed 1, otherwise Newton method may not converge to compute the gas plume thickness (Zakaria 2022-09-12)
                    averageSatColumn = std::min(averageSatColumn,1.);
                    ////

                    Scalar gasPlumeDist = calculateGasPlumeDist(averageSatColumn);

                    Scalar errorSat = 0.0;
                    Scalar errorRelPerm = 0.0;

                    // initialize the additional criteria to choose from - pressure only for capillary fringe model (Zakaria 2022-08-12 and 2022-08-30)
                    Scalar errorSatNorm = 0.0;
                    Scalar errorRelPermNorm = 0.0;
                    Scalar errorPress = 0.0;
                    Scalar errorPressNorm = 0.0;
                    Scalar errorCapPress = 0.0;
                    Scalar errorCapPressNorm = 0.0;
                    ////

                    it = mapAllColumns.lower_bound(i);

                    // store bottom pressure of column i to compute hydrostatic pressure later (Zakaria 22-08-12)
                    Element veElement = mapAllColumns.lower_bound(i)->second;
                    int eIdxGlobal = problem_.variables().index(veElement);
                    Scalar coarsePressW = problem_.variables().cellData(eIdxGlobal).pressure(wPhaseIdx);
                    Scalar bottomPressW = reconstPressureW(-deltaZ/2,gasPlumeDist,coarsePressW); //slightly different from coarsePressW when columnModel > 1 since coarsePressW is in this case the pressure at the center of the bottom cell and not the pressure at bottom
                    ////

                    for (; it != mapAllColumns.upper_bound(i); ++it)
                        {
                            Element element = it->second;
                            Scalar resSatW = problem_.spatialParams().materialLawParams(element).swr();
                            int globalIdxI = problem_.variables().index(element);
                            GlobalPosition globalPos = (element).geometry().center();
                            Scalar top = globalPos[dim - 1] + deltaZ/2.0;
                            Scalar bottom = globalPos[dim - 1] - deltaZ/2.0;

                            Scalar satW = problem_.variables().cellData(globalIdxI).saturation(wPhaseIdx);
                            Scalar krw = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element), satW);
                            // pressure in cell (Zakaria 2022-08-12)
                            Scalar pressW = problem_.variables().cellData(globalIdxI).pressure(wPhaseIdx);
                            Scalar pc = MaterialLaw::pc(problem_.spatialParams().materialLawParams(element), satW);
                            Scalar pe = problem_.spatialParams().materialLawParams(element).pe();
                            ////

                            if(veModel == 0)//calculate error for VE model
                                {
                                    if (top <= gasPlumeDist)
                                        {
                                            errorSat += std::abs(deltaZ * (satW - 1.0));
                                            errorRelPerm += std::abs(deltaZ * (krw - 1.0));

                                            // normalized errors (Zakaria 2022-08-12) 
                                            errorSatNorm += (1/satW)* std::abs(deltaZ * (satW - 1.0));
                                            errorRelPermNorm += (1/krw)* std::abs(deltaZ * (krw - 1.0));
                                            ////
                                        }
                                    else if (bottom >= gasPlumeDist)
                                        {
                                            errorSat += std::abs(deltaZ * (satW - resSatW));
                                            errorRelPerm += std::abs(deltaZ * (krw - 0.0));

                                            // normalized errors (Zakaria 2022-08-12) 
                                            errorSatNorm += (1/satW)* std::abs(deltaZ * (satW - resSatW));;
                                            errorRelPermNorm += (1/krw)* std::abs(deltaZ * (krw - 0.0));
                                            ////
                                        }
                                    else
                                        {
                                            Scalar lowerDelta = gasPlumeDist - bottom;
                                            Scalar upperDelta = top - gasPlumeDist;
                                            errorSat += std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW));
                                            errorRelPerm += std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0));

                                            // normalized errors (Zakaria 2022-08-12) 
                                            errorSatNorm += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + std::abs(upperDelta * (satW - resSatW)));
                                            errorRelPermNorm += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + std::abs(upperDelta * (krw - 0.0)));
                                            ////                                            
                                        }
                                }
                            if(veModel == 1)//calculate error for capillary fringe model
                                {
                                    if (top <= gasPlumeDist)
                                        {
                                            errorSat += std::abs(deltaZ * (satW - 1.0));
                                            errorRelPerm += std::abs(deltaZ * (krw - 1.0));

                                            // normalized errors (Zakaria 2022-08-12)
                                            errorSatNorm += (1/satW)* ( std::abs(deltaZ * (satW - 1.0)) ) ;
                                            errorRelPermNorm += (1/krw)* ( std::abs(deltaZ * (krw - 1.0)) ) ;
                                            ////
                                            
                                            // capillary pressure errors (Zakaria 2022-08-30)
                                            errorCapPress += std::abs(deltaZ * (pc - pe));
                                            errorCapPressNorm += (1/pc)* std::abs(deltaZ * (pc - pe));
                                            ////
                                        }
                                    else if (bottom >= gasPlumeDist)
                                        {
                                            Scalar eIntSat = calculateErrorSatIntegral(bottom, top, satW, gasPlumeDist);
                                            errorSat += eIntSat;

                                            Scalar eIntPerm = calculateErrorKrwIntegral(bottom, top, satW, gasPlumeDist);
                                            errorRelPerm += eIntPerm;
                        
                                            // errors in normalized variables (Zakaria 2022-08-12)
                                            errorSatNorm += (1/satW)* eIntSat;
                                            errorRelPermNorm += (1/krw)* eIntPerm;
                                            ////

                                            // capillary pressure errors (Zakaria 2022-08-30)
                                            Scalar eIntPc = calculateErrorCapPressIntegral(bottom, top, pc, gasPlumeDist);
                                            errorCapPress += eIntPc;
                                            errorCapPressNorm += (1/pc)* eIntPc;
                                            ////
                                        }
                                    else
                                        {
                                            Scalar lowerDelta = gasPlumeDist - bottom;
                                            Scalar upperDelta = top - gasPlumeDist;

                                            Scalar eIntSat = calculateErrorSatIntegral(gasPlumeDist, top, satW, gasPlumeDist);
                                            errorSat += std::abs(lowerDelta * (satW - 1.0)) + eIntSat;

                                            Scalar eIntPerm = calculateErrorKrwIntegral(gasPlumeDist, top, satW, gasPlumeDist);
                                            errorRelPerm += std::abs(lowerDelta * (krw - 1.0)) + eIntPerm;

                                            // error in normalized variables (Zakaria 2022-08-12)
                                            errorSatNorm += (1/satW)* (std::abs(lowerDelta * (satW - 1.0)) + calculateErrorSatIntegral(gasPlumeDist, top, satW, gasPlumeDist));
                                            errorRelPermNorm += (1/krw)* (std::abs(lowerDelta * (krw - 1.0)) + calculateErrorKrwIntegral(gasPlumeDist, top, satW, gasPlumeDist));
                                            ////

                                            // capillary pressure errors (Zakaria 2022-08-30)
                                            Scalar eIntPc = calculateErrorCapPressIntegral(gasPlumeDist, top, pc, gasPlumeDist);
                                            errorCapPress += std::abs(lowerDelta * (pc - pe)) + eIntPc;
                                            errorCapPressNorm += (1/pc)* (std::abs(lowerDelta * (pc - pe)) + eIntPc);
                                            ////
                                        }
                                }
                            // error in water pressure (Zakaria 2022-08-26)
                            Scalar eIntPress = calculateErrorPressIntegral(bottom, top, pressW, gasPlumeDist, bottomPressW);
                            errorPress += eIntPress;
                            errorPressNorm += (1/pressW)* eIntPress;
                            ////
                        }
                    //indicatorVector_[i] = errorSat/(domainHeight - gasPlumeDist);
                    indicatorVector_[i] = errorRelPerm/(domainHeight - gasPlumeDist);

                    // update of pressure and normalized errors (Zakaria 2022-08-12)
                    //indicatorVector_[i] =  errorSatNorm/(domainHeight - gasPlumeDist);
                    //indicatorVector_[i] =  errorRelPermNorm/(domainHeight - gasPlumeDist);
                    //indicatorVector_[i] =  errorPress/(domainHeight - gasPlumeDist);
                    //indicatorVector_[i] =  errorPressNorm/(domainHeight - gasPlumeDist);
                    //indicatorVector_[i] =  errorCapPress/(domainHeight - gasPlumeDist);
                    //indicatorVector_[i] =  errorCapPressNorm/(domainHeight - gasPlumeDist);
                    ////
                }
            Scalar absoluteError = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, GridAdapt, AbsoluteError);
            refineBound_ = 1.05 * absoluteError;
            coarsenBound_= 0.95 * absoluteError;

#if HAVE_MPI
            // communicate updated values
            typedef VectorExchange<ElementMapper, ScalarSolutionType> DataHandle;
            DataHandle dataHandle(problem_.elementMapper(), indicatorVector_);
            problem_.gridView().template communicate<DataHandle>(dataHandle,
                                                                 Dune::InteriorBorder_All_Interface,
                                                                 Dune::ForwardCommunication);

            refineBound_ = problem_.gridView().comm().max(refineBound_);
            coarsenBound_ = problem_.gridView().comm().max(coarsenBound_);

#endif
        }

        /*! \brief Calculates gasPlumeDist, distance of gas plume from bottom
         *
         * Stores minGasPlumeDist for all grid cells
         */
        Scalar calculateGasPlumeDist(Scalar satW)
        {
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);
            Scalar domainHeight = problem_.bBoxMax()[dim - 1];
            Scalar resSatW = problem_.spatialParams().materialLawParams(dummy_).swr();
            Scalar resSatN = problem_.spatialParams().materialLawParams(dummy_).snr();
            Scalar gravity = problem_.gravity().two_norm();

            Scalar gasPlumeDist = 0.0;

            if (veModel == 0) //calculate gasPlumeDist for sharp interface ve model
                {
                    gasPlumeDist = domainHeight * (satW - resSatW) / (1.0 - resSatW);
                }

            else if (veModel == 1) //calculate gasPlumeDist for capillary fringe model
                {
                    GlobalPosition globalPos = dummy_.geometry().center();
                    Scalar pRef = problem_.referencePressureAtPos(globalPos);
                    Scalar tempRef = problem_.temperatureAtPos(globalPos);
                    Scalar densityW = WettingPhase::density(tempRef, pRef);
                    Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
                    Scalar lambda = problem_.spatialParams().materialLawParams(dummy_).lambda();
                    Scalar entryP = problem_.spatialParams().materialLawParams(dummy_).pe();

                    Scalar Xi = domainHeight / 2.0; //XiStart

                    Scalar fullIntegral = 1.0 / (1.0 - lambda) * (1.0 - resSatW - resSatN) / ((densityW - densityNw) * gravity) * (std::pow(entryP, lambda)
                                                                                                                                   - std::pow(entryP, 2.0 - lambda) + std::pow((domainHeight * (densityW - densityNw) * gravity + entryP), (1.0 - lambda)));
                    //GasPlumeDist>0
                    if (fullIntegral < satW * domainHeight)
                        {
                            //solve equation for
                            for (int count = 0; count < 100; count++)
                                {
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
                    Scalar krw1 = MaterialLaw::krw(problem_.spatialParams().materialLawParams(dummy_), sat1);
                    Scalar krw2 = MaterialLaw::krw(problem_.spatialParams().materialLawParams(dummy_), sat2);
                    Scalar krw = MaterialLaw::krw(problem_.spatialParams().materialLawParams(dummy_), satW);
                    krwIntegral += std::abs((krw1 + krw2)/2.0 - krw);
                }
            krwIntegral = krwIntegral * deltaZ;

            return krwIntegral;
        }

        // Calculates integral of difference of pressure over z (Zakaria 2022-08-12)
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
            Scalar domainHeight = problem_.bBoxMax()[dim - 1];
            GlobalPosition globalPos = dummy_.geometry().center();
            Scalar pRef = problem_.referencePressureAtPos(globalPos);
            Scalar tempRef = problem_.temperatureAtPos(globalPos);
            Scalar resSatW = problem_.spatialParams().materialLawParams(dummy_).swr();
            Scalar resSatN = problem_.spatialParams().materialLawParams(dummy_).snr();
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar entryP = problem_.spatialParams().materialLawParams(dummy_).pe();
            Scalar lambda = problem_.spatialParams().materialLawParams(dummy_).lambda();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar reconstSaturation = 0.0;

            if (veModel == 0) //reconstruct phase saturation for sharp interface ve model
                {
                    reconstSaturation = 1.0;
                    if(height > gasPlumeDist)
                        {
                            reconstSaturation = resSatW;
                        }
                }
            else if (veModel == 1) //reconstruct phase saturation for capillary fringe model
                {
                    reconstSaturation = 1.0;
                    if(height > gasPlumeDist)
                        {
                            reconstSaturation = std::pow(((height - gasPlumeDist) * (densityW - densityNw) * problem_.gravity().two_norm() + entryP), (-lambda))
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
        Scalar reconstPressureW(Scalar height, Scalar gasPlumeDist, Scalar coarsePressW)
        {
            GlobalPosition globalPos = dummy_.geometry().center();
            Scalar pRef = problem_.referencePressureAtPos(globalPos);
            Scalar tempRef = problem_.temperatureAtPos(globalPos);
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar gravity = problem_.gravity().two_norm();
            int veModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, VE, VEModel);

            Scalar reconstPressure = coarsePressW;

            if(veModel == 0 && height <= gasPlumeDist)
                {
                    reconstPressure -= densityW * gravity * height;
                }
            else if (veModel == 0 && height > gasPlumeDist)
                {
                    reconstPressure -= densityW * gravity * gasPlumeDist + densityNw * gravity * (height - gasPlumeDist);
                }
            else if(veModel == 1)
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
            Scalar pRef = problem_.referencePressureAtPos(globalPos);
            Scalar tempRef = problem_.temperatureAtPos(globalPos);
            Scalar densityW = WettingPhase::density(tempRef, pRef);
            Scalar densityNw = NonWettingPhase::density(tempRef, pRef);
            Scalar gravity = problem_.gravity().two_norm();
            Scalar entryP = problem_.spatialParams().materialLawParams(dummy_).pe();
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

        //        /*! \brief Calculates the indicator used for refinement/coarsening for each column.
        //         *
        //         * This indicator is based on the velocity of the non-wetting phase. TODO: not working, implement like in test_dec2p
        //         */
        //        void calculateIndicator()
        //        {
        //            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        //            int numberOfColumns = numberOfCells[0];
        //            // prepare an indicator for refinement
        //            if(indicatorVector_.size() != numberOfColumns)
        //            {
        //                indicatorVector_.resize(numberOfColumns);
        //            }
        //            indicatorVector_ = 0.0;
        //
        //            Scalar velocityNwVerticalMaxTotal = 0.0;
        //            std::multimap<int, Element> mapAllColumns = problem_.getColumnMap();
        //            for(int i=0;i<numberOfColumns;i++)
        //            {
        //                Scalar velocityNwVerticalMax = 0.0;
        //                typename std::map<int, Element>::iterator it = mapAllColumns.lower_bound(i);
        //                for (; it != mapAllColumns.upper_bound(i); ++it)
        //                {
        //                    Element element = it->second;
        //                    int eIdx = problem_.variables().index(element);
        //                    for (const auto& intersection : intersections(problem_.gridView(), element))
        //                    {
        //                        GlobalPosition velocityNw = problem_.variables().cellData(eIdx).fluxData().velocity(nPhaseIdx, intersection.indexInInside());
        //                        GlobalPosition gravityNormalized = problem_.gravity();
        //                        gravityNormalized /= problem_.gravity().two_norm();
        //                        velocityNwVerticalMax = std::max(velocityNwVerticalMax, -(velocityNw * gravityNormalized));
        //                    }
        //                }
        //                indicatorVector_[i] = velocityNwVerticalMax;
        //                velocityNwVerticalMaxTotal = std::max(velocityNwVerticalMaxTotal, velocityNwVerticalMax);
        //            }
        //            refineBound_ = refinetol_ * velocityNwVerticalMaxTotal;
        //            coarsenBound_ = coarsentol_ * velocityNwVerticalMaxTotal;
        //
        //            #if HAVE_MPI
        //            // communicate updated values
        //            typedef VectorExchange<ElementMapper, ScalarSolutionType> DataHandle;
        //            DataHandle dataHandle(problem_.elementMapper(), indicatorVector_);
        //            problem_.gridView().template communicate<DataHandle>(dataHandle,
        //                                                 Dune::InteriorBorder_All_Interface,
        //                                                 Dune::ForwardCommunication);
        //
        //            refineBound_ = problem_.gridView().comm().max(refineBound_);
        //            coarsenBound_ = problem_.gridView().comm().max(coarsenBound_);
        //
        //            #endif
        //        }








        //        /*! \brief Calculates the indicator used for refinement/coarsening for each column.
        //         *
        //         * This indicator is based on the saturation and the position.
        //         */
        //        void calculateIndicator()
        //        {
        //            CellArray numberOfCells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, "Grid", Cells);
        //            int numberOfColumns = numberOfCells[0];
        //            // prepare an indicator for refinement
        //            if(indicatorVector_.size() != numberOfColumns)
        //            {
        //                indicatorVector_.resize(numberOfColumns);
        //            }
        //            indicatorVector_ = 0.0;
        //            std::multimap<int, Element> mapAllColumns = problem_.getColumnMap();
        //            for(int i=0;i<numberOfColumns;i++)
        //            {
        //                typename std::map<int, Element>::iterator it = mapAllColumns.lower_bound(i);
        //                for (; it != mapAllColumns.upper_bound(i); ++it)
        //                {
        //                    int globalIdxI = problem_.variables().index(it->second);
        //                    GlobalPosition globalPos = (it->second).geometry().center();
        //                    Scalar satW = problem_.variables().cellData(globalIdxI).saturation(wPhaseIdx);
        //                    //mark 2D column depending on indicator
        //                    if(problem_.timeManager().time() < 5e4)
        //                    {
        //                        if(globalPos[0] < 20.0)
        //                        {
        //                            indicatorVector_[i] = 1.0;
        //                            continue;
        //                        }
        //                        else if(satW < 0.98)
        //                        {
        //                            indicatorVector_[i] = 1.0;
        //                            continue;
        //                        }
        //                    }
        //                    else if(problem_.timeManager().time() > 5e4 && problem_.timeManager().time() < 7e4)
        //                    {
        //                        if(globalPos[0] < 40.0)
        //                        {
        //                            indicatorVector_[i] = 1.0;
        //                            continue;
        //                        }
        //                    }
        //                    else
        //                    {
        //                        if(globalPos[0] < 20.0)
        //                        {
        //                            indicatorVector_[i] = 1.0;
        //                            continue;
        //                        }
        //                    }
        //                }
        //            }
        //            refineBound_ = absoluteError;
        //            coarsenBound_ = absoluteError;
        //            #if HAVE_MPI
        //            // communicate updated values
        //            typedef VectorExchange<ElementMapper, ScalarSolutionType> DataHandle;
        //            DataHandle dataHandle(problem_.elementMapper(), indicatorVector_);
        //            problem_.gridView().template communicate<DataHandle>(dataHandle,
        //                                     Dune::InteriorBorder_All_Interface,
        //                                     Dune::ForwardCommunication);
        //
        //            refineBound_ = problem_.gridView().comm().max(refineBound_);
        //            coarsenBound_ = problem_.gridView().comm().max(coarsenBound_);
        //
        //            #endif
        //        }


        /*! \brief Indicator function for marking of grid cells for refinement
         *
         * Returns true if an element should be refined.
         *
         *  \param element A grid element
         */
        bool refine(int columnNumber)
        {
            return (indicatorVector_[columnNumber] > refineBound_);
        }

        /*! \brief Indicator function for marking of grid cells for coarsening
         *
         * Returns true if an element should be coarsened.
         *
         *  \param element A grid element
         */
        bool coarsen(int columnNumber)
        {
            return (indicatorVector_[columnNumber] < coarsenBound_);
        }

        /*! \brief Initializes the adaption indicator class*/
        void init()
        {
            refineBound_ = 0.;
            coarsenBound_ = 0.;
        };

        /*! @brief Constructs a GridAdaptionIndicator instance
         *
         *  This standard indicator is based on the saturation gradient.
         *  It checks the local gradient compared to the maximum global gradient.
         *  The indicator is compared locally to a refinement/coarsening threshold to decide whether
         *  a cell should be marked for refinement or coarsening or should not be adapted.
         *
         * \param problem The problem object
         */
        GridAdaptionIndicator2PVEFullD (Problem& problem):
            problem_(problem)
        {
            refinetol_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, GridAdapt, RefineTolerance);
            coarsentol_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, GridAdapt, CoarsenTolerance);
        }

    protected:
        Problem& problem_;
        Scalar refinetol_;
        Scalar coarsentol_;
        Scalar refineBound_;
        Scalar coarsenBound_;
        ScalarSolutionType indicatorVector_;
        static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
        Element dummy_;
    };
}

#endif
