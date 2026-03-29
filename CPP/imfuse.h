#ifndef _IMFUSE_H
#define _IMFUSE_H

#define ITK_LEGACY_FUTURE_REMOVE

#include "debugging.h"

#include "imtypes.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include "gammatransf.h"

// @defgroup ConstantsForImfuseImageGeneration
/**
 @addtogroup ConstantsForImfuseImageGeneration
 @{
*/
constexpr unsigned int neisize = 1;    /*!< The constant to set the size of the pixel drawn for each event */

/**
 @brief AdjustModes documentation
         Enumeration with the possibilities to apply a gamma transformation to obtain the output images.
         NoAdj mean not to apply any gamma transformation, Up means to apply a convex one (decrease low gray values, increase high gray values) and Down means to apply a concave one (increase low gray values, decrease high gray values)
*/
enum AdjustModes { NoAdj, Up, Down };

const double adjust_gamma_up = 0.95;           /*!< Normalized value to which the maximum value present in the background image will be maped for type of gamma transformation Up */
const double adjust_gamma_down = 0.75;         /*!< Normalized value to which the maximum value present in the background image will be maped for type of gamma transformation Down */
/**
 @}
*/

template <class TPixInput, class TPixOutput>
class GammaTransf;

using GammaFilterType = itk::UnaryFunctorImageFilter<ImageBinType, ImageBinType, GammaTransf<BinPixelType, BinPixelType>>;

using ImageCalculatorFilterType = itk::MinimumMaximumImageCalculator<ImageBinType>;

#endif
