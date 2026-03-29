#ifndef _IMPARAMS_H
#define _IMPARAMS_H

#include <iostream>
#include <vector>

#include "colors.h"

/**
 @brief Structure to account for all the parameters needed by the buildim program to generate R data and/or an image. The expression cfpar means to which parameter of the configuration file this field corresponds.
*/
typedef struct
{
 float osp[2];                            /*!< Horizontal and vertical spacing in micrometers of the original image. cfpar: Spacing */
 int ppbin;                               /*!< Number of original image pixels per bin (a bin is a pixel in the generated pseudoimage). cfpar: PixPerBin */
 unsigned long w;                         /*!< Width of the original image in pixels. */
 unsigned long h;                         /*!< Height of the original image in pixels. cpfar: Dimensions */
 bool exfov;                              /*!< True if area of interest is described as a region identifier; false if it is as a set of two points, upper left/lower right. */
 std::string idfov;                      /*!< Identifier of the region of interest, if exfov is true. Empty string otherwise. cfpar AoIID */
 float fovwin[4];                         /*!< (xu,yu) (xl,yl) Coordinates of upper left and lower right of the area of interest if exfov is false. All 0.0 otherwise. cfpar: AoIWin */
 bool onlytargets;                        /*!< True if only target genes are allowed in the list of requested genes, false otherwise. cfpar: OnlyTargets */
 int qvthres;                             /*!< Threshold on the quality value to accept a transcription event. cfpar: QvThres */
 std::vector<int> numcodes;              /*!< Vector of codewords of requested genes, if given by codeword; empty vector if not. cfpar: TargetIds */
 std::vector<std::string> targetnames; /*!< Vector of short gene names of requested genes, if given by name; empty vector if not. cfpar: TargetNames */
 bool gmarks;                              /*!< True if request to generate marks in the R files, false otherwise. cfpar: GenMarks */
 std::vector<ColortabEntry> colortab;    /*!< Vector of colors assigned to each gene, as a ColortabEntry. Empty vector if no image is to be generated. cfpar: ColorFile */
 bool checkandexit;                       /*!< Parameters for internal use by the programmer. As true, the configuration file is parsed and checked and the program exits. */
} ImParams;

#endif
