/*
 * DIP_Driver.c
 *
 *  Created on: Sep 19, 2024
 *      Author: ozand
 */

#include "DIP_Driver.h"

// Initialize the DIP Manager
void initDIPManager(DIP *dip) {

    // Pixel-Based Operations
    dip->histForm = &DIP_histogramFormation;
    dip->histEq = &DIP_histogramEqualization;
    dip->histSpec = &DIP_histogramSpecification;
    dip->histCompare = &DIP_compareHist;
    dip->negative = &DIP_negative;
    dip->powerTransform = &DIP_powerTransform;
    dip->grayscaleWatermark = &DIP_grayscaleWatermark;

    // Filtering Image in Spatial Domain
    dip->filter2D = &DIP_filter2D;
    dip->gaussianBlur = &DIP_gaussianBlur;
    dip->medianBlur = &DIP_medianBlur;

    // Filtering in Frequency Domain
    dip->fourier = &DIP_fourier;
    dip->abs = &DIP_abs;
    dip->fourierInv = &DIP_fourierInv;

    // From Pixels to Objects
    dip->sobelFilter = &DIP_sobelFilter;
    dip->canny = &DIP_canny;
    dip->regionGrowing = &DIP_regionGrowing;
    dip->grayscaleThreshold = &DIP_grayscaleThreshold;
    dip->otsuMethod = &DIP_otsuMethod;
    dip->basicSegmentation = &DIP_basicSegmentation;

    // Morphological Image Processing
    dip->dilate = &DIP_dilate;
    dip->erode = &DIP_erode;
    dip->closing = &DIP_closing;
    dip->opening = &DIP_opening;
    dip->connectedComponents = &DIP_connectedComponents;

    // Geometrical Transformations
    dip->translate = &DIP_translate;
    dip->scale = &DIP_scale;
    dip->rotate = &DIP_rotate;
    dip->shear = &DIP_shear;

    //Misc
    dip->blenLinear = &DIP_blendLinear;
    dip->templateMatching = &DIP_template_matching;
    dip->laplacian = &DIP_Laplacian;
    dip->pyrDown = &DIP_pyrDown;
    //dip->pyrUp = &DIP_pyrUp;
    dip->sqrBoxFilter = &DIP_sqrBoxFilter;
    dip->filter2Dsep = &DIP_sepFilter2D;
    dip->adaptiveThreshold = &DIP_adaptiveThreshold;
    dip->shearInv = &DIP_inverse_shear;
    dip->warpPerspective = &DIP_warpPerspective;
    dip->perspectiveTransform = &DIP_perspectiveTransform;
    dip->perspectiveTransform = &DIP_template_matching;

}
