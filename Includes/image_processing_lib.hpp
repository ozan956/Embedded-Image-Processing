#ifndef IMAGE_PROCESS_H
#define IMAGE_PROCESS_H

#include <iostream>
#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <random>
#include <complex>
#include <string>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>

const double PI = 3.14159265358979323846;
/**
 * @brief Using as PI Number.
 */
const double M_PI = 3.14159265358979323846;


/**
 * @brief Critical Image Properties Definitions
 */
#define WIDTH  50
#define HEIGHT  50
#define IMAGE_SIZE WIDTH*HEIGHT
#define MAX_INTENSITY 255

/**
 * @brief Labeling Definitions
 */
#define UNLABELED  -1;

/**
 * @brief Spatial Filter Function Parameter Definitions
 */
#define LINEAR_FILTER 0
#define NON_LINEAR_FILTER 1
#define LOW_PASS 0
#define HIGH_PASS 1
#define BAND_PASS 2
#define MEDIAN 3
#define MINIMUM 4
#define MAXIMUM 5


/**
 * @defgroup MacroDefinitions Macro Definitions & Declarations
 * @{
 */

  /**
   * @brief Returns the array based position of pixel.
   */
#define INDEX_TO_OFFSET(x, y) ((y) * WIDTH + (x))

/** @} */ // End of MacroDefinitions group




/**
 * @defgroup IntensityTransformations Intensity Transformations
 * @{
 */

 /**
  * @brief Negate the pixel intensities of an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data with negated intensities.
  */
void IMAGE_Negative(const uint8_t* imgData, uint8_t* outData);
/**
 * @brief Apply power law transformation to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with power-law transformed intensities.
 * @param[in] gamma The power-law transformation parameter.
 */
void IMAGE_powerTransform(const uint8_t* imgData, uint8_t* outData, uint8_t gamma);

/** @} */ // End of IntensityTransformations group

/**
 * @defgroup HistogramOperations Histogram Operations
 * @{
 */

 /**
  * @brief Calculate the histogram of an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] histogram Array to store the histogram.
  * @return Total number of pixels in the image.
  */
int IMAGE_histogramFormation(const uint8_t* imgData, int* histogram);

/**
 * @brief Perform histogram equalization on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with equalized intensities.
 */
void IMAGE_histogramEqualization(uint8_t* imgData, uint8_t* outData);

/**
 * @brief Perform histogram specification on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with specified histogram.
 * @param[in] targetHistogram The target histogram to match.
 */
void IMAGE_histogramSpecification(const uint8_t* imgData, uint8_t* outData, const int* targetHistogram);

/** @} */ // End of HistogramOperations group

/**
 * @defgroup Watermarking Digital Watermarking
 * @{
 */

 /**
  * @brief Embed a grayscale watermark into an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data with embedded watermark.
  * @param[in] watermark The watermark to embed.
  * @param[in] x X-coordinate of the watermark insertion point.
  * @param[in] y Y-coordinate of the watermark insertion point.
  */
void IMAGE_grayscaleWatermark(std::vector<uint8_t>& imgData, std::vector<uint8_t>& outData, const std::string& watermark, int x, int y);

/** @} */ // End of Watermarking group

/**
 * @defgroup SpatialFilters Spatial Filters
 * @{
 */

 /**
  * @brief Apply spatial filtering to an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data after spatial filtering.
  * @param[in] filterType Type of filter (LINEAR_FILTER or NON_LINEAR_FILTER).
  * @param[in] filterName Name of the filter (LOW_PASS, HIGH_PASS, etc.).
  */
void IMAGE_spatialFilter(const uint8_t* imgData, uint8_t* outData, int filterType, int filterName);

/** @} */ // End of SpatialFilters group

/**
 * @defgroup DFT 2D Discrete Fourier Transform (DFT)
 * @{
 */



void manualNormalize(cv::Mat& src, cv::Mat& dst);

void fft(std::vector<std::complex<double>>& data);

void fft2D(std::vector<std::vector<std::complex<double>>>& image);

void ifft2D(std::vector<std::vector<std::complex<double>>>& spectrum, cv::Mat& reconstructedImage);

void IMAGE_idealBPFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double minFrequency, double maxFrequency);

void IMAGE_idealFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double cutoffFrequency, int type);

void IMAGE_gaussianFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double cutoffFrequency, int filterType);



/** @} */ // End of DFT group

/**
 * @defgroup FrequencyFilters Frequency Domain Filters
 * @{
 */

 // ... Add functions for frequency domain filtering here ...

 /** @} */ // End of FrequencyFilters group

 /**
  * @defgroup Thresholding Thresholding and Segmentation
  * @{
  */

  /**
   * @brief Apply grayscale thresholding to an image.
   *
   * @param[in] imgData Input image data.
   * @param[out] outData Output binary image after thresholding.
   * @param[in] threshold Threshold value.
   */
void IMAGE_grayscaleThreshold(const uint8_t* imgData, uint8_t* outData, int threshold);

/**
 * @brief Calculate the Otsu's threshold for binarization.
 *
 * @param[in] histogram Histogram of the image.
 * @param[in] totalPixels Total number of pixels in the image.
 * @return Calculated Otsu's threshold.
 */
int IMAGE_otsuThreshold(const int* histogram, int totalPixels);

/**
 * @brief Apply Otsu's thresholding to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output binary image after Otsu's thresholding.
 */
void IMAGE_grayscaleOutsu(const uint8_t* imgData, uint8_t* outData);

/**
 * @brief Apply Sobel edge detection to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output edge-detected image.
 */
void IMAGE_sobelFilter(const uint8_t* imgData, uint8_t* outData);

/**
 * @brief Perform region growing on an image from a seed point.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output segmented image.
 * @param[in] seedX X-coordinate of the seed point.
 * @param[in] seedY Y-coordinate of the seed point.
 * @param[in] THRESHOLD Threshold for region growing.
 */
void regionGrowing(const uint8_t* imgData, uint8_t* outData, int seedX, int seedY, int THRESHOLD);



void IMAGE_basicSegmentation(const uint8_t* imgData, uint8_t* outData, uint8_t targetValue);

/** @} */ // End of Thresholding group

/**
 * @defgroup MorphologicalOps Morphological Image Processing
 * @{
 */
void IMAGE_erosion(const uint8_t* imgData, uint8_t* outData, int kernelSize);

void IMAGE_dilation(const uint8_t* imgData, uint8_t* outData, int kernelSize);

void IMAGE_opening(const uint8_t* imgData, uint8_t* outData, int kernelSize);

void IMAGE_closing(const uint8_t* imgData, uint8_t* outData, int kernelSize);

void labelConnectedComponents(const uint8_t* imgData, uint8_t* labeledData);

 /** @} */ // End of MorphologicalOps group

#endif // IMAGE_PROCESS_H
