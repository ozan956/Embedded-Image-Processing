/*
 * DIP_UTIL.h
 *
 *  Created on: May 4, 2024
 *      Author: ozand
 */

#ifndef INC_DIP_UTIL_H_
#define INC_DIP_UTIL_H_
#include "DIP.h"

/**
 * @brief Applies Sobel operator with mirror padding to the input image.
 *
 * This utility function calculates the Sobel filter values at a specific pixel
 * location (y, x) using the specified Sobel filter coefficients (fx, fy).
 * The function takes into account the mirror padding for edge pixels.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] pixelValue Pointer to the integer where the calculated pixel value will be stored.
 * @param[in] y Y-coordinate of the pixel to process.
 * @param[in] x X-coordinate of the pixel to process.
 * @param[in] fx Filter coefficient for the x-direction.
 * @param[in] fy Filter coefficient for the y-direction.
 */
void sobelUtilizedMirrorPadding(const Image *inImg, int *pixelValue, int y,
		int x, int fx, int fy);

/**
 * @brief Applies Sobel filter to the input image and calculates the gradient angle.
 *
 * This function computes the Sobel filter on the input image and outputs the
 * gradient angle in the provided angle image.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] angle Pointer to the output image structure where gradient angles will be stored.
 */
void DIP_UTIL_sobelFilter(const Image *inImg, Image *angle);

/**
 * @brief Applies mirror padding to the input image.
 *
 * This function modifies the input image by applying mirror padding at the
 * specified pixel location (i, j). The function stores the padded values in
 * the provided float arrays (r and q) based on the given condition.
 *
 * @param[in,out] inImg Pointer to the input image structure.
 * @param[in] i Row index for the pixel.
 * @param[in] j Column index for the pixel.
 * @param[out] r Pointer to the float array for storing one set of padded values.
 * @param[out] q Pointer to the float array for storing another set of padded values.
 * @param[in] condition Condition to determine the padding behavior.
 */
void mirrorPadding(Image *inImg, int i, int j, float *r, float *q,
		int condition);

/**
 * @brief Performs non-maximum suppression on the input image.
 *
 * This function applies non-maximum suppression to the input image to
 * thin out the detected edges. The result is stored in the nonmax image.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] nonmax Pointer to the output image structure for storing non-max suppressed results.
 */
void DIP_UTIL_nonMaxSupression(Image *inImg, Image *nonmax);

/**
 * @brief Applies double thresholding to the input image.
 *
 * This function performs double thresholding on the input image to classify
 * the pixels into strong, weak, or non-edges based on the specified thresholds.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[in] lowThresholdRatio Ratio for the low threshold.
 * @param[in] highThresholdRatio Ratio for the high threshold.
 * @param[out] outImg Pointer to the output image structure for storing the result.
 * @param[in] weak Value representing weak edge pixels.
 * @param[in] strong Value representing strong edge pixels.
 */
void DIP_UTIL_doubleTreshold(Image *inImg, float lowThresholdRatio,
		float highThresholdRatio, Image *outImg, int weak, int strong);

/**
 * @brief Applies hysteresis thresholding to the input image.
 *
 * This function performs hysteresis thresholding on the input image,
 * linking weak edge pixels to strong edge pixels to form continuous edges.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] outImg Pointer to the output image structure for storing the result.
 * @param[in] weak Value representing weak edge pixels.
 * @param[in] strong Value representing strong edge pixels.
 */
void DIP_UTIL_hysteresis(Image *inImg, Image *outImg, int weak, int strong);

/**
 * @brief Compares two unsigned 8-bit integers for qsort.
 *
 * This utility function is used for comparing two unsigned 8-bit integers
 * for sorting purposes.
 *
 * @param[in] a Pointer to the first unsigned 8-bit integer.
 * @param[in] b Pointer to the second unsigned 8-bit integer.
 * @return An integer less than, equal to, or greater than zero if the first
 *         argument is considered to be respectively less than, equal to, or
 *         greater than the second.
 */
int compare_uint8(const void *a, const void *b);

/**
 * @brief Computes the largest power of 2 less than or equal to the given value.
 *
 * This function takes an unsigned integer as input and returns the largest power of 2
 * that is less than or equal to that integer. If the input value is zero, the function
 * will return 0 since there is no power of 2 that satisfies this condition.
 *
 * @param[in] value The unsigned integer value for which to compute the power of 2.
 * @return The largest power of 2 less than or equal to the input value.
 *         Returns 0 if the input value is 0.
 *
 * @note This function utilizes the __builtin_clz function to count the number of leading
 *       zeros in the input value. It assumes that the input value is a non-negative integer.
 */
unsigned int orderof2(const unsigned int value);

/**
 * @brief Generates a Gabor kernel.
 *
 * This function creates a Gabor filter kernel based on the specified parameters.
 *
 * @param[out] filter 2D array to store the generated Gabor filter.
 * @param[in] size Size of the filter (should be odd).
 * @param[in] sigma Standard deviation of the Gaussian envelope.
 * @param[in] theta Orientation of the filter in radians.
 * @param[in] lambda Wavelength of the cosine factor.
 * @param[in] gamma Spatial aspect ratio.
 * @param[in] psi Phase offset.
 */
void getGaborKernel(float filter[][3], int size, double sigma, double theta, double lambda, double gamma, double psi);

/**
 * @brief Generates a Gaussian kernel.
 *
 * This function creates a Gaussian filter kernel based on the specified
 * parameters.
 *
 * @param[out] filter 2D array to store the generated Gaussian filter.
 * @param[in] size Size of the filter (should be odd).
 * @param[in] sigma Standard deviation of the Gaussian.
 */
void getGaussianKernel(float filter[][3], int size, double sigma);

/**
 * @brief Creates a Laplacian filter.
 *
 * This function generates a Laplacian filter of the specified size.
 *
 * @param[in] ksize Size of the filter.
 * @return A pointer to a 2D array representing the Laplacian filter.
 */
float **createLaplacianFilter(int ksize);

#endif /* INC_DIP_UTIL_H_ */
