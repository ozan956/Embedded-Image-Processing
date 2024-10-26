#ifndef DIP_H
#define DIP_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <limits.h>
#include <float.h>
#include "esp_log.h"

#include "esp_dsp.h"
#include <math.h>

#include "dsps_fft2r.h"

extern "C"

#define IMAGE_FORMAT_RGB565		1
#define IMAGE_FORMAT_GRAYSCALE	3

#define IMAGE_RES_QQVGA			0
#define IMAGE_RES_QVGA			1
#define IMAGE_RES_480x272		2
#define IMAGE_RES_VGA   		3

#define IMAGE_RES_VGA_Width			640
#define IMAGE_RES_VGA_Height		480
#define IMAGE_RES_QVGA_Width		320
#define IMAGE_RES_QVGA_Height		240
#define IMAGE_RES_QQVGA_Width		160
#define IMAGE_RES_QQVGA_Height		120
#define IMAGE_RES_480x272_Width		480
#define IMAGE_RES_480x272_Height	270

#define MAX_INTENSITY  255

/**
 * @brief Histogram comparison method definitions.
 */
#define HIST_COMP_CORRELATION           0
#define HIST_COMP_CHI_SQUARE            1
#define HIST_COMP_INTERSECTION          2
#define HIST_COMP_BHATTACHARYYA         3

/**
 * @brief Mathematical constants.
 */
#define pi                              acos(-1)

/**
 * @brief Thresholding types.
 */
#define THRESH_BINARY                   0
#define THRESH_BINARY_INV               1

/**
 * @brief Default border type.
 */
#define BORDER_DEFAULT                  0


// Macro to get the image size based on width
#define GET_IMAGE_SIZE(width) \
    ((width) == IMAGE_RES_QQVGA_Width ? IMAGE_RES_QQVGA : \
     (width) == IMAGE_RES_QVGA_Width ? IMAGE_RES_QVGA : \
     (width) == IMAGE_RES_480x272_Width ? IMAGE_RES_480x272 : \
     (width) == IMAGE_RES_VGA_Width ? IMAGE_RES_VGA : \
     -1)  // Return -1 if width does not match any known size


// Define a structure to represent an image
typedef struct {
    uint32_t width;        // Width of the image in pixels
    uint32_t height;       // Height of the image in pixels
    uint8_t *pixels;       // Pointer to the pixel data
    uint32_t size;         // Size of the image data in bytes
    uint8_t format;      // Number of color channels (e.g., 3 for RGB, 4 for RGBA)

    float **fourier;
}Image;

// Structure to hold filter kernel
typedef struct {
    int size;
    float **data;
    uint8_t *data_8;
} FilterKernel;

/**
 * @brief Structure to represent size in 2D space.
 */
typedef struct {
    int width;           ///< Width (number of columns).
    int height;          ///< Height (number of rows).
} Size;


/**
 * @brief Creates an image with a specified size and format.
 *
 * @param size Size of the image.
 * @param format Format of the image (e.g., grayscale, RGB).
 * @return Pointer to the created Image.
 */
Image* createImage(uint8_t size, uint8_t format);

/**
 * @brief Converts an input image to grayscale or another color format.
 *
 * This function takes an input image and converts it to either
 * grayscale or another specified color format based on the
 * toGrayscale parameter.
 *
 * @param input Pointer to the input Image structure, containing pixel data
 *              and image dimensions.
 * @param output Pointer to the output Image structure, where the converted
 *               image will be stored.
 * @param toGrayscale A flag indicating the conversion type:
 *                    - 1 (true) to convert to grayscale
 *                    - 0 (false) to perform no conversion (copy input to output)
 */
void cvtColor(Image *input, Image *output, uint8_t toGrayscale);

/**
 * @brief Resizes an input image to a specified size using nearest neighbor interpolation.
 *
 * This function takes an input image and resizes it to a square output image
 * of size `fourier_size` by calculating the nearest neighbor pixel values
 * based on the width and height ratios of the original image to the target size.
 *
 * @param inImg Pointer to the input Image structure containing the original image data.
 * @param outImg Pointer to the output Image structure where the resized image will be stored.
 * @param fourier_size The desired size of the output image (both width and height).
 */
void DIP_Resize(Image *inImg, Image *outImg, int fourier_size);

/**
 * @brief Performs a negative transformation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (negative transformed).
 */
void DIP_negative(const Image *inImg, Image *outImg);

/**
 * @brief Applies a power-law (gamma) transformation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (gamma transformed).
 * @param gamma Gamma value for transformation.
 */
void DIP_powerTransform(const Image *inImg, Image *outImg, uint8_t gamma);


/**
 * @brief Forms a histogram from the grayscale values of the input image.
 *
 * @param inImg Input image.
 * @param histogram Array to store the histogram (256 elements for 8-bit images).
 * @return 0 on success, non-zero on error.
 */
int DIP_histogramFormation(const Image *inImg, int *histogram);

/**
 * @brief Performs histogram equalization on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (equalized).
 */
void DIP_histogramEqualization(const Image *inImg, Image *outImg);

/**
 * @brief Modifies the histogram of the input image to match a target histogram.
 *
 * @param inImg Input image.
 * @param outImg Output image with histogram specification.
 * @param targetHistogram Target histogram array.
 */
void DIP_histogramSpecification(const Image *inImg, Image *outImg,
		const int *targetHistogram);

/**
 * @brief Compare two histograms.
 *
 * @param[in] hist1 First histogram (array of size MAX_INTENSITY + 1).
 * @param[in] hist2 Second histogram (array of size MAX_INTENSITY + 1).
 * @param[in] method Comparison method (e.g., HIST_COMP_CORRELATION, HIST_COMP_CHI_SQUARE).
 * @return Comparison result as a double.
 */
double DIP_compareHist(int *hist1, int *hist2, int method);

/**
 * @brief Blends two images using linear interpolation.
 *
 * This function takes a base image and an overlay image, and blends them
 * together using a specified alpha value. The result is stored in the
 * output image. The alpha value determines the transparency of the
 * overlay image: an alpha of 0 means the overlay is fully transparent,
 * and an alpha of 1 means the overlay is fully opaque.
 *
 * @param baseImg Pointer to the base Image structure, which serves as the background.
 * @param overlayImg Pointer to the overlay Image structure, which is blended on top of the base image.
 * @param outImg Pointer to the output Image structure, where the blended result will be stored.
 * @param alpha A float value representing the blend ratio, in the range [0, 1].
 */
void DIP_blendLinear(const Image *baseImg, const Image *overlayImg, Image *outImg, float alpha);

/**
 * @brief Embeds a grayscale watermark into the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image with watermark.
 * @param watermark String representing the watermark.
 * @param x X coordinate for watermark position.
 * @param y Y coordinate for watermark position.
 */
void DIP_grayscaleWatermark(const Image *inImg, Image *outImg,
		const char *watermark, int x, int y);

/**
 * @brief Applies a 2D filter to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (filtered).
 * @param size Size of the filter matrix (e.g., 3x3, 5x5).
 * @param filter 2D array representing the filter.
 */
void DIP_filter2D(const Image *inImg, Image *outImg, int size,
		const float *filter);

/**
 * @brief Applies a separable 2D filter to an input image using two kernels.
 *
 * This function performs a two-pass filtering operation on an input image.
 * The first pass applies a horizontal kernel to filter the rows of the image,
 * storing the intermediate results in a temporary image. The second pass
 * applies a vertical kernel to the temporary image, producing the final
 * filtered output image. A delta value is added after filtering to adjust the result.
 *
 * @param inImg Pointer to the input Image structure containing the original image data.
 * @param outImg Pointer to the output Image structure where the filtered result will be stored.
 * @param kernelSizeX Size of the horizontal kernel.
 * @param kernelX Pointer to the horizontal kernel array for filtering.
 * @param kernelSizeY Size of the vertical kernel.
 * @param kernelY Pointer to the vertical kernel array for filtering.
 * @param delta A double value to be added to the filtered result before normalization.
 */
void DIP_sepFilter2D(const Image *inImg, Image *outImg, int kernelSizeX, float *kernelX,
                     int kernelSizeY, float *kernelY, double delta);

/**
 * @brief Applies a square box filter to an input image, calculating the sum of squared pixel values
 * within a specified filter size neighborhood.
 *
 * This function processes each pixel in the input image by applying a square box filter, which computes
 * the sum of the squares of the pixel values in the neighborhood defined by the filter size. The result
 * can be normalized to fit within the 8-bit range [0, 255].
 *
 * @param inImg Pointer to the input Image structure containing the original image data.
 * @param outImg Pointer to the output Image structure where the filtered result will be stored.
 * @param filtersize Size of the square filter (must be an odd integer).
 * @param normalize Integer flag to indicate whether to normalize the result (1 for true, 0 for false).
 */
void DIP_sqrBoxFilter(const Image *inImg, Image *outImg, int filtersize, int normalize);

/**
 * @brief Applies a Laplacian filter to the input image to detect edges.
 *
 * This function computes the Laplacian of the input image using a specified kernel size,
 * scaling factor, and delta value for result adjustment. The filter enhances edges by
 * calculating the second derivative of the image intensity.
 *
 * @param inImg Pointer to the input Image structure containing the original image data.
 * @param outImg Pointer to the output Image structure where the filtered result will be stored.
 * @param ddepth Desired depth of the output image (not used in this implementation).
 * @param ksize Size of the Laplacian kernel (must be positive and odd).
 * @param scale Scaling factor applied to the computed values.
 * @param delta Value added to the filtered result before clamping.
 * @param borderType Type of border handling (e.g., BORDER_DEFAULT).
 */
void DIP_Laplacian(const Image *inImg, Image *outImg, int ddepth, int ksize, double scale, double delta, int borderType);

/**
 * @brief Applies a Gaussian blur to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (blurred).
 * @param filtersize Size of the Gaussian filter.
 * @param sigma Standard deviation of the Gaussian distribution.
 */
void DIP_gaussianBlur(const Image *inImg, Image *outImg, int filtersize,
		int sigma);

/**
 * @brief Applies a median blur to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (blurred).
 * @param filtersize Size of the filter.
 */
void DIP_medianBlur(const Image *inImg, Image *outImg, int filtersize);

/**
 * @brief Performs a downsampling operation on the source image using Gaussian blur.
 *
 * This function first applies a Gaussian blur to the input image to reduce aliasing,
 * then downsamples the blurred image by selecting every second pixel to create a smaller output image.
 *
 * @param src Pointer to the source Image structure (input image).
 * @param dst Pointer to the destination Image structure (output image).
 * @param dstsize Pointer to a Size structure specifying the desired output size (can be NULL).
 * @param borderType Type of border handling (not used in this implementation).
 */
void DIP_pyrDown(const Image *src, Image *dst, const Size *dstsize, int borderType);

/**
 * @brief Performs an upsampling operation on the input image.
 *
 * This function first upsamples the input image by inserting zeros between pixels
 * to create a larger image, and then applies Gaussian blur to smooth the upsampled image.
 *
 * @param inImg Pointer to the input Image structure (input image).
 * @param outImg Pointer to the output Image structure (filtered output image).
 * @param dstSize Pointer to a Size structure specifying the desired output size.
 * @param borderType Type of border handling (not used in this implementation).
 * @param filtersize Size of the Gaussian filter to apply after upsampling.
 * @param sigma Standard deviation for the Gaussian filter.
 */
void DIP_pyrUp(const Image *inImg, Image *outImg, Size *dstSize, int borderType, int filtersize, int sigma);

/**
 * @brief Computes the Fourier transform of the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (Fourier transformed).
 */
void DIP_fourier(Image *inImg, Image *outImg, int FFT_SIZE);

/**
 * @brief Gets the magnitude spectrum from a Fourier transformed image.
 *
 * @param inImg Input image (Fourier domain).
 * @param outImg Output image (magnitude spectrum).
 * @param fftshift 1 to apply fftshift, 0 otherwise.
 */
void DIP_abs(Image *inImg, Image *outImg, int FFT_SIZE, int quadrant_change);

/**
 * @brief Computes the inverse Fourier transform of the input image.
 *
 * @param inImg Input image (Fourier domain).
 * @param outImg Output image (spatial domain).
 */
void DIP_fourier_inv(Image *inImg, Image *outImg, int FFT_SIZE);

/**
 * @brief Apply grayscale thresholding to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output binary image after thresholding.
 * @param[in] threshold Threshold value.
 */
void DIP_grayscaleThreshold(const Image *inImg, Image *outImg, int threshold);

/**
 * @brief Apply adaptive grayscale thresholding to an image.
 *
 * @param[in] inImg Input image data.
 * @param[out] outImg Output binary image after thresholding.
 * @param[in] maxValue Non-zero value assigned to the pixels for which the condition is satisfied.
 * @param[in] adaptiveMethod Adaptive thresholding algorithm (currently only mean is implemented).
 * @param[in] thresholdType Type of thresholding (THRESH_BINARY or THRESH_BINARY_INV).
 * @param[in] blockSize Size of a pixel neighborhood (must be odd).
 * @param[in] C Constant subtracted from the mean.
 */
void DIP_adaptiveThreshold(const Image *inImg, Image *outImg, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C);

/**
 * @brief Applies Otsu's method to perform automatic thresholding.
 *
 * @param inImg Input image.
 * @param outImg Output image (binary result).
 */
void DIP_otsuMethod(const Image *inImg, Image *outImg);

/**
 * @brief Apply Sobel edge detection to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output edge-detected image.
 */
void DIP_sobelFilter(const Image *inImg, Image *outImg);

/**
 * @brief Performs region growing segmentation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (segmented).
 * @param y Starting Y coordinate for region growing.
 * @param x Starting X coordinate for region growing.
 * @param threshold Intensity threshold for region growing.
 */
void regionGrowing(Image *inImg, Image *outImg, int y, int x, int threshold);

/**
 * @brief Performs basic segmentation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (segmented).
 * @param targetValue Value used for segmentation.
 */
void DIP_basicSegmentation(const Image *inImg, Image *outImg,
		uint8_t targetValue);

/**
 * @brief Performs image erosion using a specified kernel size.
 *
 * @param inImg Input image.
 * @param outImg Output image (eroded).
 * @param kernelSize Size of the erosion kernel.
 */
void DIP_erode(const Image *inImg, Image *outImg, int kernelSize);

/**
 * @brief Performs image dilation using a specified kernel size.
 *
 * @param inImg Input image.
 * @param outImg Output image (dilated).
 * @param kernelSize Size of the dilation kernel.
 */
void DIP_dilate(const Image *inImg, Image *outImg, int kernelSize);

/**
 * @brief Performs morphological opening (erosion followed by dilation) on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (opened).
 * @param kernelSize Size of the kernel.
 */
void DIP_opening(const Image *inImg, Image *outImg, int kernelSize);

/**
 * @brief Performs morphological closing (dilation followed by erosion) on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (closed).
 * @param kernelSize Size of the kernel.
 */
void DIP_closing(const Image *inImg, Image *outImg, int kernelSize);

/**
 * @brief Identifies and labels connected components in a binary image.
 *
 * @param inImg Input image.
 * @param outImg Output image (labeled components).
 */
void DIP_connectedComponents(const Image *inImg, Image *outImg);

/**
 * @brief Applies the Canny edge detection algorithm to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (edges).
 * @param lowThreshold Low threshold for edge detection.
 * @param HighThreshold High threshold for edge detection.
 */
void DIP_canny(const Image *inImg, Image *outImg, int lowThreshold,
		int HighThreshold);

/**
 * @brief Translates the input image by a specified distance in the X and Y directions.
 *
 * @param inImg Input image.
 * @param outImg Output image (translated).
 * @param dx Translation distance along the X axis.
 * @param dy Translation distance along the Y axis.
 */
void translate(const Image *inImg, Image *outImg, int dx, int dy);

/**
 * @brief Scales the input image by specified factors in the X and Y directions.
 *
 * @param inImg Input image.
 * @param outImg Output image (scaled).
 * @param scaleX Scaling factor along the X axis.
 * @param scaleY Scaling factor along the Y axis.
 */
void DIP_scale(const Image *inImg, Image *outImg, float scaleX, float scaleY);

/**
 * @brief Rotates the input image by a specified angle.
 *
 * @param inImg Input image.
 * @param outImg Output image (rotated).
 * @param angle Rotation angle in degrees.
 */
void DIP_rotate(const Image *inImg, Image *outImg, float angle);

/**
 * @brief Shears the input image by specified factors along the X and Y axes.
 *
 * @param inImg Input image.
 * @param outImg Output image (sheared).
 * @param shx Shear factor along the X axis.
 * @param shy Shear factor along the Y axis.
 */
void DIP_shear(const Image *inImg, Image *outImg, float shx, float shy);

/**
 * @brief Inverse shears the input image by specified factors along the X and Y axes.
 *
 * @param inImg Input image.
 * @param outImg Output image (sheared).
 * @param shx Shear factor along the X axis.
 * @param shy Shear factor along the Y axis.
 */
void DIP_inverse_shear(const Image *inImg, Image *outImg, float shx, float shy);

/**
 * @brief Applies a perspective warp transformation to an input image.
 *
 * This function takes an input image and applies a perspective warp
 * transformation using a given 3x3 transformation matrix. The output
 * image is initialized to black (or a specified border value).
 *
 * @param[in] inImg Pointer to the input image structure.
 *                  Must not be NULL and should contain valid pixel data.
 * @param[out] outImg Pointer to the output image structure.
 *                    Must not be NULL and should be allocated with sufficient size.
 * @param[in] M 3x3 transformation matrix used for perspective warping.
 *              This matrix defines how the input image is transformed
 *              into the output image.
 *
 * @note The output image dimensions should be set prior to calling this function.
 * @warning Ensure that the input and output image pixel data is properly allocated
 *          before invoking this function. If the dimensions of the output image
 *          are larger than the input image, the additional pixels will remain black.
 */
void DIP_warpPerspective(const Image *inImg, Image *outImg, float M[3][3]);

/**
 * @brief Performs a perspective transformation on an input image.
 *
 * This function takes an input image and applies a perspective transformation
 * defined by a 3x3 transformation matrix. The output image is initialized to
 * a default value (e.g., black or transparent) before the transformation is applied.
 *
 * @param[in] inImg Pointer to the input image structure.
 *                  Must not be NULL and should contain valid pixel data.
 * @param[out] outImg Pointer to the output image structure.
 *                    Must not be NULL and should be allocated with sufficient size.
 * @param[in] h 3x3 transformation matrix used for perspective transformation.
 *              This matrix defines how each pixel in the input image is mapped
 *              to the output image.
 *
 * @note The dimensions of the output image should be set prior to calling this function.
 * @warning Ensure that the input and output image pixel data is properly allocated
 *          before invoking this function. If the new coordinates exceed the output
 *          image dimensions, the corresponding pixels will be ignored.
 *
 * @sa Image
 */
void DIP_perspectiveTransform(const Image *inImg, Image *outImg, float h[3][3]);

/**
 * @brief Performs template matching using normalized cross-correlation.
 *
 * This function applies template matching on a source image using a given
 * template image and outputs the results into a result image. The method
 * utilizes normalized cross-correlation for the matching process. The result
 * image will highlight areas where the template matches the source image
 * above a certain correlation threshold.
 *
 * @param[in] sourceImg Pointer to the source image structure.
 *                      Must not be NULL and should contain valid pixel data.
 * @param[in] templateImg Pointer to the template image structure.
 *                        Must not be NULL and should contain valid pixel data.
 * @param[out] resultImg Pointer to the result image structure.
 *                       Must not be NULL and should be allocated with sufficient size.
 *
 * @note The result image will be modified to highlight matched regions with
 *       the same pixel values as in the source image. The rest of the pixels
 *       will remain unchanged or be initialized as per the result image's
 *       original state.
 *
 * @warning Ensure that the input and output image pixel data is properly allocated
 *          before invoking this function. If the template size exceeds the source
 *          image dimensions, no matches will be found.
 */
void DIP_TemplateMatch(const Image *sourceImg, const Image *templateImg, Image *resultImg);

// DIP UTIL ************************************************************************************

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
void DIP_UTIL_sobelFilter(const Image *inImg, Image *angle, float *fourierx, float *fouriery);

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
void DIP_UTIL_nonMaxSuppression(const Image *inImg, Image *nonmax, float *GradX, float *GradY);

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

#endif /* DIP_H_ */
