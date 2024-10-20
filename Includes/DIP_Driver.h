/*
 * DIP_Driver.h
 *
 *  Created on: Sep 19, 2024
 *      Author: ozand
 */

#ifndef DIP_INCLUDE_DIP_DRIVER_H_
#define DIP_INCLUDE_DIP_DRIVER_H_

#include "DIP.h"
#include "DIP_LTDC.h"
#include "DIP_SERIAL.h"
#include "DIP_UTIL.h"

// Define the DIP Manager structure
typedef struct {

    // Pixel-Based Operations

        // Histogram of an Image

		/**
		 * @brief Calculate the histogram of an image.
		 *
		 * @param[in] imgData Input image data.
		 * @param[out] histogram Array to store the histogram.
		 * @return Total number of pixels in the image.
		 */
        int (*histForm)(const Image *inImg, int *histogram);

        /**
         * @brief Perform histogram equalization on an image.
         *
         * @param[in] imgData Input image data.
         * @param[out] outData Output image data with equalized intensities.
         */
        void (*histEq)(const Image *inImg, Image *outImg);

        /**
         * @brief Perform histogram specification on an image.
         *
         * @param[in] imgData Input image data.
         * @param[out] outData Output image data with specified histogram.
         * @param[in] targetHistogram The target histogram to match.
         */
        void (*histSpec)(const Image *inImg, Image *outImg, const int* targetHistogram);

        /**
         * @brief Compare two histograms.
         *
         * @param[in] hist1 First histogram (array of size MAX_INTENSITY + 1).
         * @param[in] hist2 Second histogram (array of size MAX_INTENSITY + 1).
         * @param[in] method Comparison method (e.g., HIST_COMP_CORRELATION, HIST_COMP_CHI_SQUARE).
         * @return Comparison result as a double.
         */
        double (*histCompare)(int *hist1, int *hist2, int method);

        // Negation Law

        /**
         * @brief Negate the pixel intensities of an image.
         *
         * @param[in] imgData Input image data.
         * @param[out] outData Output image data with negated intensities.
         */
        void (*negative)(const Image *inImg, Image *outImg);

        // Power Law

        /**
         * @brief Apply power law transformation to an image.
         *
         * @param[in] imgData Input image data.
         * @param[out] outData Output image data with power-law transformed intensities.
         * @param[in] gamma The power-law transformation parameter.
         */
        void (*powerTransform)(const Image *inImg, Image *outImg, uint8_t gamma);

        // Digital Watermarking

        /**
         * @brief Embeds a grayscale watermark into the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image with watermark.
         * @param watermark String representing the watermark.
         * @param x X coordinate for watermark position.
         * @param y Y coordinate for watermark position.
         */
        void (*grayscaleWatermark)(const Image *inImg, Image *outImg, const char* watermark, int x, int y);

        /**
         * @brief Blends two images linearly using a specified alpha value.
         *
         * @param baseImg Base image.
         * @param overlayImg Overlay image.
         * @param outImg Output image (blended).
         * @param alpha Alpha value for blending (0.0 - 1.0).
         */
        void (*blenLinear)(const Image *baseImg, const Image *overlayImg, Image *outImg, float alpha);

    // Filtering Image in Spatial Domain

        // Linear Filters

        /**
         * @brief Applies a 2D filter to the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (filtered).
         * @param size Size of the filter matrix (e.g., 3x3, 5x5).
         * @param filter 2D array representing the filter.
         */
        void (*filter2D)(const Image *inImg, Image *outImg, int size, float filter[size][size]);

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
        void (*filter2Dsep)(const Image *inImg, Image *outImg, int kernelSizeX, float kernelX[kernelSizeX],
                             int kernelSizeY, float kernelY[kernelSizeY], double delta);

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
        void (*sqrBoxFilter)(const Image *inImg, Image *outImg, int filtersize, int normalize);

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
        void (*laplacian)(const Image *inImg, Image *outImg, int ddepth, int ksize, double scale, double delta, int borderType);

        /**
         * @brief Applies a Gaussian blur to the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (blurred).
         * @param filtersize Size of the Gaussian filter.
         * @param sigma Standard deviation of the Gaussian distribution.
         */
        void (*gaussianBlur)(const Image *inImg, Image *outImg, int filtersize, int sigma);

        // Non-Linear Filters

        /**
         * @brief Applies a median blur to the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (blurred).
         * @param filtersize Size of the filter.
         */
        void (*medianBlur)(const Image *inImg, Image *outImg, int filtersize);

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
        void (*pyrDown)(const Image *src, Image *dst, const Size *dstsize, int borderType);

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
        void (*pyrUp)(const Image *inImg, Image *outImg, Size *dstSize, int borderType, int filtersize, int sigma);

    // Filtering in Frequency Domain

        /**
         * @brief Computes the Fourier transform of the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (Fourier transformed).
         */
        void (*fourier)(const Image *inImg);

        /**
         * @brief Takes the absolute value of all pixel values in the output image.
         *
         * @param outImg Output image to be processed.
         */
        void (*abs)(Image *outImg);

        /**
         * @brief Computes the inverse Fourier transform of the input image.
         *
         * @param inImg Input image (Fourier domain).
         * @param outImg Output image (spatial domain).
         */
        void (*fourierInv)(const Image *inImg, Image *outImg);

    // From Pixels to Objects

        // Edge Detection

        /**
         * @brief Apply Sobel edge detection to an image.
         *
         * @param[in] imgData Input image data.
         * @param[out] outData Output edge-detected image.
         */
        void (*sobelFilter)(const Image *inImg, Image *angle);

        /**
         * @brief Applies the Canny edge detection algorithm to the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (edges).
         * @param lowThreshold Low threshold for edge detection.
         * @param HighThreshold High threshold for edge detection.
         */
        void (*canny)(const Image *inImg, Image *outImg, int lowThreshold, int highThreshold);

        // Region Growing

        /**
         * @brief Performs region growing segmentation on the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (segmented).
         * @param y Starting Y coordinate for region growing.
         * @param x Starting X coordinate for region growing.
         * @param threshold Intensity threshold for region growing.
         */
        void (*regionGrowing)(Image *inImg, Image *outImg, int y, int x, int threshold);

        // Thresholding

        /**
         * @brief Applies a threshold to a grayscale image.
         *
         * @param inImg Input image.
         * @param outImg Output image (thresholded).
         * @param threshold Threshold value.
         */
        void (*grayscaleThreshold)(const Image *inImg, Image *outImg, int threshold);

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
        void (*adaptiveThreshold)(const Image *inImg, Image *outImg, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C);

        /**
         * @brief Applies Otsu's method to perform automatic thresholding.
         *
         * @param inImg Input image.
         * @param outImg Output image (binary result).
         */
        void (*otsuMethod)(const Image *inImg, Image *outImg);

        // Segmentation

        /**
         * @brief Performs basic segmentation on the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (segmented).
         * @param targetValue Value used for segmentation.
         */
        void (*basicSegmentation)(const Image *inImg, Image *outImg, uint8_t targetValue);

    // Morphological Image Processing

        // Dilation

        /**
         * @brief Performs image dilation using a specified kernel size.
         *
         * @param inImg Input image.
         * @param outImg Output image (dilated).
         * @param kernelSize Size of the dilation kernel.
         */
        void (*dilate)(const Image *inImg, Image *outImg, int kernelSize);

        // Erosion

        /**
         * @brief Performs image erosion using a specified kernel size.
         *
         * @param inImg Input image.
         * @param outImg Output image (eroded).
         * @param kernelSize Size of the erosion kernel.
         */
        void (*erode)(const Image *inImg, Image *outImg, int kernelSize);

        // Closing

        /**
         * @brief Performs morphological closing (dilation followed by erosion) on the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (closed).
         * @param kernelSize Size of the kernel.
         */
        void (*closing)(const Image *inImg, Image *outImg, int kernelSize);

        // Opening

        /**
         * @brief Performs morphological opening (erosion followed by dilation) on the input image.
         *
         * @param inImg Input image.
         * @param outImg Output image (opened).
         * @param kernelSize Size of the kernel.
         */
        void (*opening)(const Image *inImg, Image *outImg, int kernelSize);

        // Labeling

        /**
         * @brief Identifies and labels connected components in a binary image.
         *
         * @param inImg Input image.
         * @param outImg Output image (labeled components).
         */
        void (*connectedComponents)(const Image *inImg, Image *outImg);

    // Geometrical Transformations

        /**
         * @brief Translates the input image by a specified distance in the X and Y directions.
         *
         * @param inImg Input image.
         * @param outImg Output image (translated).
         * @param dx Translation distance along the X axis.
         * @param dy Translation distance along the Y axis.
         */
        void (*translate)(const Image *inImg, Image *outImg, int dx, int dy);

        /**
         * @brief Scales the input image by specified factors in the X and Y directions.
         *
         * @param inImg Input image.
         * @param outImg Output image (scaled).
         * @param scaleX Scaling factor along the X axis.
         * @param scaleY Scaling factor along the Y axis.
         */
        void (*scale)(const Image *inImg, Image *outImg, float scaleX, float scaleY);

        /**
         * @brief Rotates the input image by a specified angle.
         *
         * @param inImg Input image.
         * @param outImg Output image (rotated).
         * @param angle Rotation angle in degrees.
         */
        void (*rotate)(const Image *inImg, Image *outImg, float angle);

        /**
         * @brief Shears the input image by specified factors along the X and Y axes.
         *
         * @param inImg Input image.
         * @param outImg Output image (sheared).
         * @param shx Shear factor along the X axis.
         * @param shy Shear factor along the Y axis.
         */
        void (*shear)(const Image *inImg, Image *outImg, float shx, float shy);

        /**
         * @brief Inverse shears the input image by specified factors along the X and Y axes.
         *
         * @param inImg Input image.
         * @param outImg Output image (sheared).
         * @param shx Shear factor along the X axis.
         * @param shy Shear factor along the Y axis.
         */
        void (*shearInv)(const Image *inImg, Image *outImg, float shx, float shy);

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
         *
         * @sa Image
         */
        void (*warpPerspective)(const Image *inImg, Image *outImg, float M[3][3]);

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
        void (*perspectiveTransform)(const Image *inImg, Image *outImg,
                                       float h[3][3]);

    // Object Detection
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
        void (*templateMatching)(const Image *sourceImg, const Image *templateImg, Image *resultImg);


    //Miscellaneous

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
        void (*cvtColor)(Image *input, Image *output, uint8_t toGrayscale);

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
        void (*resize)(Image *inImg, Image *outImg, int fourier_size);

} DIP;

#endif /* DIP_INCLUDE_DIP_DRIVER_H_ */
