/*
 * DIP.c
 *
 *  Created on: Mar 28, 2024
 *      Author: ozand
 */
#include "DIP.h"
#include "DIP_UTIL.h"

uint32_t allocatedSize = 0;

/**
 * @brief Creates an image with a specified size and format.
 *
 * @param size Size of the image.
 * @param format Format of the image (e.g., grayscale, RGB).
 * @return Pointer to the created Image.
 */
Image* createImage(uint8_t size, uint8_t format) {

    // Declare a pointer to an Image structure and assign it to a memory address in SDRAM
    // SDRAM_BANK_ADDR and WRITE_READ_ADDR are predefined constants, likely pointing to the starting address of memory
    // allocatedSize keeps track of the current memory offset for storing new data
    Image *image = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;

    // Update the allocated size to account for the new Image structure
    allocatedSize += sizeof(Image);

    // Check if the image pointer is NULL, which would mean memory allocation failed
    if (image == NULL) {
        return NULL;  // Failed to allocate memory
    }

    // Determine the resolution (width and height) based on the provided size argument
    switch (size) {
        case IMAGE_RES_QQVGA: {
            image->width = IMAGE_RES_QQVGA_Width;  // Set QQVGA width (160 pixels)
            image->height = IMAGE_RES_QQVGA_Height;  // Set QQVGA height (120 pixels)
            break;
        }
        case IMAGE_RES_QVGA: {
            image->width = IMAGE_RES_QVGA_Width;  // Set QVGA width (320 pixels)
            image->height = IMAGE_RES_QVGA_Height;  // Set QVGA height (240 pixels)
            break;
        }
        case IMAGE_RES_480x272: {
            image->width = IMAGE_RES_480x272_Width;  // Set custom resolution width (480 pixels)
            image->height = IMAGE_RES_480x272_Height;  // Set custom resolution height (272 pixels)
            break;
        }
        case IMAGE_RES_VGA: {
            image->width = IMAGE_RES_VGA_Width;  // Set VGA width (640 pixels)
            image->height = IMAGE_RES_VGA_Height;  // Set VGA height (480 pixels)
            break;
        }
        default: {
            // Invalid resolution, return NULL to indicate an error
            return NULL;
        }
    }

    // Set the image size and format based on the format argument
    switch (format) {
        case IMAGE_FORMAT_RGB565: {
            image->size = image->width * image->height * 2;  // 2 bytes per pixel for RGB565 format
            image->format = IMAGE_FORMAT_RGB565;  // Assign the RGB565 format constant
            break;
        }
        case IMAGE_FORMAT_GRAYSCALE: {
            image->size = image->width * image->height;  // 1 byte per pixel for grayscale
            image->format = IMAGE_FORMAT_GRAYSCALE;  // Assign the Grayscale format constant
            break;
        }
        default: {
            // Invalid format, return NULL to indicate an error
            return NULL;
        }
    }

    // Assign the pointer to the pixel data, allocating memory for it after the Image structure
    image->pixels = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;

    // Update allocatedSize to account for the memory needed to store the image pixel data
    allocatedSize += image->size;

    // Return the pointer to the newly created Image structure
    return image;
}

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
void cvtColor(Image *input, Image *output, uint8_t toGrayscale) {

    // Check if the input image is valid (non-NULL and with valid pixel data)
    if (!input || !input->pixels) {
        return; // Invalid input image, return early
    }

    // Set output image dimensions to match the input image
    output->width = input->width;
    output->height = input->height;

    if (toGrayscale) {
        // Convert the image to Grayscale format

        // Set the output image format and size for Grayscale (1 byte per pixel)
        output->format = IMAGE_FORMAT_GRAYSCALE;
        output->size = output->width * output->height;
        // Allocate memory for the grayscale image in SDRAM
        output->pixels = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
        allocatedSize += output->size;

        // Check if memory allocation succeeded
        if (!output->pixels) {
            return; // Memory allocation failed
        }

        // If the input format is RGB565, perform conversion to Grayscale
        if (input->format == IMAGE_FORMAT_RGB565) {
            for (uint32_t i = 0; i < output->height * output->width; i++) {
                // Extract R, G, B components from RGB565
                uint16_t rgb565 = ((uint16_t*)input->pixels)[i];
                uint8_t r = (rgb565 >> 11) & 0x1F;  // Extract R (5 bits)
                uint8_t g = (rgb565 >> 5) & 0x3F;   // Extract G (6 bits)
                uint8_t b = rgb565 & 0x1F;          // Extract B (5 bits)

                // Scale 5-6-5 components to 8-bit values
                r = (r << 3) | (r >> 2);  // Scale R to 8 bits
                g = (g << 2) | (g >> 4);  // Scale G to 8 bits
                b = (b << 3) | (b >> 2);  // Scale B to 8 bits

                // Convert RGB to Grayscale using standard luminance formula
                output->pixels[i] = (uint8_t)(0.299f * r + 0.587f * g + 0.114f * b);
            }
        }
        // Handle conversion from YUV to Grayscale
        else if (input->format == IMAGE_FORMAT_YUV) {
            // Assuming YUV format with Y as the grayscale component
            for (uint32_t i = 0; i < output->height * output->width; i++) {
                uint8_t y = input->pixels[i]; // Use Y component for Grayscale
                output->pixels[i] = y;
            }
        }

    } else {
        // Convert from Grayscale to a different format (e.g., RGB565 or YUV)

        // Convert Grayscale to RGB565
        if (input->format == IMAGE_FORMAT_GRAYSCALE) {
            output->format = IMAGE_FORMAT_RGB565; // Set the output format to RGB565
            output->size = output->width * output->height * 2; // 2 bytes per pixel for RGB565
            // Allocate memory for the RGB565 image in SDRAM
            output->pixels = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
            allocatedSize += output->size;

            // Check if memory allocation succeeded
            if (!output->pixels) {
                return; // Memory allocation failed
            }

            // Convert Grayscale to RGB565 by copying the gray value into all RGB channels
            for (uint32_t i = 0; i < output->height * output->width; i++) {
                uint8_t gray = input->pixels[i];
                // Pack grayscale value into RGB565 format (R, G, B channels)
                uint16_t rgb565 = ((gray >> 3) << 11) | ((gray >> 2) << 5) | (gray >> 3);
                ((uint16_t*)output->pixels)[i] = rgb565; // Store as RGB565
            }
        }
        // Handle conversion from YUV to RGB565
        else if (input->format == IMAGE_FORMAT_YUV) {
            output->format = IMAGE_FORMAT_RGB565; // Set output format to RGB565
            output->size = output->width * output->height * 2; // 2 bytes per pixel for RGB565
            // Allocate memory for the RGB565 image in SDRAM
            output->pixels = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
            allocatedSize += output->size;

            // Check if memory allocation succeeded
            if (!output->pixels) {
                return; // Memory allocation failed
            }

            // Convert YUV to RGB565 by calculating the RGB values from YUV
            for (uint32_t i = 0; i < output->height * output->width; i++) {
                uint8_t y = input->pixels[i * 3];     // Y component
                uint8_t u = input->pixels[i * 3 + 1]; // U component
                uint8_t v = input->pixels[i * 3 + 2]; // V component

                // Convert YUV to RGB using the standard YUV to RGB conversion formula
                int r = y + (1.402f * (v - 128));
                int g = y - (0.344136f * (u - 128)) - (0.714136f * (v - 128));
                int b = y + (1.772f * (u - 128));

                // Clamp the RGB values to ensure they are within the valid range [0, 255]
                r = r < 0 ? 0 : (r > 255 ? 255 : r);
                g = g < 0 ? 0 : (g > 255 ? 255 : g);
                b = b < 0 ? 0 : (b > 255 ? 255 : b);

                // Pack the RGB values into RGB565 format (R, G, B channels)
                uint16_t rgb565 = ((r >> 3) << 11) | ((g >> 2) << 5) | (b >> 3);
                ((uint16_t*)output->pixels)[i] = rgb565; // Store as RGB565
            }
        }
    }
}


/**
 * @brief Allocates space for the Fourier transform data of the input image.
 *
 * This function checks the format of the input image and allocates memory
 * in a designated SDRAM area for storing the Fourier transform results.
 * If the image format is RGB565, the function exits early without
 * performing any operations.
 *
 * @param inImg Pointer to the Image structure containing the image data
 *               and relevant properties such as format and size.
 */
void addFourier(Image *inImg) {
    // Check if the input image format is RGB565
    if (inImg->format == IMAGE_FORMAT_RGB565) {
        // Exit early if the format is not suitable for Fourier transformation
        return; // Note: returning NULL is incorrect here as the return type is void
    }

    // Allocate memory for Fourier transform results
    // Note: Currently, the size allocated may exceed the required amount;
    // it should be edited to fit the actual requirements.

    // Allocate space for the Fourier transform in the x-direction
    inImg->fourierx = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
    allocatedSize += inImg->size * 8; // Update the allocated size for the x-direction

    // Allocate space for the Fourier transform in the y-direction
    inImg->fouriery = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
    allocatedSize += inImg->size * 8; // Update the allocated size for the y-direction

    // No return value is needed as the function's return type is void
    return;
}

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
void DIP_Resize(Image *inImg, Image *outImg, int fourier_size) {

	float width_ratio = (float) inImg->width / fourier_size;
	float height_ratio = (float) inImg->height / fourier_size;

	for (int y = 0; y < outImg->size; y++)
		outImg->pixels[y] = 0xff;

	for (int y = 0; y < fourier_size; y++) {
		for (int x = 0; x < fourier_size; x++) {
			int nearest_x = (int) (x * width_ratio);
			int nearest_y = (int) (y * height_ratio);

			// Get the pixel value from the nearest neighbor
			outImg->pixels[y * fourier_size + x] = inImg->pixels[nearest_y
					* inImg->width + nearest_x];
		}
	}
}

/**
 * @brief Negate the pixel intensities of an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with negated intensities.
 */
void DIP_negative(const Image *inImg, Image *outImg) {

	uint8_t *imgData = inImg->pixels;
	uint8_t *outData = outImg->pixels;

	for (int i = 0; i < inImg->size; ++i) {
		outData[i] = MAX_INTENSITY - imgData[i];
	}
}

/**
 * @brief Apply power law transformation to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with power-law transformed intensities.
 * @param[in] gamma The power-law transformation parameter.
 */
void DIP_powerTransform(const Image *inImg, Image *outImg, uint8_t gamma) {

	uint8_t *imgData = inImg->pixels;
	uint8_t *outData = outImg->pixels;

	// Calculate the scaling factor 'c' to keep pixel values in the valid intensity range
	float c = MAX_INTENSITY / pow(MAX_INTENSITY, gamma);

	// Perform the power transform
	for (int i = 0; i < inImg->size; ++i) {

		outData[i] = (uint8_t) (c * pow(imgData[i], gamma));
	}
}

/**
 * @brief Calculate the histogram of an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] histogram Array to store the histogram.
 * @return Total number of pixels in the image.
 */
int DIP_histogramFormation(const Image *inImg, int *histogram) {

	uint8_t *imgData = inImg->pixels;

	// Calculate histogram
	for (int i = 0; i < inImg->size; ++i) {
		histogram[imgData[i]]++;
	}
	// Return the total number of pixels in the image
	return inImg->size;
}

/**
 * @brief Perform histogram equalization on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with equalized intensities.
 */
void DIP_histogramEqualization(const Image *inImg, Image *outImg) {

	uint8_t *imgData = inImg->pixels;
	uint8_t *outData = outImg->pixels;

	// Calculate histogram
	int histogram[MAX_INTENSITY + 1] = { 0 };
	for (int i = 0; i < inImg->size; ++i) {
		histogram[imgData[i]]++;
	}

	// Calculate cumulative distribution function
	int cumulativeHistogram[MAX_INTENSITY + 1] = { 0 };
	cumulativeHistogram[0] = histogram[0];
	for (int i = 1; i <= MAX_INTENSITY; ++i) {
		cumulativeHistogram[i] = cumulativeHistogram[i - 1] + histogram[i];
	}

	// Calculate equalization mapping
	float scale = (float) MAX_INTENSITY / inImg->size;
	uint8_t equalizationMap[MAX_INTENSITY + 1];
	for (int i = 0; i <= MAX_INTENSITY; ++i) {
		equalizationMap[i] = (uint8_t) (round(cumulativeHistogram[i] * scale));
	}

	// Perform histogram equalization on the image
	for (int i = 0; i < inImg->size; ++i) {
		outData[i] = equalizationMap[imgData[i]];
	}
}

/**
 * @brief Perform histogram specification on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with specified histogram.
 * @param[in] targetHistogram The target histogram to match.
 */
void DIP_histogramSpecification(const Image *inImg, Image *outImg,
		const int *targetHistogram) {

	uint8_t *imgData = inImg->pixels;
	uint8_t *outData = outImg->pixels;

	int sourceHistogram[MAX_INTENSITY + 1] = { 0 };
	for (int i = 0; i < inImg->size; ++i) {
		sourceHistogram[imgData[i]]++;
	}

	// Calculate cumulative distribution function for source and target histograms
	int sourceCumulativeHistogram[MAX_INTENSITY + 1] = { 0 };
	int targetCumulativeHistogram[MAX_INTENSITY + 1] = { 0 };

	sourceCumulativeHistogram[0] = sourceHistogram[0];
	targetCumulativeHistogram[0] = targetHistogram[0];

	for (int i = 1; i <= MAX_INTENSITY; ++i) {
		sourceCumulativeHistogram[i] = sourceCumulativeHistogram[i - 1]
				+ sourceHistogram[i];
		targetCumulativeHistogram[i] = targetCumulativeHistogram[i - 1]
				+ targetHistogram[i];
	}

	// Calculate equalization mapping
	float scale = (float) MAX_INTENSITY / inImg->size;
	uint8_t equalizationMap[MAX_INTENSITY + 1];
	for (int i = 0; i <= MAX_INTENSITY; ++i) {
		int nearestMatch = 0;
		int minDiff = abs(
				sourceCumulativeHistogram[i] - targetCumulativeHistogram[0]);

		for (int j = 1; j <= MAX_INTENSITY; ++j) {
			int diff = abs(
					sourceCumulativeHistogram[i]
							- targetCumulativeHistogram[j]);
			if (diff < minDiff) {
				minDiff = diff;
				nearestMatch = j;
			}
		}

		equalizationMap[i] = nearestMatch;
	}

	// Perform histogram specification on the image
	for (int i = 0; i < inImg->size; ++i) {
		outData[i] = equalizationMap[imgData[i]];
	}
}

/**
 * @brief Compare two histograms.
 *
 * @param[in] hist1 First histogram (array of size MAX_INTENSITY + 1).
 * @param[in] hist2 Second histogram (array of size MAX_INTENSITY + 1).
 * @param[in] method Comparison method (e.g., HIST_COMP_CORRELATION, HIST_COMP_CHI_SQUARE).
 * @return Comparison result as a double.
 */
double DIP_compareHist(int *hist1, int *hist2, int method) {
    double result = 0.0;

    // Normalize histograms
    int total1 = 0, total2 = 0;
    for (int i = 0; i <= MAX_INTENSITY; ++i) {
        total1 += hist1[i];
        total2 += hist2[i];
    }

    double normHist1[MAX_INTENSITY + 1], normHist2[MAX_INTENSITY + 1];
    for (int i = 0; i <= MAX_INTENSITY; ++i) {
        normHist1[i] = (total1 > 0) ? (double)hist1[i] / total1 : 0;
        normHist2[i] = (total2 > 0) ? (double)hist2[i] / total2 : 0;
    }

    switch (method) {
        case HIST_COMP_CORRELATION: {
            double sum1 = 0.0, sum2 = 0.0, sum1Sq = 0.0, sum2Sq = 0.0, pSum = 0.0;
            for (int i = 0; i <= MAX_INTENSITY; ++i) {
                sum1 += normHist1[i];
                sum2 += normHist2[i];
                sum1Sq += normHist1[i] * normHist1[i];
                sum2Sq += normHist2[i] * normHist2[i];
                pSum += normHist1[i] * normHist2[i];
            }

            double num = pSum - (sum1 * sum2);
            double den = sqrt((sum1Sq - (sum1 * sum1)) * (sum2Sq - (sum2 * sum2)));

            result = (den == 0) ? 0 : num / den;
            break;
        }
        case HIST_COMP_CHI_SQUARE: {
            for (int i = 0; i <= MAX_INTENSITY; ++i) {
                if (normHist1[i] + normHist2[i] > 0) {
                    result += ((normHist1[i] - normHist2[i]) * (normHist1[i] - normHist2[i])) / (normHist1[i] + normHist2[i]);
                }
            }
            break;
        }
        case HIST_COMP_INTERSECTION: {
            for (int i = 0; i <= MAX_INTENSITY; ++i) {
                result += fmin(normHist1[i], normHist2[i]);
            }
            break;
        }
        case HIST_COMP_BHATTACHARYYA: {
            double pSum = 0.0;
            for (int i = 0; i <= MAX_INTENSITY; ++i) {
                pSum += sqrt(normHist1[i] * normHist2[i]);
            }
            result = sqrt(1.0 - pSum);  // Proper calculation using normalized histograms
            break;
        }
        default:
            result = -1.0; // Invalid method
            break;
    }

    return result;
}


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
void DIP_blendLinear(const Image *baseImg, const Image *overlayImg, Image *outImg, float alpha) {
    // Ensure alpha is in the range [0, 1]
    if (alpha < 0.0f) alpha = 0.0f;
    if (alpha > 1.0f) alpha = 1.0f;

    uint8_t *baseData = baseImg->pixels;
    uint8_t *overlayData = overlayImg->pixels;
    uint8_t *outData = outImg->pixels;

    // Check that the base image and overlay image have the same dimensions
    if (baseImg->width != overlayImg->width || baseImg->height != overlayImg->height) {
        printf("Base image and overlay image must have the same dimensions.\n");
        return;
    }

    // Blend the two images
    for (size_t i = 0; i < baseImg->size; ++i) {
        // Blend each pixel
        outData[i] = (uint8_t)((1.0f - alpha) * baseData[i] + alpha * overlayData[i]);
    }
}


/**
 * @brief Embed a grayscale watermark into an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with embedded watermark.
 * @param[in] watermark The watermark to embed.
 * @param[in] x X-coordinate of the watermark insertion point.
 * @param[in] y Y-coordinate of the watermark insertion point.
 */
void DIP_grayscaleWatermark(const Image *inImg, Image *outImg,
		const char *watermark, int x, int y) {

	uint8_t *imgData = inImg->pixels;
	uint8_t *outData = outImg->pixels;

	// Make sure the watermark fits within the image dimensions
	int textWidth = strlen(watermark);
	if (x < 0 || x + textWidth >= inImg->width || y < 0 || y >= inImg->height) {
		printf("Watermark does not fit within the image boundaries.\n");
		return;
	}

	// Insert the watermark into the image
	for (size_t i = 0; i < inImg->size; ++i) {

		if (i >= y * inImg->width + x
				&& i <= y * inImg->width + x + textWidth) {
			outData[i] = watermark[y];
		} else {
			outData[i] = imgData[i];
		}

	}

}

/**
 * @brief Applies a 2D filter to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (filtered).
 * @param size Size of the filter matrix (e.g., 3x3, 5x5).
 * @param filter 2D array representing the filter.
 */
void DIP_filter2D(const Image *inImg, Image *outImg, int size,
		float filter[size][size]) {

	int filtersize = size;

	for (int y = 0; y < inImg->height; y++) {
		for (int x = 0; x < inImg->width; x++) {
			int sum = 0;
			for (int fy = -filtersize / 2; fy <= filtersize / 2; ++fy) {
				for (int fx = -filtersize / 2; fx <= filtersize / 2; ++fx) {
					int pixelValue = 0;
					if (y + fy > inImg->height - 1 || x + fx > inImg->width - 1
							|| y + fy < 0 || x + fx < 0) {
						pixelValue = 0;
					} else {
						pixelValue = inImg->pixels[(y + fy) * inImg->width
								+ (x + fx)];
					}
					sum += ((float) pixelValue * filter[fy][fx]);
				}
			}

			sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum);
			outImg->pixels[y * inImg->width + x] = (uint8_t) sum; // Normalize by dividing by 9 for a 3x3 filter
		}
	}

}

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
void DIP_sepFilter2D(const Image *inImg, Image *outImg, int kernelSizeX, float kernelX[kernelSizeX],
                     int kernelSizeY, float kernelY[kernelSizeY], double delta) {

    // Temporary image to store the result after the first pass (horizontal filtering)
    Image *tempImg = (Image*) createImage(IMAGE_RES_480x272,
    		IMAGE_FORMAT_GRAYSCALE);

    // First pass: filter rows with kernelX (horizontal filtering)
    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
            float sum = 0.0;
            for (int fx = -kernelSizeX / 2; fx <= kernelSizeX / 2; ++fx) {
                int pixelValue = 0;
                int xIndex = x + fx;

                // Handle borders by ignoring pixels outside bounds (border handling)
                if (xIndex < 0 || xIndex >= inImg->width) {
                    pixelValue = 0;
                } else {
                    pixelValue = inImg->pixels[y * inImg->width + xIndex];
                }

                sum += ((float) pixelValue * kernelX[fx + kernelSizeX / 2]);
            }

            sum = sum + delta; // Add delta after the filter operation
            sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum); // Normalize to [0, 255]
            tempImg->pixels[y * tempImg->width + x] = (uint8_t)sum;
        }
    }

    // Second pass: filter columns with kernelY (vertical filtering)
    for (int y = 0; y < outImg->height; y++) {
        for (int x = 0; x < outImg->width; x++) {
            float sum = 0.0;
            for (int fy = -kernelSizeY / 2; fy <= kernelSizeY / 2; ++fy) {
                int pixelValue = 0;
                int yIndex = y + fy;

                // Handle borders by ignoring pixels outside bounds (border handling)
                if (yIndex < 0 || yIndex >= tempImg->height) {
                    pixelValue = 0;
                } else {
                    pixelValue = tempImg->pixels[yIndex * tempImg->width + x];
                }

                sum += ((float) pixelValue * kernelY[fy + kernelSizeY / 2]);
            }

            sum = sum + delta; // Add delta after the filter operation
            sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum); // Normalize to [0, 255]
            outImg->pixels[y * outImg->width + x] = (uint8_t)sum;
        }
    }

}

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
void DIP_sqrBoxFilter(const Image *inImg, Image *outImg, int filtersize, int normalize) {
    int halfSize = filtersize / 2;
    int maxPixelValue = 255;
    double maxSumSq = (filtersize * filtersize) * (maxPixelValue * maxPixelValue);

    for (int y = 0; y < inImg->height; ++y) {
        for (int x = 0; x < inImg->width; ++x) {
            double sumSq = 0.0;

            // Apply the square box filter
            for (int fy = -halfSize; fy <= halfSize; ++fy) {
                for (int fx = -halfSize; fx <= halfSize; ++fx) {
                    int ny = y + fy;
                    int nx = x + fx;

                    // Handle boundary conditions by skipping out-of-bounds indices
                    if (ny >= 0 && ny < inImg->height && nx >= 0 && nx < inImg->width) {
                        uint8_t pixel = inImg->pixels[ny * inImg->width + nx];
                        sumSq += pixel * pixel;  // Sum of squares
                    }
                }
            }

            // Normalize based on the maximum possible sum of squares
            if (normalize) {
                sumSq /= maxSumSq;
                sumSq *= 255;  // Scale to the range [0, 255]
            }

            // Clamp the result to the valid 8-bit range [0, 255]
            sumSq = sumSq < 0 ? 0 : (sumSq > 255 ? 255 : sumSq);

            // Assign the result to the output image
            outImg->pixels[y * inImg->width + x] = (uint8_t)sumSq;
        }
    }
}

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
void DIP_Laplacian(const Image *inImg, Image *outImg, int ddepth, int ksize, double scale, double delta, int borderType) {
    // Check if ksize is valid (positive and odd)
    if (ksize <= 0 || ksize % 2 == 0) {
        return; // Invalid ksize, exit the function
    }

    // Create the Laplacian filter
    float **filter = createLaplacianFilter(ksize);
    if (!filter) {
        return; // Handle memory allocation failure
    }

    int halfSize = ksize / 2;

    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
            double sum = 0.0;

            for (int fy = -halfSize; fy <= halfSize; ++fy) {
                for (int fx = -halfSize; fx <= halfSize; ++fx) {
                    int pixelValue = 0;
                    int newY = y + fy;
                    int newX = x + fx;

                    // Handle borders
                    if (newY < 0 || newY >= inImg->height || newX < 0 || newX >= inImg->width) {
                        if (borderType == BORDER_DEFAULT) {
                            pixelValue = 0; // Set to 0 for border handling
                        }
                        // You can add more border handling options here
                    } else {
                        pixelValue = inImg->pixels[newY * inImg->width + newX];
                    }

                    sum += pixelValue * filter[fy + halfSize][fx + halfSize];
                }
            }

            // Apply scale and delta, and clamp the values to [0, 255]
            sum = sum * scale + delta;
            sum = fmin(fmax(sum, 0), 255); // Clamp to [0, 255]
            outImg->pixels[y * inImg->width + x] = (uint8_t)sum;
        }
    }

    for (int i = 0; i < ksize; i++) {
        free(filter[i]);
    }
    free(filter);
}

/**
 * @brief Applies a Gaussian blur to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (blurred).
 * @param filtersize Size of the Gaussian filter.
 * @param sigma Standard deviation of the Gaussian distribution.
 */
void DIP_gaussianBlur(const Image *inImg, Image *outImg, int filtersize,
		int sigma) {

	FilterKernel kernel;
	kernel.size = filtersize;
	kernel.data = (float**) malloc(filtersize * sizeof(float*));
	double sum = 0.0;
	for (int i = 0; i < filtersize; ++i) {
		kernel.data[i] = (float*) malloc(filtersize * sizeof(float));
		for (int j = 0; j < filtersize; ++j) {
			int x = i - filtersize / 2;
			int y = j - filtersize / 2;
			kernel.data[i][j] = exp(-(x * x + y * y) / (2 * sigma * sigma))
					/ (2 * M_PI * sigma * sigma);
			sum += kernel.data[i][j];
		}
	}
	// Normalize the kernel
	for (int i = 0; i < filtersize; ++i) {
		for (int j = 0; j < filtersize; ++j) {
			kernel.data[i][j] /= sum;
		}
	}

	for (int y = 0; y < inImg->height; y++) {
		for (int x = 0; x < inImg->width; x++) {
			int sum = 0;
			for (int fy = -filtersize / 2; fy <= filtersize / 2; ++fy) {
				for (int fx = -filtersize / 2; fx <= filtersize / 2; ++fx) {
					int pixelValue = 0;
					if (y + fy > inImg->height - 1 || x + fx > inImg->width - 1
							|| y + fy < 0 || x + fx < 0) {
						pixelValue = 0;
					} else {
						pixelValue = inImg->pixels[(y + fy) * inImg->width
								+ (x + fx)];
					}
					sum += ((float) pixelValue
							* kernel.data[fy + filtersize / 2][fx
									+ filtersize / 2]);
				}
			}

			sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum);
			outImg->pixels[y * inImg->width + x] = (uint8_t) sum;

		}
	}
}

/**
 * @brief Applies a median blur to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (blurred).
 * @param filtersize Size of the filter.
 */
void DIP_medianBlur(const Image *inImg, Image *outImg, int filtersize) {

	FilterKernel kernel;
	kernel.size = filtersize;
	kernel.data_8 = (uint8_t*) malloc(
			filtersize * filtersize * sizeof(uint8_t*));

	for (int y = 0; y < inImg->height; y++) {

		for (int x = 0; x < inImg->width; x++) {

			//get filter items
			for (int fy = -filtersize / 2; fy <= filtersize / 2; ++fy) {
				for (int fx = -filtersize / 2; fx <= filtersize / 2; ++fx) {
					if (y + fy > inImg->height - 1 || x + fx > inImg->width - 1
							|| y + fy < 0 || x + fx < 0) {
						kernel.data_8[(fy + filtersize / 2) * filtersize + fx
								+ filtersize / 2] = 0;
					} else {
						kernel.data_8[(fy + filtersize / 2) * filtersize + fx
								+ filtersize / 2] = inImg->pixels[(y + fy)
								* inImg->width + (x + fx)];
					}
				}
			}

			// Apply median filtering directly using kernel_data pointer
			qsort(kernel.data_8, filtersize * filtersize, sizeof(uint8_t),
					(const void*) compare_uint8);

			// Assign the median value to the output image
			outImg->pixels[y * inImg->width + x] =
					(uint8_t) kernel.data_8[(filtersize / 2) * filtersize
							+ filtersize / 2];

		}
	}
}

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
void DIP_pyrDown(const Image *src, Image *dst, const Size *dstsize, int borderType) {
    // Determine output size
    int outputWidth = (dstsize != NULL) ? dstsize->width : (src->width + 1) / 2;
    int outputHeight = (dstsize != NULL) ? dstsize->height : (src->height + 1) / 2;

    // Apply Gaussian blur before downsampling
    Image *blurredImg = (Image*) createImage(IMAGE_RES_480x272,
    		IMAGE_FORMAT_GRAYSCALE);
    // Using your Gaussian blur function
    int filterSize = 5; // Example filter size, can be adjusted
    int sigma = 1;      // Example sigma value, can be adjusted
    DIP_gaussianBlur(src, blurredImg, filterSize, sigma);

    for(int i=0;i<dst->size;i++){
    	dst->pixels[i] = 0x00;
    }
    // Downsample the image by selecting every second pixel
    for (int y = 0; y < outputHeight; y++) {
        for (int x = 0; x < outputWidth; x++) {
            dst->pixels[y * src->width + x] = blurredImg->pixels[2 * y * src->width + 2 * x];
        }
    }

}

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
void DIP_pyrUp(const Image *inImg, Image *outImg, Size *dstSize, int borderType, int filtersize, int sigma) {
    // Default output size calculation if dstSize is not set
    if (dstSize->width == 0 || dstSize->height == 0) {
        dstSize->width = inImg->width * 2;
        dstSize->height = inImg->height * 2;
    }

    //check
    Image *processed = (Image*) createImage(IMAGE_RES_480x272,
    		IMAGE_FORMAT_GRAYSCALE);

    // Initialize output image to zero
    memset(processed->pixels, 0, processed->width * processed->height * sizeof(uint8_t));

    // Upsample by inserting zeros
    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
        	processed->pixels[(y * 2) * processed->width + (x * 2)] = inImg->pixels[y * inImg->width + x];
        }
    }

    // Apply Gaussian blur to the upsampled image using your existing function
    DIP_gaussianBlur(processed, outImg, filtersize, sigma);
}

/**
 * @brief Computes the Fourier transform of the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (Fourier transformed).
 */
void DIP_fourier(const Image *inImg, Image *outImg) {

	int imageN = 256;

	addFourier(outImg);

	float *fourier = outImg->fourierx;
	float *fourier2 = outImg->fouriery;

	for (int row = 0; row < imageN * imageN; row++) {
		fourier[2 * row] = (uint32_t) inImg->pixels[row];
		fourier[2 * row + 1] = 0x00000000;
	}

	for (int i = 0; i < imageN; i++) {
		arm_cfft_f32(&arm_cfft_sR_f32_len256, fourier + imageN * i * 2, 0, 1);
	}

	for (int k = 0; k < imageN; k++) {
		for (int j = 0; j < imageN; j++) {
			fourier2[2 * j + k * imageN * 2] = (float) fourier[j * imageN * 2
					+ k * 2];
			fourier2[2 * j + 1 + k * imageN * 2] = (float) fourier[j * imageN
					* 2 + k * 2 + 1];
		}
	}

	for (int i = 0; i < imageN; i++) {
		arm_cfft_f32(&arm_cfft_sR_f32_len256, fourier2 + imageN * i * 2, 0, 1);
	}

}

/**
 * @brief Gets the magnitude spectrum from a Fourier transformed image.
 *
 * @param inImg Input image (Fourier domain).
 * @param outImg Output image (magnitude spectrum).
 * @param fftshift 1 to apply fftshift, 0 otherwise.
 */
void DIP_abs(Image *outImg) {

	int imageN = 256;

	uint8_t *tesst = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize = allocatedSize + imageN * imageN * 4;
	int *test_arr = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize = allocatedSize + imageN * imageN * 4;

	float *fourier = outImg->fourierx;
	float *fourier2 = outImg->fouriery;

	int max = 0;

	for (int i = 0; i < imageN * imageN; i++) {
		test_arr[i] = (int) log(
				sqrt(
						fourier2[2 * i] * fourier2[2 * i]
								+ fourier2[2 * i + 1] * fourier2[2 * i + 1])
						+ 1);

		if (test_arr[i] > max)
			max = test_arr[i];
	}
	int temps = 0;
	for (int i = 0; i < imageN * imageN; i++) {
		temps = 255 * ((float) test_arr[i] / max);
		tesst[i] = temps;
	}

	int halfWidth = imageN / 2;
	int halfHeight = imageN / 2;
	int width = imageN;

	for (int y = 0; y < halfHeight; y++) {
		for (int x = 0; x < halfWidth; x++) {
			// Swap quadrants
			int quad1Index = y * width + x;
			int quad2Index = quad1Index + halfWidth;
			int quad3Index = (y + halfHeight) * width + x;
			int quad4Index = quad3Index + halfWidth;

			uint32_t temp = tesst[quad1Index];
			tesst[quad1Index] = tesst[quad4Index];
			tesst[quad4Index] = temp;

			temp = tesst[quad2Index];
			tesst[quad2Index] = tesst[quad3Index];
			tesst[quad3Index] = temp;
		}
	}

	for (int i = 0; i < outImg->size; i++)
		outImg->pixels[i] = 0xff;

	for (int i = 0; i < imageN * imageN; i++) {
		int row = i / imageN;
		int col = i % imageN;
		outImg->pixels[row * 480 + col] = tesst[i];
	}
}

/**
 * @brief Computes the inverse Fourier transform of the input image.
 *
 * @param inImg Input image (Fourier domain).
 * @param outImg Output image (spatial domain).
 */
void DIP_fourierInv(const Image *inImg, Image *outImg) {

	int imageN = 256;
	float *fourier = inImg->fourierx;
	float *fourier2 = inImg->fouriery;

	if (inImg->fourierx == NULL && inImg->fouriery == NULL) {
		return -1;
	}

	for (int i = 0; i < imageN; i++) {
		arm_cfft_f32(&arm_cfft_sR_f32_len256, fourier2 + imageN * i * 2, 1, 1);
	}

	for (int k = 0; k < imageN; k++) {
		for (int j = 0; j < imageN; j++) {
			fourier[2 * j + k * imageN * 2] = (float) fourier2[j * imageN * 2
					+ k * 2];
			fourier[2 * j + 1 + k * imageN * 2] = (float) fourier2[j * imageN
					* 2 + k * 2 + 1];
		}
	}

	for (int i = 0; i < imageN; i++) {
		arm_cfft_f32(&arm_cfft_sR_f32_len256, fourier + imageN * i * 2, 1, 1);
	}

	uint8_t *tesst = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize = allocatedSize + imageN * imageN * 4;
	int *test_arr = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize = allocatedSize + imageN * imageN * 4;

	int max;
	for (int i = 0; i < imageN * imageN; i++) {

		test_arr[i] = (int) sqrt(
				fourier[2 * i] * fourier[2 * i]
						+ fourier[2 * i + 1] * fourier[2 * i + 1]);

		if (test_arr[i] > max)
			max = test_arr[i];
	}

	for (int i = 0; i < outImg->size; i++)
		outImg->pixels[i] = 0xff;

	for (int i = 0; i < imageN * imageN; i++) {
		int row = i / imageN;
		int col = i % imageN;
		outImg->pixels[row * 480 + col] = test_arr[i];
	}

}

/**
 * @brief Apply grayscale thresholding to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output binary image after thresholding.
 * @param[in] threshold Threshold value.
 */
void DIP_grayscaleThreshold(const Image *inImg, Image *outImg, int threshold) {
	for (int i = 0; i < inImg->size; ++i) {
		outImg->pixels[i] = (inImg->pixels[i] > threshold) ? 0xFF : 0x00;
	}
}

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
void DIP_adaptiveThreshold(const Image *inImg, Image *outImg, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C) {
    int halfBlockSize = blockSize / 2;

    for (int y = 0; y < inImg->height; ++y) {
        for (int x = 0; x < inImg->width; ++x) {
            int sum = 0;
            int count = 0;

            // Calculate the mean of the neighborhood pixels
            for (int j = -halfBlockSize; j <= halfBlockSize; ++j) {
                for (int i = -halfBlockSize; i <= halfBlockSize; ++i) {
                    int neighborX = x + i;
                    int neighborY = y + j;

                    // Check if the neighbor coordinates are within bounds
                    if (neighborX >= 0 && neighborX < inImg->width &&
                        neighborY >= 0 && neighborY < inImg->height) {
                        sum += inImg->pixels[neighborY * inImg->width + neighborX];
                        count++;
                    }
                }
            }

            // Calculate the mean and apply the threshold
            int mean = (count > 0) ? (sum / count) : 0; // Avoid division by zero
            int thresholdValue = (int)(mean - C);

            // Apply the thresholding based on the specified type
            if (thresholdType == THRESH_BINARY) {
                outImg->pixels[y * inImg->width + x] = (inImg->pixels[y * inImg->width + x] > thresholdValue) ? maxValue : 0;
            } else if (thresholdType == THRESH_BINARY_INV) {
                outImg->pixels[y * inImg->width + x] = (inImg->pixels[y * inImg->width + x] > thresholdValue) ? 0 : maxValue;
            }
        }
    }
}

/**
 * @brief Applies Otsu's method to perform automatic thresholding.
 *
 * @param inImg Input image.
 * @param outImg Output image (binary result).
 */
void DIP_otsuMethod(const Image *inImg, Image *outImg) {

	int histogram[256] = { 0 };

	DIP_histogramFormation(inImg, histogram);

	double prob[256], omega[256]; /* prob of graylevels */
	double myu[256]; /* mean value for separation */
	double max_sigma, sigma[256]; /* inter-class variance */
	int i; /* Loop variable */
	int threshold; /* threshold for binarization */

	/* calculation of probability density */
	for (i = 0; i < 256; i++) {
		prob[i] = (double) histogram[i] / (inImg->size);
	}

	/* omega & myu generation */
	omega[0] = prob[0];
	myu[0] = 0.0; /* 0.0 times prob[0] equals zero */
	for (i = 1; i < 256; i++) {
		omega[i] = omega[i - 1] + prob[i];
		myu[i] = myu[i - 1] + i * prob[i];
	}

	/* sigma maximization
	 sigma stands for inter-class variance
	 and determines optimal threshold value */
	threshold = 0;
	max_sigma = 0.0;
	for (i = 0; i < 256 - 1; i++) {
		if (omega[i] != 0.0 && omega[i] != 1.0)
			sigma[i] = pow(myu[256 - 1] * omega[i] - myu[i], 2)
					/ (omega[i] * (1.0 - omega[i]));
		else
			sigma[i] = 0.0;
		if (sigma[i] > max_sigma) {
			max_sigma = sigma[i];
			threshold = i;
		}
	}

	DIP_grayscaleThreshold(inImg, outImg, threshold);
}

/**
 * @brief Apply Sobel edge detection to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output edge-detected image.
 */
// Function to apply the Sobel filter
void DIP_sobelFilter(const Image *inImg, Image *outImg) {
	int sobelX[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
	int sobelY[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };

	// Iterate over each pixel of the input image
	for (int y = 0; y < inImg->height; ++y) {
		for (int x = 0; x < inImg->width; ++x) {
			int gradientX = 0;
			int gradientY = 0;
			int pixelValue = 0;
			// Compute gradient in x and y directions

			//BORDERLAR OPENCV ILE BIREBIR GITMEK ADINA DUZELTILDI.

			for (int fy = 0; fy < 3; ++fy) {
				for (int fx = 0; fx < 3; ++fx) {

					if (x == inImg->width - 1 && fx == 2) { // right border
						pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
								+ (x + fx - 1 - 1)];
						;
					} else if (y + fy - 1 < 0 || x + fx - 1 < 0) { // top
						if (y + fy - 1 < 0 || x + fx - 1 < 0) {
							pixelValue = inImg->pixels[(y + fy) * inImg->width
									+ (x + fx)];
						} else if (y + fy - 1 < 0) {
							pixelValue = inImg->pixels[(y + fy) * inImg->width
									+ (x + fx - 1)];
						} else {
							pixelValue = inImg->pixels[(y + fy - 1)
									* inImg->width + (x + fx)];
						}

					} else if ((y + fy - 1) * inImg->width + (x + fx - 1)
							> inImg->size) { // bottom
						pixelValue = inImg->pixels[(y + fy - 1 - 1)
								* inImg->width + (x + fx - 1)];
					} else if (x == 0 && fx == 0) { // left border
						pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
								+ (x + fx)];
					} else {
						pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
								+ (x + fx - 1)];
					}

					gradientX += pixelValue * sobelX[fy][fx];
					gradientY += pixelValue * sobelY[fy][fx];
				}
			}

			// Calculate gradient magnitude
			double magnitude = sqrt(
					(double) (gradientX * gradientX + gradientY * gradientY));

			// Normalize and clamp the magnitude
			outImg->pixels[y * inImg->width + x] = (uint8_t) fmin(
					fmax(magnitude, 0), 255);
		}
	}
}

/**
 * @brief Performs region growing segmentation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (segmented).
 * @param y Starting Y coordinate for region growing.
 * @param x Starting X coordinate for region growing.
 * @param threshold Intensity threshold for region growing.
 */
void DIP_regionGrowing(Image *inImg, Image *outImg, int y, int x, int threshold) {

	Image *mask = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);

	for(int i=0;i<inImg->size;i++)mask->pixels[i]=0;

    // Define connectivity (8-connected neighborhood)
    int connectivity[8][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

    // Initialize the region with the seed point
	int *queue = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize = allocatedSize + inImg->height * inImg->width * 2;

    int queueFront = 0, queueRear = 0;
    queue[queueRear++] = x;
    queue[queueRear++] = y;

    // Start region growing
    while (queueFront < queueRear) {
        // Get the current pixel from the queue
        int cx = queue[queueFront++];
        int cy = queue[queueFront++];

        // Check 8-connected neighbors
        for (int i = 0; i < 8; i++) {
            int nx = cx + connectivity[i][0];
            int ny = cy + connectivity[i][1];

            // Check if the neighbor is within the image bounds
            if (nx >= 0 && nx < inImg->height && ny >= 0 && ny < inImg->width) {
                // Check if the neighbor pixel has not been visited and intensity within threshold
                if (mask->pixels[nx * inImg->width + ny]  != 1 && abs(inImg->pixels[nx * inImg->width + ny] - inImg->pixels[x * inImg->width + y]) <= threshold) {
                    // Add the neighbor pixel to the queue
                    queue[queueRear++] = nx;
                    queue[queueRear++] = ny;

                    // Mark the neighbor pixel as visited
                    mask->pixels[nx * inImg->width + ny] = 1;
                }
            }
        }
    }


	for(int i=0;i<inImg->size;i++)outImg->pixels[i] = inImg->pixels[i] * mask->pixels[i];
}

/**
 * @brief Performs basic segmentation on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (segmented).
 * @param targetValue Value used for segmentation.
 */
void DIP_basicSegmentation(const Image *inImg, Image *outImg,
		uint8_t targetValue) {
	for (int i = 0; i < inImg->size; ++i) {
		outImg->pixels[i] = (inImg->pixels[i] == targetValue) ? 255 : 0;
	}
}

/**
 * @brief Performs image erosion using a specified kernel size.
 *
 * @param inImg Input image.
 * @param outImg Output image (eroded).
 * @param kernelSize Size of the erosion kernel.
 */
void DIP_erode(const Image *inImg, Image *outImg, int kernelSize) {
	int kernelRadius = kernelSize / 2;

	for (int y = 0; y < inImg->height; ++y) {
		for (int x = 0; x < inImg->width; ++x) {
			uint8_t minPixelValue = MAX_INTENSITY;

			for (int ky = -kernelRadius; ky <= kernelRadius; ++ky) {
				for (int kx = -kernelRadius; kx <= kernelRadius; ++kx) {
					int offsetX = x + kx;
					int offsetY = y + ky;

					if (offsetX >= 0 && offsetX < inImg->width && offsetY >= 0
							&& offsetY < inImg->height) {
						int pixelValue = inImg->pixels[offsetY * inImg->width
								+ offsetX];
						if (pixelValue < minPixelValue) {
							minPixelValue = pixelValue;
						}
					}
				}
			}

			outImg->pixels[y * inImg->width + x] = minPixelValue;
		}
	}
}

/**
 * @brief Performs image dilation using a specified kernel size.
 *
 * @param inImg Input image.
 * @param outImg Output image (dilated).
 * @param kernelSize Size of the dilation kernel.
 */
void DIP_dilate(const Image *inImg, Image *outImg, int kernelSize) {
	int kernelRadius = kernelSize / 2;

	for (int y = 0; y < inImg->height; ++y) {
		for (int x = 0; x < inImg->width; ++x) {
			uint8_t maxPixelValue = 0;

			for (int ky = -kernelRadius; ky <= kernelRadius; ++ky) {
				for (int kx = -kernelRadius; kx <= kernelRadius; ++kx) {
					int offsetX = x + kx;
					int offsetY = y + ky;

					if (offsetX >= 0 && offsetX < inImg->width && offsetY >= 0
							&& offsetY < inImg->height) {
						int pixelValue = inImg->pixels[offsetY * inImg->width
								+ offsetX];
						if (pixelValue > maxPixelValue) {
							maxPixelValue = pixelValue;
						}
					}
				}
			}

			outImg->pixels[y * inImg->width + x] = maxPixelValue;
		}
	}
}

/**
 * @brief Performs morphological opening (erosion followed by dilation) on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (opened).
 * @param kernelSize Size of the kernel.
 */
void DIP_opening(const Image *inImg, Image *outImg, int kernelSize) {
	// Apply erosion followed by dilation
	Image *temp = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_erode(inImg, temp, kernelSize);
	DIP_dilate(temp, outImg, kernelSize);
}


/**
 * @brief Performs morphological closing (dilation followed by erosion) on the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (closed).
 * @param kernelSize Size of the kernel.
 */
void DIP_closing(const Image *inImg, Image *outImg, int kernelSize) {
	// Apply dilation followed by erosion
	Image *temp = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_dilate(inImg, temp, kernelSize);
	DIP_erode(temp, outImg, kernelSize);
}

/**
 * @brief Identifies and labels connected components in a binary image.
 *
 * @param inImg Input image.
 * @param outImg Output image (labeled components).
 */
void DIP_connectedComponents(const Image *inImg, Image *outImg) {
    int label = 1;

	Image *equivalences = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);

	for(int i=0;i<inImg->size;i++)equivalences->pixels[i]=0x00;
	for(int i=0;i<inImg->size;i++)outImg->pixels[i]=0x00;
    // First pass
    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
            int index = y * inImg->width + x;
            if (inImg->pixels[index] != 0) { // If pixel is part of an object
                int left = (x > 0) ? outImg->pixels[index - 1] : 0;
                int up = (y > 0) ? outImg->pixels[index - inImg->width] : 0;

                if (left && up) { // Both neighbors are labeled
                	outImg->pixels[index] = (left < up) ? left : up; // Assign the minimum label
                    if (left != up) {
                        // Store the equivalence between labels
                        int minLabel = (left < up) ? left : up;
                        int maxLabel = (left > up) ? left : up;
                        equivalences->pixels[maxLabel] = minLabel;
                    }
                } else if (left || up) { // One neighbor is labeled
                	outImg->pixels[index] = (left > 0) ? left : up; // Assign the labeled label
                } else {
                	outImg->pixels[index] = label++; // Assign a new label
                }
            }
        }
    }

    // Second pass to resolve equivalences
    for (int i = 0; i < inImg->width * inImg->height; i++) {
        if (outImg->pixels[i] > 0) {
            int root = outImg->pixels[i];
            while (equivalences->pixels[root] > 0) {
                root = equivalences->pixels[root];
            }
            outImg->pixels[i] = root;
        }
    }
}

/**
 * @brief Applies the Canny edge detection algorithm to the input image.
 *
 * @param inImg Input image.
 * @param outImg Output image (edges).
 * @param lowThreshold Low threshold for edge detection.
 * @param HighThreshold High threshold for edge detection.
 */
void DIP_canny(const Image *inImg, Image *outImg, int lowThreshold,
		int HighThreshold) {
	uint32_t toBeFreed = allocatedSize;
	//1- grayscale conversion
	// image here is grayscale

	//2- Noise Reduction - Gaussian Filter
	Image *noiseReduced = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_gaussianBlur(inImg, noiseReduced, 5, 1);

	// 3 -Gradient Calculation - Sobel Filter
	Image *SOBELLED_ANGLES = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	addFourier(SOBELLED_ANGLES);
	DIP_UTIL_sobelFilter(noiseReduced, SOBELLED_ANGLES);

	// 3- Non Maximum Supression
	Image *non_max_supressed = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_UTIL_nonMaxSupression(SOBELLED_ANGLES, non_max_supressed);

	// 4- Double Thresholding
	Image *double_thresholded = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_UTIL_doubleTreshold(non_max_supressed, lowThreshold, HighThreshold,
			double_thresholded, 0, 255);

	// 5- Hysteresis
	DIP_UTIL_hysteresis(double_thresholded, outImg, 0, 255);

	// free the memory that used in steps of the process.
	allocatedSize = toBeFreed;
}

/**
 * @brief Translates the input image by a specified distance in the X and Y directions.
 *
 * @param inImg Input image.
 * @param outImg Output image (translated).
 * @param dx Translation distance along the X axis.
 * @param dy Translation distance along the Y axis.
 */
void DIP_translate(const Image *inImg, Image *outImg, int dx, int dy) {

	for (int i = 0; i < inImg->size; i++)outImg->pixels[i] = 0x00;

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            int newX = j - dx;
            int newY = i - dy;

            if (newX >= 0 && newX < inImg->width && newY >= 0 && newY < inImg->height) {
                outImg->pixels[newY * inImg->width + newX] = inImg->pixels[i * inImg->width + j];
            }
        }
    }
}

/**
 * @brief Scales the input image by specified factors in the X and Y directions.
 *
 * @param inImg Input image.
 * @param outImg Output image (scaled).
 * @param scaleX Scaling factor along the X axis.
 * @param scaleY Scaling factor along the Y axis.
 */
void DIP_scale(const Image *inImg, Image *outImg, float scaleX, float scaleY) {

	for (int i = 0; i < inImg->size; i++)outImg->pixels[i] = 0x00;

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            int newX = (int)((float)j * (float)scaleX);
            int newY = (int)((float)i * (float)scaleY);

            if (newX >= 0 && newX < inImg->width && newY >= 0 && newY < inImg->height) {
            	outImg->pixels[newY * inImg->width + newX] = inImg->pixels[i * inImg->width + j];
            }
        }
    }
}

/**
 * @brief Rotates the input image by a specified angle.
 *
 * @param inImg Input image.
 * @param outImg Output image (rotated).
 * @param angle Rotation angle in degrees.
 */
void DIP_rotate(const Image *inImg, Image *outImg, float angle) {
    float cosAngle = cos(angle);
    float sinAngle = sin(angle);

	for (int i = 0; i < inImg->size; i++)outImg->pixels[i] = 0x00;

    int centerX = inImg->width / 2;
    int centerY = inImg->height / 2;

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            int newX = (int)((j - centerX) * cosAngle - (i - centerY) * sinAngle) + centerX;
            int newY = (int)((j - centerX) * sinAngle + (i - centerY) * cosAngle) + centerY;

            if (newX >= 0 && newX < inImg->width && newY >= 0 && newY < inImg->height) {
            	outImg->pixels[newY * inImg->width + newX] = inImg->pixels[i * inImg->width + j];
            }
        }
    }
}

/**
 * @brief Shears the input image by specified factors along the X and Y axes.
 *
 * @param inImg Input image.
 * @param outImg Output image (sheared).
 * @param shx Shear factor along the X axis.
 * @param shy Shear factor along the Y axis.
 */
void DIP_shear(const Image *inImg, Image *outImg, float shx, float shy) {

	for (int i = 0; i < inImg->size; i++)outImg->pixels[i] = 0x00;

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            int newX = j + (int)(shx * i);
            int newY = i + (int)(shy * j);

            if (newX >= 0 && newX < inImg->width && newY >= 0 && newY < inImg->height) {
            	outImg->pixels[newY * inImg->width + newX] = inImg->pixels[i * inImg->width + j];
            }
        }
    }
}

/**
 * @brief Inverse shears the input image by specified factors along the X and Y axes.
 *
 * @param inImg Input image.
 * @param outImg Output image (sheared).
 * @param shx Shear factor along the X axis.
 * @param shy Shear factor along the Y axis.
 */
void DIP_inverse_shear(const Image *inImg, Image *outImg, float shx, float shy) {
    // Clear the output image
    for (int i = 0; i < outImg->size; i++) {
        outImg->pixels[i] = 0x00;
    }

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            // Calculate original coordinates
            int originalX = j - (int)(shx * i);
            int originalY = i - (int)(shy * j);

            if (originalX >= 0 && originalX < inImg->width && originalY >= 0 && originalY < inImg->height) {

                outImg->pixels[originalY * inImg->width + originalX] = inImg->pixels[i * inImg->width + j];

            }
        }
    }
}

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
void DIP_warpPerspective(const Image *inImg, Image *outImg, float M[3][3]) {
    // Initialize the output image to black (or any border value)
    for (int i = 0; i < outImg->size; i++) {
        outImg->pixels[i] = 0x00;
    }

    for (int y = 0; y < outImg->height; ++y) {
        for (int x = 0; x < outImg->width; ++x) {
            // Calculate the inverse mapping
            float denom = M[2][0] * x + M[2][1] * y + M[2][2];

            if (denom == 0) {
                continue;  // Avoid division by zero
            }

            float srcX = (M[0][0] * x + M[0][1] * y + M[0][2]) / denom;
            float srcY = (M[1][0] * x + M[1][1] * y + M[1][2]) / denom;

            // Round to nearest integer
            int intSrcX = (int)roundf(srcX);
            int intSrcY = (int)roundf(srcY);

            // Ensure the coordinates are within bounds
            if (intSrcX >= 0 && intSrcX < inImg->width && intSrcY >= 0 && intSrcY < inImg->height) {
                outImg->pixels[y * outImg->width + x] = inImg->pixels[intSrcY * inImg->width + intSrcX];
            }
        }
    }
}

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
void DIP_perspectiveTransform(const Image *inImg, Image *outImg,
                               float h[3][3]) {
    // Initialize output image pixels to 0
    for (int i = 0; i < outImg->size; i++) {
        outImg->pixels[i] = 0x00; // Set to a default value (e.g., black or transparent)
    }

    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
            // Compute the perspective coordinates
            float newX = (h[0][0] * x + h[0][1] * y + h[0][2]) /
                         (h[2][0] * x + h[2][1] * y + h[2][2]);
            float newY = (h[1][0] * x + h[1][1] * y + h[1][2]) /
                         (h[2][0] * x + h[2][1] * y + h[2][2]);

            // Round the new coordinates to the nearest integer
            int intNewX = (int)floor(newX + 0.5);
            int intNewY = (int)floor(newY + 0.5);

            // Check if the new coordinates are within bounds of the output image
            if (intNewX >= 0 && intNewX < outImg->width &&
                intNewY >= 0 && intNewY < outImg->height) {
                outImg->pixels[intNewY * outImg->width + intNewX] =
                    inImg->pixels[y * inImg->width + x];
            }
        }
    }
}


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
void DIP_template_matching(const Image *sourceImg, const Image *templateImg, Image *resultImg) {

    int templateWidth = templateImg->width;
    int templateHeight = templateImg->height;

    // Create buffers for the template and source regions
    float32_t *templateBuffer =  SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
	allocatedSize += templateWidth * templateHeight *sizeof(float32_t);


    float32_t *sourceBuffer =  SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
	allocatedSize += templateWidth * templateHeight *sizeof(float32_t);

    float32_t *sourceFull =  SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
	allocatedSize += sourceImg->width * sourceImg->height *sizeof(float32_t);

    for (int ty = 0; ty < sourceImg->height; ty++) {
        for (int tx = 0; tx < sourceImg->width; tx++) {
        	sourceFull[ty * sourceImg->width + tx] = (float32_t)sourceImg->pixels[(ty) * sourceImg->width + (tx)];
        }
    }

    // Populate the template buffer
    for (int ty = 0; ty < templateHeight; ty++) {
        for (int tx = 0; tx < templateWidth; tx++) {
            templateBuffer[ty * templateWidth + tx] = (float32_t)templateImg->pixels[ty * templateWidth + tx];
        }
    }

    // Calculate template mean and stddev
    float32_t templateMean, templateStdDev;
    arm_mean_f32(templateBuffer, templateWidth * templateHeight, &templateMean);
    arm_std_f32(templateBuffer, templateWidth * templateHeight, &templateStdDev);

    for (int y = 0; y <= sourceImg->height - templateHeight; y++) {
        for (int x = 0; x <= sourceImg->width - templateWidth; x++) {
            // Populate the source buffer for the current region

            for (int ty = 0; ty < templateHeight; ty++) {
                memcpy(sourceBuffer + ty * templateWidth,
                       sourceFull + (y + ty) * sourceImg->width + x,  // Start at (x, y) in the source
                       templateWidth * sizeof(float32_t));  // Correct the size of each row in bytes
            }

            // Calculate source mean and stddev
            float32_t sourceMean, sourceStdDev;
            arm_mean_f32(sourceBuffer, templateWidth * templateHeight, &sourceMean);
            arm_std_f32(sourceBuffer, templateWidth * templateHeight, &sourceStdDev);

            // Calculate normalized cross-correlation
            float32_t correlationSum;
            arm_dot_prod_f32(sourceBuffer, templateBuffer, templateWidth * templateHeight, &correlationSum);

            float32_t normalizedCrossCorrelation = (correlationSum - (templateWidth * templateHeight * sourceMean * templateMean)) /
                                                    (templateWidth * templateHeight * sourceStdDev * templateStdDev);

            // If the normalized cross-correlation exceeds a certain threshold, mark it
            if (normalizedCrossCorrelation > 0.8) { // Adjust threshold as needed
                for (int ty = 0; ty < templateHeight; ty++) {
                    for (int tx = 0; tx < templateWidth; tx++) {
                        resultImg->pixels[(y + ty) * resultImg->width + (x + tx)] = sourceImg->pixels[(y + ty) * resultImg->width + (x + tx)]; // Highlight match
                    }
                }
            }
        }
    }

}


