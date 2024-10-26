#include "dip.h"
#include <esp32/rom/rtc.h>

/*
debugging issues
*/

const int WHT_LED =  22;
#include "Arduino.h"


Image* createImage(uint8_t size, uint8_t format) {

	//test;
	Image *image = new Image;

	if (image == NULL) {
		return NULL;  // Failed to allocate memory
	}

	switch (size) {
	case IMAGE_RES_QQVGA: {
		image->width = IMAGE_RES_QQVGA_Width;
		image->height = IMAGE_RES_QQVGA_Height;
		break;
	}
	case IMAGE_RES_QVGA: {
		image->width = IMAGE_RES_QVGA_Width;
		image->height = IMAGE_RES_QVGA_Height;
		break;
	}
	case IMAGE_RES_480x272: {
		image->width = IMAGE_RES_480x272_Width;
		image->height = IMAGE_RES_480x272_Height;
		break;
	}
	case IMAGE_RES_VGA: {
		image->width = IMAGE_RES_VGA_Width;
		image->height = IMAGE_RES_VGA_Height;
		break;
	}
	default: {
		//should not be go in here
		return NULL;
	}
	}

	switch (format) {
	case IMAGE_FORMAT_RGB565: {
		image->size = image->width * image->height * 2;
		image->format = IMAGE_FORMAT_RGB565;
		break;
	}
	case IMAGE_FORMAT_GRAYSCALE: {
		image->size = image->width * image->height;
		image->format = IMAGE_FORMAT_GRAYSCALE;
		break;
	}
	default: {
		//should not be go in here
		return NULL;
	}
	}

	image->pixels = (uint8_t *)ps_malloc(image->width*image->height * sizeof(uint8_t));

	return image;
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

void DIP_filter2D(const Image *inImg, Image *outImg, int size,
		const float *filter) {

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
					// Calculate the index for the 1D filter array
					int filterIndex = (fy + filtersize / 2) * filtersize + (fx + filtersize / 2);
					sum += ((float) pixelValue * filter[filterIndex]);
				}
			}

			sum = sum < 0 ? 0 : (sum > 255 ? 255 : sum);
			outImg->pixels[y * inImg->width + x] = (uint8_t) sum; // Normalize by dividing by 9 for a 3x3 filter
		}
	}
}


void DIP_gaussianBlur(const Image *inImg, Image *outImg, int filtersize, int sigma) {
    FilterKernel kernel;
    kernel.size = filtersize;
    kernel.data = (float**) malloc(filtersize * sizeof(float*));
    double sum = 0.0;

    // Generate Gaussian kernel
    for (int i = 0; i < filtersize; ++i) {
        kernel.data[i] = (float*) malloc(filtersize * sizeof(float));
        for (int j = 0; j < filtersize; ++j) {
            int x = i - filtersize / 2;
            int y = j - filtersize / 2;
            kernel.data[i][j] = exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma);
            sum += kernel.data[i][j];
        }
    }

    // Normalize the kernel
    for (int i = 0; i < filtersize; ++i) {
        for (int j = 0; j < filtersize; ++j) {
            kernel.data[i][j] /= sum;
        }
    }

    // Apply Gaussian filter with mirror padding
    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
            float blurredPixel = 0.0;
            for (int fy = -filtersize / 2; fy <= filtersize / 2; ++fy) {
                for (int fx = -filtersize / 2; fx <= filtersize / 2; ++fx) {
                    // Mirror padding
                    int mirroredY = y + fy < 0 ? - (y + fy) : (y + fy >= inImg->height ? 2 * inImg->height - (y + fy) - 2 : y + fy);
                    int mirroredX = x + fx < 0 ? - (x + fx) : (x + fx >= inImg->width ? 2 * inImg->width - (x + fx) - 2 : x + fx);

                    int pixelValue = inImg->pixels[mirroredY * inImg->width + mirroredX];
                    blurredPixel += pixelValue * kernel.data[fy + filtersize / 2][fx + filtersize / 2];
                }
            }
            outImg->pixels[y * inImg->width + x] = (uint8_t) (blurredPixel < 0 ? 0 : (blurredPixel > 255 ? 255 : blurredPixel));
        }
    }

    // Free allocated memory for kernel
    for (int i = 0; i < filtersize; ++i) {
        free(kernel.data[i]);
    }
    free(kernel.data);
}

//TODO OZAN
void DIP_medianBlur(const Image *inImg, Image *outImg, int filtersize) {
  
  // Allocate a temporary array for sorting the filter window
  uint8_t *kernel_data = (uint8_t*) malloc(filtersize * filtersize * sizeof(uint8_t));
  
  for (int y = 0; y < inImg->height; y++) {
    for (int x = 0; x < inImg->width; x++) {

      // Collect filter items in the neighborhood
      int index = 0;
      for (int fy = -filtersize / 2; fy <= filtersize / 2; ++fy) {
        for (int fx = -filtersize / 2; fx <= filtersize / 2; ++fx) {
          
          // Handle boundary conditions
          if (y + fy >= 0 && y + fy < inImg->height && x + fx >= 0 && x + fx < inImg->width) {
            kernel_data[index++] = inImg->pixels[(y + fy) * inImg->width + (x + fx)];
          } else {
            // Handle pixels outside of the image as 0
            kernel_data[index++] = 0;
          }
        }
      }

      // Sort the kernel data to find the median
      qsort(kernel_data, filtersize * filtersize, sizeof(uint8_t), compare_uint8);

      // Assign the median value to the output image
      outImg->pixels[y * inImg->width + x] = kernel_data[filtersize * filtersize / 2];
    }
  }
  
  // Free the allocated memory for the kernel data
  free(kernel_data);
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

void DIP_otsuMethod(const Image *inImg, Image *outImg) {

	int histogram[256] = { 0 };

	DIP_histogramFormation(inImg, histogram);

	float prob[256], omega[256]; 
	float myu[256]; 
	float max_sigma, sigma[256]; 
	int i; 
	int threshold; 

	// calculation of probability density
	for (i = 0; i < 256; i++) {
		prob[i] = (float) histogram[i] / (inImg->size);
	}

	// omega & myu generation 
	omega[0] = prob[0];
	myu[0] = 0.0; // 0.0 times prob[0] equals zero 
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



void regionGrowing(Image *inImg, Image *outImg, int y, int x, int threshold) {

	Image *mask = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);

	for(int i=0;i<inImg->size;i++)mask->pixels[i]=0;

    // Define connectivity (8-connected neighborhood)
    int connectivity[8][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}, {0, 1}, {1, -1}, {1, 0}, {1, 1}};

    // Initialize the region with the seed point
    size_t totalSize = inImg->height * inImg->width * 2;  // Size in bytes

    // Allocate memory using malloc
    int *queue = (int*) malloc(totalSize);

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



void DIP_basicSegmentation(const Image *inImg, Image *outImg,
		uint8_t targetValue) {
	for (int i = 0; i < inImg->size; ++i) {
		outImg->pixels[i] = (inImg->pixels[i] == targetValue) ? 255 : 0;
	}
}


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


//***************************13.4 Opening***************************
void DIP_opening(const Image *inImg, Image *outImg, int kernelSize) {

	// Apply erosion followed by dilation
	Image *temp = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_erode(inImg, temp, kernelSize);
	DIP_dilate(temp, outImg, kernelSize);
}

//***************************13.5 Closing************************
void DIP_closing(const Image *inImg, Image *outImg, int kernelSize) {

	// Apply dilation followed by erosion
	Image *temp = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_dilate(inImg, temp, kernelSize);
	DIP_erode(temp, outImg, kernelSize);
}


void DIP_connectedComponents(const Image *inImg, Image *outImg) {
    int label = 1;

	Image *equivalences = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);

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



void DIP_canny2(const Image *inImg, Image *outImg, int lowThreshold,
		int HighThreshold) {


	//1- grayscale conversion
	// image here is grayscale

	//2- Noise Reduction - Gaussian Filter
	Image *noiseReduced = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_gaussianBlur(inImg, noiseReduced, 5, 1);

	// 3 -Gradient Calculation - Sobel Filter
  float *GradX = new float[inImg->size];
  float *GradY = new float[inImg->size];
  Image *SOBELLED_ANGLES = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_UTIL_sobelFilter(noiseReduced, SOBELLED_ANGLES, GradX, GradY);

	// 3- Non Maximum Supression
	Image *non_max_supressed = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_UTIL_nonMaxSuppression(SOBELLED_ANGLES, non_max_supressed, GradX, GradY);

	// 4- Double Thresholding
	Image *double_thresholded = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_UTIL_doubleTreshold(non_max_supressed, lowThreshold, HighThreshold,
			double_thresholded, 0, 255);

	// 5- Hysteresis
	DIP_UTIL_hysteresis(double_thresholded, outImg, 0, 255);

}

void DIP_canny(const Image *inImg, Image *outImg, int lowThreshold,
		int HighThreshold) {


	//1- grayscale conversion
	// image here is grayscale

	//2- Noise Reduction - Gaussian Filter
	Image *noiseReduced = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_gaussianBlur(inImg, noiseReduced, 5, 1.4); // Kernel size 5, sigma 1.4

	// 3 -Gradient Calculation - Sobel Filter
  float *gradX = (float *)ps_malloc(inImg->size * sizeof(float));
  float *gradY = (float *)ps_malloc(inImg->size * sizeof(float));
  Image *gradientMagnitude = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_UTIL_sobelFilter(noiseReduced, gradientMagnitude, gradX, gradY);

  // Step 3: Non-Maximum Suppression
  Image *nonMaxSuppressed = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
  DIP_UTIL_nonMaxSuppression(gradientMagnitude, nonMaxSuppressed, gradX, gradY);

	// 4- Double Thresholding
	Image *double_thresholded = createImage(GET_IMAGE_SIZE(inImg->width), inImg->format);
	DIP_UTIL_doubleTreshold(nonMaxSuppressed, lowThreshold, HighThreshold,
			double_thresholded, 0, 255);

	// 5- Hysteresis
	DIP_UTIL_hysteresis(double_thresholded, outImg, 0, 255);

  // Free intermediate memory
  free(noiseReduced->pixels);
  free(gradX); free(gradY);
  free(gradientMagnitude->pixels);
  free(nonMaxSuppressed->pixels); 
  free(double_thresholded->pixels); 
}


// Function to perform translation
void translate(const Image *inImg, Image *outImg, int dx, int dy) {

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


void DIP_shear(const Image *inImg, Image *outImg, float shx, float shy) {

  // Calculate new dimensions
  int newWidth = (int)(inImg->width + fabs(shx) * inImg->height);
  int newHeight = (int)(inImg->height + fabs(shy) * inImg->width);

  // Update outImg dimensions
  outImg->width = newWidth;
  outImg->height = newHeight;
  outImg->size = newWidth * newHeight;
  outImg->pixels = (uint8_t *)ps_malloc(outImg->width*outImg->height * sizeof(uint8_t));

	for (int i = 0; i < outImg->size; i++)outImg->pixels[i] = 0x00;

    for (int i = 0; i < inImg->height; i++) {
        for (int j = 0; j < inImg->width; j++) {
            int newX = j + (int)round(shx * i);
            int newY = i + (int)round(shy * j);

            // Check bounds for the output image
            if (newX >= 0 && newX < outImg->width && newY >= 0 && newY < outImg->height) {
              outImg->pixels[newY * outImg->width + newX] = inImg->pixels[i * inImg->width + j];
            }
        }
    }
}

void DIP_inverse_shear(const Image *inImg, Image *outImg, float shx, float shy) {
  // Initialize outImg dimensions to match inImg dimensions
  outImg->width = (int)(inImg->width - fabs(shx) * inImg->height);
  outImg->height = (int)(inImg->height - fabs(shy) * inImg->width);
  outImg->size = outImg->width * outImg->height;
  outImg->pixels = (uint8_t *)ps_malloc(outImg->width*outImg->height * sizeof(uint8_t));

  for (int i = 0; i < outImg->size; i++)outImg->pixels[i] = 0x00;

    // Apply backward mapping for inverse shear transformation
    for (int i = 0; i < outImg->height; i++) {
        for (int j = 0; j < outImg->width; j++) {
            // Calculate the source coordinates in the sheared image
            float sourceX = j + shx * i;
            float sourceY = i + shy * j;

            int srcX = (int)round(sourceX);
            int srcY = (int)round(sourceY);

            // Check bounds for the input image
            if (srcX >= 0 && srcX < inImg->width && srcY >= 0 && srcY < inImg->height) {
                outImg->pixels[i * outImg->width + j] = inImg->pixels[srcY * inImg->width + srcX];
            }
        }
    }
}


// DIP UTIL *********************************************************************************************************

void DIP_abs(Image *inImg, Image *outImg, int FFT_SIZE, int quadrant_change){

    float **abs = (float **)ps_malloc(FFT_SIZE * sizeof(float *));
    if (abs == NULL) {
        printf("Memory allocation failed!\n");
        return;
    }

    for (int i = 0; i < FFT_SIZE; i++) {
        abs[i] = (float *)ps_malloc(FFT_SIZE * sizeof(float));
        if (abs[i] == NULL) {
            printf("Memory allocation for row %d failed!\n", i);
            // Free previously allocated rows to avoid memory leaks
            for (int j = 0; j < i; j++) {
                free(abs[j]);
            }
            free(abs);
            return;
        }
    }
  float max=0;
  for (int i = 0; i < FFT_SIZE; i++) {
      for (int j = 0; j < FFT_SIZE * 2; j += 2) {
        if(quadrant_change){
          abs[i][j/2] = log(sqrtf(outImg->fourier[i][j] * outImg->fourier[i][j] + outImg->fourier[i][j + 1] * outImg->fourier[i][j + 1]) + 1);
        }
        else{
          abs[i][j/2] = sqrtf(outImg->fourier[i][j] * outImg->fourier[i][j] + outImg->fourier[i][j + 1] * outImg->fourier[i][j + 1]);
        }

          if(max < abs[i][j/2]){
            max = abs[i][j/2];
          }
      }
  }

  //allocating fourier memory to img pointer
  uint8_t *test = (uint8_t *)ps_malloc(FFT_SIZE * FFT_SIZE * sizeof(uint8_t *));


  for(int i= 0;i<outImg->size;i++){
    outImg->pixels[i]=0xff;
  }

  for(int i=0;i<FFT_SIZE;i++){
    for(int j=0;j<FFT_SIZE;j++){
      outImg->pixels[i * outImg->width + j] = (255 * (abs[i][j] / max));
      test[i * FFT_SIZE + j] = outImg->pixels[i * outImg->width + j];
    }
  }

  if(quadrant_change){
    int halfHeight = FFT_SIZE/2;
    int halfWidth = FFT_SIZE/2;

    for (int y = 0; y < halfHeight; y++) {
      for (int x = 0; x < halfWidth; x++) {
        // Swap quadrants
        int quad1Index = y * FFT_SIZE + x;
        int quad2Index = quad1Index + halfWidth;
        int quad3Index = (y + halfHeight) * FFT_SIZE + x;
        int quad4Index = quad3Index + halfWidth;

        uint32_t temp = test[quad1Index];
        test[quad1Index] = test[quad4Index];
        test[quad4Index] = temp;

        temp = test[quad2Index];
        test[quad2Index] = test[quad3Index];
        test[quad3Index] = temp;
      }
    }

    for(int i=0;i<FFT_SIZE;i++){
      for(int j=0;j<FFT_SIZE;j++){
        outImg->pixels[i * outImg->width + j] = test[i * FFT_SIZE + j];
      }
    }
  }
}


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

void DIP_fourier(Image *inImg, Image *outImg, int FFT_SIZE){

  //allocating fourier memory to img pointer
  outImg->fourier = (float **)ps_malloc(FFT_SIZE * sizeof(float *));

  // Allocate memory for each row
  for (int i = 0; i < FFT_SIZE; i++) {
    outImg->fourier[i] = (float *)ps_malloc(FFT_SIZE * sizeof(float) * 2);
  }

  //fill the array
  for(int i=0;i<FFT_SIZE;i++){
    for(int j=0;j<FFT_SIZE*2;j++){
      if(j%2 == 0){
        outImg->fourier[i][j] = outImg->pixels[i*FFT_SIZE + j/2];
      }
      else{
        outImg->fourier[i][j] = 0.0f;
      }

    }

  }
  // Initialize the FFT library
  dsps_fft2r_init_fc32(NULL, CONFIG_DSP_MAX_FFT_SIZE);

  // Perform FFT on rows
  for (int i = 0; i < FFT_SIZE; i++) {
      dsps_fft2r_fc32(outImg->fourier[i], FFT_SIZE);  // Apply FFT on each row
      dsps_bit_rev_fc32(outImg->fourier[i], FFT_SIZE);  // Bit reversal
  }

  // Transpose the matrix to apply FFT on columns
    float **transposed_data = (float **)ps_malloc(FFT_SIZE * sizeof(float *));

    for (int i = 0; i < FFT_SIZE; i++) {
        transposed_data[i] = (float *)ps_malloc(FFT_SIZE * sizeof(float) * 2);
    }


  for (int i = 0; i < FFT_SIZE; i++) {
      for (int j = 0; j < FFT_SIZE * 2; j++) {
          transposed_data[j / 2][i * 2 + (j % 2)] = outImg->fourier[i][j];  // Transpose real and imaginary parts separately
      }
  }

  // Perform FFT on transposed columns (which are now rows)
  for (int i = 0; i < FFT_SIZE; i++) {
      dsps_fft2r_fc32(transposed_data[i], FFT_SIZE);  // Apply FFT on each row (which corresponds to columns of original)
      dsps_bit_rev_fc32(transposed_data[i], FFT_SIZE);  // Bit reversal
  }

  // Transpose back to get the final 2D FFT result
  for (int i = 0; i < FFT_SIZE; i++) {
      for (int j = 0; j < FFT_SIZE * 2; j++) {
          outImg->fourier[j / 2][i * 2 + (j % 2)] = transposed_data[i][j];  // Transpose back real and imaginary parts

      }
  }
}

void DIP_fourier_inv(Image *inImg, Image *outImg, int FFT_SIZE){

    // Step 7: Perform inverse FFT (iFFT) using the same approach
    // First, conjugate the input for iFFT
    for (int i = 0; i < FFT_SIZE; i++) {
        for (int j = 0; j < FFT_SIZE; j++) {
            inImg->fourier[i][j * 2 + 1] = -inImg->fourier[i][j * 2 + 1];  // Negate imaginary part for conjugation
        }
    }

    // Perform inverse FFT on rows
    for (int i = 0; i < FFT_SIZE; i++) {
        dsps_fft2r_fc32(inImg->fourier[i], FFT_SIZE);  // Perform FFT (used for iFFT)
        dsps_bit_rev_fc32(inImg->fourier[i], FFT_SIZE);  // Bit reversal
    }

    //allocating fourier memory to img pointer
    float **transposed_data = (float **)ps_malloc(FFT_SIZE * sizeof(float *));

    // Allocate memory for each row
    for (int i = 0; i < FFT_SIZE; i++) {
      transposed_data[i] = (float *)ps_malloc(FFT_SIZE * sizeof(float) * 2);
    }
    // Transpose the matrix for column-wise FFT
    for (int i = 0; i < FFT_SIZE; i++) {
        for (int j = 0; j < FFT_SIZE; j++) {
            transposed_data[j][i * 2] = inImg->fourier[i][j * 2];      // Transpose real parts
            transposed_data[j][i * 2 + 1] = inImg->fourier[i][j * 2 + 1];  // Transpose imaginary parts
        }
    }

    // Perform inverse FFT on transposed columns (now rows)
    for (int i = 0; i < FFT_SIZE; i++) {
        dsps_fft2r_fc32(transposed_data[i], FFT_SIZE);  // Perform FFT (used for iFFT)
        dsps_bit_rev_fc32(transposed_data[i], FFT_SIZE);  // Bit reversal
    }

    // Transpose back
    for (int i = 0; i < FFT_SIZE; i++) {
        for (int j = 0; j < FFT_SIZE; j++) {
            inImg->fourier[j][i * 2] = transposed_data[i][j * 2];      // Transpose real parts back
            inImg->fourier[j][i * 2 + 1] = transposed_data[i][j * 2 + 1];  // Transpose imaginary parts back
        }
    }

    // Conjugate again after iFFT and scale by N_SAMPLES * N_SAMPLES
    for (int i = 0; i < FFT_SIZE; i++) {
        for (int j = 0; j < FFT_SIZE * 2; j++) {
            inImg->fourier[i][j] /= (FFT_SIZE * FFT_SIZE);  // Scale the result
        }
    }
}

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

void DIP_sepFilter2D(const Image *inImg, Image *outImg, int kernelSizeX, float *kernelX,
                     int kernelSizeY, float *kernelY, double delta) {

    Image *tempImg = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
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
 * @brief Creates a Laplacian filter.
 *
 * This function generates a Laplacian filter of the specified size.
 *
 * @param[in] ksize Size of the filter.
 * @return A pointer to a 2D array representing the Laplacian filter.
 */
float **createLaplacianFilter(int ksize) {
    // Check if ksize is valid (positive and odd)
    if (ksize <= 0 || ksize % 2 == 0) {
        return NULL; // Invalid ksize, return NULL
    }

    // Allocate memory for the filter
    float **filter = (float **)ps_malloc(ksize * sizeof(float *));
    for (int i = 0; i < ksize; i++) {
        filter[i] = (float *)ps_malloc(ksize * sizeof(float));
    }

    // Create the Laplacian filter based on ksize
    int center = ksize / 2;
    for (int y = 0; y < ksize; y++) {
        for (int x = 0; x < ksize; x++) {
            if (y == center && x == center) {
                filter[y][x] = (float)(ksize * ksize - 1); // Central value
            } else {
                filter[y][x] = -1.0f; // Surrounding values
            }
        }
    }

    return filter;
}

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
            dst->pixels[y * outputWidth + x] = blurredImg->pixels[2 * y * src->width + 2 * x];
        }
    }

    dst->height = outputHeight;
    dst->width = outputWidth;
    dst->size = outputWidth* outputHeight;
}

void DIP_pyrUp(const Image *inImg, Image *outImg, Size *dstSize, int borderType, int filtersize, int sigma) {
    // Default output size calculation if dstSize is not set
    if (dstSize->width == 0 || dstSize->height == 0) {
        dstSize->width = inImg->width * 2;
        dstSize->height = inImg->height * 2;
    }

    outImg->width = dstSize->width;
    outImg->height = dstSize->height;
    outImg->size = dstSize->width *dstSize->height;
    outImg->pixels = (uint8_t *)ps_malloc(outImg->width*outImg->height * sizeof(uint8_t));

    Image *temp = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
    temp->width = dstSize->width;
    temp->height = dstSize->height;
    temp->size = dstSize->width *dstSize->height;
    temp->pixels = (uint8_t *)ps_malloc(temp->width*temp->height * sizeof(uint8_t));

    for(int i=0;i<temp->size;i++){
      temp->pixels[i]=0x00;
    }
    // Upsample by inserting zeros
    for (int y = 0; y < inImg->height; y++) {
        for (int x = 0; x < inImg->width; x++) {
        	temp->pixels[(y * 2) * outImg->width + (x * 2)] = inImg->pixels[y * inImg->width + x];
        }
    }

    // Apply Gaussian blur to the upsampled image using your existing function
    DIP_gaussianBlur(temp, outImg, filtersize, sigma);
}

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

void DIP_warpPerspective(const Image *inImg, Image *outImg, float transform[3][3]) {

    // Transform corners of the input image to find the output bounds
    float corners[4][2] = {
        {0, 0},
        {inImg->width - 1, 0},
        {0, inImg->height - 1},
        {inImg->width - 1, inImg->height - 1}
    };
    
    float minX = FLT_MAX, minY = FLT_MAX;
    float maxX = -FLT_MAX, maxY = -FLT_MAX;

    for (int i = 0; i < 4; i++) {
        float x = corners[i][0];
        float y = corners[i][1];
        float w = transform[2][0] * x + transform[2][1] * y + transform[2][2];
        float newX = (transform[0][0] * x + transform[0][1] * y + transform[0][2]) / w;
        float newY = (transform[1][0] * x + transform[1][1] * y + transform[1][2]) / w;

        if (newX < minX) minX = newX;
        if (newY < minY) minY = newY;
        if (newX > maxX) maxX = newX;
        if (newY > maxY) maxY = newY;
    }

    // Calculate new width and height for the output image
    int newWidth = (int)(maxX - minX + 120);
    int newHeight = (int)(maxY - minY + 120);

    outImg->width = newWidth;
    outImg->height = newHeight;
    outImg->size = newWidth * newHeight;
    outImg->pixels = (uint8_t *)ps_malloc(outImg->width * outImg->height * sizeof(uint8_t));

    // Initialize the output image with a background color (black)
    for (int i = 0; i < outImg->size; i++) {
        outImg->pixels[i] = 0x00;
    }

    // Adjust the transformation to center the image within the output bounds
    float offsetX = -minX;
    float offsetY = -minY;

    // Apply inverse mapping to populate the output image
    for (int y = 0; y < outImg->height; y++) {
        for (int x = 0; x < outImg->width; x++) {
            
            // Apply the inverse of the perspective transformation
            float w = transform[2][0] * (x - offsetX) + transform[2][1] * (y - offsetY) + transform[2][2];
            float srcX = (transform[0][0] * (x - offsetX) + transform[0][1] * (y - offsetY) + transform[0][2]) / w;
            float srcY = (transform[1][0] * (x - offsetX) + transform[1][1] * (y - offsetY) + transform[1][2]) / w;

            // Check if the computed source coordinates are within bounds of the input image
            if (srcX >= 0 && srcX < inImg->width && srcY >= 0 && srcY < inImg->height) {
                int srcXInt = (int)round(srcX);
                int srcYInt = (int)round(srcY);

                // Map the pixel to the output image
                outImg->pixels[y * outImg->width + x] = inImg->pixels[srcYInt * inImg->width + srcXInt];
            }
        }
    }
}

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

void DIP_template_matching2(const Image *sourceImg, const Image *templateImg, Image *resultImg) {

    int templateWidth = templateImg->width;
    int templateHeight = templateImg->height;

    // Allocate buffers for the template and source regions
    float *templateBuffer = (float *)ps_malloc(templateWidth * templateHeight * sizeof(float));

    float *sourceBuffer = (float *)ps_malloc(templateWidth * templateHeight * sizeof(float));

    float *sourceFull = (float *)ps_malloc(sourceImg->width * sourceImg->height * sizeof(float));

    // Copy source image pixels to sourceFull buffer
    for (int ty = 0; ty < sourceImg->height; ty++) {
        for (int tx = 0; tx < sourceImg->width; tx++) {
            sourceFull[ty * sourceImg->width + tx] = (float)sourceImg->pixels[ty * sourceImg->width + tx];
        }
    }

    // Populate the template buffer
    for (int ty = 0; ty < templateHeight; ty++) {
        for (int tx = 0; tx < templateWidth; tx++) {
            templateBuffer[ty * templateWidth + tx] = (float)templateImg->pixels[ty * templateWidth + tx];
        }
    }

    // Calculate template mean and stddev manually
    float templateMean = 0, templateStdDev = 0;
    for (int i = 0; i < templateWidth * templateHeight; i++) {
        templateMean += templateBuffer[i];
    }
    templateMean /= (templateWidth * templateHeight);

    for (int i = 0; i < templateWidth * templateHeight; i++) {
        templateStdDev += powf(templateBuffer[i] - templateMean, 2);
    }
    templateStdDev = sqrtf(templateStdDev / (templateWidth * templateHeight));

    // Perform template matching
    for (int y = 0; y <= sourceImg->height - templateHeight; y++) {
        for (int x = 0; x <= sourceImg->width - templateWidth; x++) {

            // Populate the source buffer for the current region
            for (int ty = 0; ty < templateHeight; ty++) {
                memcpy(sourceBuffer + ty * templateWidth,
                       sourceFull + (y + ty) * sourceImg->width + x,
                       templateWidth * sizeof(float));
            }

            // Calculate source mean and stddev manually
            float sourceMean = 0, sourceStdDev = 0;
            for (int i = 0; i < templateWidth * templateHeight; i++) {
                sourceMean += sourceBuffer[i];
            }
            sourceMean /= (templateWidth * templateHeight);

            for (int i = 0; i < templateWidth * templateHeight; i++) {
                sourceStdDev += powf(sourceBuffer[i] - sourceMean, 2);
            }
            sourceStdDev = sqrtf(sourceStdDev / (templateWidth * templateHeight));

            // Calculate normalized cross-correlation manually
            float correlationSum = 0;
            for (int i = 0; i < templateWidth * templateHeight; i++) {
                correlationSum += (sourceBuffer[i] - sourceMean) * (templateBuffer[i] - templateMean);
            }
            float normalizedCrossCorrelation = correlationSum / (templateWidth * templateHeight * sourceStdDev * templateStdDev);

            // If the normalized cross-correlation exceeds a certain threshold, mark it
            if (normalizedCrossCorrelation > 0.8) { // Adjust threshold as needed
                for (int ty = 0; ty < templateHeight; ty++) {
                    for (int tx = 0; tx < templateWidth; tx++) {
                        resultImg->pixels[(y + ty) * resultImg->width + (x + tx)] =
                            sourceImg->pixels[(y + ty) * resultImg->width + (x + tx)]; // Highlight match
                    }
                }
            }
        }
    }
}


void DIP_TemplateMatch(const Image *sourceImg, const Image *templateImg, Image *resultImg) {
    int templateWidth = templateImg->width;
    int templateHeight = templateImg->height;

    float *templateBuffer = (float *)ps_malloc(templateWidth * templateHeight * sizeof(float));
    float *sourceBuffer = (float *)ps_malloc(templateWidth * templateHeight * sizeof(float));
    float *sourceFull = (float *)ps_malloc(sourceImg->width * sourceImg->height * sizeof(float));

    if (!templateBuffer || !sourceBuffer || !sourceFull) {
        printf("Memory allocation failed.\n");
        return;
    }

    // Normalize the source image pixels
    for (int ty = 0; ty < sourceImg->height; ty++) {
        for (int tx = 0; tx < sourceImg->width; tx++) {
            sourceFull[ty * sourceImg->width + tx] = (float)sourceImg->pixels[ty * sourceImg->width + tx] / 255.0;
        }
    }

    // Populate and normalize the template buffer
    for (int ty = 0; ty < templateHeight; ty++) {
        for (int tx = 0; tx < templateWidth; tx++) {
            templateBuffer[ty * templateWidth + tx] = (float)templateImg->pixels[ty * templateWidth + tx] / 255.0;
        }
    }

    // Calculate template mean and standard deviation manually
    float templateMean = 0, templateStdDev = 0;
    for (int i = 0; i < templateWidth * templateHeight; i++) {
        templateMean += templateBuffer[i];
    }
    templateMean /= (templateWidth * templateHeight);

    for (int i = 0; i < templateWidth * templateHeight; i++) {
        templateStdDev += powf(templateBuffer[i] - templateMean, 2);
    }
    templateStdDev = sqrtf(templateStdDev / (templateWidth * templateHeight));

    // Perform template matching
    for (int y = 0; y <= sourceImg->height - templateHeight; y++) {
        for (int x = 0; x <= sourceImg->width - templateWidth; x++) {

            // Copy current region of the source image into sourceBuffer
            for (int ty = 0; ty < templateHeight; ty++) {
                memcpy(sourceBuffer + ty * templateWidth,
                       sourceFull + (y + ty) * sourceImg->width + x,
                       templateWidth * sizeof(float));
            }

            // Calculate source mean and stddev manually
            float sourceMean = 0, sourceStdDev = 0;
            for (int i = 0; i < templateWidth * templateHeight; i++) {
                sourceMean += sourceBuffer[i];
            }
            sourceMean /= (templateWidth * templateHeight);

            for (int i = 0; i < templateWidth * templateHeight; i++) {
                sourceStdDev += powf(sourceBuffer[i] - sourceMean, 2);
            }
            sourceStdDev = sqrtf(sourceStdDev / (templateWidth * templateHeight));

            // Calculate normalized cross-correlation
            float correlationSum = 0;
            for (int i = 0; i < templateWidth * templateHeight; i++) {
                correlationSum += (sourceBuffer[i] - sourceMean) * (templateBuffer[i] - templateMean);
            }
            float normalizedCrossCorrelation = correlationSum / (templateWidth * templateHeight * sourceStdDev * templateStdDev);

            // If NCC exceeds threshold, highlight the match
            if (normalizedCrossCorrelation > 0.8) {  // Adjust threshold if needed
                for (int ty = 0; ty < templateHeight; ty++) {
                    for (int tx = 0; tx < templateWidth; tx++) {
                        resultImg->pixels[(y + ty) * resultImg->width + (x + tx)] =
                            sourceImg->pixels[(y + ty) * sourceImg->width + (x + tx)];
                    }
                }
            }
        }
    }

    // Free buffers
    free(templateBuffer);
    free(sourceBuffer);
    free(sourceFull);
}




void sobelUtilizedMirrorPadding(const Image *inImg, int *pixelValue, int y,
		int x, int fx, int fy) {

	//sol ust kose
	if (y == 0 && x == 0) {
		if (fx == 0 && fy == 0) { // sol ust kose capraz
			*pixelValue = inImg->pixels[(y + fy) * inImg->width + (x + fx)];
		} else if (fx == 0) { //solda
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width + (x + fx)];
		} else if (fy == 0) {
			*pixelValue = inImg->pixels[(y + fy) * inImg->width + (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	//sol alt kose
	else if (y == inImg->height - 1 && x == 0) {
		if (fx == 0 && fy == 2) { // sol alt kose capraz
			*pixelValue = inImg->pixels[(y + fy - 1 - 1) * inImg->width
					+ (x + fx)];
		} else if (fx == 0) { //solda
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width + (x + fx)];
		} else if (fy == 2) { //asagida
			*pixelValue = inImg->pixels[(y + fy - 1 - 1) * inImg->width
					+ (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	//sag ust kose
	else if (y == 0 && x == inImg->width - 1) {
		if (fx == 2 && fy == 0) { // sag ust kose capraz
			*pixelValue = inImg->pixels[(y + fy) * inImg->width
					+ (x + fx - 1 - 1)];
		} else if (fx == 2) { //sagda
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1 - 1)];
		} else if (fy == 0) { //yukarda
			*pixelValue = inImg->pixels[(y + fy) * inImg->width + (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	//sag alt kose
	else if (y == inImg->height - 1 && x == inImg->width - 1) {
		if (fx == 2 && fy == 2) { // sag alt kose capraz
			*pixelValue = inImg->pixels[(y + fy - 1 - 1) * inImg->width
					+ (x + fx - 1 - 1)];
		} else if (fx == 2) { //sagda
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1 - 1)];
		} else if (fy == 2) { //asagida
			*pixelValue = inImg->pixels[(y + fy - 1 - 1) * inImg->width
					+ (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	// ust
	else if (y == 0) {
		if (fy == 0) {
			*pixelValue = inImg->pixels[(y + fy) * inImg->width + (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}

	}
	//sag
	else if (x == inImg->width - 1) {
		if (fx == 2) {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1 - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}

	//alt
	else if (y == inImg->height - 1) {
		if (fy == 2) {
			*pixelValue = inImg->pixels[(y + fy - 1 - 1) * inImg->width
					+ (x + fx - 1)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	//sol
	else if (x == 0) {
		if (fx == 0) {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width + (x + fx)];
		} else {
			*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width
					+ (x + fx - 1)];
		}
	}
	//normal
	else {
		*pixelValue = inImg->pixels[(y + fy - 1) * inImg->width + (x + fx - 1)];
	}

}

void DIP_UTIL_sobelFilter(const Image *inImg, Image *outImg, float *GradX, float *GradY) {
	int sobelX[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
	int sobelY[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };


	float gradMax = 0;
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
					sobelUtilizedMirrorPadding(inImg, &pixelValue, y, x, fx,
							fy);
					gradientX += pixelValue * sobelX[fy][fx];
					gradientY += pixelValue * sobelY[fy][fx];
				}
			}

			// Calculate gradient magnitude
			GradY[y * inImg->width + x] = sqrt(
					(float) (gradientX * gradientX + gradientY * gradientY));

			gradMax =
					(GradY[y * inImg->width + x] > gradMax) ?
							GradY[y * inImg->width + x] : gradMax;

			GradX[y * inImg->width + x] = (float) atan2f(gradientY,
					gradientX);
		}
	}

    // Second pass to normalize and store results in outImg
    if (gradMax > 0) { // Avoid division by zero
        for (int i = 0; i < outImg->size; i++) {
            outImg->pixels[i] = (uint8_t) ((GradY[i] / gradMax) * 255);
        }
    } else {
        // Handle case when gradMax is zero (e.g., blank image)
        memset(outImg->pixels, 0, outImg->size);
    }

}

void mirrorPadding(Image *inImg, int i, int j, float *r, float *q,
		int condition) {

	uint8_t *G = inImg->pixels;

	if (condition == 0) { //

		if (j == inImg->width - 1) {
			*r = G[i * inImg->width + (j - 1)];
			*q = G[i * inImg->width + (j)] - 1;
		} else if (j == 0) {
			*r = G[i * inImg->width + (j)] - 1;
			*q = G[i * inImg->width + (j + 1)];
		} else {
			*r = G[i * inImg->width + (j - 1)];
			*q = G[i * inImg->width + (j + 1)];
		}

	} else if (condition == 1) {

		//sol ust kose
		if (i == 0 && j == 0) {
			*r = G[(i) * inImg->width + (j + 1)];
			*q = G[(i + 1) * inImg->width + (j)] - 1;
		}
		//sol alt kose
		else if (i == inImg->height - 1 && j == 0) {
			*r = G[(i - 1) * inImg->width + (j + 1)];
			*q = G[(i) * inImg->width + (j)] - 1;
		}
		//sag ust kose
		else if (i == 0 && j == inImg->width - 1) {
			*r = G[(i) * inImg->width + (j)] - 1;
			*q = G[(i + 1) * inImg->width + (j - 1)];
			;
		}
		//sag alt kose
		else if (i == inImg->height - 1 && j == inImg->width - 1) {
			*r = G[(i - 1) * inImg->width + (j)] - 1;
			*q = G[(i) * inImg->width + (j)];
			;
		}
		// ust
		else if (i == 0) {
			*r = G[(i) * inImg->width + (j + 1)] - 1;
			*q = G[(i + 1) * inImg->width + (j - 1)];
			;
		}
		//sag
		else if (j == inImg->width - 1) {
			*r = G[(i - 1) * inImg->width + (j)] - 1;
			*q = G[(i + 1) * inImg->width + (j - 1)];
			;
		}

		//alt
		else if (i == inImg->height - 1) {
			*r = G[(i - 1) * inImg->width + (j + 1)];
			*q = G[(i) * inImg->width + (j - 1)] - 1;
		}
		//sol
		else if (j == 0) {
			*r = G[(i - 1) * inImg->width + (j + 1)];
			*q = G[(i + 1) * inImg->width + (j)] - 1;
		}
		//normal
		else {
			*r = G[(i - 1) * inImg->width + (j + 1)];
			*q = G[(i + 1) * inImg->width + (j - 1)];
		}
	} else if (condition == 2) {

		//ust
		if (i == 0) {
			*r = G[(i) * inImg->width + j] - 1;
			*q = G[(i + 1) * inImg->width + j];
		}
		//alt
		else if (i == inImg->height - 1) {
			*r = G[(i - 1) * inImg->width + j];
			*q = G[(i) * inImg->width + j] - 1;
		}
		//normal
		else {
			*r = G[(i - 1) * inImg->width + j];
			*q = G[(i + 1) * inImg->width + j];
		}

	} else {

		//sol ust kose
		if (i == 0 && j == 0) {
			*r = G[(i + 1) * inImg->width + (j + 1)];
			*q = G[(i) * inImg->width + (j)] - 1;
		}
		//sol alt kose
		else if (i == inImg->height - 1 && j == 0) {
			*r = G[(i) * inImg->width + (j + 1)];
			*q = G[(i - 1) * inImg->width + (j)] - 1;
		}
		//sag ust kose
		else if (i == 0 && j == inImg->width - 1) {
			*r = G[(i + 1) * inImg->width + (j)] - 1;
			*q = G[(i) * inImg->width + (j - 1)];
		}
		//sag alt kose
		else if (i == inImg->height - 1 && j == inImg->width - 1) {
			*r = G[(i) * inImg->width + (j)] - 1;
			*q = G[(i - 1) * inImg->width + (j - 1)];
		}
		// ust
		else if (i == 0) {
			*r = G[(i + 1) * inImg->width + (j + 1)];
			*q = G[(i) * inImg->width + (j - 1)] - 1;
		}
		//sag
		else if (j == inImg->width - 1) {
			*r = G[(i + 1) * inImg->width + (j)] - 1;
			*q = G[(i - 1) * inImg->width + (j - 1)];
		}

		//alt
		else if (i == inImg->height - 1) {
			*r = G[(i) * inImg->width + (j + 1)] - 1;
			*q = G[(i - 1) * inImg->width + (j - 1)];
		}
		//sol
		else if (j == 0) {
			*r = G[(i + 1) * inImg->width + (j + 1)];
			*q = G[(i - 1) * inImg->width + (j)] - 1;
		}
		//normal
		else {
			*r = G[(i + 1) * inImg->width + (j + 1)];
			*q = G[(i - 1) * inImg->width + (j - 1)];
		}

	}

}

void DIP_UTIL_nonMaxSuppression(const Image *inImg, Image *nonmax, float *GradX, float *GradY) {
    uint8_t *G = inImg->pixels;   // Gradient magnitudes
    float *theta = GradX;      // Gradient angles (in radians)

    // Non-maximum suppression
    for (int i = 1; i < inImg->height - 1; i++) { // Avoid borders for simplicity
        for (int j = 1; j < inImg->width - 1; j++) {
            float q = 255.0, r = 255.0;
            float angle = theta[i * inImg->width + j] * (180.0 / M_PI); // Convert to degrees

            // Map angle to 0, 45, 90, or 135 degrees
            if (angle < 0) angle += 180.0;

            // Determine neighboring pixels in gradient direction
            if ((0 <= angle && angle < 22.5) || (157.5 <= angle && angle <= 180)) { // 0 degrees
                q = G[i * inImg->width + (j + 1)];
                r = G[i * inImg->width + (j - 1)];
            } else if (22.5 <= angle && angle < 67.5) { // 45 degrees
                q = G[(i + 1) * inImg->width + (j - 1)];
                r = G[(i - 1) * inImg->width + (j + 1)];
            } else if (67.5 <= angle && angle < 112.5) { // 90 degrees
                q = G[(i + 1) * inImg->width + j];
                r = G[(i - 1) * inImg->width + j];
            } else if (112.5 <= angle && angle < 157.5) { // 135 degrees
                q = G[(i - 1) * inImg->width + (j - 1)];
                r = G[(i + 1) * inImg->width + (j + 1)];
            }

            // Suppress non-maxima
            if (G[i * inImg->width + j] >= q && G[i * inImg->width + j] >= r) {
                nonmax->pixels[i * inImg->width + j] = G[i * inImg->width + j];
            } else {
                nonmax->pixels[i * inImg->width + j] = 0;
            }
        }
    }
}

void DIP_UTIL_doubleTreshold(Image *inImg, float lowThresholdRatio,
		float highThresholdRatio, Image *outImg, int weak, int strong) {
	float highThreshold = 0;
	float lowThreshold = 0;
	weak = 25;
	strong = 255;

	for (int i = 0; i < inImg->size; i++) {
		if (inImg->pixels[i] > highThreshold)
			highThreshold = inImg->pixels[i];
	}

	highThreshold = highThresholdRatio;
	lowThreshold = lowThresholdRatio;

	// Double thresholding
	for (int i = 0; i < inImg->height; i++) {
		for (int j = 0; j < inImg->width; j++) {
			if (inImg->pixels[i * inImg->width + j] >= highThreshold) {
				outImg->pixels[i * inImg->width + j] = strong;
			} else if (inImg->pixels[i * inImg->width + j] < lowThreshold) {
				outImg->pixels[i * inImg->width + j] = 0;
			} else {
				outImg->pixels[i * inImg->width + j] = weak;
			}
		}
	}
}

void DIP_UTIL_hysteresis(Image *inImg, Image *outImg, int weak, int strong) {
	int i, j;
	weak = 25;
	strong = 255;

	for (int i = 0; i < inImg->size; i++)
		outImg->pixels[i] = inImg->pixels[i];

	// Iterate through the image pixels
	for (i = 1; i < inImg->height - 1; i++) {
		for (j = 1; j < inImg->width - 1; j++) {
			// Check if the current pixel is weak
			if (outImg->pixels[i * inImg->width + j] == weak) {

				// Check neighboring pixels for strong edges
				if (outImg->pixels[(i + 1) * inImg->width + (j - 1)] == strong
						|| outImg->pixels[(i + 1) * inImg->width + j] == strong
						|| outImg->pixels[(i + 1) * inImg->width + (j + 1)]
								== strong
						|| outImg->pixels[i * inImg->width + (j - 1)] == strong
						|| outImg->pixels[i * inImg->width + (j + 1)] == strong
						|| outImg->pixels[(i - 1) * inImg->width + (j - 1)]
								== strong
						|| outImg->pixels[(i - 1) * inImg->width + j] == strong
						|| outImg->pixels[(i - 1) * inImg->width + (j + 1)]
								== strong) {
					outImg->pixels[i * inImg->width + j] = strong;
				} else {
					outImg->pixels[i * inImg->width + j] = 0;
				}
			}
		}
	}
}

int compare_uint8(const void *a, const void *b) {
	return (*(uint8_t*) a - *(uint8_t*) b);
}

unsigned int orderof2(const unsigned int value) {
	return (unsigned int) 1
			<< ((sizeof(unsigned int) * CHAR_BIT) - __builtin_clz(value) - 1);
}

void getGaborKernel(float filter[][3], int size, double sigma, double theta, double lambda, double gamma, double psi) {
    int center = size / 2;
    double sum = 0.0;

    for (int y = -center; y <= center; ++y) {
        for (int x = -center; x <= center; ++x) {
            double xTheta = x * cos(theta) + y * sin(theta);
            double yTheta = -x * sin(theta) + y * cos(theta);

            // Gabor kernel equation
            double expPart = exp(-0.5 * (pow(xTheta, 2) + pow(gamma * yTheta, 2)) / (sigma * sigma));
            double cosPart = cos(2 * M_PI * xTheta / lambda + psi);

            filter[y + center][x + center] = (float)(expPart * cosPart);
            sum += filter[y + center][x + center];
        }
    }

    // Normalize the filter
    if (sum != 0) {
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                filter[y][x] = filter[y][x] / (sum);
            }
        }
    }
}

void getGaussianKernel(float filter[][3], int size, double sigma) {
    double sum = 0.0;
    int center = size / 2;

    for (int y = -size / 2; y <= size / 2; ++y) {
        for (int x = -size / 2; x <= size / 2; ++x) {
            double value = (1.0 / (2.0 * M_PI * sigma * sigma)) * exp(-(x * x + y * y) / (2 * sigma * sigma));
            filter[(int)(y + center)][(int)(x + center)] = (float)value;
            sum += filter[(int)(y + center)][(int)(x + center)];
        }
    }

    // Normalize the filter
    if (sum != 0) {
        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                filter[y][x] /= sum;
            }
        }
    }
}