/*
 * DIP.c
 *
 *  Created on: Mar 28, 2024
 *      Author: ozand
 */
#include "DIP.h"
#include "DIP_UTIL.h"
uint32_t allocatedSize = 0;

// Macro for comparing float values
#define COMPARE_FLOAT(a, b) (fabs((a) - (b)) < FLT_EPSILON ? 0 : ((a) < (b) ? -1 : 1))

Image* createImage(uint8_t size, uint8_t format) {

	//test;
	Image *image = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
	allocatedSize += sizeof(Image);

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

	image->pixels = SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize;
	allocatedSize += image->size;

	return image;
}

void addFourier(Image *inImg) {

	if (inImg->format == IMAGE_FORMAT_RGB565) {
		return NULL;
	}

	//should be editted. size lar gereginden falza koyuluyor right now.

	inImg->fourierx = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize += inImg->size * 8;
	inImg->fouriery = (SDRAM_BANK_ADDR + WRITE_READ_ADDR + allocatedSize);
	allocatedSize += inImg->size * 8;
	return;
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

/**
 * @brief Apply spatial filtering to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data after spatial filtering.
 * @param[in] filterType Type of filter (LINEAR_FILTER or NON_LINEAR_FILTER).
 * @param[in] filterName Name of the filter (LOW_PASS, HIGH_PASS, etc.).
 */

// Comparison function for qsort
int compare_uint8(const void *a, const void *b) {
	return (*(uint8_t*) a - *(uint8_t*) b);
}

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

// private

#include <limits.h>

unsigned int orderof2(const unsigned int value) {
	return (unsigned int) 1
			<< ((sizeof(unsigned int) * CHAR_BIT) - __builtin_clz(value) - 1);
}

// end of dft  -> add dft filtering precisely

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

void regionGrowing(Image *inImg, Image *outImg, int y, int x, int threshold) {

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



void IMAGE_basicSegmentation(const Image *inImg, Image *outImg,
		uint8_t targetValue) {
	for (int i = 0; i < inImg->size; ++i) {
		outImg->pixels[i] = (inImg->pixels[i] == targetValue) ? 255 : 0;
	}
}

/**
 * @defgroup MorphologicalOps Morphological Image Processing
 * @{
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
	Image *temp = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_erode(inImg, temp, kernelSize);
	DIP_dilate(temp, outImg, kernelSize);
}

//***************************13.5 Closing************************
void DIP_closing(const Image *inImg, Image *outImg, int kernelSize) {
	// Apply dilation followed by erosion
	Image *temp = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);
	DIP_dilate(inImg, temp, kernelSize);
	DIP_erode(temp, outImg, kernelSize);
}


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



// affine transforms

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


// Function to perform scaling
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

// Function to perform rotation
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

// Function to perform shear
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
