/*
 * DIP_UTIL.c
 *
 *  Created on: May 4, 2024
 *      Author: ozand
 */
#include "DIP_UTIL.h"

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

/**
 * @brief Applies Sobel filter to the input image and calculates the gradient angle.
 *
 * This function computes the Sobel filter on the input image and outputs the
 * gradient angle in the provided angle image.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] angle Pointer to the output image structure where gradient angles will be stored.
 */
void DIP_UTIL_sobelFilter(const Image *inImg, Image *angle) {
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
			angle->fouriery[y * inImg->width + x] = sqrt(
					(float) (gradientX * gradientX + gradientY * gradientY));

			gradMax =
					(angle->fouriery[y * inImg->width + x] > gradMax) ?
							angle->fouriery[y * inImg->width + x] : gradMax;

			angle->fourierx[y * inImg->width + x] = (float) atan2(gradientY,
					gradientX);
		}
	}

	for (int i = 0; i < angle->size; i++) {
		angle->pixels[i] = (uint8_t) ((angle->fouriery[i] / gradMax) * 255);
	}

}

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

/**
 * @brief Performs non-maximum suppression on the input image.
 *
 * This function applies non-maximum suppression to the input image to
 * thin out the detected edges. The result is stored in the nonmax image.
 *
 * @param[in] inImg Pointer to the input image structure.
 * @param[out] nonmax Pointer to the output image structure for storing non-max suppressed results.
 */
void DIP_UTIL_nonMaxSupression(Image *inImg, Image *nonmax) {
	int i, j;

	uint8_t *G = inImg->pixels;
	float *theta = inImg->fourierx;

	addFourier(nonmax);

	for (int i = 0; i < inImg->size; i++)
		nonmax->fourierx[i] = 0x0000;

	Image *angles = (Image*) createImage(IMAGE_RES_480x272,
	IMAGE_FORMAT_GRAYSCALE);

	addFourier(angles);

	// Convert angles to degrees
	for (i = 0; i < inImg->height; i++) {
		for (j = 0; j < inImg->width; j++) {
			angles->fourierx[i * inImg->width + j] = theta[i * inImg->width + j]
					* 180.0 / M_PI;
			if (angles->fourierx[i * inImg->width + j] < 0)
				angles->fourierx[i * inImg->width + j] += 180.0;
		}
	}

	// Non-maximum suppression
	for (i = 0; i < inImg->height; i++) {
		for (j = 0; j < inImg->width; j++) {
			float q, r;
			q = 255.0;
			r = 255.0;
			if ((0.0 <= angles->fourierx[i * inImg->width + j]
					&& angles->fourierx[i * inImg->width + j] < 22.5)
					|| (157.5 <= angles->fourierx[i * inImg->width + j]
							&& angles->fourierx[i * inImg->width + j] <= 180.0)) {
				mirrorPadding(inImg, i, j, &r, &q, 0);
			} else if (22.5 <= angles->fourierx[i * inImg->width + j]
					&& angles->fourierx[i * inImg->width + j] < 67.5) {
				mirrorPadding(inImg, i, j, &r, &q, 1);
			} else if (67.5 <= angles->fourierx[i * inImg->width + j]
					&& angles->fourierx[i * inImg->width + j] < 112.5) {
				mirrorPadding(inImg, i, j, &r, &q, 2);
			} else if (112.5 <= angles->fourierx[i * inImg->width + j]
					&& angles->fourierx[i * inImg->width + j] < 157.5) {
				mirrorPadding(inImg, i, j, &r, &q, 3);
			}

			if ((G[i * inImg->width + j] >= q)
					&& (G[i * inImg->width + j] >= r)) {
				nonmax->pixels[i * inImg->width + j] = G[i * inImg->width + j];
			} else {
				nonmax->pixels[i * inImg->width + j] = (uint8_t) 0;
			}
		}
	}

}

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
int compare_uint8(const void *a, const void *b) {
	return (*(uint8_t*) a - *(uint8_t*) b);
}


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
unsigned int orderof2(const unsigned int value) {
	return (unsigned int) 1
			<< ((sizeof(unsigned int) * CHAR_BIT) - __builtin_clz(value) - 1);
}


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
    float **filter = (float **)malloc(ksize * sizeof(float *));
    for (int i = 0; i < ksize; i++) {
        filter[i] = (float *)malloc(ksize * sizeof(float));
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
