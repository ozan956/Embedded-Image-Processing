/*
 * DIP_UTIL.c
 *
 *  Created on: May 4, 2024
 *      Author: ozand
 */
#include "DIP_UTIL.h"

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
