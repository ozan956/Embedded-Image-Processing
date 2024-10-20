/*
 * DIP_UTIL.h
 *
 *  Created on: May 4, 2024
 *      Author: ozand
 */

#ifndef INC_DIP_UTIL_H_
#define INC_DIP_UTIL_H_
#include "DIP.h"

void sobelUtilizedMirrorPadding(const Image *inImg, int *pixelValue, int y,
		int x, int fx, int fy);

void DIP_UTIL_sobelFilter(const Image *inImg, Image *angle);

void mirrorPadding(Image *inImg, int i, int j, float *r, float *q,
		int condition);

void DIP_UTIL_nonMaxSupression(Image *inImg, Image *nonmax);

void DIP_UTIL_doubleTreshold(Image *inImg, float lowThresholdRatio,
		float highThresholdRatio, Image *outImg, int weak, int strong);

void DIP_UTIL_hysteresis(Image *inImg, Image *outImg, int weak, int strong);

#endif /* INC_DIP_UTIL_H_ */
