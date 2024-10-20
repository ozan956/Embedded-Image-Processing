/*
 * DIP.h
 *
 *  Created on: Mar 28, 2024
 *      Author: ozand
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <complex.h>
#include "arm_math.h"
#include "arm_const_structs.h"

#ifndef DIP_H
#define DIP_H

#ifdef CAMERA
#define CAMERA_FRAME_BUFFER               0xC0000000

#endif

#define SDRAM_BANK_ADDR                 ((uint32_t)0xC0000000)
#define WRITE_READ_ADDR     			((uint32_t)0x1000)


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


//error codes

#define IMG_CREATE_ERR -1
#define IMG_NO_ERR  0

#define MAX_INTENSITY  255
//#define DBL_MAX 1.7976931348623158e+308

// FOR SPATİAL FİLTERS
#define LINEAR_FILTER 0
#define NON_LINEAR_FILTER 1
#define LOW_PASS 0
#define HIGH_PASS 1
#define BAND_PASS 2
#define MEDIAN 3
#define MINIMUM 4
#define MAXIMUM 5
#define pi acos(-1)



// Define a structure to represent an image
typedef struct {
    uint32_t width;        // Width of the image in pixels
    uint32_t height;       // Height of the image in pixels
    uint8_t *pixels;       // Pointer to the pixel data
    uint32_t size;         // Size of the image data in bytes
    uint8_t format;      // Number of color channels (e.g., 3 for RGB, 4 for RGBA)
    float *fourierx;       // Pointer to the pixel data
    float *fouriery;       // Pointer to the pixel data

}Image;

typedef struct {
    double real;
    double imag;
} Complex;

// Structure to hold filter kernel
typedef struct {
    int size;
    float **data;
    uint8_t *data_8;
} FilterKernel;


#define INDEX_TO_OFFSET(x, y,WIDTH) ((y) * WIDTH + (x))


void DIP_negative(const Image *inImg, Image *outImg);

void DIP_powerTransform(const Image *inImg, Image *outImg, uint8_t gamma);

int DIP_histogramFormation(const Image *inImg, int *histogram);

void DIP_histogramEqualization(const Image *inImg, Image *outImg);

void DIP_histogramSpecification(const Image *inImg, Image *outImg, const int* targetHistogram);

void DIP_grayscaleWatermark(const Image *inImg, Image *outImg, const char* watermark, int x, int y);

void DIP_filter2D(const Image *inImg, Image *outImg, int size,
		float filter[size][size]);

void DIP_FFT(const Image *inImg);

void DIP_getMagnitudeSpec(const Image *inImg, Image *outImg, int fftshift);

void DIP_scale(const Image *inImg, Image *outImg, float scaleX, float scaleY);

void DIP_rotate(const Image *inImg, Image *outImg, float angle);

void DIP_shear(const Image *inImg, Image *outImg, float shx, float shy);

#endif //DIP_H
