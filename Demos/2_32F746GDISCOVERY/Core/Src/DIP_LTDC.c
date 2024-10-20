/*
 * DIP_LTDC.c
 *
 *  Created on: Apr 22, 2024
 *      Author: ozand
 */
#include "DIP_LTDC.h"

uint8_t LTDC_showImage(Image *test){

	switch(test->format){
	case IMAGE_FORMAT_RGB565: {
		HAL_LTDC_SetPixelFormat(&hltdc, LTDC_PIXEL_FORMAT_RGB565, 0);break;
	}
	case IMAGE_FORMAT_GRAYSCALE: {
		HAL_LTDC_SetPixelFormat(&hltdc, LTDC_PIXEL_FORMAT_L8, 0);break;
	}
	default:{
		// should not go there
		return -1;
	}
	}

	HAL_LTDC_SetAddress(&hltdc, (uint32_t) test->pixels,
	LTDC_LAYER_1);
	HAL_LTDC_Reload(&hltdc, LTDC_RELOAD_IMMEDIATE);

	return 0;
}
