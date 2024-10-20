/*
 * DIP_SERIAL.h
 *
 *  Created on: Mar 28, 2024
 *      Author: ozand
 */
#include "stm32f7xx_hal.h"
#include "DIP.h"

#ifndef INC_DIP_SERIAL_H_
#define INC_DIP_SERIAL_H_

/**********************************/
// Peripherals and Pins
extern UART_HandleTypeDef huart1;

/**********************************/
// Functions
void SERIAL_imageCapture(Image *img);
void SERIAL_imageSend(Image *img);

#endif /* INC_DIP_SERIAL_H_ */
