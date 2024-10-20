/*
 * DIP_SERIAL.c
 *
 *  Created on: Mar 28, 2024
 *      Author: ozand
 */
#include "DIP.h"
#include "DIP_SERIAL.h"

void SERIAL_imageCapture(Image *img) {

	uint8_t request_start_sequence[3] = "STR";

	uint16_t _blocksize = 65535, _lastblocksize = 0;
	uint32_t i = 0, _blockCount = 0;

	HAL_Delay(6000);
	uint16_t sizear[3] = { img->width, img->height, img->format };
	int size = 480 * 270;

	if (size < 65536)
		_blocksize = size;

	_blockCount = size / _blocksize;
	_lastblocksize = (uint16_t) (size % _blocksize);

	HAL_UART_Transmit(&huart1, request_start_sequence, 3, 100);
	HAL_Delay(1);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->width), sizeof(uint16_t), 100);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->height), sizeof(uint16_t),
			100);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->format), sizeof(uint16_t),
			100);

	for (i = 0; i < _blockCount; i++)
		HAL_UART_Receive(&huart1, img->pixels + (i * _blocksize), _blocksize,
				10000);

	if (_lastblocksize)
		HAL_UART_Receive(&huart1, img->pixels + (i * _blocksize),
				_lastblocksize, 10000);

	return;
}

void SERIAL_imageSend(Image *img) {
	uint8_t request_start_sequence[3] = "STW";

	uint16_t _blocksize = 65535, _lastblocksize = 0;
	uint32_t i = 0, _blockCount = 0;

	if (img->size < 65536)
		_blocksize = img->size;

	_blockCount = img->size / _blocksize;
	_lastblocksize = (uint16_t) (img->size % _blocksize);

	HAL_UART_Transmit(&huart1, request_start_sequence, 3, 100);
	HAL_Delay(1);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->width), sizeof(uint16_t), 100);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->height), sizeof(uint16_t),
			100);
	HAL_UART_Transmit(&huart1, (uint8_t*) (&img->format), sizeof(uint8_t),
			100);
	HAL_Delay(200);
	for (i = 0; i < _blockCount; i++)
		HAL_UART_Transmit(&huart1, img->pixels + (i * _blocksize), _blocksize,
				10000);

	if (_lastblocksize)
		HAL_UART_Transmit(&huart1, img->pixels + (i * _blocksize),
				_lastblocksize, 10000);
	HAL_Delay(200);
}

