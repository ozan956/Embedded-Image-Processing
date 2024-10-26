#ifndef INC_DIP_SERIAL_HPP_
#define INC_DIP_SERIAL_HPP_

#include "Arduino.h"
#include "esp_camera.h"
#include "dip.h"

#define UART_BUF_SIZE 50000

void Capture(Image *img);
void Send(Image *img);
void delay_ms(uint32_t ms);

#endif /* INC_DIP_SERIAL_HPP_ */