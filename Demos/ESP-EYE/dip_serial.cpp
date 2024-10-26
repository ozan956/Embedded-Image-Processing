#include "dip_serial.h"
const int RED_LED =  21;
const int WHT_LED =  22;

void Capture(Image *img)
{
  int readbyte = 0;
    //Serial.flush();
    //Serial.begin(115200);
    uint32_t startms = millis();  
    uint8_t request_start_sequence[4] = "STR";

    uint16_t _blocksize = UART_BUF_SIZE, _lastblocksize = 0;
    uint32_t i = 0, _blockCount = 0, len = img->height*img->width;

    if (len < UART_BUF_SIZE)
    {
        _blocksize = len;
    }

    _blockCount = len / _blocksize;
    _lastblocksize = (uint16_t)(len % _blocksize);

    Serial.write(request_start_sequence, 3);
    delay_ms(1);

    uint16_t w = img->width;
    uint16_t h = img->height;
    uint16_t f = img->format;
    
    Serial.write((uint8_t *)&w, sizeof(w));Serial.flush();
    Serial.write((uint8_t *)&h, sizeof(h));Serial.flush();
    Serial.write((uint8_t *)&f, sizeof(f));

    // capture chunks of data with the legnth of UART_BUF_SIZE
    char test;
    Serial.readBytes(&test, 1);
    for (i = 0; i < _blockCount; i++)
    {
        uint32_t j = 0;
        while (j < _blocksize)
        {
            Serial.readBytes(img->pixels + (i * _blocksize) + j, 1);
            j++;readbyte++;
        }
    }
    // capture remainder data
    if (_lastblocksize)
    {
        uint32_t j = 0;
        while (j < _lastblocksize)
        {
            Serial.readBytes(img->pixels + (i * _blocksize) + j, 1);
            j++;readbyte++;
        }
    }

    if(readbyte == 76800){
    digitalWrite(WHT_LED, HIGH);
    delay_ms(400);
    digitalWrite(WHT_LED, LOW);
    delay_ms(400);
    digitalWrite(WHT_LED, HIGH);
    delay_ms(400);
    digitalWrite(WHT_LED, LOW);
    delay_ms(400);
    }
    
    delay_ms(200);

}

void Send(Image *img)
{
    //Serial.flush();
    //Serial.begin(2000000);
    uint8_t request_start_sequence[4] = "STW";
    uint8_t a = 0;

    uint16_t _blocksize = UART_BUF_SIZE, _lastblocksize = 0;
    uint32_t i = 0, _blockCount = 0, len = img->width*img->height;

    if (img->size < UART_BUF_SIZE)
        _blocksize = len;

    _blockCount = len / _blocksize;
    _lastblocksize = (uint16_t)(len % _blocksize);

    Serial.write(request_start_sequence, 3);
    delay_ms(1);


    uint16_t w = img->width;
    uint16_t h = img->height;
    uint16_t f = img->format;
    
    Serial.write((uint8_t *)&w, sizeof(w));Serial.flush();
    Serial.write((uint8_t *)&h, sizeof(h));Serial.flush();
    Serial.write((uint8_t *)&f, sizeof(f));Serial.flush();
    delay_ms(100);

    // send chunks of data with the legnth of UART_BUF_SIZE
    for (i = 0; i < _blockCount; i++)
        Serial.write(img->pixels + (i * _blocksize), _blocksize);Serial.flush();

    // send remainder data
    if (_lastblocksize)
        Serial.write(img->pixels + (i * _blocksize), _lastblocksize);Serial.flush();

    delay_ms(200);
    //Serial.flush();
    //Serial.begin(115200);
}





void delay_ms(uint32_t ms){
  uint32_t startms = millis();  
  while(millis() < (startms + ms)){
    delay(1);
  }
}