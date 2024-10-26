#include "dip_serial.h"
#include "dsp_platform.h"
#include "esp_log.h"
#include "dip.h"
#include "esp_dsp.h"
#include <math.h>
#include "DIP_Driver.h"
#include "dsps_fft2r.h"

const int PIN_BUTTON = 15;   
const int RED_LED =  21;
const int WHT_LED =  22;
int button_state = 0;

Image *img = nullptr;
#define TWO_PI 6.28318530718

Image *outimg = nullptr;
Image *fourier = nullptr;
Image *resized = nullptr;
DIP dip;

const int N_SAMPLES = 32;  // Number of samples for each dimension (must be a power of 2)
//const float PI = 3.14159265359;

SET_LOOP_TASK_STACK_SIZE(16 * 1024* 2);

void setup() {

  pinMode(PIN_BUTTON, INPUT_PULLUP);
  pinMode(RED_LED, OUTPUT);
  pinMode(WHT_LED, OUTPUT);
 
  initDIPManager(&dip);
  //Setup pins
  img = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
  outimg = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
  fourier = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
  resized = createImage(IMAGE_RES_480x272, IMAGE_FORMAT_GRAYSCALE);
  for(int i=0;i<320*240;i++){
    if(i > 320*120)
      img->pixels[i] = (uint8_t)0x00;
    else
      img->pixels[i] = (uint8_t)0xFF;
  }

  Serial.begin(250000);
    //Serial.printf("Arduino Stack was set to %d bytes", getArduinoLoopTaskStackSize());

    // Print unused stack for the task that is running setup()
   // Serial.printf("\nSetup() - Free Stack Space: %d", uxTaskGetStackHighWaterMark(NULL));


}

void loop() {
    //Serial.printf("\nSetup() - Free Stack Space: %d", uxTaskGetStackHighWaterMark(NULL));
    
  button_state = digitalRead(PIN_BUTTON);
  delay(100);
    
  if (button_state == LOW) {
    digitalWrite(RED_LED, HIGH);
    Capture(img);
    digitalWrite(RED_LED, LOW);

    delay(1000);

    digitalWrite(WHT_LED, HIGH);
    Send(img);
    digitalWrite(WHT_LED, LOW); 
/*
    DIP_negative(img,outimg);

    delay(1000);
    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 
    delay(10);


    DIP_powerTransform(img,outimg,3);

    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

    DIP_histogramEqualization(img, outimg);

    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

*/

//DENENECEK
/*
	 int istenilen_hist[MAX_INTENSITY + 1] = { 0 };
	 for (int i = 0; i < test_sdram->size; ++i) {
	 istenilen_hist[test_sdram->pixels[i]]++;
	 }

	 DIP_histogramSpecification(test, histogram_spec, istenilen_hist);

*/

/*
	char watermark[] = "test_data";
	DIP_grayscaleWatermark(img, outimg, watermark, 100, 100);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 
  

	float filter[] = {1.0/9,1.0/9,1.0/9,
	1.0/9,1.0/9,1.0/9,
	1.0/9,1.0/9,1.0/9};

	DIP_filter2D(img, outimg,3,filter);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_medianBlur(img, outimg, 5);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_grayscaleThreshold(img, outimg, 120);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_otsuMethod(img, outimg);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 


  DIP_sobelFilter(img, outimg);


  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

	regionGrowing(img, outimg, 97,199,20);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 


  IMAGE_basicSegmentation(img, outimg, 100);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

	DIP_dilate(img, outimg, 5);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_erode(img, outimg, 5);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);


	DIP_opening(img, outimg, 5);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_closing(img, outimg, 5);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);

  DIP_canny(img, outimg, 30, 70);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);


	translate(img, outimg, 10, 10);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);

	float scaleX = 1.2f;
	float scaleY = 1.2f;
	DIP_scale(img, outimg, scaleX, scaleY);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);

	DIP_rotate(img, outimg, 45);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);

	DIP_shear(img, outimg, 0.2, 0.1);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);



  int FFT_SIZE = 256;
  DIP_Resize(img,outimg,FFT_SIZE);
  DIP_fourier(img, outimg, FFT_SIZE);
  DIP_abs(outimg, outimg, FFT_SIZE, 1);
  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 
  DIP_fourier_inv(outimg, fourier, FFT_SIZE);
  DIP_abs(outimg, outimg, FFT_SIZE, 0);
  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 


    // Define Gaussian kernels (3x3 Gaussian blur in 1D)
    float kernelX[3] = {-1, 0, 1}; // Horizontal edge detection kernel
    float kernelY[3] = {1, 2, 1};  // Vertical edge detection kernel

    // Image dimensions (example size, ensure it matches the actual size of inImg)
    int width = 128;
    int height = 128;

    // Apply the separable 2D filter
    DIP_sepFilter2D(img, outimg, 3, kernelX, 3, kernelY, 0.0);
 
    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

    // Apply the separable 2D filter
    DIP_sqrBoxFilter(img, outimg, 3,1);
 
    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

    // Set Laplacian filter parameters
    int ddepth = -1;     // Depth of the output image, typically -1 to use same as input
    int ksize = 3;       // Size of the Laplacian kernel (must be positive and odd)
    double scale = 1.0;  // Scaling factor for the computed Laplacian values
    double delta = 0.0;  // Optional offset to add to the results
    int borderType = BORDER_DEFAULT;  // Border handling method

    // Apply the Laplacian filter
    DIP_Laplacian(img, outimg, ddepth, ksize, scale, delta, borderType);

    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

    // Define the desired output size (half of the original dimensions)
    Size dstSize = {img->width /2, img->height /2};


    DIP_pyrDown(img, outimg, &dstSize, BORDER_DEFAULT);

    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 

    dstSize.width *=2;
    dstSize.height *=2;

    // Perform downsampling with Gaussian blur
    DIP_pyrUp(outimg, resized, &dstSize, BORDER_DEFAULT, 3,1);

    digitalWrite(WHT_LED, HIGH);
    Send(resized);
    digitalWrite(WHT_LED, LOW); 

    // Set parameters for adaptive thresholding
    double maxValue = 255.0;
    int adaptiveMethod = 0;          // Currently, only mean is implemented
    int thresholdType = THRESH_BINARY; // Use binary thresholding
    int blockSize = 7;               // Block size (should be odd)
    double C = 10.0;                 // Constant to subtract from mean

    // Apply adaptive thresholding
    DIP_adaptiveThreshold(img, outimg, maxValue, adaptiveMethod, thresholdType, blockSize, C);

    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 



  DIP_shear(img, outimg, 0.2, 0.3);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 

  DIP_inverse_shear(outimg, resized, 0.2, 0.3);

  digitalWrite(WHT_LED, HIGH);
  Send(resized);
  digitalWrite(WHT_LED, LOW); 



    // Step 3: Define a perspective transformation matrix
    // This matrix simulates a slight tilt on the document
float M[3][3] = {
    {1.2, 0.2, -50},
    {0.1, 1.1, -40},
    {0.0004, 0.0006, 1.0}
};

    // Step 4: Apply perspective warp transformation
    DIP_warpPerspective(img, outimg, M);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW); 


float h[3][3] = {
    {1.0, 0.3, 0.0},
    {0.3, 1.0, 0.0},
    {0.001, 0.002, 1.0}
};

    // Apply perspective transformation
    DIP_perspectiveTransform(img, outimg, h);


    digitalWrite(WHT_LED, HIGH);
    Send(outimg);
    digitalWrite(WHT_LED, LOW); 



  resized->width = 80;
  resized->height = 80;
  resized->size = 80 * 80;
  resized->pixels = (uint8_t *)ps_malloc(80 * 80 * sizeof(uint8_t));

  for(int i=0;i<80;i++){
    for(int y=0;y<80;y++){
      resized->pixels[i] = img->pixels[(i+50) * img->width + (y+20)];
    }
  }


  for(int i=0;i<img->height;i++){
    for(int y=0;y<img->width;y++){
      outimg->pixels[i] = img->pixels[(i) * img->width + (y)] - 50;
    }
  }

  DIP_TemplateMatch(img, resized, outimg);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);
*/

  dip.canny(img,outimg,20,60);

  digitalWrite(WHT_LED, HIGH);
  Send(outimg);
  digitalWrite(WHT_LED, LOW);

  } else {
    digitalWrite(RED_LED, LOW);
    digitalWrite(WHT_LED, LOW);
  }

}