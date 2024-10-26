import serial
import cv2
import numpy as np
import time
import sys 

IMAGE_FORMAT_RGB565		= 1
IMAGE_FORMAT_YUV422		= 2
IMAGE_FORMAT_GRAYSCALE	= 3

IMAGE_READ = 82
IMAGE_WRITE = 87
ser = ""

width = 100
height = 100
transferType = 0
imageFormat = 1

img565 = np.zeros(1)
im = np.zeros(1)


nullWindow = np.zeros((width,height))
cv2.imshow('',nullWindow)


formatArray = {1:"RGB565", 2:"YUV422", 3:"GRAYSCALE"}
requestArray = {82: "IMAGE_READ", 87: "IMAGE_WRITE"}
 

def init(serialPort):   
    global ser 
    ser = serial.Serial(serialPort, 250000, timeout=None)
    
    ser.flush()
    print(ser.name, "Opened")
    print("")


def terminate():
    global ser
    cv2.destroyAllWindows()
    ser.close()
    print('')
    sys.exit("ESC pressed")

def serial_read():
    return (int.from_bytes(ser.read(1),byteorder='little'))

def serial_write(b):
    global ser
    if(cv2.waitKey(10) == 27):
        terminate()
    out = int(b).to_bytes(1, byteorder='little')                
    ser.write(out)  
        

def poll_transfer_request():
    global width
    global height
    global imageFormat
    global ser
    
    
    while 1:
        if(cv2.waitKey(1) == 27):
            terminate() 
        
        # read start bytes
        if(serial_read() == 0x53):
            if(serial_read() == 0x54):
            
                # read image detail bytes
                
                transferType = serial_read()
                
                width = serial_read()
                width = (serial_read()<<8 & 0xFF00) | width 
                
                height = serial_read()
                height = (serial_read()<<8 & 0xFF00) | height 
                
                imageFormat = serial_read()
                print(imageFormat)
                print("request ",  requestArray[transferType])
                print("width ", width)
                print("height ", height)
                print("format ", formatArray[imageFormat])
                
                time.sleep(0.1)
                
                return transferType

def send_image(sendImg):  
    global ser 
    global width
    global height
    global imageFormat
    img = cv2.imread(sendImg)
          
    im = cv2.resize(img,(width,height))
    
    if imageFormat == IMAGE_FORMAT_RGB565:  
        
        # convert image to BGR565 format
        img565 = cv2.cvtColor(im,cv2.COLOR_BGR2BGR565)
        img565 = img565.astype(np.uint8)        
        
        # Byte swap -> Opencv and ili9341 are not compatiable in byte order.
        imgh = img565[:,:,0].copy()
        img565[:,:,0] = img565[:,:,1].copy()
        img565[:,:,1] = imgh.copy() 
        
        # 2D to 1D conversion
        img565 = np.reshape(img565, (width*height*2))  
        
        # send image
        ser.write(img565)
        
    elif imageFormat == IMAGE_FORMAT_GRAYSCALE:
        # convert image to Grayscale format
        imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
        
        # 2D to 1D conversion
        imgray = np.reshape(imgray, (width*height))    
        print(len(imgray))
        # send image
        ser.write(imgray)
        
    elif imageFormat == IMAGE_FORMAT_YUV422:
        
        # convert image to YUV format
        img_yuv = cv2.cvtColor(im, cv2.COLOR_BGR2YUV)
        
        # Generate YUV422 from OpenCV YUV format
        y0 = np.expand_dims(img_yuv[...,0][::,::2], axis=2)
        u = np.expand_dims(img_yuv[...,1][::,::2], axis=2)
        y1 = np.expand_dims(img_yuv[...,0][::,1::2], axis=2)
        v = np.expand_dims(img_yuv[...,2][::,::2], axis=2)
        img422 = np.concatenate((y0, u, y1, v), axis=2)
        img422 = img422.reshape(img422.shape[0], img422.shape[1] * 2, int(img422.shape[2] / 2))
        
        # 2D to 1D conversion
        img422 = np.reshape(img422, (width*height*2))  
        
        # send image
        ser.write(img422)
        
            
    print("ready")
    print("")

def read_image():  
    global ser 
    global width
    global height
    global imageFormat
    global img565
    ser.timeout = 10
    
    timestamp = time.strftime('%Y_%m_%d_%H%M%S', time.localtime())  
    
    if imageFormat == IMAGE_FORMAT_RGB565:
        # read RGB565 image from serial port
        img565 = ser.read(width*height*2)
        img565 = list(img565) 
        
        # 1D to 2D conversion
        img565 = np.reshape(img565, (height, width, 2))
        img565 = img565.astype(np.uint8)
        
        # Byte swap -> Opencv and OV7670 are not compatiable in byte order.
        imgh = img565[:,:,0].copy()
        img565[:,:,0] = img565[:,:,1].copy()
        img565[:,:,1] = imgh.copy() 
        
        # convert BGR565 image to BGR fomrat
        im = cv2.cvtColor(img565,cv2.COLOR_BGR5652BGR)       
        
        # show and save image
        cv2.imshow("img565", im)        
        cv2.imwrite("C:/Users/OutImagesESP/"+ timestamp +"_rgb.png", im)        
        
    
    elif imageFormat == IMAGE_FORMAT_GRAYSCALE:   
        # read Grayscale image from serial port 
        
        imgray =  ser.read(width*height)
        print(len(imgray))
        imgray = list(imgray)      

        # 1D to 2D conversion
        imgray = np.reshape(imgray, (height, width))
        im = imgray.astype(np.uint8)          
        # show and save image
        cv2.imwrite("C:/capture/" + timestamp + "_gray.png", im)
        cv2.imshow("imgray", im)
        
    elif imageFormat == IMAGE_FORMAT_YUV422:   
        # read YUV422 image from serial port 
        img422 = ser.read(width*height*2)
        img422 = list(img422)
        
        # 1D to 2D conversion
        img422 = np.reshape(img422, (height, width, 2))
        img422 = img422.astype(np.uint8)
                
        # convert YUV422 image to BGR fomrat
        im = cv2.cvtColor(img422,cv2.COLOR_YUV2BGR_YVYU)     
        
        # show and save image
        cv2.imshow("yuv", im)              
        cv2.imwrite("C:/Users/OutImagesESP/"+ timestamp +"_yuv.png", im)        
        
    ser.timeout = 0.1
    print("ready")  
    print("")  