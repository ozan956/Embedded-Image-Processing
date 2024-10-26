import PC_Image_Transfer_Library as serialImage

serialPort = 'COM3' 
Img = "cartest.png"
#Img = "panda.png"
#Img = "barbara.bmp"

serialImage.init(serialPort)


while 1:
    
    # press ESC to terminate
    req_type = serialImage.poll_transfer_request()  
    
    if (req_type == serialImage.IMAGE_READ):
        serialImage.send_image(Img)  
        
    elif (req_type == serialImage.IMAGE_WRITE):
        serialImage.read_image()  
        
 