import serial
import time

port = "COM3"  # Change this to your UART port (e.g., "/dev/ttyUSB0" for Linux)
baudrate = 2000000  # Adjust as needed

def send_command(command, ser):
    """Send a command byte over UART and receive a response."""
    ser.write(bytes([command]))  # Send the command as a single byte
    response = ser.read(1)  # Read a single byte response
    if response:
        print(f"Received response: {response.hex()}")
    else:
        print("No response received")

def main():

    
    try:
        ser = serial.Serial(port, baudrate, timeout=0.1)
        time.sleep(2)  # Allow some time for the connection to establish
        
        while True:
            user_input = input("Enter command (on: 0, off: 1, exit: q): ")
            
            if user_input.lower() == 'q':
                print("Exiting...")
                break
            elif user_input == '0':
                send_command(0x00, ser)
            elif user_input == '1':
                send_command(0x01, ser)
            else:
                print("Invalid input. Enter 0, 1, or q to quit.")
                
    except serial.SerialException as e:
        print(f"Serial error: {e}")
    finally:
        if 'ser' in locals() and ser.is_open:
            ser.close()
            print("Serial connection closed.")

if __name__ == "__main__":
    main()
