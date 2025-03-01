from PIL import Image
import numpy as np

def convert(input, output):
    # Open the image
    img = Image.open(input).convert("RGB")
    img = img.resize((480, 270))  # Ensure the image is of size 480x270
    
    # Convert to NumPy array
    img_data = np.array(img, dtype=np.uint8)
    
    # Convert to 16-bit RGB565 format
    def rgb_to_rgb565(r, g, b):
        return ((r & 0xF8) << 8) | ((g & 0xFC) << 3) | (b >> 3)
    
    img_16bit = [rgb_to_rgb565(r, g, b) for row in img_data for r, g, b in row]
    
    # Generate C header
    with open(output, "w") as f:
        f.write("#include <stdint.h>\n")
        f.write(f"const uint16_t image_data[] = {{\n")
        
        for i, pixel in enumerate(img_16bit):
            f.write(f"  0x{pixel:04X},")
            if (i + 1) % 12 == 0:
                f.write("\n")
        
        f.write("};\n")

# Example usage
convert("panda.png", "img_data.h")
