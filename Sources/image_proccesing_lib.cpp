#include <image_proccesing_lib.hpp>

using namespace std;
using namespace cv;

/**
 * @defgroup IntensityTransformations Intensity Transformations
 * @{
 */

 /**
  * @brief Negate the pixel intensities of an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data with negated intensities.
  */
void IMAGE_Negative(const uint8_t* imgData, uint8_t* outData) {
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = MAX_INTENSITY - imgData[i];
    }
}


/**
 * @brief Apply power law transformation to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with power-law transformed intensities.
 * @param[in] gamma The power-law transformation parameter.
 */
void IMAGE_powerTransform(const uint8_t* imgData, uint8_t* outData, uint8_t gamma) {
    // Calculate the scaling factor 'c' to keep pixel values in the valid intensity range
    float c = MAX_INTENSITY / pow(MAX_INTENSITY, gamma);

    // Perform the power transform
    for (int i = 0; i < IMAGE_SIZE; ++i) {

        outData[i] = (uint8_t)(c * pow(imgData[i], gamma));
    }
}

/** @} */ // End of IntensityTransformations group

/**
 * @defgroup HistogramOperations Histogram Operations
 * @{
 */

 /**
  * @brief Calculate the histogram of an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] histogram Array to store the histogram.
  * @return Total number of pixels in the image.
  */
int IMAGE_histogramFormation(const uint8_t* imgData, int* histogram) {
    // Initialize histogram array with zeros
    fill(histogram, histogram + MAX_INTENSITY + 1, 0);

    // Calculate histogram
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        histogram[imgData[i]]++;
    }

    // Return the total number of pixels in the image
    return IMAGE_SIZE;
}

/**
 * @brief Perform histogram equalization on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with equalized intensities.
 */
void IMAGE_histogramEqualization(uint8_t* imgData, uint8_t* outData) {
    // Calculate histogram
    int histogram[MAX_INTENSITY + 1] = { 0 };
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        histogram[imgData[i]]++;
    }

    // Calculate cumulative distribution function
    int cumulativeHistogram[MAX_INTENSITY + 1] = { 0 };
    cumulativeHistogram[0] = histogram[0];
    for (int i = 1; i <= MAX_INTENSITY; ++i) {
        cumulativeHistogram[i] = cumulativeHistogram[i - 1] + histogram[i];
    }

    // Calculate equalization mapping
    float scale = static_cast<float>(MAX_INTENSITY) / IMAGE_SIZE;
    uint8_t equalizationMap[MAX_INTENSITY + 1];
    for (int i = 0; i <= MAX_INTENSITY; ++i) {
        equalizationMap[i] = static_cast<uint8_t>(round(cumulativeHistogram[i] * scale));
    }

    // Perform histogram equalization on the image
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = equalizationMap[imgData[i]];
    }
}

/**
 * @brief Perform histogram specification on an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output image data with specified histogram.
 * @param[in] targetHistogram The target histogram to match.
 */
void IMAGE_histogramSpecification(const uint8_t* imgData, uint8_t* outData, const int* targetHistogram) {
    int sourceHistogram[MAX_INTENSITY + 1] = { 0 };
    IMAGE_histogramFormation(imgData, sourceHistogram);

    // Calculate cumulative distribution function for source and target histograms
    int sourceCumulativeHistogram[MAX_INTENSITY + 1] = { 0 };
    int targetCumulativeHistogram[MAX_INTENSITY + 1] = { 0 };

    sourceCumulativeHistogram[0] = sourceHistogram[0];
    targetCumulativeHistogram[0] = targetHistogram[0];

    for (int i = 1; i <= MAX_INTENSITY; ++i) {
        sourceCumulativeHistogram[i] = sourceCumulativeHistogram[i - 1] + sourceHistogram[i];
        targetCumulativeHistogram[i] = targetCumulativeHistogram[i - 1] + targetHistogram[i];
    }

    // Calculate equalization mapping
    float scale = static_cast<float>(MAX_INTENSITY) / IMAGE_SIZE;
    uint8_t equalizationMap[MAX_INTENSITY + 1];
    for (int i = 0; i <= MAX_INTENSITY; ++i) {
        int nearestMatch = 0;
        int minDiff = abs(sourceCumulativeHistogram[i] - targetCumulativeHistogram[0]);

        for (int j = 1; j <= MAX_INTENSITY; ++j) {
            int diff = abs(sourceCumulativeHistogram[i] - targetCumulativeHistogram[j]);
            if (diff < minDiff) {
                minDiff = diff;
                nearestMatch = j;
            }
        }

        equalizationMap[i] = nearestMatch;
    }

    // Perform histogram specification on the image
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = equalizationMap[imgData[i]];
    }
}

/** @} */ // End of HistogramOperations group

/**
 * @defgroup Watermarking Digital Watermarking
 * @{
 */

 /**
  * @brief Embed a grayscale watermark into an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data with embedded watermark.
  * @param[in] watermark The watermark to embed.
  * @param[in] x X-coordinate of the watermark insertion point.
  * @param[in] y Y-coordinate of the watermark insertion point.
  */
void IMAGE_grayscaleWatermark(vector<uint8_t>& imgData, vector<uint8_t>& outData, const string& watermark, int x, int y) {
    // Make sure the watermark fits within the image dimensions
    int textWidth = watermark.length();
    if (x < 0 || x + textWidth >= WIDTH || y < 0 || y >= HEIGHT) {
        cout << "Watermark does not fit within the image boundaries." << endl;
        return;
    }

    // Convert the watermark to grayscale
    vector<uint8_t> watermarkData(textWidth, 255); // Set all watermark pixels to white (255)

    // Insert the watermark into the image
    for (size_t i = 0; i < watermarkData.size(); ++i) {
        outData[y * WIDTH + x + i] = watermarkData[i];
    }
}

/** @} */ // End of Watermarking group

/**
 * @defgroup SpatialFilters Spatial Filters
 * @{
 */

 /**
  * @brief Apply spatial filtering to an image.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output image data after spatial filtering.
  * @param[in] filterType Type of filter (LINEAR_FILTER or NON_LINEAR_FILTER).
  * @param[in] filterName Name of the filter (LOW_PASS, HIGH_PASS, etc.).
  */

void IMAGE_spatialFilter(const uint8_t* imgData, uint8_t* outData, int filterType, int filterName) {

    int filterSize = 3;

    int lowPassFilter[3][3] = {
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1}
    };

    int highPassFilter[3][3] = {
        {0, -1, 0},
        {-1, 8, -1},
        {0, -1, 0}
    };

    int bandPassFilter[3][3] = {
        {-1, -1, -1},
        {-1,  8, -1},
        {-1, -1, -1}
    };

    if (filterType == LINEAR_FILTER) {
        for (int y = 1; y < HEIGHT - 1; ++y) {
            for (int x = 1; x < WIDTH - 1; ++x) {
                int sum = 0;

                for (int fy = 0; fy < 3; ++fy) {
                    for (int fx = 0; fx < 3; ++fx) {
                        int pixelValue = imgData[(y + fy - 1) * WIDTH + (x + fx - 1)];

                        if (filterName == 0) { // Low-pass filter
                            sum += pixelValue * lowPassFilter[fy][fx];
                        }
                        else if (filterName == 1) { // High-pass filter
                            sum -= pixelValue * highPassFilter[fy][fx];
                        }
                        else if (filterName == 2) { // Band-pass filter
                            if (fy == 1 && fx == 1) {
                                sum += pixelValue * bandPassFilter[fy][fx];
                            }
                            else {
                                sum -= pixelValue * bandPassFilter[fy][fx];
                            }
                        }
                    }
                }

                outData[y * WIDTH + x] = static_cast<uint8_t>(sum / 9); // Normalize by dividing by 9 for a 3x3 filter
            }
        }
    }
    else if (filterType == NON_LINEAR_FILTER) {

        int filterRadius = filterSize / 2;
        vector<uint8_t> filterValues(filterSize * filterSize);

        for (int y = filterRadius; y < HEIGHT - filterRadius; ++y) {
            for (int x = filterRadius; x < WIDTH - filterRadius; ++x) {
                int index = INDEX_TO_OFFSET(x, y);

                int filterIndex = 0;
                for (int fy = -filterRadius; fy <= filterRadius; ++fy) {
                    for (int fx = -filterRadius; fx <= filterRadius; ++fx) {
                        filterValues[filterIndex] = imgData[INDEX_TO_OFFSET(x + fx, y + fy)];
                        filterIndex++;
                    }
                }

                if (filterName == MEDIAN) {
                    // Apply median filtering
                    sort(filterValues.begin(), filterValues.end());
                    outData[index] = filterValues[filterSize * filterSize / 2];
                }
                else if (filterName == MINIMUM) {
                    // Apply order statistics filtering (e.g., maximum or minimum)
                    sort(filterValues.begin(), filterValues.end());
                    // You can adjust this index based on the desired order statistic
                    outData[index] = filterValues[0]; // Minimum value
                }
                else if (filterName == MAXIMUM) {
                    // Apply order statistics filtering (e.g., maximum or minimum)
                    sort(filterValues.begin(), filterValues.end());
                    // You can adjust this index based on the desired order statistic
                    outData[index] = filterValues[filterSize * filterSize - 1]; // Minimum value
                }
            }
        }
    }

}

/** @} */ // End of SpatialFilters group

/**
 * @defgroup DFT 2D Discrete Fourier Transform (DFT)
 * @{
 */

void manualNormalize(cv::Mat& src, cv::Mat& dst) {
    double minVal = DBL_MAX;
    double maxVal = -DBL_MAX;
    int srcType = src.type();

    for (int i = 0; i < src.rows; i++) {
        for (int j = 0; j < src.cols; j++) {
            double value;
            if (srcType == CV_8U) {
                value = static_cast<double>(src.at<uchar>(i, j));
            }
            else if (srcType == CV_32F) {
                value = static_cast<double>(src.at<float>(i, j));
            }
            else if (srcType == CV_64F) {
                value = src.at<double>(i, j);
            }
            else {
                std::cerr << "Unsupported data type" << std::endl;
                return;
            }

            if (value < minVal) minVal = value;
            if (value > maxVal) maxVal = value;
        }
    }

    for (int i = 0; i < src.rows; i++) {
        for (int j = 0; j < src.cols; j++) {
            double value;
            if (srcType == CV_8U) {
                value = static_cast<double>(src.at<uchar>(i, j));
            }
            else if (srcType == CV_32F) {
                value = static_cast<double>(src.at<float>(i, j));
            }
            else if (srcType == CV_64F) {
                value = src.at<double>(i, j);
            }

            double normalizedValue = (value - minVal) / (maxVal - minVal);
            dst.at<double>(i, j) = normalizedValue;
        }
    }
}

void fft(std::vector<std::complex<double>>& data) {
    int n = data.size();
    if (n <= 1) {
        return;
    }

    std::vector<std::complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int i = 0; i < n / 2; ++i) {
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * i / n) * odd[i];
        data[i] = even[i] + t;
        data[i + n / 2] = even[i] - t;
    }
}

void fft2D(std::vector<std::vector<std::complex<double>>>& image) {
    int rows = image.size();
    int cols = image[0].size();

    for (int i = 0; i < rows; ++i) {
        fft(image[i]);
    }

    for (int j = 0; j < cols; ++j) {
        std::vector<std::complex<double>> column(rows);
        for (int i = 0; i < rows; ++i) {
            column[i] = image[i][j];
        }
        fft(column);
        for (int i = 0; i < rows; ++i) {
            image[i][j] = column[i];
        }
    }
}

void ifft2D(std::vector<std::vector<std::complex<double>>>& spectrum, cv::Mat& reconstructedImage) {
    int rows = spectrum.size();
    int cols = spectrum[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            spectrum[i][j] = std::conj(spectrum[i][j]);
        }
    }

    for (int i = 0; i < rows; ++i) {
        fft(spectrum[i]);
    }

    for (int j = 0; j < cols; ++j) {
        std::vector<std::complex<double>> column(rows);
        for (int i = 0; i < rows; ++i) {
            column[i] = spectrum[i][j];
        }
        fft(column);
        for (int i = 0; i < rows; ++i) {
            spectrum[i][j] = column[i];
        }
    }

    // Calculate the magnitude of the IFT result
    reconstructedImage.create(rows, cols, CV_64F);
    double maxMagnitude = -DBL_MAX;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double magnitude = std::abs(spectrum[i][j]);
            if (magnitude > maxMagnitude) {
                maxMagnitude = magnitude;
            }
        }
    }

    // Normalize the IFT result
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double magnitude = std::abs(spectrum[i][j]);
            reconstructedImage.at<double>(i, j) = magnitude / maxMagnitude;
        }
    }
}

void IMAGE_idealBPFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double minFrequency, double maxFrequency) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double wavelengthTopLeft = sqrt((i * i) + (j * j)); // Wavelength from top-left corner
            double wavelengthTopRight = sqrt((i * i) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from top-right corner
            double wavelengthBottomLeft = sqrt(((rows - i - 1) * (rows - i - 1)) + (j * j)); // Wavelength from bottom-left corner
            double wavelengthBottomRight = sqrt(((rows - i - 1) * (rows - i - 1)) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from bottom-right corner

            // Calculate frequency as the reciprocal of wavelength
            double frequencyTopLeft = 1.0 / wavelengthTopLeft;
            double frequencyTopRight = 1.0 / wavelengthTopRight;
            double frequencyBottomLeft = 1.0 / wavelengthBottomLeft;
            double frequencyBottomRight = 1.0 / wavelengthBottomRight;

            // Check if the frequency is within the desired range
            if ((frequencyTopLeft >= minFrequency && frequencyTopLeft <= maxFrequency) ||
                (frequencyTopRight >= minFrequency && frequencyTopRight <= maxFrequency) ||
                (frequencyBottomLeft >= minFrequency && frequencyBottomLeft <= maxFrequency) ||
                (frequencyBottomRight >= minFrequency && frequencyBottomRight <= maxFrequency)) {
                // Preserve pixels within the specified frequency range
                // You can modify this part to preserve the original pixel values instead of setting them to zero
                image[i][j] = std::complex<double>(0.0, 0.0);
            }
        }
    }
}

void IMAGE_idealFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double cutoffFrequency, int filterType) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double wavelengthTopLeft = sqrt((i * i) + (j * j)); // Wavelength from top-left corner
            double wavelengthTopRight = sqrt((i * i) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from top-right corner
            double wavelengthBottomLeft = sqrt(((rows - i - 1) * (rows - i - 1)) + (j * j)); // Wavelength from bottom-left corner
            double wavelengthBottomRight = sqrt(((rows - i - 1) * (rows - i - 1)) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from bottom-right corner

            // Calculate frequency as the reciprocal of wavelength
            double frequencyTopLeft = 1.0 / wavelengthTopLeft;
            double frequencyTopRight = 1.0 / wavelengthTopRight;
            double frequencyBottomLeft = 1.0 / wavelengthBottomLeft;
            double frequencyBottomRight = 1.0 / wavelengthBottomRight;

            // Check if the frequency is below the cutoff frequency
            if (frequencyTopLeft <= cutoffFrequency &&
                frequencyTopRight <= cutoffFrequency &&
                frequencyBottomLeft <= cutoffFrequency &&
                frequencyBottomRight <= cutoffFrequency) {
                // Preserve pixels with frequencies below the cutoff
                // You can modify this part to preserve the original pixel values instead of setting them to zero
                if (filterType == LOW_PASS) {
                    image[i][j] = std::complex<double>(0.0, 0.0);
                }
                else {

                }

            }
            else {
                if (filterType == LOW_PASS) {

                }
                else {
                    image[i][j] = std::complex<double>(0.0, 0.0);
                }
            }
        }
    }
}

void IMAGE_gaussianFreqFilter(std::vector<std::vector<std::complex<double>>>& image, int rows, int cols, double cutoffFrequency,int filterType) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double wavelengthTopLeft = sqrt((i * i) + (j * j)); // Wavelength from top-left corner
            double wavelengthTopRight = sqrt((i * i) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from top-right corner
            double wavelengthBottomLeft = sqrt(((rows - i - 1) * (rows - i - 1)) + (j * j)); // Wavelength from bottom-left corner
            double wavelengthBottomRight = sqrt(((rows - i - 1) * (rows - i - 1)) + ((cols - j - 1) * (cols - j - 1))); // Wavelength from bottom-right corner

            // Calculate frequency as the reciprocal of wavelength
            double frequencyTopLeft = 1.0 / wavelengthTopLeft;
            double frequencyTopRight = 1.0 / wavelengthTopRight;
            double frequencyBottomLeft = 1.0 / wavelengthBottomLeft;
            double frequencyBottomRight = 1.0 / wavelengthBottomRight;

            // Calculate the Gaussian filter response based on the cutoff frequency
            double gaussianResponse = exp(-(frequencyTopLeft * frequencyTopLeft + frequencyTopRight * frequencyTopRight +
                frequencyBottomLeft * frequencyBottomLeft + frequencyBottomRight * frequencyBottomRight) /
                (2 * cutoffFrequency * cutoffFrequency));

            // Apply the Gaussian high-pass filter response to the pixel

            if (filterType == LOW_PASS) {
                image[i][j] *= (1.0 - gaussianResponse);
            }
            else if (filterType == HIGH_PASS) {
                image[i][j] *= gaussianResponse;
            }
        }
    }
}




/** @} */ // End of DFT group

/**
 * @defgroup FrequencyFilters Frequency Domain Filters
 * @{
 */

 // ... Add functions for frequency domain filtering here ...

 /** @} */ // End of FrequencyFilters group

 /**
  * @defgroup Thresholding Thresholding and Segmentation
  * @{
  */

  /**
   * @brief Apply grayscale thresholding to an image.
   *
   * @param[in] imgData Input image data.
   * @param[out] outData Output binary image after thresholding.
   * @param[in] threshold Threshold value.
   */
void IMAGE_grayscaleThreshold(const uint8_t* imgData, uint8_t* outData, int threshold) {
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = (imgData[i] > threshold) ? 0xFF : 0x00;
    }
}

/**
 * @brief Calculate the Otsu's threshold for binarization.
 *
 * @param[in] histogram Histogram of the image.
 * @param[in] totalPixels Total number of pixels in the image.
 * @return Calculated Otsu's threshold.
 */
int IMAGE_otsuThreshold(const int* histogram, int totalPixels) {
    int bestThreshold = 0;
    double bestVariance = 0.0;

    for (int t = 0; t <= MAX_INTENSITY; ++t) {
        int bgPixels = 0;
        int fgPixels = 0;
        int bgSum = 0;
        int fgSum = 0;

        for (int i = 0; i <= t; ++i) {
            bgPixels += histogram[i];
            bgSum += i * histogram[i];
        }

        for (int i = t + 1; i <= MAX_INTENSITY; ++i) {
            fgPixels += histogram[i];
            fgSum += i * histogram[i];
        }

        if (bgPixels == 0 || fgPixels == 0) {
            continue;
        }

        double bgMean = static_cast<double>(bgSum) / bgPixels;
        double fgMean = static_cast<double>(fgSum) / fgPixels;
        double betweenClassVariance = bgPixels * fgPixels * pow(bgMean - fgMean, 2) / totalPixels;

        if (betweenClassVariance > bestVariance) {
            bestVariance = betweenClassVariance;
            bestThreshold = t;
        }
    }

    return bestThreshold;
}

/**
 * @brief Apply Otsu's thresholding to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output binary image after Otsu's thresholding.
 */
void IMAGE_grayscaleOutsu(const uint8_t* imgData, uint8_t* outData) {
    int histogram[MAX_INTENSITY + 1];
    int totalPixels = IMAGE_histogramFormation(imgData, histogram);
    int threshold = IMAGE_otsuThreshold(histogram, totalPixels);

    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = (imgData[i] > threshold) ? MAX_INTENSITY : 0;
    }
}

/**
 * @brief Apply Sobel edge detection to an image.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output edge-detected image.
 */
void IMAGE_sobelFilter(const uint8_t* imgData, uint8_t* outData) {
    int sobelX[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };

    int sobelY[3][3] = {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1}
    };

    for (int y = 1; y < HEIGHT - 1; ++y) {
        for (int x = 1; x < WIDTH - 1; ++x) {
            int gradientX = 0;
            int gradientY = 0;

            for (int fy = 0; fy < 3; ++fy) {
                for (int fx = 0; fx < 3; ++fx) {
                    int pixelValue = imgData[(y + fy - 1) * WIDTH + (x + fx - 1)];
                    gradientX += pixelValue * sobelX[fy][fx];
                    gradientY += pixelValue * sobelY[fy][fx];
                }
            }

            // Calculate the gradient magnitude
            outData[y * WIDTH + x] = static_cast<uint8_t>(sqrt(gradientX * gradientX + gradientY * gradientY));
        }
    }
}

/**
 * @brief Perform region growing on an image from a seed point.
 *
 * @param[in] imgData Input image data.
 * @param[out] outData Output segmented image.
 * @param[in] seedX X-coordinate of the seed point.
 * @param[in] seedY Y-coordinate of the seed point.
 * @param[in] THRESHOLD Threshold for region growing.
 */
 /**
  * @brief Perform region growing on an image from a seed point.
  *
  * @param[in] imgData Input image data.
  * @param[out] outData Output segmented image.
  * @param[in] seedX X-coordinate of the seed point.
  * @param[in] seedY Y-coordinate of the seed point.
  * @param[in] THRESHOLD Threshold for region growing.
  */
  // Define a structure for a pixel coordinate
typedef struct {
    int x;
    int y;
} Pixel;

// Function to perform region growing
void regionGrowing(const uint8_t* imgData, uint8_t* outData, int seedX, int seedY, int THRESHOLD) {
    // Create a queue for pixel coordinates
    Pixel* pixelQueue = (Pixel*)malloc(WIDTH * HEIGHT * sizeof(Pixel));
    int front = 0;
    int rear = 0;

    // Initialize the seed point
    Pixel seedPixel = { seedX, seedY };
    pixelQueue[rear++] = seedPixel;

    uint8_t seedValue = imgData[INDEX_TO_OFFSET(seedX, seedY)];

    while (front != rear) {
        // Dequeue a pixel
        Pixel currentPixel = pixelQueue[front++];

        int currentX = currentPixel.x;
        int currentY = currentPixel.y;

        if (currentX >= 0 && currentX < WIDTH && currentY >= 0 && currentY < HEIGHT && outData[INDEX_TO_OFFSET(currentX, currentY)] == 0) {
            uint8_t currentValue = imgData[INDEX_TO_OFFSET(currentX, currentY)];
            outData[INDEX_TO_OFFSET(currentX, currentY)] = currentValue;

            // Define neighboring pixels
            int dx[] = { 1, -1, 0, 0 };
            int dy[] = { 0, 0, 1, -1 };

            for (int i = 0; i < 4; ++i) {
                int newX = currentX + dx[i];
                int newY = currentY + dy[i];

                if (newX >= 0 && newX < WIDTH && newY >= 0 && newY < HEIGHT) {
                    uint8_t neighborValue = imgData[INDEX_TO_OFFSET(newX, newY)];

                    if (abs(seedValue - neighborValue) <= THRESHOLD) {
                        // Enqueue the neighboring pixel
                        Pixel neighborPixel = { newX, newY };
                        pixelQueue[rear++] = neighborPixel;
                    }
                }
            }
        }
    }

    free(pixelQueue);
}


void IMAGE_basicSegmentation(const uint8_t* imgData, uint8_t* outData, uint8_t targetValue) {
    for (int i = 0; i < IMAGE_SIZE; ++i) {
        outData[i] = (imgData[i] == targetValue) ? 255 : 0;
    }
}


/** @} */ // End of Thresholding group

/**
 * @defgroup MorphologicalOps Morphological Image Processing
 * @{
 */

void IMAGE_erosion(const uint8_t* imgData, uint8_t* outData, int kernelSize) {
    int kernelRadius = kernelSize / 2;

    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            uint8_t minPixelValue = MAX_INTENSITY;

            for (int ky = -kernelRadius; ky <= kernelRadius; ++ky) {
                for (int kx = -kernelRadius; kx <= kernelRadius; ++kx) {
                    int offsetX = x + kx;
                    int offsetY = y + ky;

                    if (offsetX >= 0 && offsetX < WIDTH && offsetY >= 0 && offsetY < HEIGHT) {
                        int pixelValue = imgData[offsetY * WIDTH + offsetX];
                        if (pixelValue < minPixelValue) {
                            minPixelValue = pixelValue;
                        }
                    }
                }
            }

            outData[y * WIDTH + x] = minPixelValue;
        }
    }
}

//***************************13.3 Dilation***************************
void IMAGE_dilation(const uint8_t* imgData, uint8_t* outData, int kernelSize) {
    int kernelRadius = kernelSize / 2;

    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            uint8_t maxPixelValue = 0;

            for (int ky = -kernelRadius; ky <= kernelRadius; ++ky) {
                for (int kx = -kernelRadius; kx <= kernelRadius; ++kx) {
                    int offsetX = x + kx;
                    int offsetY = y + ky;

                    if (offsetX >= 0 && offsetX < WIDTH && offsetY >= 0 && offsetY < HEIGHT) {
                        int pixelValue = imgData[offsetY * WIDTH + offsetX];
                        if (pixelValue > maxPixelValue) {
                            maxPixelValue = pixelValue;
                        }
                    }
                }
            }

            outData[y * WIDTH + x] = maxPixelValue;
        }
    }
}

//***************************13.4 Opening***************************
void IMAGE_opening(const uint8_t* imgData, uint8_t* outData, int kernelSize) {
    // Apply erosion followed by dilation
    vector<uint8_t> tempData(IMAGE_SIZE);
    IMAGE_erosion(imgData, tempData.data(), kernelSize);
    IMAGE_dilation(tempData.data(), outData, kernelSize);
}

//***************************13.5 Closing************************
void IMAGE_closing(const uint8_t* imgData, uint8_t* outData, int kernelSize) {
    // Apply dilation followed by erosion
    vector<uint8_t> tempData(IMAGE_SIZE);
    IMAGE_dilation(imgData, tempData.data(), kernelSize);
    IMAGE_erosion(tempData.data(), outData, kernelSize);
}



// Structure to represent an equivalency in labels
typedef struct {
    uint8_t parent;
    uint8_t rank;
} Equivalency;

// Function to initialize the equivalency table
void initializeEquivalencyTable(Equivalency* table, int size) {
    for (int i = 0; i < size; i++) {
        table[i].parent = i;
        table[i].rank = 0;
    }
}

// Function to find the root of an equivalency class
uint8_t findRoot(Equivalency* table, uint8_t label) {
    while (table[label].parent != label) {
        label = table[label].parent;
    }
    return label;
}

// Function to union two equivalency classes by rank
void unionLabels(Equivalency* table, uint8_t label1, uint8_t label2) {
    uint8_t root1 = findRoot(table, label1);
    uint8_t root2 = findRoot(table, label2);

    if (root1 != root2) {
        if (table[root1].rank < table[root2].rank) {
            table[root1].parent = root2;
        }
        else if (table[root1].rank > table[root2].rank) {
            table[root2].parent = root1;
        }
        else {
            table[root2].parent = root1;
            table[root1].rank++;
        }
    }
}

// Function to perform connected component labeling
void labelConnectedComponents(const uint8_t* imgData, uint8_t* labeledData) {
    Equivalency equivalencyTable[256]; // 256 is used for 8-bit labels

    initializeEquivalencyTable(equivalencyTable, 256);

    uint8_t currentLabel = 1;

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int offset = INDEX_TO_OFFSET(x, y);
            if (imgData[offset] == 0) {
                continue; // Skip background pixels
            }

            uint8_t neighbors[4] = { 0, 0, 0, 0 }; // Top-left, above, top-right, left

            if (x > 0) {
                neighbors[0] = labeledData[offset - 1];
            }

            if (y > 0) {
                neighbors[1] = labeledData[offset - WIDTH];
                if (x > 0) {
                    neighbors[2] = labeledData[offset - WIDTH - 1];
                }
            }

            if (x > 0) {
                neighbors[3] = labeledData[offset - 1];
            }

            uint8_t minNeighbor = UINT8_MAX;

            for (int i = 0; i < 4; i++) {
                if (neighbors[i] > 0 && neighbors[i] < minNeighbor) {
                    minNeighbor = neighbors[i];
                }
            }

            if (minNeighbor == UINT8_MAX) {
                labeledData[offset] = currentLabel;
                currentLabel++;
            }
            else {
                labeledData[offset] = minNeighbor;
                for (int i = 0; i < 4; i++) {
                    if (neighbors[i] > 0 && neighbors[i] != minNeighbor) {
                        unionLabels(equivalencyTable, minNeighbor, neighbors[i]);
                    }
                }
            }
        }
    }

    // Second pass to resolve equivalencies
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            int offset = INDEX_TO_OFFSET(x, y);
            if (labeledData[offset] > 0) {
                labeledData[offset] = findRoot(equivalencyTable, labeledData[offset]);
            }
        }
    }
}





 /** @} */ // End of MorphologicalOps group
