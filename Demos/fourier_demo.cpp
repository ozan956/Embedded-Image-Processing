#include <image_proccesing_lib.hpp>

using namespace std;
using namespace cv;


int main() {
    // Load the image
    Mat originalImage = imread("C:/Users/ozand/OneDrive/Masaüstü/test.png", cv::IMREAD_GRAYSCALE);

    if (originalImage.empty()) {
        cerr << "Error: Could not open or find the image." << endl;
        return -1;
    }

    resize(originalImage, originalImage, Size(512, 512));
    // Define the desired size for the resized image
    int desiredWidth = 512; // Change this to your desired width
    int desiredHeight = 512;
    printf("%d %d", originalImage.cols, originalImage.rows);
    // Calculate the padding dimensions
    int padWidth = desiredWidth - originalImage.cols;
    int padHeight = desiredHeight - originalImage.rows;

    // Create a new image with the desired size and pad the original image
    Mat paddedImage(desiredHeight, desiredWidth, CV_8U, Scalar(0));
    originalImage.copyTo(paddedImage(Rect(padWidth / 2, padHeight / 2, originalImage.cols, originalImage.rows)));

    // Perform the FFT on the padded image
    int rows = paddedImage.rows;
    int cols = paddedImage.cols;


    vector<vector<complex<double>>>
        image(rows, vector<complex<double>>(cols));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            image[i][j] = complex<double>(paddedImage.at<uchar>(i, j), 0.0);
        }
    }

    // Perform the 2D FFT
    fft2D(image);

    // Calculate the magnitude spectrum and shift it for visualization
    cv::Mat spectrum(rows, cols, CV_64F);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double magnitude = std::abs(image[i][j]);
            spectrum.at<double>((i + rows / 2) % rows, (j + cols / 2) % cols) = log(1 + magnitude);
        }
    }

    // Normalize the spectrum
    cv::Mat normalizedImage(spectrum.size(), CV_64F);
    manualNormalize(spectrum, normalizedImage);

    // Display the spectrum
    cv::imshow("Spectrum", normalizedImage);

    IMAGE_gaussianFreqFilter(image, 512, 512, 1.0 / 100, HIGH_PASS);

    // Calculate the magnitude spectrum and shift it for visualization
    cv::Mat spectrum2(rows, cols, CV_64F);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double magnitude = std::abs(image[i][j]);
            spectrum2.at<double>((i + rows / 2) % rows, (j + cols / 2) % cols) = log(1 + magnitude);
        }
    }

    // Normalize the spectrum
    cv::Mat normalizedImage2(spectrum2.size(), CV_64F);
    manualNormalize(spectrum2, normalizedImage2);

    // Display the spectrum
    cv::imshow("Spectrum2", normalizedImage2);



    // Perform the inverse Fourier Transform to reconstruct the filtered image
    Mat filteredImage;
    ifft2D(image, filteredImage);

    // Remove the padding from the filtered image
    Rect roi(padWidth / 2, padHeight / 2, originalImage.cols, originalImage.rows);
    Mat finalFilteredImage = filteredImage(roi).clone();

    // Display the filtered image
    cv::imshow("Filtered Image", finalFilteredImage);
    cv::imshow("*Image", originalImage);
    // Calculate and display the magnitude spectra
    Mat magnitudeSpectrumOriginal, magnitudeSpectrumFiltered;

    cv::waitKey(0);

    return 0;
}
