#include <image_proccesing_lib.hpp>

using namespace std;
using namespace cv;


int main() {
    // Load the image
    Mat originalImage = imread("C:/Users/ozand/OneDrive/Masaüstü/test.png", cv::IMREAD_GRAYSCALE);
    resize(originalImage, originalImage, Size(WIDTH, HEIGHT));
    if (originalImage.empty()) {
        cerr << "Error: Could not open or find the image." << endl;
        return -1;
    }
    Mat outputImage(HEIGHT, WIDTH, CV_8U, Scalar(0));
    Mat outputImage2(HEIGHT, WIDTH, CV_8U, Scalar(0));
    IMAGE_grayscaleThreshold(originalImage.data, outputImage2.data, 150);
    cv::imshow("test", outputImage2);
    labelConnectedComponents(outputImage2.data, outputImage.data);

    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            printf("%d ", outputImage.data[INDEX_TO_OFFSET(x, y)]);
            outputImage.data[INDEX_TO_OFFSET(x, y)] = outputImage.data[INDEX_TO_OFFSET(x, y)] * 255;
        }
        printf("\n");
    }

    cv::imshow("test2", outputImage);
    cv::waitKey(0);

    return 0;
}
