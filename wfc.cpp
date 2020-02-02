
#include<cstdio>
#include<string>

#include "bitmap_image.hpp"

// Load a bitmap file
bitmap_image loadImage(const std::string& filename) {
    std::cout << "LOG : Reading input bitmap from: "  << filename << std::endl;
    bitmap_image image = bitmap_image(filename);

    return image;
}



int main() {
    // Configuration settings here (eventually these should be read from a file)
    const int N = 3;
    const std::string filename = "City.bmp";
    const std::string cellSheetFilename = "cellSheet.bmp";

    std::cout << "LOG : Using " << filename << " as input, with N = " << N << " and writing cell sheet to " << cellSheetFilename << "." << std::endl;

    // Load the input bitmap
    bitmap_image inputImage = loadImage(filename);

    // Check to make sure the file was read correctly
    if(!inputImage) {
        std::cout << "ERROR : Could not load bitmap." << std::endl;
        return 1;
    }
    std::cout << "LOG : Bitmap loaded." << std::endl;

    const int height = inputImage.height();
    const int width = inputImage.width();

    // Number of (overlapping) NxN cells
    //const int numCells = height * width;

    // Write a large bitmap containing every cell type
    // The cell types are separated by a 1px black border
    // and are N x N pixels wide.
    std::cout << "LOG : Preparing cell sheet." << std::endl;
    const rgb_t borderColor = make_colour(0, 0, 255);
    const int cellSheetWidth = width * N + width + 1;
    const int cellSheetHeight = height * N + height + 1;
    bitmap_image cellSheet = bitmap_image(cellSheetWidth, cellSheetHeight);
    cellSheet.clear();

    // prepare the vertical borders
    for(int i = 0; i < width + 1; ++i)
        for(int j = 0; j < cellSheetWidth; ++j)
            cellSheet.set_pixel((N + 1) * i, j, borderColor);
    // prepare the horizontal borders
    for(int j = 0; j < height + 1; ++j)
        for(int i = 0; i < cellSheetHeight; ++i)
            cellSheet.set_pixel(i, (N + 1) * j, borderColor);
    // fill in the actual cells
    for(int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j) {
            // loop through the pixels in the (i,j)th cell
            for(int k = 0; k < N; ++k)
                for(int l = 0; l < N; ++l) {
                    rgb_t color;
                    inputImage.get_pixel((i + k) % width, (j + l) % height, color);
                    //std::cout << "LOG : Reading pixel (" << i + k << ", " << j + l << ") as color (" << +color.red << ", " << +color.green << ", " << +color.blue << ")." << std::endl;
                    // take into account the offset from the 1px borders
                    int ii = i + 1 + i * N + k;
                    int jj = j + 1 + j * N + l;
                    //std::cout << "LOG : Setting cell sheet pixel (" << ii << ", " << jj << ")." << std::endl;
                    cellSheet.set_pixel(ii, jj, color);
                }
        }
    cellSheet.save_image(cellSheetFilename);
    std::cout << "LOG : Cell sheet saved to " << cellSheetFilename << std::endl;

    return 0;
}

