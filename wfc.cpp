
#include<cstdio>
#include<string>
#include<vector>
#include<algorithm>

#include "bitmap_image.hpp"

// a convenient modulus function for wrapping out of bounds indices
int mod(int a, int b) {
    if(a >= 0)
        return a % b;
    else
        return (a % b) + b;
}

// Load a bitmap file
bitmap_image loadImage(const std::string& filename) {

    std::cout << "LOG : Reading input bitmap from: "  << filename << std::endl;
    bitmap_image image = bitmap_image(filename);

    return image;
}

// convenience class for working with 2d arrays of colors, stored in a single array
class Pattern {
    public:
    static int N;
    rgb_t* pattern;

    ~Pattern() {
        delete[] pattern;
    }

    inline bool operator== (const Pattern& rhs) {
        if(N == 0) {
            std::cout << "ERROR : N = 0 !" << std::endl;
            return false;
        }
        for(int j = 0; j < N; ++j)
            for(int k = 0; k < N; ++k)
                if(rhs.pattern[j * N + k] != this->pattern[j * N + k])
                    return false;
        return true;
    }
};

int Pattern::N;

class Cell {
    public:
    Pattern colorMatrix;
    int frequency;
    std::vector<Pattern> up, down, left, right;
};

// TODO make this return the cell data
void createCellSheet(const int width, const int height, const int N, const rgb_t borderColor, bitmap_image inputImage, const std::string cellSheetFilename) {
    // Write a large bitmap containing every cell type
    // The cell types are separated by a 1px black border
    // and are N x N pixels wide.
    std::cout << "LOG : Preparing cell sheet." << std::endl;

    const int cellSheetWidth = width * N + width + 1;
    const int cellSheetHeight = height * N + height + 1;
    bitmap_image cellSheet = bitmap_image(cellSheetWidth, cellSheetHeight);
    cellSheet.clear();

    /*
        The algorithm is the following:
        1. while reading the input bitmap file populate a rawMatrix of patterns
        2. parse this matrix to create a vector of cells
    */

    Pattern* rawMatrix = new Pattern[width * height];
    // prepare the vertical borders
    for(int i = 0; i < width + 1; ++i)
        for(int j = 0; j < cellSheetHeight; ++j)
            cellSheet.set_pixel((N + 1) * i, j, borderColor);
    // prepare the horizontal borders
    for(int j = 0; j < height + 1; ++j)
        for(int i = 0; i < cellSheetWidth; ++i)
            cellSheet.set_pixel(i, (N + 1) * j, borderColor);
    // fill in the actual cells
    for(int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j) {
            // start filling in rawMatrix's (i,j)th cell's colorMatrix
            rgb_t* colorMatrix = new rgb_t[N * N];
            // loop through the pixels in the (i,j)th cell
            for(int k = 0; k < N; ++k)
                for(int l = 0; l < N; ++l) {
                    rgb_t color;
                    inputImage.get_pixel(mod(i + k, width), mod(j + l, height), color);
                    //std::cout << "LOG : Reading pixel (" << i + k << ", " << j + l << ") as color (" << +color.red << ", " << +color.green << ", " << +color.blue << ")." << std::endl;
                    // take into account the offset from the 1px borders
                    int ii = i + 1 + i * N + k;
                    int jj = j + 1 + j * N + l;
                    //std::cout << "LOG : Setting cell sheet pixel (" << ii << ", " << jj << ")." << std::endl;
                    cellSheet.set_pixel(ii, jj, color);
                    colorMatrix[k * N + l] = color;
                }
            rawMatrix[i * width + j].pattern = colorMatrix;
        }
    cellSheet.save_image(cellSheetFilename);
    std::cout << "LOG : Cell sheet saved to " << cellSheetFilename << std::endl;

    std::cout << "LOG : Reading data into Cell data structures..." << std::endl;

    // so much for the creation of the rawMatrix of patterns
    // loop through all the patterns, creating our vector of Cells
    std::vector<Cell> cells;
    for(int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j) {
            std::cout << "LOG : Processing (" << i << ", " << j << ")." << std::endl;
            //Pattern p = rawMatrix[i * width + j];
            // if this pattern is not already represented by a cell, add it
            bool found = false;
            for(std::vector<Cell>::iterator it = cells.begin(); it != cells.end(); ++it) {
                // write manual duplicate finder here
                if(it->colorMatrix == rawMatrix[i * width + j]) {
                    found = true;
                    std::cout << "LOG : Found duplicate pattern!" << std::endl;
                    // increment the frequency with which this pattern occurs
                    ++ it->frequency;
                    // check the four directions around this rawMatrix element and add those patterns as allowable to the cell if they aren't already there
                    Pattern upPattern = rawMatrix[mod(i-1, height) * width + j];
                    if(std::find(it->up.begin(), it->up.end(), upPattern) == it->up.end())
                        it->up.push_back(upPattern);
                    Pattern downPattern = rawMatrix[mod(i + 1, height) * width + j];
                    if(std::find(it->down.begin(), it->down.end(), downPattern) == it->down.end())
                        it->down.push_back(downPattern);
                    Pattern leftPattern = rawMatrix[i * width + mod(j - 1, width)];
                    if(std::find(it->left.begin(), it->left.end(), leftPattern) == it->left.end())
                        it->left.push_back(leftPattern);
                    Pattern rightPattern = rawMatrix[i * width + mod(j + 1, width)];
                    if(std::find(it->right.begin(), it->right.end(), rightPattern) == it->right.end())
                        it->right.push_back(rightPattern);
                    // we can now break since the cells as constructed iteratively like this must have unique patterns
                    break;
                }
            }
            // if it wasn't accounted for
            if(found == false) {
                Cell c;
                c.colorMatrix = rawMatrix[i * width + j];
                c.frequency = 1;
                cells.push_back(c);
                std::cout << "LOG : Adding a new cell type." << std::endl;
            }
        }
    std::cout << "LOG : ...done." << std::endl;
    std::cout << "LOG : Found " << cells.size() << " unique cells." << std::endl;
    delete[] rawMatrix;
}

int main() {

    // cell size
    const int N = 3;
    Pattern::N = N;
    // input bitmap
    const std::string filename = "City.bmp";
    // cell sheet output
    const std::string cellSheetFilename = "cellSheet.bmp";
    // cell sheet border color
    const rgb_t borderColor = make_colour(0, 0, 255);
    // list of CellTypes (pixel data and )

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

    createCellSheet(width, height, N, borderColor, inputImage, cellSheetFilename);

    return 0;
}

