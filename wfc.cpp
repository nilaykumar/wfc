
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

        Pattern() {
            this->pattern = new rgb_t[N * N];
            for(int i = 0; i < N * N; ++i)
                (this->pattern)[i] = make_colour(0, 0, 0);
        }

        // copy constructor
        // can we just call the assignment constructor from here?
        Pattern(const Pattern &other) {
            // can we just call the constructor instead?
            this->pattern = new rgb_t[N * N];
            for(int i = 0; i < N * N; ++i)
                (this->pattern)[i] = make_colour(0, 0, 0);
            for(int j = 0; j < N; ++j)
                for(int k = 0; k < N; ++k)
                    this->pattern[j * N + k] = other.pattern[j * N + k];
        }

        // assignment operator
        Pattern& operator=(const Pattern &other) {
            // can we just call the constructor instead?
            this->pattern = new rgb_t[N * N];
            for(int i = 0; i < N * N; ++i)
                (this->pattern)[i] = make_colour(0, 0, 0);
            for(int j = 0; j < N; ++j)
                for(int k = 0; k < N; ++k)
                    this->pattern[j * N + k] = other.pattern[j * N + k];
            return *this;
        }

        ~Pattern() {
            delete[] pattern;
        }

        inline bool operator== (const Pattern& rhs) {
            for(int j = 0; j < N; ++j)
                for(int k = 0; k < N; ++k)
                    if(rhs.pattern[j * N + k] != this->pattern[j * N + k])
                        return false;
            return true;
        }
        
};

// Why do I need this...?
int Pattern::N;

std::ostream& operator<<(std::ostream& stream, const Pattern& p) {
    stream << "Pattern{ ";
    for(int i = 0; i < Pattern::N; ++i) {
        for(int j = 0; j < Pattern::N; ++j)
            stream << "(" << +(p.pattern[i * Pattern::N + j].red) << ", " << +(p.pattern[i * Pattern::N + j].green) << ", " << +(p.pattern[i * Pattern::N + j].blue) << ") ";
        if(i != Pattern::N - 1)
            stream << "  //  ";
    }
    stream << "}";
    return stream;
}

class Cell {
    public:
        Pattern colorMatrix;
        int frequency;
        std::vector<Pattern> up, down, left, right;
    
    // empty constructor
    Cell() = default;

    // copy constructor
    Cell(const Cell &other) {
        this->colorMatrix = other.colorMatrix;
        this->frequency = other.frequency;
        this->up = other.up;
        this->down = other.down;
        this->left = other.left;
        this->right = other.right;
    }

};

std::ostream& operator<<(std::ostream &stream, const Cell &cell) {
    stream << "Cell{ \n\t" << cell.colorMatrix << "\n\tFrequency{ " << cell.frequency << " }\n";
    stream << "\tUp{ " << cell.up.size() << " }\n";
    stream << "\tDown{ " << cell.down.size() << " }\n";
    stream << "\tLeft{ " << cell.left.size() << " }\n";
    stream << "\tRight{ " << cell.right.size() << " }\n";
    stream << "}";
    return stream;
}

std::vector<Cell> createCellSheet(const int width, const int height, const int N, const rgb_t borderColor, bitmap_image inputImage, const std::string cellSheetFilename) {
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
    // loop through all the patterns, creating our vector of Cells
    std::cout << "LOG : Reading data into Cell data structures..." << std::endl;

    std::vector<Cell> cells;
    for(int i = 0; i < height; ++i)
        for(int j = 0; j < width; ++j) {
            //std::cout << "LOG : Processing (" << i << ", " << j << ")." << std::endl;
            // is this pattern already represented by a cell?
            bool present = false;
            Cell* duplicate = NULL;

            for(std::vector<Cell>::iterator it = cells.begin(); it != cells.end(); ++it)
                if((*it).colorMatrix == rawMatrix[i * width + j]) {
                    present = true;
                    duplicate = &(*it);
                    break;
                }
            // if this pattern hasn't been seen before, introduce it, and its neighboring patterns
            if(present == false) {
                Cell c;
                c.colorMatrix = rawMatrix[i * width + j];
                c.frequency = 1;
                c.up.push_back(rawMatrix[mod(i - 1, height) * width + j]);
                c.down.push_back(rawMatrix[mod(i + 1, height) * width + j]);
                c.left.push_back(rawMatrix[i * width + mod(j - 1, width)]);
                c.right.push_back(rawMatrix[i * width + mod(j + 1, width)]);
                cells.push_back(c);
                std::cout << "LOG: Found new cell at (" << i << ", " << j << ")." << std::endl;
                //std::cout << c << std::endl;
            } else {
                ++ duplicate->frequency;
                // check the four directions around this rawMatrix element and add those patterns as allowable to the cell if they aren't already there
                Pattern upPattern = rawMatrix[mod(i - 1, height) * width + j];
                if(std::find(duplicate->up.begin(), duplicate->up.end(), upPattern) == duplicate->up.end())
                    duplicate->up.push_back(upPattern);
                Pattern downPattern = rawMatrix[mod(i + 1, height) * width + j];
                if(std::find(duplicate->down.begin(), duplicate->down.end(), downPattern) == duplicate->down.end())
                    duplicate->down.push_back(downPattern);
                Pattern leftPattern = rawMatrix[i * width + mod(j - 1, width)];
                if(std::find(duplicate->left.begin(), duplicate->left.end(), leftPattern) == duplicate->left.end())
                    duplicate->left.push_back(leftPattern);
                Pattern rightPattern = rawMatrix[i * width + mod(j + 1, width)];
                if(std::find(duplicate->right.begin(), duplicate->right.end(), rightPattern) == duplicate->right.end())
                    duplicate->right.push_back(rightPattern);
                //std::cout << "LOG: Found duplicate cell at (" << i << ", " << j << ")." << std::endl;
                //std::cout << *duplicate << std::endl;
            }
        }
    std::cout << "LOG : ...done." << std::endl;
    std::cout << "LOG : Found " << cells.size() << " unique cells." << std::endl;
    delete[] rawMatrix;

    return cells;
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

    std::vector<Cell> cellTypes = createCellSheet(width, height, N, borderColor, inputImage, cellSheetFilename);
    for(Cell c : cellTypes)
        std::cout << c << std::endl;

    return 0;
}

