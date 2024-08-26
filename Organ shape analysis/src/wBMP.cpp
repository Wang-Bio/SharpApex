#include "../include/wBMP.h"

#pragma pack(push, 1)
struct BMPHeader {
    uint16_t bfType;
    uint32_t bfSize;
    uint16_t bfReserved1;
    uint16_t bfReserved2;
    uint32_t bfOffBits;
};

struct BMPInfoHeader {
    uint32_t biSize;
    int32_t  biWidth;
    int32_t  biHeight;
    uint16_t biPlanes;
    uint16_t biBitCount;
    uint32_t biCompression;
    uint32_t biSizeImage;
    int32_t  biXPelsPerMeter;
    int32_t  biYPelsPerMeter;
    uint32_t biClrUsed;
    uint32_t biClrImportant;
};
#pragma pack(pop)

// Bitmap font for digits 0-9, 5x7 matrix
const uint8_t digitFont[10][7][5] = {
    { // 0
        {0,1,1,1,0},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    },
    { // 1
        {0,0,1,0,0},
        {0,1,1,0,0},
        {1,0,1,0,0},
        {0,0,1,0,0},
        {0,0,1,0,0},
        {0,0,1,0,0},
        {1,1,1,1,1}
    },
    { // 2
        {0,1,1,1,0},
        {1,0,0,0,1},
        {0,0,0,0,1},
        {0,0,0,1,0},
        {0,0,1,0,0},
        {0,1,0,0,0},
        {1,1,1,1,1}
    },
    { // 3
        {0,1,1,1,0},
        {1,0,0,0,1},
        {0,0,0,0,1},
        {0,0,1,1,0},
        {0,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    },
    { // 4
        {0,0,0,1,0},
        {0,0,1,1,0},
        {0,1,0,1,0},
        {1,0,0,1,0},
        {1,1,1,1,1},
        {0,0,0,1,0},
        {0,0,0,1,0}
    },
    { // 5
        {1,1,1,1,1},
        {1,0,0,0,0},
        {1,1,1,1,0},
        {0,0,0,0,1},
        {0,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    },
    { // 6
        {0,1,1,1,0},
        {1,0,0,0,0},
        {1,1,1,1,0},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    },
    { // 7
        {1,1,1,1,1},
        {0,0,0,0,1},
        {0,0,0,1,0},
        {0,0,1,0,0},
        {0,1,0,0,0},
        {0,1,0,0,0},
        {0,1,0,0,0}
    },
    { // 8
        {0,1,1,1,0},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    },
    { // 9
        {0,1,1,1,0},
        {1,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,1},
        {0,0,0,0,1},
        {1,0,0,0,1},
        {0,1,1,1,0}
    }
};

namespace wBMP{
    int rangeTics(int range){
        if(range<5){
            return 1;
        }
        else if(range<10){
            return 2;
        }
        else if(range<30){
            return 5;
        }
        else if(range<50){
            return 10;
        }
        else if(range<100){
            return 20;
        }else{
            return 10*rangeTics(range/10);
        }
    }

    // Function to draw a single digit at a specified position
    void drawDigit(std::vector<uint8_t>& image, int width, int height, int x, int y, int digit, uint8_t r, uint8_t g, uint8_t b) {
        if (digit < 0 || digit > 9) return; // Invalid digit

        for (int i = 0; i < 7; ++i) { // 7 rows
            for (int j = 0; j < 5; ++j) { // 5 columns
                if (digitFont[digit][i][j] == 1) {
                    int pixelX = x + j;
                    int pixelY = y + i;
                    if (pixelX >= 0 && pixelX < width && pixelY >= 0 && pixelY < height) {
                        int index = (pixelY * width + pixelX) * 3;
                        image[index + 0] = b; // Blue
                        image[index + 1] = g; // Green
                        image[index + 2] = r; // Red
                    }
                }
            }
        }
    }

    // Function to draw a negative sign at a specified position
    void drawNegativeSign(std::vector<uint8_t>& image, int width, int height, int x, int y, uint8_t r, uint8_t g, uint8_t b) {
        for (int i = 0; i < 3; ++i) { // 3 pixels long horizontal line
            int pixelX = x + i;
            int pixelY = y + 3; // Middle of the 7-pixel height
            if (pixelX >= 0 && pixelX < width && pixelY >= 0 && pixelY < height) {
                int index = (pixelY * width + pixelX) * 3;
                image[index + 0] = b; // Blue
                image[index + 1] = g; // Green
                image[index + 2] = r; // Red
            }
        }
    }

    // Function to draw a number at a specified position
    void drawNumber(std::vector<uint8_t>& image, int width, int height, int x, int y, int number, uint8_t r, uint8_t g, uint8_t b) {
        std::string numStr = std::to_string(number);
        int spacing = 6; // Space between digits
        int startX = x;

        if (number < 0) {
            drawNegativeSign(image, width, height, startX, y, r, g, b);
            startX += spacing; // Move start position to the right for digits
            numStr = numStr.substr(1); // Remove negative sign from the string
        }

        for (size_t i = 0; i < numStr.size(); ++i) {
            int digit = numStr[i] - '0';
            drawDigit(image, width, height, startX + i * spacing, y, digit, r, g, b);
        }
    }

    std::vector<int> generateTicks(int minValue, int maxValue){
                bool debugMode=0;
        std::vector<int> Ticks;
        int Tics = rangeTics(maxValue-minValue);
        //calculate the minTics
        int minTics;
        if(minValue%Tics==0){
            minTics=minValue;
        }else{
            if(minValue>0){
                minTics=minValue/Tics*Tics;
            }
            else{
                minTics=minValue/Tics*Tics-Tics;
            }
            
        }
        int maxTics;
        if(maxValue%Tics==0){
            maxTics=maxValue;
        }else{
            maxTics=maxValue/Tics*Tics+Tics;
        }
                if(debugMode==1){std::cout<<"minValue, maxValue: "<<minValue<<","<<maxValue<<"Tics, minTics, maxTics "<<Tics<<","<<minTics<<","<<maxTics<<std::endl;};
        for(int i=0;i<((maxTics-minTics)/Tics+1);i++){
            Ticks.push_back(minTics+i*Tics);
        }
            if(debugMode==1){wFile::vec_cout(Ticks);};
        return Ticks;
    }

    void drawSmoothCircle(std::vector<uint8_t>& image, int width, int height, int centerX, int centerY, int radius, uint8_t r, uint8_t g, uint8_t b) {
        for (int y = -radius; y <= radius; ++y) {
            for (int x = -radius; x <= radius; ++x) {
                float distance = std::sqrt(x * x + y * y);
                if (distance <= radius + 0.5) { // Add a small value to smooth the edges
                    int pixelX = centerX + x;
                    int pixelY = centerY + y;
                    if (pixelX >= 0 && pixelX < width && pixelY >= 0 && pixelY < height) {
                        int index = (pixelY * width + pixelX) * 3;
                        image[index + 0] = b; // Blue
                        image[index + 1] = g; // Green
                        image[index + 2] = r; // Red
                    }
                }
            }
        }
    }

    void drawLine(std::vector<uint8_t>& image, int width, int height, int x0, int y0, int x1, int y1, uint8_t r, uint8_t g, uint8_t b) {
        int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
        int dy = -abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
        int err = dx + dy, e2;

        while (true) {
            if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) {
                int index = (y0 * width + x0) * 3;
                image[index + 0] = b; // Blue
                image[index + 1] = g; // Green
                image[index + 2] = r; // Red
            }
            if (x0 == x1 && y0 == y1) break;
            e2 = 2 * err;
            if (e2 >= dy) { err += dy; x0 += sx; }
            if (e2 <= dx) { err += dx; y0 += sy; }
        }
    }

    void drawContour(std::vector<std::pair<double,double>> DotplotInput,std::string recordBMP){
                    bool debugMode=0;
            //1. normalize and map data points to the image dimensions
            double x_min = static_cast<int>(floor(wPair::findMinX(DotplotInput)));
            double x_max = static_cast<int>(ceil(wPair::findMaxX(DotplotInput)));
            double y_min = static_cast<int>(floor(wPair::findMinY(DotplotInput)));
            double y_max = static_cast<int>(ceil(wPair::findMaxY(DotplotInput)));
            double width_x = x_max-x_min;
            double height_y = y_max-y_min;
                    if(debugMode==1){std::cout<<"width_x "<<width_x<<", height_y "<<height_y<<std::endl;};
            const int width = 640;
            const int height = 480;

            std::vector<uint8_t> image(width * height * 3, 255); // White background

            for(const auto& dot: DotplotInput){
                int x = static_cast<int>((dot.first-x_min)/width_x*(width-1));
                int y = static_cast<int>((height_y-dot.second+y_min)/height_y*(height-1));
                //drawDot(image, width, x, y, 0, 100, 0); 
                drawSmoothCircle(image,width,height,x,y,5,0,100,0);
                        if(debugMode==1){std::cout<<"Dot "<<dot.first<<","<<dot.second<<" drawn at"<<x<<","<<y<<std::endl;};   
            }

             // Add axes
            
            int x_axis = static_cast<int>((-y_min) / height_y * (height - 1));
            int y_axis = static_cast<int>((-x_min) / width_x * (width - 1));
            drawLine(image, width, height, y_axis, 0, width - 1, 0, 0, 0, 0); // upper X-axis (Black color)
            drawLine(image, width, height, y_axis, height - 1, width - 1, height - 1, 0, 0, 0); // lower X-axis (Black color)
            drawLine(image, width, height, y_axis, 0, y_axis, height - 1, 0, 0, 0); // left Y-axis (Black color)
            drawLine(image, width, height, width - 1, 0, width - 1, height - 1, 0, 0, 0); // right Y-axis (Black color)

             // Add tics
                if(debugMode==1){std::cout<<"y_min "<<y_min<<" y_max "<<y_max<<" x_min "<<x_min<<" x_max "<<x_max<<std::endl;};
            std::vector<int> y_tics = generateTicks(y_min,y_max);   
            std::vector<int> x_tics = generateTicks(x_min,x_max);
                // draw xtics
            for(int i=0;i<x_tics.size();i++){
                int x_tic = i*(width / (x_tics.size()-1));
                drawLine(image, width, height, x_tic, height + 5, x_tic, height - 5, 0, 0, 0);
                int label = x_tics[i];
                drawNumber(image, width, height, x_tic, height - 10, label, 0, 0, 0);
            }
                //draw ytics
            for (int i = 0; i <= y_tics.size(); ++i) {
                int y_tic = i * (height / (y_tics.size()-1));
                drawLine(image, width, height, y_axis - 5, y_tic, y_axis + 5, y_tic, 0, 0, 0);
                int label = y_tics[y_tics.size()-1-i];
                drawNumber(image, width, height, y_axis + 10, y_tic, label, 0, 0, 0);
            }
            // Margin sizes
            const int marginSize = 40;
            const int newWidth = width + 2 * marginSize;
            const int newHeight = height + 2 * marginSize;

            std::vector<uint8_t> newImage(newWidth * newHeight * 3, 255); // White background with margin

            // Copy original image to the center of the new image with margins
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    int originalIndex = (y * width + x) * 3;
                    int newIndex = ((y + marginSize) * newWidth + (x + marginSize)) * 3;
                    newImage[newIndex + 0] = image[originalIndex + 0]; // Blue
                    newImage[newIndex + 1] = image[originalIndex + 1]; // Green
                    newImage[newIndex + 2] = image[originalIndex + 2]; // Red
                }
            }

            BMPHeader bmpHeader;
            bmpHeader.bfType = 0x4D42; // "BM"
            bmpHeader.bfSize = sizeof(BMPHeader) + sizeof(BMPInfoHeader) + newImage.size();
            bmpHeader.bfReserved1 = 0;
            bmpHeader.bfReserved2 = 0;
            bmpHeader.bfOffBits = sizeof(BMPHeader) + sizeof(BMPInfoHeader);

            BMPInfoHeader bmpInfoHeader;
            bmpInfoHeader.biSize = sizeof(BMPInfoHeader);
            bmpInfoHeader.biWidth = newWidth;
            bmpInfoHeader.biHeight = -newHeight; // Negative to indicate top-down DIB
            bmpInfoHeader.biPlanes = 1;
            bmpInfoHeader.biBitCount = 24;
            bmpInfoHeader.biCompression = 0;
            bmpInfoHeader.biSizeImage = newImage.size();
            bmpInfoHeader.biXPelsPerMeter = 0;
            bmpInfoHeader.biYPelsPerMeter = 0;
            bmpInfoHeader.biClrUsed = 0;
            bmpInfoHeader.biClrImportant = 0;

            std::ofstream file(recordBMP, std::ios::binary);
            if (!file) {
                std::cerr << "Unable to open file for writing\n";
                exit(1);
            }

            file.write(reinterpret_cast<const char*>(&bmpHeader), sizeof(bmpHeader));
            file.write(reinterpret_cast<const char*>(&bmpInfoHeader), sizeof(bmpInfoHeader));
            file.write(reinterpret_cast<const char*>(newImage.data()), newImage.size());

            if(debugMode==1){std::cout << "Image saved as "<<recordBMP<<std::endl;};
        }
        
        void drawCurvature(std::vector<double> DotplotInput,std::string recordBMP){
            std::vector<std::pair<double,double>> dotplot;
            for(int i=0;i<DotplotInput.size();i++){
                dotplot.push_back(std::make_pair(i,DotplotInput[i]));
            }
            drawContour(dotplot, recordBMP);
        }


}