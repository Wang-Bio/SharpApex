#include "../include/plot.h"
#include "../include/growth.h"

namespace plot{
    void save_contour_txt(std::string save_txt_file, Contour *ct){
        std::ofstream fout(save_txt_file);
        fout<<"index x y"<<std::endl;
        for(int i=0;i<ct->pt.size();i++){
            fout<<i<<" "<<ct->pt[i]->x<<" "<<ct->pt[i]->y<<std::endl;
        }
        fout.close();
    }

    void save_contour_png(std::string save_png_file, Contour *ct){
        FILE *pipe = popen("gnuplot -persist", "w");
        double range = ct->initial_radius * pow(1+homogeneous_growth_rate,ct->end_step) + 50;
        //double range = 50;

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png size 1000,1000\n");
            fprintf(pipe, "set output '%s'\n",save_png_file.c_str());
            fprintf(pipe, "set size ratio -1\n");
            fprintf(pipe, "set yrange [-%f:%f]\n",range,range);
            fprintf(pipe, "set xrange [-%f:%f]\n",range,range);

            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 2 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < ct->pt.size(); i++) {
                fprintf(pipe, "%f %f\n", ct->pt[i]->x, ct->pt[i]->y);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    void save_contour_png(std::string png_file_name, Contour *ct,double range_x1, double range_x2, double range_y1, double range_y2){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png size 1000,1000\n");
            fprintf(pipe, "set output '%s'\n",png_file_name.c_str());
            fprintf(pipe, "set size ratio -1\n");
            //fprintf(pipe, "set yrange [-%f:%f]\n",0,1.0);
            //fprintf(pipe, "set xrange [-%f:%f]\n",-0.5,0.5);

            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 2 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < ct->pt.size(); i++) {
                fprintf(pipe, "%f %f\n", ct->pt[i]->x, ct->pt[i]->y);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    //this function assumed that the real contour is already normalized and sampled into 200 points
    void input_real_contour(std::string real_contour_file, Contour* ct){
        std::ifstream fin(real_contour_file,std::ios::in);
        if(!fin.is_open()){
            std::cout<<"Error: missing "<<real_contour_file<<" (the real contour file) !"<<std::endl;
            exit(1);
        }
        std::string description;
        fin>>description;
        fin>>description;
        fin>>description;
        fin>>description;
        
        for(int i=0;i<200;i++){
            int index_tmp;
            fin>>index_tmp;
            Point_element *p_tmp = new Point_element;
            fin>>p_tmp->x;
            fin>>p_tmp->y;
            int z_tmp;
            fin>>z_tmp;
            ct->pt.push_back(p_tmp);
        }
    }

    void create_image_panel(std::vector<std::string> imagePaths, int nRows, int nCols,std::string outputImagePath){
        std::vector<cv::Mat> images;
        for(auto path : imagePaths){
            cv::Mat img = cv::imread(path, cv::IMREAD_COLOR);
            if (img.empty()) {
                std::cerr << "Error reading image: " << path << std::endl;
                continue; // Handle error appropriately
            }
            images.push_back(img);
        }

        int singleWidth = images[0].cols;
        int singleHeight = images[0].rows;

         cv::Mat panel(nRows * singleHeight, nCols * singleWidth, images[0].type());

    // Copy each image into its position
    for (int i = 0; i < images.size(); ++i) {
        int row = i / nCols;
        int col = i % nCols;

        // Calculate where to place the image (top-left corner)
        int startX = col * singleWidth;
        int startY = row * singleHeight;

        // Copy the image to the corresponding position
        images[i].copyTo(panel(cv::Rect(startX, startY, singleWidth, singleHeight)));
    }

    // Save the final panel
    cv::imwrite(outputImagePath, panel);

    }

}