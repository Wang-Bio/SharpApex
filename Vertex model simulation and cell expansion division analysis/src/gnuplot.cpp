#include "../include/gnuplot.h"

namespace gnuplot{
    void vv(const std::vector<Vertex> &vv){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            //fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set size ratio -1\n");

            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 0.5 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < vv.size(); i++) {
                fprintf(pipe, "%f %f\n",vv[i].loc.x, vv[i].loc.y);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

        void vv(const std::vector<Vertex> &vv,const std::string & save_file){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set size ratio -1\n");
            fprintf(pipe, "set output '%s'\n",save_file.c_str());

            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 0.5 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < vv.size(); i++) {
                fprintf(pipe, "%f %f\n",vv[i].loc.x, vv[i].loc.y);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }
    
    void vv(const std::vector<double> &vd1, const std::vector<double> & vd2, const std::string & save_file){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set size ratio -1\n");
            fprintf(pipe, "set output '%s'\n",save_file.c_str());

            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 1.0 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < vd1.size(); i++) {
                fprintf(pipe, "%f %f\n",vd1[i], vd2[i]);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }
    

    void curvature_plot(vector<double> vd){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            //fprintf(pipe, "set terminal pngcairo\n");
            //fprintf(pipe, "set output '%s'\n",save_image_fileName.c_str());
            fprintf(pipe, "set yrange [:20]\n");

            
            fprintf(pipe,"set arrow from 0,0 to 100,0 nohead dashtype 2 lw 4 lc 'grey'\n");
            // Plot with points
            fprintf(pipe, "plot '-' with points pt 7 ps 1.5 lc rgb 'dark-green' notitle \n");
            // Output data points
            for(int i = 0; i < vd.size(); i++) {
                fprintf(pipe, "%d %f\n", i, vd[i]);
            }
            
            // End of data
            fprintf(pipe, "e\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    void boxplot(const std::vector<double>&, const std::string &){

    }
}