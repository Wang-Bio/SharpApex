#include "../include/class.h"
#include "../include/initialization.h"
#include "../include/plot.h"
#include "../include/growth.h"
#include "../include/geo.h"


#include <iomanip>
#include <sstream>
#include <string>
//#include "../include/mathw.h"

std::string mode = "panel";
std::string real_contour_txt = "../input/1_0_contour_equal_distance_100.txt";

std::string doubleToString(double value, int precision) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();
}


int main(){
if(mode=="single_simulation"){
    
}
else if(mode=="simulation"){
    //initialization

    //plot::save_contour_txt("initial.txt",ct);
    //plot::save_contour_png("initial.png",ct);
    Contour *ct_real = new Contour;
    plot::input_real_contour(real_contour_txt,ct_real);
    plot::save_contour_png("real_contour.png",ct_real,-1,1,-0.5,1.5);

    std::vector<double> beta_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> ya_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> yb_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> beta_mid_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};

    for(int beta_i=0;beta_i<beta_vec.size();beta_i++){
        for(int beta_mid_i=0;beta_mid_i<beta_mid_vec.size();beta_mid_i++){
            if(beta_mid_vec[beta_mid_i]<beta_vec[beta_i]){
                //
            }
            else{
                for(int ya_i=0;ya_i<ya_vec.size();ya_i++){
                    for(int yb_i=0;yb_i<yb_vec.size();yb_i++){
                        if(ya_vec[ya_i]<yb_vec[yb_i]){

                        }
                        else{
                            Contour *ct = new Contour;
                            int element_number = 1000;
                            initial::circle(ct,element_number,5,0,0);
                            ct->end_step=200;
                            //parameter setting
                            boundary_apical=ya_vec[ya_i];
                            boundary_basal=yb_vec[yb_i];
                            heterogeneous_growth_rate_x=beta_vec[beta_i];
                            quadratic_mid_growth_rate=beta_mid_vec[beta_mid_i];
                            for(int step_i=0;step_i<ct->end_step;step_i++){
                                geo::contour_center(ct);
                                geo::contour_min_max(ct);
                                growth::contour(ct);
                                geo::continuity_check(ct);
                            }
                            std::cout<<"Simulation result for ya "<<ya_vec[ya_i]<<" yb "<<yb_vec[yb_i]<<" beta "<<beta_vec[beta_i]<<" beta_mid "<<beta_mid_vec[beta_mid_i]<<std::endl;
                            std::string save_png_file = "ya_"+doubleToString(ya_vec[ya_i],2)+"_yb_"+doubleToString(yb_vec[yb_i],2)+"_beta_"+doubleToString(beta_vec[beta_i],2)+"_beta_mid_"+doubleToString(beta_mid_vec[beta_mid_i],2)+".png";
                            plot::save_contour_png(save_png_file,ct);
                            //std::string save_txt_file = "ya_"+doubleToString(ya_vec[ya_i],2)+"_yb_"+doubleToString(yb_vec[yb_i],2)+"_beta_"+doubleToString(beta_vec[beta_i],2)+"_beta_mid_"+doubleToString(beta_mid_vec[beta_mid_i],2)+".txt";
                            //plot::save_contour_txt(save_txt_file,ct);
                            ct->~Contour();
                        }
                    }
                }
            }
        }       
    }
}
else if(mode=="panel_"){
    std::vector<std::vector<std::string>> all_subfigures;

    std::vector<double> beta_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> ya_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> yb_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};

    for(int beta_i=0;beta_i<beta_vec.size();beta_i++){
        std::vector<std::string> panel_subfigures;
        for(int ya_i=0;ya_i<ya_vec.size();ya_i++){
            for(int yb_i=0;yb_i<yb_vec.size();yb_i++){
                if(ya_vec[ya_i]<yb_vec[yb_i]){
                    panel_subfigures.push_back("blankImage.png");
                }
                else{
                    std::string subfigure = "ya_"+doubleToString(ya_vec[ya_i],2)+"_yb_"+doubleToString(yb_vec[yb_i],2)+"_beta_"+doubleToString(beta_vec[beta_i],2)+".png";
                    panel_subfigures.push_back(subfigure);
                }  
            }
        }
        all_subfigures.push_back(panel_subfigures);
        std::string save_panel = "panel_beta_"+doubleToString(beta_vec[beta_i],2)+".png";
        plot::create_image_panel(panel_subfigures,ya_vec.size(),yb_vec.size(),save_panel);
    }
}else if(mode=="panel"){
    std::vector<std::vector<std::string>> all_subfigures;

    std::vector<double> beta_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> ya_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> yb_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};
    std::vector<double> beta_mid_vec={0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.00};

    for(int beta_i=0;beta_i<beta_vec.size();beta_i++){
        std::vector<std::string> panel_subfigures;
        for(int ya_i=0;ya_i<ya_vec.size();ya_i++){
            for(int yb_i=0;yb_i<yb_vec.size();yb_i++){
                if(ya_vec[ya_i]<yb_vec[yb_i]){
                    panel_subfigures.push_back("blankImage.png");
                }
                else{
                    std::string subfigure = "ya_"+doubleToString(ya_vec[ya_i],2)+"_yb_"+doubleToString(yb_vec[yb_i],2)+"_beta_"+doubleToString(beta_vec[beta_i],2)+".png";
                    panel_subfigures.push_back(subfigure);
                }  
            }
        }
        all_subfigures.push_back(panel_subfigures);
        std::string save_panel = "panel_beta_"+doubleToString(beta_vec[beta_i],2)+".png";
        plot::create_image_panel(panel_subfigures,ya_vec.size(),yb_vec.size(),save_panel);
    }
}
else if(mode=="test"){
    cv::Mat blankImage = cv::Mat::zeros(cv::Size(1000, 1000), CV_8UC3);

    // If you want a different background color, e.g., white
     blankImage.setTo(cv::Scalar(255, 255, 255));

    // Save the image to a file
    cv::imwrite("blankImage.png", blankImage);
}
    /*
    //a single simulation of contour growth
    ct->end_step=200;
    ct->growth_rate=0.02;
    int show_step=99;
    for(int i=0;i<ct->end_step;i++){
        std::cout<<"Step "<<i<<std::endl;
        if(i%show_step==0){
                std::string save_png_file = "step_" + std::to_string(i) + ".png";
                plot::save_contour_png(save_png_file,ct);
        }
        //std::string save_txt_file = "step_" + std::to_string(i) + ".txt";
        //plot::save_contour_txt(save_txt_file,ct);
        geo::contour_center(ct);
        geo::contour_min_max(ct);
        growth::contour(ct);
        //geo::continuity_check(ct);
    }
    */
    return 0;
}