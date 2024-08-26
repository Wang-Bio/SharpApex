#include "../include/growth.h"
#include "../include/mathw.h"


std::string growth_mode="biregion_growth_quadratic";
double homogeneous_growth_rate=0.02;
double boundary_apical=1.0;
double boundary_basal=0.3;
double heterogeneous_growth_rate_x=0.01;

double quadratic_mid_growth_rate=0.8;

namespace growth{
    void contour(Contour* ct){
        for(int i=0;i<ct->pt.size();i++){
            Growth_vector gv;
            if(growth_mode=="homogeneous"){
                gv=growth::homogeneous_growth(ct->pt[i],ct);
            }
            else if(growth_mode=="biregion_growth"){
                gv=growth::biregion_growth(ct->pt[i],ct);
            }
            else if(growth_mode=="polarized_growth"){
                gv=growth::polarized_growth(ct->pt[i],ct);
            }
            else if(growth_mode=="biregion_growth_quadratic"){
                gv=growth::biregion_growth_quadratic(ct->pt[i],ct);
            }
            else{
                std::cout<<"Fatal error: No growth mode was selected ! ("<<growth_mode<<")"<<std::endl;
            }
            ct->pt[i]->x += gv.x;
            ct->pt[i]->y += gv.y;
        }
    }

    Growth_vector homogeneous_growth(Point_element *pt, Contour *ct){
        Growth_vector gv;
        gv.x = homogeneous_growth_rate * pt->x;
        gv.y = homogeneous_growth_rate * pt->y;
        return gv;
    }

    Growth_vector polarized_growth(Point_element *pt, Contour *ct){
        Growth_vector gv;
        gv.x = homogeneous_growth_rate * heterogeneous_growth_rate_x * pt->x;
        gv.y = homogeneous_growth_rate * pt->y;
        return gv;
    }

    Growth_vector biregion_growth(Point_element *pt, Contour *ct){
        Growth_vector gv;
        double relative_y = (pt->y-ct->min_y)/(ct->max_y-ct->min_y);
        if(relative_y<boundary_basal){
            //homogeneous growth in basal region
            gv=homogeneous_growth(pt,ct);
        }
        else if(relative_y>boundary_apical){
            //vertcial growth in apical region
            gv=polarized_growth(pt,ct);
        }
        else{
            //intermediate region: intermedate growth
            //(r,y_b), (r/beta,y_a)
            double slope_tmp = (homogeneous_growth_rate*heterogeneous_growth_rate_x-homogeneous_growth_rate)/(boundary_apical-boundary_basal);
            double intercept_tmp = (homogeneous_growth_rate*boundary_apical-homogeneous_growth_rate*heterogeneous_growth_rate_x*boundary_basal)/(boundary_apical-boundary_basal);
            double r_x_growth_rate = slope_tmp*relative_y+intercept_tmp;
            
            gv.x = r_x_growth_rate * pt->x;
            gv.y = homogeneous_growth_rate * pt->y;
        }
        return gv;

    }

    Growth_vector biregion_growth_quadratic(Point_element *pt, Contour *ct){
        Growth_vector gv;
        double relative_y = (pt->y-ct->min_y)/(ct->max_y-ct->min_y);
        if(relative_y<boundary_basal){
            //homogeneous growth in basal region
            gv=homogeneous_growth(pt,ct);
        }
        else if(relative_y>boundary_apical){
            //vertcial growth in apical region
            gv=polarized_growth(pt,ct);
        }
        else{
            //intermediate region: intermedate growth
            //(r,y_b), (r/beta,y_a)
            double r=homogeneous_growth_rate;
            double beta=heterogeneous_growth_rate_x;
            double rmid=homogeneous_growth_rate*quadratic_mid_growth_rate;
            double ya=boundary_apical;
            double yb=boundary_basal;
            double quadratic_a = (2.0*r+2.0*r*beta-4.0*rmid)/((ya-yb)*(ya-yb));
            double quadratic_b = (yb*r*(-1.0-3.0*beta)+ya*r*(-3.0-beta)+4.0*rmid*(yb+ya))/((ya-yb)*(ya-yb));
            double quadratic_c = (r*ya*ya+r*yb*yb*beta+r*ya*yb*(1+beta)-4.0*rmid*ya*yb)/((ya-yb)*(ya-yb));
            double r_x_growth_rate = quadratic_a*relative_y*relative_y+quadratic_b*relative_y+quadratic_c;
            
            gv.x = r_x_growth_rate * pt->x;
            gv.y = homogeneous_growth_rate * pt->y;
        }
        return gv;
    }
    /*
    Growth_vector polarized_growth(double angles, Point_element *pt, Contour *ct){
        Growth_vector gv;
        double relative_y = (pt->y-ct->min_y)/(ct->max_y-ct->min_y);
        double relative_x = (pt->x-ct->min_x)/(ct->max_x-ct->min_x);
        double k_line = pt->y/pt->x;
        if(relative_y<0.5){
            //this point belong to the basal part
            double growth_vector_x = newtonRaphson_lower(4.0,k_line,1000,1.0E-7);
            double growth_vector_y = growth_vector_x*k_line;
            gv.x = homogeneous_growth_rate * growth_vector_y;
            gv.y = homogeneous_growth_rate * growth_vector_x;
        std::cout<<"For point "<<pt->x<<","<<pt->y<<" its k is "<<k_line<<" and its growth vector is "<<growth_vector_x<<","<<growth_vector_y<<std::endl;

        }
        else{
            //this point belongs to the apical part
            //triangle setting
            double growth_vector_x = newtonRaphson_upper(4.0,k_line,1000,1.0E-7);
            double growth_vector_y = growth_vector_x*k_line;
            gv.x = homogeneous_growth_rate * growth_vector_y;
            gv.y = homogeneous_growth_rate * growth_vector_x;
        std::cout<<"For point "<<pt->x<<","<<pt->y<<" its k is "<<k_line<<" and its growth vector is "<<growth_vector_x<<","<<growth_vector_y<<std::endl;
            
        }

        return gv;
    }

    Growth_vector biregion_growth(double angles, Point_element *pt, Contour *ct){
        Growth_vector gv;
        double relative_y = (pt->y-ct->min_y)/(ct->max_y-ct->min_y);
        double relative_x = (pt->x-ct->min_x)/(ct->max_x-ct->min_x);
        //std::cout<<"relative_y "<<relative_y<<std::endl;
        if(relative_y<0.5){
            //this point belong to the basal part
            gv.x = homogeneous_growth_rate * std::cos(angles);
            gv.y = homogeneous_growth_rate * std::sin(angles);
            return gv;
        }
        else{
            //this point belongs to the apical part
            //the apical region growth vector is set to be sin(x+pi/2), we need to solve the equation sin(x+pi/2) = tan(theta)*x
            gv.x = (1.0/relative_y-1.0)*homogeneous_growth_rate * std::cos(angles);
            gv.y = (1.0/relative_y)*homogeneous_growth_rate * std::sin(angles);
            return gv;
        }
    }
    */
}