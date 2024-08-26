#include "../include/geo.h"

namespace geo{
    void contour_center(Contour *ct){
        double x_tmp=0, y_tmp=0;
        for(int i=0;i<ct->pt.size();i++){
            x_tmp+=ct->pt[i]->x;
            y_tmp+=ct->pt[i]->y;
        }
        ct->center.x = x_tmp/ct->pt.size();
        ct->center.y = y_tmp/ct->pt.size();
    }

    void contour_min_max(Contour *ct){
        double min_x=ct->pt[0]->x;
        double max_x=ct->pt[0]->x;
        double min_y=ct->pt[0]->y;
        double max_y=ct->pt[0]->y;

        for(int i=0;i<ct->pt.size();i++){
            if(min_x>ct->pt[i]->x){
                min_x=ct->pt[i]->x;
            }
            if(max_x<ct->pt[i]->x){
                max_x=ct->pt[i]->x;
            }
            if(min_y>ct->pt[i]->y){
                min_y=ct->pt[i]->y;
            }
            if(max_y<ct->pt[i]->y){
                max_y=ct->pt[i]->y;
            }
        }

        ct->min_x=min_x;
        ct->min_y=min_y;
        ct->max_x=max_x;
        ct->max_y=max_y;
        //std::cout<<"min_x "<<ct->min_x<<"; max_x "<<ct->max_x<<"; min_y "<<ct->min_y<<"; max_y "<<ct->max_y<<std::endl;
    }

    void continuity_check(Contour* ct){
        double distance_for_continuity=1.0E-1;
        for(int pi=0;pi<ct->pt.size();pi++){
            int index1=pi;
            int index2=pi+1;
            if(pi==ct->pt.size()-1){
                index2=0;
            }

            double distance = (ct->pt[index2]->y-ct->pt[index1]->y)*(ct->pt[index2]->y-ct->pt[index1]->y)+(ct->pt[index2]->x-ct->pt[index1]->x)*(ct->pt[index2]->x-ct->pt[index1]->x);
            
            if(distance>distance_for_continuity){
                //use interpolation to add new points 
                Point_element* new_pt = new Point_element;
                new_pt->x = (ct->pt[index2]->x+ct->pt[index1]->x)/2.0;
                new_pt->y = (ct->pt[index2]->y+ct->pt[index1]->y)/2.0;
                ct->pt.insert(ct->pt.begin()+index1+1,new_pt);
                //std::cout<<"new points added"<<std::endl;
            }        
        }
    }

    void normalized_contour(Contour* ct){

        for(auto p:ct->pt){
            double relative_y = (p->y-ct->min_y)/(ct->max_x-ct->min_y);
            double relative_x = (p->x-ct->center.x)/(ct->max_y-ct->min_y);
            Point_element *p_tmp = new Point_element;
            p_tmp->x= relative_x;
            p_tmp->y= relative_y;
            ct->normalized.push_back(p_tmp);
        }
    }

    void sampling_contour(Contour* ct){
        for(auto p:ct->normalized){
            
        }
    }
}