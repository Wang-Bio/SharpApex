/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/force2dv.h"
#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"

using namespace std;

namespace force{

//calculate edge elastic force and area elastic force for each vertices and execute the vertex motion

//F_i = -sigma_L * L_k (vector) / L_k (scalar)
void line_elastic_force(Organ* p_g){
    //calculate the elastic force for each edge
    
    //has two potential energy modes for elastic forces of edges: 1. simple; 2. L_std
    if(mechanics_mode=="simple"||mechanics_mode=="S_std_spatial"){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            _vec<double> relative = p_g->p_v[p_g->p_l[li]->vi[0]]->loc-p_g->p_v[p_g->p_l[li]->vi[1]]->loc;
            double dist = relative.norm();
            p_g->p_l[li]->length = dist;
                
            _vec<double> frc;
            if(p_g->p_l[li]->IsOutermost==0){
                frc = (-1.0)*sigma_L* relative/dist;
            }
            else{
                frc = (-1.0)*sigma_O* relative/dist;
            }
                
            p_g->p_v[p_g->p_l[li]->vi[0]]->frc_edge += frc;
            p_g->p_v[p_g->p_l[li]->vi[1]]->frc_edge -= frc;
            p_g->p_l[li]->edgeForce = frc.norm();

        }
    }
    
    else if(mechanics_mode == "L_std"){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            _vec<double> relative = p_g->p_v[p_g->p_l[li]->vi[0]]->loc-p_g->p_v[p_g->p_l[li]->vi[1]]->loc;
            double dist = relative.norm();
            p_g->p_l[li]->length = dist;
            _vec<double> frc;
            if(p_g->p_l[li]->IsOutermost==0){
                frc = (-1.0)*sigma_L* (dist -L_std)*relative/dist;
            }
            else{
                frc = (-1.0)*sigma_O* (dist -L_std)*relative/dist;
            }
                
            p_g->p_v[p_g->p_l[li]->vi[0]]->frc_edge += frc;
            p_g->p_v[p_g->p_l[li]->vi[1]]->frc_edge -= frc;
            p_g->p_l[li]->edgeForce = frc.norm();
        }
    }
    
    else if(mechanics_mode == "sigma_O_spatial"){
        //cout<<"p_g->p_l[0]->vi[1] "<<p_g->p_l[0]->vi[1]<<endl;
        //cout<<"p_g->p_v[p_g->p_l[0]->vi[0]]->loc.y "<<p_g->p_v[p_g->p_l[0]->vi[0]]->loc.y<<endl;
        //cout<<"p_g->p_v[0]->loc.y "<<p_g->p_v[0]->loc.y<<endl;
        vector<double> line_y_position;
        //double line_y_max=(p_g->p_v[p_g->p_l[0]->vi[0]]->loc.y+p_g->p_v[p_g->p_l[0]->vi[1]]->loc.y)/2.0;
        double line_y_max;
        double line_y_min;
        //double line_y_min=line_y_max;
        //cout<<"line y max: "<<line_y_max<<"; line y min: "<<line_y_min<<endl;
        for(int li=0;li<(int)p_g->p_l.size();li++){
            
            double line_y_tmp=0;
            line_y_tmp += p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y;
            line_y_tmp += p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y;
            line_y_tmp = line_y_tmp/2.0;
            if(line_y_max<line_y_tmp){
                line_y_max = line_y_tmp;
            }
            if(line_y_min>line_y_tmp){
                line_y_min = line_y_tmp;
            }
            //cout<<"For line "<<li<<" the line_y_position is "<<line_y_tmp<<endl;
        }

        for(int li=0;li<(int)p_g->p_l.size();li++){
            double line_y_tmp=0;
            line_y_tmp += p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y;
            line_y_tmp += p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y;
            line_y_tmp = line_y_tmp/2.0;
            line_y_tmp = (line_y_tmp-line_y_min)/(line_y_max-line_y_min);
            line_y_position.push_back(line_y_tmp);
        }

                    

        //cout<<"line y max: "<<line_y_max<<"; line y min: "<<line_y_min<<endl;
        for(int li=0;li<(int)p_g->p_l.size();li++){
            _vec<double> relative = p_g->p_v[p_g->p_l[li]->vi[0]]->loc-p_g->p_v[p_g->p_l[li]->vi[1]]->loc;
            double dist = relative.norm();
            p_g->p_l[li]->length = dist;

            double sigma_O_tmp;
            if(line_y_position[li]<sigma_O_spatial_y_boundary){
                sigma_O_tmp=sigma_O_spatial_max;
            }
            else{
                //sigma_O_tmp=(sigma_O_spatial_y_boundary-1.0)/(sigma_O_spatial_max-sigma_O_spatial_min)*line_y_position[li]+(sigma_O_spatial_max-sigma_O_spatial_y_boundary*sigma_O_spatial_min)/(sigma_O_spatial_max-sigma_O_spatial_min);
                sigma_O_tmp = (sigma_O_spatial_max-sigma_O_spatial_min)/(sigma_O_spatial_y_boundary-1.0)*line_y_position[li]+(sigma_O_spatial_min*sigma_O_spatial_y_boundary-sigma_O_spatial_max)/(sigma_O_spatial_y_boundary-1.0);
            }
            //cout<<"For line "<<li<<" the line_y_position is "<<line_y_position[li]<<"and the sigma_O_tmp is "<<sigma_O_tmp<<endl;
            _vec<double> frc;
            if(p_g->p_l[li]->IsOutermost==0){
                frc = (-1.0)*sigma_L* relative/dist;
            }
            else{
                frc = (-1.0)*sigma_O_tmp* relative/dist;
            }
            p_g->p_v[p_g->p_l[li]->vi[0]]->frc_edge += frc;
            p_g->p_v[p_g->p_l[li]->vi[1]]->frc_edge -= frc;
            p_g->p_l[li]->edgeForce = frc.norm();
        }

    }
    
    else{
        cout<<"Fatal Error: No mechanics mode is selected ! Current mechanics mode is "<<mechanics_mode<<endl;
        cout<<"End of simulation"<<endl;
        exit(1);
    }
}

//requires the assumption that vertex indices inside a cell should be arranged in anticlockwise direction
void cell_elastic_force(Organ *p_g){
    if(mechanics_mode == "simple"||mechanics_mode=="L_std"||mechanics_mode=="sigma_O_spatial"){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->area = cell_geo::cell_area(p_g,ci);
            Cell *cp = p_g->p_c[ci];

            //calculate delta S/ delta x_i
            for (int j = 0; j < (int)cp->vi.size(); j++) {
            Vertex *vp[3];
            vp[1] = p_g->p_v[cp->vi[j]];

            if (j == 0) {
                vp[0] = p_g->p_v[cp->vi[cp->vi.size() - 1]];
            }
            else {
                vp[0] = p_g->p_v[cp->vi[j - 1]];
            }

            if (j == (int)cp->vi.size() - 1) {
                vp[2] = p_g->p_v[cp->vi[0]];
            }
            else {
                vp[2] = p_g->p_v[cp->vi[j + 1]];
            }

            _vec<double> s_grad = _vec<double>(0.0, 0.0, 0.0);
            s_grad.x =  (vp[2]->loc.y - vp[0]->loc.y);
            s_grad.y =  (vp[0]->loc.x - vp[2]->loc.x);

            _vec<double> frc_tmp = (-1.0) * kappa_S * (p_g->p_c[ci]->area - S_std) * s_grad;
            vp[1]->frc_area += frc_tmp;
            }   
        }
    }

   else if(mechanics_mode=="S_std_spatial"){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->area = cell_geo::cell_area(p_g,ci);
            Cell *cp = p_g->p_c[ci];

            //calculate delta S/ delta x_i
            for (int j = 0; j < (int)cp->vi.size(); j++) {
            Vertex *vp[3];
            vp[1] = p_g->p_v[cp->vi[j]];

            if (j == 0) {
                vp[0] = p_g->p_v[cp->vi[cp->vi.size() - 1]];
            }
            else {
                vp[0] = p_g->p_v[cp->vi[j - 1]];
            }

            if (j == (int)cp->vi.size() - 1) {
                vp[2] = p_g->p_v[cp->vi[0]];
            }
            else {
                vp[2] = p_g->p_v[cp->vi[j + 1]];
            }

            _vec<double> s_grad = _vec<double>(0.0, 0.0, 0.0);
            s_grad.x =  (vp[2]->loc.y - vp[0]->loc.y);
            s_grad.y =  (vp[0]->loc.x - vp[2]->loc.x);

            _vec<double> frc_tmp = (-1.0) * kappa_S * (p_g->p_c[ci]->area - p_g->p_c[ci]->S_std) * s_grad;
            vp[1]->frc_area += frc_tmp;
            }   
        }
   }
   else{
        cout<<"Fatal Error: No mechanics mode is selected ! Current mechanics mode is "<<mechanics_mode<<endl;
        cout<<"End of simulation"<<endl;
        exit(1);
   } 
}

void force_reset(Organ *p_g){
    //force reset to zero
    for(int vi=0; vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->frc_edge = _vec<double>(0.0,0.0,0.0);
        p_g->p_v[vi]->frc_area = _vec<double>(0.0,0.0,0.0);
    }
}

void calcForceMotion(Organ* p_g){
    
    force::line_elastic_force(p_g);
    force::cell_elastic_force(p_g);
    
    //the elastic forces drive vertices to move
    for(int vi=0; vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->loc +=(p_g->p_v[vi]->frc_area+p_g->p_v[vi]->frc_edge)*delta_time;
    }
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->area>10.0 || p_g->p_c[ci]->area<0.1){
                //std::cout<<"The area of cell "<<ci<<" is abnormal: "<<p_g->p_c[ci]->area<<std::endl;
            }
    }
    
    force::force_reset(p_g);
    
}

void forceShapeInitiation(Organ *p_g, int initiation_time){
    for(int pre_step=0;pre_step<initiation_time;pre_step++){
        force::calcForceMotion(p_g);
    }
}

void S_std_spatial_geo(Organ* p_g){
    //assume that organ_geo::organ_max_min_x_y is performed before; y_min_cell and y_max_cell is already known in organ
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
        double y_cell_relative = (p_g->p_c[ci]->center.y-p_g->y_min_cell)/(p_g->y_max_cell-p_g->y_min_cell);
        p_g->p_c[ci]->S_std = (S_std_apex-S_std_base)*y_cell_relative + S_std_base;
        //std::cout<<"S_std for "<<ci<<" : "<<p_g->p_c[ci]->S_std<<std::endl;
    }
}

}