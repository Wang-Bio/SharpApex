/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#define _USE_MATH_DEFINES
#include "../include/division.h"
#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"


//const double M_PI =3.1415926;

//biased angles division parameters


namespace division_frequency{

    void no_control(Organ* p_g){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->frequency_modifier=1.0;
        }
    }

    void balance_control(Organ* p_g){
        //for temporal pattern control
        double BF_tmp = BF;
        F_modifier = BF_tmp + (BF_tmp-1)*(double)p_g->N_epi_cell/(double)p_g->N_inner_cell;
        
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->IsEpidermal==1){
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*F_modifier;
            }
        }
        std::cout<<"BF is "<< BF <<" and current F_modifier is "<<F_modifier<<std::endl;
    }

    void area_control(Organ* p_g){
        
        //sort area for each cell from the smallest to largest

        //preparation for sorting
        std::vector<double>cell_area_sort;
        std::vector<double>cell_area_rank;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(mechanics_mode=="S_std_spatial"){
                cell_area_sort.push_back(p_g->p_c[ci]->area/p_g->p_c[ci]->S_std);
            }
            else{
                cell_area_sort.push_back(p_g->p_c[ci]->area);
            }
            cell_area_rank.push_back(ci);
        }
        
        
        //get the area rank for each cell
        for(int ci=0;ci<(int)p_g->p_c.size()-1;ci++){
            for(int cj=0;cj<(int)p_g->p_c.size()-ci-1;cj++){
                if(cell_area_sort[cj]>cell_area_sort[cj+1]){
                    //change the elements
                    double sort_temp = cell_area_sort[cj+1];
                    cell_area_sort[cj+1] = cell_area_sort[cj];
                    cell_area_sort[cj] = sort_temp;
                    //change the rank
                    double rank_temp = cell_area_rank[cj+1];
                    cell_area_rank[cj+1] = cell_area_rank[cj];
                    cell_area_rank[cj] = rank_temp;
                }
            }
        }

        for(int i=0; i<(int)cell_area_sort.size();i++){
            p_g->p_c[cell_area_rank[i]]->area_rank=i;
        }

        //assign division frequency modifier to each cell : smallest 30% can not divide; modifier = 5 * (rank - 30%);
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double ci_rank_temp = p_g->p_c[ci]->area_rank/(double)p_g->p_c.size();
            if(ci_rank_temp <=area_control_lower_limit){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->area_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/(area_control_slope*(ci_rank_temp-area_control_lower_limit));
                p_g->p_c[ci]->area_modifier = 1.0/(area_control_slope*(ci_rank_temp-area_control_lower_limit));
            }
            //std::cout<<"p_g->p_c[ci]->area_modifier "<<p_g->p_c[ci]->area_modifier<<endl;
        }

        //for debug
        /*
        std::cout<<"area_control_lower_limit "<<area_control_lower_limit<<endl;
        std::cout<<"area_control_slope "<<area_control_slope<<endl;
        std::cout<<"ci ; cell area ; area_rank; area_modifier; division_frequency"<<std::endl;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            std::cout<<ci<<" ; "<<p_g->p_c[ci]->area<<" ; "<<p_g->p_c[ci]->area_rank<<" ; "<<p_g->p_c[ci]->area_modifier<<" ; "<<1/p_g->p_c[ci]->frequency_modifier<<std::endl;
        } 
        */
        
    }
    
    double calcGaussian(double x, double mu, double sigma){
        return exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
    }
    
    void Gaussian_control(Organ *p_g, double mu,double sigma){
        organ_geo::organ_center(p_g);
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }


        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,mu,sigma);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/calcGaussian(y_relative,mu,sigma);
            }
        }
    }
    
    void temporal_Gaussian_linear_mu(Organ *p_g, double mu_initial, double mu_terminal){
        organ_geo::organ_center(p_g);

        double time_by_cell_number = ((double)p_g->p_c.size()-(double)p_g->initial_cell_number)/(end_cell_number-(double)p_g->initial_cell_number);
        std::cout<<"current time by cell number is: "<<time_by_cell_number<<endl;
        double mu_a = (mu_terminal-mu_initial)/((double)end_cell_number-(double)p_g->initial_cell_number);
        double mu_b = mu_initial;

        std::cout<<"mu_a "<<mu_a<<" , mu_b "<<mu_b<<endl;

        double gau_mu_tmp=mu_a*(p_g->p_c.size()-(double)p_g->initial_cell_number)+mu_b;
        double gau_sigma_tmp=gau_sigma;
        std::cout<<"gau_mu_tmp: "<<gau_mu_tmp<<"; gau_sigma_tmp: "<<gau_sigma_tmp<<endl;
        
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,gau_mu_tmp,gau_sigma_tmp);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/calcGaussian(y_relative,gau_mu_tmp,gau_sigma_tmp);
            }
        }
    }

    void temporal_Gaussian_linear_sigma(Organ *p_g, double sigma_initial, double sigma_terminal){
        organ_geo::organ_center(p_g);

        double time_by_cell_number = ((double)p_g->p_c.size()-(double)p_g->initial_cell_number)/(end_cell_number-(double)p_g->initial_cell_number);
        std::cout<<"current time by cell number is: "<<time_by_cell_number<<endl;
        double sigma_a = (sigma_terminal-sigma_initial)/((double)end_cell_number-(double)p_g->initial_cell_number);
        double sigma_b = sigma_initial;

        double gau_sigma_tmp=sigma_a*(p_g->p_c.size()-(double)p_g->initial_cell_number)+sigma_b;
        double gau_mu_tmp=gau_mu;
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;

        std::cout<<"gau_mu_tmp: "<<gau_mu_tmp<<"; gau_sigma_tmp: "<<gau_sigma_tmp<<endl;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,gau_mu_tmp,gau_sigma_tmp);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/calcGaussian(y_relative,gau_mu_tmp,gau_sigma_tmp);
            }
        }
    }

    void temporal_even_to_Gaussian(Organ* p_g){
        double relative_time_cell_number =0.0;
        relative_time_cell_number = ((double)p_g->p_c.size()-61.0)/((double)end_cell_number-61.0);
        std::cout<<"Relative time cell number is "<<relative_time_cell_number<<endl;


        organ_geo::organ_center(p_g);
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/((1-relative_time_cell_number)+relative_time_cell_number*100);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/((1-relative_time_cell_number)+relative_time_cell_number*100);
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/((1-relative_time_cell_number)+relative_time_cell_number*calcGaussian(y_relative,gau_mu,gau_sigma));
                p_g->p_c[ci]->Gaussian_modifier = 1.0/((1-relative_time_cell_number)+relative_time_cell_number*calcGaussian(y_relative,gau_mu,gau_sigma));
            }
            //std::cout<<ci<<" "<<y_relative<<endl;
            //std::cout<<ci<<"Gaussian "<<calcGaussian(y_relative,gau_mu,gau_sigma)<<endl;
            //std::cout<<ci<<" "<<p_g->p_c[ci]->Gaussian_modifier<<endl;
        }

    }
    
    void Gaussian_xy_control(Organ* p_g, double mu_x, double sigma_x, double mu_y, double sigma_y){
        organ_geo::organ_center(p_g);
        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        double x_min=p_g->p_c[0]->center.x, x_max=p_g->p_c[0]->center.x;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }

            if(x_min>p_g->p_c[ci]->center.x){
                x_min = p_g->p_c[ci]->center.x;
            }
            
            if(x_max<p_g->p_c[ci]->center.x){
                x_max = p_g->p_c[ci]->center.x;
            }
        }

        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            double x_relative = (p_g->p_c[ci]->center.x-x_min)/(x_max-x_min);

            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,mu_y,sigma_y)*1.0/calcGaussian(y_relative,mu_x,sigma_x);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/(calcGaussian(y_relative,mu_y,sigma_y)*calcGaussian(x_relative,mu_x,sigma_x));
            }
        }
    }

    void uniform_to_Gaussian_control(Organ* p_g, double mu, double sigma, double uniform_to_Gaussian_transition_time){
        organ_geo::organ_center(p_g);

        double time_by_cell_number = ((double)p_g->p_c.size()-(double)p_g->initial_cell_number)/(end_cell_number-(double)p_g->initial_cell_number);
        std::cout<<"current time by cell number is: "<<time_by_cell_number<<endl;

        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }

        //f(y,t)=(1-at)+at*exp(-(y-\mu)^2/(2\sigma)^2), if at<1; at*exp(-(y-\mu)^2/(2\sigma)^2), if at>=1.
        if(time_by_cell_number*uniform_to_Gaussian_transition_time<1.0){
            for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                double frequency_modifier_tmp = 1.0/((1.0-time_by_cell_number*uniform_to_Gaussian_transition_time)+time_by_cell_number*uniform_to_Gaussian_transition_time*calcGaussian(y_relative,mu,sigma));
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*frequency_modifier_tmp;
                p_g->p_c[ci]->Gaussian_modifier = frequency_modifier_tmp;
            }
            }
        }
        else{
            for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_relative = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
            if(y_relative<0.01){
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*100;
                p_g->p_c[ci]->Gaussian_modifier = 100;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = p_g->p_c[ci]->frequency_modifier*1.0/calcGaussian(y_relative,mu,sigma);
                p_g->p_c[ci]->Gaussian_modifier = 1.0/calcGaussian(y_relative,mu,sigma);
            }
        }

        }
    }

    void biregion_frequency_position(Organ* p_g, double y_boundary_tmp, double relative_frequency_tmp){
        std::cout<<"y_boundary: "<<y_boundary_tmp<<"; relative_frequency: "<<relative_frequency_tmp<<endl;
        organ_geo::organ_center(p_g);
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;
        std::cout<<"y_max_tmp "<<y_max_tmp<<endl;
        std::cout<<"y_min_tmp "<<y_min_tmp<<endl;

        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
            //std::cout<<"For cell "<<ci<<" y_relative_tmp "<<y_relative_tmp<<endl;
            //if cell belongs to apical region, it will divide in a modified speed
            if(y_relative_tmp>y_boundary_tmp){
                //std::cout<<"belongs to the apical region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*1.0;
                p_g->p_c[ci]->tag = 1;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;

            }

            //else cell belongs to basal region, it will divide in normal speed
            else{
                //std::cout<<"belongs to the basal region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier/relative_frequency_tmp;
                p_g->p_c[ci]->tag = 0;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;
            }
        }

    }

    void biregion_frequency_identity(Organ* p_g, double relative_frequency_tmp){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->tag==1){
                //cell ci is originated from apical cells
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*relative_frequency_tmp;
            }
            else{
                //cell ci is originated from basal cells
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*1.0;
            }
        }
    }

    void biregion_frequency_angles_position(Organ* p_g, double y_boundary_tmp, double relative_frequency_tmp){
        std::cout<<"y_boundary: "<<y_boundary_tmp<<"; relative_frequency: "<<relative_frequency_tmp<<endl;
        organ_geo::organ_center(p_g);
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;
        std::cout<<"y_max_tmp "<<y_max_tmp<<endl;
        std::cout<<"y_min_tmp "<<y_min_tmp<<endl;

        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
            //std::cout<<"For cell "<<ci<<" y_relative_tmp "<<y_relative_tmp<<endl;
            //if cell belongs to apical region, it will divide in a modified speed
            if(y_relative_tmp>y_boundary_tmp){
                //std::cout<<"belongs to the apical region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*relative_frequency_tmp;
                p_g->p_c[ci]->tag = 1;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;

            }

            //else cell belongs to basal region, it will divide in normal speed
            else{
                //std::cout<<"belongs to the basal region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*1.0;
                p_g->p_c[ci]->tag = 0;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;
            }
        }

    }

    void tag_control(Organ *p_g){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->tag==1){
                p_g->p_c[ci]->frequency_modifier=1;
            }
        }
    }
    
    //BF_3 = (\Delta A / A_0) * ( P_0 / \Delta P)^2
    
    void balance_control_BF_3(Organ* p_g){
        /*
        int averaged_steps=20;
        
        if(p_g->step>=averaged_steps*T_division_check-1){
            //calcualte \Delta A, \Delta P
            double A_0 = p_g->geo[p_g->geo.size()-averaged_steps]->area;
            double P_0 = p_g->geo[p_g->geo.size()-averaged_steps]->perimeter;
            double delta_A = (p_g->geo[p_g->geo.size()-1]->area-A_0)/averaged_steps;
            double delta_P = (p_g->geo[p_g->geo.size()-1]->perimeter-P_0)/averaged_steps;
            double BF_3 = delta_A/A_0 * (P_0/delta_P) * (P_0/delta_P);
            
            double F_epi = 10/standard_cell_period_length;
            double F_in = BF_3*F_epi*F_epi*(1+p_g->N_epi_cell/p_g->N_inner_cell)-F_epi*p_g->N_epi_cell/p_g->N_inner_cell;
            double ratio_in_epi = F_in/F_epi;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                if(p_g->p_c[ci]->IsEpidermal==1){
                    p_g->p_c[ci]->frequency_modifier=1.0;
                }
                else{
                    p_g->p_c[ci]->frequency_modifier=1.0/ratio_in_epi;
                }
            }
            std::cout<<"For BF3 control, we have A_0 "<<A_0<<" P_0 "<<P_0<<" delta A "<<delta_A<<" delta P "<<delta_P<<" BF_3 = "<<BF_3<<" , F_in = "<<F_in<<" , ratio_in_epi = "<<ratio_in_epi<<endl;
        }
        else{
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                p_g->p_c[ci]->frequency_modifier=1.0;
            }
        }
        */
    }

    //BF_1 = (\Delta A / A_0) * ( P_0 / \Delta P)
    void balance_control_BF_1(Organ* p_g){
        /*
        int averaged_steps=20;
        if(p_g->step>=averaged_steps*T_division_check-1){
            //calcualte \Delta A, \Delta P
            double A_0 = p_g->geo[p_g->geo.size()-averaged_steps]->area;
            double P_0 = p_g->geo[p_g->geo.size()-averaged_steps]->perimeter;
            double delta_A = (p_g->geo[p_g->geo.size()-1]->area-A_0)/averaged_steps;
            double delta_P = (p_g->geo[p_g->geo.size()-1]->perimeter-P_0)/averaged_steps;
            double BF_1 = delta_A/A_0 * (P_0/delta_P);

            double ratio_in_epi = BF_1 + (BF_1-1)*p_g->N_epi_cell/p_g->N_inner_cell;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                if(p_g->p_c[ci]->IsEpidermal==1){
                    p_g->p_c[ci]->frequency_modifier=1.0;
                }
                else{
                    p_g->p_c[ci]->frequency_modifier=1.0/ratio_in_epi;
                }
            }
            std::cout<<"For BF1 control, we have A_0 "<<A_0<<" P_0 "<<P_0<<" delta A "<<delta_A<<" delta P "<<delta_P<<" BF_3 = "<<BF_1<<" , ratio_in_epi = "<<ratio_in_epi<<endl;

        }
        else{
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                p_g->p_c[ci]->frequency_modifier=1.0;
            }
        }
        */
    }

    void arrest_front(Organ* p_g, double y_arrest_front, double t_arrest_front, double k_arrest_front){
        organ_geo::organ_center(p_g);
        int current_cell_number = (int)p_g->p_c.size();
        std::cout<<"current cell number is: "<<current_cell_number<<endl;

        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }

        double arrest_front_position_tmp;
        if((current_cell_number-p_g->initial_cell_number)<t_arrest_front){
            arrest_front_position_tmp = y_arrest_front;
        }
        else{
            arrest_front_position_tmp = y_arrest_front+k_arrest_front*(current_cell_number-p_g->initial_cell_number-t_arrest_front);
        }

        std::cout<<"Current organ length is: "<<organ_geo::organ_length(p_g)<<endl;
        std::cout<<"Current arrest front position is : "<<arrest_front_position_tmp<<endl;
        if((current_cell_number-p_g->initial_cell_number)<t_arrest_front){
            for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_absolute = (p_g->p_c[ci]->center.y-y_min);
                if(y_absolute<y_arrest_front){
                    p_g->p_c[ci]->frequency_modifier = 1*p_g->p_c[ci]->frequency_modifier;
                }
                else{
                    p_g->p_c[ci]->frequency_modifier = 100000000*p_g->p_c[ci]->frequency_modifier;
                }
                //std::cout<<"Cell "<<ci<<" position is "<<y_absolute<<" and its frequency modifier is "<<p_g->p_c[ci]->frequency_modifier<<std::endl;
            }
        }
        else{
            //std::cout<<"Now we are entering the arrest front later stage"<<endl;
            for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_absolute = (p_g->p_c[ci]->center.y-y_min);
                if(y_absolute<arrest_front_position_tmp){
                    p_g->p_c[ci]->frequency_modifier = 1*p_g->p_c[ci]->frequency_modifier;
                }
                else{
                    p_g->p_c[ci]->frequency_modifier = 100000000*p_g->p_c[ci]->frequency_modifier;
                }
                std::cout<<"Cell "<<ci<<" position is "<<y_absolute<<" and its frequency modifier is "<<p_g->p_c[ci]->frequency_modifier<<std::endl;
            }
        }

        if(arrest_front_position_tmp<=0){
            std::cout<<"There is no more cell division. The simulation will go to end"<<endl;
            force::forceShapeInitiation(p_g,200000);
            output::VTK(p_g);
            return;
        }
    }

    void arrest_front_biregion(Organ* p_g, double arrest_front_biregion_x, double arrest_front_biregion_y,double arrest_front_biregion_t, double arrest_front_biregion_k){
        organ_geo::organ_center(p_g);
        int current_cell_number = (int)p_g->p_c.size();
        std::cout<<"current cell number is: "<<current_cell_number<<endl;

        double y_min=p_g->p_c[0]->center.y, y_max=p_g->p_c[0]->center.y;
        double x_min=p_g->p_c[0]->center.x, x_max=p_g->p_c[0]->center.x;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(y_min>p_g->p_c[ci]->center.y){
                y_min = p_g->p_c[ci]->center.y;
            }
            
            if(y_max<p_g->p_c[ci]->center.y){
                y_max = p_g->p_c[ci]->center.y;
            }
        }

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(x_min>p_g->p_c[ci]->center.x){
                x_min = p_g->p_c[ci]->center.x;
            }
            
            if(x_max<p_g->p_c[ci]->center.x){
                x_max = p_g->p_c[ci]->center.x;
            }
        }

        double arrest_front_biregion_tmp_x;
        double arrest_front_biregion_stage_2_t = arrest_front_biregion_x/arrest_front_biregion_k;
        if((current_cell_number-p_g->initial_cell_number)<arrest_front_biregion_t){
            arrest_front_biregion_tmp_x = 0;
        }
        else if((current_cell_number-p_g->initial_cell_number)<arrest_front_biregion_t+arrest_front_biregion_stage_2_t){
            arrest_front_biregion_tmp_x = arrest_front_biregion_k*(current_cell_number-p_g->initial_cell_number-arrest_front_biregion_t);
        }
        else{
            arrest_front_biregion_tmp_x = arrest_front_biregion_x;
        }
        std::cout<<"Current organ length is: "<<organ_geo::organ_length(p_g)<<endl;
        std::cout<<"Current arrest front position is : "<<arrest_front_biregion_y<<endl;

        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            double y_absolute = (p_g->p_c[ci]->center.y-y_min);
            double x_relative = (p_g->p_c[ci]->center.x-x_min)/(x_max-x_min);
            if(y_absolute<arrest_front_biregion_y){
                p_g->p_c[ci]->frequency_modifier = 1*p_g->p_c[ci]->frequency_modifier;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = 100000000*p_g->p_c[ci]->frequency_modifier;
            }

            if(x_relative<0.5-arrest_front_biregion_tmp_x){
                p_g->p_c[ci]->frequency_modifier = 1*p_g->p_c[ci]->frequency_modifier;
            }
            else if(x_relative>0.5+arrest_front_biregion_tmp_x){
                p_g->p_c[ci]->frequency_modifier = 1*p_g->p_c[ci]->frequency_modifier;
            }
            else{
                p_g->p_c[ci]->frequency_modifier = 100000000*p_g->p_c[ci]->frequency_modifier;
            }
            //std::cout<<"Cell "<<ci<<" position is "<<y_absolute<<" and its frequency modifier is "<<p_g->p_c[ci]->frequency_modifier<<std::endl;
        }
    }

    void meristem_position_relative_constant(Organ* p_g, double meristem_position_relative_constant_y){
        std::cout<<"meristem_position_relative_constant_y: "<<meristem_position_relative_constant_y;
        organ_geo::organ_center(p_g);
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;
        std::cout<<"; y_max_tmp "<<y_max_tmp;
        std::cout<<"; y_min_tmp "<<y_min_tmp<<endl;
    
        for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
            double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
            //std::cout<<"For cell "<<ci<<" y_relative_tmp "<<y_relative_tmp<<endl;
            //if cell belongs to apical region, it will divide in a modified speed
            if(y_relative_tmp>meristem_position_relative_constant_y){
                //std::cout<<"belongs to the apical region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*1000000000;
                p_g->p_c[ci]->tag = 1;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;

            }

            //else cell belongs to basal region, it will divide in normal speed
            else{
                //std::cout<<"belongs to the basal region"<<endl;
                //std::cout<<"frequecy_modifier_before: "<<p_g->p_c[ci]->frequency_modifier<<endl;
                p_g->p_c[ci]->frequency_modifier=p_g->p_c[ci]->frequency_modifier*1;
                p_g->p_c[ci]->tag = 0;
                //std::cout<<"frequecy_modifier_after: "<<p_g->p_c[ci]->frequency_modifier<<endl;
            }
        }
    }

}

namespace division_direction{
    _vec<double> angles(Organ* p_g, int ci){
        organ_geo::organ_center(p_g);
        _vec<double> angles_tmp;
        if(p_g->p_c[ci]->IsEpidermal==0){
            if(in_division_direction=="random"){
                division_direction::random(p_g, ci);
            }
            else if(in_division_direction=="constant_0"){
                division_direction::constant_0(p_g,ci);
            }
            else if(in_division_direction=="constant_90"){
                division_direction::constant_90(p_g,ci);
            }
            else if(in_division_direction=="yin_distribution"){
                division_direction::yin_distribution(p_g,ci);
            }
            else if(in_division_direction=="mochizuki_bias"){
                division_direction::mochizuki_bias(p_g,ci,mochizuki_bias_beta,mochizuki_bias_phi);
            }
            else if(in_division_direction=="mochizuki_bias_asymmetrical"){
                division_direction::mochizuki_bias_asymmetrical(p_g,ci,mochizuki_bias_beta_left,mochizuki_bias_beta_right,mochizuki_bias_phi_left,mochizuki_bias_phi_right);
            }
            else if(in_division_direction=="mochizuki_bias_apical_basal"){
                division_direction::mochizuki_bias_apical_basal(p_g,ci,mochizuki_bias_beta_apical,mochizuki_bias_beta_basal, mochizuki_bias_phi_apical, mochizuki_bias_phi_basal, angle_bias_y_boundary);
            }
            else if(in_division_direction=="temporal_angle_bias"){
                division_direction::temporal_angle_bias(p_g,ci,temporal_bias_beta_initial, temporal_bias_beta_terminal, temporal_bias_phi_initial, temporal_bias_phi_terminal);
            }
            else if(in_division_direction=="Gaussian_bias"){
                division_direction::Gaussian_bias(p_g,ci,Gaussian_bias_beta,Gaussian_bias_phi);
            }
            else if(in_division_direction=="biregion_angles_position"){
                division_direction::biregion_angles_position(p_g, ci, biregion_angles_position_apical_gaussian_beta, biregion_angles_position_basal_gaussian_beta, biregion_angles_position_apical_gaussian_phi, biregion_angles_position_basal_gaussian_phi, biregion_angles_position_y_boundary);
            }
            else if(in_division_direction=="biregion_angles_gradual"){
                division_direction::biregion_angles_gradual(p_g, ci, biregion_angles_gradual_apical_gaussian_beta, biregion_angles_gradual_basal_gaussian_beta, biregion_angles_gradual_apical_gaussian_phi, biregion_angles_gradual_basal_gaussian_phi, biregion_angles_gradual_y1_boundary, biregion_angles_gradual_y2_boundary);
            }
            else if(in_division_direction=="biregion_angles_identity"){
                division_direction::biregion_angles_identity(p_g, ci, biregion_angles_identity_apical_gaussian_beta, biregion_angles_identity_basal_gaussian_beta, biregion_angles_identity_apical_gaussian_phi, biregion_angles_identity_basal_gaussian_phi, biregion_angles_identity_y_boundary);
            }
            else if(in_division_direction=="biregion_frequency_angles_position"){
                division_direction::biregion_frequency_angles_position(p_g, ci, biregion_frequency_angles_position_apical_gaussian_beta, biregion_frequency_angles_position_basal_gaussian_beta, biregion_frequency_angles_position_apical_gaussian_phi, biregion_frequency_angles_position_basal_gaussian_phi, biregion_frequency_angles_position_y_boundary);
            }
            else if(in_division_direction=="Gaussian_bias_continuous_Gaussian"){
                division_direction::Gaussian_bias_continuous_Gaussian(p_g, ci, Gaussian_bias_beta_A, Gaussian_bias_beta_mu, Gaussian_bias_beta_sigma, Gaussian_bias_phi);
            }
            else if(in_division_direction=="temporal_biregion_angles"){
                division_direction::temporal_biregion_angles(p_g,ci,y_boundary_initial,y_boundary_change,y_boundary_terminal,temporal_biregion_angles_apical_gaussian_beta,temporal_biregion_angles_apical_gaussian_phi,temporal_biregion_angles_basal_gaussian_beta,temporal_biregion_angles_basal_gaussian_phi);
            }
            else if(in_division_direction=="temporal_angle_bias_Gaussian"){
                division_direction::temporal_angle_bias_Gaussian(p_g,ci,bias_sigma_initial,bias_sigma_terminal,bias_mu_initial,bias_mu_terminal,bias_t_init,bias_t_term);
            }
            else{
                std::cout<<"Fatal Error: no inner cell division direction control selected"<<endl;
                exit(-1);
            }
        }
        else if(p_g->p_c[ci]->IsEpidermal==1){
            if(epi_division_direction=="random"){
                division_direction::random(p_g,ci);
            }
            else if(epi_division_direction=="anticlinal"){
                division_direction::epi_anticlinal(p_g,ci);
                //std::cout<<"anticlinal division of cell "<<ci<<endl;
            }
            else if(epi_division_direction=="periclinal"){
                division_direction::epi_periclinal(p_g,ci);
            }
            else if(epi_division_direction=="constant_0"){
                division_direction::constant_0(p_g,ci);
            }
            else if(epi_division_direction=="constant_90"){
                division_direction::constant_90(p_g,ci);
            }
            else{
                std::cout<<"Fatal Error: no epidermal cell division direction control selected"<<endl;
                exit(-1);
            }
        }
        return p_g->p_c[ci]->axis;
    }
    
    //direction: random, anticlinal,periclinal, constant_0 (horizontal), constant_90 (vertical), and yin_distribution, mochizuki bias
    _vec<double> random(Organ* p_g, int ci){
        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution<> rand_axis(0.0,M_PI);
        double axisTheta=rand_axis(mt);
        p_g->p_c[ci]->axisTheta = axisTheta;
        p_g->p_c[ci]->axis = _vec<double>(cos(axisTheta),sin(axisTheta),0.0);

        return p_g->p_c[ci]->axis;
    }

    _vec<double> epi_anticlinal(Organ* p_g, int ci){
        if(p_g->p_c[ci]->IsEpidermal==0){
            std::cout<<"Fatal Error: inner cell is using epidermal anticlinal direction control"<<endl;
            exit(-1);
        }
        _vec<double> center_outermostEdge = (p_g->p_v[p_g->p_c[ci]->surfaceVertex[0]]->loc + p_g->p_v[p_g->p_c[ci]->surfaceVertex[1]]->loc)/2.0;
        p_g->p_c[ci]->axis = center_outermostEdge - p_g->p_c[ci]->center;
        p_g->p_c[ci]->axisTheta = atan(p_g->p_c[ci]->axis.y/p_g->p_c[ci]->axis.x);
        return p_g->p_c[ci]->axis;
    }

    _vec<double> epi_periclinal(Organ* p_g, int ci){
        if(p_g->p_c[ci]->IsEpidermal==0){
            std::cout<<"Fatal Error: inner cell is using epidermal periclinal direction control"<<endl;
            exit(-1);
        }
        p_g->p_c[ci]->axis = p_g->p_v[p_g->p_c[ci]->surfaceVertex[0]]->loc - p_g->p_v[p_g->p_c[ci]->surfaceVertex[1]]->loc;
        p_g->p_c[ci]->axisTheta = atan(p_g->p_c[ci]->axis.y/p_g->p_c[ci]->axis.x);
        return p_g->p_c[ci]->axis;
    }

    _vec<double> constant_0(Organ* p_g,int ci){
        p_g->p_c[ci]->axis = _vec<double>(0.0,1.0,0.0);
        p_g->p_c[ci]->axisTheta = 0.0;
        return p_g->p_c[ci]->axis;
    }

    _vec<double> constant_90(Organ* p_g, int ci){
        p_g->p_c[ci]->axis = _vec<double>(1.0,0.0,0.0);
        p_g->p_c[ci]->axisTheta = M_PI/2.0;
        return p_g->p_c[ci]->axis;
    }

    _vec<double> yin_distribution(Organ* p_g, int ci){
        std::uniform_real_distribution<> rand_axis_10_170(M_PI/18.0,17.0*M_PI/18.0);
        std::uniform_real_distribution<> rand_axis_0_10(0.0,M_PI/18.0);
        std::uniform_real_distribution<> rand_axis_170_180(17.0*M_PI/18.0,M_PI);
        std::random_device rnd;
        std::mt19937 mt(rnd());
        std::uniform_real_distribution<> rand_axis_0_180_or_10_170(0.0,1.0);
        double axis_inner_judge = rand_axis_0_180_or_10_170(mt);
        //std::cout<<"The axis judge is "<<axis_inner_judge<<std::endl;
        if(axis_inner_judge<0.17){
        p_g->p_c[ci]->axisTheta = rand_axis_0_10(mt);
        p_g->p_c[ci]->axis = _vec<double>(cos(p_g->p_c[ci]->axisTheta),sin(p_g->p_c[ci]->axisTheta),0.0);
        //std::cout<<"The axis theta belongs to the 0_10 group and the value is "<<p_g->p_c[ci]->axisTheta <<std::endl;
        }
        else if(axis_inner_judge>0.83){
        p_g->p_c[ci]->axisTheta = rand_axis_170_180(mt);
        p_g->p_c[ci]->axis = _vec<double>(cos(p_g->p_c[ci]->axisTheta),sin(p_g->p_c[ci]->axisTheta),0.0);
        //std::cout<<"The axis theta belongs to the 170_180 group and the value is "<<p_g->p_c[ci]->axisTheta <<std::endl;
        }
        else {
        p_g->p_c[ci]->axisTheta = rand_axis_10_170(mt);
        p_g->p_c[ci]->axis = _vec<double>(cos(p_g->p_c[ci]->axisTheta),sin(p_g->p_c[ci]->axisTheta),0.0);
        //std::cout<<"The axis theta belongs to the 10_170 group and the value is "<<p_g->p_c[ci]->axisTheta <<std::endl;
        }

        return p_g->p_c[ci]->axis;
    }

    _vec<double> mochizuki_bias(Organ* p_g, int ci, double bias_beta, double bias_phi){
        double axisTheta = wangMath::mochizuki_bias_single_random_sampling(bias_beta, bias_phi);
        p_g->p_c[ci]->axisTheta = axisTheta;
        p_g->p_c[ci]->axis = _vec<double>(cos(axisTheta),sin(axisTheta),0.0);
        return p_g->p_c[ci]->axis;
    }    

    _vec<double> mochizuki_bias_asymmetrical(Organ* p_g, int ci, double bias_beta_left, double bias_beta_right, double bias_phi_left, double bias_phi_right){
        double x_max=p_g->p_c[0]->center.x;
        double x_min=p_g->p_c[0]->center.x;

        for(int i=0; i<(int)p_g->p_c.size();i++){
            if(x_max<p_g->p_c[i]->center.x){
                x_max=p_g->p_c[i]->center.x;
            }
            if(x_min>p_g->p_c[i]->center.x){
                x_min=p_g->p_c[i]->center.x;
            }
        }

        double x_relative_tmp = (p_g->p_c[ci]->center.x-x_min)/(x_max-x_min);

        if(x_relative_tmp<0.5){
            //this cell belongs to the left part
            mochizuki_bias(p_g,ci,bias_beta_left,bias_phi_left);
        }
        else{
            //this cell belongs to the right part
            mochizuki_bias(p_g,ci,bias_beta_right,bias_phi_right);
        }

        return p_g->p_c[ci]->axis;

    }

    _vec<double> mochizuki_bias_apical_basal(Organ* p_g, int ci, double bias_beta_apical, double bias_beta_basal, double bias_phi_apical, double bias_phi_basal, double bias_y_boundary){
        double y_max=p_g->p_c[0]->center.y;
        double y_min=p_g->p_c[0]->center.y;

        for(int i=0; i<(int)p_g->p_c.size();i++){
            if(y_max<p_g->p_c[i]->center.y){
                y_max=p_g->p_c[i]->center.y;
            }
            if(y_min>p_g->p_c[i]->center.y){
                y_min=p_g->p_c[i]->center.y;
            }
        }

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
        
        if(y_relative_tmp>bias_y_boundary){
            mochizuki_bias(p_g,ci,bias_beta_apical,bias_phi_apical);
        }
        else{
            mochizuki_bias(p_g,ci,bias_beta_basal,bias_phi_basal);
        }

        return p_g->p_c[ci]->axis;
    }

    _vec<double> temporal_angle_bias(Organ* p_g, int ci, double bias_beta_initial, double bias_beta_terminal, double bias_phi_initial, double bias_phi_terminal){
        organ_geo::organ_center(p_g);

        double time_by_cell_number = ((double)p_g->p_c.size()-(double)p_g->initial_cell_number)/(end_cell_number-(double)p_g->initial_cell_number);
        std::cout<<"current time by cell number is: "<<time_by_cell_number<<endl;

        double beta_tmp = bias_beta_initial + (bias_beta_terminal-bias_beta_initial)*time_by_cell_number;
        double phi_tmp = bias_phi_initial + (bias_phi_terminal-bias_phi_initial)*time_by_cell_number;

        mochizuki_bias(p_g,ci,beta_tmp,phi_tmp);
        return p_g->p_c[ci]->axis;
    }

    _vec<double> temporal_angle_bias_Gaussian(Organ* p_g, int ci, double bias_sigma_initial, double bias_sigma_terminal, double bias_mu_initial, double bias_mu_terminal, int bias_t_init, int bias_t_term){
        organ_geo::organ_center(p_g);
        

        double time_by_cell_number = ((double)p_g->p_c.size()-(double)p_g->initial_cell_number);
        std::cout<<"current time by cell number is: "<<time_by_cell_number<<endl;
        double sigma_tmp=0, mu_tmp=0;
        if(bias_t_init==bias_t_term){
            if(time_by_cell_number<bias_t_init){
                sigma_tmp = bias_sigma_initial;
                mu_tmp = bias_mu_initial;
            }
            else{
                sigma_tmp = bias_sigma_terminal;
                mu_tmp = bias_mu_terminal;
            }
        }
        else{
            if(time_by_cell_number<bias_t_init){
                sigma_tmp = bias_sigma_initial;
                mu_tmp = bias_mu_initial;
            }
            else if(time_by_cell_number>bias_t_term){
                sigma_tmp = bias_sigma_terminal;
                mu_tmp = bias_mu_terminal;
            }
            else{
                sigma_tmp = (bias_sigma_initial)+time_by_cell_number*(bias_sigma_terminal-bias_sigma_initial)/(bias_t_term-bias_t_init);
                mu_tmp = (bias_mu_initial)+time_by_cell_number*(bias_mu_terminal-bias_mu_initial)/(bias_t_term-bias_t_init);
            }
        }

        Gaussian_bias(p_g,ci,sigma_tmp,mu_tmp);
        return p_g->p_c[ci]->axis;
    }



    _vec<double> Gaussian_bias(Organ* p_g, int ci, double bias_beta, double bias_phi){
        if(bias_beta==0){
           division_direction::random(p_g,ci);
        }
        else{
            random_device rnd;
            mt19937 mt(rnd());
            normal_distribution<> normal_axis(bias_phi,1/bias_beta);
            double axisTheta= normal_axis(mt);
            while(axisTheta<-M_PI||axisTheta>2*M_PI){
                axisTheta= normal_axis(mt);
            }
            if(axisTheta<0){
                axisTheta += M_PI;
            }
            if(axisTheta>M_PI){
                axisTheta -= M_PI;
            }

            p_g->p_c[ci]->axisTheta = axisTheta;
            p_g->p_c[ci]->axis = _vec<double>(cos(axisTheta),sin(axisTheta),0.0);
        }
        
        return p_g->p_c[ci]->axis;
    }

    _vec<double> biregion_angles_position(Organ* p_g, int ci, double bias_beta_apical, double bias_beta_basal, double bias_phi_apical, double bias_phi_basal, double bias_y_boundary){
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
        //std::cout<<" y_max = "<<y_max_tmp<<" ; y_min = "<<y_min_tmp<<" ; y position for cell ci = "<<p_g->p_c[ci]->center.y<<" ; y_relative for cell ci = "<<y_relative_tmp<<" ; y_boundary = "<<bias_y_boundary<<" ; apical bias phi = "<<biregion_angles_position_apical_gaussian_phi<<" ; apical bias beta = "<<biregion_angles_position_apical_gaussian_beta<<" ; basal bias phi = "<<biregion_angles_position_basal_gaussian_phi<<" ; basal bias beta = "<<biregion_angles_position_basal_gaussian_beta<<endl;  
        if(y_relative_tmp>bias_y_boundary){
            Gaussian_bias(p_g,ci,bias_beta_apical,bias_phi_apical);
            //std::cout<<"Cell "<<ci<<" belongs to apical region: (beta,phi)=("<<bias_beta_apical<<","<<bias_phi_apical<<")"<<endl;
        }
        else{
            Gaussian_bias(p_g,ci,bias_beta_basal,bias_phi_basal);
            //std::cout<<"Cell "<<ci<<" belongs to basal region: (beta,phi)=("<<bias_beta_basal<<","<<bias_phi_basal<<")"<<endl;
        }

        return p_g->p_c[ci]->axis;
    }
    _vec<double> biregion_angles_gradual(Organ* p_g, int ci, double bias_beta_apical, double bias_beta_basal, double bias_phi_apical, double bias_phi_basal, double bias_y1_boundary, double bias_y2_boundary){
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
        //std::cout<<" y_max = "<<y_max_tmp<<" ; y_min = "<<y_min_tmp<<" ; y position for cell ci = "<<p_g->p_c[ci]->center.y<<" ; y_relative for cell ci = "<<y_relative_tmp<<" ; y_boundary = "<<bias_y_boundary<<" ; apical bias phi = "<<biregion_angles_position_apical_gaussian_phi<<" ; apical bias beta = "<<biregion_angles_position_apical_gaussian_beta<<" ; basal bias phi = "<<biregion_angles_position_basal_gaussian_phi<<" ; basal bias beta = "<<biregion_angles_position_basal_gaussian_beta<<endl;  
        if(y_relative_tmp>bias_y1_boundary){
            Gaussian_bias(p_g,ci,bias_beta_apical,bias_phi_apical);
            //std::cout<<"Cell "<<ci<<" belongs to apical region: (beta,phi)=("<<bias_beta_apical<<","<<bias_phi_apical<<")"<<endl;
        }
        else if(y_relative_tmp<bias_y2_boundary){
            Gaussian_bias(p_g,ci,bias_beta_basal,bias_phi_basal);
            //std::cout<<"Cell "<<ci<<" belongs to basal region: (beta,phi)=("<<bias_beta_basal<<","<<bias_phi_basal<<")"<<endl;
        }
        else{
            double bias_beta_joint = bias_beta_basal+(bias_beta_apical-bias_beta_basal)/(bias_y1_boundary-bias_y2_boundary)*(y_relative_tmp-bias_y2_boundary);
            double bias_phi_joint = bias_phi_basal+(bias_phi_apical-bias_phi_basal)/(bias_y1_boundary-bias_y2_boundary)*(y_relative_tmp-bias_y2_boundary);
            Gaussian_bias(p_g,ci,bias_beta_joint,bias_phi_joint);
            std::cout<<"Cell "<<ci<<" belongs to joint region: (beta,phi)=("<<bias_beta_joint<<","<<bias_phi_joint<<")"<<endl;
        }

        return p_g->p_c[ci]->axis;
    }

    _vec<double> biregion_frequency_angles_position(Organ* p_g, int ci, double bias_beta_apical, double bias_beta_basal, double bias_phi_apical, double bias_phi_basal, double bias_y_boundary){
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
        //std::cout<<" y_max = "<<y_max_tmp<<" ; y_min = "<<y_min_tmp<<" ; y position for cell ci = "<<p_g->p_c[ci]->center.y<<" ; y_relative for cell ci = "<<y_relative_tmp<<" ; y_boundary = "<<bias_y_boundary<<" ; apical bias phi = "<<biregion_angles_position_apical_gaussian_phi<<" ; apical bias beta = "<<biregion_angles_position_apical_gaussian_beta<<" ; basal bias phi = "<<biregion_angles_position_basal_gaussian_phi<<" ; basal bias beta = "<<biregion_angles_position_basal_gaussian_beta<<endl;  
        if(y_relative_tmp>bias_y_boundary){
            Gaussian_bias(p_g,ci,bias_beta_apical,bias_phi_apical);
            //std::cout<<"Cell "<<ci<<" belongs to apical region: (beta,phi)=("<<bias_beta_apical<<","<<bias_phi_apical<<")"<<endl;
        }
        else{
            Gaussian_bias(p_g,ci,bias_beta_basal,bias_phi_basal);
            //std::cout<<"Cell "<<ci<<" belongs to basal region: (beta,phi)=("<<bias_beta_basal<<","<<bias_phi_basal<<")"<<endl;
        }

        return p_g->p_c[ci]->axis;
    }

    _vec<double> biregion_angles_identity(Organ* p_g, int ci, double bias_beta_apical, double bias_beta_basal, double bias_phi_apical, double bias_phi_basal, double bias_y_boundary){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->tag==1){
                //cell ci is originated from apical cells
                Gaussian_bias(p_g,ci,bias_beta_apical,bias_phi_apical);
            }
            else{
                //cell ci is originated from basal cells
                Gaussian_bias(p_g,ci,bias_beta_basal,bias_phi_basal);
            }
        }

        return p_g->p_c[ci]->axis;
    }

    _vec<double> Gaussian_bias_continuous_Gaussian(Organ* p_g, int ci, double bias_beta_A, double bias_beta_mu, double bias_beta_sigma, double bias_phi){
        double y_max=p_g->p_c[0]->center.y;
        double y_min=p_g->p_c[0]->center.y;

        for(int i=0; i<(int)p_g->p_c.size();i++){
            if(y_max<p_g->p_c[i]->center.y){
                y_max=p_g->p_c[i]->center.y;
            }
            if(y_min>p_g->p_c[i]->center.y){
                y_min=p_g->p_c[i]->center.y;
            }
        }

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min)/(y_max-y_min);
        
        double bias_beta_tmp = bias_beta_A*exp(-(y_relative_tmp-bias_beta_mu)*(y_relative_tmp-bias_beta_mu)/(2*bias_beta_sigma*bias_beta_sigma));
        Gaussian_bias(p_g,ci,bias_beta_tmp, bias_phi);

        return p_g->p_c[ci]->axis;
    }   

    _vec<double> temporal_biregion_angles(Organ* p_g, int ci, double y_boundary_initial, double y_boundary_change, double y_boundary_terminal, double temporal_biregion_angles_apical_beta, double temporal_biregion_angles_apical_phi, double temporal_biregion_angles_basal_beta,double temporal_biregion_angles_basal_phi){
        double y_max_tmp = p_g->y_max_cell;
        double y_min_tmp = p_g->y_min_cell;

        double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
        
        double relative_time = (p_g->p_c.size()-p_g->initial_cell_number)/(end_cell_number-p_g->initial_cell_number);

        double bias_y_boundary = y_boundary_initial + relative_time * y_boundary_change;
        if(bias_y_boundary<y_boundary_terminal){

        }
        else{
            bias_y_boundary=y_boundary_terminal;
        }
        
        if(y_relative_tmp>bias_y_boundary){
            Gaussian_bias(p_g,ci,temporal_biregion_angles_apical_beta,temporal_biregion_angles_apical_phi);
            //std::cout<<"Cell "<<ci<<" belongs to apical region: (beta,phi)=("<<temporal_biregion_angles_apical_beta<<","<<temporal_biregion_angles_apical_phi<<")"<<endl;
        }
        else{
            Gaussian_bias(p_g,ci,temporal_biregion_angles_basal_beta,temporal_biregion_angles_basal_phi);
            //std::cout<<"Cell "<<ci<<" belongs to basal region: (beta,phi)=("<<temporal_biregion_angles_basal_beta<<","<<temporal_biregion_angles_basal_phi<<")"<<endl;
        }

        return p_g->p_c[ci]->axis;
    }

}

namespace division{

    template<typename T>
    void findAndErase(std::vector<T> &vector, T search) {
            auto itr = std::find(vector.begin(), vector.end(), search);
            if(itr == vector.end()) {
            std::cerr << "Couldn't find " << search << std::endl;
            return;
            }
            vector.erase(itr);
        }    

    void Global(Organ *p_g){
        //division::Axis(p_g);
        int cell_division_count=0;
        int inner_cell_division_count=0;
        int epi_cell_division_count=0;
        
        //std::cout<<"Start division frequency judgement ---- ";

        if(division_control=="random_picking"){
            std::random_device rnd;
            std::mt19937 mt(rnd());
            std::uniform_int_distribution<> division_random_picking(0,p_g->p_c.size());

            int ci_dividing = division_random_picking(mt);
            if(p_g->p_c.size()>=end_cell_number){
                goto EndDivision2;
            }
            //std::cout<<"Inner cell "<<ci<<" dividing"<<std::endl;
            division::One(p_g,ci_dividing);
            std::cout<<"Cell "<<ci_dividing<<" divided"<<std::endl;
            cell_division_count =1;
            division::Record(p_g,ci_dividing);
        }
        else{
        //1. prepare for division frequency control
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->frequency_modifier = 1.0;
        }        

        //2. cell division frequency control
        if(division_control=="no"){
            division_frequency::no_control(p_g);
        }
        else if(division_control=="balance"){
            division_frequency::balance_control(p_g);
        }
        else if(division_control=="area"){
            division_frequency::area_control(p_g);
        }
        else if(division_control=="Gaussian_area"){
            division_frequency::area_control(p_g);
            division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
        }
        else if(division_control=="Gaussian_balance"){
            division_frequency::balance_control(p_g);
            division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
        }
        else if(division_control=="balance_BF_3"){
            division_frequency::balance_control_BF_3(p_g);
        }
        else if(division_control=="balance_BF_1"){
            division_frequency::balance_control_BF_1(p_g);
        }
        else if(division_control=="temporal_Gaussian_linear_mu"){
            division_frequency::temporal_Gaussian_linear_mu(p_g,temporal_Gaussian_initial_mu,temporal_Gaussian_terminal_mu);
            division_frequency::area_control(p_g);
        }
        else if(division_control=="temporal_Gaussian_linear_sigma"){
            division_frequency::temporal_Gaussian_linear_sigma(p_g, temporal_Gaussian_initial_sigma,temporal_Gaussian_terminal_sigma);
            division_frequency::area_control(p_g);
        }
        else if(division_control=="temporal_even_to_Gaussian"){
            division_frequency::temporal_even_to_Gaussian(p_g);
            division_frequency::area_control(p_g);
        }
        else if(division_control=="uniform_to_Gaussian"){
            division_frequency::uniform_to_Gaussian_control(p_g, gau_mu, gau_sigma, uniform_to_Gaussian_transition_time);
            division_frequency::area_control(p_g);
        }
        else if(division_control=="Gaussian_xy"){
            division_frequency::Gaussian_xy_control(p_g, gau_mu_x, gau_sigma_x, gau_mu_y, gau_sigma_y);
            division_frequency::area_control(p_g);
        }
        else if(division_control=="biregion_frequency_position"){
            division_frequency::area_control(p_g);
            division_frequency::biregion_frequency_position(p_g, biregion_frequency_position_y_boundary, biregion_frequency_position_relative_frequency);
        }
        else if(division_control=="biregion_frequency_identity"){
            division_frequency::area_control(p_g);
            division_frequency::biregion_frequency_identity(p_g,biregion_identity_relative_frequency);
        }
        else if(division_control=="biregion_frequency_angles_position"){
            division_frequency::area_control(p_g);
            division_frequency::biregion_frequency_angles_position(p_g, biregion_frequency_angles_position_y_boundary, biregion_frequency_angles_position_relative_frequency);
        }
        else if(division_control=="arrest_front"){
            division_frequency::area_control(p_g);
            division_frequency::arrest_front(p_g,y_arrest_front,t_arrest_front,k_arrest_front);
        }
        else if(division_control=="arrest_front_biregion"){
            division_frequency::area_control(p_g);
            division_frequency::arrest_front_biregion(p_g,arrest_front_biregion_x,arrest_front_biregion_y,arrest_front_biregion_t,arrest_front_biregion_k);
        }
        else if(division_control=="meristem_position_relative_constant"){
            division_frequency::area_control(p_g);
            division_frequency::meristem_position_relative_constant(p_g,meristem_position_relative_constant_y);
        }
        else{
            std::cout<<"Fatal Error: no division frequency control selected"<<endl;
            exit(-1);
        }

        //3. cell division determination
        //std::cout<<"Start cell division determination ---- ";
        
        for(int ci=0;ci<p_g->p_c.size();ci++){
            //3.1 check if end (cell number reached to termination condition)
            if(p_g->p_c.size()>=end_cell_number){
                std::cout<<"It is already "<<end_cell_number<<" cells."<<endl;
                goto EndDivision2;
            }
            //3.2 if the cell time is over threshold, then the cell will divide
            if(p_g->p_c[ci]->cellTime>standard_cell_period_length*p_g->p_c[ci]->frequency_modifier){
                random_device rnd;
                mt19937 mt(rnd());
                uniform_real_distribution<> division_judge(0,1);
                double division_judge_tmp = division_judge(mt);
                if(division_judge_tmp<F_apparent){
                    //3.3 cell division angles judgement 

                    //std::cout<<"Cell division direction judgement ---- "<<endl;
                    division_direction::angles(p_g,ci);
                    //3.4 cell division execution
                    division::One(p_g,ci);
                    //3.5 cell division information std::cout and record
                    cell_division_count++;
                    if(p_g->p_c[ci]->IsEpidermal==1){
                        std::cout<<"Epidermal cell "<<ci<<" divided"<<"; division frequency: "<<1.0/p_g->p_c[ci]->frequency_modifier<<"; division angle: "<<p_g->p_c[ci]->axisTheta<<endl;
                        epi_cell_division_count++;
                    }
                    else{
                        std::cout<<"Inner cell "<<ci<<" divided"<<"; division frequency: "<<1.0/p_g->p_c[ci]->frequency_modifier<<"; division angle: "<<p_g->p_c[ci]->axisTheta<<endl;
                        inner_cell_division_count++;
                    }
                    organ_geo::epidermal_identity(p_g);
                    division::Record(p_g,ci);
                }
            }
        }
        if(cell_division_count==0){
            std::cout<<"No cell divided in this step"<<endl;
        }
        else if(cell_division_count==1){
            std::cout<<"1 cell divided in this step"<<endl;
        }
        else{
            std::cout<<cell_division_count<<" cells divided in this step"<<endl;
        }
        std::cout<<"Now, we have "<<(int)p_g->p_c.size()-cell_division_count<<" + "<<cell_division_count<<" cells: inner cells = "<<p_g->N_inner_cell-inner_cell_division_count<<" + "<<inner_cell_division_count<<"; epidermal cells = "<<p_g->N_epi_cell-epi_cell_division_count<<" + "<<epi_cell_division_count<<endl;
        std::cout<<"Organ area: "<<p_g->area<<"; Organ perimeter: "<<p_g->perimeter<<endl;
        }
        EndDivision2:;
    }

    void One(Organ *p_g, int cellula_idx){
    //this function will add two new vertex, three new lines and one new cell, four cells need to be changed
    Cell *cp = p_g->p_c[cellula_idx];
    
    _vec<double> axis = p_g->p_c[cellula_idx]->axis; //from division::Axis
    
    _vec<double> center = cp->center; //from organ_geo::organ_center
    
    geo_basic::pLine division = std::make_pair(center, center+axis);
    std::vector<std::pair<int, _vec<double>>> crosspoint;
    std::vector<int> epoint_idx;
    
    //std::cout<<"Cell center is "<<center.x<<","<<center.y<<";Cell axisTheta is"<<p_g->p_c[cellula_idx]->axisTheta<<";Cell axis is "<<axis.x<<","<<axis.y<<std::endl;

    // sort lines counterclockwise, based on counterclockwise vertex 
    std::vector<int> anticlockwise_lidx;
    anticlockwise_lidx = cell_geo::cell_counterclock_line(p_g,cellula_idx);
    int cp_vsize = cp->vi.size();
    
        //counterclockwise sorted lines will
        for(int lidx: anticlockwise_lidx){
        int tmp_epoint_idx1=p_g->p_l[lidx]->vi[0];
        int tmp_epoint_idx2=p_g->p_l[lidx]->vi[1];

        geo_basic::pLine segment = std::make_pair(p_g->p_v[tmp_epoint_idx1]->loc,p_g->p_v[tmp_epoint_idx2]->loc);
        if(geo_basic::isIntersectSegmentLine(segment,division)){
            epoint_idx.push_back(tmp_epoint_idx1);
            epoint_idx.push_back(tmp_epoint_idx2);

            crosspoint.push_back(std::make_pair(lidx, geo_basic::crossPoint(segment, division)));
        }
        }

        if(crosspoint.size() != 2) {
            std::cout<<"Fatal bug of cross point"<<std::endl;
            //p_g->p_c[cellula_idx]->cellDivisionCount=100;
            //throw "Cannot divide cellulas";
            return;
        }
            p_g->p_c[cellula_idx]->cellDivisionCount++;
            p_g->p_c[cellula_idx]->cellTime =0.0;
            
        int lim_num1 = -1, lim_num2 = -1;

        for(int i = 0; i < cp_vsize; ++i){
            int vidx = cp->vi[i];
            if(vidx == epoint_idx[1]){
            lim_num1 = i;
            }
            if(vidx == epoint_idx[2]){
            lim_num2 = i-1;
            }
        }
        assert(lim_num1 != -1);
        assert(lim_num2 != -1);

        int v1_idx = p_g->p_v.size();
        int v2_idx = p_g->p_v.size() + 1;

        int l1_idx = p_g->p_l.size();
        int l2_idx = p_g->p_l.size() + 1;
        int l3_idx = p_g->p_l.size() + 2;

        int c1_idx = p_g->p_c.size();

        // vertex v1
        Vertex *v1 = new Vertex;
        v1->li.push_back(crosspoint[0].first);
        v1->li.push_back(l1_idx);
        v1->li.push_back(l3_idx);
        for(int c: p_g->p_l[crosspoint[0].first]->ci) {
            v1->ci.push_back(c);
        }
        v1->ci.push_back(c1_idx);
        v1->loc = crosspoint[0].second;
        v1->loc.z = 0.0;
        v1->vi = p_g->p_v.size();
        //v1->loc[1] = ???

        // vertex v2
        Vertex *v2 = new Vertex;
        v2->li.push_back(crosspoint[1].first);
        v2->li.push_back(l2_idx);
        v2->li.push_back(l3_idx);
        for(int c: p_g->p_l[crosspoint[1].first]->ci) {
            v2->ci.push_back(c);
        }
        v2->ci.push_back(c1_idx);
        v2->loc = crosspoint[1].second;
        v2->loc.z = 0.0;
        v2->vi = p_g->p_v.size()+1;

        //v2->loc[0][1] = ???

        // vertex epoint_idx[1]
        p_g->p_v[epoint_idx[1]]->li.push_back(l1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[1]]->li, crosspoint[0].first);
        p_g->p_v[epoint_idx[1]]->ci.push_back(c1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[1]]->ci, cellula_idx);

        // vertex epoint_idx[2]
        p_g->p_v[epoint_idx[2]]->li.push_back(l2_idx);
        division::findAndErase(p_g->p_v[epoint_idx[2]]->li, crosspoint[1].first);
        p_g->p_v[epoint_idx[2]]->ci.push_back(c1_idx);
        division::findAndErase(p_g->p_v[epoint_idx[2]]->ci, cellula_idx);

        for(int i = lim_num1 + 1; i <= lim_num2; ++i) {
            int vidx = cp->vi[i];
            p_g->p_v[vidx]->ci.push_back(c1_idx);
            division::findAndErase(p_g->p_v[vidx]->ci, cellula_idx);
        }

        // line l1
        Line *l1 = new Line;
        l1->vi[0] = v1_idx;
        l1->vi[1] = epoint_idx[1];
        for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
            l1->ci.push_back(cidx);
        }
        l1->ci.push_back(c1_idx);
        division::findAndErase(l1->ci, cellula_idx);
        l1->edgeForce = 0.0;
        l1->li = p_g->p_l.size();
        //l1->Balanced_Length = balanced_length;
        //l1->K_LENGTH = K_LENGTH;
        //l1->K2_LENGTH = K2_LENGTH;
        //l1->LENGTH_EQ = LENGTH_EQ;

        // line l2
        Line *l2 = new Line;
        l2->vi[0] = v2_idx;
        l2->vi[1] = epoint_idx[2];
        for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
            l2->ci.push_back(cidx);
        }
        l2->ci.push_back(c1_idx);
        division::findAndErase(l2->ci, cellula_idx);
        l2->edgeForce = 0.0;
        l2->li = p_g->p_l.size()+1;

        //l2->Balanced_Length = balanced_length;
        //l2->K_LENGTH = K_LENGTH;
        //l2->K2_LENGTH = K2_LENGTH;
        //l2->LENGTH_EQ = LENGTH_EQ;
        // line l3
        Line *l3 = new Line;
        l3->vi[0] = v1_idx;
        l3->vi[1] = v2_idx;
        l3->ci.push_back(cellula_idx);
        l3->ci.push_back(c1_idx);
        l3->edgeForce = 0.0;
        l3->li = p_g->p_l.size()+2;
        //l3->Balanced_Length = balanced_length;
        //l3->K_LENGTH = K_LENGTH;
        //l3->K2_LENGTH = K2_LENGTH;
        //l3->LENGTH_EQ = LENGTH_EQ;

        // line crosspoint[0].first
        p_g->p_l[crosspoint[0].first]->vi[1] = v1_idx;

        // line crosspoint[1].first
        p_g->p_l[crosspoint[1].first]->vi[0] = v2_idx;

        // line lim_num1-th --- lim_num2-th
        for(int i = lim_num1; i <= lim_num2; ++i) {
            int lidx = anticlockwise_lidx[i];
            p_g->p_l[lidx]->ci.push_back(c1_idx);
            division::findAndErase(p_g->p_l[lidx]->ci, cellula_idx);
        }

        // cellula c1
        Cell *c1 = new Cell;
        c1->vi.push_back(v1_idx);
        for(int i = lim_num1; i <= lim_num2 + 1; ++i) {
            int vidx = cp->vi[i];
            c1->vi.push_back(vidx);
        }
        c1->vi.push_back(v2_idx);
        c1->li.push_back(l1_idx);
        for(int i = lim_num1; i <= lim_num2; ++i) {
            int lidx = anticlockwise_lidx[i];
            c1->li.push_back(lidx);
        }
        c1->li.push_back(l2_idx);
        c1->li.push_back(l3_idx);
        c1->center = cp->center;
        
        //c1->Balanced_Area = balanced_area;
        //c1->K_AREA = K_AREA;

        c1->cellDivisionCount = p_g->p_c[cellula_idx]->cellDivisionCount;
        c1->tag = p_g->p_c[cellula_idx]->tag;
        c1->cellTime =0;
        // cellula cellula_idx
        cp->vi.erase(cp->vi.begin() + lim_num1, cp->vi.begin() + lim_num2 + 2);
        cp->vi.push_back(v1_idx);
        cp->vi.push_back(v2_idx);
        for(int i = lim_num1; i <= lim_num2; ++i) {
            //std::cerr << anticlockwise_lidx[i] << std::endl;
            division::findAndErase(cp->li, anticlockwise_lidx[i]);
        }
        cp->li.push_back(l3_idx);

        // cellula adjusting line-crosspoint[0].first (if exist)
        int adj_cidx1 = -1;
        for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
            if(cidx != cellula_idx) {
            adj_cidx1 = cidx;
            }
        }
        if(adj_cidx1 != -1) {
            p_g->p_c[adj_cidx1]->vi.push_back(v1_idx);
            p_g->p_c[adj_cidx1]->li.push_back(l1_idx);
        }

        // cellula adjusting line-crosspoint[1].first (if exist)
        int adj_cidx2 = -1;
        for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
            if(cidx != cellula_idx) {
            adj_cidx2 = cidx;
            }
        }
        if(adj_cidx2 != -1) {
            p_g->p_c[adj_cidx2]->vi.push_back(v2_idx);
            p_g->p_c[adj_cidx2]->li.push_back(l2_idx);
        }

        p_g->p_c.push_back(c1);
        p_g->p_l.push_back(l1);
        p_g->p_l.push_back(l2);
        p_g->p_l.push_back(l3);
        p_g->p_v.push_back(v1);
        p_g->p_v.push_back(v2);

        //isConsistent(p_g);
            //std::cout<<"hello 1"<<endl;
        //sortAntiClockwise(p_g);
        organ_geo::organ_vertex_counterclockwise_sort(p_g);
            //std::cout<<"hello 2"<<endl;

    }

    void Record(Organ *p_g, int cidx){
        DivisionRecord *dr = new DivisionRecord;
        dr->time = p_g->step;
        dr->cidx = cidx;
        dr->IsEpidermal = p_g->p_c[cidx]->IsEpidermal;
        dr->axisTheta = p_g->p_c[cidx]->axisTheta; 
        dr->center_x = p_g->p_c[cidx]->center.x;
        dr->center_y = p_g->p_c[cidx]->center.y;
        dr->division_count = p_g->p_c[cidx]->cellDivisionCount;
        //record frequency_modifier
        int inner_cell_number=0, peripheral_cell_number=0;
        double sum_in_frequency_modifier=0, sum_epi_frequency_modifier=0;
        for(int ci=0; ci<p_g->p_c.size(); ci++){
            if(p_g->p_c[ci]->IsEpidermal==0){
                sum_in_frequency_modifier+=p_g->p_c[ci]->frequency_modifier;
                inner_cell_number++;
            }
            else{
                sum_epi_frequency_modifier+=p_g->p_c[ci]->frequency_modifier;
                peripheral_cell_number++;
            }
        }
        double av_in_frequency_modifier = sum_in_frequency_modifier/inner_cell_number;
        double av_epi_frequency_modifier = sum_epi_frequency_modifier/peripheral_cell_number;
        dr->av_in_frequency_modifier = av_in_frequency_modifier;
        dr->av_epi_frequency_modifier = av_epi_frequency_modifier;
        p_g->d_r.push_back(dr);
    }

    void cell_time_initialization(Organ* p_g){
        if(division_control == "random_picking"){
    }
    else{
        if(division_control=="area"){
                division_frequency::area_control(p_g);
            }
            else if(division_control=="no"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="Gaussian_balance"){
                division_frequency::balance_control(p_g);
                division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
            }
            else if(division_control=="Gaussian_area"){
                division_frequency::area_control(p_g);
                division_frequency::Gaussian_control(p_g,gau_mu,gau_sigma);
            }
            else if(division_control=="balance"){
                division_frequency::balance_control(p_g);
            }
            else if(division_control=="balance_BF_3"){
                division_frequency::balance_control_BF_3(p_g);
            }
            else if(division_control=="balance_BF_1"){
                division_frequency::balance_control_BF_1(p_g);
            }
            else if(division_control=="temporal_Gaussian_linear_mu"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="temporal_Gaussian_linear_sigma"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="temporal_even_to_Gaussian"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="Gaussian_xy"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="uniform_to_Gaussian"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="biregion_frequency_position"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="biregion_frequency_angles_position"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="biregion_frequency_identity"){

                organ_geo::organ_center(p_g);
                double y_max_tmp = p_g->y_max_cell;
                double y_min_tmp = p_g->y_min_cell;
                for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
                    double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
                    if(y_relative_tmp>biregion_identity_y_boundary){
                        //cell ci belongs to the apical region
                        p_g->p_c[ci]->tag=1;
                    }
                    else{
                        //cell ci belongs to the basal region
                        p_g->p_c[ci]->tag=0;
                    }
                }

                division_frequency::biregion_frequency_identity(p_g,biregion_identity_relative_frequency);
            }
            
            else if(division_control=="arrest_front"){
                division_frequency::no_control(p_g);
                //division_frequency::area_control(p_g);
                //division_frequency::arrest_front(p_g,y_arrest_front,t_arrest_front,k_arrest_front);
            }
            else if(division_control=="arrest_front_biregion"){
                division_frequency::no_control(p_g);
            }
            else if(division_control=="meristem_position_relative_constant"){
                division_frequency::no_control(p_g);
            }
            else{
                std::cout<<"Fatal error: no cell division control is selected! (initialization)"<<endl;
                exit(-1);
            }
        
        if(in_division_direction=="biregion_angles_identity"){
            organ_geo::organ_center(p_g);
                double y_max_tmp = p_g->y_max_cell;
                double y_min_tmp = p_g->y_min_cell;
                for(int ci=0; ci<(int)p_g->p_c.size(); ci++){
                    double y_relative_tmp = (p_g->p_c[ci]->center.y-y_min_tmp)/(y_max_tmp-y_min_tmp);
                    if(y_relative_tmp>biregion_angles_identity_y_boundary){
                        //cell ci belongs to the apical region
                        std::cout<<"Cell "<<ci<<"belongs to the apical region"<<endl;
                        p_g->p_c[ci]->tag=1;
                    }
                    else{
                        //cell ci belongs to the basal region
                        std::cout<<"Cell "<<ci<<"belongs to the basal region"<<endl;
                        p_g->p_c[ci]->tag=0;
                    }
                }
        }

        random_device rnd;
        mt19937 mt(rnd());
        uniform_real_distribution <> cell_time_init(0,standard_cell_period_length);
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellTime = cell_time_init(mt)*p_g->p_c[ci]->frequency_modifier;
        }
    }
    }

    //cell time will be added for each cell
    void cellTimeAdd(Organ *p_g, int increase_index){
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellTime = p_g->p_c[ci]->cellTime+increase_index;
        }
    }

}