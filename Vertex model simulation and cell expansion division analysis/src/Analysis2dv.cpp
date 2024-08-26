#include "../include/Analysis2dv.h"

void analysis_2D::add_param_name(std::string parameter_name_1, std::string parameter_name_2){
    param_name_1 = parameter_name_1;
    param_name_2 = parameter_name_2;
}

void analysis_2D::add_param_single(double param_1,double param_2){
    param_2D.push_back(make_pair(param_1,param_2));
}

void analysis_2D::add_param_rect(std::vector<double> param_1,std::vector<double> param_2){
    for(int i=0;i<param_1.size();i++){
        for(int j=0;j<param_2.size();j++){
            param_2D.push_back(make_pair(param_1[i],param_2[j]));
        }
    }
}

void analysis_2D::add_param_rect_incremental(double param_1_init,double param_1_term,double param_1_incremental,double param_2_init,double param_2_term,double param_2_incremental){
    int param_1_size = (param_1_term-param_1_init)/param_1_incremental;
    int param_2_size = (param_2_term-param_2_init)/param_2_incremental;

    for(int i=0;i<param_1_size;i++){
        for(int j=0;j<param_2_size; j++){
            double param_1_tmp = param_1_init+i*param_1_incremental;
            double param_2_tmp = param_2_init+j*param_2_incremental;
            param_2D.push_back(make_pair(param_1_tmp,param_2_tmp));
        }
    }    
}

void analysis_2D::add_param_triang_right_up(std::vector<double> param_1,std::vector<double> param_2){
    for(int i=0;i<param_1.size();i++){
        for(int j=0;j<param_2.size();j++){
            if(param_1[i]<param_2[j]||param_1[i]==param_2[j])
            param_2D.push_back(make_pair(param_1[i],param_2[j]));
        }
    }
}

void analysis_2D::add_param_triang_right_up_incremental(double param_1_init,double param_1_term,double param_1_incremental,double param_2_init,double param_2_term,double param_2_incremental){
    cout<<"param_1_term: "<<param_1_term<<" param_1_init: "<<param_1_init<<" param_1_incremental: "<<param_1_incremental<<endl;
    double param_1_size = (param_1_term-param_1_init)/param_1_incremental+1;
    double param_2_size = (param_2_term-param_2_init)/param_2_incremental+1;
    std::cout<<"param_1_size: "<<param_1_size<<", param_2_size: "<<param_2_size<<std::endl;
    for(int i=0;i<param_1_size;i++){
        for(int j=0;j<param_2_size; j++){
            double param_1_tmp = param_1_init+i*param_1_incremental;
            double param_2_tmp = param_2_init+j*param_2_incremental;
            if(param_1_tmp>param_2_tmp||param_1_tmp==param_2_tmp)
            param_2D.push_back(make_pair(param_1_tmp,param_2_tmp));
        }
    }    
}

//file_directory
void analysis_2D::add_all_file_directory(std::string basic_directory){
    for(int i=0;i<size();i++){
        char parameter_1_file[100];
        sprintf(parameter_1_file,"_%.1f",param_2D[i].first);
        std::string parameter_1_file_s = param_name_1+parameter_1_file;

        char parameter_2_file[100];
        sprintf(parameter_2_file,"_%.1f",param_2D[i].second);
        std::string parameter_2_file_s = param_name_2+parameter_2_file;
        std::string complete_directory = basic_directory + parameter_1_file_s+"_"+parameter_2_file_s+"/";
        file_directory.push_back(complete_directory);
    }
}

void analysis_2D::add_all_file_directory_3(std::string basic_directory){
    for(int i=0;i<size();i++){
        char parameter_1_file[100];
        sprintf(parameter_1_file,"_%.0f",param_2D[i].first);
        std::string parameter_1_file_s = param_name_1+parameter_1_file;

        char parameter_2_file[100];
        sprintf(parameter_2_file,"_%.0f",param_2D[i].second);
        std::string parameter_2_file_s = param_name_2+parameter_2_file;
        std::string param_name_3 = "y_30";
        std::string parameter_3_file_s = param_name_3;
        std::string complete_directory = basic_directory + parameter_3_file_s + "_" + parameter_1_file_s + "_" + parameter_2_file_s;
        file_directory.push_back(complete_directory);
    }
}

std::string analysis_2D::get_file_directory(int index){
    return file_directory[index];
}

std::pair<std::string, std::string> analysis_2D::VTK_directory(int repeat_time_tmp, int index){
    std::string Cell_VTK, Line_VTK;

    char repeat_file[100];
    sprintf(repeat_file,"%d",repeat_time_tmp);
    std::string repeat_file_s = repeat_file;
    std::string VTK_file_path = get_file_directory(index)+repeat_file_s;

    int length_tmp = VTK_file_path.length();
    char directory_tmp[length_tmp+1];
    strcpy(directory_tmp, VTK_file_path.c_str());        
    chdir(directory_tmp);
    int file_index = file_process::getFileNum(VTK_file_path);
    cout<<"file_index: "<<file_index<<endl;
    if(file_index==0){
        cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
        exit(1);
    }
    char str_cell[100], str_line[100];
    sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-2);
    sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-2);
    cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
    Cell_VTK = VTK_file_path+"/"+str_cell;
    Line_VTK = VTK_file_path+"/"+str_line;

    return make_pair(Cell_VTK,Line_VTK);
}

std::pair<std::string, std::string> analysis_2D::VTK_directory(int index){
    std::string Cell_VTK, Line_VTK;
    std::string VTK_file_path = get_file_directory(index);

    int length_tmp = VTK_file_path.length();
    char directory_tmp[length_tmp+1];
    strcpy(directory_tmp, VTK_file_path.c_str());        
    chdir(directory_tmp);
    int file_index = file_process::getFileNum(VTK_file_path);
    cout<<"file_index: "<<file_index<<endl;
    if(file_index==0){
        cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
        exit(1);
    }
    char str_cell[100], str_line[100];
    sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-5);
    sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-5);
    cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
    Cell_VTK = VTK_file_path+"/"+str_cell;
    Line_VTK = VTK_file_path+"/"+str_line;

    return make_pair(Cell_VTK,Line_VTK);
}

std::string analysis_2D::get_parameter_s(int index){
    char parameter_1_tmp[10];
    sprintf(parameter_1_tmp,"%.2f",param_2D[index].first);
    std::string parameter_1_s = parameter_1_tmp;
    char parameter_2_tmp[10];
    sprintf(parameter_2_tmp,"%.2f",param_2D[index].second);
    std::string parameter_2_s = parameter_2_tmp;
    return param_name_1+"_"+parameter_1_s+"_"+param_name_2+"_"+parameter_2_s;
}


//basic_properties
void analysis_2D::set_repeat_time(int repeat_time_tmp){
    repeat_time = repeat_time_tmp;
}

int analysis_2D::get_repeat_time(void){
    return repeat_time;
}

int analysis_2D::size(void){
    return param_2D.size();
}

std::pair<double,double> analysis_2D::get_parameter(int index){
    return param_2D[index];
}

std::pair<std::string,std::string> analysis_2D::get_param_name(void){
    return make_pair(param_name_1,param_name_2);
}


//output
void analysis_2D::print(void){
    std::cout<<param_name_1<<" "<<param_name_2<<std::endl;
    for(int i=0;i<size();i++){
        cout<<param_2D[i].first<<" "<<param_2D[i].second<<std::endl;
    }
}

//geometric_recording
void analysis_2D::geometrics_recording(Organ* p_g,int index,int repeat_time_tmp){
    std::cout<<"**********Start geometrics recording**********"<<std::endl;
    add_curvature_analysis(p_g,index,repeat_time_tmp);
    add_similarity_analysis(p_g,index,repeat_time_tmp);
    std::cout<<"**********Geometrics recording finished**********"<<std::endl;
}

void analysis_2D::add_curvature_analysis(Organ* p_g, int index, int repeat_time_tmp){
    min_curvature_all.push_back(p_g->minimum_curvature);
    max_curvature_all.push_back(p_g->maximum_curvature);
    range_curvature_all.push_back(p_g->maximum_curvature-p_g->minimum_curvature);
    accumulated_negative_curvature_all.push_back(p_g->accumulated_negative_curvature);
    //std::cout<<"Repeat_time "<<repeat_time<<", repeat_time_tmp "<<repeat_time_tmp<<std::endl;
    if(repeat_time==1){
        min_curvature_av.push_back(p_g->minimum_curvature);
        max_curvature_av.push_back(p_g->maximum_curvature);
        range_curvature_av.push_back(p_g->maximum_curvature-p_g->minimum_curvature);
        accumulated_negative_curvature_av.push_back(p_g->accumulated_negative_curvature);
    }
    else{
        if(repeat_time_tmp==(repeat_time-1)){
            double min_curvature_av_tmp=0, max_curvature_av_tmp=0, range_curvature_av_tmp=0, accumulated_negative_curvature_av_tmp=0;
            for(int i=0;i<repeat_time;i++){
                min_curvature_av_tmp +=min_curvature_all[index*repeat_time+i]/repeat_time;
                max_curvature_av_tmp +=max_curvature_all[index*repeat_time+i]/repeat_time;
                range_curvature_av_tmp +=range_curvature_all[index*repeat_time+i]/repeat_time;
                accumulated_negative_curvature_av_tmp +=accumulated_negative_curvature_all[index*repeat_time+i]/repeat_time;
            }
            min_curvature_av.push_back(min_curvature_av_tmp);
            max_curvature_av.push_back(max_curvature_av_tmp);
            range_curvature_av.push_back(range_curvature_av_tmp);
            accumulated_negative_curvature_av.push_back(accumulated_negative_curvature_av_tmp);
        }
    }    
}

void analysis_2D::add_similarity_analysis(Organ* p_g,int index, int repeat_time_tmp){
    similarity_index_all.push_back(p_g->similarity_index);
    if(repeat_time==1){
        similarity_index_av.push_back(p_g->similarity_index);
    }
    else{
        if(repeat_time_tmp==(repeat_time-1)){
            double similarity_index_av_tmp=0;
            for(int i=0;i<repeat_time;i++){
                similarity_index_av_tmp +=similarity_index_all[index*repeat_time+i]/repeat_time;
            }
            similarity_index_av.push_back(similarity_index_av_tmp);
        }
    }
}

void analysis_2D::geometrics_output(std::string analysis_filename){
    std::string output_file_min_curvature = analysis_filename+"min_curvature.txt";
    std::ofstream fout1(output_file_min_curvature);
    for(int i=0;i<size();i++){
        fout1<<i<<" "<<param_2D[i].first<<" "<<param_2D[i].second<<" "<<min_curvature_av[i]<<endl;
    }
    fout1.close();
    
    std::string output_file_accumulated_negative_curvature = analysis_filename+"accumulated_negative_curvature.txt";
    std::ofstream fout2(output_file_accumulated_negative_curvature);
    for(int i=0;i<size();i++){
        fout2<<i<<" "<<param_2D[i].first<<" "<<param_2D[i].second<<" "<<accumulated_negative_curvature_av[i]<<endl;
    }
    fout2.close();

    std::string output_file_similarity_index = analysis_filename+"similarity_index.txt";
    std::ofstream fout3(output_file_similarity_index);
    for(int i=0;i<size();i++){
        fout3<<i<<" "<<param_2D[i].first<<" "<<param_2D[i].second<<" "<<similarity_index_av[i]<<endl;
    }
    fout3.close();
}
