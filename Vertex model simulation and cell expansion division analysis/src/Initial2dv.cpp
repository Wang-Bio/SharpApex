#include "../include/Initial2dv.h"


string real_organ_contour_imagej_txt = "/mnt/c/again_and_again/codes/git/plant_vertex_model/analysis_example/mature_ts_contour_for_similarity_index.txt";

namespace initialization{
    void organ(Organ* p_g){
        
        std::cout<<"*****Start Initialization*****"<<std::endl;
        geo::calcGeometrics(p_g);
        //cout<<"similarity_calculation_required: "<<similarity_calculation_required<<endl;
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
            //cout_fout_debug::cout_vector_vertex(real_organ_contour_processed_for_similarity_index);
        }
        p_g->step=0;
        output::VTK(p_g);
        std::cout<<"Start force initialization"<<std::endl;
        force::forceShapeInitiation(p_g,200000);
        std::cout<<"End force initialization"<<std::endl;
        division::cell_time_initialization(p_g);
        p_g->step=1;
        output::VTK(p_g);
        output::geo_initial();
        geo::calcGeometrics(p_g);
        std::cout<<"*****End Initialization*****"<<std::endl;
        //output::geo(p_g);    
    }

    void organ_continue(Organ* p_g){

        std::cout<<"******Start Initialization******"<<std::endl;
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
            //cout_fout_debug::cout_vector_vertex(real_organ_contour_processed_for_similarity_index);
        }
        geo::calcGeometrics(p_g);
        std::cout<<"*****End Initialization*****"<<std::endl;
    }
    
}