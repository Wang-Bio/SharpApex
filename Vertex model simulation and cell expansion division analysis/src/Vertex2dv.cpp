#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"

    
    void Vertex::print_Cartesian(void){
        std::cout<<"("<<loc.x<<","<<loc.y<<")"<<std::endl;
    }

    void Vertex::print_Polar(void){
        std::cout<<"("<<r<<","<<theta<<")"<<std::endl;
    }

    double Vertex::distance_from_vertex(Vertex v2){
        return sqrt((loc.x-v2.loc.x)*(loc.x-v2.loc.x)+(loc.y-v2.loc.y)*(loc.y-v2.loc.y));
    }

    void Vertex::Cartesian_to_Polar(void){
        r = sqrt(loc.x*loc.x+loc.y*loc.y);
        theta = atan(loc.y/loc.x); 
    }

    void Vertex::Polar_to_Cartesian(void){
        loc.x = r*cos(theta);
        loc.y = r*sin(theta);
    }

    bool Vertex::collinear_points(Vertex v2,Vertex v3){
        //if the slopes of l12 and l13 are the same and so does the intercepts of l12 and l13, then point p1,p2,p3 is collinear, otherwise not collinear
        double slope_12 = (loc.y-v2.loc.y)/(loc.x-v2.loc.x);
        double slope_23 = (v2.loc.y-v3.loc.y)/(v2.loc.x-v3.loc.x);

        double intercept_12 = (loc.x*v2.loc.y-v2.loc.x*loc.y)/(loc.x-v2.loc.x);
        double intercept_23 = (v2.loc.x*v3.loc.y-v3.loc.x*v2.loc.y)/(v2.loc.x-v3.loc.x);

        if(abs(slope_12-slope_23)<EPS_geo&&abs(intercept_12-intercept_23)<EPS_geo)
        {
            return 1;
        }
        else{
            return 0;
        }
    }
    
    bool Vertex::same_vertex(Vertex v2){
        if(abs(distance_from_vertex(v2))<EPS_geo){
            return 1;
        }
        else{
            return 0;
        }
    }


namespace vVertex{
    std::vector<Vertex> read_from_txt(const std::string &file, int x_column, int y_column){
        std::vector<Vertex> vv;
        std::ifstream inFile(file);
        std::string line;

        if(!inFile.is_open()){
            std::cerr << "Error opening file: "<<file<<std::endl;
            exit(1);
        }

        while(getline(inFile,line)){
            std::istringstream iss(line);
            std::string token;
            double value;
            std::vector<double> values;
            while (getline(iss, token, '\t')){
                std::istringstream converter(token);
                if(converter >> value){
                    values.push_back(value);
                }
            }

            if (x_column < values.size() && y_column < values.size()) {
                Vertex v_tmp;
                v_tmp.loc.x = values[x_column];
                v_tmp.loc.y = values[y_column];
                vv.push_back(v_tmp);
            }
        }

        inFile.close();
        return vv;
    }

    void save_txt(const std::string & file, const std::vector<Vertex> & vv){
        std::ofstream fout(file);
        for(auto v:vv){
            fout<<v.loc.x<<" "<<v.loc.y<<std::endl;
        }

        fout.close();
    }

    void print(const std::vector<Vertex> &vv){
        for(auto v_tmp : vv){
            std::cout<<v_tmp.loc.x<<" "<<v_tmp.loc.y<<std::endl;
        }
    }

    std::vector<Vertex> vertical_reflection(const std::vector<Vertex> & vv){
        std::vector<Vertex> post_vv;
        for(auto v:vv){
            Vertex v_tmp;
            v_tmp.loc.x = v.loc.x;
            v_tmp.loc.y = -v.loc.y;
            post_vv.push_back(v_tmp);
        }
        return post_vv;
    }

    std::vector<Vertex> normalization(const std::vector<Vertex> &vv){
        vector<Vertex> vv_normalized;

    //1. calculate y_max and y_min
        double y_max = vv[0].loc.y;
        double y_min = vv[0].loc.y;
        double x_center = vv[0].loc.x;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y<y_min){
                y_min=vv[vi].loc.y;
            }
            if(vv[vi].loc.y>y_max){
                y_max=vv[vi].loc.y;
                x_center=vv[vi].loc.x;
            }
        }
            //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

    //2. calculate vv_normalization
        double norm_length = y_max-y_min;
        //double x_center=0;
        //for(int vi=0;vi<vv.size();vi++){
        //    x_center+=vv[vi].loc.x/(double)vv.size();
        //}
            //std::cout<<"x_center: "<<x_center<<endl;
            //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex v_tmp;
            v_tmp.loc.x = (vv[vi].loc.x-x_center)/norm_length;
            v_tmp.loc.y = (vv[vi].loc.y-y_min)/norm_length;
                //std::cout<<"y "<<vv[vi].loc.y<<";y_norm "<<v_tmp.loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    int index_for_y(const std::vector<Vertex> & vv,double y){
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y==y){
                return vi;
            }
        }
    } 

    std::vector<Vertex> sample(const std::vector<Vertex> &vv){
        //sampled vertex points must be within the y of two Vertices, then just do the linear interpolation
        //the left part and right part are calculated differently and thus, the division between left part and right part should be made
        //the vv should be ordered in anticlockwise direction
        
        int left_part = index_for_y(vv,1.0);

        double sampling_distance=0.01;
        std::vector<Vertex> sampledVertex;
        int numSamples = 2/sampling_distance;
        for(int i=0;i<numSamples/2;++i){
            double y = i*sampling_distance;
            for(int j=0; j<left_part; j++){
                if(vv[j].loc.y<= y&&y<=vv[j+1].loc.y){
                    double t = (y-vv[j].loc.y)/(vv[j+1].loc.y-vv[j].loc.y);
                    Vertex newVertex;
                    newVertex.loc.x = vv[j].loc.x + t*(vv[j+1].loc.x-vv[j].loc.x);
                    newVertex.loc.y= y;
                    sampledVertex.push_back(newVertex);
                    break;
                }
            }
        }

        for(int i=0;i<numSamples/2;++i){
            double y=1-i*sampling_distance;
            //cout<<"y "<<y;
            for(int j=left_part; j<vv.size(); j++){
                if(vv[j].loc.y>=y&&y>=vv[j+1].loc.y){
                    double t = (y-vv[j].loc.y)/(vv[j+1].loc.y-vv[j].loc.y);
                    Vertex newVertex;
                    newVertex.loc.x = vv[j].loc.x + t*(vv[j+1].loc.x-vv[j].loc.x);
                    //cout<<" x "<<newVertex.loc.x<<std::endl;
                    newVertex.loc.y= y;
                    sampledVertex.push_back(newVertex);
                    break;
                }
            }
        }

        return sampledVertex;
    }

    //can only be used when all vv has been sampled into the same numebr and in the same order 
    std::vector<Vertex> vec_averaged(const std::vector<std::vector<Vertex>> & vvv){
        std::vector<Vertex> post_vv;
        int vv_size = vvv[0].size();
        int vvv_size = vvv.size();
        for(int i=0;i<vv_size;i++){
            Vertex v_av;
            for(auto vv:vvv){
               v_av.loc+=(vv[i].loc/(double)vvv_size);
            }
            post_vv.push_back(v_av);
        }
        return post_vv;
    }


}