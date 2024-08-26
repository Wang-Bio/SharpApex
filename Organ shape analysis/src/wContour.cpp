#include "../include/wContour.h"

double pair_distance(const std::pair<double,double>& p1, const std::pair<double,double>& p2){return std::sqrt((p1.first-p2.first)*(p1.first-p2.first)+(p1.second-p2.second)*(p1.second-p2.second));};
std::pair<double,double> pair_mid(const std::pair<double,double>& p1, const std::pair<double,double>& p2){return std::make_pair((p1.first+p2.first)/2.0,(p1.second+p2.second)/2.0);};


namespace wContour{
    std::vector<std::pair<double,double>> normalization(const std::vector<std::pair<double,double>>& rawContour){
    //normalization was based on y~(0,1)
            bool debugMode=0;
        std::vector<std::pair<double,double>> normContour;
            if(rawContour.size()==0){std::cerr<<"rawContour has no point "<<std::endl;};
        double y_max=rawContour[0].second, y_min=rawContour[0].second;
        double x_center=0;
        for(int i=0;i<rawContour.size();i++){
            if(y_max<rawContour[i].second){
                y_max=rawContour[i].second;
            }
            if(y_min>rawContour[i].second){
                y_min=rawContour[i].second;
            }
            x_center+=rawContour[i].first/(double)rawContour.size();;
        }
            if(debugMode==1){std::cout<<"y_max "<<y_max<<", y_min"<<y_min<<", x_center "<<x_center<<std::endl;};
        double organ_length_approximate = y_max-y_min;
        //normalized: x_norm = (x-x_center)/organ_length_approximate; y_norm = (y-y_min)/organ_length_approximate
        for(int i=0;i<rawContour.size();i++){
            double x_norm = (rawContour[i].first-x_center)/organ_length_approximate;
            double y_norm = (rawContour[i].second-y_min)/organ_length_approximate;
            normContour.push_back(std::make_pair(x_norm,y_norm));
                if(debugMode==1){std::cout<<"point"<<i<<"before normalization "<<rawContour[i].first<<","<<rawContour[i].second<<"; after normalization "<<x_norm<<","<<y_norm<<std::endl;};
        }
        return normContour;
    }

    //imageJ get the y posi upside down, we need to flip it back to normal
    std::vector<std::pair<double,double>> flipY(const std::vector<std::pair<double,double>>& rawContour){
        std::vector<std::pair<double,double>> flipContour;
        for(int i=0;i<rawContour.size();i++){
            double x_flip = rawContour[i].first;
            double y_flip = -rawContour[i].second;
            flipContour.push_back(std::make_pair(x_flip,y_flip));
        }
        return flipContour;
    }

std::vector<std::pair<double,double>> parameterization(const std::vector<std::pair<double,double>>& normContour, const int& boundaryPointNum){
    //0. preparation
        std::vector<std::pair<double,double>> parameterizedContour;
            bool debugMode=0;
        double contourPerimeter = calcPerimeter(normContour);
            if(debugMode==1){std::cout<<"(wContour::parameterization) Contour perimeter: "<<contourPerimeter<<std::endl;};
        double contourPerimeter_averagePointDistance = contourPerimeter/(double)boundaryPointNum;
            //double contourPerimeter_averagePointDistance = 0.05;
            if(debugMode==1){std::cout<<"(wContour::parameterization) Contour perimeter average point distance "<<contourPerimeter_averagePointDistance<<std::endl;};
    //1. the first parameterized point, just pick the first point
        parameterizedContour.push_back(normContour[0]);
            if(debugMode==1){std::cout<<"The first parameterized point (it is just the first point) "<<parameterizedContour[0].first<<","<<parameterizedContour[0].second<<std::endl;};
    //2. the i-th boundary point: incremental algorithm
        double required_parameterized_distance = contourPerimeter_averagePointDistance;
        int current_point_index=1;
        bool distance_initial_is_boundary_vertex=0;
        double current_points_distance;
        for(int i=0;i<boundaryPointNum-1;i++){
                    if(debugMode==1){std::cout<<"Searching for the param point "<<i<<std::endl;}

            required_parameterized_distance_larger_than_normContour_points_distance:;
                    if(debugMode==1){std::cout<<"required_parameterized_distance "<<required_parameterized_distance<<std::endl;};
            if(distance_initial_is_boundary_vertex==0){
                current_points_distance = pair_distance(normContour[current_point_index],parameterizedContour[i]);
                    if(debugMode==1){std::cout<<"current points distance "<<current_points_distance<<" between normContour "<<current_point_index<<", parameterizedContour "<<i<<std::endl;};
            }
            else{
                current_points_distance = pair_distance(normContour[current_point_index],normContour[current_point_index-1]);
                    if(debugMode==1){std::cout<<"Current points distance "<<current_points_distance<<" between normContour "<<current_point_index<<","<<current_point_index-1<<std::endl;};
            }
            if(required_parameterized_distance<current_points_distance){
                //the new parameterized point will be generated on the line segment between parameterized point i and normContour current_point_index
                double param_x, param_y;
                if(distance_initial_is_boundary_vertex==0){
                    param_x = parameterizedContour[i].first + (normContour[current_point_index].first-parameterizedContour[i].first)*required_parameterized_distance/current_points_distance;
                    param_y = parameterizedContour[i].second + (normContour[current_point_index].second-parameterizedContour[i].second)*required_parameterized_distance/current_points_distance;
                    parameterizedContour.push_back(std::make_pair(param_x,param_y));
                }
                else{
                    param_x = normContour[current_point_index-1].first + (normContour[current_point_index].first-normContour[current_point_index-1].first)*required_parameterized_distance/current_points_distance;
                    param_y = normContour[current_point_index-1].second + (normContour[current_point_index].second-normContour[current_point_index-1].second)*required_parameterized_distance/current_points_distance;
                    parameterizedContour.push_back(std::make_pair(param_x,param_y));
                }
                        if(debugMode==1){std::cout<<"The param point "<<i<<" posi "<<param_x<<","<<param_y<<std::endl;};
            }
            else if(required_parameterized_distance==current_points_distance){
                //the new parameterized point is just the normContour current_point_index
                parameterizedContour.push_back(normContour[current_point_index]);
                current_point_index++;
            }
            else if(required_parameterized_distance>current_points_distance){
                //the new parameterized point is not on the line segmenet between parameterized point i nad normContour current_point_index
                required_parameterized_distance -= current_points_distance;
                //the new parameterized point should still be in the later line segment
                current_point_index++;
                distance_initial_is_boundary_vertex=1;
                goto required_parameterized_distance_larger_than_normContour_points_distance;
            }
            required_parameterized_distance = contourPerimeter_averagePointDistance;
            distance_initial_is_boundary_vertex=0;
        }

    return parameterizedContour;
}

    double calcPerimeter(const std::vector<std::pair<double,double>>& contour){
            bool debugMode=0;
        double perimeter=0;
        for(int i=0;i<contour.size();i++){
            double contour_length_i;
            if(i!=contour.size()-1)
                contour_length_i = std::sqrt((contour[i].first-contour[i+1].first)*(contour[i].first-contour[i+1].first)+(contour[i].second-contour[i+1].second)*(contour[i].second-contour[i+1].second));
            else
                contour_length_i = std::sqrt((contour[i].first-contour[0].first)*(contour[i].first-contour[0].first)+(contour[i].second-contour[0].second)*(contour[i].second-contour[0].second));

                if(debugMode==1){std::cout<<i<<" contour_length_i "<<contour_length_i<<std::endl;};
            perimeter+=contour_length_i;
        }
                if(debugMode==1){std::cout<<"Perimeter "<<perimeter<<std::endl;};
        return perimeter;
    }


std::vector<double> vec_curvature_circleFittingKasa_threeBoundaryPoints(const std::vector<std::pair<double,double>>& paramContour, const int& pointsAway, const int& numAveraging){
            bool debugMode=0;
            if(debugMode==1){std::cout<<"wContour::curvature_circleFittingKasa_threeBoundaryPoints: pointsAway "<<pointsAway<<"; numAveraging "<<numAveraging<<std::endl;};
    std::vector<double> curvature, curvature_av, curvature_absolute_av;

    for(int i=0;i<(int)paramContour.size();i++){
        std::pair<double,double> p_before, p_0, p_after;
        int p_before_index, p_after_index;
        p_0 = paramContour[i];
        if(i<pointsAway)
            p_before_index = i-pointsAway+paramContour.size();   
        else
            p_before_index = i-pointsAway;
    
        if(i>(paramContour.size()-pointsAway-1)){
            p_after_index = i+pointsAway-paramContour.size();
        }
        else{
            p_after_index = i+pointsAway;
        }
                if(debugMode==1){std::cout<<"p_before "<<p_before_index<<" i "<<i<<" p_after "<<p_after_index<<";";};
        p_before = paramContour[p_before_index];
        p_after = paramContour[p_after_index];  
        double curvature_absolute = single_curvature_circleFittingKasa_threeBoundaryPoints(p_before,p_0,p_after);
        curvature_absolute_av.push_back(curvature_absolute);
        double curvature_value;
        std::pair<double,double> mid_point = pair_mid(p_before,p_after);
        if(pointInPolygon_rayCasting(mid_point,paramContour)==0){
            curvature_value = -curvature_absolute;
        }
        else{
            curvature_value = curvature_absolute;
        }
                if(debugMode==1){std::cout<<" curvature absolute "<<curvature_absolute<<"; curvature "<<curvature_value<<std::endl;};
        curvature.push_back(curvature_value);
    }
                //if(debugMode==1){vd_fout(curvature_absolute_av,fileDirectory+fileName+"_curvature_absolute.txt");};
                //if(debugMode==1){vd_fout(curvature,fileDirectory+fileName+"_curvature.txt");};
    curvature_av = vd_average(curvature,numAveraging);

    return curvature_av;
}

double single_curvature_circleFittingKasa_threeBoundaryPoints(const std::pair<double,double>& p_before, const std::pair<double,double>& p_0, const std::pair<double,double>& p_after){
//reference: 
        //1. http://www.ne.jp/asahi/paleomagnetism.rock-magnetism/basics/pmag/circ/circ1E.html
        //2. https://people.cas.uab.edu/~mosya/cl/CircleFitByKasa.cpp        
        //generating the matrix
        double Mx=0, My=0, Mz=0, Mxy=0, Mxz=0, Myz=0, Mxx=0, Myy=0;
        double B_fit, C_fit, D_fit;
        double a_fit, b_fit, R_fit;

        std::vector<std::pair<double,double>> vv_fitting = std::vector<std::pair<double,double>> {p_before,p_0,p_after};
        int n_fit = (int)vv_fitting.size();
        for(int vi=0; vi<(int)vv_fitting.size(); vi++){
            double x_tmp = vv_fitting[vi].first;
            double y_tmp = vv_fitting[vi].second;
            double z_tmp = x_tmp*x_tmp+y_tmp*y_tmp;

            Mx  += x_tmp;
            My  += y_tmp;
            Mz  += z_tmp;
            Mxy += x_tmp*y_tmp;
            Mxz += x_tmp*z_tmp;
            Myz += y_tmp*z_tmp;
            Mxx += x_tmp*x_tmp;
            Myy += y_tmp*y_tmp;
        }

        //solving the linear matrix
        // A * X = B => X
        double l11 = sqrt(Mxx);
        double l21 = Mxy/l11;
        double l22 = sqrt(Myy-l21*l21);
        double l31 = Mx/l11;
        double l32 = (My-l31*l21)/l22;
        double l33 = sqrt(n_fit-(l31*l31+l32*l32));

        double y1= - Mxz/l11;
        double y2= - (Myz+l21*y1)/l22;
        double y3= - (Mz+l31*y1+l32*y2)/l33;

        D_fit = y3/l33;
        C_fit = (y2-l32*D_fit)/l22;
        B_fit = (y1-l21*C_fit-l31*D_fit)/l11;

        a_fit = -B_fit/2.0;
        b_fit = -C_fit/2.0;
        R_fit = sqrt(a_fit*a_fit+b_fit*b_fit-D_fit);

        return 1.0/R_fit;
}

bool pointInPolygon_rayCasting(const std::pair<double,double>& mid_point, const std::vector<std::pair<double,double>>& paramContour){
    bool insideBool = 0;
    int ray_casting_index=0;
    int n = (int)paramContour.size();
            bool debugMode=1;
    for(int i=0;i<n;++i){
        int next = (i+1)%n;
        double x0_tmp = paramContour[i].first;
        double y0_tmp = paramContour[i].second;
        double x1_tmp = paramContour[next].first;
        double y1_tmp = paramContour[next].second;

        double x_point = mid_point.first;
        double y_point = mid_point.second;
        if((y0_tmp>y_point)!=(y1_tmp>y_point)){
                if(x_point<(x0_tmp+(y_point-y0_tmp)/(y1_tmp-y0_tmp)*(x1_tmp-x0_tmp))){
                    ray_casting_index++;
                }
            }
    }
    if(ray_casting_index%2==0){
            //even
            insideBool = 0;
        }
        if(ray_casting_index%2==1){
            //odd
            insideBool = 1;
        }
        //std::cout<<"ray_casting_index "<<ray_casting_index<<endl;
        return insideBool;
}

std::vector<double> vd_average(const std::vector<double>& vd, const int& numAveraging){
    std::vector<double> output_vd(vd.size());

        int halfNeighbors = (numAveraging-1)/2;

        for(int i=0;i<vd.size();++i){
            double sum=0.0;

            for(int j=i-halfNeighbors-1;j<i+halfNeighbors;++j){
                if(j>=0 && j<vd.size()){
                    sum+=vd[j];
                }
                else if(j<0){
                    sum+=vd[vd.size()-1+j];
                }
                else if(j>vd.size()){
                    sum+=vd[j-vd.size()];
                }
            }

            output_vd[i] = sum/numAveraging;
        }

        return output_vd;
}

void vd_fout(const std::vector<double>& vd, const std::string& fileName){
    if(vd.size()==0){
        std::cerr<<"No elements in vd (vd_fout for"<<fileName<<")"<<std::endl;
        exit(1);
    }
    std::ofstream fout(fileName);
    for(int i=0;i<vd.size();i++){
        fout<<i<<" "<<vd[i]<<std::endl;
    }
    fout.close();
}
}