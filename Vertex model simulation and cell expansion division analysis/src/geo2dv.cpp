/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include "../include/geo2dv.h"
#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"

bool similarity_calculation_required=0;
vector<Vertex*> real_organ_contour_processed_for_similarity_index;

namespace vertex_geo{
    //calculate the distance between vertex v1 and vertex v2
    double vertex_distance(Vertex* v1, Vertex* v2){
        double vertex_distance_tmp=0;
        vertex_distance_tmp = sqrt((v1->loc.y-v2->loc.y)*(v1->loc.y-v2->loc.y)+(v1->loc.x-v2->loc.x)*(v1->loc.x-v2->loc.x));
        return vertex_distance_tmp;
    }

    double vertex_distance(Vertex v1, Vertex v2){
        double vertex_distance_tmp=0;
        vertex_distance_tmp = sqrt((v1.loc.y-v2.loc.y)*(v1.loc.y-v2.loc.y)+(v1.loc.x-v2.loc.x)*(v1.loc.x-v2.loc.x));
        return vertex_distance_tmp;
    }

    //judge the relationship between vertex v1 and vertex v2; if return=0, they are different points; if return=1, they are so close that be seemed as the same point 
    bool vertex_relationship(Vertex* v1, Vertex* v2){
        bool index_tmp = 0;
        double epsilon_tmp =  10e-6;

        if(vertex_geo::vertex_distance(v1,v2)<epsilon_tmp){
            index_tmp = 1;
        }
        else{}

        return index_tmp;
    }

}

namespace line_geo{
    //calculate the length of line li
    double line_length(Organ* p_g, int li){
        p_g->p_l[li]->length = (p_g->p_v[p_g->p_l[li]->vi[0]]->loc - p_g->p_v[p_g->p_l[li]->vi[1]]->loc).norm();
        return p_g->p_l[li]->length;
    }


    //calculate the intersection of a line defined with slope and intercept, with cell wall li. if they are parallel _vec.z = 1.0; if they simply have no intersection _vec.z = -1.0
    _vec<double> line_cell_wall_intersection(double slope, double intercept, Organ* p_g, int li){
        _vec<double> intersection_result;
        //std::cout<<"li "<<li<<endl;
        double cell_wall_slope = (p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y-p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y)/(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x-p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
        double cell_wall_intercept = p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y-cell_wall_slope*p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x;
        //std::cout<<"cell wall slope: "<<cell_wall_slope<<" , cell_wall_intercept: "<<cell_wall_intercept<<endl;
        
        //now we have the line slope, line intercept, cell wall slope and cell wall intercept: y=a1x+b1, y=a2x+b2
        // => there are three possibilities for two lines: they could be parallel, intersected, or identical (ignore identical situation)
        //parallel: a1=a2
        //std::cout<<slope<<" "<<cell_wall_slope<<endl;
        if(abs(slope-cell_wall_slope)<1e-9){
            intersection_result = _vec<double> {0.0,0.0,1.0};
            //std::cout<<"parallel"<<endl;
        }
        //intersected
        else{
                double intersection_x = (cell_wall_intercept - intercept)/(slope - cell_wall_slope);
                double intersection_y = cell_wall_slope*intersection_x + cell_wall_intercept;
                //std::cout<<"intersection_x "<<intersection_x<<", intersection_y "<<intersection_y<<endl;
                double segment_min_y = std::min(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y);
                double segment_max_y = std::max(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y); 
                double segment_min_x = std::min(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
                double segment_max_x = std::max(p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x,p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x);
                
                intersection_result = _vec<double> {intersection_x,intersection_y,0.0};
                if(intersection_y <= segment_max_y && intersection_y >= segment_min_y &&intersection_x >= segment_min_x&&intersection_x <= segment_max_x){
                        //this intersection is correct
                        //std::cout<<"li "<<li<<"; endpoint 1 for cell wall "<<p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x<<","<<p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y<<"; endpoint 2 for cell wall "<<p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x<<","<<p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y<<"; Line slope: "<<slope<<", line intercept "<<intercept<<";cell wall slope: "<<cell_wall_slope<<" , cell_wall_intercept: "<<cell_wall_intercept<<"; intersection_x "<<intersection_x<<", intersection_y "<<intersection_y<<endl;
                        //std::cout<<" "<<std::endl;
                }
                else{
                    intersection_result.z=-1.0;
                    //this intersection is out of the line segment
                }    
    }
    return intersection_result;
}

    Distance_point_line distance_line_segment_to_vertex(Line li, Vertex p1){
        //reference 1: http://www.csharphelper.com/howtos/howto_point_segment_distance.html
        //reference 2: https://www.youtube.com/watch?v=egmZJU-1zPU
        Distance_point_line Result;
        Vertex closest_point;
        Vertex d1 = li.d1;
        Vertex d2 = li.d2;
        double dx = d2.loc.x - d1.loc.x;
        double dy = d2.loc.y - d1.loc.y;
        if((dx==0)&&(dy==0)){
            //It's a point not a line segment
            closest_point = d1;
            dx = p1.loc.x - d1.loc.x;
            dy = p1.loc.y - d1.loc.y;
            Result.distance=sqrt(dx*dx+dy*dy);
            Result.Closest_Point=d1;
            return Result;
        }

        // Calculate the t that minimizes the distance
        double t = ((p1.loc.x - d1.loc.x)*dx + (p1.loc.y-d1.loc.y)*dy)/(dx*dx+dy*dy);
        Result.t =t;
        //std::cout<<"t "<<t<<endl;
        // See if this represents one of the segment's ends points or a point in the middle
        if(t<0){
            //std::cout<<"the closest point should be the end point d1"<<endl;
            closest_point = d1;
        }
        else if(t>1){
            //std::cout<<"the closest point should be the end point d2"<<endl;
            closest_point = d2;
        }
        else{
            //std::cout<<"the closest point is between d1 and d2"<<endl;
            closest_point.loc = _vec<double>{d1.loc.x + t*dx, d1.loc.y + t*dy,0.0};
            
        }
        //std::cout<<"Closest";
        //closest.print_Cartesian();
        Result.Closest_Point = closest_point;
        Result.distance = p1.distance_from_vertex(closest_point);
        return Result;
    }

    //calculate line segment distance to vertex
    Distance_point_line distance_line_segment_to_vertex(Organ* p_g,int li,Vertex p1){
        //reference 1: http://www.csharphelper.com/howtos/howto_point_segment_distance.html
        //reference 2: https://www.youtube.com/watch?v=egmZJU-1zPU
        Distance_point_line Result;
        
        Vertex closest_point;
        Vertex d1 = *p_g->p_v[p_g->p_l[li]->vi[0]];
        Vertex d2 = *p_g->p_v[p_g->p_l[li]->vi[1]];

        double dx = d2.loc.x - d1.loc.x;
        double dy = d2.loc.y - d1.loc.y;
        if((dx==0)&&(dy==0)){
            //It's a point not a line segment
            closest_point = d1;
            dx = p1.loc.x - d1.loc.x;
            dy = p1.loc.y - d1.loc.y;
            Result.distance=sqrt(dx*dx+dy*dy);
            Result.Closest_Point=d1;
            return Result;
        }

        // Calculate the t that minimizes the distance
        double t = ((p1.loc.x - d1.loc.x)*dx + (p1.loc.y-d1.loc.y)*dy)/(dx*dx+dy*dy);
        Result.t =t;
        //std::cout<<"t "<<t<<endl;
        // See if this represents one of the segment's ends points or a point in the middle
        if(t<0){
            //std::cout<<"the closest point should be the end point d1"<<endl;
            closest_point = d1;
        }
        else if(t>1){
            //std::cout<<"the closest point should be the end point d2"<<endl;
            closest_point = d2;
        }
        else{
            //std::cout<<"the closest point is between d1 and d2"<<endl;
            closest_point.loc = _vec<double>{d1.loc.x + t*dx, d1.loc.y + t*dy,0.0};
            
        }
        //std::cout<<"Closest";
        //closest.print_Cartesian();
        Result.Closest_Point = closest_point;
        Result.distance = p1.distance_from_vertex(closest_point);
        return Result;
    }
    
    //=0, the vertex is outside of the line; =1, the vertex is inside the line
    bool line_vertex_relationship(Line l1, Vertex p1){
        //check if (p1.loc.x,p1.loc.y) satisfy the l1 eqution y = slope*x+intercept 
        double y_p1 = p1.loc.x*l1.slope+l1.intercept;
        if(abs(y_p1-p1.loc.y)<EPS_geo){
            //the vertex is inside the line
            return true;
        }else{
            //the vertex is outside the line
            return false;
        }

    }

    //=0, the vertex is outside of the linesegment; =1, the vertex is inside the line segment
    bool line_segment_vertex_relationship(Line ls1, Vertex p1){
        //check if the relationship between Vertex p1 and Line 
        if(line_vertex_relationship(ls1,p1)==0){
            return false;
        }
        else{
            //the Vertex p1 is inside the line, but maybe outside of the line segment
            //check if the (p1.x,p1.y) is inside the two endpoints
            if(p1.loc.x >= min(ls1.d1.loc.x, ls1.d2.loc.x) && p1.loc.x <= max(ls1.d1.loc.x, ls1.d2.loc.x) &&
               p1.loc.y >= min(ls1.d1.loc.y, ls1.d2.loc.y) && p1.loc.y <= max(ls1.d1.loc.y, ls1.d2.loc.y)){
                return true;
               }else{
                return false;
               }
        }
    } 

    //=0, intersected; =1, parallel; =2 perpendicular; =3, identical
    pair<int,Vertex> lines_relationship(Line l1,Line l2)
    {   
        int relationship;
        Vertex intersection;
        //in the case of inifinely large slope x=1, and x=2
        if(l1.d1.loc.x==l1.d2.loc.x&&l2.d1.loc.x==l2.d2.loc.x){
            if(l1.d1.loc.x==l2.d1.loc.x){
                relationship=3;
            }else{
                relationship=1;
            }
        }
        else if(l1.d1.loc.x==l1.d2.loc.x){
            if(l2.d1.loc.y==l2.d2.loc.y){
                relationship=2;
                intersection.loc.x=l1.d1.loc.x;
                intersection.loc.y=l2.d1.loc.y;
            }else{
                relationship=0;
                l2.calc_slope_intercept();
                intersection.loc.x=l1.d1.loc.x;
                intersection.loc.y=intersection.loc.x*l2.slope+l2.intercept;
            }
        }
        else if(l2.d1.loc.x==l2.d2.loc.x){
            if(l1.d1.loc.y==l1.d2.loc.y){
                relationship=2;
                intersection.loc.x=l2.d1.loc.x;
                intersection.loc.y=l1.d1.loc.y;
            }else{
                relationship=0;
                l1.calc_slope_intercept();
                intersection.loc.x=l2.d1.loc.x;
                intersection.loc.y=intersection.loc.x*l1.slope+l1.intercept;
            }
        }
        else
        {
            l1.calc_slope_intercept();
            l2.calc_slope_intercept();
            if(l1.slope==l2.slope){
                if(l1.intercept==l2.intercept){
                    relationship=3;
                }
                else{
                    relationship=1;
                }
            }
            else if(l1.slope*l2.slope==-1){
                relationship = 2;
                intersection.loc.x = (l2.intercept-l1.intercept)/(l1.slope-l2.slope);
                intersection.loc.y = intersection.loc.x*l1.slope + l1.intercept;
            }
            else{
                relationship=0;
                intersection.loc.x = (l2.intercept-l1.intercept)/(l1.slope-l2.slope);
                intersection.loc.y = intersection.loc.x*l1.slope + l1.intercept;
            }
        }

        return make_pair(relationship,intersection);
    }

    //=0, intersected; =1, parallel; =2, perpendicular; =3, identical; =4, overlapping; =5, touching; =6, disjoint
    pair<int,Vertex> segments_relationship(Line ls1, Line ls2)
    {   
        int segments_rela;
        Vertex segments_intersection;
        pair<int, Vertex> lines_rela = line_geo::lines_relationship(ls1,ls2);
        
        //if the two lines are intersected
        if(lines_rela.first==0){
            //judge if the intersection still existed in the segments 
            if(line_segment_vertex_relationship(ls1,lines_rela.second)){
                //the intersection existed in the segments 
                if(lines_rela.second.same_vertex(ls1.d1)||lines_rela.second.same_vertex(ls1.d2)){
                    //the intersection is one of the endpoint => touching
                    segments_rela = 5;
                    segments_intersection = lines_rela.second;
                }else{
                    segments_rela = 0;
                    segments_intersection = lines_rela.second;
                }
            }else{
                //the intersection existed outside the segments
                segments_rela = 6;
            }
        } 
        else if(lines_rela.first==1){
        //if the two lines are parallel
            segments_rela = 1;
        }
        else if(lines_rela.first==2){
            //judge if the perpendicular intersection still existed in the segments 
            if(line_segment_vertex_relationship(ls1,lines_rela.second)){
                //the intersection existed in the segments 
                if(lines_rela.second.same_vertex(ls1.d1)||lines_rela.second.same_vertex(ls1.d2)){
                    //the intersection is one of the endpoint => touching
                    segments_rela = 5;
                    segments_intersection = lines_rela.second;
                }else{
                    //the perpendicular intersection still exists
                    segments_rela = 2;
                }
            }else{
                //the intersection existed outside the segments
                segments_rela = 6;
            }
        }
        else{
            //line segment is the same
            if(ls1.same_segment(ls2)){
                segments_rela = 3;
            }
            else if((ls1.d1.loc.x>max(ls2.d1.loc.x,ls2.d2.loc.x))&&(ls1.d2.loc.x>max(ls2.d1.loc.x,ls2.d2.loc.x))||
                   (ls1.d1.loc.x<min(ls2.d1.loc.x,ls2.d2.loc.x))&&(ls1.d2.loc.x<min(ls2.d1.loc.x,ls2.d2.loc.x))){
                segments_rela = 1;
            }
            else if(ls1.d1.same_vertex(ls2.d1)||ls1.d1.same_vertex(ls2.d2)){
                segments_rela = 5;
                segments_intersection = ls1.d1;
            }
            else if(ls1.d2.same_vertex(ls2.d1)||ls1.d2.same_vertex(ls2.d2)){
                segments_rela = 5;
                segments_intersection = ls1.d2;
            }
            else{
                segments_rela = 4;
            }

        }
   
        return make_pair(segments_rela,segments_intersection);
    }

   vector<Line*> region_partition_lines_EdU(vector<Vertex*> vv){
        vector<Line*> partition_lines;
        sort(vv.begin(),vv.end(),[](Vertex* va, Vertex* vb){
            return va->loc.x > va->loc.y;
        });

        for(int i=0;i<(int)vv.size();i++){
            if(i==0){
                Line* li = new Line;
                li->slope = -1.0/((vv[0]->loc.y-vv[1]->loc.y)/(vv[0]->loc.x-vv[1]->loc.x));
                li->intercept = vv[0]->loc.y-li->slope*vv[0]->loc.x;
                partition_lines.push_back(li);
            }
            else if(i=vv.size()-1){
                Line* li = new Line;
                li->slope = -1.0/((vv[i-1]->loc.y-vv[i]->loc.y)/(vv[i-1]->loc.x-vv[i]->loc.x));
                li->intercept = vv[i]->loc.y-li->slope*vv[i]->loc.x;
                partition_lines.push_back(li);
            }
            else{
                Line* li = new Line;
                li->slope = -1.0/((vv[i-1]->loc.y-vv[i+1]->loc.y)/(vv[i-1]->loc.x-vv[i+1]->loc.x));
                li->intercept = vv[i]->loc.y-li->slope*vv[i]->loc.x;
                partition_lines.push_back(li);
            }
        }
        return partition_lines;
   } 

    vector<double> angles_normalization_changing_axis(vector<Vertex*> vv_EdU, vector<Vertex*> vv_pairs){
        vector<Line*> partition_lines;
        vector<double> axis_normalizing_angles;
        sort(vv_pairs.begin(),vv_pairs.end(),[](Vertex* va, Vertex* vb){
            return va->loc.y < vb->loc.y;
        });
        //cout_fout_debug::cout_vector_vertex(vv_pairs);

        for(int i=0;i<(int)vv_pairs.size();i++){
            if(i==0){
                Line* li = new Line;
                li->slope = -1.0/((vv_pairs[0]->loc.y-vv_pairs[1]->loc.y)/(vv_pairs[0]->loc.x-vv_pairs[1]->loc.x));
                li->intercept = vv_pairs[0]->loc.y-li->slope*vv_pairs[0]->loc.x;
                partition_lines.push_back(li);
            }
            else if(i==vv_pairs.size()-1){
                Line* li = new Line;
                li->slope = -1.0/((vv_pairs[i-1]->loc.y-vv_pairs[i]->loc.y)/(vv_pairs[i-1]->loc.x-vv_pairs[i]->loc.x));
                li->intercept = vv_pairs[i]->loc.y-li->slope*vv_pairs[i]->loc.x;
                partition_lines.push_back(li);
            }
            else{
                Line* li = new Line;
                li->slope = -1.0/((vv_pairs[i-1]->loc.y-vv_pairs[i+1]->loc.y)/(vv_pairs[i-1]->loc.x-vv_pairs[i+1]->loc.x));
                li->intercept = vv_pairs[i]->loc.y-li->slope*vv_pairs[i]->loc.x;
                partition_lines.push_back(li);
            }
            double normalizing_angle = atan(-1.0/partition_lines[i]->slope)* 180.0 / M_PI;
            if(normalizing_angle<0){
                normalizing_angle+=180.0;
            }
            axis_normalizing_angles.push_back(normalizing_angle);
        }
        //cout_fout_debug::cout_vector_line(partition_lines);
        //cout_fout_debug::cout_vector_double(axis_normalizing_angles);
        

        vector<double> EdU_normalized_angles;
        vector<double> EdU_angles;
        for(int i=0;i<(int)vv_EdU.size()/2;i++){
            double EdU_angle = atan2((vv_EdU[2*i+1]->loc-vv_EdU[2*i]->loc).y,(vv_EdU[2*i+1]->loc-vv_EdU[2*i]->loc).x)*(180.0/M_PI);
            if(EdU_angle<0){
                //EdU_angle+=180;
            }
            EdU_angles.push_back(EdU_angle);
            double pairs_center_x = (vv_EdU[2*i+1]->loc+vv_EdU[2*i]->loc).x/2.0;
            double pairs_center_y = (vv_EdU[2*i+1]->loc+vv_EdU[2*i]->loc).y/2.0;
            int above_the_line_j = -1;
            for(int j=0;j<(int)partition_lines.size();j++){
                if(partition_lines[j]->slope*pairs_center_x+partition_lines[j]->intercept<pairs_center_y){
                    above_the_line_j=j;
                }
            }
            double EdU_normalized_angle;
            if(above_the_line_j!=-1){
                //the normalizing angles is from axis_point j
                EdU_normalized_angle = 90.0+EdU_angle-axis_normalizing_angles[above_the_line_j];
                std::cout<<"For EdU pairs "<<i<<" , its EdU_angle is "<<EdU_angle<<" and its axis_angle is "<<axis_normalizing_angles[above_the_line_j]<<"; after normalization by midvein, its angle is "<<EdU_normalized_angle<<endl;
            }
            else{
                EdU_normalized_angle = 90.0+EdU_angle-axis_normalizing_angles[(int)vv_pairs.size()-1];
                std::cout<<"For EdU pairs "<<i<<" , its EdU_angle is "<<EdU_angle<<" and its axis_angle is "<<axis_normalizing_angles[(int)vv_pairs.size()-1]<<"; after normalization by midvein, its angle is "<<EdU_normalized_angle<<endl;
            }
            EdU_normalized_angles.push_back(EdU_normalized_angle);
        }

        for(int i=0; i<EdU_normalized_angles.size();i++){
            if(EdU_normalized_angles[i]<0){
                EdU_normalized_angles[i]+=360.0;
            }
            if(EdU_normalized_angles[i]>180.0){
                EdU_normalized_angles[i]-=180.0;
            }
        }
        return EdU_normalized_angles;
        //return axis_normalizing_angles;
    }
}

namespace cell_geo{

    //calculate the center of cell ci in p_g
    _vec<double> cell_center(Organ* p_g, int ci){
        _vec<double> center_tmp = _vec<double>{0.0,0.0,0.0};
        for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
            center_tmp = center_tmp + p_g->p_v[p_g->p_c[ci]->vi[vi]]->loc;
        }
        center_tmp = center_tmp/(double)p_g->p_c[ci]->vi.size();
        p_g->p_c[ci]->center = center_tmp;
        return center_tmp;
    }

    //calculate the area of cell ci in p_g
    //assumptions: vertex indices are sorted in anticlockwise direction
    double cell_area(Organ* p_g, int ci){
        double i_area =0;
        Cell *cp = p_g->p_c[ci];

        //i番目の細胞の面積を計算。点が反時計回りに格納されていることを前提にしている。
        for (int j = 0; j < (int)cp->vi.size(); j++) {
        Vertex *vp[2];
        vp[0] = p_g->p_v[cp->vi[j]];
        if (j != (int)cp->vi.size() - 1) {
            vp[1] = p_g->p_v[cp->vi[j + 1]];
        }
        else if (j == (int)cp->vi.size() - 1) {
            vp[1] = p_g->p_v[cp->vi[0]];
        }
        else {
            std::cout << "Bug.Area" << std::endl;
            exit(0);
        }
        i_area += 0.5 * (vp[0]->loc.x * vp[1]->loc.y - vp[1]->loc.x * vp[0]->loc.y);
        }
        p_g->p_c[ci]->area=i_area;
        return i_area;
    }

    //calculate the perimeter of cell ci in p_g
    double cell_perimeter(Organ* p_g, int ci){

        double perimeter_tmp=0;
        for(int li=0;li<p_g->p_c[ci]->li.size();li++){
            double length_tmp = (p_g->p_v[p_g->p_l[p_g->p_c[ci]->li[li]]->vi[0]]->loc-p_g->p_v[p_g->p_l[p_g->p_c[ci]->li[li]]->vi[1]]->loc).norm();
            perimeter_tmp += length_tmp;
        }
        p_g->p_c[ci]->perimeter = perimeter_tmp;
        //std::cout<<"Cell "<<ci<<" perimeter: "<<perimeter_tmp<<endl;
        return perimeter_tmp;
    }

    //calculate the regularity of cell ci in p_g
    //regularity defines the similarity between a defined polygon (eg., cell ci) and the regular polygon that has the same perimeter
    //regularity ranges from (0,1.0], when the polygon is a regular polygon, regularity = 1.0;
    //regularity = area (cell ci)/area (the regularity polygon that has the same perimeter) 
    //reference: https://doi.org/10.1016/j.cad.2012.07.012
    double cell_regularity(Organ* p_g, int ci){
        double area_P = p_g->p_c[ci]->area;
        double perimeter_P = p_g->p_c[ci]->perimeter;
        double area_R = perimeter_P*perimeter_P/(4.0*(double)p_g->p_c[ci]->li.size()*tan(3.145926/(double)p_g->p_c[ci]->li.size()));
        p_g->p_c[ci]->regularity = area_P/area_R;
        //std::cout<<area_R<<" "<<perimeter_P<<" "<<p_g->p_c[ci]->li.size()<<" "<<tan(3.145926/(double)p_g->p_c[ci]->li.size())<<endl;
        //std::cout<<"Cell "<<ci<<" regularity "<<p_g->p_c[ci]->regularity<<endl;
        return p_g->p_c[ci]->regularity;
    }

    //sort lines in counterclockwise direction, based on vertices in counterclockwise direction
    vector<int> cell_counterclock_line(Organ* p_g, int ci){
        std::vector<int> anticlockwise_lidx;
        Cell *cp = p_g->p_c[ci];
        int cp_vsize = cp->vi.size();
        for(int i = 0; i<cp_vsize;++i){
            for(int lidx: cp->li) {
            if(p_g->p_l[lidx]->vi[0] == cp->vi[i] && p_g->p_l[lidx]->vi[1]==cp->vi[(i+1)%cp_vsize]){
                anticlockwise_lidx.push_back(lidx);
                break;
            }else if(p_g->p_l[lidx]->vi[0] == cp->vi[(i+1)%cp_vsize] && p_g->p_l[lidx]->vi[1] == cp->vi[i]){
                std::swap(p_g->p_l[lidx]->vi[0], p_g->p_l[lidx]->vi[1]);
                anticlockwise_lidx.push_back(lidx);
            }
            }
        }
        return anticlockwise_lidx;
    }

}

namespace organ_geo{

    //calculate the length of all lines within an organ
    void organ_line_length(Organ* p_g){
        for(int li=0;li<(int)p_g->p_l.size();li++){
            line_geo::line_length(p_g,li);
        }
    }

    //calculate the perimeter of all cells within an organ
    double organ_cell_perimeter(Organ* p_g){
        double av_perimeter_tmp=0;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            av_perimeter_tmp +=cell_geo::cell_perimeter(p_g,ci);
        }
        av_perimeter_tmp = av_perimeter_tmp/(double)p_g->p_c.size();
        return av_perimeter_tmp;
    }

    //calculate the center of the organ: sum of the center of all cells
    _vec<double> organ_center(Organ* p_g){
        _vec<double> center_tmp = _vec<double>{0.0,0.0,0.0};
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            center_tmp += cell_geo::cell_center(p_g,ci);
        }
        p_g->center = center_tmp/(double)p_g->p_c.size();
        return center_tmp;
    }

    //calculate the area of the organ: sum of the area of all cells
    double organ_area(Organ* p_g){
        double area_tmp=0.0;
        double epiArea=0.0, inArea=0.0;
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            area_tmp+=cell_geo::cell_area(p_g,ci);
            if(p_g->p_c[ci]->IsEpidermal==1){
                epiArea+=p_g->p_c[ci]->area;
            }
            else{
                inArea+=p_g->p_c[ci]->area;
            }
        }
        p_g->area=area_tmp;
        p_g->area_averaged = p_g->area/(double)p_g->p_c.size();
        p_g->epiArea = epiArea;
        p_g->epiArea_averaged = epiArea/p_g->N_epi_cell;
        p_g->inArea = inArea;
        p_g->inArea_averaged = inArea/p_g->N_inner_cell; 
        return area_tmp;
    }

    //determine the epidermal identity of all cells and lines within an organ: 
    //cell->IsSurface=1, epidermal cell;=0, inner cell; Line->IsOutermost=1, epidermal line;=0, inner line; also count the epidermal cell number and inner cell number
    void epidermal_identity(Organ* p_g){

        vector<int> surface_line;
        vector<int> surface_vertex;

        for(int vi=0; vi<(int)p_g->p_v.size();vi++){
            p_g->p_v[vi]->occurrenceInCell=0;
        }


        //count the occurrence of each vertex in all cells
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                int tmp = p_g->p_c[ci]->vi[vi];
                p_g->p_v[tmp]->occurrenceInCell++;
            }
        }

        //find surface vertex: if the occurrence of vertex is less than 3, then this vertex is a surface vertex
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->occurrenceInCell<3){
                p_g->p_v[vi]->IsSurface=1;
                surface_vertex.push_back(vi);
            }
            else{
                p_g->p_v[vi]->IsSurface=0;
            }
        }
        //Debug: check the occurrence of vertex in cell and judgement of surface vertex
        /*
        std::cout<<"Debug: check the occurrence of vertex in cell and judgement of surface vertex"<<std::endl;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            std::cout<<vi<<" "<<p_g->p_v[vi]->occurrenceInCell<<" "<<p_g->p_v[vi]->IsSurface<<std::endl;
        }
        std::cout<<"End Debug: check the occurrence of vertex in cell and judgement of surface vertex"<<std::endl;
        */

    //find outermost edge: if the two vertices connecting the edge are surface vertices, then this edge is an outermost edge
    for(int li=0; li<(int)p_g->p_l.size();li++){
        if(p_g->p_v[p_g->p_l[li]->vi[0]]->IsSurface==1 && p_g->p_v[p_g->p_l[li]->vi[1]]->IsSurface==1){
            p_g->p_l[li]->IsOutermost=1;
            surface_line.push_back(li);
        }
        else{
            p_g->p_l[li]->IsOutermost=0;
        }
    }
    //Debug: check each line's IsOutermost parameter
    /*
    std::cout<<"Debug: check IsOutermost parameter for Line"<<std::endl;
    for(int li=0;li<(int)p_g->p_l.size();li++){
        std::cout<<li<<" "<<p_g->p_l[li]->IsOutermost<<std::endl;
    }
    std::cout<<"End Debug: check IsOutermost parameter for Line"<<std::endl;
        */

        //find epidermal cell: if a cell have an outermost edge, then it is an epidermal cell
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            bool IsEpidermal_tmp=0;
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                if(p_g->p_l[p_g->p_c[ci]->li[li]]->IsOutermost==1){
                    IsEpidermal_tmp=1;
                    p_g->p_c[ci]->surfaceVertex[0]=p_g->p_l[p_g->p_c[ci]->li[li]]->vi[0];
                    p_g->p_c[ci]->surfaceVertex[1]=p_g->p_l[p_g->p_c[ci]->li[li]]->vi[1];
                    p_g->p_c[ci]->outermostLength = p_g->p_l[p_g->p_c[ci]->li[li]]->length;
                }
            }
            p_g->p_c[ci]->IsEpidermal=IsEpidermal_tmp;
        }

        int inner_cell_number=0, peripheral_cell_number=0;

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->IsEpidermal==0){
                inner_cell_number++;
            }
            else{
                peripheral_cell_number++;
            }
        }
        
        //std::cout<<"Currently, we have "<<inner_cell_number<<" inner cells, and "<<peripheral_cell_number<<" peripheral cells."<<std::endl;
        p_g->N_epi_cell = peripheral_cell_number;
        p_g->N_inner_cell = inner_cell_number;
        p_g->surface_line=surface_line;
        p_g->surface_vertex=surface_vertex;
    }

    //calculate the perimeter of the organ: sum of the length of all epidermal line
    double organ_perimeter(Organ* p_g){
        double perimeter_tmp=0.0;
        for(int li=0; li<(int)p_g->p_l.size();li++){
        if(p_g->p_l[li]->IsOutermost==1){
                perimeter_tmp += p_g->p_l[li]->length;
            }
        }
        p_g->perimeter = perimeter_tmp;
        p_g->perimeter_averaged = perimeter_tmp/(double)p_g->p_c.size();
        return perimeter_tmp;
    }

    //calculate the regularity of the organ: 
    //regularity defines the similarity between a defined polygon (eg., p_g) and the regular polygon that has the same perimeter
    //regularity ranges from (0,1.0], when the polygon is a regular polygon, regularity = 1.0;
    //regularity = area (cell ci)/area (the regularity polygon that has the same perimeter) 
    //reference: https://doi.org/10.1016/j.cad.2012.07.012
    double organ_regularity(Organ* p_g){
        double regularity_averaged_tmp=0, regularity_in_tmp=0, regularity_epi_tmp=0;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->regularity = cell_geo::cell_regularity(p_g,ci);
            //std::cout<<"Cell "<<ci<<" regularity is "<<p_g->p_c[ci]->regularity<<std::endl;
            regularity_averaged_tmp += p_g->p_c[ci]->regularity;
            if(p_g->p_c[ci]->IsEpidermal==0){
                regularity_in_tmp +=p_g->p_c[ci]->regularity;
            }
            else{
                regularity_epi_tmp +=p_g->p_c[ci]->regularity;
            }
        }

        p_g->regularity_averaged = regularity_averaged_tmp/p_g->p_c.size();
        p_g->regularity_in_av=regularity_in_tmp/p_g->N_inner_cell;
        p_g->regularity_epi_av = regularity_epi_tmp/p_g->N_epi_cell;
        //std::cout<<"regularity_averaged_tmp: "<<regularity_averaged_tmp<<std::endl;
        return regularity_averaged_tmp;
    }

    //calculate the circularity of the organ: circularity quantifies the similarity between the defined polygon (eg., p_g) and a circle
    //circularity ranges from (0,1], when circularity=1, the polygon is a perfect circle; also could be refered as roundness
    //circularity = Perimeter^2/(4pi*Area)
    //reference: https://en.wikipedia.org/wiki/Roundness
    double organ_circularity(Organ* p_g){
        p_g->circularity = 4*3.1415926*p_g->area/(p_g->perimeter*p_g->perimeter);
        return p_g->circularity;
    }

    //classify cells into different layers; eg., the epidermal cells belong to the layer 0; the subepidermal cells belong to the layer 1,...... 
    void organ_cell_layer(Organ* p_g){
        int cell_layer_number=0;
    //0. clear all the vi->layer vector and previous labelling
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            p_g->p_v[vi]->layer.clear();
        }
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->layer=-1;
        }

    //1. epidermal cells layer 0, make all vertices in epidermal cells to be layer 0
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            if(p_g->p_c[ci]->IsEpidermal==1){
                p_g->p_c[ci]->layer=0;
                for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                    p_g->p_v[p_g->p_c[ci]->vi[vi]]->layer.push_back(0);
                }
            }
        }
    //2. find inner layers step by step
        int now_checking_layer = 1;
        bool checking_layer_completed = 0;
        do{
            //layer 1

            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                    //cells containing vertices that belongs to layer 0 && cells which are not epidermal 
                    std::vector<int>::iterator iter = std::find(p_g->p_v[p_g->p_c[ci]->vi[vi]]->layer.begin(),p_g->p_v[p_g->p_c[ci]->vi[vi]]->layer.end(),now_checking_layer-1);
                        if(iter != p_g->p_v[p_g->p_c[ci]->vi[vi]]->layer.end()&&p_g->p_c[ci]->layer==-1)
                        {
                            p_g->p_c[ci]->layer = now_checking_layer;
                            break;
                        }
                    }
                if(p_g->p_c[ci]->layer==now_checking_layer){
                    for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                        p_g->p_v[p_g->p_c[ci]->vi[vi]]->layer.push_back(now_checking_layer);
                    }
                }
            }
            now_checking_layer++;
            checking_layer_completed=1;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                if(p_g->p_c[ci]->layer==-1){
                    checking_layer_completed=0;
                    cell_layer_number=now_checking_layer-1;
                }
            }

        }while(checking_layer_completed==0);
        p_g->cell_layer_number=cell_layer_number;
        //debug
        /*
        std::cout<<"Layer function degbug"<<std::endl;
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            std::cout<<ci<<" "<<p_g->p_c[ci]->layer<<std::endl;
        }
        */
    }
    

    void organ_vertex_counterclockwise_sort(Organ* p_g){
        for(int cidx=0; cidx < (int)p_g->p_c.size();++cidx){
            std::vector<int> tmp_vi;
            std::vector<int> tmp_li;

            int prev_vidx = -1, curr_vidx = p_g->p_c[cidx]->vi[0];
            tmp_vi.push_back(curr_vidx);
            while(tmp_li.size() < p_g->p_c[cidx]->li.size()) {
                for(int lidx: p_g->p_c[cidx]->li) {
                    if(p_g->p_l[lidx]->vi[0] == curr_vidx && p_g->p_l[lidx]->vi[1] != prev_vidx) {
                    prev_vidx = curr_vidx;
                    curr_vidx = p_g->p_l[lidx]->vi[1];
                    tmp_li.push_back(lidx);
                    if(curr_vidx != tmp_vi[0]){
                        tmp_vi.push_back(curr_vidx);
                    }
                        }else if(p_g->p_l[lidx]->vi[0] != prev_vidx && p_g->p_l[lidx]->vi[1] == curr_vidx) {
                        prev_vidx = curr_vidx;
                        curr_vidx = p_g->p_l[lidx]->vi[0];
                        tmp_li.push_back(lidx);
                        if(curr_vidx != tmp_vi[0]){
                            tmp_vi.push_back(curr_vidx);
                        }
                    }
                }
            }

            assert(tmp_vi.size() == p_g->p_c[cidx]->vi.size());
            assert(tmp_li.size() == p_g->p_c[cidx]->li.size());

            double area = 0.0;
            for(int i = 0; i < (int)tmp_vi.size(); ++i) {
            _vec<double> r1 = p_g->p_v[tmp_vi[i]]->loc;
            _vec<double> r2 = p_g->p_v[tmp_vi[(i + 1) % tmp_vi.size()]]->loc;
            area += 0.5 * (r1 % r2).z;
            }
            if(area < 0) {
            std::reverse(tmp_vi.begin(), tmp_vi.end());
            std::reverse(tmp_li.begin(), tmp_li.end());
            }

            for(int i = 0; i < (int)tmp_vi.size(); ++i) {
            p_g->p_c[cidx]->vi[i] = tmp_vi[i];
            p_g->p_c[cidx]->li[i] = tmp_li[i];
            }
        }

    }

    //calculate the length of organ: find the longest distance for the pairs of epidermal vertices
    double organ_length(Organ* p_g){
        
        //partition of the organ contour into 10,000 points
        int partition_number_for_length = 10000;
        vector<Vertex> anticlockwise_surface_vertex = organ_geo::organ_ordered_anticlockwise_boundary(p_g).vi;
        vector<Vertex> partitioned_contour_for_length = organ_geo::organ_boundary_points_along_polygon(p_g,anticlockwise_surface_vertex,partition_number_for_length);
        
        double organ_length = geo_vv::maxmum_distance_between_points_rotating_caliphers(partitioned_contour_for_length,p_g);
        p_g->organ_length = organ_length;
        double length_axis_slope = (p_g->length_vertex.second.y-p_g->length_vertex.first.y)/(p_g->length_vertex.second.x-p_g->length_vertex.first.x);
        p_g->length_axis_slope = length_axis_slope;
        double length_axis_intercept = p_g->length_vertex.second.y-p_g->length_vertex.second.x*length_axis_slope;
        p_g->length_axis_intercept = length_axis_intercept;
        
        return organ_length;
    }
    
    double organ_width(Organ* p_g){
        double organ_width=0;
        
        //1. calculate the width axis slope; width axis should be perpendicular to the length axis
        double length_axis_slope = p_g->length_axis_slope;
        double width_axis_slope = -1.0/length_axis_slope;
        //std::cout<<"the length axis is made by ("<<p_g->length_vertex.first.x<<","<<p_g->length_vertex.first.y<<");("<<p_g->length_vertex.second.x<<","<<p_g->length_vertex.second.y<<"); length axis slope: "<<length_axis_slope<<"; width axis slope: "<<width_axis_slope<<endl;
        //2. find intersections between the line (slope known, and intercept from the top to the base) and organ's outline
        //2.1 divide the organ by length axis into 10,000 parts
        int length_axis_precision = 10000;
        for(int i=0;i<length_axis_precision;i++){
            double length_axis_y_min = min(p_g->length_vertex.second.y,p_g->length_vertex.first.y);
            double length_axis_y_max = max(p_g->length_vertex.second.y,p_g->length_vertex.first.y);
            double length_axis_x_max = max(p_g->length_vertex.second.x,p_g->length_vertex.first.x);
            double length_axis_x_min = min(p_g->length_vertex.second.x,p_g->length_vertex.first.x);
            
            double division_point_x = length_axis_x_min + (double)i/(double)length_axis_precision*(length_axis_x_max-length_axis_x_min);
            double division_point_y = length_axis_y_min + (double)i/(double)length_axis_precision*(length_axis_y_max-length_axis_y_min);

            double line_intercept = division_point_y - width_axis_slope*division_point_x;
        //std::cout<<i<<" division_point_x "<<division_point_x<<", division_point_y "<<division_point_y<<"; width_axis_intercept: "<<division_point_y - width_axis_slope*division_point_x<<endl;
            vector<_vec<double>> intersection_line_organ_outline = line_polygon_intersection(width_axis_slope, line_intercept,p_g);
            double organ_width_tmp=0;
            if(intersection_line_organ_outline.size() == 2){
                organ_width_tmp = (intersection_line_organ_outline[0]-intersection_line_organ_outline[1]).norm();
                //std::cout<<i<<", organ_width_tmp: "<<organ_width_tmp<<endl;
                
            }

            if(organ_width_tmp>organ_width){
                organ_width = organ_width_tmp;
                p_g->width_vertex=std::make_pair(intersection_line_organ_outline[0],intersection_line_organ_outline[1]);
            }
        }
        //std::cout<<"organ_width: "<<organ_width<<endl;
        p_g->organ_width = organ_width;
        p_g->width_axis_slope = (p_g->width_vertex.second.y-p_g->width_vertex.first.y)/(p_g->width_vertex.second.x-p_g->width_vertex.first.x);
        p_g->width_axis_intercept = (p_g->width_vertex.second.y-p_g->width_vertex.second.x*p_g->width_axis_slope); 
        p_g->length_width_ratio = p_g->organ_length/p_g->organ_width;
        return organ_width;
    }

    //calculate the elliptical index of organ: 1-2*sqrt(pi)*sqrt(A_ch)/C_ch   ch means convex hull, in vertex model cases, if there is no overlap, the boundary should always be a convex hull
    //Mittenentzwei, Sarah, et al. "Definition and extraction of 2D shape indices of intracranial aneurysm necks for rupture risk assessment." International Journal of Computer Assisted Radiology and Surgery 16.11 (2021): 1977-1984.
    double elliptical_index(Organ* p_g){
        double elliptical_index_tmp = 1 - 2*sqrt(M_PI*p_g->area)/p_g->perimeter;
        p_g->elliptical_index = elliptical_index_tmp;
        return elliptical_index_tmp;
        //reference: https://link.springer.com/article/10.1007/s11548-021-02469-z
    }

    //calculate the geometric entropy of leaf organ: S_L = 1.0/4.0*(P/sqrt(A))
    double geometric_entropy(Organ* p_g){
        double geometric_entropy_tmp = 0.25*p_g->perimeter/sqrt(p_g->area);
        return geometric_entropy_tmp;
    }

    //calculate the intersections between line and p_g
    vector<_vec<double>> line_polygon_intersection(double slope, double intercept, Organ* p_g){
        vector<_vec<double>> intersection_results;
        //to find crosspoints between a line and a polygon
        //we can try to find crosspoints between a line and a line segment of the polygon surface
        
        //std::cout<<"surface_line_size "<<p_g->surface_line.size()<<endl;

        for(int li=0;li<(int)p_g->surface_line.size();li++){
            _vec<double> intersection_result = line_geo::line_cell_wall_intersection(slope, intercept, p_g, p_g->surface_line[li]);
            //std::cout<<"p_g->surface_line[li]: "<<p_g->surface_line[li]<<endl;
            if(intersection_result.z==0.0){
                intersection_results.push_back(intersection_result);
            }
        }
        /*
        std::cout<<"intersection_results: "<<endl;
        for(int i=0;i<intersection_results.size();i++){
            std::cout<<i<<" "<<intersection_results[i].x<<" "<<intersection_results[i].y<<" "<<intersection_results[i].z<<endl;
        }
        std::cout<<endl;
        */
       /*
        vector<_vec<double>> intersection_results_real;
        if(intersection_results.size()!=0){
            intersection_results_real.push_back(intersection_results[0]);
            for(int i=1;i<intersection_results.size();i++){
                for(int j=0;j<intersection_results_real.size();j++){
                    if(intersection_results[i].x!=intersection_results_real[j].x&&intersection_results[i].y!=intersection_results_real[j].y){
                        intersection_results_real.push_back(intersection_results[i]);
                        goto endloop_intersection_result_real;
                    }
                    else{
                        //std::cout<<"the same intersections found!"<<endl;
                    }
                }
                endloop_intersection_result_real:;
            }
        }
        */
        /*
        std::cout<<"intersection_results_real: ";
        for(int i=0;i<intersection_results_real.size();i++){
            std::cout<<i<<" "<<intersection_results_real[i].x<<" "<<intersection_results_real[i].y<<" "<<intersection_results_real[i].z<<endl;
        }
        std::cout<<endl;
        */
        return intersection_results;
    }

    //calculate the organ_leaf_index = organ->length/organ->width
    double organ_leaf_index(Organ* p_g){
        double leaf_index_tmp=0;
        organ_geo::organ_length(p_g);
        organ_geo::organ_width(p_g);

        if(abs(p_g->length_axis_slope)>1){
            //abs(length slope) > 1 means that length axis is vertical
        }
        else{
            //abs(length slope) < 1 means that length axis is horizontal, which means it is actually the width
            double width_tmp = p_g->organ_length;
            double length_tmp = p_g->organ_width;
            p_g->organ_length = length_tmp;
            p_g->organ_width = width_tmp;
            std::cout<<"length width has been switched"<<endl;
        }
        leaf_index_tmp = p_g->organ_length/p_g->organ_width;
        p_g->leaf_index = leaf_index_tmp;
        std::cout<<"leaf_index: "<<p_g->leaf_index<<endl;
        return leaf_index_tmp;
    }

    //calculate the cell arrangement of organ
    _vec<double> organ_cell_arrangement(Organ* p_g){
        double cell_arrangement_x =0;
        double cell_arrangement_y =0;
        _vec<double> cell_arrangement = _vec<double>{0.0,0.0,0.0};

        double length_axis_slope=p_g->length_axis_slope;
        double length_axis_intercept=p_g->length_axis_intercept;
        double width_axis_slope=p_g->width_axis_slope;
        double width_axis_intercept=p_g->width_axis_intercept;
        //calculate how many cells are crossed by length_axis
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                _vec<double> intersection_tmp=line_geo::line_cell_wall_intersection(length_axis_slope,length_axis_intercept,p_g,p_g->p_c[ci]->li[li]);
                if(intersection_tmp.z==0){
                    //they do have an intersection
                    cell_arrangement_x++;
                    goto length_axis_cell_count;
                }
                
            }
            length_axis_cell_count:;
        }

        //calculate how many cells are crossed by width axis
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                _vec<double> intersection_tmp=line_geo::line_cell_wall_intersection(width_axis_slope,width_axis_intercept,p_g,p_g->p_c[ci]->li[li]);
                if(intersection_tmp.z==0){
                    //they do have an intersection
                    cell_arrangement_y++;
                    goto width_axis_cell_count;
                }
                
            }
            width_axis_cell_count:;
        }
        if(abs(p_g->length_axis_slope)>1){
            //abs(length slope) > 1 means that length axis is vertical
        }
        else{
            double cell_arrangement_x_tmp = cell_arrangement_y;
            double cell_arrangement_y_tmp = cell_arrangement_x;
            cell_arrangement_x = cell_arrangement_x_tmp;
            cell_arrangement_y = cell_arrangement_y_tmp;
        }
        p_g->cell_arrangement_x=cell_arrangement_x;
        p_g->cell_arrangement_y=cell_arrangement_y;
        p_g->cell_arrangement_ratio = cell_arrangement_x/cell_arrangement_y;

        //std::cout<<"cell_arrangement_x: "<<cell_arrangement_x<<" ;cell_arrangment_y: "<<cell_arrangement_y<<";cell_arrangment_ratio: "<<p_g->cell_arrangement_ratio<<endl;
        cell_arrangement.x = cell_arrangement_x;
        cell_arrangement.y = cell_arrangement_y;
        cell_arrangement.z = cell_arrangement_x/cell_arrangement_y;
        return cell_arrangement;
    }

    //calculate the potential energy of whole organ
    double organ_potential_energy(Organ* p_g){
        double potential_energy_all=0;
        for(int li=0;li<(int)p_g->p_l.size();li++){
            
            if(p_g->p_l[li]->IsOutermost==1){
                potential_energy_all += sigma_O*p_g->p_l[li]->length;
                //std::cout<<"potential_energy_all: "<<potential_energy_all<<endl;
            }
            else{
                potential_energy_all += sigma_L*p_g->p_l[li]->length;
            }
        }

        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            potential_energy_all += kappa_S*(p_g->p_c[ci]->area-S_std)*(p_g->p_c[ci]->area-S_std);
        }
        //std::cout<<"potential_energy_all: "<<potential_energy_all<<endl;
        p_g->organ_potential_energy = potential_energy_all;
        return potential_energy_all;
    }

    //calculate the maxmimum and minimum x,y of the organ
    //require organ_geo::organ_center()
    void organ_max_min_x_y(Organ* p_g){
        double y_min_vertex=p_g->p_v[0]->loc.y;
        double x_min_vertex=p_g->p_v[0]->loc.x;
        double y_max_vertex=p_g->p_v[0]->loc.y;
        double x_max_vertex=p_g->p_v[0]->loc.x;

        double y_min_cell=p_g->p_c[0]->center.y;
        double x_min_cell=p_g->p_c[0]->center.x;
        double y_max_cell=p_g->p_c[0]->center.y;
        double x_max_cell=p_g->p_c[0]->center.x;

        for(int vi=0; vi<(int)p_g->p_v.size(); vi++){
            if(y_min_vertex>p_g->p_v[vi]->loc.y){
                y_min_vertex = p_g->p_v[vi]->loc.y;
            }
            if(y_max_vertex<p_g->p_v[vi]->loc.y){
                y_max_vertex = p_g->p_v[vi]->loc.y;
            }
            if(x_min_vertex>p_g->p_v[vi]->loc.x){
                x_min_vertex = p_g->p_v[vi]->loc.x;
            }
            if(x_max_vertex<p_g->p_v[vi]->loc.x){
                x_max_vertex = p_g->p_v[vi]->loc.x;
            }
        }

        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            if(y_min_cell>p_g->p_c[ci]->center.y){
                y_min_cell = p_g->p_c[ci]->center.y;
            }
            if(y_max_cell<p_g->p_c[ci]->center.y){
                y_max_cell = p_g->p_c[ci]->center.y;
            }
            if(x_min_cell>p_g->p_c[ci]->center.x){
                x_min_cell = p_g->p_c[ci]->center.x;
            }
            if(x_max_cell<p_g->p_c[ci]->center.x){
                x_max_cell = p_g->p_c[ci]->center.x;
            }
        }

        p_g->y_min_vertex = y_min_vertex;
        p_g->y_max_vertex = y_max_vertex;
        p_g->x_min_vertex = x_min_vertex;
        p_g->x_max_vertex = x_max_vertex;

        p_g->y_min_cell = y_min_cell;
        p_g->y_max_cell = y_max_cell;
        p_g->x_min_cell = x_min_cell;
        p_g->x_max_cell = x_max_cell;

        //std::cout<<"y_min_vertex: "<<p_g->y_min_vertex<<"; y_max_vertex: "<<p_g->y_max_vertex<<"; x_min_vertex: "<<p_g->x_min_vertex<<"; x_max_vertex: "<<p_g->x_max_vertex<<std::endl;;
        //std::cout<<"y_min_cell: "<<p_g->y_min_cell<<"; y_max_cell: "<<p_g->y_max_cell<<"; x_min_cell: "<<p_g->x_min_cell<<"; x_max_cell: "<<p_g->x_max_cell<<std::endl;
    }

    //order all epidermal vertices (boundary vertices)/ boundary lines into a anticlockwise vector
    Ordered_boundary organ_ordered_anticlockwise_boundary(Organ* p_g){
        //std::cout<<"Start ordering surface vertices"<<endl;
        Ordered_boundary Result;
        vector<Vertex> surface_vertices;
        vector<Vertex> anticlockwise_surface_vertices;
        vector<Line> surface_line;
        vector<int> anticlockwise_surface_lines;
        
    //0. collecting all surface lines into surface_line, all surface vertices into surface_vertices
        for(int li=0;li<(int)p_g->p_l.size();li++){
            if(p_g->p_l[li]->IsOutermost==1){
                surface_line.push_back(*p_g->p_l[li]);
                //std::cout<<"p_g->p_l[li].li "<<surface_line[surface_line.size()].li<<std::endl;
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                surface_vertices.push_back(*p_g->p_v[vi]);
                //std::cout<<"p_g->p_v[vi].vi "<<surface_vertices[surface_vertices.size()].vi<<std::endl;
            }
        }
    //1. the first surface vertex: the surface vertex with lowest y position
        Vertex first_boundary_vertex;
        first_boundary_vertex = surface_vertices[0];

        for(int vi=0;vi<(int)surface_vertices.size();vi++){
            if(first_boundary_vertex.loc.y>surface_vertices[vi].loc.y){
                first_boundary_vertex=surface_vertices[vi];
            }
        }

        anticlockwise_surface_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<anticlockwise_surface_vertices[0].vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[0].loc.x<<","<<anticlockwise_surface_vertices[0].loc.y<<")."<<endl;
    //2. the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex v_tmp1;
        Vertex v_tmp2;
        int li_tmp1,li_tmp2;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)surface_line.size(); li++){
                //std::cout<<"li "<<surface_line[li].li<<" connecting vertex: "<<surface_line[li].vi[0]<<","<<surface_line[li].vi[1]<<"; size of surface line"<<surface_line.size()<<endl;
                if(surface_line[li].vi[0]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[1]]->loc.x: "<<p_g->p_v[surface_line[li].vi[1]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[surface_line[li].vi[1]];
                    li_tmp1=li;
                    v_tmp1_found=1;
                    //if(p_g->p_v[surface_line[li]->vi[1]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[1]]);
                    //}
                }
                else if(surface_line[li].vi[0]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==1){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[1]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp2 = *p_g->p_v[surface_line[li].vi[1]];
                    li_tmp2=li;
                    break;
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[0]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[surface_line[li].vi[0]];
                    li_tmp1=li;
                    v_tmp1_found=1;
                    //if(p_g->p_v[surface_line[li]->vi[0]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[0]]);
                    //}
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[0].vi&&v_tmp1_found==1){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[surface_line[li]->vi[0]]->loc.x: "<<p_g->p_v[surface_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp2 = *p_g->p_v[surface_line[li].vi[0]];
                    li_tmp2=li;
                    break;
                    //if(p_g->p_v[surface_line[li]->vi[0]]->loc.x<anticlockwise_surface_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    anticlockwise_surface_vertices.push_back(p_g->p_v[surface_line[li]->vi[0]]);
                    //}
                }
                else{

                }
            }
            if(v_tmp1.loc.x < v_tmp2.loc.x){
                anticlockwise_surface_vertices.push_back(v_tmp1);
                anticlockwise_surface_lines.push_back(surface_line[li_tmp1].li);
                //std::cout<<"v_tmp1 "<<v_tmp1<<endl;
            }
            else{
                anticlockwise_surface_vertices.push_back(v_tmp2);
                anticlockwise_surface_lines.push_back(surface_line[li_tmp2].li);
                //std::cout<<"v_tmp2 "<<v_tmp2<<endl;
            }
        //std::cout<<" The secondary boundary vertices is "<<anticlockwise_surface_vertices[1]->vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[1]->loc.x<<","<<anticlockwise_surface_vertices[1]->loc.y<<")."<<endl;
    //3. the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
        for(int vi=1; vi<(int)surface_line.size()-1;vi++){
            for(int li=0;li<(int)surface_line.size();li++){
                if(surface_line[li].vi[0]==anticlockwise_surface_vertices[vi].vi&&surface_line[li].vi[1]!=anticlockwise_surface_vertices[vi-1].vi){
                    anticlockwise_surface_vertices.push_back(*p_g->p_v[surface_line[li].vi[1]]);
                    anticlockwise_surface_lines.push_back(surface_line[li].li);
                }
                else if(surface_line[li].vi[1]==anticlockwise_surface_vertices[vi].vi&&surface_line[li].vi[0]!=anticlockwise_surface_vertices[vi-1].vi){
                    anticlockwise_surface_vertices.push_back(*p_g->p_v[surface_line[li].vi[0]]);
                    anticlockwise_surface_lines.push_back(surface_line[li].li);
                }
            }
        }
        //std::cout<<" The third boundary vertices is "<<anticlockwise_surface_vertices[2]->vi<<", and its (x,y) is ("<<anticlockwise_surface_vertices[2]->loc.x<<","<<anticlockwise_surface_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(anticlockwise_surface_vertices);
        Result.li = anticlockwise_surface_lines;
        Result.vi = anticlockwise_surface_vertices;
        return Result;
    }

    //find boundary points with equal distance on a close contour along polygon
    vector<Vertex*> organ_boundary_points_along_polygon(Organ* p_g, vector<Vertex*> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex*> boundary_points;
        double boundary_points_average_distance = p_g->perimeter/(double)boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = anticlockwise_surface_vertices[vj]->loc+(anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }

    vector<Vertex> organ_boundary_points_along_polygon(Organ* p_g, vector<Vertex> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex> boundary_points;
        double boundary_points_average_distance = p_g->perimeter/(double)boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex boundary_points_new;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new.loc = boundary_points[pi].loc+(anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new.loc = anticlockwise_surface_vertices[vj].loc+(anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }


    vector<Vertex*> organ_ordered_boundary_points_finding_pointer(Organ* p_g, int boundary_points_number){
        //std::cout<<"Start boundary points finding"<<endl;
        double boundary_points_average_distance = p_g->perimeter/boundary_points_number; 
        //std::cout<<"Organ perimeter "<<p_g->perimeter<<"; boundary_points: "<<boundary_points_number<<"; the averaged distance between boundary_points: "<<boundary_points_average_distance<<endl;

        vector<Vertex*> boundary_vertices;
        vector<Vertex*> boundary_points;
        vector<Line*> boundary_line;

            for(int li=0;li<(int)p_g->p_l.size();li++){
                if(p_g->p_l[li]->IsOutermost==1){
                    boundary_line.push_back(p_g->p_l[li]);
                }
            }

        //we need first to arrange the surface vertices into an anticlockwisely ordered array

        //the first boundary vertex: the surface vertex with lowest y position
        Vertex* first_boundary_vertex = new Vertex;

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                first_boundary_vertex = p_g->p_v[vi];
                break;
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                if(first_boundary_vertex->loc.y>p_g->p_v[vi]->loc.y){
                    first_boundary_vertex = p_g->p_v[vi];
                }
            }
        }
        
        boundary_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<boundary_vertices[0]->vi<<", and its (x,y) is ("<<boundary_vertices[0]->loc.x<<","<<boundary_vertices[0]->loc.y<<")."<<endl;

    //the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex* v_tmp1 = new Vertex;
        Vertex* v_tmp2 = new Vertex;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)boundary_line.size(); li++){
                if(boundary_line[li]->vi[0]==boundary_vertices[0]->vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[1]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp1 = p_g->p_v[boundary_line[li]->vi[1]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[1]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                    //}
                }
                else if(boundary_line[li]->vi[0]==boundary_vertices[0]->vi&&v_tmp1_found==1){
                    v_tmp2 = p_g->p_v[boundary_line[li]->vi[1]];
                    break;
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[0]->vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[0]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = p_g->p_v[boundary_line[li]->vi[0]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[0]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                    //}
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[0]->vi&&v_tmp1_found==1){
                    v_tmp2 = p_g->p_v[boundary_line[li]->vi[0]];
                    break;
                }
                else{

                }
            }
            if(v_tmp1->loc.x < v_tmp2->loc.x){
                boundary_vertices.push_back(v_tmp1);
            }
            else{
                boundary_vertices.push_back(v_tmp2);
            }
        //std::cout<<" The secondary boundary vertices is "<<boundary_vertices[1]->vi<<", and its (x,y) is ("<<boundary_vertices[1]->loc.x<<","<<boundary_vertices[1]->loc.y<<")."<<endl;

    //the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
    for(int vi=1; vi<(int)boundary_line.size()-1;vi++){
        for(int li=0;li<(int)boundary_line.size();li++){
                if(boundary_line[li]->vi[0]==boundary_vertices[vi]->vi&&boundary_line[li]->vi[1]!=boundary_vertices[vi-1]->vi){
                    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                }
                else if(boundary_line[li]->vi[1]==boundary_vertices[vi]->vi&&boundary_line[li]->vi[0]!=boundary_vertices[vi-1]->vi){
                    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                }
        }
    }
        //std::cout<<" The third boundary vertices is "<<boundary_vertices[2]->vi<<", and its (x,y) is ("<<boundary_vertices[2]->loc.x<<","<<boundary_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(boundary_vertices);

    //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(boundary_vertices[0]);
        boundary_vertices.push_back(boundary_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //the second boundary point: the point that is the left neighbor of the first point and have a distance of the boundary_averaged_distance

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
        boundary_points_distance_larger_than_surface_vertex_distance:;
        if(distance_intial_is_boundary_vertex==0) {
            surface_vertex_distance = (boundary_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
        }
        else{
            surface_vertex_distance = (boundary_vertices[vj+1]->loc-boundary_vertices[vj]->loc).norm();
        }


    //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(boundary_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = boundary_vertices[vj]->loc+(boundary_vertices[vj+1]->loc-boundary_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(boundary_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;
    }

    vector<Vertex> organ_ordered_boundary_points_finding(Organ* p_g, int boundary_points_number){
        //std::cout<<"Start boundary points finding"<<endl;
        double boundary_points_average_distance = p_g->perimeter/boundary_points_number; 
        //std::cout<<"Organ perimeter "<<p_g->perimeter<<"; boundary_points: "<<boundary_points_number<<"; the averaged distance between boundary_points: "<<boundary_points_average_distance<<endl;

        vector<Vertex> boundary_vertices;
        vector<Vertex> boundary_points;
        vector<Line> boundary_line;

            for(int li=0;li<(int)p_g->p_l.size();li++){
                if(p_g->p_l[li]->IsOutermost==1){
                    boundary_line.push_back(*p_g->p_l[li]);
                }
            }

        //we need first to arrange the surface vertices into an anticlockwisely ordered array

        //the first boundary vertex: the surface vertex with lowest y position
        Vertex first_boundary_vertex;

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                first_boundary_vertex = *p_g->p_v[vi];
                break;
            }
        }

        for(int vi=0;vi<(int)p_g->p_v.size(); vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                if(first_boundary_vertex.loc.y>p_g->p_v[vi]->loc.y){
                    first_boundary_vertex = *p_g->p_v[vi];
                }
            }
        }
        
        boundary_vertices.push_back(first_boundary_vertex);

        //std::cout<<" The first boundary vertices is "<<boundary_vertices[0]->vi<<", and its (x,y) is ("<<boundary_vertices[0]->loc.x<<","<<boundary_vertices[0]->loc.y<<")."<<endl;

    //the second boundary vertex: the surface vertex that is the neighbor of the first surface vertex (they share the same lines)
        //                              and this surface vertex should have smaller x values to make the direction anticlockwise
        Vertex v_tmp1;
        Vertex v_tmp2;
        bool v_tmp1_found=0;
            for(int li=0; li<(int)boundary_line.size(); li++){
                if(boundary_line[li].vi[0]==boundary_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[1]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[1]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[boundary_line[li].vi[1]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[1]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[1]]);
                    //}
                }
                else if(boundary_line[li].vi[0]==boundary_vertices[0].vi&&v_tmp1_found==1){
                    v_tmp2 = *p_g->p_v[boundary_line[li].vi[1]];
                    break;
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[0].vi&&v_tmp1_found==0){
                    //this outermost line connects the first boundary vertex
                    //std::cout<<"p_g->p_v[boundary_line[li]->vi[0]]->loc.x: "<<p_g->p_v[boundary_line[li]->vi[0]]->loc.x<<endl;
                    v_tmp1 = *p_g->p_v[boundary_line[li].vi[0]];
                    v_tmp1_found=1;
                    //if(p_g->p_v[boundary_line[li]->vi[0]]->loc.x<boundary_vertices[0]->loc.x){
                        //the potential secondary boundary vertex is on the anticlockwise direction of the first boundary vertex
                    //    boundary_vertices.push_back(p_g->p_v[boundary_line[li]->vi[0]]);
                    //}
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[0].vi&&v_tmp1_found==1){
                    v_tmp2 = *p_g->p_v[boundary_line[li].vi[0]];
                    break;
                }
                else{

                }
            }
            if(v_tmp1.loc.x < v_tmp2.loc.x){
                boundary_vertices.push_back(v_tmp1);
            }
            else{
                boundary_vertices.push_back(v_tmp2);
            }
        //std::cout<<" The secondary boundary vertices is "<<boundary_vertices[1].vi<<", and its (x,y) is ("<<boundary_vertices[1].loc.x<<","<<boundary_vertices[1].loc.y<<")."<<endl;

    //the tertiary and following boundary vertices could be retrieved using incremental algorithm: 
    for(int vi=1; vi<(int)boundary_line.size()-1;vi++){
        for(int li=0;li<(int)boundary_line.size();li++){
                if(boundary_line[li].vi[0]==boundary_vertices[vi].vi&&boundary_line[li].vi[1]!=boundary_vertices[vi-1].vi){
                    boundary_vertices.push_back(*p_g->p_v[boundary_line[li].vi[1]]);
                }
                else if(boundary_line[li].vi[1]==boundary_vertices[vi].vi&&boundary_line[li].vi[0]!=boundary_vertices[vi-1].vi){
                    boundary_vertices.push_back(*p_g->p_v[boundary_line[li].vi[0]]);
                }
        }
    }
        //std::cout<<" The third boundary vertices is "<<boundary_vertices[2]->vi<<", and its (x,y) is ("<<boundary_vertices[2]->loc.x<<","<<boundary_vertices[2]->loc.y<<")."<<endl;
        //cout_fout_debug::cout_vector_vertex(boundary_vertices);

    //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(boundary_vertices[0]);
        boundary_vertices.push_back(boundary_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //the second boundary point: the point that is the left neighbor of the first point and have a distance of the boundary_averaged_distance

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
        boundary_points_distance_larger_than_surface_vertex_distance:;
        if(distance_intial_is_boundary_vertex==0) {
            surface_vertex_distance = (boundary_vertices[vj+1].loc-boundary_points[pi].loc).norm();
        }
        else{
            surface_vertex_distance = (boundary_vertices[vj+1].loc-boundary_vertices[vj].loc).norm();
        }


    //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex boundary_points_new;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new.loc = boundary_points[pi].loc+(boundary_vertices[vj+1].loc-boundary_points[pi].loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new.loc = boundary_vertices[vj].loc+(boundary_vertices[vj+1].loc-boundary_vertices[vj].loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(boundary_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi].loc.x<<","<<boundary_points[pi].loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;
    }


    //judge if a point is inside a polygon (eg., p_g) by ray casting algorithm; if inside bool=1; else outside bool=0
    //reference: 1. https://en.wikipedia.org/wiki/Point_in_polygon#
    //reference: 2. https://www.youtube.com/watch?v=RSXM9bgqxJM
    bool point_in_polygon_ray_casting(Vertex* point_tmp, Organ* p_g){
        bool inside_bool=0;
        int ray_casting_index=0;
        for(int li=0; li<(int)p_g->p_l.size();li++){
            if(p_g->p_l[li]->IsOutermost==1){
                double y0_tmp = p_g->p_v[p_g->p_l[li]->vi[0]]->loc.y;
                double x0_tmp = p_g->p_v[p_g->p_l[li]->vi[0]]->loc.x;
                double y1_tmp = p_g->p_v[p_g->p_l[li]->vi[1]]->loc.y;
                double x1_tmp = p_g->p_v[p_g->p_l[li]->vi[1]]->loc.x;

                double y_point = point_tmp->loc.y;
                double x_point = point_tmp->loc.x;
                if((y0_tmp>y_point)!=(y1_tmp>y_point)){
                    if(x_point<(x0_tmp+(y_point-y0_tmp)/(y1_tmp-y0_tmp)*(x1_tmp-x0_tmp))){
                        ray_casting_index++;
                        //std::cout<<"Ray casting passed line "<<li<<"; p0 "<<x0_tmp<<","<<y0_tmp<<"; p1 "<<x1_tmp<<","<<y1_tmp<<endl;
                    }
                }
            }
        }

        if(ray_casting_index%2==0){
            //even
            inside_bool = 0;
        }
        if(ray_casting_index%2==1){
            //odd
            inside_bool = 1;
        }
        //std::cout<<"ray_casting_index "<<ray_casting_index<<endl;
        return inside_bool;
    }
   
    
   
    //do calculate the circle fitting of all boundary points by kasa
    vector<double> curvature_circle_fitting_kasa_three_boundary_points(Organ* p_g, vector<Vertex*> boundary_points, int points_away, int NumAveraging){

        vector<double> curvature;
        vector<double> curvature_averaged;

        int boundary_points_number = (int)boundary_points.size();

    //1. start circle fitting (Kasa)
    for(int vi=0;vi<(int)boundary_points.size();vi++){
        Vertex* p_minus = new Vertex;
        Vertex* p_0 = new Vertex;
        Vertex* p_plus = new Vertex;

        p_0 = boundary_points[vi];
        if(vi<points_away){
            p_minus=boundary_points[vi-points_away+boundary_points_number];
        }
        else{
            p_minus=boundary_points[vi-points_away];
        }

        if(vi>(boundary_points_number-points_away-1)){
            p_plus=boundary_points[vi+points_away-boundary_points_number];
        }
        else{
            p_plus=boundary_points[vi+points_away];
        }
        vector<Vertex*> fitting_points;
        fitting_points.push_back(p_minus);
        fitting_points.push_back(p_0);
        fitting_points.push_back(p_plus);
        Vertex* midpoint_tmp = new Vertex;
        midpoint_tmp->loc.x = (p_minus->loc.x + p_plus->loc.x)/2.0;
        midpoint_tmp->loc.y = (p_minus->loc.y + p_plus->loc.y)/2.0;
        double curvature_abstract_value = boundary_geo::surface_vertex_curvature_kasa(fitting_points);
        double curvature_value=0;
        
        if(organ_geo::point_in_polygon_ray_casting(midpoint_tmp,p_g)==0){
            curvature_value= -curvature_abstract_value;
        }
        else{
            curvature_value = curvature_abstract_value;
        }
    
        curvature.push_back(curvature_value);
    }

    //2. taking average of curvature value

        if(NumAveraging==1){
            return curvature;
        }
        else{
            curvature_averaged = geo_vv::vd_averaged(curvature,NumAveraging);
            return curvature_averaged;
        }  
    }

    //calculating curvature by circle fitting on boundary points with equal distance along polygon space: based on previous four functions
    vector<double> curvature_circle_fitting_kasa_along_polygon(Organ* p_g,int boundary_points_number, int points_away, int NumAveraging){
        vector<double> curvature_averaged;
        //std::cout<<"Start Curvature calculation based on Circle_fitting"<<endl;
        //1. order surface vertices into anticlockwise order
        vector<Vertex> anticlockwise_surface_vertices = organ_geo::organ_ordered_anticlockwise_boundary(p_g).vi;
        anticlockwise_surface_vertices = geo_vv::normalization(anticlockwise_surface_vertices);
        //std::cout<<"Ordered_bondary_generated"<<endl;
        //2. find boundary points with equal distance along the closed contour of polygon
        vector<Vertex> boundary_points = organ_geo::organ_ordered_boundary_points_finding(p_g, boundary_points_number);
        boundary_points = geo_vv::normalization(boundary_points);
        //cout_fout_debug::fout_vector_vertex(boundary_points,"boundary_points.txt");
        //std::cout<<"Boundary_points_generated"<<endl;
        //3. do the circle fitting for each points: the point and the two points that are 10 points away
        //curvature_averaged = organ_geo::curvature_circle_fitting_kasa_three_boundary_points(p_g, boundary_points, points_away);
        curvature_averaged = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(boundary_points,points_away,NumAveraging);

        p_g->minimum_curvature=geo_vv::vd_minimum(curvature_averaged);
        p_g->maximum_curvature=geo_vv::vd_maximum(curvature_averaged);
        p_g->accumulated_negative_curvature=geo_vv::accumulated_negative(curvature_averaged);

        
        /*
        for(Vertex* anti_surface_vertex : anticlockwise_surface_vertices){
            if(anti_surface_vertex != nullptr){
                delete anti_surface_vertex;
                anti_surface_vertex = nullptr;
            }
            
        }
        for(Vertex* bound_points : boundary_points){
            if(bound_points != nullptr){
                delete bound_points;
                bound_points = nullptr;
            }
        }
        boundary_points.clear();
        anticlockwise_surface_vertices.clear();
        */
        return curvature_averaged;
    }

    double from_curvature_to_tip_position(Organ* p_g,int boundary_points_number, int points_away, int NumAveraging){
        //std::cout<<"Start Curvature calculation based on Circle_fitting"<<endl;
        //1. order surface vertices into anticlockwise order
        vector<Vertex> anticlockwise_surface_vertices = organ_geo::organ_ordered_anticlockwise_boundary(p_g).vi;
        anticlockwise_surface_vertices = geo_vv::normalization(anticlockwise_surface_vertices);
        //std::cout<<"Ordered_bondary_generated"<<endl;
        //2. find boundary points with equal distance along the closed contour of polygon
        vector<Vertex> boundary_points = organ_geo::organ_ordered_boundary_points_finding(p_g, boundary_points_number);
        boundary_points = geo_vv::normalization(boundary_points);
        //cout_fout_debug::fout_vector_vertex(boundary_points,"boundary_points.txt");
        //std::cout<<"Boundary_points_generated"<<endl;
        //3. do the circle fitting for each points: the point and the two points that are 10 points away
        //curvature_averaged = organ_geo::curvature_circle_fitting_kasa_three_boundary_points(p_g, boundary_points, points_away);
        
        vector<double> curvature_averaged = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(boundary_points,points_away,NumAveraging);
        double minimum_curvature=curvature_averaged[0];
        int minimum_curvature_position;

        for(int i=0;i<curvature_averaged.size();i++){
            if(minimum_curvature>curvature_averaged[i]){
                minimum_curvature=curvature_averaged[i];
                minimum_curvature_position=i;
            }
        }

        std::cout<<"Minimum curvature: "<<minimum_curvature<<std::endl;
        std::cout<<"Minimum curvature position (boundary point index): "<<minimum_curvature_position<<std::endl;
        //from boundary point index to relative y position
        double min_c_position_relative = boundary_points[minimum_curvature_position].loc.y;
        std::cout<<"min_c position (relative y): "<<min_c_position_relative<<std::endl;
        //from relative y position to absolute y position
        double min_c_position_absolute = min_c_position_relative * (p_g->y_max_vertex-p_g->y_min_vertex); 
        std::cout<<"min_c position (absolute y): "<<min_c_position_absolute<<std::endl;
        double tip_length_min_c_position_absolute = (1-min_c_position_relative) * (p_g->y_max_vertex-p_g->y_min_vertex); 
        std::cout<<"tip_length_min_c_position_absolute "<<tip_length_min_c_position_absolute<<std::endl;

        p_g->tip_length_min_c_position_absolute = tip_length_min_c_position_absolute;
        p_g->tip_length_min_c_position_relative = 1-min_c_position_relative;

        return p_g->tip_length_min_c_position_relative;
    }

    vector<Vertex> organ_boundary_points_euclidean(Organ* p_g, double boundary_distance){
        vector<Vertex> boundary_points;

        //0. get the boundary vertex and boundary lines sorted in anticlockwise direction;
        Ordered_boundary Ordered_boundary_tmp = organ_geo::organ_ordered_anticlockwise_boundary(p_g);
        vector<Vertex> boundary_vertex = Ordered_boundary_tmp.vi;
        vector<int> boundary_line = Ordered_boundary_tmp.li;

        //1. find the first boundary points: the surface vertices with lowest y position
        //std::cout<<"********************* Finding Point 0 ***********************"<<endl;
        boundary_points.push_back(boundary_vertex[0]);
        //std::cout<<"The 0 boundary point "<<boundary_points[0]->loc.x<<","<<boundary_points[0]->loc.y<<endl;
        //std::cout<<"********************* Point 0 Found ***********************"<<endl;
        //2. find the second boundary points:
         //2.1 draw the first circle with center on the first point
         //std::cout<<"********************* Finding Point 1 ***********************"<<endl;
         Circle cir_0;
         cir_0.center = boundary_points[0].loc;
         cir_0.radius = boundary_distance;
         //2.2 find the cross point between cir_0 and the organ boundary
         int ci=0;  //now we are drawing the cir_0
         int vj=0; //now we are finding intersection between line 0 and cir_0
         
          //2.2.1 search if there is any cross point between the cir_0 and line {boundary_points[0],boundary_vertex[1]}
          Line ls_0;
          ls_0.set_endpoints(boundary_points[0],boundary_vertex[1]);
          Intersection_relationship ls_0_cir_0_intersections = circle_geo::intersection_line_segment_circle(ls_0,cir_0);
          //std::cout<<"Intersection relationship between cir_0 and ls_0 is "<<ls_0_cir_0_intersections.Relationship<<endl;
          if(ls_0_cir_0_intersections.Relationship>0)
           //std::cout<<"The intersection is "<<ls_0_cir_0_intersections.cross_points[0].loc.x<<","<<ls_0_cir_0_intersections.cross_points[0].loc.y<<endl;
          if(ls_0_cir_0_intersections.Relationship==0)
          {
            //2.2.2 if there is no crosspoint between cir_0 and ls_0 (boundary_points[0], boundary_vertex[1])
            //try to find the crosspoint between the cir_0 and ls_1,ls_2,...,ls_j, until the crosspoint is found
            for( ; vj<(int)boundary_vertex.size(); vj++){
                Line ls_j;
                ls_j.set_endpoints(boundary_vertex[vj],boundary_vertex[vj+1]);
                Intersection_relationship ls_j_cir_0_intersections = circle_geo::intersection_line_segment_circle(ls_j,cir_0);
                if(ls_j_cir_0_intersections.Relationship==0){

                }
                else{
                    boundary_points.push_back(ls_j_cir_0_intersections.cross_points[0]);
                    ci++;
                    break;
                }
                if(vj==((int)boundary_points.size()-1))
                {
                    std::cout<<"The distance between boundary points is too large, there is no qualified second boundary point"<<endl;
                    exit(1);
                }
            }
          }
          else{
            //2.2.3 there is crosspoint between cir_0 and ls_0, this should be the second boundary point
            boundary_points.push_back(ls_0_cir_0_intersections.cross_points[0]);
            ci++;
          }
        //std::cout<<"The 1 boundary point "<<boundary_points[1]->loc.x<<","<<boundary_points[1]->loc.y<<endl;
        //std::cout<<"********************* Point 1 Found ***********************"<<endl;
        //std::cout<<"********************* Finding the Point 2 "<<" ***********************"<<endl;

        //3. find the third, fourth, ... i-th, boundary points
         for(; vj<(int)boundary_vertex.size()-2;)
         {

            //3.1 draw a circle on the current boundary points (i-th)
            Circle cir_i;
            cir_i.center = boundary_points[ci].loc;
            //std::cout<<"cir_"<<ci<<" center "<<cir_i.center.x<<","<<cir_i.center.y<<endl;
            cir_i.radius = boundary_distance;
            //3.2 find the intersection between the cir_i and the remaining boundary surface (excluding the searched region)
            //3.2.1 find the intersection between the cir_i and boundary_line ls_j (the i-th boundary_point, and the vj+1 boundary vertex);
            Line ls_j;
            //std::cout<<"vj"<<vj<<endl;
            ls_j.set_endpoints(boundary_points[ci],boundary_vertex[vj+1]);
            Intersection_relationship ls_j_cir_i_intersection = circle_geo::intersection_line_segment_circle(ls_j,cir_i);
            //std::cout<<"0 Intersection relationship between cir"<<ci<<" and ls"<<vj<<" is "<<ls_j_cir_i_intersection.Relationship<<endl;
            if(ls_j_cir_i_intersection.Relationship>0){
                //std::cout<<"The intersection is "<<ls_j_cir_i_intersection.cross_points[0].loc.x<<","<<ls_j_cir_i_intersection.cross_points[0].loc.y<<endl;
            }
            //3.2.2 if there is no intersection, try to find intersection between cir_i and ls_j+1, ls_j+2,..., until the intersection is found
            if(ls_j_cir_i_intersection.Relationship==0){
                
                for( ; vj<(int)boundary_vertex.size()-2; ){
                    //3.2.2.1 ls_j+k is made by v_j+k boundary vertex and v_j+k+1 boundary vertex
                    vj++;
                    Line ls_jk;
                    ls_jk.set_endpoints(boundary_vertex[vj], boundary_vertex[vj+1]);
                    Intersection_relationship ls_jk_cir_i_intersection = circle_geo::intersection_line_segment_circle(ls_jk,cir_i);
                    //std::cout<<"1 Intersection relationship between cir"<<ci<<" and ls"<<vj<<" is "<<ls_jk_cir_i_intersection.Relationship<<endl;
                    if(ls_jk_cir_i_intersection.Relationship==0){
                        //std::cout<<"No intersection found"<<endl;
                    }
                    else{
                        Vertex current_boundary_point;
                        current_boundary_point.loc = ls_jk_cir_i_intersection.cross_points[0].loc;
                        boundary_points.push_back(current_boundary_point);
                        //cout_fout_debug::cout_vector_vertex(boundary_points);                        
                        //std::cout<<"The "<<ci+1<<"-th boundary point "<<boundary_points[ci+1]->loc.x<<","<<boundary_points[ci+1]->loc.y<<endl;
                        //std::cout<<"********************* Point "<<ci+1<< " Found ***********************"<<endl;
                        ci++;
                        //std::cout<<"********************* Finding the Point "<<ci+1<<" ***********************"<<endl;
                        break;
                    }
                }
            }
            else{
            //3.2.3 there is intersection between cir_i and ls_j, this should be the i-th boundary points
                Vertex current_boundary_point;
                current_boundary_point.loc = ls_j_cir_i_intersection.cross_points[0].loc;
                boundary_points.push_back(current_boundary_point);
                //cout_fout_debug::cout_vector_vertex(boundary_points);
                //std::cout<<"The "<<ci+1<<"-th boundary point "<<boundary_points[ci+1]->loc.x<<","<<boundary_points[ci+1]->loc.y<<endl;
                //std::cout<<"********************* Point "<<ci+1<< " Found ***********************"<<endl;
                ci++;
                //std::cout<<"********************* Finding the Point "<<ci+1<<" ***********************"<<endl;

            }

            
         }
            //4. for the last intersection: it should between vj_end and vj_0
            Circle cir_last;
            Line ls_last;
            Intersection_relationship last_Intersection_result;
            last_boundary_points_euclidean:;
            cir_last.center = boundary_points[boundary_points.size()-1].loc;
            cir_last.radius = boundary_distance;
            ls_last.set_endpoints(boundary_points[boundary_points.size()-1],boundary_points[0]);
            last_Intersection_result = circle_geo::intersection_line_segment_circle(ls_last,cir_last);
            if(last_Intersection_result.Relationship==0)
            {

            }
            else
            {
                Vertex last_boundary_point;
                last_boundary_point.loc = last_Intersection_result.cross_points[0].loc;
                boundary_points.push_back(last_boundary_point);
                //std::cout<<"The "<<ci+1<<"-th boundary point "<<boundary_points[ci+1]->loc.x<<","<<boundary_points[ci+1]->loc.y<<endl;
                //std::cout<<"********************* Point "<<ci+1<< " Found ***********************"<<endl;
                ci++;
                //std::cout<<"********************* Finding the Point "<<ci+1<<" ***********************"<<endl;
                goto last_boundary_points_euclidean;
            }
             //std::cout<<"********************* Point "<<ci<<" is the Last Point ***********************"<<endl;

            //std::cout<<"boundary_lines size "<<boundary_line.size()<<endl;
            //5. check the distance between the last boundary and the first boundary point
            double dist_last_first = vertex_geo::vertex_distance(boundary_points[0],boundary_points[boundary_points.size()-1]);
            //store the distance in 
            boundary_points[0].loc.z = dist_last_first;
            //std::cout<<"Distance between the last boundary point and the first boudary point is: "<<dist_last_first<<endl;
        
            //6. output the result to "boundary_points.txt"
            //cout_fout_debug::fout_vector_vertex(boundary_points,"boundary_points.txt");

          

        return boundary_points;

    }

    vector<tuple<double,double,double>> organ_boundary_points_euclidean_point_number(Organ* p_g, int boundary_point_number){
        vector<tuple<double,double,double>> point_number;;
        double initial_distance_guess=0.1;
        double distance_step = 0.0001;
        double terminal_distance_guess=10.0;

        vector<double> boundary_point_number_vec;

        vector<Vertex*> boundary_point_number_vv;

        for(double i=initial_distance_guess; i < terminal_distance_guess; i=i+distance_step){
            vector<Vertex> boundary_points = organ_geo::organ_boundary_points_euclidean(p_g,i);
            double boundary_points_number = boundary_points.size();
            point_number.push_back(make_tuple(i,boundary_points_number,boundary_points[0].loc.z));
            //std::cout<<"Detecting distance: "<<i<<" , the boundary point number "<<boundary_points_number<<endl;
        }      

        cout_fout_debug::fout_tuple_double(point_number,"point_number.txt");
        

        return point_number;
    }
 
    vector<double> curvature_graph_of_a_function(Organ* p_g,vector<Vertex*> boundary_points)
    {   
        //graph of a function: y=f(x)
        //k = abs(y'')/(1+y'^2)^(3/2)
        //0. organ boundary points ordered in anticlockwise direction and has the same euclidean distance 
        int n=boundary_points.size();
        vector<double> curvature_graph_function;

        //1. calculate the first derivative of f(x_i)
        vector<double> first_derivative_fx = geo_vv::first_derivative(boundary_points);

        //2. calculate the second derivative of f(x_i)
        vector<double> second_derivative_fx = geo_vv::second_derivative(boundary_points);        

        //3. calculate the curvature based on the formula: k = abs(y'')/(1+y'^2)^(3/2)
        for(int i=0;i<n;i++){
            double curvature_tmp = abs(second_derivative_fx[i])/pow((1+first_derivative_fx[i]*first_derivative_fx[i]),1.50);
            curvature_graph_function.push_back(curvature_tmp);
        }
        return curvature_graph_function;
    }

    vector<double> curvature_arc_length_parameterization(Organ* p_g, vector<Vertex*> boundary_points, double boundary_distance){
        //0. organ boundary points ordered in anticlockwise direction and had the same euclidean distance
        int n=boundary_points.size();
        vector<double> curvature_arc_length_parameters;
        boundary_points[0]->loc.z=0;

        //1. s_i for point i
        for(int i=1;i<n;i++){
            boundary_points[i]->loc.z = boundary_points[i-1]->loc.z+boundary_distance;
        }

        //2. calculate the first derivative of x_i(s) and y_i(s)
        vector<Vertex*> vv_xs, vv_ys;
         for(int i=0;i<n;i++){
            Vertex* vi = new Vertex;
            vi->loc.x = boundary_points[i]->loc.z;
            vi->loc.y = boundary_points[i]->loc.x;
            vv_xs.push_back(vi);
            Vertex* vj = new Vertex;
            vj->loc.x = boundary_points[i]->loc.z;
            vj->loc.y = boundary_points[i]->loc.y;
            vv_ys.push_back(vj);
        }
        vector<double> first_derivative_xs = geo_vv::first_derivative(vv_xs);
        vector<double> first_derivative_ys = geo_vv::first_derivative(vv_ys);

        //3. calculate the second derivative of x_i(s) and y_i(s)
        vector<double> second_derivative_xs = geo_vv::second_derivative_fourth_order(vv_xs);
        vector<double> second_derivative_ys = geo_vv::second_derivative_fourth_order(vv_ys);
        cout_fout_debug::fout_vector_double(first_derivative_xs,"xs_first.txt");
        cout_fout_debug::fout_vector_double(first_derivative_ys,"ys_first.txt");
        cout_fout_debug::fout_vector_double(second_derivative_xs,"xs_second.txt");
        cout_fout_debug::fout_vector_double(second_derivative_ys,"ys_second.txt");
        //4. calculate the curvature based on the second derivative of x_i(s) and y_i(s)
        //k(s) = sqrt(x''(s)^2 + y''(s)^2)
        for(int i=0;i<n;i++){
            double curvature_tmp = sqrt(second_derivative_xs[i]*second_derivative_xs[i]+second_derivative_ys[i]*second_derivative_ys[i]);
            curvature_arc_length_parameters.push_back(curvature_tmp);
        }

        //5. calculating the symbol of curvature: based on the formula
        // Center of curvature = gamma(s) + 1/k(s)^2*T'(s)
        // gamma(s) = (x(s),y(s))
        // k(s) = \sqrt(x''(s),y''(s))
        // T(s) = (x'(s),y'(s))
        // T'(s) = (x''(s),y''(s))
        for(int i=0;i<n;i++){
            Vertex curvature_center_i;
            curvature_center_i.loc = boundary_points[i]->loc + 1/(second_derivative_xs[i]*second_derivative_xs[i]+second_derivative_ys[i]*second_derivative_ys[i])*_vec<double>{second_derivative_xs[i],second_derivative_ys[i],0.0};
            //whether vertex curvature center_i is inside the organ or outside of the organ
            bool symbol_i = point_in_polygon_ray_casting(&curvature_center_i,p_g);
            if(symbol_i==0){
                //negative 
                curvature_arc_length_parameters[i] = 0.0-curvature_arc_length_parameters[i];
            }
            else{
                //positive
            }
        }

        return curvature_arc_length_parameters;
    }

    double complexity_index_information_theory(Organ* p_g, vector<Vertex*> boundary_points){
            double complexity_index=0;
            //find the optimal resolution of J
            //calculate the quantization error
            int initial_J=1;
            int step_J=1;
            int terminal_J=25;

            double error_max=0;
            vector<pair<double,double>> entropy_error_vec;
            for(int Ji=0; Ji<(int)terminal_J; Ji=Ji+step_J){
                pair<double,double> entropy_error = organ_geo::distance_entropy_error_information_theory(p_g,Ji, boundary_points);
                if(error_max<entropy_error.second){
                    error_max=entropy_error.second;
                }
                entropy_error_vec.push_back(entropy_error);
            }

            double cost_function_min=0;
            int cost_function_min_J=0;
            for(int Ji=0; Ji<(int)terminal_J; Ji=Ji+step_J)
            {
                double cost_function_J = entropy_error_vec[Ji].first/log2(p_g->N_epi_cell) + entropy_error_vec[Ji].second/error_max;
                if(cost_function_J>cost_function_min){
                    cost_function_min_J = Ji;
                    cost_function_min = cost_function_J;
                }
            }

            std::cout<<"cost_function_min_J "<<cost_function_min_J<<endl;
            std::cout<<"The distance entropy is "<<entropy_error_vec[cost_function_min_J].first<<" and the quantization error is "<<entropy_error_vec[cost_function_min_J].second<<endl;
        }

    pair<double,double> distance_entropy_error_information_theory(Organ* p_g, int J, vector<Vertex*> boundary_points)
        {   
            //0. preparing the boundary points (400)

            //1. calculate the normalized global distance
            vector<double> global_distance;
            double global_distance_max=0;
            vector<double> global_distance_normalized;
            for(int vi=0;vi<(int)boundary_points.size();vi++){
                double distance_vi = (boundary_points[vi]->loc - p_g->center).norm();
                global_distance.push_back(distance_vi);
                if(global_distance[vi]>global_distance_max){
                    global_distance_max=global_distance[vi];
                }
            }

            for(int vi=0;vi<(int)boundary_points.size();vi++){
                global_distance_normalized.push_back(global_distance[vi]/global_distance_max);
            }
            //cout_fout_debug::cout_vector_double(global_distance_normalized);
            
            //2. create the distance histogram, calcualte quantization values (the averaged distance in a single bin)
            double bins = pow(2,J);
            vector<int> distance_histogram(bins);
            vector<double> quantization_value(bins);
            for(int vi=0;vi<(int)global_distance_normalized.size();vi++){
                for(int ji=0;ji<bins;ji++){
                    if(global_distance_normalized[vi]>(ji/bins)&&global_distance_normalized[vi]<((ji+1)/bins)){
                        distance_histogram[ji]++;
                        quantization_value[ji] += global_distance_normalized[vi];
                    }
                }
            }
            for(int ji=0; ji<bins; ji++){
                if(distance_histogram[ji]!=0){
                    quantization_value[ji]= quantization_value[ji]/distance_histogram[ji];
                }
            }
            //cout_fout_debug::cout_vector_int(distance_histogram);
            
            //3. calculate the distance entropy
            double distance_entropy=0;
            //to get PDF from histogram bins, the bins should be normalized by the total number of bins and bin width
            //vector<double> distance_PDF = distance_histogram/(bins*bins_width) = distance_histogram;
            vector<double> distance_PDF(bins);
            for(int ji=0; ji<bins; ji++){
                distance_PDF[ji] = distance_histogram[ji]/(double)global_distance_normalized.size();
            } 
            for(int ji=0; ji<bins; ji++){
                // H = - \Sigma pj(k)*log2(pj(k))
                if(distance_PDF[ji]!=0){
                    distance_entropy += -distance_PDF[ji]*log2(distance_PDF[ji]);  
                    //std::cout<<"Distance entropy: "<<distance_entropy<<endl;
                }
            }
            
            //4. calculate quantization error
            double quantization_error=0;

            for(int vi=0; vi<(int)global_distance_normalized.size();vi++){
                double quantization_vi;
                for(int ji=0;ji<bins;ji++){
                    if(global_distance_normalized[vi]>(ji/bins)&&global_distance_normalized[vi]<((ji+1)/bins)){
                        quantization_vi = quantization_value[ji];
                    }
                }
                quantization_error += pow((global_distance_normalized[vi] - quantization_vi),2)/global_distance_normalized.size();
            }

            quantization_error = sqrt(quantization_error);
            std::cout<<"J "<<J<<endl;
            std::cout<<"Distance entropy: "<<distance_entropy<<endl;
            std::cout<<"Quantization_error "<<quantization_error<<endl;
            return make_pair(distance_entropy,quantization_error);
        }

    //minor analysis
    vector<tuple<int,double,double>> n_gons_analysis(Organ* p_g){
        vector<tuple<int,double,double>> n_gons_analysis_results;
        int minimum_lines_of_polygon=3;
        int maxmum_lines_of_polygon=11;
        
        for(int i=minimum_lines_of_polygon;i<maxmum_lines_of_polygon;i++){
            tuple<int,double,double> pi = {i,0,0}; 
            n_gons_analysis_results.push_back(pi);
        }
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            int ni=p_g->p_c[ci]->li.size();
            for(auto&t : n_gons_analysis_results){
                if(std::get<0>(t) == ni){
                    std::get<1>(t) ++;
                    break;
                }
            }
        }

        

        for(auto&t : n_gons_analysis_results){
            std::get<2>(t) = std::get<1>(t)/(double)p_g->p_c.size();
        }

        for(const auto& t: n_gons_analysis_results){
            std::cout << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << std::endl;
        }
        
        std::ofstream file("n_gons.txt",std::ios::app);
        if (!file) {
            std::cerr << "Failed to open the file for appending." << std::endl;
            exit(1);
        }
        file <<p_g->step<<" ";
        for(const auto& t: n_gons_analysis_results){
            file << std::get<1>(t) << "  " << std::get<2>(t)<<" ";
        }
        file<<std::endl;
        file.close();

        return n_gons_analysis_results;
    }

}

namespace geo{

void calcGeometrics(Organ* p_g){ 
    
    organ_geo::organ_line_length(p_g);
    organ_geo::organ_cell_perimeter(p_g);


    //1. epidermal identity
    organ_geo::epidermal_identity(p_g);
    //2. organ center
    organ_geo::organ_center(p_g);
    organ_geo::organ_max_min_x_y(p_g);
    //3. organ cell layer
    organ_geo::organ_cell_layer(p_g);
    //4. organ area
    organ_geo::organ_area(p_g);
    //5. organ perimeter
    organ_geo::organ_perimeter(p_g);
    //6. organ circularity
    organ_geo::organ_circularity(p_g);
    //7. organ regularity
    organ_geo::organ_regularity(p_g);
    //8. organ length, width, length/width ratio
    organ_geo::organ_length(p_g);
    organ_geo::organ_width(p_g);
    //9. cell arrangement
    organ_geo::organ_cell_arrangement(p_g);
    //10. organ elliptical index
    organ_geo::elliptical_index(p_g);
    //11. organ overlap index
    //organ_geo::organ_overlap(p_g);
    //12. similarity index
    if(similarity_calculation_required==1){
        boundary_geo::similarity_cal_during_simulation(p_g,real_organ_contour_processed_for_similarity_index);
    }
    //13. curvature analysis
    organ_geo::curvature_circle_fitting_kasa_along_polygon(p_g,100,10,9); 
    
    //for mechanics mode
    if(mechanics_mode=="S_std_spatial"){
        force::S_std_spatial_geo(p_g);
    }
}

void calcGeometrics_cout(Organ* p_g){ 
    std::cout<<"**********Start to calcGeometrics_cout**********"<<std::endl;
    organ_geo::organ_line_length(p_g);
    organ_geo::organ_cell_perimeter(p_g);
    organ_geo::organ_max_min_x_y(p_g);

    //1. epidermal identity
    organ_geo::epidermal_identity(p_g);
    std::cout<<"Number of epidermal cells: "<<p_g->N_epi_cell<<"; number of inner cells: "<<p_g->N_inner_cell<<std::endl;
    //2. organ center
    organ_geo::organ_center(p_g);
    std::cout<<"Organ center: ("<<p_g->center.x<<","<<p_g->center.y<<")"<<std::endl;
    //3. organ cell layer
    organ_geo::organ_cell_layer(p_g);
    std::cout<<"Organ has "<<p_g->cell_layer_number<<" layers of cells"<<std::endl;
    //4. organ area
    organ_geo::organ_area(p_g);
    std::cout<<"organ_area: "<<p_g->area<<std::endl;
    //5. organ perimeter
    organ_geo::organ_perimeter(p_g);
    std::cout<<"organ_perimeter: "<<p_g->perimeter<<std::endl;
    //6. organ circularity
    organ_geo::organ_circularity(p_g);
    std::cout<<"organ_circularity: "<<p_g->circularity<<std::endl;
    //7. organ regularity
    organ_geo::organ_regularity(p_g);
    std::cout<<"organ_regularity: "<<p_g->regularity_averaged<<std::endl;
    //8. organ length, width, length/width ratio
    organ_geo::organ_length(p_g);
    std::cout<<"organ_length: "<<p_g->organ_length<<std::endl;
    organ_geo::organ_width(p_g);
    std::cout<<"organ_width: "<<p_g->organ_width<<std::endl;
    //9. cell arrangement
    organ_geo::organ_cell_arrangement(p_g);
    std::cout<<"The x-axis cell layers is: "<<p_g->cell_arrangement_x<<", the y-axis cell layers is: "<<p_g->cell_arrangement_y<<", the cell arrangement ratio is: "<<p_g->cell_arrangement_ratio<<std::endl;
    //10. organ elliptical index
    organ_geo::elliptical_index(p_g);
    std::cout<<"organ_elliptical_index: "<<p_g->elliptical_index<<std::endl;
    //11. organ_similarity_index
    if(similarity_calculation_required==1){
        boundary_geo::similarity_cal_during_simulation(p_g,real_organ_contour_processed_for_similarity_index);
    }
    std::cout<<"organ_similarity_index: "<<p_g->similarity_index<<std::endl;
    //13. curvature analysis
    organ_geo::curvature_circle_fitting_kasa_along_polygon(p_g,100,10,3); 
    std::cout<<"organ_min_curvature: "<<p_g->minimum_curvature<<std::endl;
    //
    std::cout<<"**********Geometrics output finished**********"<<std::endl;
}

void basic_VTK_analysis(Organ *p_g){

    //organ_geo::organ_vertex_counterclockwise_sort(p_g);
    organ_geo::organ_line_length(p_g);
    organ_geo::organ_cell_perimeter(p_g);
    organ_geo::organ_center(p_g);
    organ_geo::organ_area(p_g);
    organ_geo::epidermal_identity(p_g);
    organ_geo::organ_perimeter(p_g);
    organ_geo::organ_regularity(p_g);
    organ_geo::organ_circularity(p_g);
    organ_geo::organ_cell_layer(p_g);
    //organ_geo::organ_overlap(p_g);
    organ_geo::organ_max_min_x_y(p_g);
    boundary_geo::similarity_cal_during_simulation(p_g,real_organ_contour_processed_for_similarity_index);
}

double distP(double x1, double y1, double x2, double y2){
        return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    }

bool relaP(double x1, double y1, double x2, double y2){
        double sigma = 1e-7;
        if(abs(x1-x2)<sigma&&abs(y1-y2)<sigma){
            return 1;
        }
        else{
            return 0;
        }
    }
    
//relationship between lines, =2, identical, =1, parallel, =0, intersection
int relaL(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2){
        //y=ax+b => a = (y1-y2)/(x1-x2), b = (x1y2-y1x2)/(x1-x2)
        double a1 = (yi1-yi2)/(xi1-xi2);
        double a2 = (yj1-yj2)/(xj1-xj2);
        double b1 = (xi1*yi2-yi1*xi2)/(xi1-xi2);
        double b2 = (xj1*yj2-yj1*xj2)/(xj1-xj2);

        double sigma = 1e-7;
        if(abs(a1-a2)<sigma){
            if(abs(b1-b2)<sigma){
                return 2;
            }
            else{
                return 1;
            }   
        }
        else{
            return 0;
        }

    }
    
//find intersection between lines
pair<double, double> intersectLine(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2){
        //confirm that these two lines intersect
        if(relaL(xi1,yi1, xi2, yi2, xj1, yj1, xj2, yj2)!=0){
            std::cout<<"These two lines do not intersect !"<<endl;
        }

        //y1=ax1+b, y2=ax2+b => a = (y1-y2)/(x1-x2), b = (x1y2-y1x2)/(x1-x2)
        double a1 = (yi1-yi2)/(xi1-xi2);
        double a2 = (yj1-yj2)/(xj1-xj2);
        double b1 = (xi1*yi2-yi1*xi2)/(xi1-xi2);
        double b2 = (xj1*yj2-yj1*xj2)/(xj1-xj2);
        //y=a1x+b1, y=a2x+b2 => x = (b2-b1)/(a1-a2), y = (a1b2-a2b1)/(a1-a2)
        double intersect_x = (b2-b1)/(a1-a2);
        double intersect_y = (a1*b2-a2*b1)/(a1-a2);

        pair<double, double> intersect;
        intersect.first = intersect_x;
        intersect.second = intersect_y;
        return intersect;
    }

//=0, if not intersected, =1 if intersected
int relaRS_1(double xi1, double yi1, double xj1, double yj1, double xj2, double yj2){
        double max_yj = max(yj1,yj2);
        double min_yj = min(yj1,yj2);
        if(yi1>max_yj||yi1<min_yj){
            return 0;
        }
        
        //find intersection between y=yi1(x>xi1) and y=ax+b
        double a = (yj1-yj2)/(xj1-xj2);
        double b = (xj1*yj2-yj1*xj2)/(xj1-xj2);
        //do not consider coincidence
        double intersect_x = (yi1-b)/a;
        double intersect_y = yi1;
        if(intersect_x>xi1){
            return 1;
        }
        else{
            return 0;
        }
    }

void batch(Batch *ba, Organ *p_g){
    //collect N_in, N_epi
    int inner_cell_number=0, peripheral_cell_number=0;

    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
        if(p_g->p_c[ci]->IsEpidermal==0){
            inner_cell_number++;
        }
        else{
            peripheral_cell_number++;
        }
    }
    ba->N_epi.push_back(peripheral_cell_number);
    ba->N_in.push_back(inner_cell_number);
    //collect organ perimeter and organ area
    double organ_area, organ_perimeter=0;
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
        organ_area += p_g->p_c[ci]->area;
        if(p_g->p_c[ci]->IsEpidermal==1){
            organ_perimeter += p_g->p_c[ci]->outermostLength;
        }
    
    }
    ba->organ_area.push_back(organ_area);
    ba->organ_perimeter.push_back(organ_perimeter);
    //collect averaged inner and epidermal cell area
    double in_organ_area, epi_organ_area=0;
    for(int ci=0;ci<(int)p_g->p_c.size();ci++){
        if(p_g->p_c[ci]->IsEpidermal==1){
            epi_organ_area += p_g->p_c[ci]->area;
        }
        else{
            in_organ_area += p_g->p_c[ci]->area;
        }
    }
    ba->averaged_epi_area.push_back(epi_organ_area/peripheral_cell_number);
    ba->averaged_inner_area.push_back(in_organ_area/inner_cell_number);
    //collect averaged circumference 
    ba->averaged_perimeter.push_back(organ_perimeter/peripheral_cell_number);
    
    //collect averaged cell regularity
    ba->averaged_regularity_in.push_back(p_g->regularity_in_av);
    ba->averaged_regularity_epi.push_back(p_g->regularity_epi_av);
    //collect overlaps informtion
    ba->overlap_area.push_back(p_g->overlapArea);
    ba->overlap_index.push_back(p_g->overlap_index);
    ba->real_area.push_back(p_g->realArea);
    
}




}

namespace geo_basic{
 
    double EPS = 1e-9;

    bool isIntersectSegmentLine(pLine &s, pLine &l){
        //line s and line l are not parallel to each other
        if(abs(((s.second - s.first) % (l.second - l.first )).z)< EPS){
            return false;
        }
        
        _vec<double> m1 = (l.second - l.first) % (s.first - l.first)  ;
        _vec<double> m2 = (l.second - l.first) % (s.second - l.first) ;
    
        if( m1.z * m2.z > EPS) return false;
        return true;
    }

    _vec<double> crossPoint(geo_basic::pLine &s, geo_basic::pLine &t){
        _vec<double> sv = s.second - s.first;
        _vec<double> tv = t.second - t.first;
        assert((sv % tv).norm() > geo_basic::EPS);
        double length = (tv % (t.first - s.first)).z / (tv % sv).z;
        return s .first + sv * length;
    }
}

namespace geo_vv{
    double maximum_distance_between_points_brute_force(vector<Vertex> vv, Organ* p_g){
        double maximum_distance=0;   
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=1+1;vj<(int)vv.size();vj++){
                double distance_tmp = vertex_geo::vertex_distance(vv[vi],vv[vj]);
                if(maximum_distance<distance_tmp){
                    maximum_distance = distance_tmp;
                    p_g->length_vertex= std::make_pair(vv[vi].loc,vv[vj].loc);
                }
            }
        }
        return maximum_distance;
    }

    double maxmum_distance_between_points_rotating_caliphers(vector<Vertex> vv, Organ* p_g){
        double maximum_distance=0;
        int n=vv.size();
        int j=1;
        for (int i = 0; i < n; ++i) {
        while (vertex_geo::vertex_distance(vv[i], vv[(j + 1) % n]) > vertex_geo::vertex_distance(vv[i], vv[j]))
            j = (j + 1) % n;

        maximum_distance = std::max(maximum_distance, vertex_geo::vertex_distance(vv[i], vv[j]));
        p_g->length_vertex= std::make_pair(vv[i].loc,vv[j].loc);
    }

        return maximum_distance;
    }



    vector<Vertex*> after_ImageJ_process(vector<Vertex*> vv){
        vector<Vertex*> vv_processed;

        //1. normalization
        vector<Vertex*> vv_normalized = geo_vv::normalization(vv);
        
        //2. swap from bottom to top 
        
        for(int vi=0; vi<vv_normalized.size(); vi++){
                vv_normalized[vi]->loc.y=1-vv_normalized[vi]->loc.y;
        }
        
        return vv_normalized;
        
    }

    vector<Vertex*> normalization_by_perimeter(vector<Vertex*> vv){
        vector<Vertex*> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0]->loc.y;
        double y_min = vv[0]->loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y<y_min){
                y_min=vv[vi]->loc.y;
            }
            if(vv[vi]->loc.y>y_max){
                y_max=vv[vi]->loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = perimeter_vv_boundary(vv);
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi]->loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = (vv[vi]->loc.x-x_center)/norm_length;
            v_tmp->loc.y = (vv[vi]->loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi]->loc.y<<";y_norm "<<v_tmp->loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    vector<Vertex*> normalization(vector<Vertex*> vv){
        vector<Vertex*> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0]->loc.y;
        double y_min = vv[0]->loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y<y_min){
                y_min=vv[vi]->loc.y;
            }
            if(vv[vi]->loc.y>y_max){
                y_max=vv[vi]->loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = y_max-y_min;
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi]->loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = (vv[vi]->loc.x-x_center)/norm_length;
            v_tmp->loc.y = (vv[vi]->loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi]->loc.y<<";y_norm "<<v_tmp->loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    vector<Vertex> normalization(vector<Vertex> vv){
        vector<Vertex> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0].loc.y;
        double y_min = vv[0].loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y<y_min){
                y_min=vv[vi].loc.y;
            }
            if(vv[vi].loc.y>y_max){
                y_max=vv[vi].loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = y_max-y_min;
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi].loc.x/(double)vv.size();
        }
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

    vector<Vertex> normalization_y_only(vector<Vertex> vv){
        vector<Vertex> vv_normalized;

        //calculate y_max and y_min
        double y_max = vv[0].loc.y;
        double y_min = vv[0].loc.y;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y<y_min){
                y_min=vv[vi].loc.y;
            }
            if(vv[vi].loc.y>y_max){
                y_max=vv[vi].loc.y;
            }
        }
        //std::cout<<"y_min "<<y_min<<"; y_max "<<y_max<<endl;

        //calculate vv_normalization
        double norm_length = y_max-y_min;
        double x_center=0;
        for(int vi=0;vi<vv.size();vi++){
            x_center+=vv[vi].loc.x/(double)vv.size();
        }
        //std::cout<<"x_center: "<<x_center<<endl;
        //std::cout<<"Norm length: "<<norm_length<<endl;
        for(int vi=0;vi<(int)vv.size();vi++){
            Vertex v_tmp;
            v_tmp.loc.x = (vv[vi].loc.x-x_center);
            v_tmp.loc.y = (vv[vi].loc.y-y_min)/norm_length;
            //std::cout<<"y "<<vv[vi].loc.y<<";y_norm "<<v_tmp.loc.y<<endl;
            vv_normalized.push_back(v_tmp);
        }
        return vv_normalized;
    }

    bool equal_Euclidean_distance_except_last(vector<Vertex*> vv)
    {   
        vector<double> Eu_dist;
        for(int vi=0;vi<(int)vv.size()-1;vi++){
            Eu_dist.push_back(vertex_geo::vertex_distance(vv[vi],vv[vi+1]));
        }

        //compare all element with the first element

        for(int vi=0;vi<(int)Eu_dist.size();vi++){
            if(abs(Eu_dist[vi]-Eu_dist[0])>EPS_geo){
                return false;
            }
        }
        return true;
    }

    vector<double> first_derivative(vector<Vertex*> vv){
        int vn=vv.size()-1;
        vector<double> derivatives;
        //1. use forward difference to calculate the derivative of the first point: 
        //Forward difference formula: dy_i/dx_i ≈ (y_(i+1) - y_i) / (x_(i+1) - x_i)
        //dy_0/dx_0 ≈ (y_(1) - y_0) / (x_(1) - x_0)
        double derivatives_0 = (vv[1]->loc.y-vv[0]->loc.y)/(vv[1]->loc.x-vv[0]->loc.x);
        derivatives.push_back(derivatives_0);
        //2. use central difference to calculate the derivative of the middle points:
        //Central difference formula: dy_i/dx_i ≈ (y_(i+1) - y_(i-1)) / (x_(i+1) - x_(i-1))
        for(int i=1; i<vn-1;++i){
            double derivatives_tmp = (vv[i+1]->loc.y-vv[i-1]->loc.y)/(vv[i+1]->loc.x - vv[i-1]->loc.x);
            derivatives.push_back(derivatives_tmp);
        }
        //3. use the backward difference to calculate the derivative of the last point:
        //Backward difference formula: dy_i/dx_i ≈ (y_i - y_(i-1)) / (x_i - x_(i-1))
        //dy_vn/dx_vn ≈ (y_vn - y_(vn-1)) / (x_vn - x_(vn-1))
        double derivatives_last = (vv[vn]->loc.y-vv[vn-1]->loc.y)/(vv[vn]->loc.x-vv[vn-1]->loc.x);
        derivatives.push_back(derivatives_last);

        return derivatives;
    }

    vector<double> second_derivative(vector<Vertex*> vv){
        int vn=vv.size()-1;
        vector<double> derivatives;
        //1. Forward difference for the second derivative
        //formula: d²y_i/dx_i² ≈ (y_(i+2) - 2*y_(i+1) + y_i) / ((x_(i+1) - x_i) * (x_(i+2) - x_(i+1)))
        //d²y_0/dx_0² ≈ (y_(2) - 2*y_(1) + y_0) / ((x_(1) - x_0) * (x_(2) - x_(1)))
        double derivatives_0 = (vv[2]->loc.y-2*vv[1]->loc.y+vv[0]->loc.y)/((vv[1]->loc.x-vv[0]->loc.x)*(vv[2]->loc.x-vv[1]->loc.x));
        derivatives.push_back(derivatives_0);

        //2. Central difference for the second derivative 
        //Central difference formula: d²y_i/dx_i² ≈ (y_(i+1) - 2*y_i + y_(i-1)) / ((x_(i+1) - x_i) * (x_i - x_(i-1)))
        for(int i=1; i<vn-1;++i){
            double derivatives_tmp = (vv[i+1]->loc.y-2*vv[i]->loc.y+vv[i-1]->loc.y)/((vv[i+1]->loc.x-vv[i]->loc.x)*(vv[i]->loc.x-vv[i-1]->loc.x));
            derivatives.push_back(derivatives_tmp);
        }

        //3. Backward difference for the second derivative
        //Backward difference formula: d²y_i/dx_i² ≈ (y_i - 2*y_(i-1) + y_(i-2)) / ((x_i - x_(i-1)) * (x_(i-1) - x_(i-2)))
        //d²y_vn/dx_vn² ≈ (y_vn - 2*y_(vn-1) + y_(vn-2)) / ((x_vn - x_(vn-1)) * (x_(vn-1) - x_(vn-2)))
        double derivatives_last = (vv[vn]->loc.y-2*vv[vn-1]->loc.y+vv[vn-2]->loc.y)/((vv[vn]->loc.x-vv[vn-1]->loc.x)*(vv[vn-1]->loc.x-vv[vn-2]->loc.x));
        derivatives.push_back(derivatives_last);

        return derivatives;
    }

    vector<double> second_derivative_fourth_order(vector<Vertex*> vv){
        int vn=vv.size()-1;
        
        vector<double> derivatives;
        //1. forward difference for the second derivative 
        //forward difference formua: (2*y_vn - 5*y_(vn+1) + 4*y_(vn+2) -y_(vn+3))/(dx^3)
        double dx_0 = (vv[3]->loc.x - vv[0]->loc.x)/4.0;
        double derivatives_0 = (2*vv[0]->loc.y-5*vv[1]->loc.y + 4*vv[2]->loc.y - vv[3]->loc.y)/(dx_0*dx_0*dx_0);
        derivatives.push_back(derivatives_0);
        double dx_1 = (vv[4]->loc.x - vv[1]->loc.x)/4.0;
        double derivatives_1 = (2*vv[1]->loc.y-5*vv[2]->loc.y + 4*vv[3]->loc.y - vv[4]->loc.y)/(dx_1*dx_1*dx_1);
        derivatives.push_back(derivatives_1);
        //2. central difference for the second derivative
        //formula: (-y_(vn+2)+16y_(vn+1) - 30 y_vn + 16 y_(vn-1) - y_(vn-2))/(12dx^2)
        for(int i=2;i<vn-2;++i){
            double dx_i = (vv[i+2]->loc.x-vv[i-2]->loc.x)/5.0; 
            double derivatives_tmp = (-vv[i+2]->loc.y + 16*vv[i+1]->loc.y - 30*vv[i]->loc.y + 16*vv[i-1]->loc.y - vv[i-2]->loc.y)/(12*dx_i*dx_i);
            derivatives.push_back(derivatives_tmp);
        }
        //3. backward difference for the second derivative
        //formula: (2*y_vn - 5*y_(vn-1) + 4*y_(vn-2) -y_(vn-3))/(dx^3)

        double dx_last_2 = (vv[vn-2]->loc.x-vv[vn-5]->loc.x)/4.0;
        double derivatives_last_2 = (2*vv[vn-2]->loc.y - 5*vv[vn-3]->loc.y + 4*vv[vn-4]->loc.y - vv[vn-5]->loc.y)/(dx_last_2*dx_last_2*dx_last_2);
        //derivatives.push_back(dx_last_2);
        double dx_last_1 = (vv[vn-1]->loc.x-vv[vn-4]->loc.x)/4.0;
        double derivatives_last_1 = (2*vv[vn-1]->loc.y - 5*vv[vn-2]->loc.y + 4*vv[vn-3]->loc.y - vv[vn-4]->loc.y)/(dx_last_1*dx_last_1*dx_last_1);
        //derivatives.push_back(dx_last_1);

        return derivatives;
    }

    //the vv should be ordered whether clockwise or anticlockwise
    double area_vv_boundary(vector<Vertex*> vv){
        double area =0.0;
        int n=vv.size();
        for(int i=0; i<n; ++i){
            int next = (i+1)%n;
            area += vv[i]->loc.x*vv[next]->loc.y-vv[next]->loc.x*vv[i]->loc.y;
        }
        return abs(area)/2.0;
    }

    double perimeter_vv_boundary(vector<Vertex*> vv){
        double perimeter =0.0;
        int n = vv.size();
        for(int i=0; i<n; ++i){
            int next=(i+1)%n;
            perimeter += vertex_geo::vertex_distance(vv[i],vv[next]);
            //std::cout<<"perimeter "<<i<<" "<<perimeter<<endl;
        }
        return perimeter;
    }

    double perimeter_vv_boundary(vector<Vertex> vv){
        double perimeter =0.0;
        int n = vv.size();
        for(int i=0; i<n; ++i){
            int next=(i+1)%n;
            perimeter += vertex_geo::vertex_distance(vv[i],vv[next]);
            //std::cout<<"perimeter "<<i<<" "<<perimeter<<endl;
        }
        return perimeter;
    }


    vector<Vertex> organ_boundary_points_along_polygon(vector<Vertex> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex> boundary_points;
        double perimeter_vv = geo_vv::perimeter_vv_boundary(anticlockwise_surface_vertices);
        double boundary_points_average_distance = perimeter_vv/boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex boundary_points_new;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new.loc = boundary_points[pi].loc+(anticlockwise_surface_vertices[vj+1].loc-boundary_points[pi].loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new.loc = anticlockwise_surface_vertices[vj].loc+(anticlockwise_surface_vertices[vj+1].loc-anticlockwise_surface_vertices[vj].loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }

    vector<Vertex*> organ_boundary_points_along_polygon(vector<Vertex*> anticlockwise_surface_vertices, int boundary_points_number){
        vector<Vertex*> boundary_points;
        double perimeter_vv = geo_vv::perimeter_vv_boundary(anticlockwise_surface_vertices);
        double boundary_points_average_distance = perimeter_vv/boundary_points_number;

        //now we have the ordered boundary vertices: the boundary vertex with lowest y is ranked as the first; the other boundary vertices are ordered anticlockwise
    //1. the first boundary point: just pick the first boundary vertex (with lowest y)
        boundary_points.push_back(anticlockwise_surface_vertices[0]);
        anticlockwise_surface_vertices.push_back(anticlockwise_surface_vertices[0]); //link the last boundary vertices with the first boundary vertices
        
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    //2. the i-th boundary point: incremental algorithm

        double required_boundary_points_distance = boundary_points_average_distance;
        int vj = 0; //the current boundary_vertex index
        bool distance_intial_is_boundary_vertex=0;
        double surface_vertex_distance;

    for(int pi=0; pi<boundary_points_number-1; pi++){
        
            boundary_points_distance_larger_than_surface_vertex_distance:;
            if(distance_intial_is_boundary_vertex==0) {
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc).norm();
            }
            else{
                surface_vertex_distance = (anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc).norm();
        }


        //judgement: if boundary_point_average_distance >=< surface_vertex_distance
        if(required_boundary_points_distance < surface_vertex_distance){
        // the new boundary point will be generated on the line segment between surface vertex i and surface vertex i+1
            Vertex* boundary_points_new = new Vertex;
            if(distance_intial_is_boundary_vertex==0) {
                boundary_points_new->loc = boundary_points[pi]->loc+(anticlockwise_surface_vertices[vj+1]->loc-boundary_points[pi]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }
            else{
                boundary_points_new->loc = anticlockwise_surface_vertices[vj]->loc+(anticlockwise_surface_vertices[vj+1]->loc-anticlockwise_surface_vertices[vj]->loc)*required_boundary_points_distance/surface_vertex_distance;
            }

            boundary_points.push_back(boundary_points_new);
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<","<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance == surface_vertex_distance){
        // the new boundary point will be exactly the suface vertex i+1
            boundary_points.push_back(anticlockwise_surface_vertices[vj+1]);
            vj++;
            //std::cout<<pi<<" location: ("<<boundary_points[pi]->loc.x<<","<<boundary_points[pi]->loc.y<<")."<<"; on line "<<vj<<vj+1<<"; d: "<<required_boundary_points_distance<<"; |vj+1vj|: "<<surface_vertex_distance<<endl;
        }
        else if(required_boundary_points_distance > surface_vertex_distance){
        // the new boundary point might be generated on the line segment between surface vertex i+1 and surface vertex i+2
            required_boundary_points_distance -= surface_vertex_distance;
            vj++;
            distance_intial_is_boundary_vertex=1;
            goto boundary_points_distance_larger_than_surface_vertex_distance;
        }
        
        required_boundary_points_distance = boundary_points_average_distance;
        distance_intial_is_boundary_vertex=0;
        //cout_fout_debug::cout_vector_vertex(boundary_points);

    }

        //cout_fout_debug::fout_vector_vertex(boundary_points,"50_boundary_points.txt");
        return boundary_points;

    }

    bool point_in_polygon_ray_casting(Vertex* point_tmp, vector<Vertex*> vv){
        bool inside_bool=0;
        int ray_casting_index=0;
        int n = vv.size();
        for(int i=0; i<n; ++i){
            int next = (i+1)%n;
            double y0_tmp = vv[i]->loc.y;
            double x0_tmp = vv[i]->loc.x;
            double y1_tmp = vv[next]->loc.y;
            double x1_tmp = vv[next]->loc.x;
            double y_point = point_tmp->loc.y;
            double x_point = point_tmp->loc.x;
            if((y0_tmp>y_point)!=(y1_tmp>y_point)){
                if(x_point<(x0_tmp+(y_point-y0_tmp)/(y1_tmp-y0_tmp)*(x1_tmp-x0_tmp))){
                    ray_casting_index++;
                    //std::cout<<"Ray casting passed line "<<li<<"; p0 "<<x0_tmp<<","<<y0_tmp<<"; p1 "<<x1_tmp<<","<<y1_tmp<<endl;
                }
            }

        }

        if(ray_casting_index%2==0){
            //even
            inside_bool = 0;
        }
        if(ray_casting_index%2==1){
            //odd
            inside_bool = 1;
        }
        //std::cout<<"ray_casting_index "<<ray_casting_index<<endl;
        return inside_bool;
    }
    
    bool point_in_polygon_ray_casting(Vertex point_tmp, vector<Vertex> vv){
        bool inside_bool=0;
        int ray_casting_index=0;
        int n = vv.size();
        for(int i=0; i<n; ++i){
            int next = (i+1)%n;
            double y0_tmp = vv[i].loc.y;
            double x0_tmp = vv[i].loc.x;
            double y1_tmp = vv[next].loc.y;
            double x1_tmp = vv[next].loc.x;
            double y_point = point_tmp.loc.y;
            double x_point = point_tmp.loc.x;
            if((y0_tmp>y_point)!=(y1_tmp>y_point)){
                if(x_point<(x0_tmp+(y_point-y0_tmp)/(y1_tmp-y0_tmp)*(x1_tmp-x0_tmp))){
                    ray_casting_index++;
                    //std::cout<<"Ray casting passed line "<<li<<"; p0 "<<x0_tmp<<","<<y0_tmp<<"; p1 "<<x1_tmp<<","<<y1_tmp<<endl;
                }
            }

        }

        if(ray_casting_index%2==0){
            //even
            inside_bool = 0;
        }
        if(ray_casting_index%2==1){
            //odd
            inside_bool = 1;
        }
        //std::cout<<"ray_casting_index "<<ray_casting_index<<endl;
        return inside_bool;
    }
    
    //do calculate the circle fitting of all boundary points by kasa
    vector<double> curvature_circle_fitting_kasa_three_boundary_points(vector<Vertex*> boundary_points, int points_away, int NumAveraging){

        vector<double> curvature;
        vector<double> curvature_averaged;

        int boundary_points_number = (int)boundary_points.size();

    //1. start circle fitting (Kasa)
    for(int vi=0;vi<(int)boundary_points.size();vi++){
        Vertex* p_minus = new Vertex;
        Vertex* p_0 = new Vertex;
        Vertex* p_plus = new Vertex;

        p_0 = boundary_points[vi];
        if(vi<points_away){
            p_minus=boundary_points[vi-points_away+boundary_points_number];
        }
        else{
            p_minus=boundary_points[vi-points_away];
        }

        if(vi>(boundary_points_number-points_away-1)){
            p_plus=boundary_points[vi+points_away-boundary_points_number];
        }
        else{
            p_plus=boundary_points[vi+points_away];
        }
        vector<Vertex*> fitting_points;
        fitting_points.push_back(p_minus);
        fitting_points.push_back(p_0);
        fitting_points.push_back(p_plus);
        Vertex* midpoint_tmp = new Vertex;
        midpoint_tmp->loc.x = (p_minus->loc.x + p_plus->loc.x)/2.0;
        midpoint_tmp->loc.y = (p_minus->loc.y + p_plus->loc.y)/2.0;
        double curvature_abstract_value = boundary_geo::surface_vertex_curvature_kasa(fitting_points);
        double curvature_value=0;
        
        if(geo_vv::point_in_polygon_ray_casting(midpoint_tmp,boundary_points)==0){
            curvature_value= -curvature_abstract_value;
        }
        else{
            curvature_value = curvature_abstract_value;
        }
    
        curvature.push_back(curvature_value);
    }

    //2. taking average of curvature value
        if(NumAveraging==1){
            return curvature;
        }
        else{
            curvature_averaged = geo_vv::vd_averaged(curvature,NumAveraging);
            return curvature_averaged;
        }    
    }

    vector<double> curvature_circle_fitting_kasa_three_boundary_points(vector<Vertex> boundary_points, int points_away, int NumAveraging){

        vector<double> curvature;
        vector<double> curvature_averaged;

        int boundary_points_number = (int)boundary_points.size();

    //1. start circle fitting (Kasa)
    for(int vi=0;vi<(int)boundary_points.size();vi++){
        Vertex p_minus;
        Vertex p_0;
        Vertex p_plus;

        p_0 = boundary_points[vi];
        if(vi<points_away){
            p_minus=boundary_points[vi-points_away+boundary_points_number];
        }
        else{
            p_minus=boundary_points[vi-points_away];
        }

        if(vi>(boundary_points_number-points_away-1)){
            p_plus=boundary_points[vi+points_away-boundary_points_number];
        }
        else{
            p_plus=boundary_points[vi+points_away];
        }
        vector<Vertex> fitting_points;
        fitting_points.push_back(p_minus);
        fitting_points.push_back(p_0);
        fitting_points.push_back(p_plus);
        Vertex midpoint_tmp;
        midpoint_tmp.loc.x = (p_minus.loc.x + p_plus.loc.x)/2.0;
        midpoint_tmp.loc.y = (p_minus.loc.y + p_plus.loc.y)/2.0;
        double curvature_abstract_value = boundary_geo::surface_vertex_curvature_kasa(fitting_points);
        double curvature_value=0;
        
        if(geo_vv::point_in_polygon_ray_casting(midpoint_tmp,boundary_points)==0){
            curvature_value= -curvature_abstract_value;
        }
        else{
            curvature_value = curvature_abstract_value;
        }
    
        curvature.push_back(curvature_value);
    }

    //2. taking average of curvature value
        if(NumAveraging==1){
            return curvature;
        }
        else{
            curvature_averaged = geo_vv::vd_averaged(curvature,NumAveraging);
            return curvature_averaged;
        }    
    }


    double vd_minimum(vector<double> vd){
        double min_vd=vd[0];
        for(int i=1; i<vd.size();i++){
            if(min_vd>vd[i]){
                min_vd=vd[i];
                //std::cout<<i<<" "<<min_vd<<endl;
            }
        }
        return min_vd;
    }

    double vd_maximum(vector<double> vd){
        double max_vd=vd[0];
        for(int i=1; i<vd.size();i++){
            if(max_vd<vd[i]){
                max_vd=vd[i];
                //std::cout<<i<<" "<<max_vd<<endl;
            }
        }
        return max_vd;
    }

    double accumulated_negative(vector<double> vd){
        double accumulated=0;
        for(int i=0; i<vd.size(); i++){
            if(vd[i]<0){
                accumulated+=vd[i]/(double)vd.size();
            }
        }
        return accumulated;
    }

    vector<double> vd_averaged(vector<double> vd, int numNeighbors){
        vector<double> output_vd(vd.size());

        int halfNeighbors = (numNeighbors-1)/2;

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

            output_vd[i] = sum/numNeighbors;
        }

        return output_vd;
    }

    vector<Vertex*> vector_vertex_sampling(vector<Vertex*> vv,double sampling_distance){
        vector<Vertex*> vv_sampled;
        //divide the organ into left side and right side based on the theta 
        //when sampling_distance = 0.01, if theta_0<theta_100, left side, else right side
        vector<double> theta;
        vector<Vertex*> left_side_vertex;
        vector<Vertex*> right_side_vertex;
        //find center of outline vertices
        double center_x=0,center_y=0;
        for(int i=0;i<(int)vv.size();i++){
            center_x = vv[i]->loc.x/(double)vv.size();
            center_y = vv[i]->loc.y/(double)vv.size();
        }
        //std::cout<<"center_x "<<center_x<<" center_y "<<center_y<<endl;
        for(int i=0;i<(int)vv.size();i++){
            double x_tmp=vv[i]->loc.x-center_x;
            double y_tmp=vv[i]->loc.y-center_y;
            double theta_tmp = atan2(y_tmp,x_tmp);
            if(theta_tmp<0){

            }
            theta.push_back(theta_tmp);
        }

        double x_0,x_100,theta_0,theta_100;
        int vi_0,vi_100;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi]->loc.y==0.00){
                x_0 = vv[vi]->loc.x;
                theta_0 = theta[vi];
                vi_0=vi;
            }
            else if(vv[vi]->loc.y==1.00){
                x_100 = vv[vi]->loc.x;
                theta_100 = theta[vi];
                vi_100=vi;
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                left_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0||vi==vi_100){
                
            }
            else if(theta[vi]>theta_0&&theta[vi]<theta_100){
                left_side_vertex.push_back(vv[vi]);
            }
            else{
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                right_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                left_side_vertex.push_back(vv[vi]);
            }
        }

        //sorting left_side and right_side by y value
        left_side_vertex = sort_vector_vertex_ascend(left_side_vertex);
        right_side_vertex = sort_vector_vertex_descend(right_side_vertex);

        Vertex* v_0 = new Vertex;
        Vertex* v_100 = new Vertex;
        v_0->loc.x = x_0;
        v_0->loc.y =0;
        v_100->loc.x=x_100;
        v_100->loc.y=1.0;

        vv_sampled.push_back(v_0);
        for(int i=1;i<1/sampling_distance;i++){
            for(int vi=0;vi<left_side_vertex.size();vi++){
                //std::cout<<(double)i*sampling_distance<<" "<<left_side_vertex[vi]->loc.y<<" "<<left_side_vertex[vi+1]->loc.y<<endl;
                if(left_side_vertex[vi]->loc.y==(double)sampling_distance*i){
                    vv_sampled.push_back(left_side_vertex[vi]);
                }
                else if(left_side_vertex[vi]->loc.y<(double)i*sampling_distance&&left_side_vertex[vi+1]->loc.y>(double)i*sampling_distance){
                    //std::cout<<"fit"<<endl;
                    Vertex* vt_tmp = new Vertex;
                    vt_tmp->loc.y = i*sampling_distance;
                    vt_tmp->loc.x = linear_fitting(left_side_vertex[vi],left_side_vertex[vi+1],i*sampling_distance);
                    //vt_tmp->loc.x = i*sampling_distance;
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        vv_sampled.push_back(v_100);

        for(int i=1/sampling_distance+1;i<2/sampling_distance;i++){
            for(int vi=0;vi<right_side_vertex.size();vi++){
                if(right_side_vertex[vi]->loc.y==(2-sampling_distance*i)){
                    vv_sampled.push_back(right_side_vertex[vi]);
                }
                else if(right_side_vertex[vi]->loc.y>(2-i*sampling_distance)&&right_side_vertex[vi+1]->loc.y<(2-i*sampling_distance)){
                    Vertex* vt_tmp = new Vertex;
                    vt_tmp->loc.y = 2-i*sampling_distance;
                    vt_tmp->loc.x = linear_fitting(right_side_vertex[vi],right_side_vertex[vi+1],2-i*sampling_distance);
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        /*
        ofstream fout3("sampled.txt");
        fout3<<"x y"<<endl;
        for(int i=0;i<vv_sampled.size();i++){
            fout3<<vv_sampled[i]->loc.x<<" "<<vv_sampled[i]->loc.y<<endl;
        }
        fout3.close();
        
        
        ofstream fout("sampled_left.txt");
        fout<<"x y"<<endl;
        for(int i=0;i<left_side_vertex.size();i++){
            fout<<left_side_vertex[i]->loc.x<<" "<<left_side_vertex[i]->loc.y<<endl;
        }
        fout.close();

        ofstream fout2("sampled_right.txt");
        fout2<<"x y"<<endl;
        for(int i=0;i<right_side_vertex.size();i++){
            fout2<<right_side_vertex[i]->loc.x<<" "<<right_side_vertex[i]->loc.y<<endl;
        }
        fout2.close();
        */
        return vv_sampled;
    }

     vector<Vertex> vector_vertex_sampling__(vector<Vertex> vv,double sampling_distance){
        vector<Vertex> vv_sampled;
        //divide the organ into left side and right side based on the theta 
        //when sampling_distance = 0.01, if theta_0<theta_100, left side, else right side
        vector<double> theta;
        vector<Vertex> left_side_vertex;
        vector<Vertex> right_side_vertex;
        //find center of outline vertices
        double center_x=0,center_y=0;
        for(int i=0;i<(int)vv.size();i++){
            center_x = vv[i].loc.x/(double)vv.size();
            center_y = vv[i].loc.y/(double)vv.size();
        }
        //std::cout<<"center_x "<<center_x<<" center_y "<<center_y<<endl;
        for(int i=0;i<(int)vv.size();i++){
            double x_tmp=vv[i].loc.x-center_x;
            double y_tmp=vv[i].loc.y-center_y;
            double theta_tmp = atan2(y_tmp,x_tmp);
            if(theta_tmp<0){

            }
            theta.push_back(theta_tmp);
        }

        double x_0,x_100,theta_0,theta_100;
        int vi_0,vi_100;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y==0.00){
                x_0 = vv[vi].loc.x;
                theta_0 = theta[vi];
                vi_0=vi;
            }
            else if(vv[vi].loc.y==1.00){
                x_100 = vv[vi].loc.x;
                theta_100 = theta[vi];
                vi_100=vi;
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                left_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0||vi==vi_100){
                
            }
            else if(theta[vi]>theta_0&&theta[vi]<theta_100){
                left_side_vertex.push_back(vv[vi]);
            }
            else{
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                right_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                left_side_vertex.push_back(vv[vi]);
            }
        }

        //sorting left_side and right_side by y value
        left_side_vertex = sort_vector_vertex_ascend(left_side_vertex);
        right_side_vertex = sort_vector_vertex_descend(right_side_vertex);

        Vertex v_0;
        Vertex v_100;
        v_0.loc.x = x_0;
        v_0.loc.y =0;
        v_100.loc.x=x_100;
        v_100.loc.y=1.0;

        vv_sampled.push_back(v_0);
        for(int i=1;i<100;i++){
            for(int vi=0;vi<left_side_vertex.size();vi++){
                //std::cout<<(double)i*0.01<<" "<<left_side_vertex[vi]->loc.y<<" "<<left_side_vertex[vi+1]->loc.y<<endl;
                if(left_side_vertex[vi].loc.y==(double)0.01*i){
                    vv_sampled.push_back(left_side_vertex[vi]);
                }
                else if(left_side_vertex[vi].loc.y<(double)i*0.01&&left_side_vertex[vi+1].loc.y>(double)i*0.01){
                    //std::cout<<"fit"<<endl;
                    Vertex vt_tmp;
                    vt_tmp.loc.y = i*0.01;
                    vt_tmp.loc.x = linear_fitting(left_side_vertex[vi],left_side_vertex[vi+1],i*0.01);
                    //vt_tmp->loc.x = i*0.01;
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        vv_sampled.push_back(v_100);

        for(int i=101;i<200;i++){
            for(int vi=0;vi<right_side_vertex.size();vi++){
                if(right_side_vertex[vi].loc.y==(2-0.01*i)){
                    vv_sampled.push_back(right_side_vertex[vi]);
                }
                else if(right_side_vertex[vi].loc.y>(2-i*0.01)&&right_side_vertex[vi+1].loc.y<(2-i*0.01)){
                    Vertex vt_tmp;
                    vt_tmp.loc.y = 2-i*0.01;
                    vt_tmp.loc.x = linear_fitting(right_side_vertex[vi],right_side_vertex[vi+1],2-i*0.01);
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        /*
        ofstream fout3("sampled.txt");
        fout3<<"x y"<<endl;
        for(int i=0;i<vv_sampled.size();i++){
            fout3<<vv_sampled[i]->loc.x<<" "<<vv_sampled[i]->loc.y<<endl;
        }
        fout3.close();
        
        
        ofstream fout("sampled_left.txt");
        fout<<"x y"<<endl;
        for(int i=0;i<left_side_vertex.size();i++){
            fout<<left_side_vertex[i]->loc.x<<" "<<left_side_vertex[i]->loc.y<<endl;
        }
        fout.close();

        ofstream fout2("sampled_right.txt");
        fout2<<"x y"<<endl;
        for(int i=0;i<right_side_vertex.size();i++){
            fout2<<right_side_vertex[i]->loc.x<<" "<<right_side_vertex[i]->loc.y<<endl;
        }
        fout2.close();
        */
        return vv_sampled;
    }

    vector<Vertex> vector_vertex_sampling(vector<Vertex> vv,double sampling_distance){
        std::vector<Vertex> sampledVertex;
        int numSamples = 2/sampling_distance;
        for(int i=0;i<numSamples/2;++i){
            double y = i*sampling_distance;
            for( size_t j=0; j<vv.size()/2; j++){
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
            for( size_t j=vv.size()/2; j<vv.size(); j++){
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

    vector<Vertex> vector_vertex_sampling_(vector<Vertex> vv,double sampling_distance){
        vector<Vertex> vv_sampled;
        //divide the organ into left side and right side based on the theta 
        //when sampling_distance = sampling_distance, if theta_0<theta_100, left side, else right side
        vector<double> theta;
        vector<Vertex> left_side_vertex;
        vector<Vertex> right_side_vertex;
        //find center of outline vertices
        double center_x=0,center_y=0;
        for(int i=0;i<(int)vv.size();i++){
            center_x = vv[i].loc.x/(double)vv.size();
            center_y = vv[i].loc.y/(double)vv.size();
        }
        //std::cout<<"center_x "<<center_x<<" center_y "<<center_y<<endl;
        for(int i=0;i<(int)vv.size();i++){
            double x_tmp=vv[i].loc.x-center_x;
            double y_tmp=vv[i].loc.y-center_y;
            double theta_tmp = atan2(y_tmp,x_tmp);
            if(theta_tmp<0){

            }
            theta.push_back(theta_tmp);
        }

        double x_0,x_100,theta_0,theta_100;
        int vi_0,vi_100;
        for(int vi=0;vi<vv.size();vi++){
            if(vv[vi].loc.y==0.00){
                x_0 = vv[vi].loc.x;
                theta_0 = theta[vi];
                vi_0=vi;
            }
            else if(vv[vi].loc.y==1.00){
                x_100 = vv[vi].loc.x;
                theta_100 = theta[vi];
                vi_100=vi;
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                left_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0||vi==vi_100){
                
            }
            else if(theta[vi]>theta_0&&theta[vi]<theta_100){
                left_side_vertex.push_back(vv[vi]);
            }
            else{
                right_side_vertex.push_back(vv[vi]);
            }
        }

        for(int vi=0;vi<vv.size();vi++){
            if(vi==vi_0){
                right_side_vertex.push_back(vv[vi]);
            }
            else if(vi==vi_100){
                left_side_vertex.push_back(vv[vi]);
            }
        }

        //sorting left_side and right_side by y value
        left_side_vertex = sort_vector_vertex_ascend(left_side_vertex);
        right_side_vertex = sort_vector_vertex_descend(right_side_vertex);

        Vertex v_0;
        Vertex v_100;
        v_0.loc.x = x_0;
        v_0.loc.y =0;
        v_100.loc.x=x_100;
        v_100.loc.y=1.0;

        vv_sampled.push_back(v_0);
        for(int i=1;i<1/sampling_distance;i++){
            for(int vi=0;vi<left_side_vertex.size();vi++){
                //std::cout<<(double)i*sampling_distance<<" "<<left_side_vertex[vi]->loc.y<<" "<<left_side_vertex[vi+1]->loc.y<<endl;
                if(left_side_vertex[vi].loc.y==(double)sampling_distance*i){
                    vv_sampled.push_back(left_side_vertex[vi]);
                }
                else if(left_side_vertex[vi].loc.y<(double)i*sampling_distance&&left_side_vertex[vi+1].loc.y>(double)i*sampling_distance){
                    //std::cout<<"fit"<<endl;
                    Vertex vt_tmp;
                    vt_tmp.loc.y = i*sampling_distance;
                    vt_tmp.loc.x = linear_fitting(left_side_vertex[vi],left_side_vertex[vi+1],i*sampling_distance);
                    //vt_tmp->loc.x = i*sampling_distance;
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        vv_sampled.push_back(v_100);

        for(int i=1/sampling_distance+1;i<2/sampling_distance;i++){
            for(int vi=0;vi<right_side_vertex.size();vi++){
                if(right_side_vertex[vi].loc.y==(2-sampling_distance*i)){
                    vv_sampled.push_back(right_side_vertex[vi]);
                }
                else if(right_side_vertex[vi].loc.y>(2-i*sampling_distance)&&right_side_vertex[vi+1].loc.y<(2-i*sampling_distance)){
                    Vertex vt_tmp;
                    vt_tmp.loc.y = 2-i*sampling_distance;
                    vt_tmp.loc.x = linear_fitting(right_side_vertex[vi],right_side_vertex[vi+1],2-i*sampling_distance);
                    vv_sampled.push_back(vt_tmp);
                }
                else{
                    //std::cout<<"Error in finding fitted simulated outline"<<endl;
                }
            }
        }
        /*
        ofstream fout3("sampled.txt");
        fout3<<"x y"<<endl;
        for(int i=0;i<vv_sampled.size();i++){
            fout3<<vv_sampled[i]->loc.x<<" "<<vv_sampled[i]->loc.y<<endl;
        }
        fout3.close();
        
        
        ofstream fout("sampled_left.txt");
        fout<<"x y"<<endl;
        for(int i=0;i<left_side_vertex.size();i++){
            fout<<left_side_vertex[i]->loc.x<<" "<<left_side_vertex[i]->loc.y<<endl;
        }
        fout.close();

        ofstream fout2("sampled_right.txt");
        fout2<<"x y"<<endl;
        for(int i=0;i<right_side_vertex.size();i++){
            fout2<<right_side_vertex[i]->loc.x<<" "<<right_side_vertex[i]->loc.y<<endl;
        }
        fout2.close();
        */
        return vv_sampled;
    }
    
    vector<Vertex*> sort_vector_vertex_ascend(vector<Vertex*> vv){
        vector<Vertex*> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi]->loc.y);
        }

        sort(sort_y.begin(),sort_y.end());
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj]->loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex> sort_vector_vertex_ascend(vector<Vertex> vv){
        vector<Vertex> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi].loc.y);
        }

        sort(sort_y.begin(),sort_y.end());
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj].loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex*> sort_vector_vertex_descend(vector<Vertex*> vv){
        vector<Vertex*> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi]->loc.y);
        }

        sort(sort_y.begin(),sort_y.end(),comp_descend);
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj]->loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    vector<Vertex> sort_vector_vertex_descend(vector<Vertex> vv){
        vector<Vertex> sort_vertex;

        vector<double> sort_y;
        for(int vi=0;vi<(int)vv.size();vi++){
            sort_y.push_back(vv[vi].loc.y);
        }

        sort(sort_y.begin(),sort_y.end(),comp_descend);
        for(int vi=0;vi<(int)vv.size();vi++){
            for(int vj=0;vj<(int)vv.size();vj++){
                if(sort_y[vi]==vv[vj].loc.y){
                    sort_vertex.push_back(vv[vj]);
                }
            }
        }

        return sort_vertex;
    }

    bool comp_descend(double y1,double y2){
        return (y1>y2);
    }

    double linear_fitting(Vertex* v1, Vertex* v2, double y_i){
        double x_i;
        if(v1->loc.x == v2->loc.x){
            x_i = v1->loc.x;
        }
        else{
            double a_tmp = (v1->loc.y-v2->loc.y)/(v1->loc.x-v2->loc.x);
            double b_tmp = (v1->loc.x*v2->loc.y-v2->loc.x*v1->loc.y)/(v1->loc.x-v2->loc.x);
            x_i = (y_i-b_tmp)/a_tmp;
        }
        //std::cout<<"a_tmp "<<a_tmp<<" b_tmp "<<b_tmp<<" x_i "<<x_i<<endl;
        return x_i;
    }

    double linear_fitting(Vertex v1, Vertex v2, double y_i){
        double x_i;
        if(v1.loc.x == v2.loc.x){
            x_i = v1.loc.x;
        }
        else{
            double a_tmp = (v1.loc.y-v2.loc.y)/(v1.loc.x-v2.loc.x);
            double b_tmp = (v1.loc.x*v2.loc.y-v2.loc.x*v1.loc.y)/(v1.loc.x-v2.loc.x);
            x_i = (y_i-b_tmp)/a_tmp;
        }
        //std::cout<<"a_tmp "<<a_tmp<<" b_tmp "<<b_tmp<<" x_i "<<x_i<<endl;
        return x_i;
    }

    vector<Vertex*> vv_x_swap(vector<Vertex*> vv){
        vector<Vertex*> vv_swaped;
        for(int i=0;i<vv.size();i++){
            Vertex* v_tmp = new Vertex;
            v_tmp->loc.x = -vv[i]->loc.x;
            v_tmp->loc.y = vv[i]->loc.y;
            vv_swaped.push_back(v_tmp);
        }
        return vv_swaped;
    }

    vector<Vertex> vv_x_swap(vector<Vertex> vv){
        vector<Vertex> vv_swaped;
        for(int i=0;i<vv.size();i++){
            Vertex v_tmp;
            v_tmp.loc.x = -vv[i].loc.x;
            v_tmp.loc.y = vv[i].loc.y;
            vv_swaped.push_back(v_tmp);
        }
        return vv_swaped;
    }

    
}

namespace boundary_geo{
    //calcualte circle fitting of a vector of points
    double surface_vertex_curvature_kasa(vector<Vertex*> vv_fitting){
        //reference: 
        //1. http://www.ne.jp/asahi/paleomagnetism.rock-magnetism/basics/pmag/circ/circ1E.html
        //2. https://people.cas.uab.edu/~mosya/cl/CircleFitByKasa.cpp

        /*
        Vertex* p1_tmp = new Vertex;
        p1_tmp->loc.x = 3.0;
        p1_tmp->loc.y = 3.0;
        Vertex* p2_tmp = new Vertex;
        p2_tmp->loc.x = 2.0;
        
        p2_tmp->loc.y = 4.0;
        Vertex* p3_tmp = new Vertex;
        p3_tmp->loc.x = 2+sqrt(2.0)/2.0;
        p3_tmp->loc.y = 3+sqrt(2.0)/2.0;
        vv_fitting.push_back(p1_tmp);
        vv_fitting.push_back(p2_tmp);
        vv_fitting.push_back(p3_tmp);
        */
        

        //generating the matrix
        double Mx=0, My=0, Mz=0, Mxy=0, Mxz=0, Myz=0, Mxx=0, Myy=0;
        double B_fit, C_fit, D_fit;
        double a_fit, b_fit, R_fit;
        int n_fit = (int)vv_fitting.size();

        for(int vi=0; vi<(int)vv_fitting.size(); vi++){
            double x_tmp = vv_fitting[vi]->loc.x;
            double y_tmp = vv_fitting[vi]->loc.y;
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
        /*
        for(int vi=0;vi<(int)vv_fitting.size();vi++){
            std::cout<<"fitting points "<<vi<<" ("<<vv_fitting[vi]->loc.x<<","<<vv_fitting[vi]->loc.y<<")."<<endl;
        }

        std::cout<<"fitted circle center is ("<<a_fit<<","<<b_fit<<") and fitted radius is "<<R_fit<<endl;
        */
        return 1.0/R_fit;
    }

    double surface_vertex_curvature_kasa(vector<Vertex> vv_fitting){
        //reference: 
        //1. http://www.ne.jp/asahi/paleomagnetism.rock-magnetism/basics/pmag/circ/circ1E.html
        //2. https://people.cas.uab.edu/~mosya/cl/CircleFitByKasa.cpp

        /*
        Vertex* p1_tmp = new Vertex;
        p1_tmp->loc.x = 3.0;
        p1_tmp->loc.y = 3.0;
        Vertex* p2_tmp = new Vertex;
        p2_tmp->loc.x = 2.0;
        p2_tmp->loc.y = 4.0;
        Vertex* p3_tmp = new Vertex;
        p3_tmp->loc.x = 2+sqrt(2.0)/2.0;
        p3_tmp->loc.y = 3+sqrt(2.0)/2.0;
        vv_fitting.push_back(p1_tmp);
        vv_fitting.push_back(p2_tmp);
        vv_fitting.push_back(p3_tmp);
        */
        

        //generating the matrix
        double Mx=0, My=0, Mz=0, Mxy=0, Mxz=0, Myz=0, Mxx=0, Myy=0;
        double B_fit, C_fit, D_fit;
        double a_fit, b_fit, R_fit;
        int n_fit = (int)vv_fitting.size();

        for(int vi=0; vi<(int)vv_fitting.size(); vi++){
            double x_tmp = vv_fitting[vi].loc.x;
            double y_tmp = vv_fitting[vi].loc.y;
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
        /*
        for(int vi=0;vi<(int)vv_fitting.size();vi++){
            std::cout<<"fitting points "<<vi<<" ("<<vv_fitting[vi]->loc.x<<","<<vv_fitting[vi]->loc.y<<")."<<endl;
        }

        std::cout<<"fitted circle center is ("<<a_fit<<","<<b_fit<<") and fitted radius is "<<R_fit<<endl;
        */
        return 1.0/R_fit;
    }


    double similarity_Index_1(vector<Vertex*> vv_simu,vector<Vertex*> vv_real, double sampling_distance){
        //adjust the real y_position
        /*
        for(int vi=0;vi<(int)vv_real.size();vi++){
            vv_real[vi]->loc.y = abs(1-vv_real[vi]->loc.y);
        }
        vector<Vertex*> vv_real_align;
        for(int i=0;i<101;i++){
            vv_real_align.push_back(vv_real[100-i]);
        }
        for(int i=101;i<200;i++){
            vv_real_align.push_back(vv_real[i-100]);
        }
        //output::vv_foutXY(vv_real_align,"real_aligned.txt");
        */
        double error_sum=0;
        for(int vi=0;vi<(int)vv_simu.size();vi++){
            double x_error = abs(vv_simu[vi]->loc.x-vv_real[vi]->loc.x);
            error_sum += x_error;
            //std::cout<<vv_simu[vi]->loc.x<<" "<<vv_real[vi]->loc.x<<" "<<x_error<<endl;
        }
        double error_area = error_sum*sampling_distance;
        return error_area;
    }

    double similarity_Index_1(vector<Vertex> vv_simu,vector<Vertex*> vv_real, double sampling_distance){
        //adjust the real y_position
        /*
        for(int vi=0;vi<(int)vv_real.size();vi++){
            vv_real[vi]->loc.y = abs(1-vv_real[vi]->loc.y);
        }
        vector<Vertex*> vv_real_align;
        for(int i=0;i<101;i++){
            vv_real_align.push_back(vv_real[100-i]);
        }
        for(int i=101;i<200;i++){
            vv_real_align.push_back(vv_real[i-100]);
        }
        //output::vv_foutXY(vv_real_align,"real_aligned.txt");
        */
        double error_sum=0;
        for(int vi=0;vi<(int)vv_simu.size();vi++){
            double x_error = abs(vv_simu[vi].loc.x-vv_real[vi]->loc.x);
            error_sum += x_error;
            //std::cout<<vv_simu[vi]->loc.x<<" "<<vv_real[vi]->loc.x<<" "<<x_error<<endl;
        }
        double error_area = error_sum*sampling_distance;
        return error_area;
    }
    
    double similarity_Index_1(vector<Vertex> vv_simu,vector<Vertex> vv_real, double sampling_distance){
        //adjust the real y_position
        /*
        for(int vi=0;vi<(int)vv_real.size();vi++){
            vv_real[vi]->loc.y = abs(1-vv_real[vi]->loc.y);
        }
        vector<Vertex*> vv_real_align;
        for(int i=0;i<101;i++){
            vv_real_align.push_back(vv_real[100-i]);
        }
        for(int i=101;i<200;i++){
            vv_real_align.push_back(vv_real[i-100]);
        }
        //output::vv_foutXY(vv_real_align,"real_aligned.txt");
        */
        double error_sum=0;
        for(int vi=0;vi<(int)vv_simu.size();vi++){
            double x_error = abs(vv_simu[vi].loc.x-vv_real[vi].loc.x);
            error_sum += x_error;
            //std::cout<<vv_simu[vi]->loc.x<<" "<<vv_real[vi]->loc.x<<" "<<x_error<<endl;
        }
        double error_area = error_sum*sampling_distance;
        return error_area;
    }

    

    double similarity_Index_2(string real_contour_file, string simulated_contour_file, double sampling_distance){
        vector<Vertex*> real_contour = readV::xy_txt_to_vertexXY(real_contour_file);
        vector<Vertex*> real_contour_normalized = geo_vv::after_ImageJ_process(real_contour);
        vector<Vertex*> real_contour_sampled = geo_vv::vector_vertex_sampling(real_contour_normalized,sampling_distance);
        //gnu_plot::organ_contour_plot(real_contour_normalized);

        
        vector<Vertex*> simulated_contour = readV::read_vv(simulated_contour_file);
        vector<Vertex*> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
        vector<Vertex*> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
        vector<Vertex*> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
        vector<Vertex*> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
        //gnu_plot::organ_contour_plot(simulated_contour_normalized);
        cout_fout_debug::fout_vector_vertex(simulated_outline_sampled,"simulated_petal.txt");

        double similarity_index_tmp1 = boundary_geo::similarity_Index_1(simulated_outline_sampled,real_contour_sampled,sampling_distance);
        double similarity_index_tmp2 = boundary_geo::similarity_Index_1(simulated_outline_swapped_sampled,real_contour_sampled,sampling_distance);
        
        double similarity_index;

        if(similarity_index_tmp1>similarity_index_tmp2){
            similarity_index = similarity_index_tmp2;
        }
        else{
            similarity_index = similarity_index_tmp1;
        }
        std::cout<<"similarity_Index1 "<<similarity_index_tmp1<<" similarity_index2 "<<similarity_index_tmp2 <<endl;
        std::cout<<"similarity_index "<<similarity_index<<endl;

        return similarity_index;
    }

    vector<Vertex*> read_and_process_real_organ_contour_imagej(string real_contour_file){
        double sampling_distance=0.01;
        vector<Vertex*> real_contour = readV::xy_txt_to_vertexXY(real_contour_file);
        vector<Vertex*> real_contour_normalized = geo_vv::after_ImageJ_process(real_contour);
        vector<Vertex*> real_contour_sampled = geo_vv::vector_vertex_sampling(real_contour_normalized,sampling_distance);
        //std::cout<<"real_contour_sampled size: "<<real_contour_sampled.size()<<std::endl;
        //cout_fout_debug::cout_vector_vertex(real_contour_sampled);
        return real_contour_sampled;
    }

    double similarity_cal_during_simulation(Organ* p_g,vector<Vertex*> real_contour_sampled){
        //std::cout<<"real_contour_sampled size: "<<real_contour_sampled.size()<<std::endl;
        //cout_fout_debug::cout_vector_vertex(real_contour_sampled);
        double similarity_index;
        vector<Vertex> simulated_contour;
        double sampling_distance=0.01;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            if(p_g->p_v[vi]->IsSurface==1){
                simulated_contour.push_back(*p_g->p_v[vi]);
            }
        }
        vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
        vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
        vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
        vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
        //gnu_plot::organ_contour_plot(simulated_contour_normalized);
        //cout_fout_debug::fout_vector_vertex(simulated_outline_sampled,"simulated_petal.txt");
        //cout_fout_debug::cout_vector_vertex(simulated_outline_swapped_sampled);

        double similarity_index_tmp1 = boundary_geo::similarity_Index_1(simulated_outline_sampled,real_contour_sampled,sampling_distance);
        double similarity_index_tmp2 = boundary_geo::similarity_Index_1(simulated_outline_swapped_sampled,real_contour_sampled,sampling_distance);

        if(similarity_index_tmp1>similarity_index_tmp2){
            similarity_index = similarity_index_tmp2;
        }
        else{
            similarity_index = similarity_index_tmp1;
        }
        //std::cout<<"similarity_Index1 "<<similarity_index_tmp1<<" similarity_index2 "<<similarity_index_tmp2 <<endl;
        //std::cout<<"similarity_index "<<similarity_index<<endl;
        p_g->similarity_index=similarity_index;
        return similarity_index;

    }

    double similarity_for_contour_contour(std::string filename1,std::string filename2, double sampling_distance){
        

        std::vector<Vertex> vv1 = readV::read_vv_(filename1);
        std::vector<Vertex> vv2 = readV::read_vv_(filename2); 

        vector<Vertex> vv1_normalized = geo_vv::normalization(vv1);
        vector<Vertex> vv1_sampled = geo_vv::vector_vertex_sampling(vv1_normalized,sampling_distance);
        vector<Vertex> vv2_normalized = geo_vv::normalization(vv2);
        vector<Vertex> vv2_sampled = geo_vv::vector_vertex_sampling(vv2_normalized,sampling_distance);
        //gnu_plot::organ_contour_plot(simulated_contour_normalized);
        //std::cout<<"***************| vv1_normlized | ******************"<<std::endl;
        //cout_fout_debug::cout_vector_vertex(vv1_normalized);
        //std::cout<<"***************| vv1_sampled | ******************"<<std::endl;
        //cout_fout_debug::cout_vector_vertex(vv1_sampled);
        //std::cout<<"***************| vv2_normlized | ******************"<<std::endl;
        //cout_fout_debug::fout_vector_vertex(vv2, "vv2_contour.txt");
        //cout_fout_debug::fout_vector_vertex(vv2_normalized,"vv2_normalized.txt");
        //cout_fout_debug::fout_vector_vertex(vv1, "vv1_contour.txt");
        //cout_fout_debug::fout_vector_vertex(vv1_normalized,"vv1_normalized.txt");
        //std::cout<<"***************| vv2_sampled | ******************"<<std::endl;
        //cout_fout_debug::fout_vector_vertex(vv2_sampled, "vv2_sampled.txt");
        //cout_fout_debug::fout_vector_vertex(vv1_sampled, "vv1_sampled.txt");

        double similarity_index_tmp = boundary_geo::similarity_Index_1(vv1_sampled,vv2_sampled,sampling_distance);
        //std::cout<<"Difference index: "<<similarity_index_tmp<<std::endl;
        return similarity_index_tmp;
    }
}

namespace circle_geo{
  

    //0. no intersection; 1. one intersection; 2. two intersections
    //not tested
    Intersection_relationship intersection_line_circle(Organ* p_g, int li, Circle cir){
        Intersection_relationship Intersection_Result;
        //the equation for circle: (x-h)^2 + (y-m)^2 = r^2
        //the equation for line: y=mx+b
        //To combine these two equations, we have:
        //(𝑥−ℎ)^2+(𝑚𝑥+𝑏−𝑘)^2=𝑟^2
        //⇒(1+𝑚^2 ) 𝑥^2+(−2ℎ+2𝑚𝑏−2𝑘𝑚)𝑥+(ℎ^2+𝑏^2+𝑘^2−2𝑏𝑘−𝑟^2 )=0
        //This is a simple quadratic equation: 𝐴𝑥^2+𝐵𝑥+𝐶=0
        //The solution is 𝑥_1,2=(−𝐵±√(𝐵^2−4𝐴𝐶))/2𝐴, 𝑖𝑓 √(𝐵^2−4𝐴𝐶)≥0
        p_g->p_l[li]->calc_slope_intercept(*p_g);
        double m_tmp = p_g->p_l[li]->slope;
        double b_tmp = p_g->p_l[li]->intercept; 
        
        double A_tmp = 1+m_tmp*m_tmp;
        double B_tmp = -2*cir.center.x+2*m_tmp*b_tmp-2*cir.center.y*m_tmp;
        double C_tmp = cir.center.x*cir.center.x + b_tmp*b_tmp + cir.center.y*cir.center.y - 2*b_tmp*cir.center.y -cir.radius*cir.radius;
        Quadratic_Solution Result_Quadratic = wangMath::Quadratic_Equation_Solve(A_tmp,B_tmp,C_tmp);
        if(Result_Quadratic.delta>0){
            //there are two intersections 
            Intersection_Result.Relationship = 2;
            Vertex p1_intersection;
            Vertex p2_intersection;
            p1_intersection.loc.x = Result_Quadratic.x1;
            p1_intersection.loc.y = p1_intersection.loc.x*m_tmp+b_tmp;
            p2_intersection.loc.x = Result_Quadratic.x2;
            p2_intersection.loc.y = p2_intersection.loc.x*m_tmp+b_tmp;

            Intersection_Result.cross_points.push_back(p1_intersection);
            Intersection_Result.cross_points.push_back(p2_intersection);
        }
        else if(Result_Quadratic.delta==0){
            //there is only one intersection
            Intersection_Result.Relationship = 1;
            Vertex p1_intersection;
            p1_intersection.loc.x = Result_Quadratic.x1;
            p1_intersection.loc.y = p1_intersection.loc.x*m_tmp+b_tmp;
            Intersection_Result.cross_points.push_back(p1_intersection);
        }
        else{
            //there is no intersection
            Intersection_Result.Relationship = 0;
        }

        return Intersection_Result;
    }
    
    //0. no intersection; 1. one intersection; 2. two intersections
    //not tested
    Intersection_relationship intersection_line_segment_circle(Organ* p_g, int li, Circle cir){
        Intersection_relationship Intersection_Result_line = circle_geo::intersection_line_circle(p_g,li,cir);
        Intersection_relationship Intersection_Result_line_segment;
        Intersection_Result_line_segment.Relationship=0;
        if(Intersection_Result_line.Relationship==2){
            //check if these two intersection lie on the line segment or outside of the line segment
            if(abs(line_geo::distance_line_segment_to_vertex(p_g,li,Intersection_Result_line.cross_points[0]).distance)<EPS_geo){
                //the cross point 0 is indeed inside on the line segment
                Intersection_Result_line_segment.Relationship++;
                Intersection_Result_line_segment.cross_points.push_back(Intersection_Result_line.cross_points[0]);
            }
            if(abs(line_geo::distance_line_segment_to_vertex(p_g,li,Intersection_Result_line.cross_points[1]).distance)<EPS_geo){
                //the cross point 1 is indeed inside on the line segment
                Intersection_Result_line_segment.Relationship++;
                Intersection_Result_line_segment.cross_points.push_back(Intersection_Result_line.cross_points[1]);
            }
        }
        else if(Intersection_Result_line.Relationship==1){
            //check if this intersection lies on the line segment or outside of the line segment
            if(abs(line_geo::distance_line_segment_to_vertex(p_g,li,Intersection_Result_line.cross_points[0]).distance)<EPS_geo){
                //the cross point 0 is indeed inside on the line segment
                Intersection_Result_line_segment.Relationship++;
                Intersection_Result_line_segment.cross_points.push_back(Intersection_Result_line.cross_points[0]);
            }

        }else{
            //no intersection
            Intersection_Result_line_segment.Relationship=0;
        }
        return Intersection_Result_line_segment;
    }
    /*
    //there are seven different relationships for a line segment and a circle:
    //1. cut twice, where the segment and the circle has two distinct intersections
    //2. cut once, where the segment and the circle has one intersection
    //3. tangent, where the segment and the circle only have one intersection
    //4. miss outside, where the segment is outside of the circle and they have no intersection
    //5. miss inside, where the segment is inside of the circle and they have no intersection
    //6. end cut once (outside), where the segment and the circle has one intersection, the intersection is at the end of the line segment, the other part of line segment is outside of the cirlce 
    //7. end cut once (inside), where the segment and the circle has one intersection, the intersection is at the end of the line segment, the other part of line segment is inside of the circle
    //please notice that 3,6 actually has the overlap situation, where a tangent cross point happens to be the end cut, in that case, it is defined as case 6
    //the function returns a class which has an integer index to identify the case, an a vector<Vertex> to store the cross points

    Intersection_relationship intersection_line_segment_circle(Organ* p_g, int li, Circle cir)
    {
        Intersection_relationship Result;
        Vertex d1 = *p_g->p_v[p_g->p_l[li]->vi[0]];
        Vertex d2 = *p_g->p_v[p_g->p_l[li]->vi[1]];
        
        //calculate the distance between the circle center and the line segment
        Vertex cir_center;
        cir_center.loc = cir.center;
        Distance_point_line dist = line_geo::distance_line_segment_to_vertex(p_g,li,cir_center);
        if(abs(dist.distance-cir.radius)<EPS_geo){
            //distance same with radius: 3.tangent or 6. End cut once (outside)
            //the cross point is the closest point
            Result.cross_points.push_back(dist.Closest_Point);
            if(dist.t<0||dist.t>1){
                Result.Relationship=6;
            }
            else{
                Result.Relationship=3;
            }
        }
        else if(dist.distance>cir.radius)
        {
            //distance longer than radius: 4. miss outside
            Result.Relationship=4;
        }
        else{
            //distance shorter than radius 1.cut twice, 2. cut once, 5. miss inside, 7 end cut once (inside)
            //The distance between endpoints and circle center:
            double dist1 = cir_center.distance_from_vertex(d1);
            double dist2 = cir_center.distance_from_vertex(d2);

            //d1,d2>radius  => case 1
            if(dist1>cir.radius&&dist2>cir.radius)
            {
                Result.Relationship=1;

            }
            //d1>radius>d2  => case 2
            if((dist1>cir.radius&&dist2<cir.radius)||(dist1<cir.radius&&dist2>cir.radius))
            {
                Result.Relationship=2;
            }
            //d1,d2<radius  => case 5
            if(dist1<cir.radius&&dist2<cir.radius){
                Result.Relationship=5;
            }
            //d1=radius, d2<radius => case 7
            if(dist.t<0||dist.t>1){
                Result.Relationship=7;
            }

        }

        return Result;
    }
    */

    Intersection_relationship intersection_line_segment_circle(Line ls ,Circle cir){
        Intersection_relationship Intersection_Result_line;
        Intersection_relationship Intersection_Result_segment;
        ls.calc_slope_intercept();
        double m_t = ls.slope;
        double b_t = ls.intercept;
        double A_t = 1+m_t*m_t;
        double B_t = -2*cir.center.x+2*m_t*b_t-2*cir.center.y*m_t;
        double C_t = cir.center.x*cir.center.x + b_t*b_t + cir.center.y*cir.center.y - 2*b_t*cir.center.y -cir.radius*cir.radius;
        Quadratic_Solution Result_Quadratic = wangMath::Quadratic_Equation_Solve(A_t,B_t,C_t);
        if(Result_Quadratic.delta>0){
            //there are two intersections 
            Intersection_Result_line.Relationship = 2;
            Vertex p1_intersection;
            Vertex p2_intersection;
            p1_intersection.loc.x = Result_Quadratic.x1;
            p1_intersection.loc.y = p1_intersection.loc.x*m_t+b_t;
            p2_intersection.loc.x = Result_Quadratic.x2;
            p2_intersection.loc.y = p2_intersection.loc.x*m_t+b_t;

            Intersection_Result_line.cross_points.push_back(p1_intersection);
            Intersection_Result_line.cross_points.push_back(p2_intersection);
        }
        else if(Result_Quadratic.delta==0){
            //there is only one intersection
            Intersection_Result_line.Relationship = 1;
            Vertex p1_intersection;
            p1_intersection.loc.x = Result_Quadratic.x1;
            p1_intersection.loc.y = p1_intersection.loc.x*m_t+b_t;
            Intersection_Result_line.cross_points.push_back(p1_intersection);
        }
        else{
            //there is no intersection
            Intersection_Result_line.Relationship = 0;
        }

        //check if the intersections between the line and circle, lies inside the line segment, if outside they will be ignored
        Intersection_Result_segment.Relationship=0;

        if(Intersection_Result_line.Relationship==2){
            if(ls.distance_from_point(Intersection_Result_line.cross_points[0])<EPS_geo){
                //the crosspoint 0 is indeed inside the line segment
                Intersection_Result_segment.Relationship++;
                Intersection_Result_segment.cross_points.push_back(Intersection_Result_line.cross_points[0]);
            }
            if(ls.distance_from_point(Intersection_Result_line.cross_points[1])<EPS_geo){
                //the crosspoint 1 is indeed inside the line segment
                Intersection_Result_segment.Relationship++;
                Intersection_Result_segment.cross_points.push_back(Intersection_Result_line.cross_points[1]);
            }
        }
        else if(Intersection_Result_line.Relationship==1){
            if(ls.distance_from_point(Intersection_Result_line.cross_points[0])<EPS_geo){
                //the crosspoint 0 is indeed inside the line segment
                Intersection_Result_segment.Relationship++;
                Intersection_Result_segment.cross_points.push_back(Intersection_Result_line.cross_points[0]);
            }
        }
        else{
            //no intersection
        }

        return Intersection_Result_segment;
    }

}


namespace geo_analysis{
    void geo_plot(Geometrics_analysis_class ga, std::string x_axis, std::string y_axis, std::string save_plot_txt){
        std::vector<double> x_axis_value;
        std::vector<double> y_axis_value;
        
        if(x_axis=="cell_number"){
            std::vector<double> inner_cell_number = ga.getColumnData("inner_cell_number");
            std::vector<double> epi_cell_number = ga.getColumnData("epidermal_cell_number");

            for(int i=0;i<inner_cell_number.size();i++){
                double value_i = inner_cell_number[i]+epi_cell_number[i];
                x_axis_value.push_back(value_i);
            }
        }
        else{
            x_axis_value = ga.getColumnData(x_axis);
        }
        
        y_axis_value = ga.getColumnData(y_axis);

        cout_fout_debug::fout_vector_x_y(x_axis_value, y_axis_value, save_plot_txt);
    }  
}