/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#include "../include/IO2dv.h"
#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"
#include "../include/gnuplot.h"

namespace readV{

    //read_organ_txt 
    //reads organ information (vertex position, line connectivity and cell connectivity) from the designated "initialOrganFile".txt file
    //the "initialOrganFile" requires a specified format: 
    //description_of_file
    //number_of_vertex
    //vertex_index, x,y,z position of vertex vi
    //i x_i y_i z_i
    //number_of_line
    //line index, the two vertex indices connected by line lj
    //j l(j,1) l(j,2)
    //number of hexagon
    //hegaxon index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4) cl(k,5) cl(k,6)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4) cv(k,5) cv(k,6)
    //number of pentagon
    //pentagon index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4) cl(k,5)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4) cv(k,5)
    //number of quandrangle
    //quandrangle index, the six line indices and vertex indices (vertex indices must be in anticlockwise direction)
    //k cl(k,1) cl(k,2) cl(k,3) cl(k,4)
    //k cv(k,1) cv(k,2) cv(k,3) cv(k,4)
    void read_organ_txt(Organ* p_g,int fixed_index){
        cout<<"*************|Initial Shape (primordia) Information ("<<initialOrganFile<<") |***************"<<endl;
        //extern string initialOrganFile;
        //cout<<"initialOrganFile: "<<initialOrganFile<<endl;
        ifstream fin(initialOrganFile,ios::in);
        if(!fin.is_open()){
            initialOrganFile = "../"+initialOrganFile;
            ifstream fin(initialOrganFile,ios::in);
        }
        if(!fin.is_open()){
            cout<<"Error: missing "<<initialOrganFile <<" (the initial state file!)"<<endl;        
            exit(1);
        }
        string description;
        int vertex_number,line_number,hexagon_number,pentagon_number,quadrangle_number;
        fin>>description;
        fin>>vertex_number;

        std::cout<<"Description: "<<description<<std::endl;
        std::cout<<"Number of vertices in initial organ: "<<vertex_number<<std::endl;

        for(int i=0; i<vertex_number;i++)
        {   
            Vertex *tmp = new Vertex;
            int vertex_index;
            //fin>>vertex_index;
            fin>>tmp->loc.x;
            fin>>tmp->loc.y;
            fin>>tmp->loc.z;
            tmp->loc.z =0.0;
            tmp->vi=i;

            p_g->p_v.push_back(tmp);
        }


        //Debug: use initialOrgan to check the force amd motion calculation
        /*
        std::cout<<"Debug: use initialOrgan to check the force amd motion calculation"<<std::endl;
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
        p_g->p_v[vi]->loc[0].x = 1.1*p_g->p_v[vi]->loc[0].x;
        p_g->p_v[vi]->loc[0].y = 1.1*p_g->p_v[vi]->loc[0].y;
        }
        */
        //debug: check vertex input
        /*
        std::cout<<"Debug: check vertex input"<<std::endl;
        for(int i=0; i<vertex_number;i++)
        {
            std::cout<<i<<" "<<p_g->p_v[i]->loc[0].x<<" "<<p_g->p_v[i]->loc[0].y<<std::endl;
        }
        std::cout<<"End Debug: check vertex input" <<std::endl;
        */

        fin>>line_number;
        if(fixed_index==0)
    std::cout<<"Number of lines in initial state: "<<line_number<<std::endl;

    for(int i=0; i<line_number;i++)
    {
        int line_index;
        Line *tmp = new Line;
        
        fin>>line_index;
        fin>>tmp->vi[0];
        fin>>tmp->vi[1];
        tmp->li=i;

        p_g->p_l.push_back(tmp);
    }

        //debug: check line input
        /*
        std::cout<<"Debug: check line input"<<std::endl;
        for(int i=0; i<line_number;i++){
            std::cout<<i<<" "<<p_g->p_l[i]->vi[0]<<" "<<p_g->p_l[i]->vi[1]<<std::endl;
        }
        std::cout<<"End Debug: check line input"<<std::endl;
        */

        for(int i=0; i<line_number;i++){
        p_g->p_v[p_g->p_l[i]->vi[0]]->li.push_back(i);
        p_g->p_v[p_g->p_l[i]->vi[1]]->li.push_back(i);
    }

    //debug: check which lines vertex belongs to
    /*
    std::cout<<"Debug: check which lines vertex belongs to."<<std::endl;
    for(int i=0; i<vertex_number; i++){
        std::cout<<"For vertex "<<i<<", it belongs to line ";
        for(int j=0; j<p_g->p_v[i]->li.size();j++){
            std::cout<<p_g->p_v[i]->li[j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<"End Debug: check which lines vertex belongs to."<<std::endl;
        */


        fin>>hexagon_number;
        
        
        for(int i=0; i<hexagon_number; i++){
            Cell *tmp = new Cell;
            int hexagon_index,line_index,vertex_index;
            fin>>hexagon_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);

            fin>>hexagon_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            
            tmp->n_edges=6;
            p_g->p_c.push_back(tmp);       
        }

        fin>>pentagon_number;
        
        for(int i=0; i<pentagon_number; i++){
            Cell *tmp = new Cell;
            int pentagon_index, line_index, vertex_index;
            fin>>pentagon_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);

            fin>>pentagon_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);

            tmp->n_edges=5;
            p_g->p_c.push_back(tmp);
        }
        fin>>quadrangle_number;
        int cell_number = hexagon_number+pentagon_number+quadrangle_number;
        
        if(fixed_index==0){
            std::cout<<"Number of hexagons in initial state : "<<hexagon_number<<std::endl;
            std::cout<<"Number of pentagons in initial state : "<<pentagon_number<<std::endl;
            std::cout<<"Number of quadrangles in initial state : "<<quadrangle_number<<std::endl;
            std::cout<<"Number of all cells in initial state : "<<cell_number<<std::endl;
        }
        for(int i=0; i<quadrangle_number; i++){
            Cell *tmp = new Cell;
            int quadrangle_index,  line_index, vertex_index;
            fin>>quadrangle_index;
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            fin>>line_index;
            tmp->li.push_back(line_index);
            

            fin>>quadrangle_index;
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            fin>>vertex_index;
            tmp->vi.push_back(vertex_index);
            
            tmp->n_edges=4;
            p_g->p_c.push_back(tmp);
        }

        //Debug: check the vertex index and line index in each cell
        /*
        std::cout<<"Debug: check the vertex index and line index in each cell"<<std::endl;
        for(int i=0; i<p_g->p_c.size(); i++){
            std::cout<<"Cell "<<i<<" has line ";
            for(int j=0; j<p_g->p_c[i]->n_edges; j++){
                std::cout<<p_g->p_c[i]->li[j]<<" ";
            }
            std::cout<<std::endl;
            
            std::cout<<"Cell "<<i<<" has vertex ";
            for(int j=0; j<p_g->p_c[i]->n_edges; j++){
                std::cout<<p_g->p_c[i]->vi[j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"End Debug: check the vertex index and line index in each cell"<<std::endl;
        */

        for(int i=0; i<cell_number; i++){
            for(int j=0; j<p_g->p_c[i]->n_edges;j++){
                p_g->p_l[p_g->p_c[i]->li[j]]->ci.push_back(i);
                p_g->p_v[p_g->p_c[i]->vi[j]]->ci.push_back(i);
            }
        }

        //initialize the cellDivisionCount
        for(int ci=0; ci<(int)p_g->p_c.size();ci++){
            p_g->p_c[ci]->cellDivisionCount=0;
        }

        //debug: check which cells vertex belongs to
        /*
        std::cout<<"Debug: check which cells each vertex belongs to."<<std::endl;
        for(int i=0; i<vertex_number; i++){
        std::cout<<"For vertex "<<i<<", it belongs to cell ";
        for(int j=0; j<p_g->p_v[i]->ci.size();j++){
            std::cout<<p_g->p_v[i]->ci[j]<<" ";
        }
        std::cout<<std::endl;
        }
        std::cout<<"End Debug: check which cells each vertex belongs to."<<std::endl;
        */

        //debug: check which cells line belongs to
        /*
        std::cout<<"Debug: check which cells each line belongs to."<<std::endl;
        for(int i=0; i<line_number; i++){
        std::cout<<"For line "<<i<<", it belongs to cell ";
        for(int j=0; j<p_g->p_l[i]->ci.size();j++){
            std::cout<<p_g->p_l[i]->ci[j]<<" ";
        }
        std::cout<<std::endl;
        }
        std::cout<<"End Debug: check which cells each line belongs to."<<std::endl;
        */
       p_g->initial_cell_number = p_g->p_c.size();
       fin.close();
       cout<<"********************************************************************"<<endl;
    }

//read vertex position and line-vertex, cell-vertex connectivity from cell.vtk and line.vtk
    void oneCell(Organ* p_g,string FileName){
        ifstream fin(FileName,ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing cell.vtk ("<<FileName<<")"<<endl;
            exit(1);
        }
        char parameter_name[100];
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        fin>>parameter_name;
        int vertex_number;
        fin>>vertex_number;
        fin>>parameter_name;

        //read vertex position information
        for(int i=0; i<vertex_number;i++){
            Vertex *v_tmp = new Vertex;
            fin>>v_tmp->loc.x;
            fin>>v_tmp->loc.y;
            fin>>v_tmp->loc.z;
            p_g->p_v.push_back(v_tmp);
        }
        //debug: vertex position information
        /*
        for(int i=0; i<vertex_number;i++){
            std::cout<<"Vertex "<< i<<" position: ("<<p_g->p_v[i]->loc.x<<","<<p_g->p_v[i]->loc.y<<","<<p_g->p_v[i]->loc.z<<")."<<std::endl;
        }
        */
    
        //read cell-vertex connection informtation
        fin>>parameter_name;
        int cell_number;
        fin>>cell_number;
        fin>>parameter_name;
        for(int ci=0;ci<cell_number;ci++){
            Cell *c_tmp = new Cell;
            int cell_vertex_number;
            fin>>cell_vertex_number;
            for(int vi=0;vi<cell_vertex_number;vi++){
                int vertex_index;
                fin>>vertex_index;
                c_tmp->vi.push_back(vertex_index);
            }
            p_g->p_c.push_back(c_tmp);
        }
        if(minor_mode=="continue"&&major_mode=="simulation"){
            fin>>parameter_name;
            int index_tmp;
            fin>>index_tmp;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                fin>>index_tmp;
            }
            fin>>parameter_name;
            fin>>index_tmp;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                fin>>p_g->p_c[ci]->cellDivisionCount;
            }
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                fin>>p_g->p_c[ci]->cellTime;
            }
            //output::VTK_debug(p_g,"test");
            //exit(1);
        }


        fin.close();
    }

    void oneLine(Organ* p_g,string FileName){
        ifstream fin(FileName,ios::in);
            if(!fin.is_open()){
            cout<<"Error: missing Line.vtk ("<<FileName<<")"<<endl;
                return;
            }
            char parameter_name[100];
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            fin>>parameter_name;
            int vertex_number;
            fin>>vertex_number;
            fin>>parameter_name;

            //read vertex position information
            for(int i=0; i<vertex_number;i++){
                Vertex *v_tmp = new Vertex;
                fin>>v_tmp->loc.x;
                fin>>v_tmp->loc.y;
                fin>>v_tmp->loc.z;
                //p_g->p_v.push_back(v_tmp);
            }

            //read line-vertex connectivity information
            fin>>parameter_name;
            int line_number;
            fin>>line_number;
            fin>>parameter_name;
            for(int li=0;li<line_number;li++){
                int line_index;
                fin>>line_index;
                Line *l_tmp = new Line;
                fin>>l_tmp->vi[0];
                fin>>l_tmp->vi[1];
                p_g->p_l.push_back(l_tmp);
            }

            //debug: line-vertex connectivity information
            /*
            for(int li=0;li<(int)p_g->p_l.size();li++){
                std::cout<<li<<" "<<p_g->p_l[li]->vi[0]<<" "<<p_g->p_l[li]->vi[1]<<std::endl;
            }
            */
            fin.close();
    }

    //line-vertex, cell-vertex connectivity => vertex-cell, vertex-line, cell-line, line-cell connectivity
    void oneCellLine(Organ* p_g){
        //add line index to cell
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            for(int vi=0;vi<(int)p_g->p_c[ci]->vi.size();vi++){
                for(int li=0; li<(int)p_g->p_l.size();li++){
                    if(p_g->p_l[li]->vi[0]==p_g->p_c[ci]->vi[vi]){
                        for(int vii=0;vii<(int)p_g->p_c[ci]->vi.size();vii++){
                            if(p_g->p_l[li]->vi[1]==p_g->p_c[ci]->vi[vii]){
                                p_g->p_c[ci]->li.push_back(li);
                            }
                        }
                    }
                }
            }
        }

        //debug: cell-line connnectivity information
        /*
        for(int ci=0;ci<(int)p_g->p_c.size();ci++){
            std::cout<<ci;
            for(int li=0;li<(int)p_g->p_c[ci]->li.size();li++){
                std::cout<<" "<<p_g->p_c[ci]->li[li];
            }
            std::cout<<std::endl;
        }
        */

        //add cell index to line and vertex
        for(int i=0; i<(int)p_g->p_c.size(); i++){
            for(int j=0; j<p_g->p_c[i]->vi.size();j++){
                p_g->p_l[p_g->p_c[i]->li[j]]->ci.push_back(i);
                p_g->p_v[p_g->p_c[i]->vi[j]]->ci.push_back(i);
            }
        }

        //add line index to vertex 
        for(int i=0; i<(int)p_g->p_l.size();i++){
            //a line can only connect two vertices
            p_g->p_v[p_g->p_l[i]->vi[0]]->li.push_back(i);
            p_g->p_v[p_g->p_l[i]->vi[1]]->li.push_back(i);
        }

        for(int li=0;li<(int)p_g->p_l.size();li++){
            p_g->p_l[li]->li = li;
        }
        for(int vi=0;vi<(int)p_g->p_v.size();vi++){
            p_g->p_v[vi]->vi = vi;
        }
    }

    //all information of an Organ in one time step
    void oneVTK(Organ* p_g,string CellVTK,string LineVTK, int simulation_steps){
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
            //cout_fout_debug::cout_vector_vertex(real_organ_contour_processed_for_similarity_index);
        }
        oneCell(p_g,CellVTK);
        oneLine(p_g,LineVTK);
        oneCellLine(p_g);
        p_g->step=simulation_steps;
    }

    void oneVTK(Organ* p_g,string CellVTK,string LineVTK){
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
            //cout_fout_debug::cout_vector_vertex(real_organ_contour_processed_for_similarity_index);
        }
        oneCell(p_g,CellVTK);
        oneLine(p_g,LineVTK);
        oneCellLine(p_g);
    }

    void oneVTK(Organ* p_g,std::pair<std::string,std::string> VTK_dir){
        if(similarity_calculation_required==1){
            real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);
            //cout_fout_debug::cout_vector_vertex(real_organ_contour_processed_for_similarity_index);
        }
        oneCell(p_g,VTK_dir.first);
        oneLine(p_g,VTK_dir.second);
        oneCellLine(p_g);
    }

    //read the final VTK inside a filefolder
    void final_VTK(Organ* p_g,string VTK_file_path){
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
        //cout<<"file_index "<<file_index<<endl;
        char str_cell[100], str_line[100];
        sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-5);
        sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-5);
        if(minor_mode=="continue"&&major_mode=="simulation"){
            p_g->step = file_index/2-5;
            std::cout<<"p_g->step: "<<p_g->step<<std::endl;
        }
        cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
        string CellVTK = VTK_file_path+"/"+str_cell;
        string LineVTK = VTK_file_path+"/"+str_line;

        readV::oneVTK(p_g, CellVTK, LineVTK);
    }

    //read the outline information from imageJ extracted xy txt file
    vector<Vertex*> xy_txt_to_vertexXY(string FileName){
        vector<Vertex*> vv;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        while(!fin.eof()){
            Vertex * v_tmp = new Vertex;
            int tmp;
            fin>>tmp;
            fin>>v_tmp->loc.x;
            fin>>v_tmp->loc.y;
            vv.push_back(v_tmp);
        }
        return vv;
    }

    vector<double> read_vd(string FileName){
        vector<double> vd;
        string nothing;
        
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        fin>>nothing;
        fin>>nothing;
        while(!fin.eof()){
            double tmp;
            fin>>tmp;
            fin>>tmp;
            vd.push_back(tmp);
        }
        return vd;
    }

    vector<Vertex*> read_vv(string FileName){
        vector<Vertex*> vv;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        string s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        while(!fin.eof()){
            double tmp;
            Vertex* v_tmp = new Vertex;
            fin>>tmp;
            fin>>tmp;
            v_tmp->loc.x = tmp;
            fin>>tmp;
            v_tmp->loc.y = tmp;
            fin>>tmp;
            v_tmp->loc.z = tmp;
            vv.push_back(v_tmp);
        }
        if(vv[vv.size()-1]->loc.x==0&&vv[vv.size()-1]->loc.y==0)
            vv.pop_back();

        return vv;
    }

    vector<Vertex> read_vv_(string FileName){
        vector<Vertex> vv;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        string s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        fin>>s_tmp;
        while(!fin.eof()){
            double tmp;
            Vertex v_tmp;
            fin>>tmp;
            fin>>tmp;
            v_tmp.loc.x = tmp;
            fin>>tmp;
            v_tmp.loc.y = tmp;
            fin>>tmp;
            v_tmp.loc.z = tmp;
            vv.push_back(v_tmp);
        }
        return vv;
    }

    Geometrics_analysis_class read_geo_output_txt(string FileName){
        Geometrics_analysis_class ga;
        ifstream fin(FileName, ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing "<<FileName<<endl;
            exit(1);
        }
        
        // Reading first line of variable names
        std::string line;
        std::getline(fin,line);
        std::stringstream ss(line);
        std::vector<std::string> headers;
        std::string header;

        while( ss >>header){
            headers.push_back(header);
        }
        ga.variable_name=headers;
        // Reading the following lines of variable values
        std::vector<std::vector<double>> all_values;

        while(std::getline(fin,line)){
            ss.clear();
            ss.str(line);
            std::vector<double> values;
            double value;

            while(ss>>value){
                values.push_back(value);
            }

            all_values.push_back(values);
        }

        ga.value = all_values;

        fin.close();
        return ga;
    }

}

namespace output{

void VTK(Organ* p_g){
    {
        char fname[100];
        sprintf(fname,"2dv_line%010u.vtk", p_g->step);
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 3.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
        fout<<std::endl;
        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        fout<<std::endl;
        //線の要素を書き込む
        fout << "CELLS " << p_g->p_l.size() << " " << p_g->p_l.size() * 3 << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        Line *lp = p_g->p_l[i];
        fout << "2 " << lp->vi[0] << " " << lp->vi[1] << std::endl;
        }
        fout<<std::endl;
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_l.size() << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        fout << "3" << std::endl;
        }
        fout<<std::endl;
        fout << "CELL_DATA " <<  p_g->p_l.size() << std::endl;
        fout << "SCALARS current_Length float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->length<<std::endl;
        }
        fout<<std::endl;
        fout << "SCALARS edge_force float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->edgeForce<<std::endl;
        }
        fout.close();
    }

    {
        char fname[100];
        sprintf(fname,"2dv_cell%010u.vtk", p_g->step);
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 2.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        int cells_size = 0;

        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        cells_size += cp->vi.size() + 1;
        }

        fout << "CELLS " << p_g->p_c.size() << " " << cells_size << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        fout << cp->vi.size();
        for (int j = 0; j < (int)cp->vi.size(); j++) {
            fout << " " << cp->vi[j];
        }
        fout << std::endl;
        }
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_c.size() << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        fout << "7" << std::endl;
        }
        
        //Field 
        fout << "CELL_DATA " <<  p_g->p_c.size() << std::endl;
        //output the cell division count
        fout << "SCALARS cell_division_count float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellDivisionCount << std::endl;
        }
        fout << "SCALARS cell_time float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellTime << std::endl;
        }
        
        fout << "SCALARS cell_area float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->area << std::endl;
        }
        fout << "SCALARS cell_layer float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->layer << std::endl;
        }
        fout << "SCALARS cell_perimeter float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->perimeter << std::endl;
        }
        fout << "SCALARS cell_frequency_modifier float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->frequency_modifier<< std::endl;
        }
        if(division_control=="Gaussian_area"||division_control=="area"){
            fout << "SCALARS cell_area_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->area_modifier<< std::endl;
            }
        }
        if(division_control=="Gaussian_area"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="uniform_to_Gaussian"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="temporal_Gaussian_linear_mu"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="temporal_Gaussian_linear_sigma"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="biregion_frequency_position"||division_control=="biregion_frequency_identity"||in_division_direction=="biregion_angles_identity"){
            fout << "SCALARS apical_basal_tag float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->tag<< std::endl;
            }
        }
        fout.close();
    }

}
/*
void VTK_surface(Organ* p_g){
    {
        char fname[100];
        sprintf(fname,"2dv_surface_line%010u.vtk", p_g->step);
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 3.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
        fout<<std::endl;
        //点の位置を書き込む
        fout << "POINTS " << .size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        fout<<std::endl;
        //線の要素を書き込む
        fout << "CELLS " << p_g->p_l.size() << " " << p_g->p_l.size() * 3 << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        Line *lp = p_g->p_l[i];
        fout << "2 " << lp->vi[0] << " " << lp->vi[1] << std::endl;
        }
        fout<<std::endl;
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_l.size() << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        fout << "3" << std::endl;
        }
        fout<<std::endl;

        fout << "CELL_DATA " <<  p_g->p_l.size() << std::endl;
        fout << "SCALARS current_Length float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->length<<std::endl;
        }
        fout<<std::endl;
        fout << "SCALARS edge_force float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->edgeForce<<std::endl;
        }
        fout.close();
    }
}
*/

void VTK_debug(Organ* p_g,std::string output_fileName){
    {
        char fname[100];
        sprintf(fname,"2dv_line%010u.vtk", p_g->step);
        strcat(fname, output_fileName.c_str());
        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 3.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
        fout<<std::endl;
        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        fout<<std::endl;
        //線の要素を書き込む
        fout << "CELLS " << p_g->p_l.size() << " " << p_g->p_l.size() * 3 << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        Line *lp = p_g->p_l[i];
        fout << "2 " << lp->vi[0] << " " << lp->vi[1] << std::endl;
        }
        fout<<std::endl;
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_l.size() << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); i++) {
        fout << "3" << std::endl;
        }
        fout<<std::endl;

        fout << "CELL_DATA " <<  p_g->p_l.size() << std::endl;
        fout << "SCALARS current_Length float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->length<<std::endl;
        }
        fout<<std::endl;
        fout << "SCALARS edge_force float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
            fout <<p_g->p_l[i]->edgeForce<<std::endl;
        }
        fout.close();
    }

    {
        char fname[100];
        sprintf(fname,"2dv_cell%010u.vtk", p_g->step);
        strcat(fname, output_fileName.c_str());

        std::ofstream fout(fname);

        //vtkファイルのヘッダ
        fout << "# vtk DataFile Version 2.0" << std::endl;
        fout << "2D-vertex" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

        //点の位置を書き込む
        fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
        for (int i = 0; i < (int)p_g->p_v.size(); i++) {
        Vertex *vp = p_g->p_v[i];
        fout << vp->loc.x << " ";
        fout << vp->loc.y << " ";
        fout << vp->loc.z << std::endl;
        }
        int cells_size = 0;

        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        cells_size += cp->vi.size() + 1;
        }

        fout << "CELLS " << p_g->p_c.size() << " " << cells_size << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        Cell *cp = p_g->p_c[i];
        fout << cp->vi.size();
        for (int j = 0; j < (int)cp->vi.size(); j++) {
            fout << " " << cp->vi[j];
        }
        fout << std::endl;
        }
        //CELL_TYPES
        fout << "CELL_TYPES " << p_g->p_c.size() << std::endl;
        for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
        fout << "7" << std::endl;
        }
        
        //Field 
        fout << "CELL_DATA " <<  p_g->p_c.size() << std::endl;
        //output the cell division count
        fout << "SCALARS cell_division_count float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellDivisionCount << std::endl;
        }
        fout << "SCALARS cell_time float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->cellTime << std::endl;
        }
        
        fout << "SCALARS cell_area float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->area << std::endl;
        }
        fout << "SCALARS cell_layer float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->layer << std::endl;
        }
        fout << "SCALARS cell_perimeter float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->perimeter << std::endl;
        }
        fout << "SCALARS cell_frequency_modifier float" << std::endl;
        fout << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < (int)p_g->p_c.size(); i++) {
            fout << p_g->p_c[i]->frequency_modifier<< std::endl;
        }
        if(division_control=="Gaussian_area"||division_control=="area"){
            fout << "SCALARS cell_area_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->area_modifier<< std::endl;
            }
        }
        if(division_control=="Gaussian_area"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="uniform_to_Gaussian"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="temporal_Gaussian_linear_mu"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="temporal_Gaussian_linear_sigma"){
            fout << "SCALARS cell_Gaussian_modifier float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->Gaussian_modifier<< std::endl;
            }
        }
        if(division_control=="biregion_frequency_position"||division_control=="biregion_frequency_identity"||in_division_direction=="biregion_angles_identity"){
            fout << "SCALARS apical_basal_tag float" << std::endl;
            fout << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < (int)p_g->p_c.size(); i++) {
                fout << p_g->p_c[i]->tag<< std::endl;
            }
        }
        fout.close();
    }

}


void geo_initial(void){
    string fname = "geometrics_record.txt";
    std::ofstream fout(fname);
    fout<<"Step inner_cell_number epidermal_cell_number center.x center.y cell_layer_number area epiArea inArea epiArea_averaged inArea_averaged area_averaged perimeter perimeter_averaged circularity regularity_av regularity_in_av regularity_epi_av organ_length organ_width length_width_ratio cell_arrangement_x cell_arrangement_y cell_arrangmenet_ratio elliptical_indx overlap_index overlapArea realArea similarity_index maximum_curvature minimum_curvature accumulated_negative_curvature"<<endl;
    fout.close();
}

void geo(Organ* p_g){
    string fname = "geometrics_record.txt";
    ofstream outfile(fname, std::ios::app);
    if(!outfile){
        std::cerr<<"Failed to open "<<fname<<std::endl;
        exit(1);
    }
    outfile<<p_g->step<<" "<<p_g->N_inner_cell<<" "<<p_g->N_epi_cell<<" ";
    outfile<<p_g->center.x<<" "<<p_g->center.y<<" ";
    outfile<<p_g->cell_layer_number<<" ";
    outfile<<p_g->area<<" "<<p_g->epiArea<<" "<<p_g->inArea<<" "<<p_g->epiArea_averaged<<" "<<p_g->inArea_averaged<<" "<<p_g->area_averaged<<" ";
    outfile<<p_g->perimeter<<" "<<p_g->perimeter_averaged<<" ";
    outfile<<p_g->circularity<<" ";
    outfile<<p_g->regularity_averaged<<" "<<p_g->regularity_in_av<<" "<<p_g->regularity_epi_av<<" ";
    outfile<<p_g->organ_length<<" "<<p_g->organ_width<<" "<<p_g->length_width_ratio<<" ";
    outfile<<p_g->cell_arrangement_x<<" "<<p_g->cell_arrangement_y<<" "<<p_g->cell_arrangement_ratio<<" ";
    outfile<<p_g->elliptical_index<<" ";
    outfile<<p_g->overlap_index<<" "<<p_g->realArea<<" "<<p_g->overlapArea<<" ";
    outfile<<p_g->similarity_index<<" ";
    outfile<<p_g->maximum_curvature<<" "<<p_g->minimum_curvature<<" "<<p_g->accumulated_negative_curvature<<endl;

    outfile.close();
}

void division(Organ *p_g){
    //char fname[100] ="divisionRecord.txt";
    char fname[100];
    sprintf(fname,"divisionRecord.txt");
    std::ofstream fout(fname);
    fout<<"division_index "<<"time_step "<<"divided_cell_index "<<"divided_cell_IsEpidermal "<<"division_axisTheta "<<"divided_cell_x "<<"divided_cell_y "<<"cell_division_count "<<"av_in_frequency_modifier "<<"av_epi_frequency_modifier "<<std::endl;
    for(int di=0;di<(int)p_g->d_r.size();di++){
        fout<<di<<" "<<p_g->d_r[di]->time<<" "<<p_g->d_r[di]->cidx<<" "<<p_g->d_r[di]->IsEpidermal<<" "<<p_g->d_r[di]->axisTheta<<" "<<p_g->d_r[di]->center_x<<" "<<p_g->d_r[di]->center_y<<" "<<p_g->d_r[di]->division_count<<" "<<p_g->d_r[di]->av_in_frequency_modifier<<" "<<p_g->d_r[di]->av_epi_frequency_modifier<<std::endl;
    }
    fout.close();
}

double averaged_vec_double(vector<double> vec){
    double sum = 0;
    for(int i=0;i<vec.size();i++){
        sum+=vec[i];
    }
    return sum/vec.size();
};

double averaged_vec_int(vector<int> vec){
    double sum = 0;
    for(int i=0;i<vec.size();i++){
        sum+=vec[i];
    }
    return sum/vec.size();
};

void batch(Batch *ba){
    char fname[100];
    sprintf(fname,"batch_information.txt");
    std::ofstream fout(fname);

    fout<<"simulation_index sigma_O sigma_L sigma_O/sigma_L N_in N_epi organ_area organ_perimeter circularity regularity_in regularity_epi overlap_index"<<endl;
    
    for(int i=0;i<(int)ba->N_in.size();i++){
        fout<<i<<" "<<sigma_O<<" "<<sigma_L<<" "<<sigma_O/sigma_L<<" "<<ba->N_in[i]<<" "<<ba->N_epi[i]<<" "<<ba->organ_area[i]<<" "<<ba->organ_perimeter[i]<<" "<<4*3.1415926*ba->organ_area[i]/(ba->organ_perimeter[i]*ba->organ_perimeter[i])<<" "<<ba->averaged_regularity_in[i]<<" "<<ba->averaged_regularity_epi[i]<<" "<<ba->overlap_index[i]<<endl;
    }

    /*
    fout<<"simulation_index N_in N_epi N_all organ_area organ_perimeter averaged_in_area averaged_epi_area averaged_perimeter";
    fout<<"regularity_in_av regularity_epi_av real_area overlap_area overlap_index"<<std::endl;
    for(int i=0;i<(int)ba->N_in.size();i++){
        fout<<i<<" "<<ba->N_in[i]<<" "<<ba->N_epi[i]<<" "<<ba->N_in[i]+ba->N_epi[i]<<" "<<ba->organ_area[i]<<" "<<ba->organ_perimeter[i]<<" "<<ba->averaged_inner_area[i]<<" "<<ba->averaged_epi_area[i]<<" "<<ba->averaged_perimeter[i];
        fout<<" "<<ba->averaged_regularity_in[i]<<" "<<ba->averaged_regularity_epi[i]<<" "<<ba->real_area[i]<<" "<<ba->overlap_area[i]<<" "<<ba->overlap_index[i] <<std::endl;
    }
    */
    //output the averaged value and variance for each parameter
    fout<<" "<<endl;
    fout<<"averaged value"<<endl;
    fout<<"N_in N_epi N_all organ_area organ_perimeter averaged_in_area averaged_epi_area averaged_perimeter ";
    fout<<"regularity_in_av regularity_epi_av real_area overlap_area overlap_index";
    fout<<"sigma_O sigma_L sigma_O/sigma_L"<<std::endl;
    fout<<averaged_vec_int(ba->N_in)<<" "<<averaged_vec_int(ba->N_epi)<<" "<<averaged_vec_int(ba->N_in)+averaged_vec_int(ba->N_epi)<<" "<<averaged_vec_double(ba->organ_area)<<" "<<averaged_vec_double(ba->organ_perimeter)<<" "<<averaged_vec_double(ba->averaged_inner_area)<<" "<<averaged_vec_double(ba->averaged_epi_area)<<" "<<averaged_vec_double(ba->averaged_perimeter);
    fout<<" "<<averaged_vec_double(ba->averaged_regularity_in)<<" "<<averaged_vec_double(ba->averaged_regularity_epi)<<" "<<averaged_vec_double(ba->real_area)<<" "<<averaged_vec_double(ba->overlap_area)<<" "<<averaged_vec_double(ba->overlap_index);
    fout<<" "<<sigma_O<<" "<<sigma_L<<" "<<sigma_O/sigma_L<<std::endl;
    fout.close();

}

void batch_final_analysis(vector<double> parameter_1, vector<double> parameter_2,vector<double> vec_result,string output_filename){
    std::ofstream fout(output_filename);
    for(int i=0;i<parameter_1.size();i++){
        for(int j=0;j<parameter_2.size();j++){
            fout<<i*parameter_2.size()+j<<" "<<parameter_1[i]<<" "<<parameter_2[j]<<" "<<vec_result[i*parameter_2.size()+j]<<endl;
        }
        fout<<endl;
    }
    fout.close();
}

void simulation_log(vector<time_t> initial_time, vector<time_t> terminal_time){
    char fname[100];
    sprintf(fname,"simulationLog.txt");
    std::ofstream fout(fname);

    tm * initial_tm;
    char initial_time_string [100];
    initial_tm = localtime(&initial_time[0]);
    strftime(initial_time_string, 100, "The whole simulation started at %T, %Y %B %d", initial_tm);

    tm * terminal_tm;
    char terminal_time_string [100];
    terminal_tm = localtime(&terminal_time[terminal_time.size()-1]);
    strftime(terminal_time_string, 100, "The whole simulation terminated at %T, %Y %B %d", terminal_tm);
    

    vector<time_t> simulation_time;
    for(int i=0; i<initial_time.size();i++){
        if(i==0){
            time_t simulation_time_tmp = difftime(terminal_time[terminal_time.size()-1],initial_time[0]);
            simulation_time.push_back(simulation_time_tmp);
        }
        else{
            time_t simulation_time_tmp = difftime(terminal_time[i-1],initial_time[i]);
            simulation_time.push_back(simulation_time_tmp);
        }
    }

    vector<int> simulation_time_h;
    vector<int> simulation_time_m;
    vector<int> simulation_time_s;
    for(int i=0;i<simulation_time.size();i++){
        simulation_time_h.push_back(simulation_time[i]/3600);
        simulation_time_m.push_back((simulation_time[i]%3600)/60);
        simulation_time_s.push_back(simulation_time[i]%60);
    }

    cout<<"The whole simulation (or analysis) calculation took "<<simulation_time_h[0]<<"h "<<simulation_time_m[0]<<"m "<<simulation_time_s[0]<<"s."<<endl;
    fout<<initial_time_string<<endl;
    fout<<terminal_time_string<<endl;
    fout<<"The whole simulation (or analysis) calculation took "<<simulation_time_h[0]<<"h "<<simulation_time_m[0]<<"m "<<simulation_time_s[0]<<"s."<<endl;

    if(simulation_time.size()>1){
        for(int i=1;i<simulation_time.size();i++){
            fout<<"For the "<<i<<"th batch, it took "<<simulation_time_h[i]<<"h "<<simulation_time_m[i]<<"m "<<simulation_time_s[i]<<"s."<<endl;
        }
    }

    fout.close();
}

void simulated_contour(Organ* p_g, std::string filename){
    std::vector<Vertex> simulated_contour;
    for(int vi=0;vi<(int)p_g->p_v.size();vi++){
        if(p_g->p_v[vi]->IsSurface==1){
            Vertex contour_point;
            contour_point.loc.x = p_g->p_v[vi]->loc.x;
            contour_point.loc.y = p_g->p_v[vi]->loc.y;
            simulated_contour.push_back(contour_point);
        }
    }
    double sampling_distance=0.01;
    vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
    vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
    vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
    vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
    cout_fout_debug::fout_vector_vertex(simulated_outline_swapped_sampled,filename);
    //gnuplot::vv(simulated_outline_swapped_sampled);
}

void simulated_contour(Organ* p_g, std::string txt_filename, std::string png_filename){
    std::vector<Vertex> simulated_contour;
    for(int vi=0;vi<(int)p_g->p_v.size();vi++){
        if(p_g->p_v[vi]->IsSurface==1){
            Vertex contour_point;
            contour_point.loc.x = p_g->p_v[vi]->loc.x;
            contour_point.loc.y = p_g->p_v[vi]->loc.y;
            simulated_contour.push_back(contour_point);
        }
    }
    double sampling_distance=0.01;
    vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
    vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
    vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
    vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
    cout_fout_debug::fout_vector_vertex(simulated_outline_swapped_sampled,txt_filename);
    gnuplot::vv(simulated_outline_swapped_sampled,png_filename);
}

void simulated_contour(Organ* p_g){
    char fname[100];
    sprintf(fname,"2dv_contour%010u.vtk", p_g->step);
    std::vector<Vertex> simulated_contour;
    for(int vi=0;vi<(int)p_g->p_v.size();vi++){
        if(p_g->p_v[vi]->IsSurface==1){
            Vertex contour_point;
            contour_point.loc.x = p_g->p_v[vi]->loc.x;
            contour_point.loc.y = p_g->p_v[vi]->loc.y;
            simulated_contour.push_back(contour_point);
        }
    }
    double sampling_distance=0.01;
    vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
    vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
    vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);
    vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
    cout_fout_debug::fout_vector_vertex(simulated_outline_swapped_sampled,std::string(fname));
}

}

namespace file_process{

int getFileNum(const string &path) {   //需要用到<dirent.h>头文件
    int fileNum=0;
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(path.c_str())))
        return fileNum;
    while((ptr=readdir(pDir))!=0){
        if(strcmp(ptr->d_name,".")!=0&&strcmp(ptr->d_name,"..")!=0 )
           fileNum++;
    }
    closedir(pDir);
    return fileNum;
}


}

namespace cout_fout_debug{
    void cout_vector_vertex(vector<Vertex*> vv){
        cout<<"index x y"<<endl;
        for(int vi=0; vi<vv.size();vi++){
            cout<<vi<<" "<<vv[vi]->loc.x<<" "<<vv[vi]->loc.y<<endl;
        }
    }

    void cout_vector_vertex(vector<Vertex> vv){
        cout<<"index x y"<<endl;
        for(int vi=0; vi<vv.size();vi++){
            cout<<vi<<" "<<vv[vi].loc.x<<" "<<vv[vi].loc.y<<endl;
        }
    }

    void fout_vector_vertex(vector<Vertex*> vv, string FileName){
        std::cout<<"vv fout to "<<FileName<<std::endl;
        std::ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        for(int vi=0;vi<vv.size();vi++){
            fout<<vi<<" "<<vv[vi]->loc.x<<" "<<vv[vi]->loc.y<<" "<<vv[vi]->loc.z<<endl;
        }
        fout.close();
    }

    void fout_vector_vertex(vector<Vertex> vv, string FileName){
        std::ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        for(int vi=0;vi<vv.size();vi++){
            fout<<vi<<" "<<vv[vi].loc.x<<" "<<vv[vi].loc.y<<" "<<vv[vi].loc.z<<endl;
        }
        fout.close();
    }

    void cout_vector_double(vector<double> vd){
        cout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            cout<<di<<" "<<vd[di]<<endl;
        }
    }

    void cout_vector_int(vector<int> vi){
        cout<<"index value"<<endl;
        for(int i=0; i<vi.size();i++){
            cout<<i<<" "<<vi[i]<<endl;
        }
    }

    void cout_vector_string(vector<string> vs){
        cout<<"index string"<<endl;
        for(int i=0; i<vs.size();i++){
            cout<<i<<" "<<vs[i]<<endl;
        }
    }

    void cout_vector_line(vector<Line*> vl){
        cout<<"index slope intercept"<<endl;
        for(int li=0; li<vl.size();li++){
            cout<<li<<" "<<vl[li]->slope<<" "<<vl[li]->intercept<<endl;
        }
    }

    void cout_vector_x_y(vector<double> vx,vector<double> vy){
        std::cout<<"v_x v_y"<<std::endl;
        for(int i=0;i<vx.size();i++){
            std::cout<<vx[i]<<" "<<vy[i]<<std::endl;
        }
    }


    void fout_vector_double(vector<double> vd, string FileName){
        std::ofstream fout(FileName);
        fout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            fout<<di<<" "<<vd[di]<<endl;
        }
        fout.close();
    }

    void fout_vector_int_double(vector<int> v, vector<double> vd, string FileName){
        std::ofstream fout(FileName);
        fout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            fout<<di<<" "<<v[di]<<" "<<vd[di]<<endl;
        }
        fout.close();
    }

    void fout_vector_double(vector<double> vd1,vector<double> vd2,vector<double> vd3, string FileName,string vd1_name, string vd2_name, string vd3_name){
        std::ofstream fout(FileName);
        fout<<"index "<<vd1_name<<" "<<vd2_name<<" "<<vd3_name<<endl;
        for(int di=0;di<vd1.size();di++){
            fout<<di<<" "<<vd1[di]<<" "<<vd2[di]<<" "<<vd3[di]<<endl;
        }
        fout.close();
    }

    void fout_vector_int(vector<int> vd, string FileName){
        std::ofstream fout(FileName);
        fout<<"index value"<<endl;
        for(int di=0;di<vd.size();di++){
            fout<<di<<" "<<vd[di]<<endl;
        }
        fout.close();
    }

    void fout_vector_x_y(vector<double> vx,vector<double> vy,string FileName){
        ofstream fout(FileName);
        fout<<"x y"<<std::endl;
        for(int i=0;i<vx.size();i++){
            fout<<vx[i]<<" "<<vy[i]<<std::endl;
        }
        fout.close();
    }


    void fout_pair_double(vector<pair<double,double>> pdd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y"<<endl;
        for(int pi=0;pi<pdd.size();pi++){
            fout<<pi<<" "<<pdd[pi].first<<" "<<pdd[pi].second<<endl;
        }
        fout.close();
    }

    void fout_tuple_double(vector<tuple<double,double,double>> tddd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        int index=0;
        for(const auto& ti : tddd){
            fout<<index<<" "<<get<0>(ti)<<" "<<get<1>(ti)<<" "<<get<2>(ti)<<endl;
            index++;
        }
        fout.close();
    }

    void fout_tuple_double(vector<tuple<int,double,double>> tddd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        int index=0;
        for(const auto& ti : tddd){
            fout<<index<<" "<<get<0>(ti)<<" "<<get<1>(ti)<<" "<<get<2>(ti)<<endl;
            index++;
        }
        fout.close();
    }

        void fout_tuple_double(vector<tuple<int,int,double>> tddd,string FileName)
    {
        ofstream fout(FileName);
        fout<<"index x y z"<<endl;
        int index=0;
        for(const auto& ti : tddd){
            fout<<index<<" "<<get<0>(ti)<<" "<<get<1>(ti)<<" "<<get<2>(ti)<<endl;
            index++;
        }
        fout.close();
    }

    void fout_heatmap_data(std::vector<int> xv,std::vector<int> yv,std::vector<double> value,std::string filename){
        ofstream fout(filename);
        fout<<"index x y z"<<endl;
        for(int i=0 ; i<xv.size();i++){
            for(int j=0;j<yv.size();j++){
                fout<<j*xv.size()+i<<" "<<xv[i]<<" "<<yv[j]<<" "<<value[j*xv.size()+i]<<endl;
            }
        }
        fout.close();
    }

        void fout_heatmap_data(std::vector<double> xv,std::vector<double> yv,std::vector<double> value,std::string filename){
        ofstream fout(filename);
        fout<<"index x y z"<<endl;
        for(int i=0 ; i<xv.size();i++){
            for(int j=0;j<yv.size();j++){
                fout<<j*xv.size()+i<<" "<<xv[i]<<" "<<yv[j]<<" "<<value[j*xv.size()+i]<<endl;
            }
        }
        fout.close();
    }

    std::vector<std::pair<double,double>> fin_vec_pair_double(std::string input_file){
        std::vector<std::pair<double,double>> vdd;
        ifstream fin(input_file, ios::in);
        if(!fin.is_open()){            
            cout<<"Error: missing "<<input_file<<endl;
            exit(1);
        }
        
        while(!fin.eof()){
            double tmp;
            fin>>tmp;
            double vd1,vd2;
            fin>>vd1;
            fin>>vd2;
            vdd.push_back(std::make_pair(vd1,vd2));
        }
        return vdd;
    }
}

namespace gnu_plot{
    void dotplot(vector<double> vd, int point_size, string image_name){
        FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe != NULL) {
        // Set the terminal to png and output to a file
        fprintf(pipe, "set terminal png\n");
        fprintf(pipe, "set output '%s'\n",image_name.c_str());
        //fprintf(pipe, "set yrange [-4:8]\n");

        // Plot with points
        fprintf(pipe, "plot '-' with points pt 7 ps %d lc rgb 'dark-green' notitle \n",point_size);
        // Output data points
        for(int i = 0; i < vd.size(); i++) {
            fprintf(pipe, "%d %f\n", i, vd[i]);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
    }
    }

    void organ_contour_plot(vector<Vertex*> vv, int point_size, string save_image_name){
        FILE *pipe = popen("gnuplot -persist", "w");

        if(pipe != NULL){
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            fprintf(pipe, "set size ratio -1\n");
            
        }
        // Plot with points
        fprintf(pipe,"set xrange [-0.8:0.8]\n");
        fprintf(pipe,"set yrange [0:1]\n");
        fprintf(pipe,"set size ratio 0.625\n");
        fprintf(pipe,"unset tics\n");
        fprintf(pipe,"unset border\n");
        fprintf(pipe, "plot '-' with points pt 7 ps %d lc 'dark-green' notitle \n",point_size);
        // Output data points
        for(int i = 0; i < vv.size(); i++) {
            fprintf(pipe, "%f %f\n", vv[i]->loc.x, vv[i]->loc.y);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);

    }

    void organ_contour_plot(vector<Vertex*> vv){
        FILE *pipe = popen("gnuplot -persist", "w");

        if(pipe != NULL){
            // Set the terminal to png and output to a file
            //fprintf(pipe, "set terminal png\n");
            //fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            fprintf(pipe, "set size ratio -1\n");
            
        }
        // Plot with points
        fprintf(pipe, "plot '-' with points pt 7 notitle \n");
        // Output data points
        for(int i = 0; i < vv.size(); i++) {
            fprintf(pipe, "%f %f\n", vv[i]->loc.x, vv[i]->loc.y);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);

    }

    void curvature_plot(vector<double> vd, string save_image_fileName){
        FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe != NULL) {
        // Set the terminal to png and output to a file
        fprintf(pipe, "set terminal pngcairo\n");
        fprintf(pipe, "set output '%s'\n",save_image_fileName.c_str());
        //fprintf(pipe, "set yrange [-10:10]\n");
        fprintf(pipe,"set arrow from 0,0 to 100,0 nohead dashtype 2 lw 4 lc 'grey'\n");
        // Plot with points
        fprintf(pipe, "plot '-' with points pt 7 ps 1.5 lc rgb 'dark-green' notitle \n");
        // Output data points
        for(int i = 0; i < vd.size(); i++) {
            fprintf(pipe, "%d %f\n", i, vd[i]);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
    }
    }

    void panel_plot(vector<string> plot_filenames, int row, int column, string save_panel_name){
        
        cout<<"Now we are going to arrange plots into panel"<<endl;
        cout_fout_debug::cout_vector_string(plot_filenames);
        vector<cv::Mat> images;
        for (const auto& filename : plot_filenames) {
            images.push_back(cv::imread(filename));
        }

        // Create the panel image
        cv::Mat panel;
        
        // Arrange the images into a panel
        int images_per_row = row;
        int images_per_column = column;
        
        for (int i = 0; i < images_per_column; ++i) {
            cv::Mat row;
            for (int j = 0; j < images_per_row; ++j) {
                if (j == 0) {
                    row = images[i * images_per_row + j];
                } else {
                    cv::hconcat(row, images[i * images_per_row + j], row);
                }
                std::cout<<"successfully add the "<<i<<","<<j<<" images"<<std::endl;
            }
            if (i == 0) {
                panel = row;
            } else {
                cv::vconcat(panel, row, panel);
            }
        }

        // Save the panel image
        cv::imwrite(save_panel_name, panel);
    }

    void panel_plot_triangle(vector<string> plot_filenames, int row, int column, string save_panel_name){
        
        cout<<"Now we are going to arrange plots into panel"<<endl;
        cout_fout_debug::cout_vector_string(plot_filenames);
        vector<cv::Mat> images;
        for (const auto& filename : plot_filenames) {
            images.push_back(cv::imread(filename));
        }

        // Create the panel image
        cv::Mat panel;
        cv::Mat blank;
        for(int i=0;i<plot_filenames.size();i++){
            if(images[i].empty()!=0){
                cv::Mat blank = cv::Mat::zeros(images[i].rows, images[i].cols, images[i].type());
            }
        }
        
        // Arrange the images into a panel
        int images_per_row = row;
        int images_per_column = column;
        
        for (int i = 0; i < images_per_column; ++i) {
            cv::Mat row;
            for (int j = 0; j < images_per_row; ++j) {
                if(images[i * images_per_row + j].empty()==0){
                    cv::hconcat(row, blank, row);
                    std::cout<<i<<","<<j<<" has no image"<<std::endl;
                }
                else{
                    if (j == 0) {
                        row = images[i * images_per_row + j];
                    } else {
                        cv::hconcat(row, images[i * images_per_row + j], row);
                    }
                    std::cout<<"successfully add the "<<plot_filenames[i * images_per_row + j]<<" images"<<std::endl;
                }
            }
            if (i == 0) {
                panel = row;
            } else {
                cv::vconcat(panel, row, panel);
            }
        }

        // Save the panel image
        cv::imwrite(save_panel_name, panel);
    }

    void gaussian_plot(double mu,double sigma, string save_image_name){
        FILE *pipe = popen("gnuplot -persist", "w");

        if(pipe != NULL){
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            fprintf(pipe, "set xrange [0:180]\n");
            fprintf(pipe, "set yrange [0:1.2]\n");
            fprintf(pipe, "set samples 1000\n");
            fprintf(pipe, "set xtics font \"Arial,18\"\n");
            fprintf(pipe, "set ytics font \"Arial,18\"\n");
            fprintf(pipe, "set title 'mu=%.2f,sigma=%.2f' font \"Arial,24\"\n",mu,sigma);
            
        }
        // Plot with points
        fprintf(pipe, "plot exp(-(x-%2.f)**2/(2*(%.2f)**2)) lw 4 lc rgb 'dark-green' notitle\n",mu,sigma);
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
    }

    void histogram(std::vector<double> vd_data, std::string save_image_name){
            double min_data = vd_data[0];
            double max_data = vd_data[0];
            for(int i=0;i<vd_data.size();i++){
                if(vd_data[i]<min_data) min_data=vd_data[i];
                else if(vd_data[i]>max_data) max_data=vd_data[i];
            }
            std::cout<<"Maximum cell area: "<<max_data<<"; minimum cell area: "<<min_data<<std::endl;

            double binData = (max_data-min_data)/20.0;
            double minData = min_data;
            double maxData = max_data;
            double binCount = (maxData-minData)/binData;
            vector<int> binDatas(binCount,0);
            for(double Data : vd_data){
                int binIndex = static_cast<int>((Data-minData)/binData);
                if(binIndex>=0&&binIndex<binCount)
                    ++binDatas[binIndex];
            }
            cout_fout_debug::cout_vector_int(binDatas);
                FILE *pipe = popen("gnuplot -persist", "w");

            if (pipe != NULL) {
                // Set the terminal to png and output to a file
                fprintf(pipe, "set terminal png\n");
                fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
                fprintf(pipe, "set style fill solid border lc rgb 'black'\n");
                //fprintf(pipe, "set yrange [0:20]\n");
                fprintf(pipe, "set xtics '%.2f' \n",2*binData);

                // Plot with points
                fprintf(pipe, "plot '-' with boxes lw 2 lc rgb 'light-cyan' notitle \n");
                // Output data points
                for(int i = 0; i < binDatas.size(); i++) {
                    double binMiddle = minData + i*binData + binData/2;
                    fprintf(pipe, "%.2f %d\n", binMiddle, binDatas[i]);
                }
                
                // End of data
                fprintf(pipe, "e\n");
                
                // Flush the pipe
                fflush(pipe);

                // Close the pipe
                pclose(pipe);
        }
    }

    void multiple_time_lapse_contour_txt_plot(std::vector<std::string> input_files, std::string save_image_name, int cell_number){
        /*
        std::vector<std::vector<std::pair<double,double>>> contours;
        for(int i=0;i<input_files.size();i++){
            std::vector<std::pair<double,double>> contour;
            std::ifstream file(input_files[i]);
            if (!file.is_open()) {
                std::cerr << "Unable to open file ("<<input_files[i]<<")." << std::endl;
            }
            std::string line;
            file>>line;
            file>>line;
            file>>line;
            file>>line;
            double index,value1,value2,value3;
            // Read file line by line
            while (file >> index >> value1 >> value2 >> value3) {
                contour.push_back(std::make_pair(value1,value2));
            }
            file.close();
        }
        */
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            fprintf(pipe, "set yrange [-0.1:1.1]\n");
            fprintf(pipe, "set xrange [-0.6:0.6]\n");
            fprintf(pipe, "set title 'contours of cell number %d'\n",cell_number);
            // Plot with points
            fprintf(pipe, "plot ");
            // Output data points
            for(size_t repeat_i = 0; repeat_i<input_files.size(); repeat_i++){
                fprintf(pipe, " '%s' using 2:3 with points pt 7 ps 1 notitle", input_files[repeat_i].c_str());
                if(repeat_i<input_files.size()-1){
                    fprintf(pipe,", ");
                }
            }
            
            // End of data
            fprintf(pipe, "\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    void organ_contour_plot(std::string input_file,std::string save_image_name, int cell_number){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            fprintf(pipe, "set yrange [-0.1:1.1]\n");
            fprintf(pipe, "set xrange [-0.6:0.6]\n");
            fprintf(pipe, "set size ratio 1\n");
            fprintf(pipe, "set title 'contours of cell number %d'\n",cell_number);
            // Plot with points
            fprintf(pipe, "plot ");
            // Output data points
            fprintf(pipe, " '%s' using 2:3 with points pt 7 ps 1 notitle", input_file.c_str());

            
            // End of data
            fprintf(pipe, "\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    void organ_width_plot(std::string input_file,std::string save_image_name, int cell_number){
        FILE *pipe = popen("gnuplot -persist", "w");

        if (pipe != NULL) {
            // Set the terminal to png and output to a file
            fprintf(pipe, "set terminal png\n");
            fprintf(pipe, "set output '%s'\n",save_image_name.c_str());
            //fprintf(pipe, "set yrange [-0.1:1.1]\n");
            //fprintf(pipe, "unset xrange\n");
            fprintf(pipe, "set title 'contours of cell number %d'\n",cell_number);
            // Plot with points
            fprintf(pipe, "plot ");
            // Output data points
            fprintf(pipe, " '%s' using 2 with points pt 7 ps 1 notitle", input_file.c_str());

            
            // End of data
            fprintf(pipe, "\n");
            
            // Flush the pipe
            fflush(pipe);

            // Close the pipe
            pclose(pipe);
        }
    }

    void heatmap(std::string input_txt_file, std::string save_image_name, int cell_number){
        FILE *pipe = popen("gnuplot -persist", "w"); // Open a pipe to gnuplot
        if (pipe == NULL) {
            std::cerr << "Could not open pipe" << std::endl;
            exit(1);
        }

        fprintf(pipe, "set term png\n"); // Set the output terminal to PNG
        fprintf(pipe, "set output '%s'\n",save_image_name.c_str()); // Set the output file name
        fprintf(pipe, "set view map\n"); // Set the view to 2D heatmap
        fprintf(pipe, "set size square\n"); // Set aspect ratio to equal
        fprintf(pipe, "set palette rgbformulae 10,13,35\n");   
        //fprintf(pipe, "set logscale cb\n");
        fprintf(pipe, "set cbrange [0:0.5]\n");
        // Send the data to be plotted
        fprintf(pipe, "set title 'dynamics difference index of cell number %d'\n",cell_number);
        fprintf(pipe, "splot '%s' using 2:3:4 with image notitle\n",input_txt_file.c_str());
        fflush(pipe);
        pclose(pipe);
    }

    void heatmap(std::vector<int> x_axis, std::vector<int> y_axis, std::vector<double> values, std::string save_image_name){
        FILE *pipe = popen("gnuplot -persist", "w"); // Open a pipe to gnuplot
        if (pipe == NULL) {
            std::cerr << "Could not open pipe" << std::endl;
            exit(1);
        }

        fprintf(pipe, "set term png\n"); // Set the output terminal to PNG
        fprintf(pipe, "set output '%s'\n",save_image_name.c_str()); // Set the output file name
        fprintf(pipe, "set view map\n"); // Set the view to 2D heatmap
        fprintf(pipe, "set size square\n"); // Set aspect ratio to equal
        fprintf(pipe, "set palette rgbformulae 10,13,35\n");   
        //fprintf(pipe, "set logscale cb\n");
        fprintf(pipe, "set cbrange [0:0.5]\n");
        // Send the data to be plotted
        fprintf(pipe, "splot '-' matrix with image\n");

        int cols = x_axis.size();
        int rows = y_axis.size();
        double data[cols][rows];

        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                data[i][j] = values[i*rows+j]; // Replace this with actual data
                fprintf(pipe, "%d %d %.3f ",x_axis[i], y_axis[j], data[i][j]);
            }
            fprintf(pipe, "\n");
        }

        fflush(pipe);
        pclose(pipe);
    }
}