#include "../include/image2dv.h"



namespace Image_2dv{
    void read_image(std::string image_file, Image2dv* image){
        cv::Mat img = cv::imread(image_file, cv::IMREAD_GRAYSCALE);
        if(img.empty()){
            std::cerr<<"Fatal error: Could not read the image("<<image_file<<")"<<std::endl;
            exit(1);
        }

        for(int i=0;i<img.rows;++i){
            std::vector<Point2dv> v_pt;
            for(int j=0;j<img.cols;++j){
                double pixelValue = static_cast<double>(img.at<uchar>(j, i));
                v_pt.emplace_back(i, j, pixelValue);
            }
            image->pt.push_back(v_pt);
        }
        image->height = img.rows;
        image->width = img.cols;
    }

    void draw_image_from_points(Image2dv* image) {
    // Assuming the points fit within the dimensions of img
        cv::Mat img(image->width,image->height,CV_8UC1,cv::Scalar(0));

        for (auto v_point : image->pt) {
            for(auto point : v_point){
                // Cast the coordinates to int for pixel access
                int ix = static_cast<int>(point.x);
                int iy = static_cast<int>(point.y);

                // Set the pixel value, assuming a single-channel (grayscale) image.
                // If your image has more channels, you need to set a value for each channel.
                img.at<uchar>(iy, ix) = static_cast<uchar>(point.pixel_value);
            }
            
        }
        cv::imshow("Image",img);
        cv::waitKey(0);
    }

    Point2dv start_point_in_polygon(Image2dv* image){
        for(int i=0;i<image->height;i++){
            for(int j=0;j<image->width;j++){
                if(image->pt[i][j].checked==0&&image->pt[i][j].pixel_value>220){
                    image->pt[i][j].pixel_value=255;
                    //std::cout<<"the start point in polygon: "<<i<<","<<j<<std::endl;
                    return image->pt[i][j];
                }
            }
        }
        return Point2dv(-1,-1,-1);
    }

    void add_white_pixel_to_polygon(Image2dv* image,int pt_x,int pt_y, int pg_i){
        image->pt[pt_x][pt_y].checked=1;
        image->pt[pt_x][pt_y].polygon_index=pg_i;
        //add upper left point
        if(pt_x!=0&&pt_y!=0){
            if(image->pt[pt_x-1][pt_y-1].pixel_value>220&&image->pt[pt_x-1][pt_y-1].checked==0&&image->pt[pt_x-1][pt_y-1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x-1][pt_y-1]);
                image->pt[pt_x-1][pt_y-1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x-1<<","<<pt_y-1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_x!=0){
            //add upper point
            if(image->pt[pt_x-1][pt_y].pixel_value>220&&image->pt[pt_x-1][pt_y].checked==0&&image->pt[pt_x-1][pt_y].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x-1][pt_y]);
                image->pt[pt_x-1][pt_y].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x-1<<","<<pt_y<<" is added to the polygon "<<pg_i<<std::endl;;
            }
        }
        if(pt_x!=0&&pt_y!=image->width-1){
            //add upper right point
            if(image->pt[pt_x-1][pt_y+1].pixel_value>220&&image->pt[pt_x-1][pt_y+1].checked==0&&image->pt[pt_x-1][pt_y+1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x-1][pt_y+1]);
                image->pt[pt_x-1][pt_y+1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x-1<<","<<pt_y+1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_y!=0){
            //add left point
            if(image->pt[pt_x][pt_y-1].pixel_value>220&&image->pt[pt_x][pt_y-1].checked==0&&image->pt[pt_x][pt_y-1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x][pt_y-1]);
                image->pt[pt_x][pt_y-1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x<<","<<pt_y-1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_x!=image->height-1&&pt_y!=0){
            //add lower left point
            if(image->pt[pt_x+1][pt_y-1].pixel_value>220&&image->pt[pt_x+1][pt_y-1].checked==0&&image->pt[pt_x+1][pt_y-1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x+1][pt_y-1]);
                image->pt[pt_x+1][pt_y-1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x+1<<","<<pt_y-1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_x!=image->height-1&&pt_y!=image->width-1){
            //add lower right point
            if(image->pt[pt_x+1][pt_y+1].pixel_value>220&&image->pt[pt_x+1][pt_y+1].checked==0&&image->pt[pt_x+1][pt_y+1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x+1][pt_y+1]);
                image->pt[pt_x+1][pt_y+1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x+1<<","<<pt_y+1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_x!=image->height-1){
            //add lower point
            if(image->pt[pt_x+1][pt_y].pixel_value>220&&image->pt[pt_x+1][pt_y].checked==0&&image->pt[pt_x+1][pt_y].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x+1][pt_y]);
                image->pt[pt_x+1][pt_y].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x+1<<","<<pt_y<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
        if(pt_y!=image->width-1){
            //add right point
            if(image->pt[pt_x][pt_y+1].pixel_value>220&&image->pt[pt_x][pt_y+1].checked==0&&image->pt[pt_x][pt_y+1].polygon_index==-1){
                image->pg[pg_i].pt.push_back(image->pt[pt_x][pt_y+1]);
                image->pt[pt_x][pt_y+1].polygon_index=pg_i;
                //std::cout<<"Point "<<pt_x<<","<<pt_y+1<<" is added to the polygon "<<pg_i<<std::endl;
            }
        }
    }

    bool check_contour_finding_end(Image2dv* image){
        bool check_end=1;
        for(int i=0;i<image->height;i++){
            for(int j=0;j<image->width;j++){
                if(image->pt[i][j].polygon_index==-1&&image->pt[i][j].pixel_value>220){
                    check_end=0;
                }
            }
        }

        return check_end;
    }

    void find_polygon_from_image(Image2dv* image){
        int polygon_index=0;
        bool check_end=0;
        do{
            Polygon2dv polygon_i;
            Point2dv start_point=start_point_in_polygon(image);
            polygon_i.pt.push_back(start_point);
            image->pg.push_back(polygon_i);
            for(int i=0;i<image->pg[polygon_index].pt.size();i++){
                if(image->pt[image->pg[polygon_index].pt[i].x][image->pg[polygon_index].pt[i].y].checked==0&&image->pt[image->pg[polygon_index].pt[i].x][image->pg[polygon_index].pt[i].y].pixel_value>220)
                    add_white_pixel_to_polygon(image,image->pg[polygon_index].pt[i].x,image->pg[polygon_index].pt[i].y,polygon_index);
            }
            
            for(int i=0;i<image->height;i++){
                for(int j=0;j<image->width;j++){
                    if(image->pt[i][j].polygon_index==polygon_index){
                        image->pt[i][j].pixel_value=0;
                    }
                }
            }
            
            Image_2dv::draw_image_from_points(image);
            
            polygon_index++;
        }while(!check_contour_finding_end(image));

        
        std::cout<<"Number of cells found "<<image->pg.size()-1<<std::endl;
    }

    void find_polygon_from_image_(Image2dv* image){


        //1. start from the upper left white pixel to find the background polygon
        Polygon2dv background;
        Point2dv background_point = start_point_in_polygon(image);
        
        background.pt.push_back(background_point);
        image->pg.push_back(background);
        //2. add the first white pixel from the eight neighbor point to the background polygon
        add_white_pixel_to_polygon(image,background_point.x,background_point.y,0);
        //add_white_pixel_to_polygon(image,1,0,0);
        //std::cout<<"image->pg[0].pt.size() "<<image->pg[0].pt.size()<<std::endl;
        //3. add all the first white pixel from the unchecked neighbor points to the polygon
        for(int i=0;i<image->pg[0].pt.size();i++){
            if(image->pt[image->pg[0].pt[i].x][image->pg[0].pt[i].y].checked==0&&image->pt[image->pg[0].pt[i].x][image->pg[0].pt[i].y].pixel_value>220){
                add_white_pixel_to_polygon(image,image->pg[0].pt[i].x,image->pg[0].pt[i].y,0);
                //std::cout<<"image,image->pg[0].pt[i].x,image->pg[0].pt[i].y,0); "<<image->pg[0].pt[i].x<<","<<image->pg[0].pt[i].y<<std::endl;
            }
        }
        //4. find the start point of second polygon
        
        Point2dv start_point = start_point_in_polygon(image);
        Polygon2dv polygon_1;
        polygon_1.pt.push_back(start_point);
        image->pg.push_back(polygon_1);
        for(int i=0;i<image->pg[1].pt.size();i++){
            if(image->pt[image->pg[1].pt[i].x][image->pg[1].pt[i].y].checked==0&&image->pt[image->pg[1].pt[i].x][image->pg[1].pt[i].y].pixel_value>220){
                add_white_pixel_to_polygon(image,image->pg[1].pt[i].x,image->pg[1].pt[i].y,1);
                std::cout<<"image,image->pg[1].pt[i].x,image->pg[1].pt[i].y,0); "<<image->pg[1].pt[i].x<<","<<image->pg[1].pt[i].y<<std::endl;
            }
        }

        //debug check if the background are detected correctly
        /*
        for(int i=0;i<image->height;i++){
            for(int j=0;j<image->width;j++){
                if(image->pt[i][j].polygon_index==1){
                    image->pt[i][j].pixel_value=0;
                }
            }
        }
        
        Image_2dv::draw_image_from_points(image);
        */
    }

     void statistics_of_polygons_(Image2dv* image){
        /*
        std::vector<double> cell_area_vec;
        std::vector<double> cell_length_vec;
        std::vector<double> cell_width_vec;
        std::vector<double> cell_expansion_degree_vec;
        std::vector<double> cell_expansion_direction_vec;
        //calculate the area of polygons

        for(int pi=0;pi<image->pg.size();pi++){
            double cell_area = image->pg[pi].pt.size();
            cell_area_vec.push_back(cell_area);
            double cell_length=0;
            double cell_expansion_direction, cell_length_slope;
            double length_vertex_first_x, length_vertex_first_y, length_vertex_second_x,length_vertex_second_y;
            for(int vi=0;vi<image->pg[pi].pt.size()){
                for(int vj=0;vj<image->pg[pi].pt.size()){
                    double distance=sqrt((image->pg[pi].pt[vi].x-image->pg[pi].pt[vj].x)*(image->pg[pi].pt[vi].x-image->pg[pi].pt[vj].x)+(image->pg[pi].pt[vi].y-image->pg[pi].pt[vj].y)*(image->pg[pi].pt[vi].y-image->pg[pi].pt[vj].y));
                    cell_length=std::max(cell_length,distance);
                    length_vertex_first_x = image->pg[pi].pt[vi].x;
                    length_vertex_first_y = image->pg[pi].pt[vi].y;
                    length_vertex_second_x = image->pg[pi].pt[vj].x;
                    length_vertex_second_y = image->pg[pi].pt[vj].y;
                }
            }
            cell_expansion_direction=atan2(length_vertex_second_y-length_vertex_first_y,length_vertex_second_x-length_vertex_first_x);
            cell_length_slope = (length_vertex_second_y-length_vertex_first_y)/(length_vertex_second_x-length_vertex_first_x);
            cell_length_vec.push_back(cell_length);
            cell_expansion_direction_vec.push_back(cell_expansion_direction);

            double cell_width_slope = -1.0/cell_length_slope;
            int width_axis_precision = 0.1;
            double cell_width=0;
            bool width_found=0;
            do
            {
                
            } while ();
            
        //std::cout<<"organ_width: "<<organ_width<<endl;
            for(int vi=0;vi<image->pg[pi].pt.size()){
                for(int vj=0;vj<image->pg[pi].pt.size()){
                    double distance=sqrt((image->pg[pi].pt[vi].x-image->pg[pi].pt[vj].x)*(image->pg[pi].pt[vi].x-image->pg[pi].pt[vj].x)+(image->pg[pi].pt[vi].y-image->pg[pi].pt[vj].y)*(image->pg[pi].pt[vi].y-image->pg[pi].pt[vj].y));
                    double slope = (image->pg[pi].pt[vi].y-image->pg[pi].pt[vj].y)/(image->pg[pi].pt[vi].x-image->pg[pi].pt[vj].x)
                    if(abs(slope-cell_width_slope)<0.1){
                        cell_width=std::max(cell_width,distance);
                    }
                }
            }

            double cell_expansion_degree = cell_length/cell_width;
            cell_expansion_degree_vec.push_back(cell_expansion_degree);
        }

        //cout all results
        for(int pi=0;pi<image->pg.size();pi++){
            std::cout<<"Cell "<<pi<<" area "<<cell_area_vec[pi]<<" expansion direction "<<cell_expansion_direction_vec[pi]<<" expansion degree "<<cell_expansion_degree_vec[pi]<<std::endl;
        }
        */
     }

    void removeSmallerObjects(cv::Mat& src, cv::Mat& dst) {
        // Convert the image to grayscale
        cv::Mat gray;
        cvtColor(src, gray, cv::COLOR_BGR2GRAY);

        // Threshold the image to binary
        cv::Mat binary;
        threshold(gray, binary, 1, 255, cv::THRESH_BINARY);

        // Find all contours in the binary image
        vector<vector<cv::Point>> contours;
        findContours(binary, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

        // Find the largest contour
        int largest_contour_index = -1;
        double largest_area = 0;
        for (size_t i = 0; i < contours.size(); i++) {
            double area = cv::contourArea(contours[i]);
            if (area > largest_area) {
                largest_area = area;
                largest_contour_index = i;
            }
        }

        // Create a mask for the largest contour
        cv::Mat mask = cv::Mat::zeros(src.size(), CV_8UC1);
        if (largest_contour_index != -1) {
            drawContours(mask, contours, largest_contour_index, cv::Scalar(255), cv::FILLED);
        }

        // Use the mask to keep only the largest object
        dst = cv::Mat::zeros(src.size(), src.type());
        src.copyTo(dst, mask);
    }

}

