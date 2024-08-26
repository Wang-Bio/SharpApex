# 4. Geometric Analysis of Vertices, Lines, Cells and Organ

In this section, I'd like to summarize functions in [geo.h](../../include/geo.h), [geo.cpp](../../src/geo.cpp) and member functions in [class.h](../../include/class.h)

## 4.1 vertex_geo namespace and class Vertex
Now we are able to calculate the 
* distance between two vertices
* relationship between two vertices (=0, different; =1, identical)
* the collinear relationship between three vertices (=0, not collinear; =1, collinear), 
* conversion between Cartesian coordinates and Polar coordinates

The namespace vertex_geo includes functions:
* double vertex_distance(Vertex* v1, Vertex* v2): returns the euclidean distance between two vertices.
* bool vertex_relationship(Vertex* v1, Vertex* v2): returns 0, if they locate at different positions, returns 1, if they locate at the same position

The member function of class Vertex includes: 
* void print_Cartesian(void):       print out the (x,y) coordinates of the Vertex
* void print_Polar(void):              print out the (r,θ) of the Vertex
* double distance_from_vertex(Vertex): calculate the distance between the Vertex and another Vertex
* void Cartesian_to_Polar(void): calculate the r,θ from the x,y
* void Polar_to_Cartesian(void): calculate the x,y from the r,θ
* bool collinear_points(Vertex,Vertex):
* bool same_vertex(Vertex): if the distance between two vertices is smaller than 1E-5, they are defined as the same vertex

## 4.2 line_geo and class Line
The line (inifinite in both directions), ray (inifinite in one direction), and line segment (has two end points) are all categorized as line in this code.  
Now we are able to calculate the:  
* relationship between two lines (=0, intersected; =1, parallel; =2, perpendicular; =3, identical)
* (unfinished) relationship between two line_segments (=0, intersected; =1, parallel; =2, perpendicular; =3, identical; =4, overlapping; =5, touching; =6, disjoint)
* distance between two lines, two line_segments
* finding the crosspoints between two lines, two line_segments (if they are intersected, or the overlapping length if two line segments are overlapping or identical)

The namespace line_geo includes functions:
* double line_length(Organ* p_g, int li): returns the length of line li in organ p_g.
* _vec<double> line_cell_wall_intersection(double slope, double intercept, Organ* p_g, int li): returns the position of their intersection, if _vec.z=1.0, they are parallel, if _vec.z=-1.0, they simply have no intersection

The namespace segment_geo includes functions:
* Distance_point_line distance_line_segment_to_vertex(Organ*,int,Vertex): return the distance and the closest vertex on the line_segment

The member function of class Line includes: 
* double cal_length(Organ p_g): returns the length of the line 
* void calc_slope_intercept(): calculate the slope and intercept 



## 4.3 cell_geo namespace and class Cell
Now we are able to calculate the
* area and perimeter of the cell
* center of the cell
* regularity of the cell
and do the following process to the cell
* order the vertices and lines of a cell into anticlockwise direction

The namespace cell_geo includes functions:
* _vec<double> cell_center(Organ* p_g, int ci): returns the center of cell ci in organ p_g.
* double cell_area(Organ* p_g, int ci)
* double cell_perimeter(Organ* p_g, int ci)
* vector<int> cell_counterclock_line(Organ* p_g, int ci)
* 

## 4.4 organ_geo namespace and class Organ
Now we are able to calculate the characteristics of organ
* length of all lines of the organ and the organ perimeter
* area of all cells of the organ and the organ area
* center of organ
* regularity, circularity and overlaps of the organ
* length, width, and leaf index of the organ
* potential energy of the organ
and identify
* the epidermal identity of each cell (=1, epidermal cells; =0, inner cells)

## 4.5 boundary_geo and class Boundary
Now we are able to calculate the characteristics of organ boundary: 
* calculate the curvature by circle fitting

The boundary_geo namespace includes the function:
* double surface_vertex_curvature_kasa(vector<Vertex*> vv_fitting): returns the curvature value made by circle fitting based on Kasa method. 



## 4.6 circle_geo namespace and class Circle
Now we are able to calculate the characteristics and relationships of circles:
* relationships and intersections (if existed) between circle and line, circle and line segment


The circle_geo namespace includes the functions:
* Intersection_relationship intersection_line_circle(Organ*,int,Circle);
* Intersection_relationship intersection_line_segment_circle(Organ*,int,Circle);
* Intersection_relationship intersection_line_segment_circle(Line_Segment,Circle);