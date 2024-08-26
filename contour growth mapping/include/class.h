#ifndef _CLASS_H
#define _CLASS_H

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>


class Point_element;
class Growth_vector;
class Contour;

class Point_element{
 public:
    double x;
    double y;
};

class Contour{
 public:

    double initial_radius;
    double growth_rate;
    int end_step;

    Point_element center;
    
    std::vector<Point_element*> pt;
    std::vector<Point_element*> normalized;
    std::vector<Point_element*> sampled;


    double min_x;
    double min_y;
    double max_x;
    double max_y;
    // Destructor
    ~Contour() {
        for (auto p : pt) { // Iterate over all pointers in the vector
            delete p; // Delete the object pointed to by the pointer
        }
        for (auto p:normalized){
            delete p;
        }
        for (auto p:sampled){
            delete p;
        }
        pt.clear(); // Optional: Clear the vector, not strictly necessary since the object is being destroyed
        normalized.clear();
        sampled.clear();
    }

};



class Growth_vector{
 public:
    double x;
    double y;
};

#endif