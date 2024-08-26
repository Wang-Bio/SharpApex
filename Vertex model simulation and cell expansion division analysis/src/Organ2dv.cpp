#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"
#include "../include/class2dv.h"

Organ::~Organ(){
        for(Vertex* ptr_v : p_v){
                delete ptr_v;
                ptr_v = nullptr;
        }
        for(Line* ptr_l : p_l) {
                delete ptr_l;
                ptr_l=nullptr;
        }
        for(Cell* ptr_c : p_c) {
                delete ptr_c;
                ptr_c=nullptr;
        }
        for(DivisionRecord* ptr_dr : d_r){
                delete ptr_dr;
                ptr_dr=nullptr;
        }
}

