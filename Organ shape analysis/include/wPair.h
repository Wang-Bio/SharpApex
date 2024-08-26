#ifndef WPAIR_H
#define WPAIR_H

#include <string>
#include <vector>
#include <iostream>

namespace wPair{
    template <typename T>
    T findMaxX(std::vector<std::pair<T,T>> vecPair){
        if(vecPair.size()==0){
            std::cerr<<"There is no element in vecPair (function findMaxX)"<<std::endl;
            exit(1);
        }
        T max_tmp=vecPair[0].first;
        for(const auto& pair_tmp:vecPair){
            if(max_tmp<pair_tmp.first){
                max_tmp=pair_tmp.first;
            }
        }
        return max_tmp;
    }

    template <typename T>
    T findMinX(std::vector<std::pair<T,T>> vecPair){
        if(vecPair.size()==0){
            std::cerr<<"There is no element in vecPair (function findMaxX)"<<std::endl;
            exit(1);
        }
        T min_tmp=vecPair[0].first;
        for(const auto& pair_tmp:vecPair){
            if(min_tmp>pair_tmp.first){
                min_tmp=pair_tmp.first;
            }
        }
        return min_tmp;
    }

    template <typename T>
    T findMaxY(std::vector<std::pair<T,T>> vecPair){
        if(vecPair.size()==0){
            std::cerr<<"There is no element in vecPair (function findMaxX)"<<std::endl;
            exit(1);
        }
        T max_tmp=vecPair[0].second;
        for(const auto& pair_tmp:vecPair){
            if(max_tmp<pair_tmp.second){
                max_tmp=pair_tmp.second;
            }
        }
        return max_tmp;
    }

    template <typename T>
    T findMinY(std::vector<std::pair<T,T>> vecPair){
        if(vecPair.size()==0){
            std::cerr<<"There is no element in vecPair (function findMaxX)"<<std::endl;
            exit(1);
        }
        T min_tmp=vecPair[0].second;
        for(const auto& pair_tmp:vecPair){
            if(min_tmp>pair_tmp.second){
                min_tmp=pair_tmp.second;
            }
        }
        return min_tmp;
    }
}

#endif
