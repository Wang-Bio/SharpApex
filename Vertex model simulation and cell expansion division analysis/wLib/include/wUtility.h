#ifndef _WUTILITY_H
#define _WUTILITY_H

#include <utility>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

namespace wVec{
    template <typename T>
    std::vector<T> v_av(std::vector<T> v,int NumAv){
        std::vector<double> output_v(v.size());

        int halfNeighbors = (NumAv-1)/2;
        for(int i=0;i<v.size();i++){
            double sum=0.0;
            for(int j=1-halfNeighbors-1;j<i+halfNeighbors;++j){
                if(j>=0&&j<v.size()){
                    sum+=v[j];
                }
                else if(j<0){
                    sum+=v[v.size()-1+j];
                }
                else if(j>v.size()){
                    sum+=v[j-v.size()];
                }
            }
            output_v[i]=sum/NumAv;
        }

        return output_v;
    }

    template <typename T>
    T v_min(std::vector<T> v){
        T v_min=v[0];
        for(auto vi:v){
            if(v_min>vi){
                v_min=vi;
            }
        }
        return v_min;
    }

    template <typename T>
    T v_max(std::vector<T> v){
        T v_max=v[0];
        for(auto vi:v){
            if(v_max<vi){
                v_max=vi;
            }
        }
        return v_max;
    }

    template <typename T>
    void v_print(std::vector<T> v){
        for(size_t vi=0;vi<v.size();vi++){
            std::cout<<vi<<" "<<v[vi]<<std::endl;
        }
    }

    template <typename T>
    T v_av(std::vector<T> v){
        T v_sum = v[0];
        for(size_t vi=1;vi<v.size();vi++){
            v_sum+=v[vi];
        }
        return v_sum/v.size();
    }

}

namespace wPair{
    template <typename T>
    double distance(const std::pair<T, T>& point1, const std::pair<T, T>& point2) {
        T xDiff = point2.first - point1.first;   // Difference in x-coordinates
        T yDiff = point2.second - point1.second; // Difference in y-coordinates

        return std::sqrt(xDiff * xDiff + yDiff * yDiff);
    }

    template <typename T>
    std::pair<T,T> addition(const std::pair<T,T>& point1, const std::pair<T,T>& point2){
        T xAdd = point1.first+point2.first;
        T yAdd = point1.second+point2.second;

        return std::make_pair(xAdd,yAdd);
    }

    template <typename T>
    std::pair<T,T> subtraction(const std::pair<T,T>& point1, const std::pair<T,T>& point2){
        T xSub = point1.first-point2.first;
        T ySub = point1.second-point2.second;

        return std::make_pair(xSub,ySub);
    }

    template <typename T>
    std::pair<T,T> multiplication(const std::pair<T,T>& point1, const T r){
        T xMul = point1.first*r;
        T yMul = point1.second*r;

        return std::make_pair(xMul,yMul);
    }
}


namespace wVecPair{

    //print vec_pair
    template <typename T>
    void v_print(std::vector<std::pair<T,T>> v){
        for(size_t vi=0;vi<v.size();vi++){
            std::cout<<vi<<" "<<v[vi].first<<" "<<v[vi].second<<std::endl;
        }
    }

    //read vec_pair from txt, with designated column
    template <typename T>
    std::vector<std::pair<T,T>> read_from_txt(const std::string& file,int pair_column1,int pair_column2){
        std::vector<std::pair<T,T>> pairs;
        std::ifstream inFile(file);
        std::string line;

        if(!inFile.is_open()){
            std::cerr << "Error opening file: "<<file<<std::endl;
            exit(1);
        }

        while(getline(inFile,line)){
            std::istringstream iss(line);
            std::string token;
            T value;
            std::vector<T> values;
            while (getline(iss, token, '\t')){
                std::istringstream converter(token);
                if(converter >> value){
                    values.push_back(value);
                }
            }

            if (pair_column1 < values.size() && pair_column2 < values.size()) {
                pairs.emplace_back(values[pair_column1], values[pair_column2]);
            }
        }

        inFile.close();
        return pairs;
    }

    //calculate the maximum of the second element of vec_pairs
    template <typename T>
    T max_second(const std::vector<std::pair<T,T>> & vec_pairs){
        T maxT = vec_pairs[0].second; 
        for(auto vec_pair : vec_pairs){
            if(maxT < vec_pair.second){
                maxT = vec_pair.second;
            }
        }
        return maxT;
    }

    //calculate the minimum of the second element of vec_pairs
    template <typename T>
    T min_second(const std::vector<std::pair<T,T>> & vec_pairs){
        T minT = vec_pairs[0].second; 
        for(auto vec_pair : vec_pairs){
            if(minT > vec_pair.second){
                minT = vec_pair.second;
            }
        }
        return minT;
    }

    //calculate the normalized vec_pair
    template <typename T>
    std::vector<std::pair<T,T>> normalization(const std::vector<std::pair<T,T>>& vec_pairs){
        std::vector<std::pair<T,T>> post_vec_pairs;
        T maxT = max_second(vec_pairs);
        T minT = min_second(vec_pairs);
        T normT = maxT-minT;
            //std::cout<<"maxT "<<maxT<<" minT "<<minT<<" normT "<<normT<<std::endl;
        for(auto vec_pair:vec_pairs){
            std::pair<T,T> post_vec_pair = wPair::multiplication(vec_pair,1.0/normT);
            post_vec_pairs.push_back(post_vec_pair);
        }

        return post_vec_pairs;
    }


}



#endif