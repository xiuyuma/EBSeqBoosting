#pragma once

#include <vector>
#include <algorithm>

namespace EBS
{
    struct partition
    {
        static std::vector<std::vector<int> > Part(int n)
        {
            
            std::vector<std::vector<int> > start(1);
            start[0].push_back(1);
            if(n == 1){
                return start;
            }
            
            for(int i = 1;i < n; ++i){
                std::vector<std::vector<int> > new_p;
                size_t L = start.size();
                for(int j = 0;j < L; ++j){
                    int M = *std::max_element(start[j].begin(),start[j].end());
                    for(int k = 0;k < M; ++k){
                        start[j].push_back(k+1);
                        new_p.push_back(start[j]);
                        start[j].pop_back();
                    }
                    start[j].push_back(M+1);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                }
                start = new_p;
            }
            return start;
        }
        
        static std::vector<std::vector<int> > monoPart(int n)
        {
            std::vector<std::vector<int> > start(1);
            start[0].push_back(1);
            if(n == 1){
                return start;
            }
            
            for(int i = 1;i < n; ++i){
                std::vector<std::vector<int> > new_p;
                size_t L = start.size();
                for(int j = 0;j < L; ++j){
                    int M = *std::max_element(start[j].begin(),start[j].end());
                    //push last position again
                    start[j].push_back(start[j][i - 1]);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                    
                    start[j].push_back(M+1);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                }
                start = new_p;
            }
            return start;
        }
        
        // K - 1 "<" or "=", 1 for "<" and 0 for "="
        template <size_t K>
        static std::bitset<K-1> mapToBit(vector<int>& part)
        {
            static_asser(part.size() == K);
            
            std::bitset<K-1> res;
            
            // all one
            for(size_t i = 0; i < K - 1; i++)
            {
                res[i] = 0;
                
                if(part[i] != part[i + 1]){res[i] = 1;}
                
            }
            
            return res;
            
        }
    };

};


