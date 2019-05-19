#pragma once

#include <vector>
#include <algorithm>

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

};


