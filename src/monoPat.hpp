#pragma once

#include<vector>


std::vector<std::vector<int> > monoPart(const int& n)
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
