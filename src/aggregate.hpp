#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>

struct helper
{
    static std::vector<int> which(std::vector<int> clus, int i)
    {
        if(clus.size() < 1 || std::find(clus.begin(),clus.end(),i) == clus.end())
            return std::vector<int>();
        
        std::vector<int> res;
        
        for(size_t t = 0; t < clus.size(); t++)
        {
            if(clus[t] == i)
            {
                res.push_back(t);
            }
        }
        
        return res;
    }
    
    // clus should be array of 1 - K
    static bool clusCheck(std::vector<int> clus)
    {
        int smallest = *std::min_element(clus.begin(),clus.end());
        
        if(smallest != 1)
            return false;
        
        int largest = *std::max_element(clus.begin(),clus.end());
        
        for(int t = 1; t < largest + 1; t++)
        {
            if(std::find(clus.begin(),clus.end(),t) == clus.end())
                return false;
        }
        
        return true;
    }
};

Eigen::Matrixi aggr(Eigen::Matrixi counts, vector<int> clus)
{
    
}
