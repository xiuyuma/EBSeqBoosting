#pragma once

#include <Eigen/Dense>

struct helper
{
    static std::vector<int> which(std::vector<int> clus, int i)
    {
        if(clus.size < 1 || std::find(clus.begin(),clus.end(),i) == clus.end())
            return std::vector<int>;
        
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
    
    static bool clusCheck(vector<int> clus)
    {
        
    }
}

Eigen::Matrixi aggr(Eigen::Matrixi counts, vector<int> clus)
{
    
}
