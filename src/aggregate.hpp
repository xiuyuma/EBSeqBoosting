#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "helper.hpp"

struct aggregate
{
    using COUNTS = Eigen::MatrixXi;
    
    static COUNTS aggregate(COUNTS& counts, CLUSINFO& clusInfo)
    {
        int K = clusInfo.size.size();
        
        COUNTS res(counts.rows(),K);
        
        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += counts.col(s);
            }
        }
        
        return res;
    }
    
};
