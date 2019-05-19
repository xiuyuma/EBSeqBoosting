#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "helper.hpp"

struct aggregate
{
    using COUNTS = Eigen::MatrixXf;
    
    static COUNTS sum(COUNTS& counts, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS res(counts.rows(),K);
        
        for(auto i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += counts.col(s);
            }
        }
        
        return res;
    }
    
    static COUNTS groupMean(COUNTS& counts, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS _sum = sum(counts, clusInfo);
        
        for(auto i = 0; i < K; i++)
        {
            size_t _size = clusInfo.index[i].size();
            
            _sum.col(i) /= _size;
        }
        
        return _sum;
    }
    
    template <typename VAL>
    static VAL square(VAL x)
    {
        return x * x;
    }
    
    static COUNTS groupVar(COUNTS& counts, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS _mean = groupMean(counts, clusInfo);
        
        COUNTS res(counts.rows(),K);
        
        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += (counts.col(s) - _mean.col(i)).unaryExpr(&square<float>);
            }
            
            res.col(i) /= clusInfo.index[i].size();
        }
        
        return res;
    }
    
};
