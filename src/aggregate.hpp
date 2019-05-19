#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "helper.hpp"
#include "float.hpp"
#include "counts.hpp"

namespace EBS
{
    

struct aggregate
{
    
    static COUNTS sum(COUNTS& counts, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS res(counts.rows(),K);
        
        res.fill(0);
        
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
        
        res.fill(0);

        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += (counts.col(s) - _mean.col(i)).unaryExpr(&square<Float>);
            }

            res.col(i) /= clusInfo.index[i].size();
        }

        return res;
    }
    
};

};
