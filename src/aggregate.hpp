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
    
    static Eigen::MatrixXd groupMean(COUNTS& counts, CLUSINFO& clusInfo)
    {
        int K = clusInfo.size.size();
        
        Eigen::MatrixXd _sum = aggregate(counts, clusInfo);
        
        for(int i = 0; i < K; i++)
        {
            int _size = clusInfo.index[i].size();
            
            _sum.col(i) /= _size;
        }
        
        return _sum;
    }
    
    static Eigen::MatrixXd groupVar(COUNTS& counts, CLUSINFO& clusInfo)
    {
        int K = clusInfo.size.size();
        
        Eigen::MatrixXd _mean = groupMean(counts, clusInfo);
        
        Eigen::MatrixXd res(counts.rows(),K);
        
        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += (counts.col(s) - _mean.col(i)).unaryExpr(&[](double a){return a * a;});
            }
            
            res.col(i) /= clusInfo.index.[i].size();
        }
        
        return res;
    }
    
};
