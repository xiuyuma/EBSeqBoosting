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
    static COUNTS sum(Eigen::VectorXd& vec, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS res(1,K);
        
        res.fill(0);
        
        for(auto i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res(0,i) += vec(s);
            }
        }
        
        return res;
    }
    
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
    
    static COUNTS sum(COUNTS& counts, CLUSINFO& clusInfo, Eigen::VectorXd& sizeFactor)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS res(counts.rows(),K);
        
        res.fill(0);
        
        for(auto i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += counts.col(s) / sizeFactor(s);
            }
        }
        
        return res;
    }
    
    static COUNTS groupMean(COUNTS SUM, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        for(auto i = 0; i < K; i++)
        {
            size_t _size = clusInfo.index[i].size();
            
            SUM.col(i) /= _size;
        }
        
        return SUM;
    }
    
    template <typename VAL>
    static VAL square(VAL x)
    {
        return x * x;
    }
    
    static COUNTS groupVar(COUNTS& counts, COUNTS& MEAN, CLUSINFO& clusInfo)
    {
        size_t K = clusInfo.size.size();
        
        COUNTS res(counts.rows(),K);
        
        res.fill(0);
        
        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += (counts.col(s) - MEAN.col(i)).unaryExpr(&square<Float>);
            }
            
            res.col(i) /= clusInfo.index[i].size();
        }
        
        return res;
    }
    
    static COUNTS groupVar(COUNTS& counts, COUNTS& MEAN, CLUSINFO& clusInfo, Eigen::VectorXd& sizeFactor)
    {
        size_t K = clusInfo.size.size();

        COUNTS res(counts.rows(),K);
        
        res.fill(0);

        for(int i = 0; i < K; i++)
        {
            for(auto s:clusInfo.index[i])
            {
                res.col(i) += (counts.col(s) - MEAN.col(i) * sizeFactor(s)).unaryExpr(&square<Float>) / sizeFactor(s);
            }

            res.col(i) /= clusInfo.index[i].size();
        }

        return res;
    }
    
};

};
