#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <numeric>

namespace EBS
{
    
struct CLUSINFO
{
    std::vector<std::vector<int> > index;
    std::vector<int> size;
};

struct helper
{
    
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
    
    
    static CLUSINFO clusInfo(std::vector<int>& clus)
    {
        int K = *std::max_element(clus.begin(),clus.end());
        
        std::vector<std::vector<int> > index;
        
        std::vector<int> size;
        
        for(int i = 1; i < K + 1; i++)
        {
            index.push_back(which(clus,i));
            
            size.push_back(index[i - 1].size());
        }
        
        CLUSINFO res;
        
        res.index = index;
        
        res.size = size;
        
        return res;
    }
    
    // expected type: row of eigen matrix
    template<typename T>
    static std::vector<size_t> sortIndexes(T ROW)
    {
        size_t K = ROW.size();
        
        std::vector<size_t> tmp(K);
        
        std::iota(tmp.begin(),tmp.end(),0);
        
        //res = tmp;
        
        std::sort(tmp.begin(),tmp.end(),
                  [&ROW](size_t i1, size_t i2) {return ROW[i1] < ROW[i2];});
        
        //std::sort(res.begin(),res.end(),
        //          [&tmp](size_t i1, size_t i2){return tmp[i1] < tmp[i2];});
        
        return tmp;
    }
    
    template<typename T>
    static std::vector<size_t> sortIndexes2(T ROW)
    {
        size_t K = ROW.size();
        
        std::vector<size_t> tmp(K),res;
        
        std::iota(tmp.begin(),tmp.end(),0);
        
        res = tmp;
        
        std::sort(tmp.begin(),tmp.end(),
                  [&ROW](size_t i1, size_t i2) {return ROW[i1] < ROW[i2];});
        
        std::sort(res.begin(),res.end(),
                  [&tmp](size_t i1, size_t i2){return tmp[i1] < tmp[i2];});
        
        return res;
    }
    
};

};

