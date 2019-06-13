#pragma once

#include "aggregate.hpp"
#include <assert.h>
#include "partition.hpp"
#include "helper.hpp"

namespace EBS
{
    
    Float wrapperFunc(const std::vector<Float> &param, std::vector<Float> &grad, void *object)
    {
        
        std::cout << "called\n";
        
        Eigen::VectorXd p(param.size() + 1);
        
        Float last = 1;
        
        for(size_t i = 0; i < param.size(); i++)
        {
            p(i) = param[i];
            
            last -= p(i);
            
        }
        
        p(param.size()) = last;
        
        if(last < 0)
            return  -INT_MAX;
        
        return static_cast<NB*>(object)->OBJ(p);
        
        
    }
    
//    Float myconstraint(const std::vector<Float> &param, std::vector<Float> &grad, void *data)
//    {
//        Float res = 0;
//        for(auto x:param)
//            res += x;
//
//        return res - 1;
//    }

};
