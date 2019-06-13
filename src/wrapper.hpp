#pragma once

#include "aggregate.hpp"
#include <assert.h>
#include "partition.hpp"
#include "helper.hpp"

namespace EBS
{
    
    Float wrapperFunc(const std::vector<Float> &param, std::vector<Float> &grad, void *object)
    {
        Eigen::VectorXd p(param.size());
        
        for(size_t i = 0; i < param.size(); i++)
            p(i) = param[i];
        
        p = p / p.sum();
        
        Float res;

        res = static_cast<NB*>(object)->OBJ(p);
        
        
        std::cout << "function value " << res << "\n";
        
        return res;
        
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
