#pragma once

#include "aggregate.hpp"
#include <assert.h>
#include "partition.hpp"
#include "helper.hpp"

namespace EBS
{
    
    Float wrapperFunc(const std::vector<Float> &param, std::vector<Float> &grad, void *object)
    {
        Float alpha = param[0];
        
        Eigen::VectorXd beta(param.size() - 1);
        
        for(size_t i = 1; i < param.size(); i++)
            beta(i - 1) = param[i];
        
        return static_cast<NB*>(object)->OBJ(alpha,beta);
    }

};
