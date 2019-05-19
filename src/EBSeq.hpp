#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "helper.hpp"
#include "aggregate.hpp"

namespace EBS
{
    
    class EBSeq
    {

    public:
        virtual float kernel(vector<Float> geneCounts, vector<Float> hyperParam)
        {
            return 0;
        }
        
        
    };

};
