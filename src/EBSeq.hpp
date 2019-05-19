#pragma once

#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "helper.hpp"

class EBSeq
{

public:
    virtual float kernel(vector<float> geneCounts, vector<float> hyperParam)
    {
        return 0;
    }
    
    
    
};
