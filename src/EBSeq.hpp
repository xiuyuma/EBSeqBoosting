#pragma once

#include "aggregate.hpp"

namespace EBS
{
    
    class EBSeq
    {

    public:
        virtual Float kernel(vector<Float> geneCounts, vector<Float> hyperParam)
        {
            return 0;
        }
        
        virtual vector<Float> derivative(COUNTS counts, vector<Float> hyperParam)
        {
            return vector<Float>();
        }
    };

};
