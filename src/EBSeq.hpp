#pragma once

#include "aggregate.hpp"
#include <assert.h>

namespace EBS
{
    
    class EBSeq
    {
        
    public:
        //on log scale
        virtual Float kernel(COUNTS& sum, std::vector<Float>& hyperParam)
        {
            return 0;
        }
        
        virtual COUNTS kernelDerivative(COUNTS& sum, COUNTS& hyperParam)
        {
            return COUNTS();
        }
        
        virtual void gradientAscent(COUNTS& sum, COUNTS& hyperParam, std::vector<Float>& lrate){};
        
        
        
    protected:
        
        std::vector<Float> _hp;
        
        COUNTS _sum;
        
    };

};
