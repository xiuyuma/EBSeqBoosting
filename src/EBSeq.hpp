#pragma once

#include "aggregate.hpp"
#include <assert.h>

namespace EBS
{
    
    class EBSeq
    {
        
    public:
        //on log scale
        virtual Float kernel(COUNTS& _sum, std::vector<Float>& hyperParam)
        {
            return 0;
        }
        
        virtual COUNTS kernelDerivative(COUNTS& _sum, COUNTS& hyperParam)
        {
            return COUNTS();
        }
        
        virtual void gradientAscent(COUNTS& _sum, COUNTS& hyperParam, std::vector<Float>& lrate){};
        
        
        
    protected:
        
        std::vector<Float> _hp;
        
        COUNTS _sum;
        
    };

};
