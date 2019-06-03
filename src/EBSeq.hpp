#pragma once

#include "aggregate.hpp"
#include <assert.h>

namespace EBS
{
    
    class EBSeq
    {
        
    public:
        //on log scale
        virtual Float kernel(std::vector<int>& pat)
        {
            return 0;
        }
        
        virtual COUNTS kernelDerivative()
        {
            return COUNTS();
        }
        
        virtual void gradientAscent(){};
         
        
        
    protected:
        
        std::vector<Float> _hp;
        
        COUNTS _sum;
        
    };

};
