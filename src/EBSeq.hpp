#pragma once

#include "aggregate.hpp"
#include <assert.h>

namespace EBS
{
    
    class EBSeq
    {
        
    public:
        //on log scale
        virtual Float kernel()
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
