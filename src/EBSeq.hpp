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
        
        void gradientAscent(COUNTS& _sum, COUNTS& hyperParam, std::vector<Float>& lrate)
        {
            size_t _n = hyperParam.size();
            
            assert(lrate.size() == _n);
            
            COUNTS drv = kernelDerivative(_sum, hyperParam);
            
            for(size_t i = 0; i < _n; i++)
            {
                
            }
            
            
        };
        
        Float likelihood()
        {
            size_t G = _sum.rows();
            
            Float LL = 0;
            
            for(size_t i = 0; i < G; i++)
            {
                LL += kernel(_sum, _hp);
            }
            
            return LL;
        };
        
    protected:
        
        std::vector<Float> _hp;
        
        COUNTS _sum;
    
        
    };

};
