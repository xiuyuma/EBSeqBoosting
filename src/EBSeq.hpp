#pragma once

#include "aggregate.hpp"
#include <assert.h>

namespace EBS
{
    
    class EBSeq
    {

    public:
        //on log scale
        virtual Float kernel(VEC geneCounts, std::vector<Float> hyperParam)
        {
            return 0;
        }
        
        virtual VEC kernelDerivative(COUNTS counts, std::vector<Float> hyperParam)
        {
            return VEC();
        }
        
        void gradientAscent(COUNTS counts, std::vector<Float> hyperParam, std::vector<Float> lrate)
        {
            size_t _n = hyperParam.size();
            
            assert(lrate.size() == _n);
            
            VEC drv = kernelDerivative(counts, hyperParam);
            
            for(size_t i = 0; i < _n; i++)
            {
                
            }
            
            
        };
        
        Float likelihood(COUNTS counts, std::vector<Float> hyperParam)
        {
            size_t G = counts.rows();
            
            Float LL;
            
            for(size_t i = 0; i < G; i++)
            {
                LL += kernel(counts.row(i), _hp);
            }
            
            return LL;
        };
        
    protected:
        
        std::vector<Float> _hp;
        
    
        
    };

};
