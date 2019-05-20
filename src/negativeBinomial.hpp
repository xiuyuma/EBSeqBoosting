#pragma once


#include "EBSeq.hpp"

namespace EBS
{
  
    class NB:public EBSeq
    {
    public:
        
        void init(COUNTS& counts, CLUSINFO& clusInfo, std::vector<Float> hyperParam)
        {
            _sum = aggregate::sum(counts, clusInfo);
            
            _hp = hyperParam;
        }
        
        
        Float kernel(COUNTS& _sum, std::vector<Float>& hyperParam)
        {
            
        }
        
        
    private:
        
        COUNTS _r;
        
    };
    
};
