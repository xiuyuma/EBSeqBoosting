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
            
            COUNTS _var = aggregate::groupVar(counts, clusInfo);
            
            COUNTS var;
            
            size_t G = counts.rows();
            
            var = _var.rowwise().mean();
            
            COUNTS q(G,1), I(G,1), mn;
            
            I.fill(1);
            
            mn = counts.rowwise().mean();
            
            for(int i = 0; i < G; i++){
                
                if(abs(var(i,0) - 0) < 0.0001)
                    var(i,0) = 1;
                
                if(var(i,0) <= mn(i,0))
                    q(i,0) = 0.99;
                else
                    q(i,0) = mn(i,0) / var(i,0);
            }
            
            
            
        }
        
        
        Float kernel(COUNTS& _sum, std::vector<Float>& hyperParam)
        {
            return 0;
        }
        
        
    private:
        
        COUNTS _r;
        
    };
    
};
