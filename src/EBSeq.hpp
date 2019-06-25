#pragma once

#include "aggregate.hpp"
#include <assert.h>
#include "partition.hpp"
#include "helper.hpp"

namespace EBS
{
    
    class EBSeq
    {
        
    public:
        virtual void gradientAscent(){};
         
        EBSeq(COUNTS& scRNAexpMatrix, std::vector<int>& cellCluster)
        {
            _clusinfo = helper::clusInfo(cellCluster);
            
            _sum = aggregate::sum(scRNAexpMatrix, _clusinfo);
            
            _mean = aggregate::groupMean(scRNAexpMatrix, _clusinfo);
        }
        
        
        
        COUNTS getSUM()
        {
            return _sum;
        }
        
        COUNTS getMEAN()
        {
            return _mean;
        }
        
        // get the posterior probability
        virtual COUNTS getPOSP() = 0;
        
        void EM(size_t max_iteration, Float changeThre)
        {
            size_t Iter = 0;
            
            Float changeRatio = 10;
            
            Float lastOBJ;
            
            bool first = true;
            
            while(Iter < max_iteration && changeRatio > changeThre)
            {
                Estep();
                
                Mstep();
                
                if(first)
                {
                    lastOBJ = getOBJ();
                    
                    Iter++;
                    
                    first = false;
                    
                    continue;
                }
                
                changeRatio = (getOBJ() - lastOBJ) / lastOBJ;
                
                changeRatio = abs(changeRatio);
                
                Iter++;
            }
        }
        
        
    protected:
        
        // matrix for prior predictive scores under different DE patterns
        virtual void kernel() = 0;
        
        // function to find more probable DE patterns
        virtual void DEpat() = 0;
        
        // E - step of EM
        virtual void Estep() = 0;
        
        // M - step of EM
        virtual void Mstep() = 0;
        
        // get value of objective function
        virtual Float getOBJ()=0;
        
        
    protected:
        
        COUNTS _sum;
        
        CLUSINFO _clusinfo;
        
        COUNTS _mean;
        
    };

};
