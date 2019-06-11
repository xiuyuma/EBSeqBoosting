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
        
        virtual Float LogLikelihood()
        {
            return 0;
        }
        
    protected:
        //on log scale
        virtual COUNTS kernel(COUNTS& p)
        {
            return COUNTS();
        }
        
        virtual COUNTS kernelDerivative()
        {
            return COUNTS();
        }
        
    protected:
        
        COUNTS _sum;
        
        CLUSINFO _clusinfo;
        
        COUNTS _mean;
    };

};
