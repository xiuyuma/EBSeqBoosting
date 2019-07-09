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
         
        EBSeq(COUNTS& scRNAexpMatrix, std::vector<int>& cellCluster, Eigen::VectorXd& sizeFactor)
        {
            _clusinfo = helper::clusInfo(cellCluster);
            
            _sum = aggregate::sum(scRNAexpMatrix, _clusinfo, sizeFactor);
            
            _mean = aggregate::groupMean(_sum, _clusinfo);
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
        virtual COUNTS getPOST() = 0;
        
        void EM(size_t max_iteration, Float changeThre)
        {
            size_t Iter = 0;
            
            Float changeRatio = 10;
            
            Float lastOBJ;
            
            bool first = true;
            
            while(Iter < max_iteration && changeRatio > changeThre)
            {
                
                if(first)
                {
                    
                    Estep();
                    
                    Mstep();
                    
                    lastOBJ = getOBJ();
                    
                    Iter++;
                    
                    first = false;
                    
                    std::cout << "OBJ " << lastOBJ << "\n";
                    
                    std::cout << "size " << DEPsize() << "\n";
                    
                    continue;
                }
                
                shrinkage();
                
                Estep();
                
                Mstep();
                
                changeRatio = (getOBJ() - lastOBJ) / lastOBJ;
                
                changeRatio = abs(changeRatio);
                
                lastOBJ = getOBJ();
                
                std::cout << "OBJ " << lastOBJ << "\n";
                
                std::cout << "size " << DEPsize() << "\n";
                
                Iter++;
            }
        }
        
        
    protected:
        
        // for one step EM
        virtual void gradientAscent(){};
        
        // matrix for prior predictive scores under different DE patterns
        virtual void kernel() = 0;
        
        // function to find more probable DE patterns
        virtual void DEpat() = 0;
        
        // E - step of EM
        virtual void Estep() = 0;
        
        // M - step of EM
        virtual void Mstep() = 0;
        
        // shrinkage
        virtual void shrinkage() = 0;
        
        // get value of objective function
        virtual Float getOBJ() = 0;
        
        // get number of DE patterns
        virtual size_t DEPsize() = 0;
        
        
        
    protected:
        
        COUNTS _sum;
        
        CLUSINFO _clusinfo;
        
        COUNTS _mean;
        
    };

};
