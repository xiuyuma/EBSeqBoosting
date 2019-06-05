#pragma once


#include "EBSeq.hpp"

namespace EBS
{
  
    class NB:public EBSeq
    {
    public:
        
        
        NB(COUNTS& scRNAexpMatrix, std::vector<int>& cellCluster) : EBSeq(scRNAexpMatrix, cellCluster)
        {
            COUNTS _var = aggregate::groupVar(scRNAexpMatrix, _clusinfo);
            
            COUNTS var;
            
            size_t G = scRNAexpMatrix.rows();
            
            var = _var.rowwise().mean();
            
            COUNTS q(G,1), I(G,1), mn;
            
            I.fill(1);
            
            mn = scRNAexpMatrix.rowwise().mean();
            
            for(int i = 0; i < G; i++){
                
                if(abs(var(i,0) - 0) < 0.0001)
                    var(i,0) = 1;
                
                if(var(i,0) <= mn(i,0))
                    q(i,0) = 0.99;
                else
                    q(i,0) = mn(i,0) / var(i,0);
            }
            
            _r = (mn.cwiseProduct(q)).array() / (I - q).array();
        }
        
        void init(std::vector<Float> hyperParam, std::vector<Float>& lrate, int UC)
        {
            assert(UC < _sum.cols());
            
            _hp = hyperParam;
            
            _lrate = lrate;
            
            _uncertainty = UC;
        }
        
        // only to be called in init
        void DEpat()
        {
            size_t G = _mean.rows();
            
            for(size_t i = 0; i < G; i++)
            {
                _order.push_back(helper::sortIndexes<ROW>(_mean.row(i)));
            }
        }
        
        Float kernel(std::vector<int>& pat)
        {
            return 0;
            
            
        }
        
        void gradientAscent()
        {
            
        }
        
        
        
        
        
    private:
        
        typedef decltype(_mean.row(0)) ROW;
        
        COUNTS _r;
        
        std::vector<Float> _lrate;
        
        // prop of each nonzero pattern
        Eigen::VectorXd _p;
        
        std::unordered_map<int, std::vector<int>> _hash;
        
        int _uncertainty;
        
        std::vector<std::vector<size_t>> _order;
        
    };
    
};
