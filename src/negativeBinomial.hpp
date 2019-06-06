#pragma once


#include "EBSeq.hpp"
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>

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
        
        void init(std::vector<Float> hyperParam, std::vector<Float> lrate, int UC)
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
            
            size_t K = _sum.cols();
            
            std::vector<Float> abslogRatio(K);
            
            std::vector<int> baseClus(K);
            
            for(size_t i = 0; i < G; i++)
            {
                auto ord = helper::sortIndexes<ROW>(_mean.row(i));
                
                _order.push_back(ord);
                
                baseClus[0] = 1;
                
                for(size_t j = 1; j < K; j++)
                {
                    Float s1 = _sum(i,ord[j - 1]);
                    
                    Float s2 = _sum(i,ord[j]);
                    
                    Float r1 = _clusinfo.size[ord[j - 1]] * _r(i,0);
                    
                    Float r2 = _clusinfo.size[ord[j]] * _r(i,0);
                    
                    Float tmp = kernel2case(s1,s2,r1,r2);
                    
                    abslogRatio.push_back(abs(tmp));
                    
                    //  more favorable for equal mean
                    if(tmp > 0)
                    {
                        baseClus[j] = baseClus[j - 1];
                    }
                    else
                    {
                        // DE start a new cluster
                        baseClus[j] = baseClus[j - 1] + 1;
                    }
                    
                }
                
                auto tmpOrd = helper::sortIndexes<std::vector<Float>>(abslogRatio);
                
                std::cout << "G " << i << "\n";
                
                std::cout << baseClus[0] << "," << baseClus[1] << "," << baseClus[2] << "\n";
                
                std::cout << tmpOrd[0] << "," << tmpOrd[1] << "," << tmpOrd[2] << "\n";
            }
            
            
            
            
        }
        
        Float kernel(std::vector<int>& pat)
        {
            return 0;
            
            
        }
        
        void gradientAscent()
        {
            
        }
        
        
        Float kernel2case(Float& s1, Float& s2, Float& r1, Float& r2)
        {
            Float alpha = _hp[0];
            
            Float beta = _hp[1];
            
            Float res = lbeta(alpha + r1 + r2, beta + s1 + s2) + lbeta(alpha, beta) - lbeta(alpha + r1, beta + s1) - lbeta(alpha + r2, beta + s2);
            
            return res;
        }
        
        inline Float lbeta(Float x,Float y)
        {
            return boost::math::lgamma(x) + boost::math::lgamma(x) - boost::math::lgamma(x + y);
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
