#pragma once


#include "EBSeq.hpp"
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <set>

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
            
            for(size_t i = 0; i < G; i++){
                
                if(abs(var(i,0) - 0) < 0.0001)
                    var(i,0) = 1;
                
                if(var(i,0) <= mn(i,0))
                    q(i,0) = 0.99;
                else
                    q(i,0) = mn(i,0) / var(i,0);
            }
            
            _r = (mn.cwiseProduct(q)).array() / (I - q).array();
            
            auto tmp = _clusinfo.size;
            
            _csize.resize(1,tmp.size());
            
            for(size_t i = 0; i < tmp.size(); i++)
            {
                _csize(0,i) = tmp[i];
            }
        }
        
        
        void init(Float alpha, Eigen::VectorXd beta, std::vector<Float> lrate, int UC, Float thre, Float filter)
        {
            assert(UC < _sum.cols());
            
            _alpha = alpha;
            
            _beta = beta;
            
            _lrate = lrate;
            
            _uncertainty = UC;
            
            _threshold = thre;
            
            _filter = filter;
            
            DEpat();
            
            size_t n = _dep.size();
            
            assert(n > 0);
            
            _p.resize(n);
            
            _p.fill(1.0 / n);
            
            kernel();
        }
        
        
        size_t DEPsize()
        {
            return _dep.size();
        }
        
        std::vector<std::vector<int>> getDEP()
        {
            return _dep;
        }
        
        std::vector<size_t> getGUC()
        {
            return _guc;
        }
        
        std::vector<COUNTS> getPAT()
        {
            return _pat;
        }
        
        
        inline Float lbeta(Float x,Float y)
        {
            return boost::math::lgamma(x) + boost::math::lgamma(y) - boost::math::lgamma(x + y);
        }
        
        inline COUNTS lbeta(COUNTS& A, COUNTS& B)
        {
            return A.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) + B.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) - (A + B).unaryExpr<Float(*)(Float)>(&boost::math::lgamma);
        }
    
        void setAlphaBeta(Float alpha, Eigen::VectorXd& beta)
        {
            _alpha = alpha;
            
            _beta = beta;
        }
        
        void setP(Eigen::VectorXd& p)
        {
            _p = p;
        }
        
        
        void kernel()
        {
            _kernel.resize(_sum.rows(),_pat.size());
            
            for(size_t i = 0; i < _pat.size(); i++)
            {
                COUNTS _csum = _sum * _pat[i];
                
                COUNTS _rsum = _r * _csize * _pat[i];
                
                COUNTS A = (_rsum.array() + _alpha).matrix();
                
                COUNTS B = _csum.colwise() + _beta;
                
                COUNTS res = lbeta(A,B);
                
                res = (res.array() - boost::math::lgamma(_alpha)).matrix();
                
                res =  res.colwise() - (_beta.unaryExpr<Float(*)(Float)>(&boost::math::lgamma) + (_alpha + _beta.array()).matrix().unaryExpr<Float(*)(Float)>(&boost::math::lgamma));
                
                _kernel.col(i) = res.rowwise().sum();
                
            }
            
        }
        
        
        
        
        Float OBJ(Eigen::VectorXd& p)
        {
            setP(p);
            
            return (_kernel * _p).sum();
        }
        
        COUNTS posterior()
        {
            
            kernel();
            
            assert(abs(_p.sum() - 1) < 0.0001);
            
            Eigen::VectorXd M = _kernel.rowwise().maxCoeff();
            
            auto posp = _kernel.colwise() - M;
            
            posp.unaryExpr<Float(*)(Float)>(& exp);
            
            Eigen::VectorXd total = posp * _p;
            
            total = (1 / total.array()).matrix();
            
            //outer product of total and p
            COUNTS div = total * _p.transpose();
            
            return (posp.array() * div.array()).matrix();
        }
        
    private:
        // only to be called in init
        void DEpat()
        {
            size_t G = _mean.rows();
            
            size_t K = _sum.cols();
            
            std::vector<Float> abslogRatio(K - 1);
            
            std::vector<int> baseClus(K);
            
            // DE patterns to be considered
            std::set<std::string> dep;
            
            for(size_t i = 0; i < G; i++)
            {
                auto ord = helper::sortIndexes<ROW>(_mean.row(i));
                
                auto ord2 = helper::sortIndexes2<ROW>(_mean.row(i));
                
                baseClus[0] = 1;
                
                for(size_t j = 1; j < K; j++)
                {
                    Float s1 = _sum(i,ord[j - 1]);
                    
                    Float s2 = _sum(i,ord[j]);
                    
                    int n1 = _clusinfo.size[ord[j - 1]];
                    
                    int n2 = _clusinfo.size[ord[j]];
                    
                    Float r1 = n1 * _r(i);
                    
                    Float r2 = n2 * _r(i);
                    
                    Float tmp = kernel2case(s1,s2,r1,r2,n1,n2,i);
                    
                    //                    if(tmp < 1)
                    //                    {
                    //                        std::cout << "i " << i <<"," << s1 << "," << s2 << "," << r1 << "," << r2 << "\n";
                    //                    }
                    
                    abslogRatio[j - 1] = abs(tmp);
                    
                    
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
                
                auto baseBit = partition::mapToBit(baseClus);
                
                int localUC = 0;
                
                while(abslogRatio[localUC] < _threshold)
                {
                    localUC++;
                }
                
                _guc.push_back(localUC);
                
                localUC = std::min(localUC , _uncertainty);
                
                auto pBit = partition::genBit(localUC);
                
                // get promising DE pattern
                for(auto x:pBit)
                {
                    
                    auto newBit = baseBit;
                    
                    for(int t = 0; t < _uncertainty; t++)
                    {
                        auto tmpJ = tmpOrd[t];
                        
                        newBit[tmpJ] = x[t];
                    }
                    
                    auto newClus = partition::bitToPart(newBit);
                    
                    auto newClusOr = newClus;
                    
                    for(size_t iter = 0; iter < ord2.size(); iter++)
                    {
                        newClusOr[iter] = newClus[ord2[iter]];
                    }
                    
                    auto newClusOrd = partition::reorder(newClusOr);
                    
                    auto sClus = partition::toString<std::vector<int>>(newClusOrd);
                    
                    if(dep.find(sClus) == dep.end())
                    {
                        dep.insert(sClus);
                        
                        _dep.push_back(newClusOrd);
                        
                        _pat.push_back(partition::toMatrix(newClusOrd));
                    }
                    
                }
                
            }
            
        }
        
        
        
        Float kernel2case(Float& s1, Float& s2, Float& r1, Float& r2, int n1, int n2, size_t i)
        {
            // if too small mean, assume they are the same
            if(s1 / n1 < _filter && s2 / n2 <_filter )
            {
                return INT_MAX;
            }
            
            Float alpha = _alpha;
            
            Float beta = _beta(i);
            
            Float res = lbeta(alpha + r1 + r2, beta + s1 + s2) + lbeta(alpha, beta) - lbeta(alpha + r1, beta + s1) - lbeta(alpha + r2, beta + s2);
            
            return res;
        }
        
        
        void gradientAscent()
        {
            
        }
        
        
    private:
        
        typedef decltype(_mean.row(0)) ROW;
        
        // hyper parameter r can be estimated by MM
        COUNTS _r;
        
        // alpha for the beta prior shared by genome, estimated by EM
        Float _alpha;
        
        // beta for the beta prior specified by each gene, estimated by EM
        Eigen::VectorXd _beta;
        
        COUNTS _csize;
        
        std::vector<Float> _lrate;
        
        // prop of each nonzero pattern
        Eigen::VectorXd _p;
        
        // upper bound of unsure "<" and "=", controlling number of patterns DE
        int _uncertainty;
        
        // at least one group mean should be no smaller than this value to do DE comparison
        Float _filter;
        
        // positve threshold to decide how many uncertain patterns
        Float _threshold;
        
        // gene level uncertainty
        std::vector<size_t> _guc;
        
        // matrix for DE pattern
        std::vector<COUNTS> _pat;
        
        // vector for DE pattern
        std::vector<std::vector<int>> _dep;
        
        // kernel matrix
        COUNTS _kernel;
        
        
        
    };
    
};
