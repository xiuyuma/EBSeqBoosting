#pragma once

#include "float.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <boost/math/special_functions/gamma.hpp>

namespace EBS
{
    
    class ALGO
    {
    public:
    	struct Node
		{
    		Float rs;
    		Float cs;
            Float distToNext;
            int sz;
    		std::vector<int> indexSet;
            Node *prev, *next;
		};
        
        template<typename ROW>
        static Node* createNode(ROW& csum, ROW& rsum,std::vector<Float>& logRatio, int pos, int size)
        {
            Node* res = new Node();
            
            assert(pos < rsum.size());
            
            res->rs = rsum(pos);
            res->cs = csum(pos);
            res->sz = size;
            if(pos == rsum.size() - 1){res->distToNext = 0;}
            else{res->distToNext = logRatio[pos];}
            res->indexSet.push_back(pos);
            res->prev = nullptr;
            res->next = nullptr;
            
            
            return res;
        }
        
        template<typename ROW>
        static Node* createNodeList(ROW& csum, ROW& rsum,std::vector<Float>& logRatio, int start, int end, std::vector<int>& sizes)
        {
            Node* head = createNode<ROW>(csum,rsum,logRatio,start,sizes[0]);
            Node* prev = head;
            for(int i = start + 1; i < end + 1; i++)
            {
                Node* tmp = createNode<ROW>(csum,rsum,logRatio,i,sizes[i - start]);
                prev->next = tmp;
                tmp->prev = prev;
                prev = tmp;
            }
            
            return head;
        }
    	
    	template<typename ROW>
        static void hclust(ROW& csum, ROW& rsum, std::vector<bool>& baseBit, 
                           std::vector<Float>& logRatio, int start, int end, Float alpha, Float beta, Float thre1, Float thre2, std::vector<int>& sizes)
		{
            std::list<Node> clus;
            auto head = createNodeList<ROW>(csum,rsum,logRatio,start,end,sizes);
            
            int counter = end - start;
            
            Float minDist = -INT_MAX;
            
            Node* minDistNode = nullptr;
            
            while(counter > 0)
            {
                auto tmpNode = head;
                for(size_t i = 0; i < counter; i++)
                {
                    if(tmpNode->distToNext > minDist)
                    {
                        minDist = tmpNode->distToNext;
                        minDistNode = tmpNode;
                    }
                    
                    tmpNode = tmpNode->next;
                }
                
                if(minDist > thre1 && minDistNode != nullptr)
                {
                    merge(minDistNode,alpha,beta,thre2);
                    counter--;
                }
                if(minDist < thre1)
                {
                    break;
                }
            }
            if(counter == 0)
            {
                return;
            }
            
            Node * tmpNode = head->next;
            
            while(tmpNode != nullptr)
            {
                baseBit[tmpNode->indexSet[0] - 1] = 1;
                tmpNode = tmpNode->next;
            }
            return;
		}
        
        //merge two nodes and delete right node
        static void merge(Node* left, Float alpha, Float beta, Float filter)
        {
            auto right = left->next;
            
            left->rs += right->rs;
            left->cs += right->cs;
            left->sz += right->sz;
            
            left->distToNext = kernel2(left->cs,right->next->cs,
                                       left->rs,right->next->rs,
                                       alpha,beta,left->sz,
                                       right->next->sz,filter);
            
            for(auto s:(right->indexSet))
            {
                left->indexSet.push_back(s);
            }
            
            left->next = right->next;
            
            delete right;
            
        }
        
        
        static inline Float kernel2(Float& cs1, Float& cs2, Float& rs1, Float& rs2, Float alpha, Float beta, int n1, int n2, Float filter)
        {
            // if too small mean, assume they are the same
            if(cs1 / n1 < filter && cs2 / n2 < filter )
            {
                return INT_MAX;
            }
            
            
            Float res = lbeta(alpha + rs1 + rs2, beta + cs1 + cs2) + lbeta(alpha, beta) - lbeta(alpha + rs1, beta + cs1) - lbeta(alpha + rs2, beta + cs2);
            
            return res;
        }
        
        static inline Float lbeta(Float x,Float y)
        {
            return boost::math::lgamma(x) + boost::math::lgamma(y) - boost::math::lgamma(x + y);
        }
        
    };
};
