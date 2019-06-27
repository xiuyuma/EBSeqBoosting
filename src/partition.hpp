#pragma once

#include <vector>
#include <algorithm>
#include <unordered_map>
#include "counts.hpp"
#include <unordered_set>
#include <string>

namespace EBS
{
    struct partition
    {
        static std::vector<std::vector<int> > Part(int n)
        {
            
            std::vector<std::vector<int> > start(1);
            start[0].push_back(1);
            if(n == 1){
                return start;
            }
            
            for(int i = 1;i < n; ++i){
                std::vector<std::vector<int> > new_p;
                size_t L = start.size();
                for(int j = 0;j < L; ++j){
                    int M = *std::max_element(start[j].begin(),start[j].end());
                    for(int k = 0;k < M; ++k){
                        start[j].push_back(k+1);
                        new_p.push_back(start[j]);
                        start[j].pop_back();
                    }
                    start[j].push_back(M+1);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                }
                start = new_p;
            }
            return start;
        }
        
        static std::vector<std::vector<int> > monoPart(int n)
        {
            std::vector<std::vector<int> > start(1);
            start[0].push_back(1);
            if(n == 1){
                return start;
            }
            
            for(int i = 1;i < n; ++i){
                std::vector<std::vector<int> > new_p;
                size_t L = start.size();
                for(int j = 0;j < L; ++j){
                    int M = *std::max_element(start[j].begin(),start[j].end());
                    //push last position again
                    start[j].push_back(start[j][i - 1]);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                    
                    start[j].push_back(M+1);
                    new_p.push_back(start[j]);
                    start[j].pop_back();
                }
                start = new_p;
            }
            return start;
        }
        
        // K - 1 "<" or "=", 1 for "<" and 0 for "="
        static std::vector<bool> mapToBit(std::vector<int>& part)
        {
            size_t K = part.size();
            
            std::vector<bool> res(K - 1, 0);
            
            for(size_t i = 0; i < K - 1; i++)
            {
                if(part[i] != part[i + 1]){res[i] = 1;}
            }
            
            return res;
            
        }
        
        static std::vector<int> bitToPart(std::vector<bool>& bits)
        {
            size_t K = bits.size();
            
            std::vector<int> res(K + 1);
            
            res[0] = 1;
            
            for(size_t i = 0; i < K; i++)
            {
                if(bits[i]){res[i + 1] = res[i] + 1;}
                else{res[i + 1] = res[i];}
            }
            
            return res;
        }
        
        static std::vector<std::vector<bool>> genBit(int n)
        {
            std::vector<std::vector<bool>> res;
            
            if(n == 0)
            {
                return res;
            }
            
            res.push_back(std::vector<bool>(1,true));
            
            res.push_back(std::vector<bool>(1,false));
            
            for(int i = 1; i < n; i++)
            {
                std::vector<std::vector<bool>> tmp;
                
                for(auto x:res)
                {
                    
                    auto y = x;
                    
                    x.push_back(true);
                    
                    tmp.push_back(x);
                    
                    y.push_back(false);
                
                    tmp.push_back(y);
                }
                
                res = tmp;
            }
            
            return res;
        }
        
        static std::vector<int> reorder(std::vector<int>& part)
        {
            size_t K = part.size();
            
            std::vector<int> res(K);
            
            std::unordered_set<int> seen;
            
            size_t i = 0;
            
            int label = 1;
            
            while(seen.size() < K)
            {
                int tmpLabel = part[i];
                
                for(int j = 0; j < K; j++)
                {
                    if(part[j] == tmpLabel)
                    {
                        res[j] = label;
                        
                        seen.insert(j);
                    }
                }
                
                label++;
                
                for(size_t t = i; t < K; t++)
                {
                    if(seen.find(t) == seen.end())
                    {
                        i = t;
                        
                        break;
                    }
                }
                
            }
            
            return res;
        }
        
        template <typename T>
        static std::string toString(T& part)
        {
            std::string res = "";
            
            for(size_t i = 0; i < part.size(); i++)
            {
                res += std::to_string(part[i]);
            }
            
            return  res;
        }
        
        // convert string representation of partition to matrix
        static COUNTS toMatrix (std::vector<int>& part)
        {
            int submax = *std::max_element(part.begin(),part.end());
            
            size_t K = part.size();

            COUNTS res(K,submax);
            
            res.fill(0);
            
            for(int i = 0; i < submax; i++)
            {
                for(int j = 0; j < K; j++)
                {
                    if(part[j] == i + 1)
                    {
                        res(j,i) = 1;
                    }
                }
            }
            
            return res;
        }
        
        
        static std::unordered_map<std::vector<bool>, size_t> buildHash(std::vector<std::vector<int> >& parts)
        {
            
            std::unordered_map<std::vector<bool>, size_t> res;
            
            for(size_t i = 0; i < parts.size(); i++)
            {
                auto bit = mapToBit(parts[i]);
                
                res[bit] = i;
            }
            
            return res;
        }
        
        
        
        
    };

};


