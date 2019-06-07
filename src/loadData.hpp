#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace EBS
{
    
    // very coarse function, for test
    // need number of rows and columns and input are file of integers
    
    COUNTS readData(std::string path, char delim, int G, int n)
    {
        std::ifstream File;
        
        File.open(path);
        
        std::string element;
        
        std::string line;
        
        COUNTS res(G,n);
        
        size_t i = 0,j = 0;
        
        if(File.is_open())
        {
            
            
            while(std::getline(File,line))
            {
                
                std::stringstream tmp(line);
                
                while(std::getline(tmp,element,delim))
                {
                    res(i,j) = std::stoi(element);
                    
                    j++;
                    
                    j %= n;
                    
                    if(j == 0){i++;}
                }
                
                
            }
            
        }
        
        return res;
        
    }
    
};
