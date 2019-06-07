#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <fstream>

namespace EBS
{
    struct loadInfo
    {
        COUNTS data;
        
        vector<int> clus;
    }
    
    loadInfo readData(string path)
    {
        std::ifstream File;
        
        File.open(path);
        
        
    }
    
};
