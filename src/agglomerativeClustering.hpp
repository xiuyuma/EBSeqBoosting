#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Float.hpp"

namespace EBS
{
    
    class ALGO
    {
    	struct Node
		{
    		Float rs;
    		Float cs;
    		
    		std::vector<int> SET;
		};
    	
    	template<typename ROW>
        static void hclust(ROW& csum, ROW& rsum, std::vector<bool>& baseBit, 
        		std::vector<Float>& logRatio, int start, int end, Float alpha, Float beta)
		{
        	
		}
    }
};
