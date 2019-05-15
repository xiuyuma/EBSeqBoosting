#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include "pat.hpp"

Rcpp::List pat(int K);
RcppExport SEXP pat(SEXP K)
{
    BEGIN_RCPP
    int k=as<int>(K);
    std::vector<std::vector<int> > res = partition::Part(k);
    Eigen::MatrixXi res_m(res.size(),k);
    for(int i=0;i<res.size();i++)
        for(int j=0;j<k;j++)
            res_m(i,j) = res[i][j];
    return List::create(Named("part") = res_m);
    END_RCPP
}
