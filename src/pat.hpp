#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>

std::vector<std::vector<int> > partition(const int& n){

    std::vector<std::vector<int> > start(1);
    start[0].push_back(1);
    if(n==1){
        return start;
    }
    
    for(int i=1;i<n;++i){
        std::vector<std::vector<int> > new_p;
        size_t L = start.size();
        for(int j = 0;j < L; ++j){
            int M = *max_element(start[j].begin(),start[j].end());
            for(int k = 0;k < M; ++k){
                start[j].push_back(k+1);
                new_p.push_back(start[j]);
                start[j].pop_back();
            }
            start[j].push_back(M+1);
            new_p.push_back(start[j]);
            start[j].pop_back();
        }
        start=new_p;
    }
    return start;
}

Rcpp::List pat(int K);
RcppExport SEXP pat(SEXP K){
    BEGIN_RCPP
    int k=as<int>(K);
    std::vector<std::vector<int> > res=partition(k);
    Eigen::MatrixXi res_m(res.size(),k);
    for(int i=0;i<res.size();i++)
        for(int j=0;j<k;j++)
            res_m(i,j) = res[i][j];
    return List::create(Named("part") = res_m);
    END_RCPP
}
