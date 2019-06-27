#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include "partition.hpp"
#include "helper.hpp"
#include "aggregate.hpp"
#include "EBSeq.hpp"
#include "float.hpp"
#include "counts.hpp"
#include "negativeBinomial.hpp"
#include "loadData.hpp"
#include "wrapper.hpp"

using namespace Rcpp;
using namespace Eigen;
using namespace std;

Rcpp::List EBSeq(Rcpp::NumericMatrix scExpMatrix, Rcpp::IntegerVector groupLabel, int iter, double alpha, Rcpp::NumericVector beta, double step1, double step2, int uc, double thre, double sthre, double filter, double stopthre);

RcppExport SEXP EBSeq(SEXP scExpMatrix, SEXP groupLabel, SEXP iter, SEXP alpha, SEXP beta, SEXP step1, SEXP step2, SEXP uc, SEXP thre, SEXP sthre, SEXP filter, SEXP stopthre) {
    // param scExpMatrix: scRNA seq transcripts matrix (normalized counts required)
    // param groupLabel: group label for each cell
    // param iter: number of max iteration in EM
    // param alpha: start point of hyper parameter alpha
    // param beta: start point of hyper parameter beta
    // param step1: stepsize for gradient ascent of alpha
    // param step2: stepsize for gradietn ascent of beta
    // param uc: number of unceratin relations between means of subtypes for each gene level
    // param thre: threshold for determining whether a relation is sure or uncertain
    // param sthre: shrinkage threshold for iterative pruning space of DE patterns
    // param filter: filterthreshold for low expression gene for DE analysis
    // param stopthre: stopping threshold for EM
    
    // return: a list containing considered DE patterns and their posterior probability
    // values for alpha and beta
    
    
    
    
    
    // Rcpp wrapper
    BEGIN_RCPP
    
    // config S pointer to c++ data structure
    int itr = as<int>(iter);
    int UC = as<int>(uc);
    
    double alp = as<double>(alpha);
    double stepsizeAlp = as<double>(step1);
    double stepsizeBt = as<double>(step2);
    double threshold = as<double>(thre);
    double sthreshold = as<double>(sthre);
    double filterThre = as<double>(filter);
    double stopThre = as<double>(stopthre);
    
    NumericMatrix scExpM(scExpMatrix);
    
    IntegerVector cluster(groupLabel);
    
    NumericVector bta(beta);
    
    const int ng = scExpM.rows();
    const int nc = scExpM.cols();
    
    EBS::COUNTS data(ng,nc);
    
    std::copy(scExpM.begin(),scExpM.end(),data.data());
    
    std::vector<int> conditions(nc);
    
    std::copy(cluster.begin(),cluster.end(),conditions.begin());
    
    Eigen::VectorXd bt(ng);
    
    std::copy(bta.begin(),bta.end(),bt.data());
    
    std::vector<EBS::Float> lrate;
    
    lrate.push_back(stepsizeAlp);
    
    lrate.push_back(stepsizeBt);
    
    // create and initialize NB class object
    EBS::NB X = EBS::NB(data,conditions);
    
    X.init(alp, bt, lrate, UC, threshold, sthreshold, filterThre);
    
    // EM
    X.EM(itr, stopThre);
    
    // results to be returned
    // posterior prob
    auto POSP = X.getPOST();
    
    // DE patterns to be considered
    auto DEP = X.getDEP();
    
    // convert to R acceptable object
    Eigen::MatrixXi mDep(DEP.size(),DEP[0].size());
    
    for (size_t ri = 0; ri < mDep.rows(); ri++)
        for(size_t ci = 0; ci < mDep.cols(); ci++)
            mDep(ri,ci) = DEP[ri][ci];
    
    return Rcpp::List::create(Named("DEpattern") = mDep, Named("Posterior") = POSP, Named("Alpha") = X.getALP(), Named("Beta") = X.getBETA());
    
    END_RCPP
    
}
    
    
    
    
    
    

