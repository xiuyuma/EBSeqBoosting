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

Rcpp::List EBSeq(Rcpp::NumericMatrix scExpMatrix, Rcpp::IntegerVector groupLabel, Rcpp::IntegerVector isoLabel, Rcpp::NumericVector sizeFactor, int iter, double alpha, Rcpp::NumericVector beta, double step1, double step2, int uc, double thre, double sthre, double filter, double stopthre);

RcppExport SEXP EBSeq(SEXP scExpMatrix, SEXP groupLabel, SEXP isoLabel, SEXP sizeFactor, SEXP iter, SEXP alpha, SEXP beta, SEXP step1, SEXP step2, SEXP uc, SEXP thre, SEXP sthre, SEXP filter, SEXP stopthre)
{
    // param scExpMatrix: scRNA seq transcripts matrix (normalized counts required)
    // param groupLabel: group label for each cell
    // param isoLabel: for isoform case need to share beta within same isoform label
    // param sizeFactor: normalizing factor for raw counts (for old EBSeq), 1 for normalized counts
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
    
    IntegerVector iL(isoLabel);
    
    NumericVector szf(sizeFactor);
    
    NumericVector bta(beta);
    
    const int ng = scExpM.rows();
    const int nc = scExpM.cols();
    
    EBS::COUNTS data(ng,nc);
    std::copy(scExpM.begin(),scExpM.end(),data.data());
    
    std::vector<int> conditions(nc);
    std::copy(cluster.begin(),cluster.end(),conditions.begin());
    
    std::vector<int> iLabel(ng);
    std::copy(iL.begin(),iL.end(),iLabel.begin());
    
    Eigen::VectorXd sf(nc);
    std::copy(szf.begin(),szf.end(),sf.data());
    
    Eigen::VectorXd bt(ng);
    std::copy(bta.begin(),bta.end(),bt.data());
    
    std::vector<EBS::Float> lrate;
    lrate.push_back(stepsizeAlp);
    lrate.push_back(stepsizeBt);
    
    // create and initialize NB class object
    EBS::NB X = EBS::NB(data,conditions,sf);
    
    X.init(alp, bt, iLabel, lrate, UC, threshold, sthreshold, filterThre);
    
    // EM
    X.EM(itr, stopThre);
    
    // results to be returned
    // posterior prob
    auto POSP = X.getPOST();
    
    // DE patterns to be considered
    auto DEP = X.getDEP();
    
    // q to be returned
    auto QQ = X.getQ();
    
    // mean to be returned
    auto mm = X.getMEAN();
    
    // var to be returned
    auto var = X.getVar();
    
    // convert to R acceptable object
    Eigen::MatrixXi mDep(DEP.size(),DEP[0].size());
    
    for (size_t ri = 0; ri < mDep.rows(); ri++)
        for(size_t ci = 0; ci < mDep.cols(); ci++)
            mDep(ri,ci) = DEP[ri][ci];
    
    return Rcpp::List::create(Named("DEpattern") = mDep, Named("Posterior") = POSP, Named("Alpha") = X.getALP(), Named("Beta") = X.getBETA(), Named("q") = QQ, Named("mean") = mm, Named("var") = var);
    
    
    
}
    
    
    
// function to be called when there is no replicate for each condition
Rcpp::List EBSeqWQ(Rcpp::NumericMatrix scExpMatrix, Rcpp::IntegerVector groupLabel, Rcpp::NumericVector r, Rcpp::IntegerVector isoLabel, Rcpp::NumericVector sizeFactor, int iter, double alpha, Rcpp::NumericVector beta, double step1, double step2, int uc, double thre, double sthre, double filter, double stopthre);

RcppExport SEXP EBSeqWQ(SEXP scExpMatrix, SEXP groupLabel, SEXP r, SEXP isoLabel, SEXP sizeFactor, SEXP iter, SEXP alpha, SEXP beta, SEXP step1, SEXP step2, SEXP uc, SEXP thre, SEXP sthre, SEXP filter, SEXP stopthre)
{
    // param scExpMatrix: scRNA seq transcripts matrix (normalized counts required)
    // param groupLabel: group label for each cell
    // param r: parameter r of NB
    // param isoLabel: for isoform case need to share beta within same isoform label
    // param sizeFactor: normalizing factor for raw counts (for old EBSeq), 1 for normalized counts
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
    
    NumericVector R(r);
    
    IntegerVector iL(isoLabel);
    
    NumericVector szf(sizeFactor);
    
    NumericVector bta(beta);
    
    const int ng = scExpM.rows();
    const int nc = scExpM.cols();
    
    EBS::COUNTS data(ng,nc);
    std::copy(scExpM.begin(),scExpM.end(),data.data());
    
    std::vector<int> conditions(nc);
    std::copy(cluster.begin(),cluster.end(),conditions.begin());
    
    Eigen::VectorXd rR(ng);
    std::copy(R.begin(),R.end(),rR.data());
    
    std::vector<int> iLabel(ng);
    std::copy(iL.begin(),iL.end(),iLabel.begin());
    
    Eigen::VectorXd sf(nc);
    std::copy(szf.begin(),szf.end(),sf.data());
    
    Eigen::VectorXd bt(ng);
    std::copy(bta.begin(),bta.end(),bt.data());
    
    std::vector<EBS::Float> lrate;
    lrate.push_back(stepsizeAlp);
    lrate.push_back(stepsizeBt);
    
    // create and initialize NB class object
    EBS::NB X = EBS::NB(data,conditions,sf,rR);
    
    X.init(alp, bt, iLabel, lrate, UC, threshold, sthreshold, filterThre);
    
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
    
   
    
}
    

