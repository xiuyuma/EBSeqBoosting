# Iterative pruning for empirical bayes

### Speed up EBSeq (http://bioconductor.org/packages/release/bioc/html/EBSeq.html) under comparision of more groups

EBSeq provides an empirical bayesian way to handle means comparison of multiple groups. It struggled when the number of groups(K) become big, specifically, the time and space complexity of the algorithm increased exponentially with K. 

We developed an iterative information sharing scheme to efficiently pruning the space of differential mean(DM) patterns and select those patterns with higher marginal densities. 

The algorithm can be summaried in the following:
1) At each gene, find the pattern maximizing the prior predictive function (PT*) and construct a neighbour covering PT* based on local bayesian factors

2) Take the union of each gene specific neighbour of patterns as the pool of patterns to be considered in the EBSeq scheme, using EM algorithm to find out hyper parameters and marginal density of each patterns

3) Remove patterns with small marginal densities, repeat step 2 until convergence of the EM algorithm.


p.s.
It may be helpful using tbb for parallel computing at the first pruning. However the benefits of speeding up is very limited and asking for tbb dependency on server would create more overhead.

