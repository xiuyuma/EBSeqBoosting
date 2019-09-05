GetMultiFC=function(EBMultiOut,SmallNum=.01){
    if(ncol(EBMultiOut$Mean) < 3)
    stop("The input doesn't seem like an output from EBMultiTest")

    
    
	OutNames = rownames(EBMultiOut$PPMat)
    
	NumCondition = ncol(EBMultiOut$Mean)
    
	ConditionNames = colnames(EBMultiOut$AllParti)
    
	CondMeans = EBMultiOut$Mean
    
	colnames(CondMeans)=ConditionNames
    
	CondMeansPlus = CondMeans+SmallNum

	GeneRealMean = rowMeans(CondMeans)
    
	GeneR = EBMultiOut$RList
    
  GeneR[GeneR<=0 | is.na(GeneR)]=GeneRealMean[GeneR<=0 | is.na(GeneR)]*.99/.01

  GeneAlpha = EBMultiOut$Alpha
  
  GeneBeta = EBMultiOut$Beta
  

	FCMat=PostFCMat=matrix(0,ncol=choose(NumCondition,2),nrow=length(OutNames))
    
	rownames(FCMat) = rownames(PostFCMat) = OutNames
    
	k = 1
    
	ColNames=rep(NA,choose(NumCondition,2))
    
	for(i in 1:(NumCondition-1)){
		for(j in (i+1):NumCondition)
		{
		ColNames[k]=paste(ConditionNames[i],"Over",ConditionNames[j],sep="")
		FCMat[,k]=CondMeansPlus[,i]/CondMeansPlus[,j]


		nC1=sum(EBMultiOut$Condition==ConditionNames[i])
        
		nC2=sum(EBMultiOut$Condition==ConditionNames[j])
        
		GenePostAlphaC1 = GeneAlpha + nC1 * GeneR
        
		GenePostAlphaC2 = GeneAlpha + nC2 * GeneR
        
		GenePostBetaC1 = GeneBeta + nC1 * CondMeans[,i]
		GenePostBetaC2 = GeneBeta + nC2 * CondMeans[,j]
		GenePostQC1 = GenePostAlphaC1/(GenePostAlphaC1+GenePostBetaC1)
		GenePostQC2 = GenePostAlphaC2/(GenePostAlphaC2+GenePostBetaC2)

		GenePostFC=((1-GenePostQC1)/(1-GenePostQC2))*(GenePostQC2/GenePostQC1)
		PostFCMat[,k]= GenePostFC

		k=k+1
		}
	}
	colnames(FCMat)=colnames(PostFCMat)=ColNames
	Log2FCMat=log2(FCMat)
	Log2PostFCMat=log2(PostFCMat)
	Out=list(FCMat=FCMat,Log2FCMat=Log2FCMat,
					 PostFCMat=PostFCMat, Log2PostFCMat=Log2PostFCMat,
					 CondMeans=CondMeans, 
					 ConditionOrder=EBMultiOut$Condition)
}

