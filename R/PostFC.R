PostFC=function(EBoutput, SmallNum=.01) {
	if(ncol(EBoutput$Mean) != 2)
		stop("The input doesn't seem like an output from EBTest")
	GeneRealMeanC1= EBoutput$Mean[,1]
	GeneRealMeanC2= EBoutput$Mean[,2]
	GeneRealMeanC1Plus=GeneRealMeanC1+SmallNum
	GeneRealMeanC2Plus=GeneRealMeanC2+SmallNum
	GeneRealMean=(GeneRealMeanC1+GeneRealMeanC2)/2

	GeneRealFC=GeneRealMeanC1Plus/GeneRealMeanC2Plus

	GeneR=unlist(EBoutput$RList)
	GeneR[GeneR<=0 | is.na(GeneR)]=GeneRealMean[GeneR<=0 | is.na(GeneR)]*.99/.01

	GeneAlpha = EBoutput$Alpha
    
	GeneBeta = EBoutput$Beta
	  
		# Post alpha P_a_C1= alpha + r_C1 * n_C1
	  # Post beta P_b_C1= beta + Mean_C1 * n_C1
	  # P_q_C1= P_a_C1/ (P_a_C1 + P_b_C1)
	  # Post FC = ((1-P_q_C1)/P_q_c1) /( (1-P_q_c2)/P_q_c2)

	nC1 = sum(EBoutput$Conditions==levels(EBoutput$Conditions)[1])
    
	nC2 = sum(EBoutput$Conditions==levels(EBoutput$Conditions)[2])
    
	GenePostAlphaC1 = GeneAlpha + nC1 * GeneR
    
	GenePostAlphaC2 = GeneAlpha + nC2 * GeneR
    
	GenePostBetaC1 = GeneBeta + nC1 * GeneRealMeanC1
    
	GenePostBetaC2 = GeneBeta + nC2 * GeneRealMeanC2
    
	GenePostQC1 = GenePostAlphaC1 / (GenePostAlphaC1 + GenePostBetaC1)
    
	GenePostQC2 = GenePostAlphaC2 / (GenePostAlphaC2 + GenePostBetaC2)

	GenePostFC = ((1-GenePostQC1)/(1-GenePostQC2))*(GenePostQC2/GenePostQC1)
    
	Out=list(PostFC=GenePostFC, RealFC=GeneRealFC,
					 Direction=paste(levels(EBoutput$Conditions)[1],"Over", levels(EBoutput$Conditions)[2])
					 )

}
