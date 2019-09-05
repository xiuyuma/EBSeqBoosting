EBTest <-
function(Data,NgVector=NULL,Conditions, sizeFactors, uc, Alpha=NULL, Beta=NULL, Qtrm=1, QtrmCut=0
    ,maxround = 50, step1 = 1e-6,step2 = 0.01, thre = log(2), sthre = 0.001, filter = 10, stopthre = 1e-3)
{
	
    ## validity check
    expect_is(sizeFactors, c("numeric","integer"))
	expect_is(maxround,  c("numeric","integer"))
	if(!is.factor(Conditions))Conditions=as.factor(Conditions)
	if(is.null(rownames(Data)))stop("Please add gene/isoform names to the data matrix")

	if(!is.matrix(Data))stop("The input Data is not a matrix")
	if(length(Conditions)!=ncol(Data))stop("The number of conditions is not the same as the number of samples! ")
	if(nlevels(Conditions)>2)stop("More than 2 conditions! Please use EBMultiTest() function")
	if(nlevels(Conditions)<2)stop("Less than 2 conditions - Please check your input")
	if(length(sizeFactors)!=ncol(Data))
		stop("The number of library size factors is not the same as the number of samples!")		
	
	
    

	#Normalized and filtered
	DataNorm=GetNormalizedMat(Data, sizeFactors)
	expect_is(DataNorm, "matrix")
	Levels=levels(as.factor(Conditions))

	QuantileFor0=apply(DataNorm,1,function(i)quantile(i,Qtrm))
	AllZeroNames=which(QuantileFor0<=QtrmCut)
	NotAllZeroNames=which(QuantileFor0>QtrmCut)
    
	if(length(AllZeroNames)>0)
					    cat(paste0("Removing transcripts with ",Qtrm*100,
							    " th quantile < = ",QtrmCut," \n",
									length(NotAllZeroNames)," transcripts will be tested\n"))
                                    
	if(length(NotAllZeroNames)==0)stop("0 transcript passed")
    
	Data=Data[NotAllZeroNames,]
    
    if(!is.null(NgVector))
    {
        if(length(NgVector) != nrow(DataNorm))
        {
            stop("NgVector should have same size as number of genes")
        }
        NgVector = NgVector[NotAllZeroNames]
        NgVector = as.factor(NgVector)
        levels(NgVector) = 1:length(levels(NgVector))
    }
    
    
    if(is.null(NgVector)){NgVector = 1}

    if(is.null(Alpha))
    {
        Alpha = 0.4
    }
    
    if(is.null(Beta))
    {
        Beta = 0
    }
    
    MeanList=rowMeans(DataNorm)
    
    VarList=apply(DataNorm, 1, var)
    
    cd = Conditions
    
    levels(cd) = 1:length(levels(cd))
    
  
    res = EBSeqTest(Data,cd,uc, iLabel = NgVector,sizefactor = sizeFactors,
    iter = maxround,alpha = Alpha, beta = Beta, step1 = step1,step2 = step2,
    thre = thre, sthre = sthre, filter = filter, stopthre = stopthre)
    
    if(nrow(res$DEpattern) != 2){
        stop("too few DE patterns, try reducing sthre, increasing thre")
    }
    
    if(res$DEpattern[1,1] == res$DEpattern[1,2])
    {
        Allequal = 1
        Alldiff = 2
    }
    else
    {
        Allequal = 2
        Alldiff = 1
    }
    
    Mat = res$Posterior
    
    Matwith0 = matrix(NA,nrow = nrow(DataNorm), ncol = 2)

    rownames(Matwith0) = rownames(DataNorm)
    
    Matwith0[NotAllZeroNames,] = Mat
    
    
    rownames(Mat) = rownames(Data)
    colnames(Mat) = c(1,2)
    colnames(Mat)[Allequal] = "PPEE"
    colnames(Mat)[Alldiff] = "PPDE"
    
    parti = res$DEpattern
    
    colnames(parti) = levels(Conditions)
    
    Result=list(Alpha=res$Alpha,Beta=res$Beta,P=res$prop,
    RList=res$r, MeanList=MeanList,
    VarList=VarList, QList = res$q,
    Mean = res$mean,Var = res$var, PoolVar=res$poolVar,
    DataNorm = DataNorm, Iso = NgVector,
    AllZeroIndex = AllZeroNames,
    PPMat = Mat,AllParti = parti, PPMatWith0 = Matwith0,
    Conditions=Conditions)

}

