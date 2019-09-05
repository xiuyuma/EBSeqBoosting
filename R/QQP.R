QQP <-
function(EBOut, GeneLevel=F){
    if(!"Alpha"%in%names(EBOut))stop("The input doesn't seem like an output from EBTest/EBMultiTest")
        AlphaResult=EBOut$Alpha
        BetaResult=EBOut$Beta
	    
	# Multi
		
    QList=EBOut$QList
    
    if(length(table(EBOut$Iso)) == 1)
    {
        for(j in 1:ncol(QList)){
            if(GeneLevel==F)Main=paste("Ig","C",j)
            if(GeneLevel==T)Main=paste("Gene","C",j)
            
            tmQ = QList[,j]
            
            tmpSize = length(tmQ[tmQ<1 & !is.na(tmQ)])
            
            rdpts = rbeta(tmpSize,AlphaResult,mean(BetaResult))
            
            qqplot(tmQ[tmQ<1],
            rdpts,xlab="estimated q's", ylab="simulated q's from fitted beta",
            main=Main,
            xlim=c(0,1),ylim=c(0,1))
            fit=lm(sort(rdpts)~sort(tmQ[tmQ<1  & !is.na(tmQ)]))
            abline(fit,col="red")
            
        }
        
    }else
    {
        
        Iso = EBOut$Iso
        I = max(as.numeric(Iso))
        
        for(i in 1:I)
        {
            tmp = which(Iso == i)
            
            for(j in 1:ncol(QList)){
                if(GeneLevel==F)Main=paste("Ig",i,"C",j)
                if(GeneLevel==T)Main=paste("Gene","C",j)
                
                tmQ = QList[tmp,j]
                
                tmpSize = length(tmQ[tmQ<1 & !is.na(tmQ)])
                
                rdpts = rbeta(tmpSize,AlphaResult,mean(BetaResult[tmp]))
                
                qqplot(tmQ[tmQ<1],
                rdpts,xlab="estimated q's", ylab="simulated q's from fitted beta",
                main=Main,
                xlim=c(0,1),ylim=c(0,1))
                fit=lm(sort(rdpts)~sort(tmQ[tmQ<1  & !is.na(tmQ)]))
                abline(fit,col="red")
                
                
                
            }
        }
    }

}

