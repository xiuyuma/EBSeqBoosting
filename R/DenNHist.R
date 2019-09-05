DenNHist <-
function(EBOut, GeneLevel=F)
{	
	if(!"Alpha"%in%names(EBOut))stop("The input doesn't seem like an output from EBTest/EBMultiTest")
	
	Alpha=EBOut$Alpha
	Beta=EBOut$Beta
    
	# Multi
	
    QList=EBOut$QList
    
    if(length(EBOut$Iso) == 1)
    {
        for(j in 1:ncol(QList)){
            if(GeneLevel==F)Main=paste("Ig","C",j)
            if(GeneLevel==T)Main=paste("Gene","C",j)
            
            tmQ = QList[,j]
            
        hist(tmQ[tmQ < .98 & tmQ > 0],
            prob=T,col="blue",breaks=100,
            main=Main,
            xlim=c(0,1),xlab=paste("Q alpha=",round(Alpha,2),
            " beta=",round(mean(Beta),2),sep=""))
            
            tmpseq = seq(0.001,1,length=1000)
            
            ll = tmpseq
            
            lines(ll,dbeta(ll,Alpha,mean(Beta)),col="green",lwd=2)
            
            legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2)
        }
        
    }else
    {
        
        Iso = EBOut$Iso
        I = max(Iso)
        
        for(i in 1:I)
        {
            tmp = which(Iso == i)
            
            for(j in 1:ncol(QList)){
                if(GeneLevel==F)Main=paste("Ig",i,"C",j)
                if(GeneLevel==T)Main=paste("Gene","C",j)
                
                tmQ = QList[tmp,j]
                hist(tmQ[tmQ < .98 & tmQ > 0],
                prob=T,col="blue",breaks=100,
                main=Main,
                xlim=c(0,1),xlab=paste("Q alpha=",round(Alpha,2),
                " beta=",round(mean(Beta),2),sep=""))
            
            tmpseq = seq(0.001,1,length=1000)
            
            ll = tmpseq
            
            lines(ll,dbeta(ll,Alpha,mean(Beta[tmp])),col="green",lwd=2)
            
            legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2)
            
            
            }
        }
    }
    
    
	


	
	
}

