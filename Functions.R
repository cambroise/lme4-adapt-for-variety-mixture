##############################################
#### Function for formating entries
##############################################

GetGenotype <- function(DataList,GenoColumns){
    unique(unlist(lapply(1:length(DataList), function(nn){
        unique(as.vector(DataList[[nn]][,GenoColumns[[nn]]]))
    })))    
}

GetCrossedGenotype <- function(DataList,GenoColumns){
    CG <- unique(unlist(lapply(1:length(DataList), function(nn){
        if (length(GenoColumns[[nn]])>1){
            Comb <- combn(x = GenoColumns[[nn]],simplify = T,m=2)
            CrossList <- lapply(1:ncol(Comb), function(u){
                unique(paste(DataList[[nn]][,Comb[1,u]],DataList[[nn]][,Comb[2,u]],sep='_'))  
            })     
        } else {CrossList=NULL}
    })))
    G <- GetGenotype(DataList,GenoColumns)
    CG <- c(CG,paste(G,G,sep='_'))
}


GetPhenotype <- function(DataList,Plist){
    lapply(Plist, function(pp){
        Reduce('c',lapply(DataList, function(d){ d[[pp]] }))
    })
}

 

BuildX <- function(DataList,ListX){
    if (is.null(ListX)){
        X <- rep(1, sum(sapply(DataList, nrow)))
    } else {
        Var <- Reduce('rbind',lapply(Data, function(d) d[,colnames(d)%in%ListX,drop=F]))
        Fake <- data.frame(Y=rnorm(nrow(Var)),Var)
        X <- model.matrix(lm(Y~.,Fake))
    }
}



BuildZ_Genotype <- function(DataList,G,Gcol){
    Reduce('rbind',lapply(1:length(DataList), function(nn){
        MixtSize <- length(Gcol[[nn]])
        Reduce('+',lapply(1:MixtSize, function(mm){
            sapply(G, function(gg){
                Data[[nn]][,Gcol[[nn]][mm]]==gg
            })
        }))/MixtSize
    }))
}



BuildZ_CrossedGenotype <- function(DataList,CG,Gcol){
    Reduce('rbind',lapply(1:length(DataList), function(nn){
        MixtSize <- length(Gcol[[nn]])
        if (MixtSize>1){
            #Crossed terms (i != j)
            Comb <- combn(x = Gcol[[nn]],simplify = T,m=2)
            tmp <- Reduce('+',lapply(1:ncol(Comb), function(cc){
                sapply(CG, function(gg){
                    paste(DataList[[nn]][,Comb[1,cc]],DataList[[nn]][,Comb[2,cc]],sep='_')==gg
                })
            }))*(2/MixtSize**2)
            
            #Crossed terms (i =j)
            tmp <- tmp + Reduce('+',lapply(1:MixtSize, function(mm){
                sapply(CG, function(gg){
                    paste(DataList[[nn]][,Gcol[[nn]][mm]], DataList[[nn]][,Gcol[[nn]][mm]],sep='_')==gg
                })
            }))/(MixtSize**2)
        } else {
            #No mixture
            tmp <- sapply(CG, function(gg){
                paste(DataList[[nn]][,Gcol[[nn]]], DataList[[nn]][,Gcol[[nn]]],sep='_')==gg
            })
        }
    }))
    
}