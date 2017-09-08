rm(list=ls())
library(nnet)
source('./Functions.R')
source('./lme4-non-binary-Z.R')

##############################################
#### Import data
##############################################
## Data
#Data file list
FileList <- c('melanges','pures')
#Phenotype list
PhenoList <- c('Rendement','NbPltesM','EpiMai','NbEpisM','Rendement','NbEpPlte')
#Fixed explanatory variable list
FixedVarList <- c('Rep')

##Import
Data <- lapply(FileList, function(name){
  tmp <- paste0("Data/",name,".csv")
  d <- read.csv (tmp, sep =";", as.is = T, header = TRUE, dec=",")
  d$Rep <- as.character(d$Rep)
  return(d)
})
NbData <- length(FileList)
NbObs <- sum(sapply(Data,nrow))
##Check
lapply(Data,dim)
lapply(Data,head)


##############################################
#### Get genotype values
##############################################
##Find which columns contain the genoyptic information is in each dataset
GenoCol <- lapply(Data, function(d) grep(pattern = 'Genotype',x = colnames(d)))
##Get the list of genotypes
Genotype <- GetGenotype(Data,GenoCol)
##Get the list of crossed genotypes
CrossedGenotype <- GetCrossedGenotype(Data,GenoCol)

##############################################
#### Get phenotypes
##############################################
Phenotype <- GetPhenotype(Data,PhenoList)

##############################################
#### Build X, Zg and Zc
##############################################
##Incidence matrix for fixed effects
X <- BuildX(Data,FixedVarList)
##Incidence matrix for random genotype effects
Zg <- BuildZ_Genotype(Data,Genotype,GenoCol)
##Incidence matrix for random crossed genotype effects
Zc <- BuildZ_CrossedGenotype(Data,CrossedGenotype,GenoCol)
##Check
apply(Zg,1,sum)
apply(Zc,1,sum)
table(as.vector(Zg))
table(as.vector(Zc))

##############################################
#### Perform analysis for all phenotypes
##############################################
ListZ <- list(GMA=Zg,SMA=Zc)

## Estimation of two different model (with and without GMA)
mymod<-mylmer(Response =  Phenotype[[1]],X, ListZ,REML=FALSE) 
mymod.gma<-mylmer(Response =  Phenotype[[1]],X, ListZ[1],REML=FALSE) 

## Classical outputs
summary(mymod)
sigma(mymod)

library(lattice)
dotplot(myranef(mymod))


pr01 <- profile(mymod)
xyplot(pr01, aspect = 1.3)

confint(pr01,level=0.9)
splom(pr01)  

## Testing between models
anova(mymod.gma,mymod)
