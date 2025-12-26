
################################################################################
# Step1:General tasks
################################################################################
rm(list=ls())
getwd()
dir.create("C:/Users/Win10/Documents/GWAS2GP"); dir.create("C:/Users/Win10/Documents/GWAS2GP/data")   #Change user name
setwd("C:/Users/Win10/Documents/GWAS2GP")

download.file(sprintf("https://tassel.bitbucket.io/docs/TASSELTutorialData5.zip"), "data/TASSELTutorialData5.zip") #to unzip
zip_file <- "data/TASSELTutorialData5.zip"
unzip(zip_file, exdir = "data/unzipd")

pheno <- read.table("data/unzipd/TASSELTutorialData5/mdp_traits.txt",head=TRUE)

library(rMVP)
MVP.Data(fileHMP="data/unzipd/TASSELTutorialData5/mdp_genotype.hmp.txt", #Genotypic data in HapMap format
         filePhe="data/unzipd/TASSELTutorialData5/mdp_traits.txt",
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=T,
         out="mvp.hmp")
################################################################################
# Step2: import formatted data (MVP.Data) from working directory 
################################################################################

genotypic_dat <- attach.big.matrix("mvp.hmp.geno.desc")
phenotype_dat <- read.table("mvp.hmp.phe",head=TRUE)
map_info <- read.table("mvp.hmp.geno.map" , head = TRUE)
#popstr <- read.table("mdp_population_structure.txt", header = T, skip = 1)
Kinship <- MVP.K.VanRaden(genotypic_dat, verbose = T)

################################################################################
# GWAS run 
################################################################################
################################################################################
# PCA and Kinship as covariates
################################################################################

dir.create("C:/Users/Win10/Documents/GWAS2GP/PcaKinship")
setwd("C:/Users/Win10/Documents/GWAS2GP/PcaKinship")

GWAS_PCA_Kin_mvp <- MVP(
  phe = phenotype_dat, 
  geno = genotypic_dat, 
  map = map_info,
  method =  c("GLM", "MLM", "FarmCPU"),
  nPC.GLM = 5,
  nPC.MLM = 5, 
  nPC.FarmCPU = 5,
  K = Kinship
)

################################################################################
# hmp file to 0 1 2 format #####################################################
################################################################################

setwd("C:/Users/Win10/Documents/GWAS2GP")
getwd()
library(rMVP)
library(bigmemory)
library(data.table)

geno_numeric_list <- list(
  
  # Big matrix genotype (0/1/2, additive)
  geno = attach.big.matrix("mvp.hmp.geno.desc"),
  
  # SNP map
  map = fread("mvp.hmp.geno.map", data.table = FALSE),
  
  # Individual IDs
  ind = fread("mvp.hmp.geno.ind", header = FALSE, data.table = FALSE)
)


geno <- {
  # library(bigmemory)
  # library(data.table)
  
  # Accessing objects from the geno_numeric_list
  bm  <- geno_numeric_list$geno
  ind <- geno_numeric_list$ind
  map <- geno_numeric_list$map
  
  # Converting the big matrix into a regular matrix
  m   <- as.matrix(bm)
  
  # Assigning dimension names based on the Individual IDs and SNP Map
  dimnames(m) <- list(ind[1:nrow(m), 1], map[1:ncol(m), 1])
  
  m
}

# get genotype row names (To sync geno pheno
geno_ids <- trimws(rownames(geno))

# get phenotype IDs (first column)
pheno_ids <- trimws(pheno[[1]])

# keep only matching rows in pheno
pheno_new <- pheno[pheno_ids %in% geno_ids, ]
geno_new <- geno[geno_ids %in% pheno_ids, ]

################################################################################
### Genomic Prediction
################################################################################

library(BGLR) #import package

Y <- pheno_new 
X <- geno_new
# A<-wheat.A #pedigree relatioship matrix from BGLR
y<-Y[,2] # EarHT

#Setting the linear predictor
ETA<-list(list(X=X, model='BRR')) #Gaussian prior
ETA<-list(list(X=X, model='BL')) #Double exponential


#Fitting the model
fm<-BGLR(y=y,ETA=ETA, nIter=5000, burnIn=1000, thin = 5) #1000 interactions in total
cor(fm$y,fm$yHat) #calculate correlation
plot(fm) #plot graph pred x observed

# Extracting results from the model
bHat <- fm$ETA[[1]]$b # Estimated Marker Effects
SD.bHat <- fm$ETA[[1]]$SD.b
plot(bHat, ylab='Estimated Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_1_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda))
abline(h=fm$ETA[[1]]$lambda,col=4,lwd=2)
abline(v=fm$burnIn/fm$thin,col=4)

unlink("*.dat") #delete files created by BGLR
################################################################################
rm(list=ls())

Y <- pheno_new 
X <- geno_new
# A<-wheat.A #pedigree relatioship matrix from BGLR
y<-Y[,2] # EarHT
n<-nrow(X) # lines
p<-ncol(X) # markers

#Testing set
yNA<-y
set.seed(123)
tst<-sample(1:n,size=100,replace=FALSE)
yNA[tst]<-NA
cbind(y,yNA) #which phenotype removed

#For G
X<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(X)/p #VanRaden (2008)

#Generating heatmaps
heatmap(G) #heatmap for genomic relationship matrix
#heatmap(A) #heatmap for pedigree relationship matrix


#Fitting the G-BLUP model
ETA<-list(list(K=G,model='RKHS'))
fm<-BGLR(y=yNA,ETA=ETA, nIter=5000, burnIn=2000, thin = 5)
unlink("*.dat")

#Correlation graph between pred x observed
plot(fm$yHat,y,xlab="Phenotype",
     ylab="Pred. Gen. Value" ,cex=.8,bty="L")

points(x=y[tst],y=fm$yHat[tst],col=2,cex=.8,pch=19)
legend("topleft", legend=c("training","testing"),
       bty="n",pch=c(1,19), col=c("black","red"))

cbind(fm$y,fm$yHat)# NA have value - GEBV

# correlation comparison in training (TRN) and testing (TST) data sets
cor(fm$yHat[tst],y[tst]) #TEST
cor(fm$yHat[-tst],y[-tst]) #TRAIN

#Extracting parameters from the model
fm$mu #intercept
fm$varE #residual variance
fm$SD.varE #standard eror of residual variance
fm$ETA[[1]]$varU #genetic variance
fm$ETA[[1]]$SD.varU #Standard error of genetic variance
################################################################################
