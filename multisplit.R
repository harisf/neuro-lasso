library(Matrix)
library(parallel)
library(hdi)
library(glmnet)

##################
# INPUT PARAMETERS
##################
mouseTag = "ANM210862"
sessionTag = "20130626"
trialType = c("goodTrials","leftTrials","rightTrials")[1]
binSize = c(0.01,0.001)[1] # normally, we run multisplit only for binSize = 0.01
##################

modelmatrixDirectory = paste("/global/work/harisf/mette/modelmatrix/",mouseTag,"/",sessionTag,"/",trialType,sep="")
modelmatrixFiles = list.files(modelmatrixDirectory,pattern=paste("_b",binSize*1000,"ms.rds",sep=""))

saveMultisplitIn = paste("/global/work/harisf/mette/multisplit/",mouseTag,"/",sessionTag,"/",trialType,sep="")
if(!dir.exists(saveMultisplitIn))
  dir.create(saveMultisplitIn,recursive = TRUE)

lasso.cv.lambda.min = function (x, y, nfolds = 10, grouped = nrow(x) > 3 * nfolds,...){
  suppressMessages(library(doMC))
  registerDoMC(cores=10)
  fit.cv <- cv.glmnet(x, y, nfolds = nfolds, grouped = grouped, parallel = TRUE,
                      ...)
  sel <- predict(fit.cv, type = "nonzero", s = "lambda.min")
  sel[[1]]
}

glm.pval.x.as.matrix = function (x, y, family = "binomial", verbose = FALSE, ...){
  fit.glm <- glm(y ~ as.matrix(x), family = family, ...)
  fit.summary <- summary(fit.glm)
  if (!fit.glm$converged & verbose) {
    #print(fit.summary)
    cat(" glm.fit: algorithm did not converge.\n")
  }
  pval.sel <- coef(fit.summary)[-1, 4]
  names(pval.sel) <- colnames(x)
  pval.sel
}

for(fileName in modelmatrixFiles){
  modelMatrix = readRDS(file.path(modelmatrixDirectory, fileName))
  
  x = modelMatrix[,-1]
  y = modelMatrix[,1]
  y[which(y > 1)] = 1 # since family = "binomial". Certainly for binSize = 0.01 we get that some y > 1
  
  
  startTime = Sys.time()
  fit <- multi.split(x,y, ci = FALSE, B = 50,
                     classical.fit = glm.pval.x.as.matrix, #args.classical.fit = list(verbose = TRUE),
                     model.selector = lasso.cv.lambda.min, args.model.selector = list(family = "binomial"),
                     parallel = TRUE, ncores = 10,
                     return.selmodels = FALSE, verbose = FALSE)
  endTime = Sys.time() - startTime
  saveRDS(fit,file.path(saveMultisplitIn, fileName))
  cat("Multisplit done for neuron ",substr(fileName,start = 2,stop = regexpr("_",fileName)[1]-1),". (Time used: ",as.numeric(endTime)," ",attr(endTime,"units"),"). \n",sep="")
}










