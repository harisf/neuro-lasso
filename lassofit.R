library(glmnet)
library(doMC)
library(Matrix)

##################
# INPUT PARAMETERS
##################
mouseTag = "ANM210862"
sessionTag = "20130626"
trialType = c("goodTrials","leftTrials","rightTrials")[1]
binSize = c(0.01,0.001)[2] # normally, we fit lasso models only for binSize = 0.001
##################

modelmatrixDirectory = paste("/global/work/harisf/mette/modelmatrix/",mouseTag,"/",sessionTag,"/",trialType,sep="")
modelmatrixFiles = list.files(modelmatrixDirectory,pattern=paste("_b",binSize*1000,"ms.rds",sep=""))

saveLassofitIn = paste("/global/work/harisf/mette/lassofit/",mouseTag,"/",sessionTag,"/",trialType,sep="")
if(!dir.exists(saveLassofitIn))
  dir.create(saveLassofitIn,recursive = TRUE)

for(fileName in modelmatrixFiles){
  modelMatrix = readRDS(file.path(modelmatrixDirectory, fileName))
  
  x = modelMatrix[,-1]
  y = modelMatrix[,1]
  
  y[which(y > 1)] = 1 # since family = "binomial". Certainly for binSize = 0.01 we get that some y > 1
  
  startTime = Sys.time()
  # fit cv.glmnet model
  registerDoMC(cores = 10)
  model_lasso_cv = cv.glmnet(x,y,
                             family = "binomial",alpha = 1, nfolds = 10,
                             parallel = TRUE)
  endTime = Sys.time() - startTime
  saveRDS(model_lasso_cv,file.path(saveLassofitIn, fileName))
  cat("Lasso model fitted for neuron ",substr(fileName,start = 2,stop = regexpr("_",fileName)[1]-1),". (Time used: ",as.numeric(endTime)," ",attr(endTime,"units"),"). \n",sep="")
}
