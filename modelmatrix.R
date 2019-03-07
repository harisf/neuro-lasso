library(R.matlab)
library(foreach)
library(doParallel)
library(parallel)
library(Matrix)

##################
# INPUT PARAMETERS
##################
mouseTag = "ANM210862"
sessionTag = "20130626"
trialType = c("goodTrials","leftTrials","rightTrials")[1]
binSize = c(0.01,0.001)[2] # should only choose 0.001 or 0.01 (see the function weightedSpikeData)
deg = 5 # degree of lickOnset polynomial
nBases_history = 10
nBases_connectivity = 4
##################

# read data
data = readMat(paste("/global/work/harisf/mette/data/data_structure_",mouseTag,"/data_structure_",mouseTag,"_",sessionTag,".mat",sep=""))
totalNumberOfNeurons = length(data$obj[[12]][[1]])

# get trials that should be analyzed
getTrials = function(trialType){
  if(trialType == "goodTrials"){
    trials.good = which(data$obj[[9]][[3]][[4]][[1]] == 1) # trials where mice are performing (should be tested)
    trials.photostimConfig = which(is.nan(data$obj[[9]][[3]][[5]][[1]])) # trials where photostimulation configuration is tested (should NOT be tested)
    trials.good = trials.good[!is.element(trials.good,trials.photostimConfig)] # trials where mice are performing, AND we've taken out trials where photostimulation configuration is tested
    trials = trials.good
  }
  if(trialType == "rightTrials"){
    trials_correctR = which(data$obj[[8]][1,] == 1)
    trials_correctR = trials_correctR[is.element(trials_correctR,trials.good)]
    trials = trials_correctR
  }
  if(trialType == "leftTrials"){
    trials_correctL = which(data$obj[[8]][2,] == 1)
    trials_correctL = trials_correctL[is.element(trials_correctL,trials.good)]
    trials = trials_correctL
  }
  return(trials)
}
trials = getTrials(trialType)

discretizeSpikeData = function(neuron,TRIALS,binSize){
  eventTrials = data$obj[[12]][[3]][[neuron]][[1]][[3]]
  eventTimes = data$obj[[12]][[3]][[neuron]][[1]][[2]]
  lickOnset = mean(c(data$obj[[9]][[3]][[3]][[1]][TRIALS]),na.rm=TRUE)
  
  timeInterval = seq(0,5,binSize) # each trial lasts 5 seconds (this is true in all 3 sessions of  mouse ANM210861, but is it true for all other mice?)
  
  mat_j = matrix(NA,ncol=3,nrow=length(timeInterval)-1)
  
  registerDoParallel(cores = detectCores()-1)
  mat = foreach(trial_j = TRIALS,.combine = rbind) %dopar% {
    trialStartTime_j = data$obj[[7]][1,trial_j]
    eventTimes_j = eventTimes[which(eventTrials == trial_j)] - trialStartTime_j
    
    mat_j[,1] = rep(trial_j,length(timeInterval)-1)
    mat_j[,2] = timeInterval[2:length(timeInterval)]-lickOnset
    mat_j[,3] = as.vector(table(cut(eventTimes_j,breaks=timeInterval)))
    
    mat_j
  }
  stopImplicitCluster()
  
  colnames(mat) = c("trialId","lickOnset",paste("spikeCountj",neuron,sep=""))
  mat = as.data.frame(mat)
  return(mat)
}

discretizeAndAlignSpikeData = function(mainNeuron,TRIALS,binSize){
  spikeData_mainNeuron = discretizeSpikeData(mainNeuron,TRIALS,binSize)
  
  otherNeurons = setdiff(seq(1,totalNumberOfNeurons),mainNeuron)
  
  timeInterval = seq(0,5,binSize)
  
  registerDoParallel(cores = detectCores()-1)
  spikeData_otherNeurons = foreach(i = seq(1,totalNumberOfNeurons-1),.combine = cbind) %dopar% {
    neuron = otherNeurons[i]
    eventTimes = data$obj[[12]][[3]][[neuron]][[1]][[2]]
    eventTrials = data$obj[[12]][[3]][[neuron]][[1]][[3]]

    registerDoParallel(cores = detectCores()-1)
    spikesInTrial = foreach(trial_j = unique(spikeData_mainNeuron$trialId), .combine = list) %dopar% {
      trialStartTime_j = data$obj[[7]][1,trial_j]
      eventTimes_j = eventTimes[which(eventTrials == trial_j)] - trialStartTime_j
      spikevector = as.vector(table(cut(eventTimes_j,breaks=timeInterval)))
    }
    stopImplicitCluster()
    
    spikeData_neuron = unlist(spikesInTrial)
  }
  stopImplicitCluster()
  
  colnamesTxt = NULL
  for(i in seq(1,totalNumberOfNeurons-1))
    colnamesTxt = c(colnamesTxt,paste("spikeCountj",otherNeurons[i],sep = ""))
  colnames(spikeData_otherNeurons) = colnamesTxt
  
  return(as.data.frame(cbind(spikeData_mainNeuron,spikeData_otherNeurons)))
  
}

getBasis = function(nBases,binSize){
  b = binSize*nBases
  peaks = c(binSize,binSize*10*nBases)
  
  # nonlinearity for stretching x axis (and its inverse)
  nlin = function(x){log(x+1e-20)}
  invnl = function(x){exp(x)-1e-20}
  
  # Generate basis of raised cosines
  yrange = nlin(peaks+b)
  db = diff(yrange)/(nBases-1)
  centers = seq(yrange[1],yrange[2],db)
  maxt = invnl(yrange[2]+2*db)-b
  iht = seq(binSize,maxt,binSize)
  nt = length(iht)
  
  raisedCosineBasis = function(x,c,dc){
    (cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2
  }
  
  ihbasis = matrix(NA,nrow = nt,ncol = nBases)
  for(i in seq(1,nt)){
    for(j in seq(1,length(centers))){
      ihbasis[i,j] = raisedCosineBasis(nlin(iht+b)[i],centers[j],db)
    }
  }
  
  #matplot(ihbasis,type="b",pch=seq(1,5))
  # for plotting model coefficients
  lags = invnl(centers)-b
  
  library(pracma)
  ihbas = orth(ihbasis) # orthogonal bases
  
  return(list(bas=ihbasis,bas_orth=ihbas,lags=lags,tau_N=maxt))
}

getMedianBasis = function(nBases,binSize_original,binSize_median){
  bas_original = getBasis(nBases,binSize = binSize_original)
  
  pointsToCombine = binSize_median / binSize_original
  totalCombinedPoints = floor(dim(bas_original$bas)[1] / pointsToCombine)
  
  indexes = list()
  for(i in seq(0,totalCombinedPoints - 1))
    indexes = c(indexes,list(seq(1,pointsToCombine) + pointsToCombine * i))
  
  binSize_median = list(bas = matrix(NA,ncol = dim(bas_original$bas)[2],nrow = length(indexes)), 
                        bas_orth = matrix(NA,ncol = dim(bas_original$bas)[2],nrow = length(indexes)),
                        tau_N = bas_original$tau_N)
  
  for(i in seq(1,length(indexes))){
    binSize_median$bas[i,] = bas_original$bas[floor(median(indexes[[i]])),]
    binSize_median$bas_orth[i,] = bas_original$bas_orth[floor(median(indexes[[i]])),]
  }
  
  binSize_median
}

weightedSpikeData = function(spikeData,nBases_history,nBases_connectivity,binSize){
  if(binSize == 0.001){
    bas_hist = getBasis(nBases_history,binSize)$bas_orth
    bas_connect = getBasis(nBases_connectivity,binSize)$bas_orth
  }
  if(binSize == 0.01){
    bas_hist = getMedianBasis(nBases_history,binSize_original = 0.001,binSize_median = binSize)$bas_orth
    bas_connect = getMedianBasis(nBases_connectivity,binSize_original = 0.001,binSize_median = binSize)$bas_orth
  }
  
  #history
  registerDoParallel(cores = nBases_history)
  basisWeights = foreach(k = seq(1,nBases_history),.combine = cbind) %dopar%{
    if(binSize == 0.01)
      bsWght = convolve(c(0,spikeData[,3]),rev(bas_hist[,k]),type="open")[2:dim(spikeData)[1]]
    if(binSize == 0.001)
      bsWght = convolve(c(0,spikeData[,3]),rev(bas_hist[1:160,k]),type="open")[2:dim(spikeData)[1]] # ideally, should convolve using bas_hist[,k], but to shorten computation time we use the first 160 rows
    bsWght
  }
  stopImplicitCluster()
  
  txt=NULL
  for(k in seq(1,nBases_history)){
    txt = c(txt,paste(sub("spikeCount","",colnames(spikeData)[3]),".k",k,sep=""))
  }
  colnames(basisWeights) = txt
  
  spikeData_basis = cbind(spikeData[-1,seq(1,3)],basisWeights) #when using convolution
  colnames(spikeData_basis)[1:3] = colnames(spikeData)[1:3]
  spikeData_basis = as.data.frame(spikeData_basis)
  
  #connectivity
  registerDoParallel(cores = detectCores() - 1)
  spikeData_basis_connectivity = foreach(j = seq(4,dim(spikeData)[2]),.combine = cbind) %dopar% {
    registerDoParallel(cores = nBases_connectivity)
    basisWeights = foreach(k = seq(1,nBases_connectivity),.combine = cbind) %dopar% {
      if(binSize == 0.01)
        bsWght = convolve(c(0,spikeData[,j]),rev(bas_connect[,k]),type="open")[2:dim(spikeData)[1]]
      if(binSize == 0.001)
        bsWght = convolve(c(0,spikeData[,j]),rev(bas_connect[1:160,k]),type="open")[2:dim(spikeData)[1]] # # ideally, should convolve using bas_connect[,k], but to shorten computation time we use the first 160 rows
      bsWght
    }
    stopImplicitCluster()
    
    txt=NULL
    for(k in seq(1,nBases_connectivity)){
      txt = c(txt,paste(sub("spikeCount","",colnames(spikeData)[j]),".k",k,sep=""))
    }
    colnames(basisWeights) = txt
    
    basisWeights
  }
  stopImplicitCluster()
  
  spikeData_basis = cbind(spikeData_basis,spikeData_basis_connectivity)
  
  return(spikeData_basis)
  
}

saveModelmatrixIn = paste("/global/work/harisf/mette/modelmatrix/",mouseTag,"/",sessionTag,"/",trialType,sep="")
if(!dir.exists(saveModelmatrixIn))
  dir.create(saveModelmatrixIn,recursive = TRUE)

for(responseNeuron in seq(1,totalNumberOfNeurons)){
  start = Sys.time()
  spikeData <- discretizeAndAlignSpikeData(responseNeuron,trials,binSize)
  
  spikeData_basis <- weightedSpikeData(spikeData,nBases_history,nBases_connectivity,binSize)
  
  predMat = cbind(spikeData_basis[,-c(1,2,3)],poly(spikeData_basis[,2],degree = deg))
  txtPolyname = NULL
  for(i in seq(1,deg))
    txtPolyname = c(txtPolyname,paste("poly(lickOnset)",i,sep = ""))
  colnames(predMat)[seq(dim(predMat)[2]-deg+1,dim(predMat)[2])] = txtPolyname
  
  # save data
  evaluatedData = cbind(spikeData_basis[,3],predMat)
  colnames(evaluatedData)[1] = paste("spike.j",responseNeuron,sep="")
  tol = 1e-10
  for(i in seq(2,dim(evaluatedData)[2])){
    evaluatedData[which(abs(evaluatedData[,i]) < tol),i] = 0
  }
  evaluatedData = Matrix(as.matrix(evaluatedData),sparse = TRUE)
  
  saveRDS(evaluatedData,paste(saveModelmatrixIn,"/n",responseNeuron,"_b",binSize*1000,"ms.rds",sep=""))
  end = Sys.time() - start
  cat("Model matrix saved for neuron ",responseNeuron,". (Time used: ",as.numeric(end)," ",attr(end,"units"),"). \n",sep="")
}








