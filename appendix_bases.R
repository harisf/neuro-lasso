getBasis = function(nBases,binSize){
  b = binSize*nBases # IN THESIS TEXT: c <- b
  peaks = c(binSize,binSize*10*nBases)
  
  # nonlinearity for stretching x axis (and its inverse)
  nlin = function(x){log(x+1e-20)}
  invnl = function(x){exp(x)-1e-20}
  
  # Generate basis of raised cosines
  yrange = nlin(peaks+b)
  db = diff(yrange)/(nBases-1)
  centers = seq(yrange[1],yrange[2],db) # IN THESIS TEXT: phi_j <- pi*c[j] / (2*db)
  maxt = invnl(yrange[2]+2*db)-b
  iht = seq(binSize,maxt,binSize)
  nt = length(iht)
  
  raisedCosineBasis = function(x,c,dc){
    (cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2 # IN THESIS TEXT: a <- pi / (2*db)
  }
  
  ihbasis = matrix(NA,nrow = nt,ncol = nBases)
  for(i in seq(1,nt)){
    for(j in seq(1,length(centers))){
      ihbasis[i,j] = raisedCosineBasis(nlin(iht+b)[i],centers[j],db)
    }
  }
  
  
  # orthogonal bases
  library(pracma)
  ihbas = orth(ihbasis)
  
  return(list(bas=ihbasis,bas_orth=ihbas,tau_N=maxt))
}