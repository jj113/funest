# univariate FPCA via principal analysis by conditional estimation(PACE)
uPACE = function(testData, domain, predData=NULL, nbasis = 10, pve = 0.95, npc = NULL){
  library(MFPCA)
  tmp = funData(domain, testData)
  if(is.null(predData)){
    tmp2 = NULL
  }else{
    tmp2 = funData(domain, predData)
  }
  
  
  res = PACE(tmp, tmp2, pve=pve, npc= npc, nbasis=nbasis)
  
  return(res)
} 

# multivariate FPCA based on results from uPACE
mFPCA = function(Xi, phi, p , L, nsurv,predXi=NULL){
  
  # eigenanalysis on matrix M
  M = t(Xi)%*%Xi/(nsurv-1)
  eigen.M = eigen(M)
  values = eigen.M$values
  pve = cumsum(values)/sum(values)
  Cms = eigen.M$vectors
  index = unlist(lapply(1:length(L), function(x) rep(x, L[x])))
  
  # MFPCA score
  if(is.null(predXi)){
    predXi = Xi
  }
  rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
  for(i in 1:nrow(predXi)){
    for(m in 1:dim(Cms)[2]){
      rho[i,m] = predXi[i,]%*%Cms[,m]
    }
  }
  
  # MFPCA eigenfunction
  psis = NULL
  for(j in 1:p){
    psi = NULL
    for(m in 1:dim(Cms)[2]){
      psi = cbind(psi, phi[[j]]%*%Cms[which(index==j),m])
    }
    psis[[j]] = psi
  }
  
  out = list(eigenvalue = values, Cms = Cms, pve = pve, index=index, rho = rho, psis=psis)
  
  return(out)
}

gen_tv_list = function(n){
  tv_list = list()
  for(i in 1:n){
    tname = paste("Y", as.character(i), sep = '')
    tv_list[[tname]] = list()
  }
  return(tv_list)
}

predictSurvProb.ranger <- function (object, newdata, times, ...) {
  ptemp <- ranger:::predict.ranger(object, data = newdata, importance = "none")$survival
  pos <- prodlim::sindex(jump.times = object$unique.death.times,
                         eval.times = times)
  p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
               NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ",
               NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
  p
}


cond.prob = function(pred.mod, newdata, Tstart, Tpred){
 
  T1 = which(pred.mod$unique.death.times <= Tstart)
  T2 = which(pred.mod$unique.death.times <= Tpred)
  
  T1 = max(T1)
  T2 = max(T2)
  
  risk.Tstart = pred.mod$survival[,T1]
  risk.Tpred = pred.mod$survival[,T2]
  
  return(risk.Tpred/risk.Tstart)
}

gen_arg_list = function(n_of_tv, nofp){
  arg_list = list()
  for(i in 1:n_of_tv){
    arg_list[i] = list(list(type = "uFPCA",nbasis = nofp))
  }
  return(arg_list)
}
