funest_fit = function(long_train, surv_train, 
                      noftree = 500, nofcov = 2, split_rule = "maxstat",
                      tv_names, fv_names, nbasis = 3, 
                      nofp = 3, t_star, t_pred, ...){
  # find number of unique individuals in training sample and test sample
  max_time = max(surv_train$time)
  if(t_star >= t_pred){
    stop("t_star must be smaller than t_pred!")
  }
  if(t_pred > max_time){
    stop("t_pred must be smaller than maximum observation time in the dataset!")
  }
  train_patID = surv_train$id
  train_npat = length(train_patID)
  n_of_tv = length(tv_names)
  
  # transfer longitudinal outcomes from long to wide
  multivar.train = array(NA, c(nrow(surv_train), length(unique(long_train$obstime)), n_of_tv))
  sums = 1
  for(i in train_patID){

    visits = long_train$visit[long_train$id == i]
    multivar.train[sums,visits, 1] = long_train$Y1[long_train$id == i]
    multivar.train[sums,visits, 2] = long_train$Y2[long_train$id == i]
    multivar.train[sums,visits, 3]= long_train$Y3[long_train$id == i]
    
    sums = sums + 1
  }
  
  #---- user note: you may want to standardize your time-varying variables 
  # univariate FPCA via PACE
  Xi.train = L = phi.train = meanFun.train = varFun.train =  NULL
  
  for(p in 1:n_of_tv){
    tmp.ufpca = uPACE(multivar.train[,,p], unique(long_train$obstime), nbasis=nbasis)
    Xi.train = cbind(Xi.train, tmp.ufpca$scores) 
    L = c(L, dim(tmp.ufpca$scores)[2])
    phi.train[[p]] = t(tmp.ufpca$functions@X) 
    meanFun.train[[p]] = tmp.ufpca$mu@X 
    varFun.train[[p]] = tmp.ufpca$estVar@X 
  }
  
  # multivariate FPCA
  mFPCA.train = mFPCA(Xi=Xi.train, phi=phi.train, p=n_of_tv, L=L, nsurv = nrow(surv_train) )
  rho.train = mFPCA.train$rho 
  pve = mFPCA.train$pve
  psi = mFPCA.train$psi
  Cms = mFPCA.train$Cms
  
  train_scores = mFPCA.train$rho[,c(1:nbasis)] 
  
  train_covs = cbind(train_scores, surv_train[,unlist(fv_names)])
  d_delta = ncol(train_covs)
  
  score_names = c()
  for(q in 1:(d_delta)){
    tname = paste("score", as.character(q), sep = "")
    score_names = c(score_names, tname)
  }
  colnames(train_covs) = score_names
  
  train_data.sub = cbind(surv_train, train_covs)
  
  fmla = stats::as.formula(paste("Surv(time,event) ~ ",
                                 paste(score_names, collapse= "+")))
  
  
  rg = ranger::ranger(fmla, data = train_data.sub, num.trees = noftree,
                      mtry = nofcov, splitrule = split_rule, ...)
  
  misc = list(
    long_train = long_train,
    surv_train = surv_train,
    fmla = fmla,
    score_names = score_names,
    nofp = nofp,
    train_data.sub = train_data.sub,
    phi.train = phi.train,
    Xi.train = Xi.train,
    nbasis = nbasis,
    L = L,
    multivar.train = multivar.train
  )
  
  return(list(misc = misc,
              rg = rg
  ))
}
