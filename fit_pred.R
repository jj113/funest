funest_pred = function(funest.fit, long_test, surv_test,
                       tv_names, fv_names,
                       t_star, t_pred){
  # combine data together
  max_time = max(surv_test$time)
  if(t_star >= t_pred){
    stop("t_star must be smaller than t_pred!")
  }
  if(t_pred > max_time){
    stop("t_pred must be smaller than maximum observation time in the dataset!")
  }
  model = funest.fit[[1]]
  rg = funest.fit[[2]]
  
  long_train = model[[1]]
  surv_train = model[[2]]
  fmla = model[[3]]
  score_names = model[[4]]
  nofp = model[[5]]
  train_data.sub = model[[6]]
  phi.train = model[[7]]
  Xi.train = model[[8]]
  nbasis = model[[9]]
  L = model[[10]]
  multivar.train = model[[11]]
  
  test_patID = surv_test$id
  test_npat = length(test_patID)
  n_of_tv = length(tv_names)
  
  # transfer longitudinal outcomes from long to wide
  multivar.test = array(NA, c(nrow(surv_test), length(unique(long_train$obstime)), n_of_tv))
  sums = 1
  for(i in test_patID){
    
    visits = long_test$visit[long_test$id == i]
    multivar.test[sums,visits, 1] = long_test$Y1[long_test$id == i]
    multivar.test[sums,visits, 2] = long_test$Y2[long_test$id == i]
    multivar.test[sums,visits, 3]= long_test$Y3[long_test$id == i]
    
    sums = sums + 1
  }
  
  # univariate FPC  
  Xi.test = NULL
  for(p in 1:n_of_tv){
    tmp.ufpca = uPACE(multivar.train[,,p], unique(long_train$obstime), multivar.test[,,p], nbasis=nbasis)
    Xi.test = cbind(Xi.test, tmp.ufpca$scores) 
  }
  
  # estimate MFPC scores for test subjects
  mFPCA.test = mFPCA(Xi=Xi.train, phi=phi.train, p=n_of_tv, L=L, predXi=Xi.test, nsurv = nrow(surv_test))
  test_scores = mFPCA.test$rho
  
  test_covs = cbind(test_scores, surv_test[, unlist(fv_names)])
  colnames(test_covs) = score_names
  
  test_data.sub = cbind(surv_test, test_covs)
  
  pred.mod = predict(rg, data = test_data.sub, importance = "none")
  
  predictSurvProb.ranger <- function (object, newdata, times, ...) {
    #ptemp <- ranger:::predict.ranger(object, data = newdata, importance = "none")$survival
    ptemp <- predict(object, data = newdata, importance = "none")$survival
    pos <- prodlim::sindex(jump.times = object$unique.death.times,
                           eval.times = times)
    p <- cbind(1, ptemp)[, pos + 1, drop = FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
                 NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ",
                 NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    p
  }
  
  dummy = predictSurvProb.ranger(rg, test_data.sub, test_data.sub$time)
  
  SB = pec(object = rg, formula = fmla,
           traindata = train_data.sub,
           data=test_data.sub,
           start = t_star)
  
  upper = which(SB$time > t_pred)[1]
  lower = which(SB$time <= t_pred)[length(which(SB$time <= t_pred))]
  
  btw = c(lower, upper)
  
  idx = which(abs(c(SB$time[lower], SB$time[upper]) - t_pred) == min(abs(c(SB$time[lower], SB$time[upper]) - t_pred)))
  
  surv_time = pred.mod$unique.death.times
  sb_time = SB$time[btw[idx]]
  surv_ind = max(which(surv_time <= sb_time))
  
  pred_pb = pred.mod$survival[,surv_ind]
  
  surv_pb = data.frame(ID = surv_test$id,
                       pred_pb = pred_pb)
  
  rf.sb = SB$AppErr$ranger[btw[idx]]
  
  timeEvent = test_data.sub[, c("time", "event")]
  #pred.mod = predict(rg, data = test_data.sub, importance = "none")
  
  
  roc = tdROC::tdROC(
    X = 1- cond.prob(pred.mod = pred.mod, newdata = test_data.sub, Tstart = t_star, Tpred = t_pred),
    Y= test_data.sub$time,
    delta = test_data.sub$event,
    tau = t_pred, span = 0.05,
    nboot = 0, alpha = 0.05,
    n.grid = 1000, cut.off = 0.5)
  
  auc.m = roc$AUC$value
  
  return(list(pred_pb = surv_pb, bs = rf.sb, AUC = auc.m))
}

