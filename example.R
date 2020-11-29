#please use this version instead of the R package version (Feb 19, 2020)

source("fit_fun.R")
source("fit_pred.R")
source("utilities.R")

data(long_train)
data(surv_train)

data(long_test)
data(surv_test)

predictSurvProb.ranger = predictor_loader()
w = funest_fit(long_train, surv_train, tv_names = list("Y1", "Y2", "Y3"), noftree = 10, nbasis = 3,
               fv_names = list("W"), t_star = 5.5, t_pred = 11)

# given the trained model, we can assess the prediction accuracy on the testing dataset
# note that trajectories post t_star have already been removed in the validation set. Users need to pre-process data prior to using this function
pred = funest_pred(w, long_test, surv_test, tv_names = list("Y1", "Y2", "Y3"),
                   fv_names = list("W"), t_star = 5.5, t_pred = 11)

# brier score
pred$bs 
# AUC
pred$AUC 
