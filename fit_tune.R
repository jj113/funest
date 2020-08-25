library('mlr')
# Define task and learner
task <- makeSurvTask(id = "train_data.sub",
                     data = train_data.sub,
                     target = c("time", "event"))

learner <- makeLearner("surv.ranger")

# Choose resampling strategy and define grid
rdesc <- makeResampleDesc("CV", iters = 5)
ps <- makeParamSet(makeIntegerParam("mtry", 3, 4), # self-specified choices for nofp
                   # self-specified choices for noftree
                   makeDiscreteParam("num.trees", 200), 
                   # self-specified choices for min.node.size
                   makeIntegerParam("min.node.size", 5)) 

# Tune
res = tuneParams(learner, task, rdesc, par.set = ps,
                 control = makeTuneControlGrid())

lrn = setHyperPars(makeLearner("surv.ranger"), par.vals = res$x)

#Train on entire training dataset (using best hyperparameters):
pars = lrn$par.vals