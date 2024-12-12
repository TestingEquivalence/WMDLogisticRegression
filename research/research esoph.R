source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")

# prepare data
df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)

df=aggregate(cbind(ncases,ncontrols) ~ alcgp+tobgp+agegp, df, sum)
df$n=df$ncases+df$ncontrols
df$p=df$ncases/(df$ncases+df$ncontrols)

frm="p ~ agegp + alcgp+ tobgp"


# fitting the model and perform a single equivalence tests
###########################################################

# using logit regression
lr <- glm(frm, df, family = binomial("logit"), weights = n)
v=(df$p-lr$fitted.values)*weights(lr)/sum(weights(lr))
lr$min.distance=sqrt(sum(v*v))
write.result(lr,"lr.csv")

# using minimum distance regression
set.seed(01012021)
mdr = min_dst_logit(frm,df,weights=df$n,test = none, nSimulation = 1000, alpha=0.05,
                    fixIntercept = TRUE)
write.result(mdr,"mdr.csv")

# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities

lr = glm(frm,df, family = binomial("logit"), weights =n)

mdr = min_dst_logit(frm,df,weights=df$n,test = none, alpha=0.05, nSimulation = 0, fixIntercept = TRUE)

# compute distribution using logit regression 
res=simulatePowerLR(p=lr$fitted.values,nSimulation=1000,lr=lr, mdr)
write.results(res,"estimation_lr_power_lr.csv")

res=simulatePowerLR(p=mdr$fitted,nSimulation=1000,lr=lr, mdr)
write.results(res,"estimation_mdr_power_lr.csv")

res=simulatePowerLR(p=df$p,nSimulation=1000,lr=lr, mdr)
write.results(res,"data_power_lr.csv")

# compute distribution using minimum distance regression 
res=simulatePowerMDR(p=lr$fitted.values,nSimulation=1000, mdr)
write.results(res,"estimation_lr_power_mdr.csv")

res=simulatePowerMDR(p=mdr$fitted,nSimulation=1000, mdr)
write.results(res,"estimation_mdr_power_mdr.csv")

res=simulatePowerMDR(p=df$p,nSimulation=1000, mdr)
write.results(res,"data_set_power_mdr.csv")

# compute test power at the fitted model 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptoticBootstrapVariance, 
                    nSimulation = 200)

res=simulatePowerMDR(p=fitted(mdr),nSimulation = 1000,mdr = mdr)
write.results(res,"size_mdr.csv")

# compute test power at the random boundary points 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

res= simulatePowerAtBoundary(p=mdr$fitted,mdr, nSimulation=1000, eps=0.008)
write.csv(res,"power_mdr.csv")


