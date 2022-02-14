source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")
library(forcats)

# prepare data
df=esoph
df$agegp=factor(df$agegp,ordered = FALSE)
df$alcgp=factor(df$alcgp,ordered = FALSE)
df$tobgp=factor(df$tobgp,ordered = FALSE)

#df$alcgp=fct_collapse(df$alcgp,over80=c("80-119","120+"))
df$agegp=fct_collapse(df$agegp,i25_44=c("25-34","35-44"), over65=c("65-74","75+"))
df$tobgp=fct_collapse(df$tobgp,under10=c("0-9g/day"), over10=c("20-29","30+","10-19"))


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
mdr = min_dst_logit(frm,df,weights=df$n,test = tPercentileBootstrap, nSimulation = 1000)
write.result(mdr,"mdr.csv")



# compute distribution of the estimated regression parameters
# at the linear model (so if the linear model were true)
###########################################################

#fit two models to obtain model probabilities

lr = glm(frm,df, family = binomial("logit"), weights =n)

mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

# compute distribution using logit regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=lr$fitted.values,
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_lr_power_lr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=mdr$fitted,
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_lr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=df$p,
                         lr=lr,
                         updateLR =updateLogitModel,nSimulation=1000)
write.results(res,"data_power_lr.csv")

# compute distribution using minimum distance regression 
res=simulatePowerAtModel(df,n=df$n,
                         p=lr$fitted.values,
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_lr_power_mdr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=mdr$fitted,
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"estimation_mdr_power_mdr.csv")

res=simulatePowerAtModel(df,n=df$n,
                         p=df$p,
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"data_set_power_mdr.csv")

# compute test power at the fitted model 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptoticBootstrapVariance, nSimulation = 1000)

res=simulatePowerAtModel(df,
                         n=df$n,
                         p=fitted(mdr),
                         lr=mdr,
                         updateLR =updateMinDistanceModel,nSimulation=1000)
write.results(res,"size_mdr.csv")

# compute test power at the random boundary points 
###########################################################

# obtain minimum distance model for technical and simulate the test power
mdr = min_dst_logit(frm,df,weights=df$n,test = asymptotic)

res= simulatePowerAtBoundary(p=df$p,mdr, nSimulation=1000, eps=0.3)
write.csv(res,"power_mdr.csv")
