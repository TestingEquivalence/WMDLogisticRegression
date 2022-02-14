source("dataSets.R")
source("mdLogitRegression.R")
source("size.R")
source("power.R")

# data 

df=as.data.frame(Titanic)
df$Survived=ifelse(df$Survived=="Yes", df$Freq, 0)
df$Deceased=ifelse(df$Survived==0, df$Freq, 0)
df=aggregate(cbind(Survived,Deceased) ~ Class+Sex, df, sum)
df$n=df$Survived+df$Deceased
df$p=df$Survived/df$n

frm="p ~ Class+Sex"


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
mdr = min_dst_logit(frm,df,weights=df$n,test = tPercentileBootstrap, nSimulation = 1000)

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

res= simulatePowerAtBoundary(p=df$p,mdr, nSimulation=1000, eps=0.35)
write.csv(res,"power_mdr.csv")
