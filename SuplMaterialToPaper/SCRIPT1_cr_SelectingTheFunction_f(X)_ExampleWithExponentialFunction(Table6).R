################################################################################
#                                                                              #
# Supplementary material to Pavel et al. (2025)                                #
#                                                                              #
# SCRIPT 1 - select the variables to be used in the f(x) function              #
#                                                                              #
# This script reads data from file Pb.hbc.RData (Pine trials in Portugal)      # 
# and test models for crown ratio (cr) using the logistic function             #
#                                                                              #
# the procedure for other functions is similar                                 #
# the procedure for the variable crown length (cl) is also similar             #
#                                                                              #                                                                             #
################################################################################

## Cleaning the workspace
rm(list = ls())

# Load the libraries into R workspace.
library("tidyverse")
library("readxl")
library("openxlsx")
library("nlstools")

# Read the dataset
load("Pb.hbc.Rdata")
summary(Pb.hbc)

attach(Pb.hbc)

################################################################################
#	Fiting a crown ratio model based on the logistic function                    #

# building the base model with a subset of tree variables in f(x)
x <- invt  #test all the variables and select the one that reduces RSE the most
           # in this case it was invt
summary(nls (cr ~ (1-exp(-(a0+ax*x))), Pb.hbc, start=list(a0=0.8497,ax=-0.016)))

x <- d_h   # then select additional tree variables until no additional tree 
           # variable was significant
summary(nls (cr ~ (1-exp(-(a0+ait*invt+ax*x))), Pb.hbc, 
             start=list(a0=-0.016,ait=12.5,ax=0)))

x <- invh
summary(nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+ax*x))), Pb.hbc, 
             start=list(a0=-0.073,ait=15.9,adh=0.25,ax=0)))

# no other tree variable was significant

# test adding stand variables (SV) until no additional SV is significant
# using the best model with tree variables as the base

x <- invG   # add the stand variable that reduces the RSE the most
            # in this case it was invG (1/G)
summary(nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+ax*x))), Pb.hbc, 
             start=list(a0=-0.073,ait=15.9,adh=0.25,aih=1.86,ax=0)))

x <- invhdom
summary(nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+ax*x))),Pb.hbc, 
             start=list(a0=-0.28,ait=6.8,adh=0.22,aih=0.73,aiG=4.3,ax=0)))

x <- SDI
summary(nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
             aihdom*invhdom+ax*x))), Pb.hbc, 
             start=list(a0=-0.28,ait=7.6,adh=0.11,aih=0.83,aiG=3.11,
             aihdom=0.14,ax=0)))

## if adding SDI or Fw they come with a wrong sign, not added

# no other stand variable was significant

# adding a competition index (CI)
# using the best model with tree and stand variables as the base

# then select additional competition indices until no additional CI is significant
x <- d_dg # add the competition index that reduces the RSE the most
          # in this case it was d_dg (d/dg)

summary(nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
             aihdom*invhdom+ax*x))), Pb.hbc, start=list(a0=-0.32,ait=6.2,
             adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,ax=0)))

################################################################################
# calculate fitting and prediction statistic for the selected model
mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h++aih*invh+aiG*invG+aihdom*invhdom+
            ad_dg*d_dg))), Pb.hbc, start=list(a0=-0.32,ait=6.2,adh=0.08,
            aih=1.05,aiG=4.3,aihdom=0.43,ad_dg=0.184))

nparam=length(coef(mod))
summary(mod)
overview(mod)
y<- cr
yest <- predict(mod)                        # predicted values
res <- residuals(mod)                       # ordinary residuals
res_perc <- res/y*100
res_stand <- nlsResiduals(mod)$resi2[,2]    # standardized residuals
res_stand2 <- resid(mod,type="pearson")
mean(res_stand)
mean(res_stand2)
### press residuals
sm <- summary(mod)
W <- mod$m$gradient()                  # W matrix
sm$cov.unscaled
H <- W %*% sm$cov.unscaled%*% t(W)  # sm$cov.unscaled is (W'W)^-1
lev <- diag(H)
rpress <- res/(1-lev)

print(RSS <- sqrt(SSres/(length(y)-nparam)))
print(SSres <- sum(res**2))
print(r2 <- 1-sum(res**2)/(var(y)*(length(y)-1)))
print(adj.r2 <- 1-sum(res**2)/(length(y)-nparam)/var(y))
print(aic <- AIC(mod))
print(PRESS <- sum(rpress**2))
print(r2p <- 1-sum(rpress**2)/(var(y)*(length(y)-1)))
print(bias <- sum(rpress)/length(y))
print(precision <- sum(abs(rpress))/length(y))

##################################################################################																		                                         #
# Plotting Observed versus Predicted values  

df1 <- data.frame(x1=0:1,y1=0:1)
ggplot(Pb.hbc,aes(x=yest,y=cr)) +
  geom_point()+
  labs(title = "Maritime Pine (Portugal)") +
  xlim(0.,1.)+
  ylim(0.,1.)+
  labs(x="Estimated crown ratio", y="Observed crown ratio")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,color="1:1 line"))

##################################################################################																		                                         #
# Testing regression assumptions                     

# Normality
ggplot(Pb.hbc, aes(sample = res_stand)) +  #Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "black")+
  labs(x="normal quantiles", y="standardized residuals")+
  labs(title = "Logistic function")+
  theme_minimal()

# Heterocedasticity
df1 <- data.frame(x1=0:1,y1=0)

ggplot(Pb.hbc,aes(x=yest,y=res_stand)) +
  geom_point()+
  labs(title = "Maritime Pine (Portugal)") +
  xlim(0.,1.)+
  ylim(-7,7.)+
  labs(x="Estimated crown ratio", y="Standardized residuals")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,col="red"))

# Alternative way to make the plots to check normality assumptions

res_all <- nlsResiduals(mod)
plot(res_all, which=0)

##################################################################################																		                                         #
# Writing the results in an EXCEL file with a nice format   

names <- c("RSS","SSres","r2","adj.r2","aic","r2p","bias","precision")
eval <- c(RSS,SSres,r2,adj.r2,aic,r2p,bias,precision)
evaluation <- data.frame(names,eval)
params <- data.frame(coef(mod))
variables <- (row.names(params))
params <- cbind(variables,params)
variables <- c("m/w","a0","at","ait","adh","aih","aiG","aihdom","aG_di","ad_dg")
var_all <- data.frame(variables)
params <- merge(var_all,params,by="variables",all=TRUE)

outfile <- createWorkbook()
addWorksheet(outfile,"params")
addWorksheet(outfile,"evaluation")
writeData(outfile,"params",params)
writeData(outfile,"evaluation",evaluation)
saveWorkbook(outfile,file="cr_exponential.xlsx",overwrite=TRUE)

