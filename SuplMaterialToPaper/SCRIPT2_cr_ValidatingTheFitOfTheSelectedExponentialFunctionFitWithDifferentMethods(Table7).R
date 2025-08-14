################################################################################
#                                                                              #
# Supplementary material to Pavel et al. (2025)                                #
#                                                                              #
# SCRIPT 2 - compare different methods to fit the selected                     #
#            exponential function (EXP)                                        #
#        A - Ordinary least squares (OLS)                                      #
#        B - Ordinary least squares with an autoregressive structure (OLS&AR1) #
#        C - Mixed models techniques (MIX)                                     #
#        D - Mixed models techniques and an autoregressive structure(MIX&AR1)  #
#                                                                              #
# This script reads data from R datafile Pb.hbc.RData (Pine trials in Portugal)# 
# and validates the model with a set of 28 plots using the                     #
# leave on out technique                                                       #
#                                                                              #
# the procedure for the variable crown length (cl) is similar                  #
#                                                                              #                                                                             #
################################################################################


##cleaning the last history
rm(list = ls()) 

################################################################################
#Load the required libraries into R workspace
library("tidyverse")
library("readxl")
library("writexl")
library("openxlsx")
library("nlstools")
library("nlme")

################################################################################
# Read the dataset
load("Pb.hbc.Rdata")
summary(Pb.hbc)

################################################################################
# create the output files
# file with the parameters
outfile1 <- createWorkbook()
addWorksheet(outfile1,"OLS")
addWorksheet(outfile1,"OLS&AR1")
addWorksheet(outfile1,"MIX")
addWorksheet(outfile1,"MIX&AR1")

# file with evaluation statistics
outfile2 <- createWorkbook()
addWorksheet(outfile2,"OLS")
addWorksheet(outfile2,"OLS&AR1")
addWorksheet(outfile2,"MIX")
addWorksheet(outfile2,"MIX&AR1")
statistics <- c("bias_cr","bias_cl","bias_hbc",
                "precision_cr","precision_cl","precision_hbc",
                "r2_pred_cr","r2_pred_cl","r2_pred_hbc")

# file with the validation for each tree of the 28 plots
outfile3 <- createWorkbook()
addWorksheet(outfile3,"OLS")
addWorksheet(outfile3,"OLS&AR1")
addWorksheet(outfile3,"MIX")
addWorksheet(outfile3,"MIX&AR1")

################################################################################
# CHANGE HERE THE METHOD
METHOD <- "OLS" # "OLS","OLS&AR1","MIX","MIX&AR1"
TITLE <- paste("Maritime Pine (Portugal) - exponential -",METHOD)

################################################################################
#	Fiting the selected crown ratio model based on the exponential function 
# with the required METHOD
#_______________________________________________________________________________
if (METHOD=="OLS") {
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg))), Pb.hbc, start=list(a0=-0.32,ait=6.2,
              adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,ad_dg=0.184))
  params <- data.frame(coef(mod))
  variables <- row.names(params)
  params <- cbind(variables,params)
  res_stand <- nlsResiduals(mod)$resi2[,2]    # standardized residuals

}
#_______________________________________________________________________________
if (METHOD=="OLS&AR1") {
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg))), Pb.hbc, start=list(a0=-0.32,ait=6.2,
              adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,ad_dg=0.184))

  #### adding autocorrelation structure
  Pb.hbc$res <- residuals(mod)
  Pb.hbc$res_1 <- lag(Pb.hbc$res)
  Pb.hbc$res_1[1] <- 0
  dummy <- ifelse(Cod_Par_Med == lag(Cod_Par_Med),1,0)
  dummy[1] <-0
  
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg+ar1*res_1*dummy))), Pb.hbc, 
              start=list(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,
                       aihdom=0.43,ad_dg=0.184,ar1=0))
  params <- data.frame(coef(mod))
  variables <- row.names(params)
  params <- cbind(variables,params)
  res_stand <- nlsResiduals(mod)$resi2[,2]    # standardized residuals
  
}
#_______________________________________________________________________________
if (METHOD=="MIX") {
  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg))), 
               data=Pb.hbc,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,
               ad_dg=0.184))
  coef <- fixed.effects(mod)
  params <- data.frame(coef) 
  variables <- row.names(params)
  params <- cbind(variables,params)
  res_stand <- residuals(mod,type="normalized")    # standardized residuals  
  
}
#_______________________________________________________________________________
if (METHOD=="MIX&AR1") {
  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg))), 
               data=Pb.hbc,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,
               ad_dg=0.184))

  #### adding autocorrelation structure
  Pb.hbc$res <- residuals(mod)
  Pb.hbc$res_1 <- lag(Pb.hbc$res)
  Pb.hbc$res_1[1] <- 0
  dummy <- ifelse(Cod_Par_Med == lag(Cod_Par_Med),1,0)
  dummy[1] <-0

  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg+ar1*res_1*dummy))),data=Pb.hbc,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg+ar1 ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,ad_dg=0.184,ar1=0))
  coef <- fixed.effects(mod)
  params <- data.frame(coef)
  variables <- row.names(params)
  params <- cbind(variables,params)
  res_stand <- residuals(mod,type="normalized")    # standardized residuals  
  
}

summary(mod)
writeData(outfile1,METHOD,params)

################################################################################																		                                         #
# Plotting observed versus estimated values
yest <- predict(mod)                        # predicted values
df1 <- data.frame(x1=0:1,y1=0:1)
ggplot(Pb.hbc,aes(x=yest,y=cr)) +
  geom_point()+
  labs(title = TITLE) +
  xlim(0.,1.)+
  ylim(0.,1.)+
  labs(x="Estimated crown ratio", y="Observed crown ratio")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,color="1:1 line"))

################################################################################																		                                         #
# Testing regression assumptions                     

# Normality
ggplot(Pb.hbc, aes(sample = res_stand)) +  #Create QQplot with ggplot2 package
  stat_qq() +
  stat_qq_line(col = "red")+
  labs(x="normal quantiles", y="standardized residuals")+
  xlim(-5.,5.)+
  ylim(-5,5.)+
  labs(title = TITLE)+
  theme_minimal()

# Heterocedasticity
df1 <- data.frame(x1=0:1,y1=0)
ggplot(Pb.hbc,aes(x=yest,y=res_stand)) +
  geom_point()+
  labs(title = TITLE) +
  xlim(0.,1.)+
  ylim(-7,7.)+
  labs(x="Estimated crown length", y="Standardized residuals")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,color="red"))

################################################################################
# Validation with re-sampling
listofplots <- list("ASB103","ASB213","ASB108",
                    "COB101","COB104","COB204","COB202",
                    "CPB102","CPB203","CPB201","CPB104",
                    "GOB102","GOB112",
                    "L1B141","L1B132","L1B144","L1B138",
                    "L2B102","L2B108","L2B118","L2B127",
                    "LOB105","LOB202",
                    "PCB102","PCB108",
                    "SSB103","SSB208","SSB311")

niter <- length(listofplots)

for (i in 1:niter){
  Pb.fit <- filter(Pb.hbc,Cod_Par!=listofplots[[i]])
  Pb.val <- filter(Pb.hbc,Cod_Par==listofplots[[i]]) 
  print(length(Pb.fit$Cod_Par))
  print(length(Pb.val$Cod_Par))
#_______________________________________________________________________________
if (METHOD=="OLS") {
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg))), Pb.hbc, start=list(a0=-0.32,ait=6.2,
              adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,ad_dg=0.184))

  a0 <- coefficients(mod)[1]
  ait <- coefficients(mod)[2]
  adh <- coefficients(mod)[3]  
  aih <- coefficients(mod)[4]
  aiG <- coefficients(mod)[5]  
  aihdom <- coefficients(mod)[6]  
  ad_dg <- coefficients(mod)[7] 
 
  Pb.val$cr_est <- with(Pb.val,(1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
                                          aihdom*invhdom+ad_dg*d_dg))))
  Pb.val$cl_est <- with(Pb.val,h*cr_est)
  Pb.val$hbc_est <- with(Pb.val,h-cl_est)
}
#_______________________________________________________________________________  
if (METHOD=="OLS&AR1") {
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg))), 
              Pb.fit,start=list(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,
              aihdom=0.43,ad_dg=0.184))
  
  #### adding autocorrelation structure
  res<-residuals(mod)
  res_1 <- lag(res)
  res_1[1] <- 0
  Pb.aux <- data.frame(res,res_1)
  Pb.fit <- cbind(Pb.fit,Pb.aux)
  Pb.fit$dummy <- with(Pb.fit,ifelse(Cod_Par_Med == lag(Cod_Par_Med),1,0))
  Pb.fit$dummy[1] <-0
  mod <- nls (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
              aihdom*invhdom+ad_dg*d_dg+ar1*res_1*dummy))), Pb.fit, 
              start=list(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,
              aihdom=0.43,ad_dg=0.184,ar1=0))

  a0 <- coefficients(mod)[1]
  ait <- coefficients(mod)[2]
  adh <- coefficients(mod)[3]  
  aih <- coefficients(mod)[4]
  aiG <- coefficients(mod)[5]  
  aihdom <- coefficients(mod)[6]  
  ad_dg <- coefficients(mod)[7] 
  
  Pb.val$cr_est <- with(Pb.val,(1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
                                          aihdom*invhdom+ad_dg*d_dg))))
  Pb.val$cl_est <- with(Pb.val,h*cr_est)
  Pb.val$hbc_est <- with(Pb.val,h-cl_est)
  
} 
#_______________________________________________________________________________  
if(METHOD=="MIX") {
  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg))), 
               data=Pb.fit,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,
               ad_dg=0.184))
  
  a0 <- fixed.effects(mod)[1]
  ait <- fixed.effects(mod)[2]
  adh <- fixed.effects(mod)[3]  
  aih <- fixed.effects(mod)[4]
  aiG <- fixed.effects(mod)[5]  
  aihdom <- fixed.effects(mod)[6]  
  ad_dg <- fixed.effects(mod)[7]
  
  Pb.val$cr_est <- with(Pb.val,(1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
                                          aihdom*invhdom+ad_dg*d_dg))))
  Pb.val$cl_est <- with(Pb.val,h*cr_est)
  Pb.val$hbc_est <- with(Pb.val,h-cl_est)
  
}
#_______________________________________________________________________________
if (METHOD=="MIX&AR1") {
  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg))), 
               data=Pb.fit,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,
               ad_dg=0.184))
  
  #### adding autocorrelation structure
  Pb.fit$res <- residuals(mod)
  Pb.fit$res_1 <- lag(Pb.fit$res)
  Pb.fit$res_1[1] <- 0
  dummy <- ifelse(Pb.fit$Cod_Par_Med == lag(Pb.fit$Cod_Par_Med),1,0)
  dummy[1] <-0
  
  mod <- nlme (cr ~ (1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
               aihdom*invhdom+ad_dg*d_dg+ar1*res_1*dummy))),data=Pb.fit,
               fixed=a0+ait+adh+aih+aiG+aihdom+ad_dg+ar1 ~1,
               random=a0 ~1,
               groups= ~ Cod_Par_Med,
               start=c(a0=-0.32,ait=6.2,adh=0.08,aih=1.05,aiG=4.3,aihdom=0.43,
               ad_dg=0.184,ar1=0))
  
  a0 <- fixed.effects(mod)[1]
  ait <- fixed.effects(mod)[2]
  adh <- fixed.effects(mod)[3]  
  aih <- fixed.effects(mod)[4]
  aiG <- fixed.effects(mod)[5]  
  aihdom <- fixed.effects(mod)[6]  
  ad_dg <- fixed.effects(mod)[7]
  
  Pb.val$cr_est <- with(Pb.val,(1-exp(-(a0+ait*invt+adh*d_h+aih*invh+aiG*invG+
                                          aihdom*invhdom+ad_dg*d_dg))))
  Pb.val$cl_est <- with(Pb.val,h*cr_est)
  Pb.val$hbc_est <- with(Pb.val,h-cl_est)

}
  
  if (i==1){Pb.val_all <- Pb.val}
  if (i>1) {Pb.val_all <- rbind(Pb.val_all,Pb.val)}
  
}

Pb.val_all<-select(Pb.val_all,c(Cod_Par_Med,Id_Trat,Id_Arv,d,h,
                                cr,cr_est,cl,cl_est,hbc,hbc_est,hdom,N,G,S,t))

Pb.val_all$trial <- substr(Pb.val_all$Cod_Par_Med,1,2)
writeData(outfile3,METHOD,Pb.val_all)

res_cr<-with(Pb.val_all,cr-cr_est)
abs_res_cr<-with(Pb.val_all,abs(res_cr))
res_cl<-with(Pb.val_all,cl-cl_est)
abs_res_cl<-with(Pb.val_all,abs(res_cl))
res_hbc<-with(Pb.val_all,hbc-hbc_est)
abs_res_hbc<-with(Pb.val_all,abs(res_hbc))

bias_cr<-mean(res_cr)
precision_cr<-mean(abs_res_cr)
r2_pred_cr<-with(Pb.val_all,1-sum(res_cr**2)/(var(cr)*(length(cr)-1)))

bias_cl<-mean(res_cl)
precision_cl<-mean(abs_res_cl)
r2_pred_cl<-with(Pb.val_all,1-sum(res_cl**2)/(var(cl)*(length(cl)-1)))

bias_hbc<-mean(res_hbc)
precision_hbc<-mean(abs_res_hbc)
r2_pred_hbc<-with(Pb.val_all,1-sum(res_hbc**2)/(var(hbc)*(length(hbc)-1)))
#_______________________________________________________________________________
if (METHOD=="OLS") {
cr_ols <- c(bias_cr,bias_cl,bias_hbc,
            precision_cr,precision_cl,precision_hbc,
            r2_pred_cr,r2_pred_cl,r2_pred_hbc)
cr_ols <- cbind(statistics,cr_ols)
print(cr_ols)
writeData(outfile2,"OLS",cr_ols)

}
#_______________________________________________________________________________
if (METHOD=="OLS&AR1") {
  cr_ols_ar1 <- c(bias_cr,bias_cl,bias_hbc,
              precision_cr,precision_cl,precision_hbc,
              r2_pred_cr,r2_pred_cl,r2_pred_hbc)
  cr_ols_ar1 <- cbind(statistics,cr_ols_ar1)  
  print(cr_ols_ar1)
  writeData(outfile2,"OLS&AR1",cr_ols_ar1)  

}
#_______________________________________________________________________________
if (METHOD=="MIX") {
  cr_mix <- c(bias_cr,bias_cl,bias_hbc,
                  precision_cr,precision_cl,precision_hbc,
                  r2_pred_cr,r2_pred_cl,r2_pred_hbc)
  cr_mix <- cbind(statistics,cr_mix)  
  print(cr_mix)
  writeData(outfile2,"MIX",cr_mix)  

}
#_______________________________________________________________________________
if (METHOD=="MIX&AR1") {
  cr_mix_ar1 <- c(bias_cr,bias_cl,bias_hbc,
              precision_cr,precision_cl,precision_hbc,
              r2_pred_cr,r2_pred_cl,r2_pred_hbc)
  cr_mix_ar1 <- cbind(statistics,cr_mix_ar1)  
  print(cr_mix_ar1)
  writeData(outfile2,"MIX&AR1",cr_mix_ar1)  
  
}

################################################################################
# Plot observed versus estimated values

df1 <- data.frame(x1=0:1,y1=0:1)
ggplot(Pb.val_all,aes(x=cr,y=cr_est)) +
  geom_point()+
  labs(title = TITLE) +
  labs(x="Crown ratio", y="Estimated crown ratio")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,color="red"))

df2 <- data.frame(x1=0:12.5,y1=0:12.5)
ggplot(Pb.val_all,aes(x=cl,y=cl_est)) +
  geom_point()+
  labs(title = TITLE) +
  labs(x="Crown length", y="Estimated crown length")+
  theme_bw()+
  geom_line(data=df2,aes(x=x1,y=y1,color="red"))

df3 <- data.frame(x1=0:25,y1=0:25)
ggplot(Pb.val_all,aes(x=hbc,y=hbc_est)) +
  geom_point()+
  labs(title = TITLE) +
  labs(x="Height to the crown base", y="Estimated height to the crown base")+
  theme_bw()+
  geom_line(data=df3,aes(x=x1,y=y1,color="red"))

################################################################################
# analyse box-plots of residuals over classes of S, N and t
# to analize possible bias

Pb.val_all$Cl_S<-round(Pb.val_all$S)
Pb.val_all$Cl_N<-round(Pb.val_all$N/500)*500
Pb.val_all$Cl_t<-round(Pb.val_all$t/5)*5

TITLE <- paste("crown length residuals (m) - ",METHOD)
boxplot(res_cl~Cl_S,data=Pb.val_all,
        xlab="site index classes",
        ylab=TITLE)
boxplot(res_cl~Cl_N,data=Pb.val_all,
        xlab="stand density classes",
        ylab=TITLE)
boxplot(res_cl~Cl_t,data=Pb.val_all,
        xlab="age classes",
        ylab=TITLE)

################################################################################
# save the output files in EXCEL files
saveWorkbook(outfile1,file="cr_exponential_params.xlsx",overwrite=TRUE)
saveWorkbook(outfile2,file="cr_exponential_evaluation.xlsx",overwrite=TRUE)
saveWorkbook(outfile3,file="cr_exponential_data_val.xlsx",overwrite=TRUE)
