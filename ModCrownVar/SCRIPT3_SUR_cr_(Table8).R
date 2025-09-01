################################################################################
#                                                                              #
# Supplementary material to Pavel et al. (2025)                                #
#                                                                              #
# SCRIPT 3 - fits the system based on cr with SUR                              #
#            with Ordinary least squares (OLS)                                 #
#                                                                              #
# This script reads data from file Pb.hbc.RData (Pine trials in Portugal)      # 
# and validates the model with a set of 28 plots using the                     #
# leave on out technique                                                       #
#                                                                              #
# the procedure for the system based on cl is similar                          #
#                                                                              #
################################################################################

##cleaning the last history
rm(list = ls()) 

################################################################################
#Load the library into R workspace
library("tidyverse")
library("readxl")
library("writexl")
library("openxlsx")
library("nlstools")
library("systemfit")

################################################################################
# Read the dataset
load("Pb.hbc.Rdata")
summary(Pb.hbc)

################################################################################
# create the output files
# file with the parameters
outfile1 <- createWorkbook()
addWorksheet(outfile1,"SUR")

# file with evaluation statistics
outfile2 <- createWorkbook()
statistics <- c("bias_cr","bias_cl","bias_hbc",
                "precision_cr","precision_cl","precision_hbc",
                "r2_pred_cr","r2_pred_cl","r2_pred_hbc")
addWorksheet(outfile2,"SUR")

# file with the validation for each tree of the 28 plots
outfile3 <- createWorkbook()
addWorksheet(outfile3,"SUR")

################################################################################
# Simultaneous fitting of non-linear system of equations with SUR
# Model definition for each variable (cr & cl)
mod_cr <- cr ~(1-exp(-(a0+ait*invt+adh*d_h+aiG*invG+
                       aihdom*invhdom+ad_dg*d_dg)))

mod_cl <- cl ~ h*(1-exp(-(a0+ait*invt+adh*d_h+aiG*invG+
                              aihdom*invhdom+ad_dg*d_dg)))

# Simultaneous fitting
labels        <- list("cl","cr") # label for each model
start.values  <- c(a0=-0.28,ait=7.4,adh=0.11,
                            aiG=3.12,aihdom=1.12,ad_dg=0.08)
system.crown <- list(mod_cl,mod_cr)

model_SUR <- nlsystemfit("SUR", system.crown, start.values, data=Pb.hbc, 
                         eqnlabels=labels)

coeff <- data.frame(model_SUR$b)
variables <- row.names(coeff)
se <- data.frame(model_SUR$se)
t <- data.frame(model_SUR$t)
Pvalue <- data.frame(model_SUR$p)
ci_coeff <- cbind(variables,coeff,se,t,Pvalue)
writeData(outfile1,"SUR",ci_coeff)

print(ci_coeff)

## R2 of the fitted models
sapply(model_SUR$eq, function(x){(x)$r2}) #r2
sapply(model_SUR$eq, function(x){(x)$adjr2}) #adjr2

## Extract fitted and residuals from SUR model
fit_SUR <- data.frame(sapply(model_SUR$eq,function(x){(x)$predicted}))
names(fit_SUR) <- c('fit_cl', 'fit_cr') 
head(fit_SUR)
res_SUR <- data.frame(sapply(model_SUR$eq,function(x){(x)$res}))
names(res_SUR) <- c('res_cl','res_cr')
head(res_SUR)

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
  
  print(paste("niter=",i))
  print(length(Pb.fit$Cod_Par))
  print(length(Pb.val$Cod_Par))
  
  model_SUR <- nlsystemfit("SUR", system.crown, start.values, data=Pb.fit, eqnlabels=labels)
  
  coeff <- data.frame(model_SUR$b)
  pvalues <- data.frame(model_SUR$p)
  coeff <- cbind(coeff,pvalues)

  a0 <- coeff[1,1]
  ait <- coeff[2,1]
  adh <- coeff[3,1]  
  aiG <- coeff[4,1]  
  aihdom <- coeff[5,1]  
  ad_dg <- coeff[6,1] 
  
  Pb.val$cr_est <- with(Pb.val,(1-exp(-(a0+ait*invt+adh*d_h+aiG*invG+
                                          aihdom*invhdom+ad_dg*d_dg))))
  Pb.val$cl_est <- with(Pb.val,h*cr_est)
  Pb.val$hbc_est <- with(Pb.val,h-cl_est)
  if (i==1){Pb.val_all <- Pb.val}
  if (i>1) {Pb.val_all <- rbind(Pb.val_all,Pb.val)}
  
}

Pb.val_all<-select(Pb.val_all,c(Cod_Par_Med,Id_Trat,Id_Arv,d,h,
                                cr,cr_est,cl,cl_est,hbc,hbc_est,hdom,N,G,S,t))
Pb.val_all$trial <- substr(Pb.val_all$Cod_Par_Med,1,2)
writeData(outfile3,"SUR",Pb.val_all)

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


val_SUR <- c(bias_cr,bias_cl,bias_hbc,
            precision_cr,precision_cl,precision_hbc,
            r2_pred_cr,r2_pred_cl,r2_pred_hbc)
val_SUR <- cbind(statistics,val_SUR)
val_SUR

writeData(outfile2,"SUR",val_SUR)

################################################################################
# Plot observed versus estimated values

df1 <- data.frame(x1=0:1,y1=0:1)
ggplot(Pb.val_all,aes(x=cr,y=cr_est)) +
  geom_point()+
  labs(title = "Maritime Pine (Portugal)") +
  labs(x="Crown ratio", y="Estimated crown ratio")+
  theme_bw()+
  geom_line(data=df1,aes(x=x1,y=y1,color="red"))

df2 <- data.frame(x1=0:12.5,y1=0:12.5)
ggplot(Pb.val_all,aes(x=cl,y=cl_est)) +
  geom_point()+
  labs(title = "Maritime Pine (Portugal)") +
  labs(x="Crown length", y="Estimated crown length")+
  theme_bw()+
  geom_line(data=df2,aes(x=x1,y=y1,color="red"))

df3 <- data.frame(x1=0:25,y1=0:25)
ggplot(Pb.val_all,aes(x=hbc,y=hbc_est)) +
  geom_point()+
  labs(title = "Maritime Pine (Portugal)") +
  labs(x="Height to the crown base", y="Estimated height to the crown base")+
  theme_bw()+
  geom_line(data=df2,aes(x=x1,y=y1,color="red"))


################################################################################
# analyse box-plots of residuals over classes of S, N and t
# to analyze possible bias

Pb.val_all$Cl_S<-round(Pb.val_all$S)
Pb.val_all$Cl_N<-round(Pb.val_all$N/500)*500
Pb.val_all$Cl_t<-round(Pb.val_all$t/5)*5

boxplot(res_cl~Cl_S,data=Pb.val_all)
boxplot(res_cl~Cl_N,data=Pb.val_all)
boxplot(res_cl~Cl_t,data=Pb.val_all)

################################################################################
# save the output file in EXCEL files
saveWorkbook(outfile1,file="cr_exponential_SUR_params.xlsx",overwrite=TRUE)
saveWorkbook(outfile2,file="cr_exponential_SUR_evaluation.xlsx",overwrite=TRUE)
saveWorkbook(outfile3,file="cr_exponential_SUR_data_val.xlsx",overwrite=TRUE)
