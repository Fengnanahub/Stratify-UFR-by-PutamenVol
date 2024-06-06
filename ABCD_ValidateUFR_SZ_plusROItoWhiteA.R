#############1. Validation in ABCD ##############
## Add Regional GM vol（AAL3 ROI 77-78 and cluster 1， mm3）among 4253 sub's  t1 images(CAT12+8mm) of ABCD_v3 2year， named 
# ~\ABCDV3subareaStriatumVOL.csv 

## and add cat12 TIVetc data, named 
# ~\ABCDV3cat12TIVetc.csv

# to the fist selected 2year ABCD data

# read ABCDV3subareaStriatumVOL.csv；
library(data.table)
library(dplyr)
library(readxl)
# # read csv files
rm(list=ls())
CAT12ROIvol <- read.table('~/ABCDV3subareaStriatumVOL.csv', head = TRUE, sep=",")# 4253*8
CAT12ROIvolhead <- colnames(CAT12ROIvol) # 1*8

setwd('~/fu4_2yAnalysis')
load('~/fu4_2yAnalysis/fu4_2y_datasetV5.rda') # 3086*1274

# 1. add Cat12 striatum ROIvol to 3086*1274 across the 3086 subject
subID1 = fu4_2y_noNApsy_noNaPRSPCAgu1WhiteA$src_subject_id # 3086

colnames(CAT12ROIvol)[c(1)]<- c('src_subject_id')# revise subID as src_subject_id;
subID2 = CAT12ROIvol$src_subject_id # 4253
rowindex1in2 <- match(subID1, subID2) # ROIvol index, 8078, 1396 na;

rm(ROIvol)
ROIvol <- CAT12ROIvol[rowindex1in2,]; # 3086*8

library(dplyr)
inter <-  c('src_subject_id')
PRSPCAguWhiteA_CAT12ROIvol <- fu4_2y_noNApsy_noNaPRSPCAgu1WhiteA %>% full_join(ROIvol, by=inter)# 3086*1281

## 2. add cat12 TIV data to above PRSPCAguWhiteA_CAT12ROIvol;
CAT12TIV <- read.table('~/ABCDV3cat12TIVetc.csv', head = TRUE, sep=",")# 4255*6
CAT12TIVhead <- colnames(CAT12TIV) # 1*6

subID1X = PRSPCAguWhiteA_CAT12ROIvol1$src_subject_id # 3086

colnames(CAT12TIV)[c(1)]<- c('src_subject_id')# revise subID as src_subject_id;
subID2X = CAT12TIV$src_subject_id # 4255
rowindex1in2X <- match(subID1X, subID2X) # TIV index, 3086, 1395 na;

rm(TIVselect)
TIVselect <- CAT12TIV[rowindex1in2X,]; # 3086*6

library(dplyr)
inter <-  c('src_subject_id')
PRSPCAguWhiteA_CAT12ROIvolTIV <- PRSPCAguWhiteA_CAT12ROIvol1 %>% full_join(TIVselect, by=inter)# 3086*1286
##############################################end add 

############# remove the[V3cat12 ROI vol ==NA]&[t1 qc==NA]&[t1qc!=1]--> remain 1680
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PRSpcaGuWhiteA3086_V3cat12ROIvolTIV.rda') # 3086

# remove added cat12 vol NA row 
listNA <- which(is.na(PRSPCAguWhiteA_CAT12ROIvolTIV1$Lcau75volmm3)) # find the CAT12 VOL row contain all NA,1396;
PRSPCAguWhiteA_CAT12ROIvolTIV2 <- PRSPCAguWhiteA_CAT12ROIvolTIV1[-listNA,] # get the data of no all NA row, 1690

# remove the t1qc NA;
listNA_t1qc <- which(is.na(PRSPCAguWhiteA_CAT12ROIvolTIV2[, c('iqc_t1_1_qc_score')])) # find the row NA of t1_qc, 2
hafu_NA <- PRSPCAguWhiteA_CAT12ROIvolTIV2[listNA_t1qc,]#extract row data of NA row, 76;
PRSPCAguWhiteA_CAT12ROIvolTIV3 <- PRSPCAguWhiteA_CAT12ROIvolTIV2[-listNA_t1qc,]# 1688

# t1_qc==1
PRSPCAguWhiteA_CAT12ROIvolTIV4 <- PRSPCAguWhiteA_CAT12ROIvolTIV3[PRSPCAguWhiteA_CAT12ROIvolTIV3$iqc_t1_1_qc_score==1, ]# 1680
#############################################end remove

############extract all interest stress, vol, PSY, PRS, COV ######
colName <- colnames(PRSPCAguWhiteA_CAT12ROIvolTIV4)

match("PRS0001", colName)# 1253
match("PC10", colName)# 1269
ROIvolROIvarAll <- PRSPCAguWhiteA_CAT12ROIvolTIV4[, c('src_subject_id','eventname','interview_age', 'pubertal_sex_p', 'demo_prnt_ed_v2', 'demo_comb_income_v2', 'cat12TIV', 'race_ethnicity', 'site_id_l', 'rel_family_id', 'pps_y_ss_bother_sum', 'pps_y_ss_severity_score', 'pps_ss_mean_severity', 'SUMp_ksads_ptsd','ple_p_ss_total_bad', 'ple_p_ss_affected_bad_sum', 'ple_p_ss_affected_bad_mean', 'ple_y_ss_total_bad','ple_y_ss_affected_bad_sum', 'ple_y_ss_affected_bad_mean', colName[1253:1269], 'Lcau75volmm3', 'Rcau76volmm3', 'Lput77volmm3', 'Rput78colmm3', 'LNAcc157volmm3', 'RNAcc158volmm3', 'GHRcluster1Lcaudatemm3', 'cat12GM')] #all interest psy and stress score, PRS,ROI vol, catTIV.1680*43, +ID, eventName-->1680*45

colnames(ROIvolROIvarAll)[c(3:10)] <- c('age', 'sex', 'edu', 'incom', 'catTIV', 'race', 'siteID', 'FamilyID')
colnames(ROIvolROIvarAll)[c(28:37)] <- c('genePC1', 'genePC2', 'genePC3', 'genePC4', 'genePC5', 'genePC6', 'genePC7', 'genePC8', 'genePC9', 'genePC10')

# trans category var as factor
ROIvolROIvarAll$sex = factor(ROIvolROIvarAll$sex)
ROIvolROIvarAll$race = factor(ROIvolROIvarAll$race)
ROIvolROIvarAll$siteID = factor(ROIvolROIvarAll$siteID)
ROIvolROIvarAll$FamilyID = factor(ROIvolROIvarAll$FamilyID)
######################## end extract (1680*45) ################

## ###### Ana 2 again: using all 1680 subjects;
# age,sex, description 
meanAge <- mean(ROIvolROIvarAll$age) # 142.998
sdAge <- sd(ROIvolROIvarAll$age) # 7.451

meanEdu <- mean(ROIvolROIvarAll$edu) # 18.015
sdEdu <- sd(ROIvolROIvarAll$edu) # 18.632

unique(ROIvolROIvarAll$sex) # 1,2
nrow(ROIvolROIvarAll[ROIvolROIvarAll$sex==1, ]) #933
nrow(ROIvolROIvarAll[ROIvolROIvarAll$sex==2, ]) #747

###### Q4 subjects: Using Cat12 ROI vol median, PRS-SZ median+1sd
PRSin <- match('PRS0001', colnames(ROIvolROIvarAll)) # sz-PRS begins at which col;21
volcutoff <- median(ROIvolROIvarAll$Rput78colmm3) # 4770.1723
PRSmedian <- median(ROIvolROIvarAll$PRS0001) # PRS0001, -0.002217855
PRSmean <- mean(ROIvolROIvarAll$PRS0001) # PRS0001, -0.0022764
PRSSD <- sd(ROIvolROIvarAll$PRS0001) # PRS0001 sd, 0.0003834

## ###### PRS cutoff 1；median+SD
PRScutoff <- PRSmedian+PRSSD# +1SD, -0.001834;

# Binarize column ksads_ptsd
ROIvolROIvarAll$SUMp_ksads_ptsdBin <- ifelse(ROIvolROIvarAll$SUMp_ksads_ptsd <= 0, 0, 1) # 1680*46

dataQudrant1 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll$PRS0001>=PRScutoff), ] # 79
dataQudrant2 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll$PRS0001<PRScutoff), ] # 761
dataQudrant3 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3<volcutoff) & (ROIvolROIvarAll$PRS0001<PRScutoff), ] # 755
dataQudrant4 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3<volcutoff) & (ROIvolROIvarAll$PRS0001>=PRScutoff), ] # 85

####  ROI vol -LESN 相关性；Q4 subjects;
stressin <- match('SUMp_ksads_ptsd', colnames(dataQudrant4)) # LESN 

### (1). ROI（'GHRcluster1Lcaudatemm3' or 'Lcau75volmm3'） vol 与stress相关性# -控制age,sex; Q4被试；
library(Matrix)
result_vol1 <- c()
rm(fitvol1)

for (i in stressin){
  ## ROI vol~stressi+...
  fitvol1 <- lm(GHRcluster1Lcaudatemm3 ~ dataQudrant4[, i]+age + sex + catTIV + siteID, data = dataQudrant4)
  #calculate 95% CI for regression coefficient beta--> add to fitvol1
  fitvol1CI <- cbind(coef(summary(fitvol1)), confint(fitvol1, conf.level = 0.95))
  
  result_vol1<-rbind(result_vol1, c(colnames(dataQudrant4)[i], fitvol1CI[2, ]))
}
result_vol1 <- as.data.frame(result_vol1) 
colnames(result_vol1)[5]<- c('pvalue')


### (2) ROI（'GHRcluster1Lcaudatemm3'） vol 与 binary stress性; Q4被试；
stressin <- match('SUMp_ksads_ptsd', colnames(dataQudrant4copy))
result_vol1 <- c()
rm(fitvol1)
for (i in stressin){
  # Binarize column i
  dataQudrant4copy[, i]<- ifelse(dataQudrant4copy[, i] <= 0, 0, 1)
  
  ## ROI vol~stressi+...
  fitvol1 <- lm(GHRcluster1Lcaudatemm3 ~ dataQudrant4copy[, i]+age + sex +catTIV + siteID, data = dataQudrant4copy)
  #calculate 95% CI for regression coefficient beta--> add to fitvol1
  fitvol1CI <- cbind(coef(summary(fitvol1)), confint(fitvol1, conf.level = 0.95))

  result_vol1<-rbind(result_vol1, c(colnames(dataQudrant4copy)[i], fitvol1CI[2, ]))
}
result_vol1 <- as.data.frame(result_vol1) 
colnames(result_vol1)[5]<- c('pvalue')

### (3) 1680 sub, GHRcluster1Lcaudatemm3~ stress * R_putamenVol * PRS+age+sex+TIV+siteID
rm(datause)
datause <- ROIvolROIvarAll
colnames(ROIvolROIvarAll)

## GHRcluster1 ~ stressi*R_putamenVol*PRS0001
library(Matrix)
rm(fit1)
rm(fit1CI)
result1 <- c()
for (i in stressin){
  fit1 <- lm(GHRcluster1Lcaudatemm3 ~ datause[, i]* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
  fit1coef<-coef(summary(fit1))
  #calculate 95% CI for regression coefficient beta--> add to fit1
  fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
  
  result1<-rbind(result1, c(colnames(datause)[i], fit1CI[nrow(fit1coef), ]))
}
result1 <- as.data.frame(result1) 
colnames(result1)[5]<- c('pvalue')

# (3 again bin) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
# ------ stress为binary SUMp_ksads_ptsd------
datause[, 14]<- ifelse(datause[, 14] <= 0, 0, 1) # binary
fit1 <- lm(GHRcluster1Lcaudatemm3 ~ SUMp_ksads_ptsd* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
fit1coef<-coef(summary(fit1))
#calculate 95% CI for regression coefficient beta--> add to fit1
fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
fit1CI


####### 1680-->remove severity score>6 individuals)
# study Psychosis Risk Screening with the Prodromal Questionnaire – Brief version (PQ-B).Schizophr Res. 2011. 

##### PRS cutoff1（median+1SD）, R_putamen median separated Q4 Ana

# remove severity score>6
length(which(ROIvolROIvarAll$pps_y_ss_severity_score>6)) #200
indexL0 <- which(ROIvolROIvarAll$pps_y_ss_severity_score>6)

ROIvolROIvarAll_Repsy6 <- ROIvolROIvarAll[-indexL0, ]# 1480*46, using below;

volcutoff <- median(ROIvolROIvarAll_Repsy6$Rput78colmm3) # 4774.9701
PRSmedian <- median(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001, -0.002217855
PRSmean <- mean(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001, -0.0022795
PRSSD <- sd(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001 sd, 0.000384

## ---1480 PRS cutoff 1:median+SD;volcutoff:median Q4;
PRScutoff <- PRSmedian+PRSSD# +1SD, -0.001833;

dataQudrant1 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001>=PRScutoff), ] # 65
dataQudrant2 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001<PRScutoff), ] # 675
dataQudrant3 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3<volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001<PRScutoff), ] # 665
dataQudrant4 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3<volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001>=PRScutoff), ] # 75

#### ---->(1) ROI vol -LESN association；1480 Q4 subjects;
stressin <- match('SUMp_ksads_ptsd', colnames(dataQudrant4)) # LESN begins at which col; 14
library(Matrix)
datause <-dataQudrant4
result_vol1 <- c()
rm(fitvol1)
for (i in stressin){
  ## ROI vol~stressi+...
  fitvol1 <- lm(GHRcluster1Lcaudatemm3 ~ datause[, i]+age + sex+ catTIV + siteID, data = datause)
  #calculate 95% CI for regression coefficient beta--> add to fitvol1
  fitvol1CI <- cbind(coef(summary(fitvol1)), confint(fitvol1, conf.level = 0.95))
  
  result_vol1<-rbind(result_vol1, c(colnames(datause)[i], fitvol1CI[2, ]))
}
result_vol1 <- as.data.frame(result_vol1) 
colnames(result_vol1)[5]<- c('pvalue')

# ----> (2)LESN~ binary stress
# (1 again) ROI（'GHRcluster1Lcaudatemm3'） vol 与stress; Q4被试；
stressin <- match('SUMp_ksads_ptsd', colnames(datauseB))
result_vol1 <- c()
rm(fitvol1)
for (i in stressin){
  # Binarize column i
  datauseB[, i]<- ifelse(datauseB[, i] <= 0, 0, 1)
  
  ## ROI vol~stressi+...
  fitvol1 <- lm(GHRcluster1Lcaudatemm3 ~ datauseB[, i]+age + sex +catTIV + siteID, data = datauseB)
  #calculate 95% CI for regression coefficient beta--> add to fitvol1
  fitvol1CI <- cbind(coef(summary(fitvol1)), confint(fitvol1, conf.level = 0.95))
  
  result_vol1<-rbind(result_vol1, c(colnames(datauseB)[i], fitvol1CI[2, ]))
}
result_vol1 <- as.data.frame(result_vol1) 
colnames(result_vol1)[5]<- c('pvalue')

## --->(3)  ROI vol ~ stressi*R_putamenVol*PRS0001 -------
# (1) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
library(Matrix)
rm(fit1)
rm(fit1CI)
result1 <- c()
for (i in stressin){
  fit1 <- lm(GHRcluster1Lcaudatemm3 ~ datause[, i]* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
  fit1coef<-coef(summary(fit1))
  #calculate 95% CI for regression coefficient beta--> add to fit1
  fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
  
  result1<-rbind(result1, c(colnames(datause)[i], fit1CI[nrow(fit1coef), ]))
}
result1 <- as.data.frame(result1) # 7*7
colnames(result1)[5]<- c('pvalue')

# --- (3 again bin) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
# ------ sub 1480 , stress binary SUMp_ksads_ptsd------
datause[, 14]<- ifelse(datause[, 14] <= 0, 0, 1) # binary

fit1 <- lm(GHRcluster1Lcaudatemm3 ~ SUMp_ksads_ptsd* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
fit1coef<-coef(summary(fit1))
#calculate 95% CI for regression coefficient beta--> add to fit1
fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
fit1CI



