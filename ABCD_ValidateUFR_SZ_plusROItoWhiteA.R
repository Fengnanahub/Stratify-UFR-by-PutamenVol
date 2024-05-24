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
PRSPCAguWhiteA_CAT12ROIvol <- fu4_2y_noNApsy_noNaPRSPCAgu1WhiteA %>% full_join(ROIvol, by=inter)# 4482*1281

# remove all NA row (including src_subject_id)
listNA <- which(rowSums(is.na(PRSPCAguWhiteA_CAT12ROIvol[, c(1:1281)]))>=1280) # find the row contain all NA,1396;

PRSPCAguWhiteA_CAT12ROIvol1 <- PRSPCAguWhiteA_CAT12ROIvol[-listNA,] # get the data of no all NA row, 3086*1281

save(PRSPCAguWhiteA_CAT12ROIvol1, file = '~/PRSpcaGuWhiteA3086_V3cat12ROIvol.rda')## add ABCDv3 CAT12 ROI vol; 3086*1281;


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
PRSPCAguWhiteA_CAT12ROIvolTIV <- PRSPCAguWhiteA_CAT12ROIvol1 %>% full_join(TIVselect, by=inter)# 4481*1286

# remove all NA row (including src_subject_id)
listNAx <- which(rowSums(is.na(PRSPCAguWhiteA_CAT12ROIvolTIV[, c(1:1286)]))>=1285) # find the row contain all NA,1395;

PRSPCAguWhiteA_CAT12ROIvolTIV1 <- PRSPCAguWhiteA_CAT12ROIvolTIV[-listNAx,] # get the data of no all NA row, 3086*1286

save(PRSPCAguWhiteA_CAT12ROIvolTIV1, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PRSpcaGuWhiteA3086_V3cat12ROIvolTIV.rda')## add ABCDv3 CAT12 TIV; 3086*1286;
##############################################end add 



############# remove the[V3cat12 ROI vol ==NA]&[t1 qc==NA]&[t1qc!=1]--> remain 1680*1281
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PRSpcaGuWhiteA3086_V3cat12ROIvolTIV.rda') # 3086*1286

# remove added cat12 vol NA row 
listNA <- which(is.na(PRSPCAguWhiteA_CAT12ROIvolTIV1$Lcau75volmm3)) # find the CAT12 VOL row contain all NA,1396;
PRSPCAguWhiteA_CAT12ROIvolTIV2 <- PRSPCAguWhiteA_CAT12ROIvolTIV1[-listNA,] # get the data of no all NA row, 1690*1286

# remove the t1qc NA;
listNA_t1qc <- which(is.na(PRSPCAguWhiteA_CAT12ROIvolTIV2[, c('iqc_t1_1_qc_score')])) # find the row NA of t1_qc, 2
hafu_NA <- PRSPCAguWhiteA_CAT12ROIvolTIV2[listNA_t1qc,]#extract row data of NA row, 76;
PRSPCAguWhiteA_CAT12ROIvolTIV3 <- PRSPCAguWhiteA_CAT12ROIvolTIV2[-listNA_t1qc,]# 1688*1286

# t1_qc==1
PRSPCAguWhiteA_CAT12ROIvolTIV4 <- PRSPCAguWhiteA_CAT12ROIvolTIV3[PRSPCAguWhiteA_CAT12ROIvolTIV3$iqc_t1_1_qc_score==1, ]# 1680*1286

save(PRSPCAguWhiteA_CAT12ROIvolTIV4, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680.rda')
#############################################end remove



############extract all interest stress, vol, PSY, PRS, COV ######
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680.rda')# 1680*1286

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
save(ROIvolROIvarAll, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda')# 1680*45

write.csv(ROIvolROIvarAll, '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.csv', row.names=FALSE)# 1680*45
######################## end extract (1680*45) ################


# find the NA num in each col of 1860*43
rm(dataneedPro)

dataneedPro <- ROIvolROIvarAll
fu4_2y_ColNA1680 <- matrix(nrow = 2, ncol = ncol(dataneedPro))
for (j in 1:ncol(dataneedPro)){
  fu4_2y_ColNA1680[1, j] <- sum(is.na(dataneedPro[, j]))
  fu4_2y_ColNA1680[2, j] <- nrow(dataneedPro)- fu4_2y_ColNA1680[1, j]
}
colnames(fu4_2y_ColNA1680) <- c(colnames(dataneedPro))
rownames(fu4_2y_ColNA1680) <- c('NAnumeachCol', 'noNAnumofeachCol')

save(fu4_2y_ColNA1680, '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVarNAsum.rda', row.names=FALSE)# 1680*45


## ###### Ana 2 again: using all 1680 subjects;
load('D:/DATA_Fudan/ABCDdata/New_ABCDanaAgain20231115/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda')# 1680*45
# age,sex, description 
meanAge <- mean(ROIvolROIvarAll$age) # 142.998
sdAge <- sd(ROIvolROIvarAll$age) # 7.451

meanEdu <- mean(ROIvolROIvarAll$edu) # 18.015
sdEdu <- sd(ROIvolROIvarAll$edu) # 18.632

unique(ROIvolROIvarAll$sex) # 1,2
nrow(ROIvolROIvarAll[ROIvolROIvarAll$sex==1, ]) #933
nrow(ROIvolROIvarAll[ROIvolROIvarAll$sex==2, ]) #747

# Plot histograms for each sym, stress, vol
rm(data_hist)
rm(colN)

data_hist <- ROIvolROIvarAll[, c( 'Lput77volmm3','Rput78colmm3', 'GHRcluster1Lcaudatemm3')]
colN <- colnames(data_hist)

# layout for multiple hist plots
par(mfrow = c(2, 4))
for (i in 1:length(colN)) {
  hist(data_hist[, i], main = colN[i], xlab=paste("no_NA_num",fu4_2y_ColNA1680[2, c(colN[i])]))
}


#################### Ana6 Again##########################################
###### 6.5. Q4 subjects: Using Cat12 ROI vol 20231217, 以其他阈值划分PRS
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda') # 1680*45

PRSin <- match('PRS0001', colnames(ROIvolROIvarAll)) # sz-PRS begins at which col;21
PRSfinal <- match('PRS05', colnames(ROIvolROIvarAll))#27

volcutoff <- median(ROIvolROIvarAll$Rput78colmm3) # 4770.1723
PRSmedian <- median(ROIvolROIvarAll$PRS0001) # PRS0001, -0.002217855
PRSmean <- mean(ROIvolROIvarAll$PRS0001) # PRS0001, -0.0022764
PRSSD <- sd(ROIvolROIvarAll$PRS0001) # PRS0001 sd, 0.0003834

## ###### PRS cutoff 1；median+SD
PRScutoff <- PRSmedian+PRSSD# +1SD, -0.001834;

# Binarize column ksads_ptsd
ROIvolROIvarAll$SUMp_ksads_ptsdBin <- ifelse(ROIvolROIvarAll$SUMp_ksads_ptsd <= 0, 0, 1) # 1680*46

dataQudrant1 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll$PRS0001>=PRScutoff), ] # 79*46
dataQudrant2 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll$PRS0001<PRScutoff), ] # 761*46
dataQudrant3 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3<volcutoff) & (ROIvolROIvarAll$PRS0001<PRScutoff), ] # 755*46
dataQudrant4 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3<volcutoff) & (ROIvolROIvarAll$PRS0001>=PRScutoff), ] # 85*46
save(dataQudrant1, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/1680PRS0001cat12RputaVol_dataQ4x/t1qc1680_median1SD_Q1_s79ROIdata.rda')#
save(dataQudrant2, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/1680PRS0001cat12RputaVol_dataQ4x/t1qc1680_median1SD_Q2_s761_ROIdata.rda')#
save(dataQudrant3, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/1680PRS0001cat12RputaVol_dataQ4x/t1qc1680_median1SD_Q3_s755_ROIdata.rda')#
save(dataQudrant4, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/1680PRS0001cat12RputaVol_dataQ4x/t1qc1680_median1SD_Q4_s85ROIdata.rda')#

# plot GHRclu1vol by SUMp_ksds_ptsd in each Qudrant
## same PIC for 4 qudrant to manuscript;--- final--20240508
library(tidyverse)
# Q1
q1.manu <-dataQudrant1 |>
  group_by(SUMp_ksads_ptsdBin) |>
  summarise(ym = mean(GHRcluster1Lcaudatemm3/1000), 
            ystd = sd(GHRcluster1Lcaudatemm3/1000),
            yrad = 1.96 * ystd / sqrt(n()))
pointQ1 <- ggplot(data=q1.manu) + theme_bw() + theme(panel.grid=element_blank()) + geom_jitter(data=dataQudrant1, mapping=aes(
  x = SUMp_ksads_ptsdBin, y = GHRcluster1Lcaudatemm3/1000), colour = '#ED6F6E',alpha = 0.8, size = 1.8, width=0.25)+geom_errorbar(data = q1.manu, mapping = aes(x = SUMp_ksads_ptsdBin, ymin = ym - yrad, ymax = ym + yrad), colour = '#4C9AC9', width=0.3, size=0.2) +
  geom_point(mapping = aes(x = SUMp_ksads_ptsdBin, y = ym ),alpha = 0.8, size=1.8, color="#3371B3") + xlim(-0.5,1.5)+ylab("Volume of cluster caudate") + theme(axis.title.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+theme(axis.text.x = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+theme(axis.text.y = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
# Q2
q2.manu <-dataQudrant2 |>
  group_by(SUMp_ksads_ptsdBin) |>
  summarise(ym = mean(GHRcluster1Lcaudatemm3/1000), 
            ystd = sd(GHRcluster1Lcaudatemm3/1000),
            yrad = 1.96 * ystd / sqrt(n()))
pointQ2 <- ggplot(data=q2.manu) + theme_bw() + theme(panel.grid=element_blank()) + geom_jitter(data=dataQudrant2, mapping=aes(x = SUMp_ksads_ptsdBin, y = GHRcluster1Lcaudatemm3/1000), colour = '#ED6F6E',alpha = 0.75, size = 1.8, width=0.25)+geom_errorbar(data = q2.manu, mapping = aes(x = SUMp_ksads_ptsdBin, ymin = ym - yrad, ymax = ym + yrad), colour = '#4C9AC9', width=0.3,size=0.2) +
  geom_point(mapping = aes(x = SUMp_ksads_ptsdBin, y = ym ),alpha = 0.75, size=1.8, color="#3371B3") + xlim(-0.5,1.5)+ylab("Volume of cluster caudate") + theme(axis.title.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+theme(axis.text.x = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+theme(axis.text.y = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
# Q3
q3.manu <-dataQudrant3 |>
  group_by(SUMp_ksads_ptsdBin) |>
  summarise(ym = mean(GHRcluster1Lcaudatemm3/1000), 
            ystd = sd(GHRcluster1Lcaudatemm3/1000),
            yrad = 1.96 * ystd / sqrt(n()))
pointQ3 <- ggplot(data=q3.manu) + theme_bw() + theme(panel.grid=element_blank()) + geom_jitter(data=dataQudrant3, mapping=aes(x = SUMp_ksads_ptsdBin, y = GHRcluster1Lcaudatemm3/1000), colour = '#ED6F6E', alpha = 0.75, size = 1.8, width=0.25)+geom_errorbar(data = q3.manu, mapping = aes(x = SUMp_ksads_ptsdBin, ymin = ym - yrad, ymax = ym + yrad), colour = '#4C9AC9', width=0.3, size=0.2) +
  geom_point(mapping = aes(x = SUMp_ksads_ptsdBin, y = ym ), alpha = 0.75,size=1.8, color="#3371B3") + xlim(-0.5,1.5)+ylab("Volume of cluster caudate") + theme(axis.title.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+theme(axis.text.x = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+theme(axis.text.y = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))
# Q4
q4.manu <-dataQudrant4 |>
  group_by(SUMp_ksads_ptsdBin) |>
  summarise(ym = mean(GHRcluster1Lcaudatemm3/1000), 
            ystd = sd(GHRcluster1Lcaudatemm3/1000),
            yrad = 1.96 * ystd / sqrt(n()))
pointQ4 <- ggplot(data=q4.manu) + theme_bw() + theme(panel.grid=element_blank())  +
  geom_jitter(data=dataQudrant4, mapping=aes(
    x = SUMp_ksads_ptsdBin, y = GHRcluster1Lcaudatemm3/1000), colour = '#ED6F6E', size = 1.8, alpha = 0.8,width=0.25)+geom_errorbar(data = q4.manu, mapping = aes(x = SUMp_ksads_ptsdBin, ymin = ym - yrad, ymax = ym + yrad), colour = '#4C9AC9', width=0.3, size=0.2) +
  geom_point(mapping = aes(
    x = SUMp_ksads_ptsdBin, y = ym),alpha = 0.8, size=1.8, color="#3371B3")+ xlim(-0.5,1.5)+ylab("Volume of cluster caudate") + theme(axis.title.x = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 10, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 90))+theme(axis.text.x = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))+theme(axis.text.y = element_text(size = 9, color = "black", face = "bold", vjust = 0.5, hjust = 0.5, angle = 0))

library(patchwork)
(pointQ2+pointQ1)/(pointQ3+pointQ4)
##－－－－－－－－－－－－－－－－－－－－－－－－－－－－


##### PRS cutoff1, R_putamen median separated Q4 Ana
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda') # 1680*45

dataQudrant4 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colmm3<volcutoff) & (ROIvolROIvarAll$PRS0001>=PRScutoff), ] # 85*45

#### 6.5.2 ROI vol -症状或LESN 相关性；Q4 subjects;
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
result_vol1 <- as.data.frame(result_vol1) # 7*7
colnames(result_vol1)[5]<- c('pvalue')
result_vol1$pvalue.Adj<- p.adjust(result_vol1$pvalue, method = 'fdr', n=length(result_vol1$pvalue)) ## 7*8


# (again, ptsd binary). ROI（'GHRcluster1Lcaudatemm3'） vol 与stress性; Q4被试；
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
result_vol1 <- as.data.frame(result_vol1) # 5*7
colnames(result_vol1)[5]<- c('pvalue')
result_vol1$pvalue.Adj<- p.adjust(result_vol1$pvalue, method = 'fdr', n=length(result_vol1$pvalue)) ## 5*8


##############################################################################
### Ana 7 again(1680 sub),  GHRcluster1Lcaudatemm3~ stress * R_putamenVol * PRS+age+sex+TIV+siteID
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda')# 1680*45

stressin <- match('SUMp_ksads_ptsd', colnames(ROIvolROIvarAll)) # LESN begins at which col;

rm(datause)
datause <- ROIvolROIvarAll
colnames(ROIvolROIvarAll)

## ### 7.1 again, GHRcluster1 ~ stressi*R_putamenVol*PRS0001
library(Matrix)
rm(fit1)
rm(fit1CI)
result1 <- c()
# stress loop，12
for (i in stressin){
  fit1 <- lm(GHRcluster1Lcaudatemm3 ~ datause[, i]* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
  fit1coef<-coef(summary(fit1))
  #calculate 95% CI for regression coefficient beta--> add to fit1
  fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
  
  result1<-rbind(result1, c(colnames(datause)[i], fit1CI[nrow(fit1coef), ]))
}
result1 <- as.data.frame(result1) # 7*7
colnames(result1)[5]<- c('pvalue')
result1$pvalue.Adj<- p.adjust(result1$pvalue, method = 'fdr', n=length(result1$pvalue)) ## 7*8, no FDR sign
GHRclu1_LoopStressRputaPRS <- result1
save(GHRclu1_LoopStressRputaPRS, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Cauda_SressPutamenPRS1680/V3Cat12ROIvolTIV_t1qc1680_GHRclu1_lm_LoopStressRputaPRS.rda') 



## ### 7.1 again again, binary stress score for samel above 7.1 again ana;
rm(list=ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda')# 1680*45

rm(datause)
datause <- subset(ROIvolROIvarAll, select = -c(ple_p_ss_affected_bad_mean, ple_y_ss_affected_bad_mean))

# (3 again bin) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
# ------ 取stress为binary SUMp_ksads_ptsd------
datause[, 14]<- ifelse(datause[, 14] <= 0, 0, 1) # binary

fit1 <- lm(GHRcluster1Lcaudatemm3 ~ SUMp_ksads_ptsd* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
fit1coef<-coef(summary(fit1))
#calculate 95% CI for regression coefficient beta--> add to fit1
fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
fit1CI



## ----- Analysis 8 -------
####### 240429--(补充分析，1680被试去除severity score>6 被试)-------
# 依据study Psychosis Risk Screening with the Prodromal Questionnaire – Brief version (PQ-B).Schizophr Res. 2011. 
##### PRS cutoff1（median+1SD）, R_putamen median separated Q4 Ana
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar.rda') # 1680*45

# Binarize column ksads_ptsd
ROIvolROIvarAll$SUMp_ksads_ptsdBin <- ifelse(ROIvolROIvarAll$SUMp_ksads_ptsd <= 0, 0, 1) # 1680*46

# 移除severity score>6
length(which(ROIvolROIvarAll$pps_y_ss_severity_score>6)) #200
indexL0 <- which(ROIvolROIvarAll$pps_y_ss_severity_score>6)

ROIvolROIvarAll_Repsy6 <- ROIvolROIvarAll[-indexL0, ]# 1480*46, using below;
save(ROIvolROIvarAll_Repsy6, file= '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar1480.rda')


volcutoff <- median(ROIvolROIvarAll_Repsy6$Rput78colmm3) # 4774.9701
PRSmedian <- median(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001, -0.002217855
PRSmean <- mean(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001, -0.0022795
PRSSD <- sd(ROIvolROIvarAll_Repsy6$PRS0001) # PRS0001 sd, 0.0003843


## ---- (8.1) 1480 PRS cutoff 1:median+SD;volcutoff:median Q4;
# ------ ---------- 8.1 Ana start -----------------
PRScutoff <- PRSmedian+PRSSD# +1SD, -0.001833;

dataQudrant1 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001>=PRScutoff), ] # 65*46
dataQudrant2 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3>=volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001<PRScutoff), ] # 675*46
dataQudrant3 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3<volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001<PRScutoff), ] # 665*46
dataQudrant4 <- ROIvolROIvarAll_Repsy6[(ROIvolROIvarAll_Repsy6$Rput78colmm3<volcutoff) & (ROIvolROIvarAll_Repsy6$PRS0001>=PRScutoff), ] # 75*46
save(dataQudrant1, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/t1qc1480_median1SD_Q1_s65ROIdata.rda')#
save(dataQudrant2, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/t1qc1480_median1SD_Q2_s675ROIdata.rda')#
save(dataQudrant3, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/t1qc1480_median1SD_Q3_s665ROIdata.rda')#
save(dataQudrant4, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/t1qc1480_median1SD_Q4_s75ROIdata.rda')#

#### ---- (8.1.1) ROI vol -LESN Bin 相关性；1480 Q4 subjects;
##---- Fit a linear regression model 
# ---------------with group variable and covariates
rm(fitvol4CI)
model4 <- lm(GHRcluster1Lcaudatemm3 ~ SUMp_ksads_ptsdBin + age +sex+catTIV+siteID, data = dataQudrant4)
#calculate 95% CI for regression coefficient beta--> add to fitvol1
fitvol4CI <- cbind(coef(summary(model4)), confint(model4, conf.level = 0.95))
## same PIC for 4 qudrant to manuscript;--- final--20240508



#### ---- (8.1.2) ROI vol -LESN 相关性；1480 划分的Q4 subjects;
stressin <- match('SUMp_ksads_ptsd', colnames(dataQudrant4)) # LESN begins at which col; 14

# (1). ROI（'GHRcluster1Lcaudatemm3'） vol 与tress相关性# -Q4被试；
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
result_vol1 <- as.data.frame(result_vol1) # 7*7
colnames(result_vol1)[5]<- c('pvalue')
result_vol1$pvalue.Adj<- p.adjust(result_vol1$pvalue, method = 'fdr', n=length(result_vol1$pvalue)) ## 7*8

save(result_vol1, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/SoOn_t1qc1480_Q4median1SD_lm75_LESN_GHRclu1Vol.rda')#


# --- (8.1.2 again, stress binary). 
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
result_vol1 <- as.data.frame(result_vol1) # 5*7
colnames(result_vol1)[5]<- c('pvalue')
result_vol1$pvalue.Adj<- p.adjust(result_vol1$pvalue, method = 'fdr', n=length(result_vol1$pvalue)) ## 5*8

save(result_vol1, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/1480_PRS0001cat12RputaVol_dataQ4x/SoOn_t1qc1480_Q4median1SD_lm75_LESNbin_GHRclu1Vol.rda')#
###### ------------ End Ana 8.1 (1480 Q4 分析) ------------


###　--------- Ana 8.2：interact start -------------
#load 1480*46 (while remove 200 sub(whose severityScore >6) from sub1680)
rm(list = ls())
load('~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/PCAgu_whiteA3086_V3Cat12ROIvolTIV_t1qc1680_interstVar1480.rda') # 1480*46

stressin <- match('SUMp_ksads_ptsd', colnames(ROIvolROIvarAll_Repsy6)) # LESN 
rm(datause)
datause <- ROIvolROIvarAll_Repsy6

## --- 8.2.1: ROI vol ~ stressi*R_putamenVol*PRS0001 -------
# (1) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
library(Matrix)
rm(fit1)
rm(fit1CI)
result1 <- c()
# stress loop，14-20
for (i in stressin){
  fit1 <- lm(GHRcluster1Lcaudatemm3 ~ datause[, i]* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
  fit1coef<-coef(summary(fit1))
  #calculate 95% CI for regression coefficient beta--> add to fit1
  fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
  
  result1<-rbind(result1, c(colnames(datause)[i], fit1CI[nrow(fit1coef), ]))
}
result1 <- as.data.frame(result1) # 7*7
colnames(result1)[5]<- c('pvalue')
result1$pvalue.Adj<- p.adjust(result1$pvalue, method = 'fdr', n=length(result1$pvalue)) ## 7*8, no FDR sign

GHRclu1_LoopStressRputaPRS <- result1
save(GHRclu1_LoopStressRputaPRS, file = '~/fu4_2yAnalysis/fu2year_V3Cat12regionalVolTIV/Further_removePSY6_1480_interestVar/ROIvol_SressPutamenPRS1480/SoOn_1480_GHRclu1_lm_LoopStressRputaPRS.rda') 


# --- (1 again bin) GHRcluster1 ~ stressi*R_putamenVol*PRS0001
# ------ sub 1480 取stress为binary SUMp_ksads_ptsd------
# datause, 1480*44
datause[, 14]<- ifelse(datause[, 14] <= 0, 0, 1) # binary

fit1 <- lm(GHRcluster1Lcaudatemm3 ~ SUMp_ksads_ptsd* Rput78colmm3 *PRS0001 + age + sex + catTIV+ siteID+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data = datause)
fit1coef<-coef(summary(fit1))
#calculate 95% CI for regression coefficient beta--> add to fit1
fit1CI <- cbind(fit1coef, confint(fit1, conf.level = 0.95))
fit1CI



