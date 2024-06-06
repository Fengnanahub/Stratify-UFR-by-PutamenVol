
# install.packages("ggplot2")
require(UsingR)
require(ggplot2)


# read csv files
UKBIDlistweb <- read.table('~/UKB_FieldID_SubsetDown.csv', head = TRUE, sep=",")
head(UKBIDlistweb)
UKBIDneedmatch <- read.table('~/dataNeedbySelf1.csv', head = TRUE, sep=",")
# remove nan row
colnames(UKBIDneedmatch)

dimweb <- dim(UKBIDlistweb) # 25583*2
dimneedmatch <- dim(UKBIDneedmatch) # 1204*4

## ----------- match period  --------------
# eg. only 23321-2.0, 23321-3.0
for (i in 1:dimneedmatch[1]){
  newIDneedmatch <- paste(UKBIDneedmatch[i, 2], '-', sep = '')
  rm(rowfind)
  rowfind = grep(as.character(newIDneedmatch), UKBIDlistweb[, 1])
  if (length(rowfind)==0){
    print('Rowfind is null')
  }else{
    for (j in 1: length(rowfind)){
      UKBIDneedmatch[i, 3] <- UKBIDlistweb[rowfind[1], 2] # subset number;
      if (UKBIDlistweb[rowfind[j], 1] == as.character(paste(newIDneedmatch, '0.0', sep = ''))){
        UKBIDneedmatch[i, 4] <- UKBIDlistweb[rowfind[j], 1]
      } else if (UKBIDlistweb[rowfind[j], 1] == as.character(paste(newIDneedmatch, '1.0', sep = ''))){
        UKBIDneedmatch[i, 5] <- UKBIDlistweb[rowfind[j], 1]
      }else if (UKBIDlistweb[rowfind[j], 1] == as.character(paste(newIDneedmatch, '2.0', sep = ''))){
        UKBIDneedmatch[i, 6] <- UKBIDlistweb[rowfind[j], 1]
      }else if (UKBIDlistweb[rowfind[j], 1] == as.character(paste(newIDneedmatch, '3.0', sep = ''))){
        UKBIDneedmatch[i, 7] <- UKBIDlistweb[rowfind[j], 1]
      }
    }
  }
}
colnames(UKBIDneedmatch)[4:7]<- c('VarIDBase', 'VarIDfu1', 'VarIDfu2', 'VarIDfu3') # 1204*7
UKBIDneedmatch_period <- UKBIDneedmatch
write.csv(UKBIDneedmatch_period, file = "~/dataNeedbySelf1_matchWebPeriod.csv", row.names =FALSE)


##### --------------------- combine ROI var data start -------------------------
###### -----------1. select interest Var data in UKB subset csv -----------------
# load VarID need extracting data; UKB_subset_5
library(stringr)
UKB_ROIvarID <- read.table('~/UFR_CluCaudateValidate_Ana/ROIVarID.csv', head = TRUE, sep=",") 

# obtain varID data for each subset
subsetN <- unique(UKB_ROIvarID$subset) # 5,1,3,7,6,10,8,4

## ------ subset5 ---------- 
subsetID5 <- UKB_ROIvarID[UKB_ROIvarID$subset==5, ] # varID need extract in subset num 
UKBsubset5<- read.csv("~/UFR_CluCaudateValidate_Ana/UKB_subset_5.csv") # 502409*2501

# X22009.0.1, X22009.0.2,...X22009.0.20
PC_name = paste0("X22009.0.", 1:20)

# select VarIDbase col, ID, age, race, gene pc20; 
dataset5ROI <- UKBsubset5[, c('eid', 'X21022.0.0', 'X21000.0.0', PC_name)] 

colnames(dataset5ROI)[c(2:ncol(dataset5ROI))] <- c('age', 'race', 'BMI', paste0("genePC", 1:20)) # 502409

## similiarly select data: 
# ID, edu (subset3)
# ID, standard PRS-SZ (subset10)
# ID, ICD-10 (subset8)
# ID, traumatic life events (subset4)， replace -818,-121 as NA;
# ID, sex(0, Female; 1, male) (subset 1)...


## ------- 2. combine the subset Var data extracted above -----------
library(dplyr)
inter <-  c('eid')
UKB_ROIdata <- dataset7ROI %>% full_join(dataset6ROI, by=inter)# 502409
UKB_ROIdata <- dataset4ROInoNA %>% full_join(UKB_ROIdata, by=inter)# 502409
UKB_ROIdata <- dataset8ROI %>% full_join(UKB_ROIdata, by=inter)# 502409
UKB_ROIdata <- dataset10ROI %>% full_join(UKB_ROIdata, by=inter)# 502409
UKB_ROIdata <- dataset3ROI %>% full_join(UKB_ROIdata, by=inter)# 502409
UKB_ROIdata <- dataset1ROI %>% full_join(UKB_ROIdata, by=inter)# 502409
UKB_ROIdata <- dataset5ROI %>% full_join(UKB_ROIdata, by=inter)# 502409
save(UKB_ROIdata, file = '~/UFR_CluCaudateValidate_Ana/datasetROI_CombineALL.rds') 


#############################   start  ############################
## ------ 3. add cat12 processed TIV and volume by merge -----########
# (3.1) add TIV 
inter <-  c('eid')
UKB_ROIdata_CAT1 <- merge(UKB_ROIdata,TIV, by=inter) #39664*87

# (3.2) add striatum subarea volume 
UKB_ROIdata_CAT2 <- merge(UKB_ROIdata_CAT1, volume_unique, by='eid') #39657*94


##  -------------     select     --------------------------------
## 1. selecting sub to validating ana after Combine all,using 39657*94
# (0) Remove all in demo (eid, age,race,BMI,sex,edu,income,site,handness) 
listNA <- which(rowSums(is.na(UKB_ROIdata_CAT2[, c(1:4, 25:29)]))>8) # empty;

# (1) remove all NA in column of traumatic events--> 27350*94
# "ChildtraumnaTotal" "AdultTraumnaTotal" "TraumnaTotalall";
listNA <- which(rowSums(is.na(UKB_ROIdata[, c(61:63)]))>2) # find the row contain all NA,12307;
UKB_ROIdata_noNAtraum <- UKB_ROIdata[-listNA,] # get the data of no all NA row, 27350*94

# (2) further remove NA in colume of PRS-SZ --> 26661*94
# "stanPRSSZ" 
listNA <- which(rowSums(is.na(UKB_ROIdata_noNAtraum[, c(31:32)]))>1) # find the row contain all NA,689;
UKB_ROIdata_noNAtraum_noNAPRS <- UKB_ROIdata_noNAtraum[-listNA,] # get the data of no all NA row, 26661

# (3) further remove all NA in colume of CAT12 volume(mm3)-->26586
# "TIVml"   "GMml"   "WMml"   "CSFml"  
listNA <- which(rowSums(is.na(UKB_ROIdata_noNAtraum_noNAPRS[, c(84:87)]))>3) # find the row contain all NA, 0;
UKB_ROIdata_noNAtraum_PRS_catvol <- UKB_ROIdata_noNAtraum_noNAPRS # 26661

# (4) further remove NA in colume of CAT12 volume(mm3)
# Lput77volml"  "Rput78colml" "GHRcluster1Lcaudateml" 
listNA <- which(rowSums(is.na(UKB_ROIdata_noNAtraum_PRS_catvol[, c(88:94)]))>2) # find the row contain all NA, 0;

## ---------- factor sex,income...
UKB_ROIdata_noNAtraum_PRS_catvol$sex = factor(UKB_ROIdata_noNAtraum_PRS_catvol$sex)
UKB_ROIdata_noNAtraum_PRS_catvol$income = factor(UKB_ROIdata_noNAtraum_PRS_catvol$income)
UKB_ROIdata_noNAtraum_PRS_catvol$site = factor(UKB_ROIdata_noNAtraum_PRS_catvol$site)

# ----------- race<0 and handness<0 as NA;
UKB_ROIdata_noNAtraum_PRS_catvol$race <- ifelse(UKB_ROIdata_noNAtraum_PRS_catvol$race < 0, NA, UKB_ROIdata_noNAtraum_PRS_catvol$race)

UKB_ROIdata_noNAtraum_PRS_catvol$handness <- ifelse(UKB_ROIdata_noNAtraum_PRS_catvol$handness < 0, NA, UKB_ROIdata_noNAtraum_PRS_catvol$handness)

## race, handness as factor
UKB_ROIdata_noNAtraum_PRS_catvol$race = factor(UKB_ROIdata_noNAtraum_PRS_catvol$race)
UKB_ROIdata_noNAtraum_PRS_catvol$handness = factor(UKB_ROIdata_noNAtraum_PRS_catvol$handness)
## ------------ combine- org-select END--------------------




#############   Ana: merge combine-select3 26661-->24469 white  british ###############
## select white british
UKB_ROIdata_noNAtraum_PRS_catvol_white <- UKB_ROIdata_noNAtraum_PRS_catvol[UKB_ROIdata_noNAtraum_PRS_catvol$race %in% c(1001),] # 24469*94，British write;

# --- add SZ diagnosis ------ 
## search scz diagnosis in ICD10
library(readxl)
ICD10coding19 <- read_excel("~/UFR_CluCaudateValidate_AnaAgain/ICD10coding19.xlsx")
names_ICD10coding<-ICD10coding19$coding

names_ICD10_SCZ<-names_ICD10coding[2870:2879]

library(magrittr)
library(dplyr)
# get SCZ diagnosis
UKB_ROIdata_SCZ<- UKB_ROIdata_noNAtraum_PRS_catvol_white_TrauEvent_lifetime %>% 
  rowwise() %>% mutate(SCZ_diagnosis=
                         ifelse(sum(as.numeric(paste0("'",names_ICD10_SCZ,"'") %in% unlist(strsplit(gsub("\\[|\\]", "", ICD10),", "))))>0,TRUE,FALSE)) # 24469

## ICD10 NA，SCZ diag return NA；
UKB_ROIdata_SCZ1$SCZ_diagnosis <- ifelse(UKB_ROIdata_SCZ1[, c('ICD10')]=="", NA, UKB_ROIdata_SCZ1$SCZ_diagnosis) 
b<- UKB_ROIdata_SCZ1[, c('ICD10', 'SCZ_diagnosis')]


## --- 以standard PRS-SZ median+1sd, 以及Rputamen median 划分四象限----- ###
## select individuals without SZ diagnosis
UKB_ROIdata_SCZ1noSCZ <- UKB_ROIdata_SCZ1[UKB_ROIdata_SCZ1$SCZ_diagnosis==FALSE & !is.na(UKB_ROIdata_SCZ1$SCZ_diagnosis), ] # 20547*113, 40-70y;  age mean[sd] 55.38[7.43]; 9387 male, 11160 female, 45.69% male ;

volcutoff <- median(UKB_ROIdata_SCZ1noSCZ$Rput78colml) # 3.952535
PRSmedian <- median(UKB_ROIdata_SCZ1noSCZ$stanPRSSZ) # stanPRS-sz, -0.363557
PRSSD <- sd(UKB_ROIdata_SCZ1noSCZ$stanPRSSZ) # stanPRS-sz sd, 0.9910816

## ###### PRS cutoff 1；median+SD
PRScutoff <- PRSmedian+PRSSD# +1SD, 0.6275246;

ROIvolROIvarAll <- UKB_ROIdata_SCZ1noSCZ
dataQudrant1 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colml>=volcutoff) & (ROIvolROIvarAll$stanPRSSZ>=PRScutoff), ] # 1535
dataQudrant2 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colml>=volcutoff) & (ROIvolROIvarAll$stanPRSSZ<PRScutoff), ] # 8739
dataQudrant3 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colml<volcutoff) & (ROIvolROIvarAll$stanPRSSZ<PRScutoff), ] # 8642
dataQudrant4 <- ROIvolROIvarAll[(ROIvolROIvarAll$Rput78colml<volcutoff) & (ROIvolROIvarAll$stanPRSSZ>=PRScutoff), ] # 1631


### (1). ROI（'GHRcluster1Lcaudatemm3'） vol - adult stress；
library(Matrix)
result_vol1 <- c()
rm(fitvol1)
dataQudrant <- dataQudrant4
for (i in c(62)){
  ## ROI vol~stressi+...
  fitvol1 <- lm(GHRcluster1Lcaudateml ~ as.numeric(unlist(dataQudrant[, i]))+ age + sex + TIVml+ site, data = dataQudrant)
  #calculate 95% CI for regression coefficient beta--> add to fitvol1
  fitvol1CI <- cbind(coef(summary(fitvol1)), confint(fitvol1, conf.level = 0.95))
  
  result_vol1<-rbind(result_vol1, c(colnames(dataQudrant)[i], fitvol1CI[2, ]))
}
result_vol1 <- as.data.frame(result_vol1) # 1*7
colnames(result_vol1)[5]<- c('pvalue')
result_vol1$pvalue.Adj<- p.adjust(result_vol1$pvalue, method = 'fdr', n=length(result_vol1$pvalue)) ## 1*8


## PIC1---- GHRcluster1Lcaudateml~AdultTraumnaTotal
library(ggplot2)
PIC1 <- ggplot() + geom_point(data = dataQudrant1, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), color = "#b6D7E9", fill = "white", shape = 1, size = 2.8, stroke = 0.7) + geom_smooth(data = dataQudrant1, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), method = "lm", size = 0.9, color = "#82A7d1", linetype = "dotted") + labs(x="AdultTraumnaTotal", y = "GHRcluster1Lcaudateml")  +theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(labels = scales::scientific) 
PIC2 <- ggplot() + geom_point(data = dataQudrant2, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), color = "#b6D7E9", fill = "white", shape = 1, size = 2.8, stroke = 0.7) + geom_smooth(data = dataQudrant2, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), method = "lm", size = 0.9, color = "#82A7d1", linetype = "dotted") + labs(x="AdultTraumnaTotal", y = "GHRcluster1Lcaudateml")  +theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(labels = scales::scientific) 
PIC3 <- ggplot() + geom_point(data = dataQudrant3, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), color = "#b6D7E9", fill = "white", shape = 1, size = 2.8, stroke = 0.7) + geom_smooth(data = dataQudrant3, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), method = "lm", size = 0.9, color = "#82A7d1", linetype = "dotted") + labs(x="AdultTraumnaTotal", y = "GHRcluster1Lcaudateml")  +theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(labels = scales::scientific) 
PIC4 <- ggplot() + geom_point(data = dataQudrant4, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), color = "#b6D7E9", fill = "white", shape = 1, size = 2.8, stroke = 0.7) + geom_smooth(data = dataQudrant4, aes(x=AdultTraumnaTotal, y=GHRcluster1Lcaudateml), method = "lm", size = 0.9, color = "#82A7d1", linetype = "dotted") + labs(x="AdultTraumnaTotal", y = "GHRcluster1Lcaudateml")  +theme_bw() + theme(panel.grid = element_blank()) + scale_x_continuous(labels = scales::scientific)

library(patchwork)
(PIC2+PIC1)/(PIC3+PIC4)
