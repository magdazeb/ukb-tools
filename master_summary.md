# Summary of Master 6

```r
library(dplyr)
master6 = readRDS('Jinhee Code/master_6.rds')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')
dim(master6)
```
> [1] 4234599       9

```r
head(master6)
```
```
      eid Disease_Code ICD9.ICD10 Disease_Date   content   birthday
1 1000010        M7916      ICD10   2014-10-31 Secondary 1941-07-01
2 1000010         Z861      ICD10   2013-07-24 Secondary 1941-07-01
3 1000010         R680      ICD10   2013-11-02      Main 1941-07-01
4 1000010         Z433      ICD10   2001-09-12      Main 1941-07-01
5 1000010         R060      ICD10   2014-10-31 Secondary 1941-07-01
6 1000010         J181      ICD10   2014-06-07      Main 1941-07-01
  Diagnosed_age
1      73.38356
2      72.11233
3      72.38904
4      60.24110
5      73.38356
6      72.98356
                                                               meaning gender
1                                           M79.16 Myalgia (Lower leg) Female
2          Z86.1 Personal history of infectious and parasitic diseases Female
3 R68.0 Hypothermia, not associated with low environmental temperature Female
4                                         Z43.3 Attention to colostomy Female
5                                                       R06.0 Dyspnoea Female
6                                   J18.1 Lobar pneumonia, unspecified Female
```

```r
unique(master6$eid) %>% length
```
> [1] 457982

Split data by ICD version

```r
master6_icd9 = subset(master6,ICD9.ICD10=='ICD9')
master6_icd10 = subset(master6,ICD9.ICD10=='ICD10')

unique(master6_icd9$Disease_Code) %>% length
unique(master6_icd10$Disease_Code) %>% length
```
> [1] 3337
> [1] 11726

```r
unique(master6_icd9$eid) %>% length
unique(master6_icd10$eid) %>% length
```
> [1] 20303
> [1] 410336

```r
dis_10_freq = table(master6_icd10$Disease_Code) %>% as.data.frame
dis_10_freq_0 = subset(dis_10_freq,Freq>0)
dis_10_freq_200 = subset(dis_10_freq, Freq>200)

png('fig/hist_icd10_dis.png')
hist(dis_10_freq_0$Freq, main="ICD-10 disease records", breaks=20, xlab="Record number", ylab="Disease number")
dev.off()

png('fig/hist_icd10_dis_200.png')
hist(dis_10_freq_200$Freq, main="ICD-10 disease records (>200)", breaks=20, xlab="Record number", ylab="Disease number")
dev.off()
```


```R
dis_9_freq = table(master6_icd9$Disease_Code) %>% as.data.frame
dis_9_freq_0 = subset(dis_9_freq, Freq>0)
dis_9_freq_200 = subset(dis_9_freq, Freq>200)

png('fig/hist_icd9_dis.png')
hist(dis_9_freq_0$Freq, main="ICD-9 disease records", breaks=20, xlab="Record number", ylab="Disease number")
dev.off()

png('fig/hist_icd9_dis_200.png')
hist(dis_9_freq_200$Freq, main="ICD-9 disease records (>200)", breaks=20, xlab="Record number", ylab="Disease number")
dev.off()
```

Filter diseases only have >200 records

```R
master6_icd9_200 = subset(master6_icd9, Disease_Code %in% dis_9_freq_200$Var1)
master6_icd10_200 = subset(master6_icd10, Disease_Code %in% dis_10_freq_200$Var1)

unique(master6_icd9_200$eid) %>% length
unique(master6_icd10_200$eid) %>% length
```
> [1] 13064
> [1] 407519


Age-of-onset distribution: records and patients

Example 1. Alzheimer's disease (G30)

```R
#icd9_alz = dis_freq(master6_icd9,'290.11',200)

source('src/ard.r')
out_dir = 'fig/icd10_person'

icd10_alz = original_person(
	master6_icd10, code_ann=ICD10_Code_ann, code_num='G30', freq=200)

plot_person(icd10_alz, age_min=40, age_max=80, out_dir)
```


