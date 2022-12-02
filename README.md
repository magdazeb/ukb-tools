# I. UK Biobank tools

How to use:

* Run help command `Rscript ukb-exe.r -h`

```cmd
ukb-exe v0.1 (2020-01-12)

Usage: Rscript ukb-exe.r [functions] [options] <-b> <-t>

Functions:
    --ukb_prune | -pr  This is a first function for automated pruning for UKB data.
    --ukb_excld | -ex  This is a second function for select eid for further 
      		       analysis, such as correlation table and linear regression
                       model.
  
Global arguments:
    --base      | -b   Base RDS files are mendatory.
    --target    | -t   Target_dir or target_file is mendatory.
    --verbose   | -v   Rich description for debugging. Default is FALSE.
  
Running functions without arguments prints usage information for [functions].
```



## 1. ukb_prune

* Run help command `Rscript ukb-exe.r --ukb_prune -h`

```cmd
ukb-exe v0.1 (2020-01-12)

Usage: Rscript ukb-exe.r --ukb_prune [options] <--base|-b> <--target|-t>
    --ukb_prune | -pr  This is a first function for automated pruning for UKB data.

Required arguments:
    -b <1.rds 2.rds 3.rds>
                       Three base RDS files are needed: 
                       * 1.rds: A RDS file for UKB phenotypes. This is a
                           directly downloaded file from the UK biobank.
                       * 2.rds: A RDS file for UKB field ID description.
                       * 3.rds: A RDS file for UKB field ID code information.
                       Please keep the file order!
    -t <target_dir target_file.rds>
                       Two variables are needed for target directory/RDS file:
                       * target_dir: A target directory address for save figures
                       * target_file.rds: A target file for result phenotype table.
                       Please keep the dir/file order!

Either one of these are required:
    --fids <fid1 fid2 ...>
                       Fid numbers, a subset of field IDs to extract reliable
                       subjects.
    --anns <Reproductive_Aging / Longevity / Confounders>
                       Either Reproductive_Aging/Longevity/Confounders,
                       an annotation for a subset of field IDs in UKB field
                       description file (e.g., db_UKB/ukb_id_info.csv.rds),
                       including Reproductive_Aging, Longevity, and Confounders

Optional arguments:
    --expt             Fid numbers, a subset of field IDs as exception for filtering.
```

Usage examples for ukb_prune:

1. Longevity

```cmd
#Rscript ukb-exe.r -pr -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_lon data_UKB/ukb_prune_lon.rds -anns Longevity -expt 1845 2946 40007 > data_UKB/ukb_prune_lon.rds_log.txt
Rscript ukb-exe.r -pr -v \
	-b db_menopause/ukb44928.csv.rds db_menopause/ukb_fid_info.csv.rds db_menopause/ukb_fid_code.csv.rds \
	-t fig_menopause/prune_lon data_menopause/ukb_prune_lon.rds \
	-anns Longevity \
	-expt 1845 2946 40007 \
> data_menopause/ukb_prune_lon.rds_log.txt
```

* input RDS files
  * `db_menopause/ukb44928.csv.rds`: ukb phenotype data (downloaded from UKB and converted to RDS file)
  * `db_menopause/ukb_fid_info.csv.rds`: ukb field ID information file (manually generated)
  * `db_menopause/ukb_fid_code.csv.rds`: ukb field ID code information file (manually generated)
* `-anns Longevity`: Select a fid preset, 8 phenotypes
* `-expt 1845 2946 40007`: exceptionally fids
  * 1845: Mother's age, values are keep increasing since they are alive, Just select latest data
  * 2946: Father's age, values are keep increasing since they are alive, Just select latest data
  * 40007: Age at death, answer difference distribution is not following normal distribution. And the differences are less than -0.05~0.05 years range, which is likely very small.
* `> data_UKB/ukb_prune_lon.rds_log.txt` save messages to a designated file
* Job done: 2021-03-30 02:52:08 for 2.9 min

2. Reproductive aging

```cmd
#Rscript ukb-exe.r -pr -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_rep data_UKB/ukb_prune_rep.rds -anns Reproductive_Aging > data_UKB/ukb_prune_rep.rds_log.txt
Rscript ukb-exe.r -pr -v \
	-b db_menopause/ukb44928.csv.rds db_menopause/ukb_fid_info.csv.rds db_menopause/ukb_fid_code.csv.rds \
	-t fig_menopause/prune_rep data_menopause/ukb_prune_rep.rds \
	-anns Reproductive_Aging \
> data_menopause/ukb_prune_rep.rds_log.txt
```

* `-anns Reproductive_Aging`: Select a fid preset, 14 phenotypes
* `> data_UKB/ukb_prune_rep.rds_log.txt`: save messages to a designated file
* Job done: 2021-03-30 03:10:02 for 4.3 min

3. Confounders

```cmd
#Rscript ukb-exe.r -pr -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_con data_UKB/ukb_prune_con.rds -anns Confounders > data_UKB/ukb_prune_con.rds_log.txt
Rscript ukb-exe.r -pr -v \
	-b db_menopause/ukb44928.csv.rds db_menopause/ukb_fid_info.csv.rds db_menopause/ukb_fid_code.csv.rds \
	-t fig_menopause/prune_con data_menopause/ukb_prune_con.rds \
	-anns Confounders \
> data_menopause/ukb_prune_con.rds_log.txt
```

* `-anns Confounders`: Select a fid preset, 77 phenotypes
* `> data_UKB/ukb_prune_con.rds_log.txt`: save messages to a designated file
* Job done: 2021-03-30 03:39:34 for 25.9 min

4. Additional

* fid1767 Adopted as a child

```cmd
#Rscript ukb-exe.r -pr -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_con2 data_UKB/ukb_prune_con2.rds --fids 1767 > data_UKB/ukb_prune_con2.rds_log.txt
Rscript ukb-exe.r -pr -v \
	-b db_menopause/ukb44928.csv.rds db_menopause/ukb_fid_info.csv.rds db_menopause/ukb_fid_code.csv.rds \
	-t fig_menopause/prune_con2 data_menopause/ukb_prune_con2.rds \
	--fids 1767 \
> data_menopause/ukb_prune_con2.rds_log.txt
```

* fid21022 Age at recruitment
   
```cmd
#Rscript ukb-exe.r -pr -v -b db_UKB/ukb25461.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_con2 data_UKB/ukb_prune_con2_21022.rds --fids 21022 > data_UKB/ukb_prune_con2_21022.rds_log.txt
Rscript ukb-exe.r -pr -v \
	-b db_menopause/ukb25461.csv.rds db_menopause/ukb_fid_info.csv.rds db_menopause/ukb_fid_code.csv.rds \ #<- ukb25461?
	-t fig_menopause/prune_con2 data_menopause/ukb_prune_con2_21022.rds \
	--fids 21022 \
> data_menopause/ukb_prune_con2_21022.rds_log.txt
```



## 2. ukb_excld

* Run help command `Rscript ukb-exe.r --ukb_prune -h`

```
ukb-exe v0.1 (2020-01-12)

Usage: Rscript ukb-exe.r --ukb_excld [options] <-b> <-t>
    --ukb_excld | -ex  This is a second function for select eid for further analysis,
                       such as correlation table and linear regression model.

Required arguments:
    -b <1.rds>         One base RDS file is needed:
                       * 1.rds: A RDS file for UKB phenotypes. This is a result file
                           from the ukb_prune function.
    -t <target_dir target_file.rds>
                       Two variables are needed for target directory/RDS file.
    --fids <fid1 fid2 fid3 ...>
                       Fid numbers, a subset of field IDs to extract subjects.
                       If there are three fids, this function will automatically
                       select eids having the no[0] answers from the first fid (data
                       type should be categorical). Then select union eids from
                       the second and third fids, which data types should be 
                       continuous or integer.
                       * ex1) setdiff(Oophorectomy_no, union(Oophorectomy_yes,
                           Oophorectomy_age))
                       * ex2) setdiff(HRT_no, union(HRT_yes, HRT_start, HRT_end))
                       * In case of multiple categorial data is not considered yet.

Optional arguments:
	--cutoff <num1 num2>
	                   If data type is continuous or integer, setting cutoff criteria
	                   such as age at menarche: <9 or >20 years.
	--cat_no           True/False argument for select "No" answers from categorial
	                   data.
```

Usage examples for ukb_excld:

1. Filtering Menarche

```cmd
#Rscript ukb-exe.r --ukb_excld -v -b data_UKB/ukb_prune_rep.rds --fids fid2714 --cutoff 9 20 -t data_UKB/eid_fid2714.csv > data_UKB/eid_fid2714.csv_log.txt
Rscript ukb-exe.r --ukb_excld -v \
   -b data_menopause/ukb_prune_rep.rds \
   --fids fid2714 \
   --cutoff 9 20 \
   -t data_menopause/eid_fid2714.csv \
> data_menopause/eid_fid2714.csv_log.txt
```

* fid2714 Age when periods started (menarche)
* pruned and filtered outliers by <9 or >20 years old
* messages saved at `data_UKB/eid_fid2714.csv_log.txt`

2. Filtering Menopause

```cmd
#Rscript ukb-exe.r --ukb_excld -v -b data_UKB/ukb_prune_rep.rds --fids fid3581 --cutoff 40 60 -t data_UKB/eid_fid3581.csv > data_UKB/eid_fid3581.csv_log.txt
Rscript ukb-exe.r --ukb_excld -v \
   -b data_menopause/ukb_prune_rep.rds \
   --fids fid3581 \
   --cutoff 40 60 \
   -t data_menopause/eid_fid3581.csv \
> data_menopause/eid_fid3581.csv_log.txt
```

* fid3581 Age at menopause (last menstrual period)
* pruned and filtered outliers by <40 or >60 years old
* messages saved at `data_UKB/eid_fid3581.csv_log.txt`

3. setdiff Hysterectomy

```cmd
#Rscript ukb-exe.r --ukb_excld -v -b data_UKB/ukb_prune_rep.rds --fids fid3591 fid2824 --cat_no -t data_UKB/eid_fid3591.csv > data_UKB/eid_fid3591.csv_log.txt
Rscript ukb-exe.r --ukb_excld -v \
   -b data_menopause/ukb_prune_rep.rds \
   --fids fid3591 fid2824 \
   --cat_no \
   -t data_menopause/eid_fid3591.csv \
> data_menopause/eid_fid3591.csv_log.txt
```

* fid3591 Ever had hysterectomy (womb removed): Yes or No
  * fid2824 Age at hysterectomy: Integer, years
  * setdiff(fid3591 No, union(fid3591 Yes, fid2824) )
* messages saved at `data_UKB/eid_fid3591.csv_log.txt`

4. setdiff Bilateral oophorectomy

```cmd
#Rscript ukb-exe.r --ukb_excld -v -b data_UKB/ukb_prune_rep.rds --fids fid2834 fid3882 --cat_no -t data_UKB/eid_fid2834.csv > data_UKB/eid_fid2834.csv_log.txt
Rscript ukb-exe.r --ukb_excld -v \
   -b data_menopause/ukb_prune_rep.rds \
   --fids fid2834 fid3882 \
   --cat_no \
   -t data_menopause/eid_fid2834.csv \
> data_menopause/eid_fid2834.csv_log.txt
```

* fid2834 Bilateral oophorectomy (both ovaries removed): Yes or No
  * fid3882 Age at bilateral oophorectomy (both ovaries removed): Integer, years
  * setdiff(fid2834 No, union(fid2834 Yes, fid3882) )
* messages saved at `data_UKB/eid_fid2834.csv_log.txt`

5. setdiff Hormone-replacement therapy (HRT)

```cmd
#Rscript ukb-exe.r --ukb_excld -v -b data_UKB/ukb_prune_rep.rds --fids fid2814 fid3536 fid3546 --cat_no -t data_UKB/eid_fid2814.csv > data_UKB/eid_fid2814.csv_log.txt
Rscript ukb-exe.r --ukb_excld -v \
	-b data_menopause/ukb_prune_rep.rds \
   	--fids fid2814 fid3536 fid3546 \
   	--cat_no \
   	-t data_menopause/eid_fid2814.csv \
> data_menopause/eid_fid2814.csv_log.txt
```
   
* fid2814 Ever used hormone-replacement therapy (HRT)
  * fid3536, Age started hormone-replacement therapy (HRT)
  * fid3546, Age last used hormone-replacement therapy (HRT)
* messages saved at `data_UKB/eid_fid2814.csv_log.txt`
   
To identify the subjects who don't had hysterectomy, bilateral oophorectomy, and hormone-replacement therapy (HRT), exclude answers of ever had question is Yes as well as data existing in ages at hyterectomy, bilateral oophorectomy, and HRT.



## 3. Venn analysis

Run below codes in R.

eid list 1 - No Surgery/HRT

```R
library(dplyr)
dir = 'data_menopause'
eid_hyterectomy = read.csv(paste0(dir,'/eid_fid3591.csv')) %>% unlist
eid_oophorectomy = read.csv(paste0(dir,'/eid_fid2834.csv')) %>% unlist
eid_hormone = read.csv(paste0(dir,'/eid_fid2814.csv')) %>% unlist
eid_ectomy_li = list(eid_hyterectomy,eid_oophorectomy,eid_hormone)
names(eid_ectomy_li) = c('No hysterectomy','No oophorectomy','No HRT')
```

Draw a venn diagram

```R
source('src/venn_analysis.r')
venn_ectomy = venn_analysis(grouplist=eid_ectomy_li, dir='fig_menopause')
```

eid list 2 - Replicative lifespan

```R
eid_ectomy_inter = Reduce(intersect, eid_ectomy_li)
eid_menarche_flt = read.csv(paste0(dir,'/eid_fid2714.csv')) %>% unlist
eid_menopause_flt = read.csv(paste0(dir,'/eid_fid3581.csv')) %>% unlist

eid_rep_li = list(
    eid_menarche_flt,  # excluding age at menarche <9 or >20 years
    eid_menopause_flt, # excluding age at menopause <40 or >60 years
    eid_ectomy_inter
)
names(eid_rep_li) = c('Age at menarche','Age at menopause','No surgery/HRT')
```

Draw a venn diagram

```R
source('src/venn_analysis.r')
venn_rep = venn_analysis(grouplist=eid_rep_li, dir='fig_menopause')
```



### 4. Calculate replicative lifespan

Run below codes in R.

```R
library(dplyr)
library(ggplot2)
```

```R
ukb_rep  = readRDS('data_menopause/ukb_prune_rep.rds')
dim(ukb_rep)%>%print

eid_rep_lif = Reduce(intersect,eid_rep_li)
ukb_rep_lif_tmp = ukb_rep %>% 
    select(eid,starts_with('fid2714'),starts_with('fid3581') ) %>%
    filter(eid %in% eid_rep_lif)
dim(ukb_rep_lif_tmp) %>% print

ukb_rep_lifespan = data.frame(
    eid = ukb_rep_lif_tmp$eid,
    Rep_lifespan = ukb_rep_lif_tmp[,3]-ukb_rep_lif_tmp[,2]
)
summary(ukb_rep_lifespan$Rep_lifespan)
write.csv(ukb_rep_lifespan,'data_menopause/ukb_rep_lifespan.csv',row.names=F,quote=F)
```

Draw a distribution histogram

```R
ukb_rep_lifespan = read.csv('data_menopause/ukb_rep_lifespan.csv')
dim(ukb_rep_lifespan) %>% print
summary(ukb_rep_lifespan) %>% print

p_title = paste0('Replicative lifespan from ',nrow(ukb_rep_lifespan),' subjects.')
f_name  = 'fig_menopause/Rep_lifespan.png'
ukb_rep_lif = ukb_rep_lifespan
colnames(ukb_rep_lif) = c('eid','fid')

source('src/ukb_filter.r')
ukb_histo(
    dat    = ukb_rep_lif,  # Input data for histogram
    unit   = 'Replicative lifespan (years)',   # Unit for x-axis
    type   = 'Integer', # Data types*
    title  = p_title,   # Plot title
    f_name = f_name,    # Save file as name
    pruned = F          # Add decile information? default=T
)
# * Data types = Categorical single/Categorical multiple, Date, Integer, and Continuous
```







# II. Age-related diseases: ICD-10 code

## 1. Calculating ICD-10 background incidence

Preparing total incidence of ICD-10 data as background. Save the result as RDS file for later use.

```R
library(dplyr)

# Read data
master6 <- readRDS("Jinhee Code/master_6.rds")
master6_icd10 = subset(master6,ICD9.ICD10=='ICD10')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')

# Calculate the total incidence of ICD-10
source('src/ard.r')
total_icd10_freq = original(
    master6_icd10, 
    code_num=NULL, 
    ICD10_Code_ann, 
    totalage_freq=NULL, 
    subgroup=FALSE)[[1]]

# Save as RDS file
saveRDS(total_icd10_freq,'total_icd10_freq.rds')

# Create barplot of total ICD-10 disease incidences
total_icd10_freq_vec = total_icd10_freq$Frequency
names(total_icd10_freq_vec) = total_icd10_freq$`Age at Diagnosis`
png('fig/total_icd10_freq.png')
barplot(total_icd10_freq_vec,main='Total ICD-10 incidences',xlab='Age at Diagnosis')
dev.off()
```

ICD10 codes have 4.1 M records with age range 40-80 years old.

![](fig/icd10/total_icd10_freq.png)

### Calculating ICD-9 background incidence

```R
source('src/ard.r')
master6_icd9 = subset(master6,ICD9.ICD10=='ICD9')
total_icd9_freq = original(
    master6_icd9, 
    code_num=NULL, 
    code_ann=NULL,
    totalage_freq=NULL, 
    subgroup=FALSE)[[1]]
saveRDS(total_icd9_freq,'total_icd9_freq.rds')

# Create barplot of total ICD-10 disease incidences
total_icd9_freq_vec = total_icd9_freq$Frequency
names(total_icd9_freq_vec) = total_icd9_freq$`Age at Diagnosis`
png('fig/total_icd9_freq.png')
barplot(total_icd9_freq_vec, main='Total ICD-9 incidences', xlab='Age at Diagnosis')
dev.off()
```

ICD9 codes have 58 K (58,699) records with age range 27-57 years old.

![](fig/total_icd9_freq.png)

### Example: G35 Multiple sclerosis / E10 Type 1 diabetes

```R
library(dplyr)

# Read data
master6 <- readRDS("Jinhee Code/master_6.rds")
master6_icd10 = subset(master6,ICD9.ICD10=='ICD10')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')
total_icd10_freq = readRDS('total_icd10_freq.rds')

# Calculate incidence rates
source('src/ard.r')
dis = c("G35", "E10")
out_dir = 'fig/icd10'

for(i in 1:length(dis)) {
    # Prepare data
    inc_rate = original(master6,dis[i],ICD10_Code_ann,200,total_icd10_freq,TRUE)
    inc_norm = normalized(inc_rate)

    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=40,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=40,age_max=80,out_dir)
}
```

| Disease                                 | Incidence rate                        | Normalized incidence rate               |
| --------------------------------------- | ------------------------------------- | --------------------------------------- |
| G35 Multiple sclerosis                  | ![](fig/icd10/G35_original_log10.png) | ![](fig/icd10/G35_normalized_log10.png) |
| E10 Insulin-dependent diabetes mellitus | ![](fig/icd10/E10_original_log10.png) | ![](fig/icd10/E10_normalized_log10.png) |



### Examples: Age-related diseases

![](fig/ARD.png)

This is the plot from UKB self-reported data.

Draw incidence rate plot for major age-related diseases and its subgroups.

```R
source('src/ard.r')
dis = c("G30", "E11", "E14", "I21", "I22", "I50", "I64", "J44")
out_dir = 'fig/icd10'

for(i in 1:length(dis)) {
    # Prepare data
    inc_rate = original(master6,dis[i],ICD10_Code_ann,200,total_icd10_freq,TRUE)
    inc_norm = normalized(inc_rate)

    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=40,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=40,age_max=80,out_dir)
}
```

| Disease                                                      | Incidence rate                        | Normalized incidence rate               |
| ------------------------------------------------------------ | ------------------------------------- | --------------------------------------- |
| G30 Alzheimer's disease (AD)                                 | ![](fig/icd10/G30_original_log10.png) | ![](fig/icd10/G30_normalized_log10.png) |
| I50 Heart failure (CHF)                                      | ![](fig/icd10/I50_original_log10.png) | ![](fig/icd10/I50_normalized_log10.png) |
| J44 Other chronic obstructive pulmonary disease (COPD)       | ![](fig/icd10/J44_original_log10.png) | ![](fig/icd10/J44_normalized_log10.png) |
| I21 Acute myocardial infarction (Acute MI)                   | ![](fig/icd10/I21_original_log10.png) | ![](fig/icd10/I21_normalized_log10.png) |
| I22 Subsequent myocardial infarction (Subsequent MI)         | ![](fig/icd10/I22_original_log10.png) | ![](fig/icd10/I22_normalized_log10.png) |
| I64 Stroke, not specified as haemorrhage or infarction       | ![](fig/icd10/I64_original_log10.png) | ![](fig/icd10/I64_normalized_log10.png) |
| E11 Non-insulin-dependent diabetes mellitus (Diabetes Type II) | ![](fig/icd10/E11_original_log10.png) | ![](fig/icd10/E11_normalized_log10.png) |
| E14 Unspecified diabetes mellitus (Unspecified Diabetes)     | ![](fig/icd10/E14_original_log10.png) | ![](fig/icd10/E14_normalized_log10.png) |



## 2. Preprocessing for ICD-10

Calculating incidence of each ICD10 code.

```R
library(dplyr)

# Read data
#totalage_freq = readRDS('Jinhee Code/totalage_freq.rds')
master6          = readRDS("Jinhee Code/master_6.rds")
master6_icd10    = subset(master6,ICD9.ICD10=='ICD10')
ICD10_Code_ann   = read.delim('ICD10_DataCoding_41270.tsv')
total_icd10_freq = readRDS('total_icd10_freq.rds')
icd10 = master6_icd10$Disease_Code %>% as.character %>% unique %>% sort

# Calculate incidence rate of ICD10
## Minimum incidence criteria: 200
source('src/ard.r')
clust_result_ir = clustering_preprocess(
    master         = master6_icd10,
    code_num       = icd10,
    code_ann       = ICD10_Code_ann,
    totalage_freq  = total_icd10_freq,
    normalized     = FALSE)

saveRDS(clust_result_ir,'clust_result_icd10_ir.rds')
```

> Processing 11726 iterations:
>
> head(split)
>   100/11726 A392 Job process: 33 sec
>   200/11726 B000 Job process: 1.1 min
> ...
>   11600/11726 Z851 Job process: 1.1 hr
>   11700/11726 Z961 Job process: 1.1 hr
> List = 2089 (freq >200) ; Merging data = 41  2090 -> done
>
> Job done: 2021-02-22 02:07:25 for 1.1 hr

Calculate normalized incidence rate of ICD10

```R
# Run clustering_preprocess
## Minimum incidence criteria: 200
source('src/ard.r')
clust_result_nir = clustering_preprocess(
    master         = master6_icd10,
    code_num       = icd10,
    code_ann       = ICD10_Code_ann,
    totalage_freq  = total_icd10_freq,
    normalized     = TRUE)

saveRDS(clust_result_nir,'clust_result_icd10_nir.rds')
```

> 



## 3. Data pruning: Draw heatmap of ICD-10 normalized incidences

* [See R color palette](https://kbroman.files.wordpress.com/2014/05/crayons.png)
* `library(colorspace)` [ref](https://www.r-bloggers.com/2019/01/colorspace-new-tools-for-colors-and-palettes/), [fig](https://i2.wp.com/eeecon.uibk.ac.at/~zeileis/assets/posts/2019-01-14-colorspace/hcl-palettes-1.png?ssl=1)

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Read data
clust_result = readRDS('clust_result_icd10.rds')

# Prepare sub-sets
## clust_result:  41 2089; original data
## clust_result2: 41  873; peak age over 60
## clust_result3: 36  784; peak age over 60 & remove top 1 and bottom 2 rows (< 1,500 incidences)
source('src/ard.r')
clust_result2 = subsetting_cluster_result(clust_result, 0,0, 60)
clust_result3 = subsetting_cluster_result(clust_result, 4,2, 60)
col_fun = colorRamp2(c(0,1,2,3,4,7,10), c("white","Sky Blue","yellow Green","yellow","red","purple","black"))
```

Draw total 2,089 ICD-10 diseases heatmap

```R
png('fig/icd10/clust_result.png', width=25,height=8, units='in', res=150)
Heatmap(clust_result, cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/icd10/clust_result.png)

Draw 873 diseases filtered by peak age after 60 years

```R
png('fig/icd10/clust_result_peak60.png', width=25,height=8, units='in', res=150)
Heatmap(clust_result2, cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/icd10/clust_result_peak60.png)

Draw 784 diseases with peak age after 60 years and removed top 4 and bottom 2 rows (< 1,500 incidences)

```R
# set configurations
i = 3
wh = c(25,8)
col_split = c(8,12,20)
f_name = paste0('fig/icd10/clust_result_peak60_2,4trimed_',col_split[i],'.png')

# Draw heatmaps
png(f_name, width=25,height=8, units='in', res=150)
Heatmap(clust_result3, cluster_rows = FALSE, col = col_fun, column_split = col_split[i], column_dend_height=unit(1.5,'in'))
dev.off()

# Kill all graphics objects
graphics.off()
```

![](fig/icd10/clust_result_peak60_2,4trimed_20.png)



## 4. Extract clustered data & draw cluster plots

Load libraries & data

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

master6          = readRDS("Jinhee Code/master_6.rds")
master6_icd10    = subset(master6,ICD9.ICD10=='ICD10')
ICD10_Code_ann   = read.delim('ICD10_DataCoding_41270.tsv')
total_icd10_freq = readRDS('total_icd10_freq.rds')
```

Split normalized incidence rate data by hclust

```R
k = 20
f_name = paste0('fig/icd10/clust_result_peak60_2,4trimed_',k,'-re.png')
code_nm = 'ICD-10'
split = data.frame(cutree(hclust(dist(t(clust_result3))), k = k))
split[,2] = rownames(split)
colnames(split) = c('cluster',code_nm)

png(f_name, width=25,height=7, units='in', res=150)
Heatmap(clust_result3, column_split=split[,1], cluster_rows=FALSE, col=col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/icd10/clust_result_peak60_2,4trimed_20-re.png)

Draw correlation matrix

```R
library(corrplot)
cor_mat = cor(clust_result3, method="spearman")
f_name = paste0('fig/icd10/clust_result_peak60_2,4trimed_',k,'-corr.png')

png(f_name, width=25,height=25, units='in', res=150)
corrplot(cor_mat, method='square', tl.col="black", order="hclust", tl.cex=0.7)
dev.off()
```

![](fig/icd10/clust_result_peak60_2,4trimed_20-corr.png)

Prepare normalized incidence rate data list by cluster

```R
source('src/ard.r')
dis_cluster = normalized_by_cat(
    master         = master6_icd10,
    cat_code       = split,
    code_ann       = ICD10_Code_ann,
    totalage_freq  = total_icd10_freq,
    subgroup       = FALSE)
```

> ** Run normalized_by_cat **
>
> Processing 20 iterations:
>   1/20 1 = 317 -> dim  10306 4; Job process: 2 min
>   2/20 2 = 117 -> dim  3490 4; Job process: 2.7 min
>   3/20 3 = 2 -> dim  88 4; Job process: 2.7 min
>   4/20 4 = 41 -> dim  1335 4; Job process: 2.9 min
>   5/20 5 = 2 -> dim  101 4; Job process: 3 min
>   6/20 6 = 59 -> dim  2010 4; Job process: 3.3 min
>   7/20 7 = 145 -> dim  4662 4; Job process: 4.2 min
>   8/20 8 = 21 -> dim  683 4; Job process: 4.3 min
>   9/20 9 = 9 -> dim  314 4; Job process: 4.4 min
>   10/20 10 = 36 -> dim  1017 4; Job process: 4.6 min
>   11/20 11 = 8 -> dim  214 4; Job process: 4.7 min
>   12/20 12 = 6 -> dim  245 4; Job process: 4.7 min
>   13/20 13 = 2 -> dim  82 4; Job process: 4.7 min
>   14/20 14 = 5 -> dim  157 4; Job process: 4.8 min
>   15/20 15 = 3 -> dim  127 4; Job process: 4.8 min
>   16/20 16 = 3 -> dim  132 4; Job process: 4.8 min
>   17/20 17 = 3 -> dim  95 4; Job process: 4.8 min
>   18/20 18 = 2 -> dim  76 4; Job process: 4.9 min
>   19/20 19 = 2 -> dim  107 4; Job process: 4.9 min
>   20/20 20 = 1 -> dim  69 4; Job process: 4.9 min

Draw normalized incidence rate as line plot by cluster

```R
source('src/ard.r')
for (i in 1:length(dis_cluster)) {
    dataplotting_multi(
        dis_cluster[[i]],
        cut_top    = 4,
        cut_bottom = 2,
        out_dir = 'fig/icd10_cluster_20'
    )
}
```

| ![](fig/icd10_cluster_20/1_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/2_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/3_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/4_normalized_cut_4,2_plot.png)  |
| -------------------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------------- | -------------------------------------------------------- |
| ![](fig/icd10_cluster_20/5_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/6_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/7_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/8_normalized_cut_4,2_plot.png)  |
| ![](fig/icd10_cluster_20/9_normalized_cut_4,2_plot.png)  | ![](fig/icd10_cluster_20/10_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/11_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/12_normalized_cut_4,2_plot.png) |
| ![](fig/icd10_cluster_20/13_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/14_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/15_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/16_normalized_cut_4,2_plot.png) |
| ![](fig/icd10_cluster_20/17_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/18_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/19_normalized_cut_4,2_plot.png) | ![](fig/icd10_cluster_20/20_normalized_cut_4,2_plot.png) |



# III. Age-related diseases: Phecode

## 1. Calculating phecode background incidence

Preparing total incidence of Phecode data as background. Save the result as RDS file for later use.

```R
library(dplyr)
suppressMessages(library(stringr))

# Read data
master6 <- readRDS("Jinhee Code/master_6.rds")
phecodes = read.delim("phecode-ukbb-agg.tsv",stringsAsFactors=F)
phecodes$phecode <- sprintf("%.2f", phecodes$phecode)

# Calculate the total incidence of Phecode
source('src/ard.r')
total_phecode_freq = original_phecode(
    master6,
    code_num      = NULL, 
    code_ann      = phecodes, 
    totalage_freq = NULL,
    subgroup=FALSE)[[1]]
)

# Save as RDS file
saveRDS(total_phecode_freq,'total_phecode_freq.rds')

# Create barplot of total phecode disease incidences
total_phecode_freq_vec = total_phecode_freq$Frequency
names(total_phecode_freq_vec) = total_phecode_freq$`Age at Diagnosis`
png('fig/phecode/total_phecode_freq.png')
barplot(total_phecode_freq_vec,main='Total phecode incidences',xlab='Age at Diagnosis')
dev.off()
```

> dim(master6)
> [1] 4234599       9
> dim(master6_phecode)
> [1] 3388553       9

Phecode covers 3.4M records (80% of total 4.2M records) with age range 27-80 years old.

![](fig/phecode/total_phecode_freq.png)

### Examples: Age-related diseases (ARDs)

Draw incidence rate plot for major age-related diseases and its subgroups.

* 290	Delirium dementia and amnestic and other cognitive disorders
  * 290.1 Dementias
  * **290.11 Alzheimer's disease**
  * 290.12 Dementia with cerebral degenerations
  * 290.16 Vascular dementia
  * 290.2 Delirium due to conditions classified elsewhere
  * 290.3 Other persistent mental disorders due to conditions classified elsewhere
* 428	Congestive heart failure; nonhypertensive
  * 428.1 Congestive heart failure (CHF) NOS
  * 428.2 Heart failure NOS
* No - Chronic obstructive pulmonary disease (COPD)
* 411	Ischemic Heart Disease
  * 411.1 Unstable angina (intermediate coronary syndrome)
  * 411.2 Myocardial infarction
  * 411.3 Angina pectoris
  * 411.4 Coronary atherosclerosis
  * 411.41 Aneurysm and dissection of heart
  * 411.8 Other chronic ischemic heart disease, unspecified
  * 411.9 Other acute and subacute forms of ischemic heart disease
* Stroke, haemorrage or infarction
* 433	Cerebrovascular disease
  * 433.1 Occlusion and stenosis of precerebral arteries
  * 433.11 Occlusion of cerebral arteries, with cerebral infarction
  * 433.12 Cerebral atherosclerosis
  * 433.2 Occlusion of cerebral arteries
  * 433.21 Cerebral artery occlusion, with cerebral infarction
  * 433.3 Cerebral ischemia
  * 433.31 Transient cerebral ischemia
  * 433.32 Moyamoya disease
  * 433.5 Cerebral aneurysm
  * 433.8 Late effects of cerebrovascular disease
* 250.2	Type 2 diabetes
  * 250.21	Type 2 diabetes with ketoacidosis
  * 250.22	Type 2 diabetes with renal manifestations
  * 250.23	Type 2 diabetes with ophthalmic manifestations
  * 250.24	Type 2 diabetes with neurological manifestations
* 250.1	Type 1 diabetes

```R
library(dplyr)

# Read data
master6 <- readRDS("Jinhee Code/master_6.rds")
phecodes = read.delim("phecode-ukbb-agg.tsv",stringsAsFactors=F)
phecodes$phecode <- sprintf("%.2f", phecodes$phecode)
total_phecode_freq = readRDS('total_phecode_freq.rds')

phecode_icds = c(phecodes$icd9_ukb,phecodes$icd10_ukb)
phecode_icd = strsplit(phecode_icds,'\\, ') %>% unlist %>% unique
phecode_icd = str_replace_all(phecode_icd, "[^[:alnum:]]", "")
master6_phecode = subset(master6, Disease_Code %in% phecode_icd)

# Set configurations
source('src/ard.r')
ards_phecodes = c("290","428","411","433","250.1","250.2")
out_dir = 'fig/phecode'

for(i in 1:length(ards_phecodes)) {
    # Prepare data
    inc_rate = original_phecode(
        master6_phecode, 
        code_num=ards_phecodes[i], 
        code_ann=phecodes, 
        freq=200,
        totalage_freq=total_phecode_freq)
    inc_norm = normalized(inc_rate)
    
    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=27,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=27,age_max=80,out_dir)
}
```

> ** Run original_phecode **
>
> Processing 5 iterations:
> 250.10 Type 1 diabetes, N = 4000 (25001,25011,E100,E106,E107,E108,E109)
> 250.11 Type 1 diabetes with ketoacidosis, N = 378 (25011,E101)
> 250.12 Type 1 diabetes with renal manifestations, N = 130 (E102)
> 250.13 Type 1 diabetes with ophthalmic manifestations, N = 711 (E103)
> 250.14 Type 1 diabetes with neurological manifestations, N = 253 (E104)
>
> ** Run original_phecode **
>
> Processing 5 iterations:
> 250.20 Type 2 diabetes, N = 32951 (25000,25010,E110,E116,E117,E118,E119,E135,E136,E137,E138,E139,E149)
> 250.21 Type 2 diabetes with ketoacidosis, N = 236 (2501,E111,E131)
> 250.22 Type 2 diabetes with renal manifestations, N = 394 (E112)
> 250.23 Type 2 diabetes with ophthalmic manifestations, N = 2572 (E103,E113)
> 250.24 Type 2 diabetes with neurological manifestations, N = 1399 (E104,E114,E134,G590)
>
> ...

| Disease                                                      | Incidence rate                            | Normalized incidence rate                   |
| ------------------------------------------------------------ | ----------------------------------------- | ------------------------------------------- |
| 250.2 Type 2 diabetes                                        | ![](fig/phecode/250.2_original_log10.png) | ![](fig/phecode/250.2_normalized_log10.png) |
| 290 Delirium dementia and amnestic and other cognitive disorders | ![](fig/phecode/290_original_log10.png)   | ![](fig/phecode/290_normalized_log10.png)   |
| 411 Ischemic Heart Disease                                   | ![](fig/phecode/411_original_log10.png)   | ![](fig/phecode/411_normalized_log10.png)   |
| 428 Congestive heart failure; nonhypertensive                | ![](fig/phecode/428_original_log10.png)   | ![](fig/phecode/428_normalized_log10.png)   |
| 433 Cerebrovascular disease                                  | ![](fig/phecode/433_original_log10.png)   | ![](fig/phecode/433_normalized_log10.png)   |



### Examples: non-ARDs

```R
source('src/ard.r')
ards_phecodes = c("250.1","335")
out_dir = 'fig/phecode'

for(i in 2:length(ards_phecodes)) {
    # Prepare data
    inc_rate = original_phecode(
        master6_phecode, 
        code_num=ards_phecodes[i], 
        code_ann=phecodes, 
        freq=200,
        totalage_freq=total_phecode_freq)
    inc_norm = normalized(inc_rate)
    
    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=27,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=27,age_max=80,out_dir)
}
```

| Disease                | Incidence rate                            | Normalized incidence rate                   |
| ---------------------- | ----------------------------------------- | ------------------------------------------- |
| 250.1 Type 1 diabetes  | ![](fig/phecode/250.1_original_log10.png) | ![](fig/phecode/250.1_normalized_log10.png) |
| 335 Multiple sclerosis | ![](fig/phecode/335_original_log10.png)   | ![](fig/phecode/335_normalized_log10.png)   |



## 2. Preprocessing for Phecode

Calculating incidence rate by Phecode.

```R
library(dplyr)

# Read data
master6 <- readRDS("Jinhee Code/master_6.rds")
phecodes = read.delim("phecode-ukbb-agg.tsv",stringsAsFactors=F)
phecode_codes = sprintf("%.2f", phecodes$phecode)
phecodes$phecode = phecode_codes
total_phecode_freq = readRDS('total_phecode_freq.rds')

phecode_icds = c(phecodes$icd9_ukb,phecodes$icd10_ukb)
phecode_icd = strsplit(phecode_icds,'\\, ') %>% unlist %>% unique
phecode_icd = str_replace_all(phecode_icd, "[^[:alnum:]]", "") %>% sort
master6_phecode = subset(master6, Disease_Code %in% phecode_icd)

# Run clustering_process for calculating incidence rate (IR)
source('src/ard.r')
clust_result_phecode = clustering_preprocess_phecode(
    master         = master6_phecode,
    code_num       = phecode_codes,
    code_ann       = phecodes,
    freq           = 200,
    totalage_freq  = total_phecode_freq,
	normalized     = TRUE,
    subgroup       = FALSE
)

saveRDS(clust_result_phecode,'clust_result_phecode_ir.rds')
```



Calculating normalized incidence rate by Phecode.

```R
# Run clustering_preprocess for calculating normalized IR
## Minimum incidence criteria: 200
source('src/ard.r')
clust_result_phecode = clustering_preprocess_phecode(
    master         = master6_phecode,
    code_num       = phecode_codes,
    code_ann       = phecodes,
    freq           = 200,
    totalage_freq  = total_phecode_freq,
	normalized     = TRUE,
    subgroup       = FALSE
)

saveRDS(clust_result_phecode,'clust_result_phecode_nir.rds')
```

> ** Run clustering_preprocess_phecode **
>
> Processing 1674 iterations:
>
>   100/1674 175.00 Job process: 7.2 min
>   200/1674 244.10 Job process: 12.8 min
>   300/1674 272.13 Job process: 15.2 min
>   400/1674 290.20 Job process: 18 min
>   500/1674 342.00 Job process: 22.5 min
>   600/1674 371.33 Job process: 27.1 min
>   700/1674 425.12 Job process: 31.8 min
>   800/1674 458.00 Job process: 35.4 min
>   900/1674 525.00 Job process: 39.6 min
>   1000/1674 571.60 Job process: 44.3 min
>   1100/1674 601.12 Job process: 48.9 min
>   1200/1674 637.00 Job process: 54 min
>   1300/1674 703.00 Job process: 58 min
>   1400/1674 733.40 Job process: 1 hr
>   1500/1674 758.10 Job process: 1.2 hr
>   1600/1674 913.00 Job process: 1.3 hr
> Merging data = 54  1049 -> done
>
> Job done: 2021-02-21 10:09:12 for 1.4 hr



## 3. Data pruning: Draw heatmap of Phecode

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Read data
clust_result_phecode = readRDS('clust_result_phecode.rds')

# Prepare sub-sets
## clust_result_phecode:  54 1,040; original data
## clust_result_phecode2: 54   341; peak age over 60
## clust_result_phecode3: 52   290; peak age over 60 & remove bottom 2 rows (< 1,500 incidences)
source('src/ard.r')
clust_result_phecode2 = subsetting_cluster_result(clust_result_phecode, 0,0, 60)
clust_result_phecode3 = subsetting_cluster_result(clust_result_phecode, 0,2, 60)
col_fun = colorRamp2(c(0,1,2,3,4,7,10), c("white","Sky Blue","yellow Green","yellow","red","purple","black"))
```

Draw total 1,040 Phecode diseases heatmap

```R
png('fig/phecode/clust_result_phecode.png', width=25,height=8, units='in', res=150)
Heatmap(clust_result_phecode, cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/phecode/clust_result_phecode.png)

Draw 365 diseases filtered by peak age after 60 years

```R
png('fig/phecode/clust_result_phecode_peak60.png',
    	width=25,height=8, units='in', res=150)
Heatmap(clust_result_phecode2, cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/phecode/clust_result_phecode_peak60.png)

Draw 316 diseases with peak age after 60 years and removed bottom 2 rows (< 1,500 incidences)

```R
# set configurations
i = 3
wh = c(25,8)
col_split = c(8,12,20)
f_name = paste0('fig/phecode/clust_result_phecode_peak60_bt2trimed_',col_split[i],'.png')

# Draw heatmaps
png(f_name, width=25,height=8, units='in', res=150)
Heatmap(clust_result_phecode3, cluster_rows = FALSE, col = col_fun, column_split = col_split[i], column_dend_height=unit(1.5,'in'))
dev.off()

# Kill all graphics objects
graphics.off()
```

![](fig/phecode/clust_result_phecode_peak60_bt2trimed_20.png)



## 4. Extract clustered data & draw cluster plots

Load libraries & data

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

master6 = readRDS("Jinhee Code/master_6.rds")
phecodes = read.delim("phecode-ukbb-agg.tsv",stringsAsFactors=F)
phecodes$phecode = sprintf("%.2f", phecodes$phecode)
total_phecode_freq = readRDS('total_phecode_freq.rds')
```

Split normalized incidence rate data by hclust

```R
k = 20
f_name = paste0('fig/phecode/clust_result_phecode_peak60_bt2trimed_',k,'-re.png')
code_nm = 'Phecode'
split_phecode = data.frame(cutree(hclust(dist(t(clust_result_phecode3))), k = k))
split_phecode[,2] = rownames(split_phecode)
colnames(split_phecode) = c('cluster',code_nm)

png(f_name, width=26,height=8, units='in', res=150)
Heatmap(clust_result_phecode3, column_split=split_phecode[,1], cluster_rows=FALSE, col=col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/phecode/clust_result_phecode_peak60_bt2trimed_20-re.png)

Prepare normalized incidence rate data list by cluster

```R
source('src/ard.r')
dis_cluster_phecode = normalized_by_cat_phecode(
    master         = master6_phecode,
    cat_code       = split_phecode,
    code_ann       = phecodes,
    totalage_freq  = total_phecode_freq,
    subgroup       = FALSE)
```

> ** Run normalized_by_cat **
>
> Processing 20 iterations:
>   1/20 1 = 78 -> dim  2752 4; Job process: 3.1 min
>   2/20 2 = 41 -> dim  1436 4; Job process: 4.2 min
>   3/20 3 = 60 -> dim  2201 4; Job process: 6.6 min
>   4/20 4 = 1 -> dim  86 4; Job process: 6.7 min
>   5/20 5 = 23 -> dim  788 4; Job process: 7.1 min
>   6/20 6 = 45 -> dim  1759 4; Job process: 9.2 min
>   7/20 7 = 5 -> dim  208 4; Job process: 9.3 min
>   8/20 8 = 15 -> dim  542 4; Job process: 9.6 min
>   9/20 9 = 2 -> dim  130 4; Job process: 9.7 min
>   10/20 10 = 2 -> dim  130 4; Job process: 9.8 min
>   11/20 11 = 2 -> dim  122 4; Job process: 9.9 min
>   12/20 12 = 5 -> dim  211 4; Job process: 10 min
>   13/20 13 = 1 -> dim  75 4; Job process: 10 min
>   14/20 14 = 1 -> dim  87 4; Job process: 10.1 min
>   15/20 15 = 1 -> dim  85 4; Job process: 10.1 min
>   16/20 16 = 2 -> dim  108 4; Job process: 10.1 min
>   17/20 17 = 2 -> dim  133 4; Job process: 10.2 min
>   18/20 18 = 2 -> dim  119 4; Job process: 10.3 min
>   19/20 19 = 1 -> dim  98 4; Job process: 10.5 min
>   20/20 20 = 1 -> dim  87 4; Job process: 10.5 min
>
> Process done.

Draw normalized incidence rate as line plot by cluster

```R
source('src/ard.r')
for (i in 1:length(dis_cluster_phecode)) {
    dataplotting_multi(
        dis_cluster_phecode[[i]],
        cut_top    = 0,
        cut_bottom = 2,
        out_dir = 'fig/phecode_cluster_20'
    )
}
```

| ![](fig/phecode_cluster_20/1_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/2_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/3_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/4_normalized_cut_0,2_plot.png)  |
| ---------------------------------------------------------- | ---------------------------------------------------------- | ---------------------------------------------------------- | ---------------------------------------------------------- |
| ![](fig/phecode_cluster_20/5_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/6_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/7_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/8_normalized_cut_0,2_plot.png)  |
| ![](fig/phecode_cluster_20/9_normalized_cut_0,2_plot.png)  | ![](fig/phecode_cluster_20/10_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/11_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/12_normalized_cut_0,2_plot.png) |
| ![](fig/phecode_cluster_20/13_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/14_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/15_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/16_normalized_cut_0,2_plot.png) |
| ![](fig/phecode_cluster_20/17_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/18_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/19_normalized_cut_0,2_plot.png) | ![](fig/phecode_cluster_20/20_normalized_cut_0,2_plot.png) |



# IV. Death (ICD-10)

## 1. Calculating Death background incidence

Preparing total incidence of Death data as background. Save the result as RDS file for later use.

```R
library(dplyr)

# Read data
death = readRDS("Jinhee Code/Deathfile.rds")
colnames(death) = c('eid','Diagnosed_age','Cause_of_death','Disease_Code','meaning')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')

# Calculate the total death count
source('src/ard.r')
total_death_freq = original(
    death, 
    code_num = NULL, 
    code_ann = ICD10_Code_ann, 
    totalage_freq = NULL, 
    subgroup = FALSE)[[1]]

# Save as RDS file
saveRDS(total_death_freq,'total_death_freq.rds')

# Create barplot of total ICD-10 disease incidences
total_death_freq_vec = total_death_freq$Frequency
names(total_death_freq_vec) = total_death_freq$`Age at Diagnosis`
png('fig/death/total_death_freq.png')
barplot(total_death_freq_vec,main='Total death incidences',xlab='Age at Death')
dev.off()
```

Death (ICD10) have > 47K records with age range 40-81 years old.

![](fig/death/total_death_freq.png)

### Examples: Age-related diseases

Draw incidence rate plot for major age-related diseases and its subgroups.

```R
library(dplyr)

# Read data
death = readRDS("Jinhee Code/Deathfile.rds")
colnames(death) = c('eid','Diagnosed_age','Cause_of_death','Disease_Code','meaning')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')
total_death_freq = readRDS('total_death_freq.rds')

# Calculate incidence rates
source('src/ard.r')
dis = c("G30", "E11", "E14", "I21", "I22", "I50", "I64", "J44")
out_dir = 'fig/death'

for(i in 1:length(dis)) {
    # Prepare data
    inc_rate = original(death,dis[i],ICD10_Code_ann,150,total_death_freq,TRUE)
    inc_norm = normalized(inc_rate)

    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=40,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=40,age_max=80,out_dir)
}
```

| Disease                                                      | Incidence rate                        | Normalized incidence rate               |
| ------------------------------------------------------------ | ------------------------------------- | --------------------------------------- |
| G30 Alzheimer's disease (AD)                                 | ![](fig/death/G30_original_log10.png) | ![](fig/death/G30_normalized_log10.png) |
| J44 Other chronic obstructive pulmonary disease (COPD)       | ![](fig/death/J44_original_log10.png) | ![](fig/death/J44_normalized_log10.png) |
| I21 Acute myocardial infarction (Acute MI)                   | ![](fig/death/I21_original_log10.png) | ![](fig/death/I21_normalized_log10.png) |
| I22 Subsequent myocardial infarction (Subsequent MI)         | <100 records                          | <100 records                            |
| I64 Stroke, not specified as haemorrhage or infarction       | ![](fig/death/I64_original_log10.png) | ![](fig/death/I64_normalized_log10.png) |
| E11 Non-insulin-dependent diabetes mellitus (Diabetes Type II) | ![](fig/death/E11_original_log10.png) | ![](fig/death/E11_normalized_log10.png) |
| E14 Unspecified diabetes mellitus (Unspecified Diabetes)     | ![](fig/death/E14_original_log10.png) | ![](fig/death/E14_normalized_log10.png) |

### Example: non-ARDs

```R
# Calculate incidence rates
source('src/ard.r')
dis = c("G35")
out_dir = 'fig/death'

for(i in 1:length(dis)) {
    # Prepare data
    inc_rate = original(death,dis[i],ICD10_Code_ann,0,total_death_freq,TRUE)
    inc_norm = normalized(inc_rate)

    # Draw plots
    dataplotting(inc_rate,logbase=10,age_min=40,age_max=80,out_dir)
    dataplotting(inc_norm,logbase=10,age_min=40,age_max=80,out_dir)
}
```



## 2. Preprocessing for Death

Calculating incidence of each ICD10 code of Death.

```R
library(dplyr)

# Read data
death = readRDS("Jinhee Code/Deathfile.rds")
colnames(death) = c('eid','Diagnosed_age','Cause_of_death','Disease_Code','meaning')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')
total_death_freq = readRDS('total_death_freq.rds')
icd10_death = death$Disease_Code %>% as.character %>% unique %>% sort

# Calculate incidence rate of ICD10
## Minimum incidence criteria: 100
source('src/ard.r')
clust_result_death_ir = clustering_preprocess(
    master         = death,
    code_num       = icd10_death,
    code_ann       = ICD10_Code_ann,
    freq           = 100,
    totalage_freq  = total_death_freq,
    normalized     = FALSE)

saveRDS(clust_result_death_ir,'clust_result_death-ir.rds')
```

> ** Run clustering_preprocess **
>
> Processing 1694 iterations:
>   100/1694 C119 Job process: 0.6 sec
>   200/1694 C509 Job process: 1.2 sec
>   300/1694 C911 Job process: 1.6 sec
>   400/1694 D758 Job process: 2 sec
>   500/1694 F329 Job process: 2.3 sec
>   600/1694 G969 Job process: 2.6 sec
>   700/1694 I600 Job process: 3.1 sec
>   800/1694 J123 Job process: 3.5 sec
>   900/1694 K274 Job process: 3.9 sec
>   1000/1694 K822 Job process: 4.2 sec
>   1100/1694 M628 Job process: 4.5 sec
>   1200/1694 Q676 Job process: 4.9 sec
>   1300/1694 S127 Job process: 5.2 sec
>   1400/1694 T66 Job process: 5.6 sec
>   1500/1694 W158 Job process: 6 sec
>   1600/1694 X73 Job process: 6.3 sec
> List = 96 (freq >100) ; Merging data = 42  97 -> done
>
> Job done: 2021-02-21 21:57:28 for 6.6 sec

Calculating normalized incidence rate of ICD10

```R
# Run clustering_preprocess
## Minimum incidence criteria: 100
source('src/ard.r')
clust_result_death_nir = clustering_preprocess(
    master         = death,
    code_num       = icd10_death,
    code_ann       = ICD10_Code_ann,
    freq           = 100,
    totalage_freq  = total_death_freq,
    normalized     = TRUE)

saveRDS(clust_result_death_nir,'clust_result_death-nir.rds')
```

> ** Run clustering_preprocess **
>
> Processing 1694 iterations:
>   100/1694 C119 Job process: 0.6 sec
>   200/1694 C509 Job process: 1.1 sec
>   300/1694 C911 Job process: 1.7 sec
>   400/1694 D758 Job process: 2.1 sec
>   500/1694 F329 Job process: 2.4 sec
>   600/1694 G969 Job process: 2.8 sec
>   700/1694 I600 Job process: 3.3 sec
>   800/1694 J123 Job process: 3.7 sec
>   900/1694 K274 Job process: 4.2 sec
>   1000/1694 K822 Job process: 4.6 sec
>   1100/1694 M628 Job process: 5 sec
>   1200/1694 Q676 Job process: 5.4 sec
>   1300/1694 S127 Job process: 5.7 sec
>   1400/1694 T66 Job process: 6 sec
>   1500/1694 W158 Job process: 6.3 sec
>   1600/1694 X73 Job process: 6.7 sec
> List = 96 (freq >100) ; Merging data = 42  97 -> done
>
> Job done: 2021-02-21 21:51:04 for 7.1 sec

## 3. Data pruning: Draw heatmap of Death

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Read data
clust_result_death_nir = readRDS('clust_result_death.rds')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')

# Prepare sub-sets
## clust_result_death_ir:   42 96: incidence rate
## clust_result_death_ir3:  40 27: incidence rate

## clust_result_death_nir:  42 96; normalized incidence rate
## clust_result_death_nir2: 42 35; peak age over 60
## clust_result_death_nir3: 40 27; peak age over 60 & remove top 1 and bottom 2 rows (< 200 death)
source('src/ard.r')
clust_result_death_nir2 = subsetting_cluster_result(clust_result_death_nir, 0,0, 60)
clust_result_death_nir3 = subsetting_cluster_result(clust_result_death_nir, 0,2, 60)
col_fun = colorRamp2(c(0,1,2,3,4,7,10), c("white","Sky Blue","yellow Green","yellow","red","purple","black"))

# Extract pruned list from IR data
col_nms = colnames(clust_result_death_ir)
select_cols3 = colnames(clust_result_death_nir3)
clust_result_death_ir3 = clust_result_death_ir[, col_nms[col_nms %in% select_cols3]]
col_fun2 = colorRamp2(c(min(clust_result_death_ir3),max(clust_result_death_ir3)), c("white","royalblue"))

# Prepare disease name annotations
select_cols = colnames(clust_result_death_nir)
dis_names = subset(ICD10_Code_ann, coding %in% select_cols) %>%
	select(coding, meaning)
dis_names = dis_names[match(select_cols, dis_names$coding), ]
dis_ha = dis_names$meaning
#names(dis_ha) = dis_names$coding
ha = HeatmapAnnotation(Name=anno_text(dis_ha))
```

Draw NIR of total 96 Phecode diseases heatmap

```R
png('fig/death/clust_result_death_merge.png', width=20,height=18, units='in', res=150)

H1 = Heatmap(clust_result_death_nir, name="Norm.\nIR", cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
H2 = Heatmap(clust_result_death_ir, name="Incidence\nRate", cluster_rows = FALSE, col = col_fun2, column_dend_height=unit(1.5,'in'), bottom_annotation = ha, show_column_names = FALSE)
draw(H1 %v% H2)
dev.off()
```

![](fig/death/clust_result_death_merge.png)

Draw 35 diseases filtered by peak age after 60 years

```R
png('fig/death/clust_result_death-peak60.png',
    	width=10,height=8, units='in', res=150)
Heatmap(clust_result_death2, name="Norm.\nIR", cluster_rows = FALSE, col = col_fun, column_dend_height=unit(1.5,'in'))
dev.off()
```

![](fig/death/clust_result_death-peak60.png)

Draw 27 diseases with peak age after 60 years and removed bottom 2 rows (< 200 death)

```R
# set configurations
i = 2
wh = c(10,8)
col_split = c(5,8)
f_name = paste0('fig/death/clust_result_death-peak60_bt2trimed_',col_split[i],'.png')

# Draw heatmaps
png(f_name, width=wh[1],height=wh[2], units='in', res=150)
Heatmap(clust_result_death3, name="Norm.\nIR", cluster_rows = FALSE, col = col_fun, column_split = col_split[i], column_dend_height=unit(1.5,'in'))
dev.off()

# Kill all graphics objects
graphics.off()
```

![](fig/death/clust_result_death-peak60_bt2trimed_8.png)



## 4. Extract clustered data & draw cluster plots

Load libraries and data

```R
library(dplyr)
library(ComplexHeatmap)
library(circlize)

death = readRDS("Jinhee Code/Deathfile.rds")
colnames(death) = c('eid','Diagnosed_age','Cause_of_death','Disease_Code','meaning')
ICD10_Code_ann = read.delim('ICD10_DataCoding_41270.tsv')
total_death_freq = readRDS('total_death_freq.rds')

# Prepare disease name annotations
select_cols3 = colnames(clust_result_death_nir3)
dis_names3 = subset(ICD10_Code_ann, coding %in% select_cols3) %>%
	select(coding, meaning)
dis_names3 = dis_names3[match(select_cols3, dis_names3$coding), ]
dis_ha3 = dis_names3$meaning
ha3 = HeatmapAnnotation(Name=anno_text(dis_ha3))
```

Split normalized incidence rate data by hclust

```R
k = 8
wh = c(10,17)
f_name = paste0('fig/death/clust_result_death-peak60_bt2trimed_',k,'-merge.png')
code_nm = 'Death'
split_death = data.frame(cutree(hclust(dist(t(clust_result_death_nir3))), k = k))
split_death[,2] = rownames(split_death)
colnames(split_death) = c('cluster',code_nm)

# Draw heatmap
png(f_name, width=wh[1],height=wh[2], units='in', res=150)
H1 = Heatmap(clust_result_death_nir3, name="Norm.\nIR", column_split=split_death[,1], cluster_rows=FALSE, col=col_fun, column_dend_height=unit(1.5,'in'))
H2 = Heatmap(clust_result_death_ir3, name="Incidence\nRate", cluster_rows = FALSE, col = col_fun2, column_dend_height=unit(1.5,'in'), bottom_annotation = ha3, show_column_names = FALSE)
draw(H1 %v% H2)
dev.off()
```

![](fig/death/clust_result_death-peak60_bt2trimed_8-merge.png)

Prepare normalized incidence rate data list by cluster

```R
# Minimum death count is >100
source('src/ard.r')
death_cluster = normalized_by_cat(
    master         = death,
    cat_code       = split_death,
    code_ann       = ICD10_Code_ann,
    freq           = 100,
    totalage_freq  = total_death_freq,
    subgroup       = FALSE)
```

> ** Run normalized_by_cat **
>
> Processing 8 iterations:
>   1/8 1 = 3 -> dim  106 4; Job process: 0.3 sec
>   2/8 2 = 11 -> dim  286 4; Job process: 0.5 sec
>   3/8 3 = 1 -> dim  72 4; Job process: 0.6 sec
>   4/8 4 = 5 -> dim  185 4; Job process: 0.8 sec
>   5/8 5 = 2 -> dim  85 4; Job process: 0.9 sec
>   6/8 6 = 2 -> dim  102 4; Job process: 1.1 sec
>   7/8 7 = 1 -> dim  0 0; Job process: 1.2 sec
>   8/8 8 = 2 -> dim  70 4; Job process: 1.3 sec
>
> Job done: 2021-02-21 20:30:04 for 1.3 sec

Draw normalized incidence rate as line plot by cluster

```R
source('src/ard.r')
for (i in 1:length(death_cluster)) {
    dataplotting_multi(
        death_cluster[[i]],
        cut_top    = 0,
        cut_bottom = 2,
        out_dir = 'fig/death_cluster_8'
    )
}
```

| ![](fig/death_cluster_8/1_normalized_cut_0,2_plot.png) | ![](fig/death_cluster_8/2_normalized_cut_0,2_plot.png) |
| ------------------------------------------------------ | ------------------------------------------------------ |
| ![](fig/death_cluster_8/3_normalized_cut_0,2_plot.png) | ![](fig/death_cluster_8/4_normalized_cut_0,2_plot.png) |
| ![](fig/death_cluster_8/5_normalized_cut_0,2_plot.png) | ![](fig/death_cluster_8/6_normalized_cut_0,2_plot.png) |
| ![](fig/death_cluster_8/7_normalized_cut_0,2_plot.png) | ![](fig/death_cluster_8/8_normalized_cut_0,2_plot.png) |

