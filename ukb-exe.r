#!/usr/bin/env Rscript
#
# 2019-12-27 First version in UKB-menopause analysis.
# usage examples:
# 1. Longevity
## Rscript ukb-exe.r --ukb_prune -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_lon data_UKB/ukb_prune_lon.rds -anns Longevity -expt 1845 2946 40007 > data_UKB/ukb_prune_lon.rds_log.txt
# 2. Reproductive aging
## Rscript ukb-exe.r --ukb_prune -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_rep data_UKB/ukb_prune_rep.rds -anns Reproductive_Aging > data_UKB/ukb_prune_rep.rds_log.txt
# 3. Confounders 
## Rscript ukb-exe.r --ukb_prune -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_con data_UKB/ukb_prune_con.rds -anns Confounders -expt 1845 2946 40007 > data_UKB/ukb_prune_con.rds_log.txt
# 4. Additional
## Rscript ukb-exe.r --ukb_prune -v -b db_UKB/ukb38640.csv.rds db_UKB/ukb_fid_info.csv.rds db_UKB/ukb_fid_code.csv.rds -t fig_UKB/prune_con2 data_UKB/ukb_prune_con2.rds --fids 1767 > data_UKB/ukb_prune_con2.rds_log.txt

# Help Messages --------------------
help_message = c('
Version: 2020-01-06

Usage: Rscript ukb-exe.r [functions: --ukb_prune ...] [options: --verbose ...]
           <-b base file(s)> <-t target dir/file(s)>
','
Functions:','
    --ukb_prune | -pr     This is a first function for automated pruning for UKB data.','
        -b 1.rds 2.rds 3.rds
                          Three base RDS files are needed: 
                          * 1.rds: A RDS file for UKB phenotypes. This is a
                              directly downloaded file from the UK biobank.
                          * 2.rds: A RDS file for UKB field ID description.
                          * 3.rds: A RDS file for UKB field ID code information.
                          Please keep the file order!
        -t target_dir target_file.rds
                          Two variables are needed for target directory/RDS file:
                          * target_dir: A target directory address for save figures
                          * target_file.rds: A target file for result phenotype table.
                          Please keep the dir/file order!
        --fids
                          Fid numbers, a subset of field IDs to extract reliable
                          subjects.
        --anns
                          Either Reproductive_Aging/Longevity/Confounders,
                          an annotation for a subset of field IDs in UKB field
                          description file (e.g., db_UKB/ukb_id_info.csv.rds),
                          including Reproductive_Aging, Longevity, and Confounders
        --expt
                          Fid numbers, a subset of field IDs as exception for filtering.
','
    --ukb_excld | -ex     This is a second function for select eid for further analysis,
                          such as correlation table and linear regression model.','
        -b 1.rds
                          One base RDS file is needed:
                          * 1.rds: A RDS file for UKB phenotypes. This is a result file
                              from the ukb_prune function.
        -t target_dir target_file.rds
                          Two variables are needed for target directory/RDS file.
        --fids fid1 fid2 fid3 ...
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
','
Global options:
    -b                    Base RDS files are mendatory.
    -t                    
    --verbose    | -v     Rich description for debugging. Default is FALSE.
')

# Command line arguments --------------------
## See the Help Messages section for more information.
source('src/command.r')
args = command(commandArgs(trailingOnly=T), help_message)

# Libraries here --------------------
suppressMessages(library(dplyr))

# Invoke function and running --------------------
## Set general error handler
errors_n = 0 # error count
if(!'b' %in% names(args)) {
    paste0('\n[Error] There is no base files.') %>% cat
    errors_n = errors_n+1 }
if(!'t' %in% names(args)) {
    paste0('\n[Error] There is no target dir/files.' %>% cat)
    errors_n = errors_n+1 }
## Error STOP1
if(errors_n>0) { paste0('\n** ',errors_n,' errors are found, QUIT.') %>% cat; quit() }

## Set general option to run
source('src/pdtime.r'); t0=Sys.time()
### -v|--verbose
if('v' %in% names(args)) { verbose = TRUE
} else if('verbose' %in% names(args)) { verbose = TRUE
} else verbose = FALSE

if(names(args)[1] %in% c('ukb_prune','pr')) {
    ## Function ukb_filter|fl: src/ukb_filter.r
    paste0('\n** Preparing to run function: ',names(args)[1],'\n') %>% cat
    if(length(args$b)==3) {
        ukb_field = readRDS(args$b[1])
        paste0('  Read ',args$b[1],'; col= ',dim(ukb_field)[1],
               ', row= ',dim(ukb_field)[2],'\n') %>% cat
        fid_info  = readRDS(args$b[2])
        paste0('  Read ',args$b[2],'; col= ',dim(fid_info)[1],
               ', row= ',dim(fid_info)[2],'\n') %>% cat
        fid_code  = readRDS(args$b[3])
        paste0('  Read ',args$b[3],'; col= ',dim(fid_code)[1],
               ', row= ',dim(fid_code)[2],'\n') %>% cat
    } else { # Error message
        paste0('\n[Error] Input file number is not three.',
            ' Please see the --help|-h description.') %>% cat
        errors_n = errors_n+1 }
    
    if(length(args$t)==2) {
        dir    = args$t[1]
        f_name = args$t[2]
    } else { # Error message
        paste0('\n[Error] Target dir/file is not two.',
            ' Please see the --help|-h description.') %>% cat
        errors_n = errors_n+1 }
    
    ## Error STOP2
    if(errors_n>0) {
        paste0('\n** ',errors_n,' errors are found, QUIT.') %>% cat; quit() }
    
    if('fids' %in% names(args)) {
        fids = args$fids # Running for specific fids
        paste0('- Subset of fids = ',length(fids),'\n') %>% cat
    } else if('anns' %in% names(args)) {
        fids = fid_info %>% filter(Anno1 %in% args$anns) %>% .$FieldID
        paste0('- Subset of fids = ',length(fids),'\n') %>% cat
    } else {
        fids = fid_info$FieldID
        paste0('- No subset of fids are submitted.\n',
               '  Total ',length(fids),' field IDs in UKB field ID description file\n',
               '  (e.g., ',args$b[2],')\n',
               '  will undergo to filter process.\n') %>% cat
    }
        
    if('expt' %in% names(args)) {
        except_fid = args$expt
        paste0('- ',length(except_fid),
               ' field IDs as exception for filtering:') %>% cat
        paste0(args$expt,collapse=' ') %>% cat
    } else except_fid = NULL
    
    paste0('\n\n** Now running ',names(args)[1],': src/ukb_filter.r\n') %>% cat
    source('src/ukb_filter.r')
    n = length(fids)
    ukb_li = lapply(c(1:n),function(i) {
        t1=Sys.time()
        paste0('\n(',i,'/',n,') fid',fids[i],'\n') %>% cat
        out = ukb_prune(
            dir        = dir,                   # directory for figrues
            ukb_pheno  = ukb_field,             # download file from ukb
            fid        = as.character(fids[i]), # interesting field id
            except_fid = except_fid,            # exceptional fids. No filtering
            fid_info   = fid_info,              # UKB field id information table
            fid_code   = fid_code,              # UKB field code information table
            verbose    = verbose )              # Debugging mode
        if(is.null(out)) { 'Null data\n' %>% cat
        } else summary(out)
        paste0('Sub-',pdtime(t1,2)) %>% cat
        paste0('Total ',pdtime(t0,2),'\n') %>% cat
        return(out)
    })
    merge_li = function(li1,li2) {
        return(merge(li1,li2,by='eid',all=T)) }
    ukb_df = Reduce(merge_li,ukb_li)
    paste0('** Merge process done. col= ',dim(ukb_df)[1],
           ', row= ',dim(ukb_df)[2],'\n') %>% cat
    saveRDS(ukb_df,f_name)
    paste0('** Save RDS: ',f_name,'\n') %>% cat
    pdtime(t0,1) %>% cat

} else if(names(args)[1] %in% c('ukb_exclude','ex')) {
    ## Function ukb_selec|fl: src/ukb_filter.r
    paste0('\n** Preparing to run function: ',names(args)[1],'\n') %>% cat
    
    
} else {
    paste0('\n[Error] There is no such a function, ',names(args)[1],'.',
        ' Please see the --help|-h description.') %>% cat
    errors_n = errors_n+1
    if(errors_n>0) paste0('\n** ',errors_n,' errors are found, QUIT.\n') %>% cat
    quit()
}