# UK Biobank tools

How to use:

* Run command `Rscript ukb-exe.r -h`

* ```
  Version: 2020-01-06
  
  Usage: Rscript ukb-exe.r [functions: --ukb_prune ...] [options: --verbose ...]
             <-b base file(s)> <-t target dir/file(s)>
  
  Functions:
      --ukb_prune | -pr     This is a first function for automated pruning for UKB data.
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
  
      --ukb_excld | -ex     This is a second function for select eid for further analysis,
                            such as correlation table and linear regression model.
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
  
  Global options:
      -b                    Base RDS files are mendatory.
      -t
      --verbose    | -v     Rich description for debugging. Default is FALSE.
  ```

