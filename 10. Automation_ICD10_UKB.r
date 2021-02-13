### libraries ###

library(dplyr)
library(tidyr)
library(future.apply)


### functions ###

##Function that creates incidence rate of subset disease, total disease and total age distribution
# code_num: Disease code from ICD10
original <- function(
    master = NULL,
    code_num = NULL,
    ICD10_Code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    subgroup = TRUE
) {
    if(length(code_num)>0) {
        master_sub <- master[grep(code_num, master$Disease_Code), ]
        code_meaning_df2 <- subset(ICD10_Code_ann, ICD10_Code_ann$coding == code_num)
        code_meaning2 = code_meaning_df2$meaning
    } else {
        paste0('\n[NOTE] No input code_num -> Counting total code frequency.\n') %>% cat
        code_num = "total"
        totalage_freq = NULL
        subgroup = FALSE
        master_sub = master
        code_meaning2 = "Total Diseases"
    }
    col_nm = c("Age at Diagnosis", "Frequency", "Incidence Rate", "category", "Disease ID")

    # Count incidence of subgroups
    if(subgroup) {
        subset_name <- master_sub$Disease_Code %>% as.character %>% unique
        subset_name = subset_name[!subset_name==code_num] # exclude query code_num
        
        if (nrow(master_sub) < freq) return(NULL)

        subset_list = lapply(c(1:length(subset_name)), function(i) {
            master_subset <- master_sub[grep(subset_name[i], master_sub$Disease_Code), ]
            age_freq <- as.data.frame(table(floor(master_subset$Diagnosed_age)))
            if(nrow(age_freq)>0) {
                colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
            } else age_freq = NULL

            if (sum(age_freq$"Frequency") > freq) {
                # Add column: Incidence Rate
                age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)

                # Add columns: category, Disease ID
                code_meaning_df = subset(master_subset, Disease_Code==subset_name[i])
                code_meaning = code_meaning_df$meaning %>% as.character %>% unique
                code_meaning = code_meaning[1]
                if(length(code_meaning)==0) code_meaning = subset_name[i]
                age_freq = data.frame(age_freq,
                    category=paste0(code_meaning, " (n=",nrow(master_subset), ")"),
                    `Disease ID`=subset_name[i])
            } else age_freq = NULL
            return(age_freq)
        })
        subset_li = subset_list[!sapply(subset_list,is.null)]
        if(length(subset_li)>0) {
            totaldataframe <- data.table::rbindlist(subset_li)
            colnames(totaldataframe) = col_nm
        } else return(NULL)
    } else totaldataframe = NULL
    
    # Count incidence of code_num (total)
    age_freq <- as.data.frame(table(floor(master_sub$Diagnosed_age)))
    m = nrow(age_freq)
    if(m>0) {
        colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
        age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)
        age_freq = data.frame(age_freq,
            category = paste0(code_meaning2, " (n=", nrow(master_sub), ")"),
            `Disease ID` = code_num
        )
        colnames(age_freq) = col_nm
    } else return(NULL)

    # Add "Disease ID" column to totalage_freq
    if(is.null(totalage_freq)) { k=0 }
    else k = ncol(totalage_freq)
    if(k==3) {
        totalage_freq = data.frame(totalage_freq,`Disease ID` = "total")
        colnames(totalage_freq) = col_nm
    }
    
    total_res <- rbind(age_freq, totaldataframe, totalage_freq)
    total_res$`Age at Diagnosis` <- as.character(total_res$`Age at Diagnosis`) %>% as.numeric
    total_res$category <- as.factor(total_res$category)
  
    return(list(total_res, code_num, "original"))
}


##Function that normalize the incidence rate based on total age incidence rate
normalized <- function(
    master = NULL,
    code_num = NULL,
    ICD10_Code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    subgroup = TRUE
) {
    result <- master[grep(code_num, master$Disease_Code), ]
    if (nrow(result) < freq) return(NULL)
  
    #normalizeddf <- merge(x = original_ICD9(code_num)[[1]], y = totalage_freq, by = "Age at Diagnosis", all = TRUE)
    original_res = original(
        master,
        code_num,
        ICD10_Code_ann,
        totalage_freq = totalage_freq,
        subgroup = subgroup)[[1]]
    if(!is.null(original_res)) {
        normalized_merge <- merge(
            x = original_res,
            y = totalage_freq[,c(1,3)],
            by = "Age at Diagnosis",
            all = TRUE)

        normalized_merge <- na.omit(normalized_merge)
        normalized_merge$`Normalized Incidence Rate` = normalized_merge$`Incidence Rate.x`/normalized_merge$`Incidence Rate.y`
        normalized_res <- normalized_merge %>% select(`Age at Diagnosis`,`category`,`Normalized Incidence Rate`,`Disease ID`)
        #colnames(normalized_df) <- c("Age at Diagnosis", "category", "Incidence Rate")
    } else return(NULL)
  
    return(list(normalized_res, code_num, "normalized"))
}


clustering_preprocess <- function(master, code_num, ICD10_Code_ann, totalage_freq, subgroup = TRUE) {
    source('src/pdtime.r')
    t0=Sys.time()
    
    paste0('\n** Run clustering_preprocess **\n\n') %>% cat
    n = length(code_num); m = 100
    paste0('Processing ',n,' iterations: \n') %>% cat
    norm_result <- future_lapply(c(1:n), function(i) {
        #paste0(i,'/',n,' ',code_num[i],'\n') %>% cat # for debug
        if(i%%m==0) {paste0('  ',i,'/',n,' ',code_num[i],' ',pdtime(t0,2),'\n') %>% cat}

        norm_df <- normalized(master, code_num[i],
            ICD10_Code_ann, 200, totalage_freq, subgroup)[[1]]
        if(!is.null(norm_df)) {
            #tmp_df <- subset(tmp_list[[1]],tmp_list[[1]]$category != "Total Disease(n=3003823)")
            #colnames(temp_df)[3] <- tmp_list[[2]]
            #subset(norm_df,`Disease ID` != "total")
            return(norm_df)
        } else return(NULL)
    })
    norm_result_df = data.table::rbindlist(norm_result) %>% unique

    paste0('Merging data = ') %>% cat
    if(nrow(norm_result_df)>0) {
        clust_result = norm_result_df %>%
            select(-category) %>%
            spread(key=`Disease ID`,value=`Normalized Incidence Rate`) %>%
            select(-total)
        paste0(dim(clust_result),' ') %>% cat
    } else {
        dim(norm_result_df) %>% print
        return(NULL)
    }

    row_nms = clust_result$`Age at Diagnosis`
    clust_result$`Age at Diagnosis` <- NULL
    clust_result <- as.matrix(clust_result)
    rownames(clust_result) = row_nms
    clust_result[is.na(clust_result)] <- 0
    paste0('-> done\n') %>% cat
    
    paste0('\n',pdtime(t0,1),'\n') %>% cat
    return(clust_result)
}


### End ###