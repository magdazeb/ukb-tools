### libraries ###

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


### functions ###

##Function that creates incidence rate of subset disease, total disease and total age distribution
# code_num: Disease code from ICD-10
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
            subgroups_df <- data.table::rbindlist(subset_li)
            colnames(subgroups_df) = col_nm
            subgroups_df$`Age at Diagnosis` <- as.character(subgroups_df$`Age at Diagnosis`) %>% as.numeric
            subgroups_df$category <- as.factor(subgroups_df$category)
        } else return(NULL)
    } else subgroups_df = NULL
    
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
        age_freq$`Age at Diagnosis` <- as.character(age_freq$`Age at Diagnosis`) %>% as.numeric
        age_freq$category <- as.factor(age_freq$category)
    } else return(NULL)

    # Add "Disease ID" column to totalage_freq
    if(is.null(totalage_freq)) { k=0 }
    else k = ncol(totalage_freq)
    if(k==3) {
        totalage_freq = data.frame(totalage_freq,`Disease ID` = "total")
        colnames(totalage_freq) = col_nm
    }

    total_res <- rbind(age_freq, subgroups_df, totalage_freq)
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
        colnames(normalized_res) <- c("Age at Diagnosis", "category", "Incidence Rate", "Disease ID")
    } else return(NULL)
  
    return(list(normalized_res, code_num, "normalized"))
}


multi_normalized = function(
    master = NULL,
    code_num = NULL,
    ICD10_Code_ann = NULL, 
    totalage_freq = NULL, 
    subgroup = TRUE
) {
    library(future.apply)

    #paste0('\n** Run multi_normalized **\n\n') %>% cat
    n = length(code_num); m = 100
    paste0('Processing ',n,' iterations:\n') %>% cat
    norm_result <- future_lapply(c(1:n), function(i) {
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
    
    return(norm_result)
}


normalized_by_cat = function(
    master = NULL,
    cat_code = NULL,
    ICD10_Code_ann = NULL, 
    totalage_freq = NULL, 
    subgroup = FALSE
) {
    library(future.apply)
    source('src/pdtime.r')
    t0=Sys.time()

    paste0('\n** Run normalized_by_cat **\n\n') %>% cat
    categories = cat_code[,1] %>% unique
    n = length(categories);
    paste0('Processing ',n,' iterations:\n') %>% cat
    cat_res_li = lapply(c(1:n),function(i) {
        paste0('  ',i,'/',n,' ',categories[i],' = ') %>% cat
        code_num = cat_code[which(cat_code[,1] %in% categories[i]),2] %>% unique

        m = length(code_num)
        paste0(m,' -> ') %>% cat
        norm_li = lapply(c(1:m), function(j) {
            norm_df <- normalized(master, code_num[j],
                ICD10_Code_ann, 200, totalage_freq, subgroup)[[1]]
        })
        norm_rbind = data.table::rbindlist(norm_li) %>% unique
        paste0(dim(norm_rbind),collapse=' ') %>% cat
        paste0('; ',pdtime(t0,2)) %>% cat

        return(list(norm_rbind, categories[i], "normalized"))
    })
}


clustering_preprocess <- function(
    master = NULL, 
    code_num = NULL, 
    ICD10_Code_ann = NULL, 
    totalage_freq = NULL, 
    subgroup = TRUE
) {
    source('src/pdtime.r')
    t0=Sys.time()
    
    paste0('\n** Run clustering_preprocess **\n\n') %>% cat
    norm_result <- multi_normalized(
        master, code_num, ICD10_Code_ann, totalage_freq, subgroup)
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
    
    paste0('\n',pdtime(t0,1)) %>% cat
    return(clust_result)
}


dataplotting <- function(
    list_exp = NULL, 
    logbase = 10, 
    age_min = 27, 
    age_max = 80,
    f_dir = NULL
) {
    suppressMessages(library(ggplot2))
    if (is.null(list_exp)) return(NULL)
    
    p <- ggplot(data = list_exp[[1]],
        aes(x = `Age at Diagnosis`,
            y = `Incidence Rate`, 
            color = category,
            group = category
        ), las = 3, size = 5)+
        geom_smooth(alpha = .15, aes(fill = category))+
        scale_x_continuous(limits = c(age_min, age_max), n.breaks = 10)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        labs(y = "Incidence Rate", x = "Age at Diagnosis")+
        theme_bw()+
        theme(legend.position = "right", 
            legend.direction = "vertical",
            legend.background = element_rect(color = "steelblue", linetype = "solid"),
            legend.text = element_text(size = 8))
  
    if (!is.null(logbase)) {
        p <- p+aes(y = log(`Incidence Rate`, logbase))+
            labs(y = paste0("log", logbase, "(Incidence Rate)"))
    }
    print(list_exp[[2]])
    
    f_name = paste0(f_dir,'/',list_exp[[2]], "_", list_exp[[3]], "_", logbase, ".png")
    ggsave(f_name,p,width=10,height=6,units='in')
}


dataplotting_multi =  function(
    list_exp = NULL,
    cut_top = NULL,
    cut_bottom = NULL,
    out = 'fig'
) {
  
    if (is.null(list_exp)) return(NULL)
    if(!is.null(cut_top) | is.null(cut_bottom)) {
        if(is.null(cut_top)) cut_top = 0
        if(is.null(cut_bottom)) cut_bottom = 0
        age = list_exp[[1]]$`Age at Diagnosis` %>% unique %>% sort
        k = 1+cut_top
        l = length(age)-cut_bottom
        age_sub = age[c(k:l)]
        tb = subset(list_exp[[1]],`Age at Diagnosis` %in% age_sub)
    } else tb = list_exp[[1]]

    dis_num = 
    plotdata <- ggplot(data = tb, 
        aes(x = `Age at Diagnosis`,
            y = `Incidence Rate`,
            group = category),
        las=1, size=3)+
        geom_line(color="gray40", alpha = .5)+
        scale_x_continuous(n.breaks = 10)+ #limits = c(age_min, age_max)
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        geom_hline(yintercept=1, linetype="dashed",color="black")+
        #geom_smooth(method="loess")+
        labs(y = "Incidence Rate", x = "Age at Diagnosis")+
        theme_bw()
  
    print(list_exp[[2]])
  
    f_name = paste0(out,'/',list_exp[[2]],"_",list_exp[[3]],
        "_cut_",cut_top,',',cut_bottom,"_plot.png")
    png(f_name, width=5,height=4, unit='in', res = 300)
  
    plotdata = plotdata + theme(legend.position = "none")
    print(plotdata)
    dev.off() 
  
}


subsetting_cluster_result <- function(
    clust_result = NULL, 
    over = 0, # cut row from top to remove outliers
    less = 0, # cut row from bottom to remove outliers
    age_cut = NULL # find peak age over this criteria
) {
    clust_result <- clust_result[over:(nrow(clust_result)-less),]
    if (!is.null(age_cut)) {
        max_age <- apply(clust_result, 2, which.max)
        clust_result <- 
            clust_result[, as.numeric(rownames(clust_result)[max_age]) >= age_cut]
    }
    return(clust_result)
}


draw_hm = function(
    hm_mat = NULL,
    f_name = NULL, 
    column_split = NULL
) {
    suppressMessages(library(ComplexHeatmap))
    suppressMessages(library(circlize))

    col_fun = colorRamp2(c(0,1,2,3,4), c("white","Sky Blue","yellow Green","yellow","red"))

    png(f_name, width=wh[1], height=wh[2], units='in', res=150)
    Heatmap(hm_mat, cluster_rows = FALSE, col = col_fun, column_split = column_split[i])
    paste0('Draw plot: ',f_name,'\n') %>% cat
    dev.off()
    #graphics.off()
}


### End ###