### libraries ###

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))


### functions ###

##Function that creates incidence rate of subset disease, total disease and total age distribution
# code_num: Disease code from ICD-10
original <- function(
    master = NULL,
    code_num = NULL,
    code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    subgroup = TRUE,
    group_id = NULL,
    codeset = 'icd10'
) {
  suppressMessages(library(stringr))

  if (codeset == 'icd10' & length(code_num) == 1) {
    code_num = str_replace_all(code_num, "[^[:alnum:]]", "") # remove all special characters from string
    master_sub <- master[grep(code_num, master$Disease_Code),]
    code_meaning_df2 <- subset(code_ann, coding == code_num)
    code_meaning2 = code_meaning_df2$meaning
  } else if (codeset == 'phecode' & length(code_num) > 0) {
    # Extract icd9
    n = length(code_num$icd9)
    if (n > 0) {
      master_sub_icd9_li = lapply(c(1:n), function(i) {
        code_num_icd9 = str_replace_all(code_num$icd9[i], "[^[:alnum:]]", "")
        master_icd9 = subset(master, ICD9.ICD10 == "ICD9")
        master_icd9 %>% filter(str_detect(
                    Disease_Code,
                    paste0('^', code_num_icd9)))
      })
    } else master_sub_icd9_li = NULL
    master_sub_icd9 = data.table::rbindlist(master_sub_icd9_li) %>% unique
    # Extract icd10
    m = length(code_num$icd10)
    if (m > 0) {
      master_sub_icd10_li = lapply(c(1:m), function(j) {
        code_num_icd10 = str_replace_all(code_num$icd10[j], "[^[:alnum:]]", "")
        master_icd10 = subset(master, ICD9.ICD10 == "ICD10")
        #res = master_icd10[grep(code_num_icd10, master$Disease_Code), ]
        master_icd10 %>% filter(str_detect(
                    Disease_Code,
                    paste0('^', code_num_icd10)
                ))
      })
    } else master_sub_icd10_li = NULL
    master_sub_icd10 = data.table::rbindlist(master_sub_icd10_li) %>% unique

    # Row binding icd9 and icd10
    master_sub = rbind(master_sub_icd9, master_sub_icd10) %>% unique
    group_sub = subset(code_ann, phecode == group_id)
    code_meaning2 = paste0(group_id, ' ', group_sub$phenotype)
  } else if (is.null(code_num)) {
    paste0('\n[NOTE] No input code_num -> Counting total code frequency.\n') %>% cat
    code_num = "total"
    totalage_freq = NULL
    subgroup = FALSE
    master_sub = master
    code_meaning2 = "Total Diseases"
  } else {
    paste0('\n[ERROR] Input data is something wrong. Please check.\n') %>% cat
    stop()
  }
  col_nm = c("Age at Diagnosis", "Frequency", "Incidence Rate", "category", "Disease ID")

  # Count incidence of subgroups
  #code_nums = paste0(code_num %>% unlist, collapse=',')
  #paste0('\n',code_meaning2, ', N = ', nrow(master_sub),' (', code_nums, ')') %>% cat
  #table(master_sub$Disease_Code %>% as.character) %>% print

  if (nrow(master_sub) < freq) return(NULL)
  if (subgroup) {
    subset_name <- master_sub$Disease_Code %>%
            as.character %>% unique
    subset_name = subset_name[!subset_name == code_num] # exclude query code_num
    if (length(subset_name) > 0) {
      subset_list = lapply(c(1:length(subset_name)), function(i) {
        master_subset <- master_sub[grep(subset_name[i], master_sub$Disease_Code),]
        age_freq <- as.data.frame(table(floor(master_subset$Diagnosed_age)))
        if (nrow(age_freq) > 0) {
          colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
        } else age_freq = NULL

        if (sum(age_freq$"Frequency") > freq) {
          # Add column: Incidence Rate
          age_freq$`Incidence Rate` <- age_freq$Frequency / sum(age_freq$Frequency)

          # Add columns: category, Disease ID
          code_meaning_df = subset(master_subset, Disease_Code == subset_name[i])
          code_meaning = code_meaning_df$meaning %>% as.character %>% unique
          code_meaning = code_meaning[1]
          if (length(code_meaning) == 0) code_meaning = subset_name[i]
          age_freq = data.frame(age_freq,
                        category = paste0(code_meaning, " (n=", nrow(master_subset), ")"),
                        `Disease ID` = subset_name[i])
        } else age_freq = NULL
        return(age_freq)
      })
      subset_li = subset_list[!sapply(subset_list, is.null)]
      if (length(subset_li) > 0) {
        subgroups_df <- data.table::rbindlist(subset_li)
        colnames(subgroups_df) = col_nm
        subgroups_df$`Age at Diagnosis` <- as.character(subgroups_df$`Age at Diagnosis`) %>% as.numeric
        subgroups_df$category <- as.factor(subgroups_df$category)
      } else return(NULL)
    } else subgroups_df = NULL
  } else subgroups_df = NULL

  # Count incidence of code_num (total)
  age_freq <- as.data.frame(table(floor(master_sub$Diagnosed_age)))
  m = nrow(age_freq)
  if (m > 0) {
    colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
    age_freq$`Incidence Rate` <- age_freq$Frequency / sum(age_freq$Frequency)
    if (is.null(group_id)) {
      dis_name = code_num
    } else dis_name = group_id
    age_freq = data.frame(
        age_freq,
        category = paste0(code_meaning2, " (n=", nrow(master_sub), ")"),
        `Disease ID` = dis_name
    )
    colnames(age_freq) = col_nm
    age_freq$`Age at Diagnosis` <- as.character(age_freq$`Age at Diagnosis`) %>% as.numeric
    age_freq$category <- as.factor(age_freq$category)
  } else return(NULL)

  # Add "Disease ID" column to totalage_freq
  if (is.null(totalage_freq)) { k = 0 }
  else k = ncol(totalage_freq)
  if (k == 3) {
    totalage_freq = data.frame(totalage_freq, `Disease ID` = "total")
    colnames(totalage_freq) = col_nm
  }

  total_res <- rbind(age_freq, subgroups_df, totalage_freq)
  return(list(total_res, code_num, "original"))
}


original_phecode = function(
    master = NULL,
    code_num = NULL,
    code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    subgroup = FALSE
) {
  paste0('\n** Run original_phecode **\n\n') %>% cat
  suppressMessages(library(stringr))
  suppressMessages(library(future.apply))

  # Convert a phecode to icd9/10
  #code_ann_sub = code_ann[grep(code_num, code_ann$phecode), ]
  if (!is.null(code_num)) {
    code_ann_sub = code_ann %>% filter(str_detect(phecode, paste0('^', code_num)))
    code_num_li = lapply(c(1:nrow(code_ann_sub)), function(i) {
      row = code_ann_sub[i,]
      res_icd9 = strsplit(row$icd9_ukb, '\\, ') %>% unlist %>% unique
      res_icd10 = strsplit(row$icd10_ukb, '\\, ') %>% unlist %>% unique
      res_icd9 = res_icd9[!is.na(res_icd9)]
      res_icd10 = res_icd10[!is.na(res_icd10)]
      list(icd9 = res_icd9, icd10 = res_icd10)
    })
    names(code_num_li) = code_ann_sub$phecode
    k = length(code_num_li)
  } else {
    # Extract total phecode-icd codes
    phecode_icd9s = phecodes$icd9_ukb %>% unique
    phecode_icd10s = phecodes$icd10_ukb %>% unique

    phecode_icd9 = strsplit(phecode_icd9s, '\\, ') %>% unlist %>% unique
    phecode_icd10 = strsplit(phecode_icd10s, '\\, ') %>% unlist %>% unique
    phecode_icd9 = phecode_icd9[!is.na(phecode_icd9)]
    phecode_icd10 = phecode_icd10[!is.na(phecode_icd10)]

    # Extract data by ICD codes
    n = length(phecode_icd9)
    master_sub_icd9_li = future_lapply(c(1:n), function(i) {
      code_num_icd9 = str_replace_all(phecode_icd9[i], "[^[:alnum:]]", "")
      master_icd9 = subset(master, ICD9.ICD10 == "ICD9")
      master_icd9 %>% filter(str_detect(Disease_Code, paste0('^', code_num_icd9)))
    })
    master_sub_icd9 = data.table::rbindlist(master_sub_icd9_li) %>% unique

    m = length(phecode_icd10)
    master_sub_icd10_li = future_lapply(c(1:m), function(j) {
      code_num_icd10 = str_replace_all(phecode_icd10[j], "[^[:alnum:]]", "")
      master_icd10 = subset(master, ICD9.ICD10 == "ICD10")
      master_icd10[grep(code_num_icd10, master$Disease_Code),]
    })
    master_sub_icd10 = data.table::rbindlist(master_sub_icd10_li) %>% unique
    master_sub = rbind(master_sub_icd9, master_sub_icd10) %>% unique

    # Configuration to run original function
    master = master_sub # switch input to result
    k = 1;
    code_num = NULL;
    group_id = NULL;
    totalage_freq = NULL
  }
  paste0('Processing ', length(code_num_li), ' iterations:\n') %>% cat

  ori_res_li = lapply(c(1:k), function(i) {
    if (!is.null(code_num)) {
      code_num = code_num_li[[i]]
      group_id = names(code_num_li[i])
    }

    res_df = original(
            master,
            code_num = code_num,
            code_ann = code_ann,
            freq = freq,
            totalage_freq = totalage_freq,
            subgroup = subgroup,
            group_id = group_id,
            codeset = 'phecode'
        )[[1]]
  })
  ori_res_df = data.table::rbindlist(ori_res_li) %>% unique

  return(list(ori_res_df, code_num, "original"))
}


##Function that normalize the incidence rate based on total age incidence rate
normalized <- function(
    original_result = NULL
) {
  #paste0('\n** Run normalized **\n\n') %>% cat

  if (!is.null(original_result)) {
    totalage_freq = original_result[[1]] %>%
            filter(`Disease ID` == "total")

    original_res = original_result[[1]]
    normalized_merge <- merge(
            x = original_res,
            y = totalage_freq[, c(1, 3)],
            by = "Age at Diagnosis",
            all = TRUE)

    normalized_merge <- na.omit(normalized_merge)
    normalized_merge$`Normalized Incidence Rate` = normalized_merge$`Incidence Rate.x` / normalized_merge$`Incidence Rate.y`
    normalized_res <- normalized_merge %>% select(`Age at Diagnosis`, `category`, `Normalized Incidence Rate`, `Disease ID`)
    colnames(normalized_res) <- c("Age at Diagnosis", "category", "Incidence Rate", "Disease ID")
  } else {
    #paste0('[ERROR] NULL input detected at normlized function!\n') %>% cat
    return(NULL)
  }

  return(list(normalized_res, original_result[[2]], "normalized"))
}


normalized_by_cat = function(
    master = NULL,
    cat_code = NULL,
    code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    subgroup = FALSE
) {
  suppressMessages(library(future.apply))
  source('src/pdtime.r')
  t0 = Sys.time()

  paste0('\n** Run normalized_by_cat **\n\n') %>% cat
  categories = cat_code[, 1] %>% unique
  n = length(categories);
  paste0('Processing ', n, ' iterations:\n') %>% cat
  cat_res_li = lapply(c(1:n), function(i) {
    paste0('  ', i, '/', n, ' ', categories[i], ' = ') %>% cat
    code_num = cat_code[which(cat_code[, 1] %in% categories[i]), 2] %>% unique

    m = length(code_num)
    paste0(m, ' -> ') %>% cat
    norm_li = future_lapply(c(1:m), function(j) {
      original_res = original(master, code_num[j],
                code_ann, freq, totalage_freq, subgroup)
      norm_df <- normalized(original_res)[[1]]
    })
    norm_rbind = data.table::rbindlist(norm_li) %>% unique
    paste0(c('dim ', dim(norm_rbind)), collapse = ' ') %>% cat
    paste0('; ', pdtime(t0, 2)) %>% cat

    return(list(norm_rbind, categories[i], "normalized"))
  })
  paste0('\n', pdtime(t0, 1)) %>% cat
  return(cat_res_li)
}


normalized_by_cat_phecode = function(
    master = NULL,
    cat_code = NULL,
    code_ann = NULL,
    totalage_freq = NULL,
    subgroup = FALSE
) {
  suppressMessages(library(future.apply))
  source('src/pdtime.r')
  t0 = Sys.time()

  paste0('\n** Run normalized_by_cat_phecode **\n\n') %>% cat
  categories = cat_code[, 1] %>% unique
  n = length(categories);
  paste0('Processing ', n, ' iterations:\n') %>% cat
  cat_res_li = future_lapply(c(1:n), function(i) {
    paste0('  ', i, '/', n, ' ', categories[i], ' = ') %>% cat
    code_num = cat_code[which(cat_code[, 1] %in% categories[i]), 2] %>% unique

    m = length(code_num)
    paste0(m, '\t-> ') %>% cat
    norm_li = lapply(c(1:m), function(j) {
      original_res = original_phecode(master, code_num[j],
                code_ann, 200, totalage_freq, subgroup)
      norm_df <- normalized(original_res)[[1]]
    })
    norm_rbind = data.table::rbindlist(norm_li) %>% unique
    paste0(c('dim ', dim(norm_rbind)), collapse = ' ') %>% cat
    paste0('; ', pdtime(t0, 2)) %>% cat

    return(list(norm_rbind, categories[i], "normalized"))
  })
  paste0('\n', pdtime(t0, 1)) %>% cat
  return(cat_res_li)
}


clustering_preprocess <- function(
    master = NULL,
    code_num = NULL,
    code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    normalized = TRUE
) {
  source('src/pdtime.r')
  t0 = Sys.time()
  suppressMessages(library(future.apply))

  paste0('\n** Run clustering_preprocess **\n\n') %>% cat
  n = length(code_num);
  m = 100
  paste0('Processing ', n, ' iterations:\n') %>% cat
  norm_result_li <- future_lapply(c(1:n), function(i) {
    if (i %% m == 0) { paste0('  ', i, '/', n, ' ', code_num[i], ' ', pdtime(t0, 2)) %>% cat }

    original_result = original(master, code_num[i],
            code_ann, freq, totalage_freq, subgroup = FALSE)
    if (normalized) {
      normalized(original_result)[[1]]
    } else return(original_result[[1]])
  })
  norm_result_li[sapply(norm_result_li, is.null)] <- NULL # Remove NULL elements
  paste0('List = ', length(norm_result_li), ' (freq >', freq, ') ') %>% cat

  norm_result_df = data.table::rbindlist(norm_result_li) %>% unique
  #colnames(norm_result_df) = c('Age at Diagnosis','category','Normalized Incidence Rate','Disease ID')

  paste0('; Merging data = ') %>% cat
  if (nrow(norm_result_df) > 0) {
    clust_result = norm_result_df %>%
            select(-category, - Frequency) %>%
            spread(key = `Disease ID`, value = `Incidence Rate`) %>%
            select(-total)
    paste0(dim(clust_result), ' ') %>% cat
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

  paste0('\n', pdtime(t0, 1)) %>% cat
  return(clust_result)
}


clustering_preprocess_phecode = function(
    master = NULL,
    code_num = NULL,
    code_ann = NULL,
    freq = 200,
    totalage_freq = NULL,
    normalized = TRUE
) {
  suppressMessages(library(future.apply))

  source('src/pdtime.r')
  t0 = Sys.time()

  paste0('\n** Run clustering_preprocess_phecode **\n\n') %>% cat
  n = length(code_num);
  m = 100
  paste0('Processing ', n, ' iterations:\n') %>% cat

  norm_result_li = future_lapply(c(1:n), function(i) {
    if (i %% m == 0) { paste0('  ', i, '/', n, ' ', code_num[i], ' ', pdtime(t0, 2)) %>% cat }
    ori_result = original_phecode(master, code_num[i],
                code_ann, freq, totalage_freq, subgroup = FALSE)
    if (normalized) {
      normalized(ori_result)[[1]]
    } else return(ori_result)

  })
  norm_result_li[sapply(norm_result_li, is.null)] <- NULL # Remove NULL elements
  norm_result_df = data.table::rbindlist(norm_result_li) %>% unique

  paste0('Merging data = ') %>% cat
  if (nrow(norm_result_df) > 0) {
    clust_result = norm_result_df %>%
            select(-category) %>%
            spread(key = `Disease ID`, value = `Incidence Rate`) %>%
            select(-total)
    paste0(dim(clust_result), ' ') %>% cat
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

  paste0('\n', pdtime(t0, 1)) %>% cat
  return(clust_result)
}


dataplotting <- function(
    list_exp = NULL,
    logbase = 10,
    age_min = 27,
    age_max = 80,
    out_dir = NULL
) {
  suppressMessages(library(ggplot2))

  if (is.null(list_exp)) return(NULL)

  dis_total = list_exp[[1]] %>% filter(`Disease ID` == "total")
  dis = list_exp[[1]] %>% filter(`Disease ID` != "total")

  if (list_exp[[3]] == "normalized") {
    ylab_nm = "Norm. Incidence Rate"
  } else ylab_nm = list_exp[[3]]

  if (!is.null(logbase)) {
    dis_total$`Incidence Rate` = log(dis_total$`Incidence Rate`)
    dis$`Incidence Rate` = log(dis$`Incidence Rate`)
    ylab_nm = paste0("log", logbase, "(", ylab_nm, ")")
  }

  p = ggplot(data = dis_total, aes(
            x = `Age at Diagnosis`,
            y = `Incidence Rate`
        ), las = 3, size = 1) +
        geom_line(size = 1, color = "black", alpha = .3) +
        geom_smooth(color = "black", linetype = "dashed", alpha = .15) +
        scale_x_continuous(limits = c(age_min, age_max), n.breaks = 10) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(y = ylab_nm, x = "Age at Diagnosis") +
        theme_bw() +
        theme(legend.position = "right",
            legend.direction = "vertical",
            legend.background = element_rect(color = "steelblue", linetype = "solid"),
            legend.text = element_text(size = 8))

  p = p + geom_line(data = dis, aes(
            x = `Age at Diagnosis`,
            y = `Incidence Rate`,
            color = category,
            group = category
        ), alpha = .3, size = 1) +
        geom_smooth(data = dis, aes(
            x = `Age at Diagnosis`,
            y = `Incidence Rate`,
            color = category,
            fill = category,
            group = category
        ), alpha = .15)

  f_name = paste0(out_dir, '/', list_exp[[2]], "_",
        list_exp[[3]], "_log", logbase, ".png")
  ggsave(f_name, p, width = 8, height = 3, units = 'in')
  paste0('Draw plot: ', f_name, '\n') %>% cat
}


dataplotting_multi = function(
    list_exp = NULL,
    cut_top = NULL,
    cut_bottom = NULL,
    out_dir = 'fig'
) {

  if (is.null(list_exp)) return(NULL)
  if (!is.null(cut_top) | is.null(cut_bottom)) {
    if (is.null(cut_top)) cut_top = 0
    if (is.null(cut_bottom)) cut_bottom = 0
    age = list_exp[[1]]$`Age at Diagnosis` %>% unique %>% sort
    k = 1 + cut_top
    l = length(age) - cut_bottom
    age_sub = age[c(k:l)]
    tb = subset(list_exp[[1]], `Age at Diagnosis` %in% age_sub)
  } else tb = list_exp[[1]]

  age_min = min(tb$`Age at Diagnosis`)
  age_max = max(tb$`Age at Diagnosis`)
  tb = tb %>% filter(`Disease ID` != 'total') # remove total data
  dis_num = tb$`Disease ID` %>% unique %>% length
  plotdata <- ggplot(data = tb,
    aes(x = `Age at Diagnosis`,
        y = `Incidence Rate`,
        group = category),
    size = 1) +
    geom_line(aes(group = category), color = "gray40", alpha = .5) +
    scale_x_continuous(limits = c(age_min, age_max), n.breaks = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    stat_summary(aes(color = "mean", shape = "mean", group = 1),
        fun = mean, color = "red", geom = "line", size = 2) +
    labs(y = "Incidence Rate", x = "Age at Diagnosis") +
    ggtitle(paste0('Disease = ', dis_num)) +
    theme_bw()

  print(list_exp[[2]])

  f_name = paste0(out_dir, '/', list_exp[[2]], "_", list_exp[[3]],
        "_cut_", cut_top, ',', cut_bottom, "_plot.png")
  png(f_name, width = 4, height = 3, unit = 'in', res = 300)

  plotdata = plotdata +
        theme(legend.position = "none")
  print(plotdata)
  dev.off()
}


subsetting_cluster_result <- function(
    clust_result = NULL,
    over = 0, # cut row from top to remove outliers
    less = 0, # cut row from bottom to remove outliers
    age_cut = NULL # find peak age over this criteria
) {
  clust_result <- clust_result[over:(nrow(clust_result) - less),]
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

  col_fun = colorRamp2(c(0, 1, 2, 3, 4), c("white", "Sky Blue", "yellow Green", "yellow", "red"))

  png(f_name, width = wh[1], height = wh[2], units = 'in', res = 150)
  Heatmap(hm_mat, cluster_rows = FALSE, col = col_fun, column_split = column_split[i])
  paste0('Draw plot: ', f_name, '\n') %>% cat
  dev.off()
  #graphics.off()
}


subset_master = function( # Need to fix function
    master = NULL,
    code_num = NULL,
    freq = 200
) {
  n = length(code_num)
  paste0('\n',n,' codes were input. Extracting ') %>% cat

  master_sub_li = future_lapply(c(1:n),function(i) {
    master %>% filter(str_detect(Disease_Code, paste0('^', code_num[i]) ))
  })
  master_sub = data.table::rbindlist(master_sub_li) %>% unique
  dim(master_sub) %>% print

  if(nrow(master_sub)<freq) return(NULL)
  return(master_sub)
}


original_person = function(
    master = NULL,
    code_num = NULL,
    code_ann = NULL,
    group_nm = 'No input code group name',
    freq = 200
) {
  suppressMessages(library(stringr))
  suppressMessages(library(future.apply))

  if (!is.null(code_num)) {
    code_num = str_replace_all(code_num, "[^[:alnum:]]", "") # Remove special characters
    #master_sub = master[grep(code_num, master$Disease_Code), ]
    master_sub = subset_master(master,code_num,freq=200)
    if(length(code_num)==1) {
      code_ann_meaning = subset(code_ann, coding == code_num)$meaning
    } else code_ann_meaning = group_nm
  } else {
    # Run total diseases
    master_sub = master
    code_ann_meaning = "total diseases"
  }

  # Select patients having multiple records
  patients_freq = table(master_sub$eid) %>% as.data.frame
  patients_freq = subset(patients_freq, Freq > 1)

  # Remove earliest record per person
  patients = patients_freq$Var1 %>% unique
  n = length(patients)
  master_p_1_li = future_lapply(c(1:n), function(i) {
    patient = patients[i]
    master_p = subset(master_sub, eid == patient)
    master_p = master_p[order(master_p$Diagnosed_age),]
    master_p[1,]
  })
  master_p_1 = data.table::rbindlist(master_p_1_li)
  `%notin%` = Negate(`%in%`)
  master_p_2 = subset(master_sub, eid %notin% patients)
  master_p = rbind(master_p_1, master_p_2)

  # Check freq
  if (nrow(master_p) < freq) {
    paste0('[NOTICE] ', code_num, ' has ', nrow(master_p), ' record (freq <', freq, ').\n') %>% cat
    return(NULL)
  }

  # Count numbers by age
  age_freq = floor(master_p$Diagnosed_age) %>% table %>% as.data.frame
  m = nrow(age_freq)
  if (m > 0) {
    colnames(age_freq) = c('Age', 'Freq')
  } else return(NULL)
  age_freq$Disease_onset_prop = age_freq$Freq / sum(age_freq$Freq)
  age_freq$category = paste0(code_ann_meaning, ' (n=', nrow(master_p), ')')
  if(length(code_num)==1) {
    age_freq$`Disease ID` = code_num
  } else age_freq$`Disease ID` = group_nm

  return(list(age_freq, code_num, "Disease onset prop"))
}


plot_person <- function(
    list_exp = NULL,
    age_min = 27,
    age_max = 80,
    out_dir = NULL
) {
  suppressMessages(library(ggplot2))

  if (is.null(list_exp)) return(NULL)

  #dis_total = list_exp[[1]] %>% filter(`Disease ID` == "total")
  dis = list_exp[[1]]
  dis$Age = as.character(dis$Age) %>% as.numeric
  ylab_nm = list_exp[[3]]

  p = ggplot(data = dis, aes(
        x = Age,
        y = Disease_onset_prop
    ), las = 3, size = 1) +
    geom_line(aes(
        color = category,
        group = category
    ), alpha = .3, size = 1) +
    geom_smooth(aes(
        color = category,
        fill = category,
        group = category
    ), alpha = .15) +
    scale_x_continuous(limits = c(age_min, age_max), n.breaks = 10) +
    geom_vline(xintercept = 40, linetype = "dotted", color = "grey") +
    geom_vline(xintercept = 60, linetype = "dotted", color = "black") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(y = ylab_nm, x = "Age-of-onset") +
    theme_bw() +
    theme(legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_rect(color = "steelblue", linetype = "solid"),
        legend.text = element_text(size = 8))

  f_name = paste0(out_dir, '/', list_exp[[2]], "_",
        list_exp[[3]], ".png")
  ggsave(f_name, p, width = 8, height = 3, units = 'in')
  paste0('Draw plot: ', f_name, '\n') %>% cat
}


plot_gbd = function(
    gbd_data = NULL,
    ylab_nm = NULL,
    out_dir = NULL
) {
  suppressMessages(library(ggplot2))
  if (is.null(gbd_data)) return(NULL)
  age_order = c("1 to 4", "5 to 9", "10 to 14", "15 to 19", "20 to 24", "25 to 29", "30 to 34", "35 to 39", "40 to 44", "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74", "75 to 79", "85 to 89", "90 to 94")
  gbd_data$age_name = factor(gbd_data$age_name, levels = age_order)

  p = ggplot(data = gbd_data, aes(
      x = age_name, y = val
    ), las = 3, size = 1) +
    geom_line(aes(
        color = cause_name,
        group = cause_name
    ), alpha = .3, size = 1) +
    geom_smooth(aes(
        color = cause_name,
        fill = cause_name,
        group = cause_name
    ), alpha = .15) +
    labs(y = ylab_nm, x = "Age range") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, vjust = .5, hjust = .5),
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_rect(color = "steelblue", linetype = "solid"),
        legend.text = element_text(size = 8))

  f_name = paste0(out_dir, '/gbd_', ylab_nm, '.png')
  ggsave(f_name, p, width = 8, height = 3.5, units = 'in')
  paste0('Draw plot: ', f_name, '\n') %>% cat
}


### End ###
