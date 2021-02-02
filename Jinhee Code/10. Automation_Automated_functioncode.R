library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)

## Create a function that plots data as we want
dataplotting <- function(list_exp, logbase = 10, age_min = 27, age_max = 80) {
  
  if (is.null(list_exp)) {
    return(NULL)
  }
  
  plotdata <- ggplot(data = list_exp[[1]], aes(x = `Age at Diagnosis`,
                                               y = `Incidence Rate`, 
                                               color = category
                                               #group = category
  ), las = 3)+
    geom_line(alpha = .15
              #, aes(fill = category)
    )+
    scale_x_continuous(limits = c(age_min, age_max), n.breaks = 10)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(y = "Incidence Rate", x = "Age at Diagnosis")+
    theme_bw()
  #theme(legend.position = "right", 
  #legend.direction = "vertical",
  #legend.background = element_rect(color = "steelblue", linetype = "solid"),
  #legend.text = element_text(size = 8))
  if (!is.null(logbase)) {
    plotdata <- plotdata+
      aes(y = log(`Incidence Rate`, logbase))+
      labs(y = paste0("log", logbase, "(Incidence Rate)"))
  }
  print(list_exp[[2]])
  
  png(paste0(list_exp[[2]], "_", list_exp[[3]], "_", logbase, "_plot.png"),
      res = 80,
      width = 1000, height = 600)
  
  plotdata = plotdata + theme(legend.position = "none")
  print(plotdata)
  dev.off() 
  
}

##Function that creates incidence rate of subset disease, total disease and total age distribution
# code_num: Disease code from ICD10
original <- function(code_num, freq = 200) {
  
  result <- master6[grep(code_num, master6$Disease_Code), ]
  subset_name <- result %>% distinct(Disease_Code)
  
  if (nrow(result) < freq) {
    return(NULL)
  }
  
  subset_list <- list()
  for (i in 1:nrow(subset_name)) {
    result_subset <- master6[grep(subset_name[i,], master6$Disease_Code), ]
    age_freq <- as.data.frame(table(floor(result_subset$Diagnosed_age)))
    colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
    
    if (sum(age_freq$"Frequency") > freq) {
      age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)
      age_freq$category <- rep(paste0(result_subset$meaning[1], "(n=",
                                      nrow(result_subset), ")"), nrow(age_freq))
      subset_list[[i]] <- age_freq
    }
  }
  
  totaldataframe <- rbindlist(subset_list)
  
  age_freq <- as.data.frame(table(floor(result$Diagnosed_age)))
  colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
  age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)
  code_meaning <- subset(ICD10_Code, ICD10_Code$coding == code_num)
  age_freq$category <- rep(paste0(code_meaning$meaning, "(n=", nrow(result), ")"),
                           nrow(age_freq))
  
  totaldataframe <- rbind(totaldataframe, age_freq)
  totaldataframe <- rbind(totaldataframe, totalage_freq)
  totaldataframe$`Age at Diagnosis` <- as.numeric(as.character(totaldataframe$`Age at Diagnosis`))
  totaldataframe$category <- as.factor(totaldataframe$category)
  
  return(list(totaldataframe, code_num, "original"))
  
}

##Function that creates incidence rate of subset disease, total disease and total age distribution
# code_num: Disease code from ICD10
original_ICD9 <- function(code_num, freq = 200) {
  
  result1 <- master6[grep(code_num, master6$Disease_Code), ]
  
  geneatlas <- read.table("geneatlas.txt", header = T)
  which_ID <- which(geneatlas$ICD10 == code_num)
  
  if (length(which_ID) > 0) {
    
    ICD9code_num <- geneatlas$ICD9[which_ID]
    result2 <- master6 %>% filter(str_detect(Disease_Code, paste0("^", ICD9code_num)))
    result <- rbind(result1, result2)
    
  } else {
    result <- result1
  }
  
  if (nrow(result) < freq) {
    return(NULL)
  }
  
  age_freq <- as.data.frame(table(floor(result$Diagnosed_age)))
  colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
  age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)
  
  code_meaning <- subset(ICD10_Code, ICD10_Code$coding == code_num)
  age_freq$category <- rep(paste0(code_meaning$meaning, "(n=", nrow(result), ")"),
                           nrow(age_freq))

  totaldataframe <- rbind(totalage_freq, age_freq)
  totaldataframe$`Age at Diagnosis` <- as.numeric(as.character(totaldataframe$`Age at Diagnosis`))
  totaldataframe$category <- as.factor(totaldataframe$category)
  
  return(list(totaldataframe, code_num, "original_w.ICD9"))
  
}

##Function that normalize the incidence rate based on total age incidence rate
normalized <- function(code_num, freq = 200, totalagedf = totalage_freq) {
  
  result <- master6[grep(code_num, master6$Disease_Code), ]
  
  if (nrow(result) < freq) {
    return(NULL)
  }
  
  normalizeddf <- merge(x = original_ICD9(code_num)[[1]], y = totalagedf, by = "Age at Diagnosis", all = TRUE)
  normalizeddf <- na.omit(normalizeddf)
  normalizeddf$`Normalized Incidence Rate` = normalizeddf$`Incidence Rate.x`/normalizeddf$`Incidence Rate.y`
  normalizeddf <- normalizeddf %>% select(1,4,8)
  colnames(normalizeddf) <- c("Age at Diagnosis", "category", "Incidence Rate")
  
  return(list(normalizeddf, code_num, "normalized"))
}

normalized_phecode <- function(code_num, freq = 200) {
  
  code_meaning <- subset(phecode_ukb, phecode_ukb$phecode == code_num)
  icd9_ukb <- unlist(strsplit(code_meaning$icd9_ukb, ", "))
  icd9_ukb <- str_replace_all(icd9_ukb, "[^[:alnum:]]", "")
  
  icd10_ukb <- unlist(strsplit(code_meaning$icd10_ukb, ", "))
  icd10_ukb <- str_replace_all(icd10_ukb, "[^[:alnum:]]", "")
  
  result0 = list()
  for (i in 1:length(icd10_ukb)) {
    result0[[i]] <- master6[grep(icd10_ukb[i], master6$Disease_Code), ]
  }
  df <- do.call(rbind.data.frame, result0)
  
  result1 = list()
  for (i in 1:length(icd9_ukb)) {
    result1[[i]] <- master6 %>% filter(str_detect(Disease_Code, paste0("^", icd9_ukb[i])))
  }
  df1 <- do.call(rbind.data.frame, result1)
  
  result <- rbind(df, df1)
  result <- result[!duplicated(result),]
  
  if (nrow(result) < freq) {
    return(NULL)
  }
  
  age_freq <- as.data.frame(table(floor(result$Diagnosed_age)))
  colnames(age_freq) <- c("Age at Diagnosis", "Frequency")
  age_freq$`Incidence Rate` <- age_freq$Frequency/sum(age_freq$Frequency)
  age_freq$category <- rep(paste0(code_meaning$phenotype, "(n=", nrow(result), ")"),
                           nrow(age_freq))
  
  totaldataframe <- rbind(totalage_freq, age_freq)
  totaldataframe$`Age at Diagnosis` <- as.numeric(as.character(totaldataframe$`Age at Diagnosis`))
  totaldataframe$category <- as.factor(totaldataframe$category)
  
  normalizeddf <- merge(x = totaldataframe, y = totalage_freq, by = "Age at Diagnosis", 
                        all = TRUE)
  normalizeddf <- na.omit(normalizeddf)
  normalizeddf$`Normalized Incidence Rate` = normalizeddf$`Incidence Rate.x`/normalizeddf$`Incidence Rate.y`
  normalizeddf <- normalizeddf %>% select(1,4,8)
  colnames(normalizeddf) <- c("Age at Diagnosis", "category", "Incidence Rate")
  
  return(list(normalizeddf, code_num, "normalized_phecode"))
}

multi_normalized_phecode <- function(code_nums, category_name, freq = 200) {
  
  multi_result <- lapply(c(1:length(code_nums)), function(i) {
    normalized_phecode(code_nums[i], freq)[[1]]
  })
  
  multi_resultdf <- data.table::rbindlist(multi_result)
  
  return(list(multi_resultdf, category_name, "normalized_phecode"))
}