library(tidyr)
library(dplyr)
library(data.table)
#We are going to use "Ageatrecruitment" for future use (only 11 discrepancy)
Ageatrecruitment <- read.csv("21022.csv", header = TRUE)
Ageatattending <- read.csv("21003.csv", header = TRUE)
all(Ageatrecruitment$X21022.0.0 == Ageatattending$X21003.0.0) #False
which(Ageatrecruitment$X21022.0.0 != Ageatattending$X21003.0.0) #11 total
agefile <- merge(x = Ageatrecruitment, y = Ageatattending, by = 'eid', all = TRUE)
View(agefile[which(agefile$X21022.0.0 != agefile$X21003.0.0),]) #11 total

#make complete Selfreportedcases file - start with disease code
Noncancercode <- read.csv("20002.csv", header = TRUE)
Cancercode <- read.csv("20001.csv", header = TRUE)

#make Cancercode_1
checkna <- Cancercode[rowSums(is.na(Cancercode)) != 24, ] #45224 total
Cancercode_1 <- Cancercode %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Cancercode")
Cancercode_1 <- Cancercode_1[complete.cases(Cancercode_1), ] #52469 total
length(unique(Cancercode_1$eid)) #45224 total
Cancercode_1 <- Cancercode_1  %>% select(1,3)

#make Noncancercode_1
checkna_1 <- Noncancercode[rowSums(is.na(Noncancercode)) != 136, ] #384906 total 
Noncancercode_1 <- Noncancercode %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Noncancercode")
Noncancercode_1 <- Noncancercode_1[complete.cases(Noncancercode_1), ] #1104882 total
length(unique(Noncancercode_1$eid)) #384906 total
Noncancercode_1 <- Noncancercode_1 %>% separate(fid, c("fid", "Nothing","Num"))
Noncancercode_1$NewID <- paste0(
  Noncancercode_1$eid,".",Noncancercode_1$Nothing,".",Noncancercode_1$Num)
Noncancercode_1 <- Noncancercode_1  %>% select(6,1,5)

#make complete Selfreportedcases file - start with year occurred
Noncancerage_year <- read.csv("20009.csv", header = TRUE)
Cancerage_year <- read.csv("20007.csv", header = TRUE)

#merge Cancercode_1 and Cancerage
View(Cancerage_year[rowSums(is.na(Cancerage_year)) != 24, ]) #45224 total
Cancerage <- Cancerage_year %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Cancerage")
Cancerage <- Cancerage[complete.cases(Cancerage), ] #52469 total
Cancerage <- Cancerage  %>% select(1,3)
Cancerage[] <- lapply(Cancerage, as.integer)
all(Cancercode_1$eid == Cancerage$eid) #True

Cancercode_2 <- cbind(Cancercode_1, Cancerage)
all(Cancercode_2[1] == Cancercode_2[3]) #True
names(Cancercode_2)[3] <- "same"
all(Cancercode_2$eid == Cancercode_2$same) #True
Cancercode_2 <- Cancercode_2  %>% select(1,2,4)
rows <- unique(c(rownames(Cancercode_1), rownames(Cancerage)))
CHECK <- cbind(Cancercode_1[rows ,], Cancerage[rows ,])
all(CHECK[1] == CHECK[3]) #True

#merge Noncancercode_1 and Noncancerage
checkna_2 <- Noncancerage_year[rowSums(is.na(Noncancerage_year)) != 136, ] #384905 total
Noncancerage <- Noncancerage_year %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Noncancerage")
Noncancerage <- Noncancerage[complete.cases(Noncancerage), ] #1104871 total
length(unique(Noncancerage$eid)) #384905 total
Noncancerage <- Noncancerage %>% separate(fid, c("fid", "Nothing","Num"))
Noncancerage$NewID <- paste0(
  Noncancerage$eid,".",Noncancerage$Nothing,".",Noncancerage$Num)
Noncancerage <- Noncancerage  %>% select(6,1,5)
all(Noncancercode_1$NewID == Noncancerage$NewID) #False

Noncancercode_2 <- merge(x = Noncancercode_1, y = Noncancerage, by = 'NewID', all = TRUE)
View(Noncancercode_2[which(is.na(Noncancercode_2$Noncancerage)), ]) #11 total
Noncancercode_2 <- Noncancercode_2  %>% select(2,3,5)
colnames(Noncancercode_2) <- c("eid", "Noncancercode", "Noncancerage")
Noncancercode_2[] <- lapply(Noncancercode_2, as.integer)

#Unused code - cbind failed
#Noncancerage <- dplyr::add_count(Noncancerage, eid, sort = FALSE, name = "n")
#Noncancercode_1 <- dplyr::add_count(Noncancercode_1, eid, sort = FALSE, name = "n")
#rows <- unique(c(rownames(Noncancercode_1), rownames(Noncancerage)))
#CHECK <- cbind(Noncancercode_1[rows ,], Noncancerage[rows ,])
#names(CHECK)[4] <- "same"
#names(CHECK)[6] <- "X"
#Noncancercode_2 <- CHECK %>% select(1,2,5)

#Attach code meaning to the Noncancercode_2 and Cancercode_2
Noncancercode_meaning <- read.csv("coding6.tsv", sep = "\t")
Cancercode_meaning <- read.csv("coding3.tsv", sep = "\t")
names(Noncancercode_meaning)[1] <- "Noncancercode"
names(Cancercode_meaning)[1] <- "Cancercode"

Noncancercode_3 <- merge(x = Noncancercode_2, y = Noncancercode_meaning[1:2], by = 'Noncancercode', all = TRUE)
Noncancercode_3 <- Noncancercode_3  %>% select(2,1,3,4)
Noncancercode_3 <- Noncancercode_3[complete.cases(Noncancercode_3$eid), ]
any(duplicated(Noncancercode_3)) #True

Cancercode_3 <- merge(x = Cancercode_2, y = Cancercode_meaning[1:2], by = 'Cancercode', all = TRUE)
Cancercode_3 <- Cancercode_3  %>% select(2,1,3,4)
Cancercode_3 <- Cancercode_3[complete.cases(Cancercode_3), ]
names(Noncancercode_3)[4] <- "Noncancermeaning"
names(Cancercode_3)[4] <- "Cancermeaning"

#find duplicated data entry and erase it
duplicatedcancer <- which(duplicated(Cancercode_3))
View(Cancercode_3[duplicatedcancer,]) #1720 duplicated
duplicatednoncancer <- which(duplicated(Noncancercode_3))
View(Noncancercode_3[duplicatednoncancer,]) #13185 duplicated
Cancercode_3 <- Cancercode_3[!duplicated(Cancercode_3), ] #50749 total
Noncancercode_3 <- Noncancercode_3[!duplicated(Noncancercode_3), ] #1091697 total

#erase all the ages in the Noncancercode_3 and Cancercode_3 that are -1, -3 or NA
View(Cancercode_3[Cancercode_3$Cancerage == -1, ]) #106 total (date unknown)
View(Cancercode_3[Cancercode_3$Cancerage == -3, ]) #none
View(Cancercode_3[which(is.na(Cancercode_3$Cancerage)), ]) #none
View(Cancercode_3[Cancercode_3$Cancerage != -1, ]) #50643 total
Cancercode_3 <- Cancercode_3[Cancercode_3$Cancerage != -1, ]
any(Cancercode_3$Cancerage < 0) #none
length(which(Cancercode_3$Cancerage == 0)) #17 total
View(Cancercode_3[Cancercode_3$Cancerage == 0, ]) #17 total
saveRDS(Cancercode_3, file = "SelfreportedCancer.rds")

length(which(Noncancercode_3$Noncancerage == -1)) #36348 total (date unknown)
length(which(Noncancercode_3$Noncancerage == -3)) #235 total (prefer not to answer)
View(Noncancercode_3[which(is.na(Noncancercode_3$Noncancerage)), ]) #11 total
length(which(is.na(Noncancercode_3$Noncancerage))) #11 total
View(Noncancercode_3[Noncancercode_3$Noncancerage != -1, ]) #1055349 total
Noncancercode_4 <- Noncancercode_3[Noncancercode_3$Noncancerage != -1, ]
View(Noncancercode_4[Noncancercode_4$Noncancerage != -3, ]) #1055114 total
Noncancercode_4 <- Noncancercode_4[Noncancercode_4$Noncancerage != -3, ]
which(is.na(Noncancercode_4$Noncancerage)) #11 total
View(Noncancercode_4[complete.cases(Noncancercode_4), ]) #1055103 total
Noncancercode_4 <- Noncancercode_4[complete.cases(Noncancercode_4), ]
any(Noncancercode_4$Noncancerage < 0) #none
length(which(Noncancercode_4$Noncancerage == 0)) #6981 total
View(Noncancercode_4[Noncancercode_4$Noncancerage== 0, ]) #6981 total
saveRDS(Noncancercode_4, file = "SelfreportedNoncancer.rds")

#Selfreportedcases <- merge(x = Noncancercode_4, y = Cancercode_4, by = 'eid', all = TRUE)
#any(duplicated(Selfreportedcases))
#saveRDS(Selfreportedcases, file = "Selfreportedcases.rds")