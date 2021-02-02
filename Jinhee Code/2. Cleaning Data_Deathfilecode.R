library(tidyr)
library(dplyr)
library(data.table)

#check master4 age at death column
master4 <- readRDS("master_4.rds")
prevageatdeath <- master4 %>% select(1,12)
prevageatdeath <- prevageatdeath[complete.cases(prevageatdeath), ]
View(prevageatdeath[!duplicated(prevageatdeath), ]) #19909 total

#make complete deathfile - start
Ageatdeath <- read.csv("40007.csv", header = TRUE)
Primarycauseofdeath <- read.csv("40001.csv", header = TRUE)
Secondarycauseofdeath <- read.csv("40002.csv", header = TRUE)
Descriptioncauseofdeath <- read.csv("40010.csv", header = TRUE)

#deathfile: first column - Ageatdeath (death1)
View(Ageatdeath[complete.cases(Ageatdeath$X40007.0.0), ]) #20465 total
View(Ageatdeath[rowSums(is.na(Ageatdeath)) != 3, ]) #20465 total
Ageatdeath <- Ageatdeath %>% select(1,2)
View(Ageatdeath[rowSums(is.na(Ageatdeath)) != 1, ]) #20465 total

death1 <- Ageatdeath %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Ageatdeath")
View(death1[complete.cases(death1), ]) #20465 total
death1 <- death1[complete.cases(death1), ]
death1 <- death1 %>% select(1,3)
death1[] <- lapply(death1, as.integer)

#deathfile: second column - Primarycauseofdeath (death2)
Primarycauseofdeath <- mutate_all(Primarycauseofdeath, list(~na_if(.,"")))
View(Primarycauseofdeath[complete.cases(Primarycauseofdeath$X40001.0.0), ]) #20432 total
View(Primarycauseofdeath[rowSums(is.na(Primarycauseofdeath)) != 2, ]) #20432 total
checkduplicate <- Primarycauseofdeath[rowSums(is.na(Primarycauseofdeath)) != 2, ]
which(is.na(checkduplicate$X40001.0.0)) #none
View(checkduplicate[rowSums(is.na(checkduplicate)) != 1, ]) #60 total
View(checkduplicate[which(as.character(checkduplicate$X40001.0.0) == 
                       as.character(checkduplicate$X40001.1.0)), ]) #19 duplicate

death2 <- Primarycauseofdeath %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Primarycauseofdeath")
View(death2[complete.cases(death2), ]) #20492 total
death2 <- death2[complete.cases(death2), ]
death2 <- death2 %>% select(1,3)
View(death2[!duplicated(death2), ]) #20473 total
death2 <- death2[!duplicated(death2), ]
uniquedeath2 <- death2[!duplicated(death2$eid), ] #20432 total

#deathfile: third column - Secondarycauseofdeath (death3)
Secondarycauseofdeath <- mutate_all(Secondarycauseofdeath, list(~na_if(.,"")))
checkna <- Secondarycauseofdeath[rowSums(is.na(Secondarycauseofdeath)) != 28, ] #12660 total
which(!is.na(checkna$X40002.0.13)) #one ID
which(!is.na(checkna$X40002.1.8))  #one ID
which(!is.na(checkna$X40002.1.9))  #none
View(Secondarycauseofdeath %>% select(1:24))
Secondarycauseofdeath <- Secondarycauseofdeath %>% select(1:24)
View(Secondarycauseofdeath[rowSums(is.na(Secondarycauseofdeath)) != 23, ]) #12660 total
checkna <- Secondarycauseofdeath[rowSums(is.na(Secondarycauseofdeath)) != 23, ]
View(checkna[which(is.na(checkna$X40002.0.0)), ]) #24 total

death3 <- Secondarycauseofdeath %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Secondarycauseofdeath")
View(death3[complete.cases(death3), ]) #27182 total
death3 <- death3[complete.cases(death3), ] 
View(death3[!duplicated(death3$eid), ]) #12660 total
death3 <- death3 %>% select(1,3)
View(death3[which(duplicated(death3)), ]) #65 duplicated
View(death3[!duplicated(death3), ]) #27117 total
death3 <- death3[!duplicated(death3), ]
uniquedeath3 <- death3[!duplicated(death3$eid), ] #12660 total

#deathfile: fourth column - Descriptioncauseofdeath (death4)
Descriptioncauseofdeath <- mutate_all(Descriptioncauseofdeath, list(~na_if(.,"")))
View(Descriptioncauseofdeath[complete.cases(Descriptioncauseofdeath), ]) #17147 total
any(duplicated(Descriptioncauseofdeath$eid)) #none

death4 <- Descriptioncauseofdeath %>% 
  pivot_longer(cols = -c(eid), names_to = "fid", values_to = "Descriptioncauseofdeath")
View(death4[complete.cases(death4), ]) #17147 total
death4 <- death4 %>% select(1,3)
death4 <- death4[complete.cases(death4), ]

#merge death1, death2, death3, death4
deathfile <- merge(x = death2, y = death3, by = 'eid', all = TRUE)
colnames(deathfile) <- c("eid", "Primary", "Secondary")
deathfile <- deathfile %>% 
  pivot_longer(cols = -c(eid), names_to = "Causeofdeath", values_to = "Diseasecode")
uniquedeath <-  merge(x = uniquedeath2, y = uniquedeath3, by = 'eid', all = TRUE)
all(is.na(uniquedeath$Primarycauseofdeath) == FALSE) #True
length(which(is.na(uniquedeath$Secondarycauseofdeath))) #7772 total
View(deathfile[!duplicated(deathfile), ]) #55362 total
deathfile <- deathfile[!duplicated(deathfile), ]
length(which(is.na(deathfile$Diseasecode))) #7772 total
View(deathfile[complete.cases(deathfile), ]) #47590 total
deathfile <- deathfile[complete.cases(deathfile), ]

deathfile <- merge(x = deathfile, y = death1, by = 'eid', all = TRUE)
length(which(is.na(deathfile$Causeofdeath))) #33 total
length(which(is.na(deathfile$Diseasecode))) #33 total

deathfile <- merge(x = deathfile, y = death4, by = 'eid', all = TRUE)
deathfile <- deathfile %>% select(1,4,2,3,5) #47623 total
length(unique(deathfile$eid)) #20465 total
saveRDS(deathfile, file = "Deathfile.rds")
deathfile <- readRDS("Deathfile.rds")

#add new ageatdeath to master4 to make master5
master5 <- merge(x = master4[1:11], y = death1, by = 'eid', all = TRUE)
afterageatdeath <- master5 %>% select(1,12)
afterageatdeath <- afterageatdeath[complete.cases(afterageatdeath), ]
View(afterageatdeath[!duplicated(afterageatdeath), ]) #20465 total

any(duplicated(master5)) #this code failed :(
length(unique(master4$eid)) #413589 total
length(unique(master5$eid)) #414145 total
uniquemaster4 <- master4[!duplicated(master4$eid), ]
uniquemaster5 <- master5[!duplicated(master5$eid), ]
uniquemaster <-  merge(x = uniquemaster5, y = uniquemaster4, by = 'eid', all = TRUE)
#which(uniquemaster$Disease_Code.x != uniquemaster$Disease_Code.y)
any(is.na(master5$Disease_Code)) #True
View(master5[which(is.na(master5$Disease_Code)), ]) #556 total
saveRDS(master5, file = "master_5.rds")

#df <- data.table::rbindlist(db1)
#df <- df %>% select(1:5)
#saveRDS(df, file = "Deathfile.rds")