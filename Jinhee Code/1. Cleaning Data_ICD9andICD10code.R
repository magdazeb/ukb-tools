library(tidyr)
library(dplyr)

ICD10 <- read.csv("41270.csv", header = TRUE)
any(duplicated(ICD10$eid))
ICD10_Date <- read.csv("41280.csv", header = TRUE)
any(duplicated(ICD10_Date$eid))
View(ICD10)
View(ICD10_Date)
ICD10 <- ICD10 %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD10")
ICD10 <- ICD10[!(ICD10$ICD10==""),]
ICD10 <- ICD10 %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD10[,3])
dim(table(ICD10[,3])) == 1
ICD10$Nothing <- NULL
NewID <- paste0(ICD10$eid,".",ICD10$Num)
ICD10 <- data.frame(NewID = c(NewID), ICD10)
ICD10$Num <- NULL
ICD10_Date <- ICD10_Date %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "date")
ICD10_Date <- ICD10_Date[!(ICD10_Date$date==""),]
ICD10_Date <- ICD10_Date %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD10_Date[,3])
dim(table(ICD10_Date[,3])) == 1
ICD10_Date$Nothing <- NULL
NewID <- paste0(ICD10_Date$eid,".",ICD10_Date$Num)
ICD10_Date <- data.frame(NewID = c(NewID), ICD10_Date)
ICD10_Date$Num <- NULL
ICD10FullData <- merge(x = ICD10, y = ICD10_Date, by = 'NewID', all = TRUE)
View(ICD10FullData)
ICD10FullData <- na.omit(ICD10FullData)
all(ICD10FullData$eid.x == ICD10FullData$eid.y)
ICD10FullData$eid.y <- NULL

ICD10_Main <- read.csv("41202.csv", header = TRUE)
any(duplicated(ICD10_Main$eid))
ICD10_MainDate <- read.csv("41262.csv", header = TRUE)
any(duplicated(ICD10_MainDate$eid))
View(ICD10_Main)
View(ICD10_MainDate)
ICD10_Main <- ICD10_Main %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD10_Main")
ICD10_Main <- ICD10_Main[!(ICD10_Main$ICD10_Main==""),]
ICD10_Main <- ICD10_Main %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD10_Main[,3])
dim(table(ICD10_Main[,3])) == 1
ICD10_Main$Nothing <- NULL
NewID <- paste0(ICD10_Main$eid,".",ICD10_Main$Num)
ICD10_Main <- data.frame(NewID = c(NewID), ICD10_Main)
ICD10_Main$Num <- NULL
ICD10_MainDate <- ICD10_MainDate %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "date")
ICD10_MainDate <- ICD10_MainDate[!(ICD10_MainDate$date==""),]
ICD10_MainDate <- ICD10_MainDate %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD10_MainDate[,3])
dim(table(ICD10_MainDate[,3])) == 1
ICD10_MainDate$Nothing <- NULL
NewID <- paste0(ICD10_MainDate$eid,".",ICD10_MainDate$Num)
ICD10_MainDate <- data.frame(NewID = c(NewID), ICD10_MainDate)
ICD10_MainDate$Num <- NULL
ICD10MainFullData <- merge(x = ICD10_Main, y = ICD10_MainDate, by = 'NewID', all = TRUE)
View(ICD10MainFullData)
ICD10MainFullData <- na.omit(ICD10MainFullData)
all(ICD10MainFullData$eid.x == ICD10MainFullData$eid.y)
ICD10MainFullData$eid.y <- NULL

ICD10FullData$NewID <- paste0(ICD10FullData$eid.x,".",ICD10FullData$ICD10,".",ICD10FullData$date)
any(duplicated(ICD10FullData$NewID))
ICD10MainFullData$NewID <- paste0(ICD10MainFullData$eid.x,".",ICD10MainFullData$ICD10_Main,".",ICD10MainFullData$date)
any(duplicated(ICD10MainFullData$NewID))
ICD10Data <- merge(x = ICD10FullData, y = ICD10MainFullData, by = 'NewID', all = TRUE)
View(ICD10Data)

ICD10DataCheck <- which(is.na(ICD10Data$eid.x.x))
ICD10DataCheck <- ICD10Data[ICD10DataCheck,]
length(unique(ICD10DataCheck$eid.x.y))
ICD10DataCheck_Main <- which(is.na(ICD10Data$eid.x.y))
ICD10DataCheck_Main <- ICD10Data[ICD10DataCheck_Main,]
ICD10DataCheck$NewID0 <- paste0(ICD10DataCheck$eid.x.y,".",ICD10DataCheck$ICD10_Main)
ICD10DataCheck_Main$NewID0 <- paste0(ICD10DataCheck_Main$eid.x.x,".",ICD10DataCheck_Main$ICD10)
ICD10DataCheck <- merge(x = ICD10DataCheck_Main, y = ICD10DataCheck, by = 'NewID0', all = TRUE)
View(ICD10DataCheck)
ICD10DataCheck <- ICD10DataCheck[which(!(is.na(ICD10DataCheck$date.y.y))),]
ICD10DataCheck <- ICD10DataCheck[, -c(8:12)]
ICD10DataCheck <- ICD10DataCheck[, -c(9:13)]
all(ICD10DataCheck$NewID.x < ICD10DataCheck$NewID.y)

sum(is.na(ICD10Data$eid.x.x))
ICD10Data <- ICD10Data[!(is.na(ICD10Data$eid.x.x)),]
ICD10Data$content <- ifelse(ICD10Data$NewID %in% ICD10DataCheck$NewID.x == TRUE, "Main", "Secondary")
length(which(ICD10Data$content == "Main"))
ICD10Data$content[ICD10Data$eid.x.y >= 0] <- "Main"
length(which(ICD10Data$content == "Main"))
ICD10Data <- ICD10Data[, -c(7:11)]

ICD10_Secondary <- read.csv("41204.csv", header = TRUE)
View(ICD10_Secondary)
ICD10_Secondary <- ICD10_Secondary %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD10_Secondary")
ICD10_Secondary <- ICD10_Secondary[!(ICD10_Secondary$ICD10_Secondary==""),]
ICD10_Secondary <- ICD10_Secondary %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD10_Secondary[,3])
dim(table(ICD10_Secondary[,3])) == 1
ICD10_Secondary$Nothing <- NULL
NewID <- paste0(ICD10_Secondary$eid,".",ICD10_Secondary$Num)
ICD10_Secondary <- data.frame(NewID = c(NewID), ICD10_Secondary)
ICD10_Secondary$Num <- NULL

ICD9 <- read.csv("41271.csv", header = TRUE)
any(duplicated(ICD9$eid))
ICD9_Date <- read.csv("41281.csv", header = TRUE)
any(duplicated(ICD9_Date$eid))
View(ICD9)
View(ICD9_Date)
ICD9 <- ICD9 %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD9", values_ptypes = list(ICD9 = 'character'))
ICD9 <- ICD9[!(ICD9$ICD9==""),]
ICD9 <- na.omit(ICD9)
ICD9 <- ICD9 %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD9[,3])
dim(table(ICD9[,3])) == 1
ICD9$Nothing <- NULL
NewID <- paste0(ICD9$eid,".",ICD9$Num)
ICD9 <- data.frame(NewID = c(NewID), ICD9)
ICD9$Num <- NULL
ICD9_Date <- ICD9_Date %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "date")
ICD9_Date <- ICD9_Date[!(ICD9_Date$date==""),]
ICD9_Date <- ICD9_Date %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD9_Date[,3])
dim(table(ICD9_Date[,3])) == 1
ICD9_Date$Nothing <- NULL
NewID <- paste0(ICD9_Date$eid,".",ICD9_Date$Num)
ICD9_Date <- data.frame(NewID = c(NewID), ICD9_Date)
ICD9_Date$Num <- NULL
ICD9FullData <- merge(x = ICD9, y = ICD9_Date, by = 'NewID', all = TRUE)
View(ICD9FullData)
any(is.na(ICD9FullData))
all(ICD9FullData$eid.x == ICD9FullData$eid.y)
ICD9FullData$eid.y <- NULL

ICD9_Main <- read.csv("41203.csv", header = TRUE)
any(duplicated(ICD9_Main$eid))
ICD9_MainDate <- read.csv("41263.csv", header = TRUE)
any(duplicated(ICD9_MainDate$eid))
View(ICD9_Main)
View(ICD9_MainDate)
ICD9_Main <- ICD9_Main %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD9_Main", values_ptypes = list(ICD9_Main = 'character'))
ICD9_Main <- ICD9_Main[!(ICD9_Main$ICD9_Main==""),]
ICD9_Main <- na.omit(ICD9_Main)
ICD9_Main <- ICD9_Main %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD9_Main[,3])
dim(table(ICD9_Main[,3])) == 1
ICD9_Main$Nothing <- NULL
NewID <- paste0(ICD9_Main$eid,".",ICD9_Main$Num)
ICD9_Main <- data.frame(NewID = c(NewID), ICD9_Main)
ICD9_Main$Num <- NULL
ICD9_MainDate <- ICD9_MainDate %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "date")
ICD9_MainDate <- ICD9_MainDate[!(ICD9_MainDate$date==""),]
ICD9_MainDate <- ICD9_MainDate %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD9_MainDate[,3])
dim(table(ICD9_MainDate[,3])) == 1
ICD9_MainDate$Nothing <- NULL
NewID <- paste0(ICD9_MainDate$eid,".",ICD9_MainDate$Num)
ICD9_MainDate <- data.frame(NewID = c(NewID), ICD9_MainDate)
ICD9_MainDate$Num <- NULL
ICD9MainFullData <- merge(x = ICD9_Main, y = ICD9_MainDate, by = 'NewID', all = TRUE)
View(ICD9MainFullData)
any(is.na(ICD9MainFullData))
all(ICD9MainFullData$eid.x == ICD9MainFullData$eid.y)
ICD9MainFullData$eid.y <- NULL

ICD9FullData$NewID <- paste0(ICD9FullData$eid.x,".",ICD9FullData$ICD9,".",ICD9FullData$date)
any(duplicated(ICD9FullData$NewID))
ICD9MainFullData$NewID <- paste0(ICD9MainFullData$eid.x,".",ICD9MainFullData$ICD9_Main,".",ICD9MainFullData$date)
any(duplicated(ICD9MainFullData$NewID))
ICD9Data <- merge(x = ICD9FullData, y = ICD9MainFullData, by = 'NewID', all = TRUE)
View(ICD9Data)

ICD9DataCheck <- which(is.na(ICD9Data$eid.x.x))
ICD9DataCheck <- ICD9Data[ICD9DataCheck,]
length(unique(ICD9DataCheck$eid.x.y))
ICD9DataCheck_Main <- which(is.na(ICD9Data$eid.x.y))
ICD9DataCheck_Main <- ICD9Data[ICD9DataCheck_Main,]
ICD9DataCheck$NewID0 <- paste0(ICD9DataCheck$eid.x.y,".",ICD9DataCheck$ICD9_Main)
ICD9DataCheck_Main$NewID0 <- paste0(ICD9DataCheck_Main$eid.x.x,".",ICD9DataCheck_Main$ICD9)
ICD9DataCheck <- merge(x = ICD9DataCheck_Main, y = ICD9DataCheck, by = 'NewID0', all = TRUE)
View(ICD9DataCheck)
ICD9DataCheck <- ICD9DataCheck[which(!(is.na(ICD9DataCheck$date.y.y))),]
ICD9DataCheck <- ICD9DataCheck[, -c(8:12)]
ICD9DataCheck <- ICD9DataCheck[, -c(9:13)]
all(ICD9DataCheck$NewID.x < ICD9DataCheck$NewID.y)

sum(is.na(ICD9Data$eid.x.x))
ICD9Data <- ICD9Data[!(is.na(ICD9Data$eid.x.x)),]
ICD9Data$content <- ifelse(ICD9Data$NewID %in% ICD9DataCheck$NewID.x == TRUE, "Main", "Secondary")
length(which(ICD9Data$content == "Main"))
ICD9Data$content[ICD9Data$eid.x.y >= 0] <- "Main"
length(which(ICD9Data$content == "Main"))
ICD9Data <- ICD9Data[, -c(7:11)]

ICD9_Secondary <- read.csv("41205.csv", header = TRUE)
View(ICD9_Secondary)
ICD9_Secondary <- ICD9_Secondary %>% pivot_longer(cols = -c(eid), names_to = "fid", values_to = "ICD9_Secondary", values_ptypes = list(ICD9_Secondary = 'character'))
ICD9_Secondary <- ICD9_Secondary[!(ICD9_Secondary$ICD9_Secondary==""),]
ICD9_Secondary <- na.omit(ICD9_Secondary)
ICD9_Secondary <- ICD9_Secondary %>% separate(fid, c("fid", "Nothing","Num"))
table(ICD9_Secondary[,3])
dim(table(ICD9_Secondary[,3])) == 1
ICD9_Secondary$Nothing <- NULL
NewID <- paste0(ICD9_Secondary$eid,".",ICD9_Secondary$Num)
ICD9_Secondary <- data.frame(NewID = c(NewID), ICD9_Secondary)
ICD9_Secondary$Num <- NULL
ICD9_Secondary <- ICD9_Secondary[order(ICD9_Secondary$NewID),]

ICD9Data$NewID0 <- paste0(ICD9Data$eid.x.x,".",ICD9Data$date.x)
ICD10Data$NewID0 <- paste0(ICD10Data$eid.x.x,".",ICD10Data$date.x)
duplicatedata <- merge(x = ICD10Data, y = ICD9Data, by = 'NewID0', all = TRUE)
View(duplicatedata)
sum(is.na(duplicatedata$eid.x.x.x))
duplicatedata <- duplicatedata[!(is.na(duplicatedata$eid.x.x.x)),]
duplicatedata <- duplicatedata[!(is.na(duplicatedata$eid.x.x.y)),]
View(duplicatedata)

PatientData <- bind_rows("ICD10" = ICD10Data, "ICD9" = ICD9Data, .id = "ICD9/ICD10")
length(unique(PatientData$NewID))
PatientData <- PatientData %>% separate(NewID, c("Patient_ID", "Disease_Code","Disease_Date"), extra = "merge", fill = "left")
sum(is.na(PatientData$Disease_Code))
all(PatientData$Patient_ID == PatientData$eid.x.x)
all(PatientData$Disease_Date == PatientData$date.x)

PatientData <- select(PatientData, c(2,1,3,4,10))
PatientData <- PatientData[order(PatientData$Patient_ID),]
dim(PatientData)
View(PatientData)

library(data.table)
ICD9Data <- select(ICD9Data, c(2:7))
ICD9Data <- setnames(ICD9Data, old = c('eid.x.x','fid.x.x', 'fid.y.x', 'date.x'), 
                     new = c('patientID','fid_DiseaseCode', 'fid_Date', 'date'))
View(ICD9Data)

ICD10Data <- select(ICD10Data, c(2:7))
ICD10Data <- setnames(ICD10Data, old = c('eid.x.x','fid.x.x', 'fid.y.x', 'date.x'), 
                      new = c('patientID','fid_DiseaseCode', 'fid_Date', 'date'))
View(ICD10Data)

write.csv(ICD9Data, row.names=FALSE, file = 'ICD9Data.csv')
write.csv(ICD10Data, row.names=FALSE, file = 'ICD10Data.csv')
write.csv(PatientData, row.names=FALSE, file = 'PatientData.csv')