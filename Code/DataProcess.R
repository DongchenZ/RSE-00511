Read.csv.zdc <- function(folder.pattern, colnames){
  current.path <- getwd()
  folder.name <- list.files(pattern = folder.pattern)
  setwd(folder.name)
  file.name <- list.files(pattern="*.csv")
  Data.list <- list()
  for (i in 1:length(file.name)) {
    Data <- read.csv(file.name[i])#Read file
    Data.list[[i]] <- list(file.name = strsplit(file.name[i],".csv")[[1]], Data = as.matrix(Data[,colnames]))
  }
  setwd(current.path)
  return(Data.list)
}
folder.pattern <- "509_Landsat"
colnames <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7")
Data.list <- Read.csv.zdc(folder.pattern = folder.pattern, colnames = colnames)

Sort.name.by.date <- function(string){
  L <- length(string)
  Date <- as.Date(string)
  sorted.Date <- sort(Date)
  Index <- rep(-1,L)
  for (i in 1:L) {
    Index[i] <- which(Date == sorted.Date[[i]])
  }
  return(Index)
}
Date <- c("2019-5-3", "2019-3-7", "2019-9-10", "2020-3-5", "2019-1-1", "2019-12-1")
Index <- Sort.name.by.date(Date)
