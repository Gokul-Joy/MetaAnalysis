
library(readxl)

data <- read_excel("GC_DFS_NMA.xlsx",na="NA")
summary(data)

colSums(is.na(data))
