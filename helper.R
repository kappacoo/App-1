exampleTable<-data.frame("Sample_ID"=c(1,2,3,4,5,6,7,8), "Group"=c(rep("A",4),rep("B",4)), "Citrate"=runif(8, 1000,12000), "Citrate_13C1"=runif(8, 200, 300), "Citrate_13C2"=c(runif(4, 2000,2100), runif(4,3000,3100)),"Citrate_13C3"=c(runif(4, 1000,1100), runif(4,2000,2100)),"Citrate_13C4"=c(runif(4, 1000,1100), runif(4,2000,2100)),"Citrate_13C5"=c(runif(4, 1000,1100), runif(4,2000,2100)),"Citrate_13C6"=c(runif(4, 1000,1100), runif(4,200,2100)), "LAC"=runif(8, 5000,55000), "LAC_13C1"=runif(8, 200, 3000), "LAC_13C2"=c(runif(4, 200,210), runif(4,20,210)),"LAC_13C3"=c(runif(4, 1000,11000), runif(4, 2000, 21000)))

meta<- read.csv("Molecular Formulas.csv")
meta[,1]<- gsub("-", ".", meta[,1])
meta[,1]<- gsub("/",".", meta[,1])


