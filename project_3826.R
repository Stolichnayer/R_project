# # # # # # # # # # # # # # # # # #
#                                 #
#   Alexandros Perrakis csd3826   #
#         Final Project           #
#                                 #
# # # # # # # # # # # # # # # # # #


#1. Read data from the cleaned file
raw.data <- read.table("C:\\Users\\Alex\\Desktop\\GDS4879.clean", sep="\t", header=TRUE)
dim(raw.data)

#The first two columns contain the ID Ref and the identifie, so we remove them
data <- raw.data[,-c(1,2)]

#Lusi sto provlima me ta NA kai ta factor-values (apo to elearn) opws proteinan merika paidia
data = data.frame(sapply(data, function(x){
  if(is.factor(x)){
    as.numeric(as.character(x))
  }else{
    x
  }
}))

#2. Creating boxplot
boxplot(data)

data_log2 <- lapply(data, function(x) log2(as.numeric(as.character(x))) )

boxplot(data_log2)
#Creating another plot with average differences
difAverages <- function(v, female.alc.index, female.con.index){
  mean(v[female.alc.index]) - mean(v[female.con.index])
}

average.differences <- apply(data, 1, difAverages, 1:6, 7:12)

### The greatest difference
greatest.dif.index <- which.max(average.differences)
labels <- as.factor(c(rep("female alcoholics", 6), rep("female control", 6)))


data.plot <- data.frame(labels=labels,values=as.numeric(data[greatest.dif.index,]))
plot(data.plot$values~data.plot$labels)



#3.
#1:6	female alcoholic
#7:12	female control
#13:26 male alcoholic
#26:39 male control

#Opote: alcoholics = 1:6,  13:26
#       controls   = 7:12, 26:39

geneAverageDataNotContinuous <- function(gene_data, index1, index2){
   (mean(gene_data[index1]) +  mean(gene_data[index2])) / 2
}

geneAverageDataContinuous <- function(gene_data, index){
  mean(gene_data[index])
}

averageAlcoholics <- apply(data, 1, geneAverageDataNotContinuous, 1:6, 13:26)
averageControl <- apply(data, 1, geneAverageDataNotContinuous, 7:12, 27:39)

averageFemale  <- apply(data, 1, geneAverageDataContinuous, 1:12)
averageMale  <- apply(data, 1, geneAverageDataContinuous, 13:39)

#remove NA values
averageAlcoholics<-averageAlcoholics[!is.na(averageAlcoholics)]
averageControl<-averageControl[!is.na(averageControl)]
averageFemale<-averageFemale[!is.na(averageFemale)]
averageMale<-averageMale[!is.na(averageMale)]

# add to list genes that are higher in alcoholics
counter = 1
my_list <- list();
for(i in 1:length(averageAlcoholics)){
  if (averageAlcoholics[i] < averageControl[i]){
      my_list[counter] <- raw.data[i:i, 1:1]
      counter = counter + 1
  }
    
}

# add to list genes that are lower in Males
counter = 1
my_list2 <- list();
for(i in 1:length(averageFemale)){
  if (averageFemale[i] > averageMale[i]){
    my_list2[counter] <- raw.data[i:i, 1:1]
    counter = counter + 1
  }
  
}

pvalues_1 <- sapply(1:length(averageAlcoholics), function(x){
  t.test(averageAlcoholics, averageControl,"greater")$p.value
})

pvalues_2 <- sapply(1:length(averageFemale), function(x){
  t.test(averageFemale, averageMale,"greater")$p.value
})

#5. FDR correction
fdrpvalues = p.adjust(pvalues_1,method="BH")
fdr <- which(fdrpvalues < 0.1)
fdr

#Bonferroni correction
bonferronipvalues=p.adjust(pvalues_2,method="bonferroni")
bonferroni <- which(bonferronipvalues < 0.1)
bonferroni

#6.
gprofiler(c("NameOfGenes"), organism = "hsapiens")

