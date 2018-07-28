## Loading the package and input file
library(vegan)
data <- read.csv('No42_Pinaceae_Li2013_mt.csv', header = TRUE)
Region <- data[,1]
dat <- data[,-1]

## Count the raw number of haplotype (Nhaps)
Nhaps <- round(specnumber(dat), 2)

## Rarefaction in Haplotype
HAP_R <- round(rarefy(dat, 2), 2)

## Calculating gene diversity(haplotype diversity) H
F <-function(d) {
  f <- rowSums(d) / (rowSums(d)-1)
  Simpson <- diversity(d, 'simpson')
  H1 <- round (f * Simpson, 2)
}

H <- F(dat)

## Prepare the outfile and write in the csv
outfile <- data.frame(Region, Nhaps, HAP_R, H)
write.csv(outfile, 'No42_mt.csv')

#################################################
## Change to Gengepop format
# Extract the matrix we need
HAP <- data[1,-1]
# Build the genepop sequencing
N <- 1100 + 100 * length(HAP) - 100
gene <- seq(1100, N, 100)

# change to the genepop format
mydata2 <- list()
for(i in 1:nrow(dat)){
  mydata1 <- list()
  for(j in 1:ncol(dat)){
    value <- dat[i,j]
    if (value == 0) next
    x <- rep(paste(colnames(dat[j]), gene[j], sep = ','), value)
    mydata1[[j]] <- x}
  mydata2[[i]] <- c('pop', mydata1)}
y <- c('Title:', 'Locus_1',unlist(mydata2))
write.table(y, 'No42_mt.gpp', quote = FALSE, sep = ' ', row.names = FALSE, col.names = FALSE)