#analysis for the mob score

rm(list= ls())
setwd("~/Desktop/Dist_matrix_TE_div")

spec1 <- "Human"
if(spec1 == "Mouse"){UCSCspec = "mm10"}
if(spec1 == "Human"){UCSCspec = "hg19"}
if(spec1 == "Opossum"){UCSCspec = "monDom5"}


yn <- c("no")
for(i in yn){
print(i)
# file selection 
keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = i

file_name <- paste(spec1, "_R.cor.cons", "_NGenes.", keep.NGenes ,"_NG4s.", keep.NG4s, "_NCpGI.", keep.NCpGI, "_CpGBP.", keep.CpGBP, "_GC.", keep.GC, "_SCALE.", SCALE,"_no_conservation.txt" ,sep="")

# read in the table and seperate out conservation, TE dist, and species comparison
# make sure each one segregate with the binID

# maybe get a list of usefule species too 


all.species1 <- read.table(paste("r_class_table/", file_name, sep = ""), header = TRUE)
# if we do a heatmap with all the TE info and all the scores in it
scores <- all.species1[,grep("dist", colnames(all.species1))]
TE <- all.species1[,5:18]


all <- c(scores[,1], scores[,2], scores[,3], scores[,4])

if(SCALE == "no"){
	pca <- prcomp(scale(TE))
}
if(SCALE == "yes"){
	pca <- prcomp((TE))
}


score.bov <- scores$Horse_dist_cor

A <- boxplot(score.bov)

cols <- rep(1, length(score.bov))
cols[score.bov < A$stats[2,]] <- 2
cols[score.bov > A$stats[4,]] <- 3

plot(pca$x[,1], pca$x[,2], col = cols )
assign(i, cols)

}









