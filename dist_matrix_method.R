# get the S bins and calculate the distance matricies

# all options need to be at the top

# this script just does S1 to S2

rm(list= ls())
setwd("~/Desktop/Human_and_Mouse_TE.div/")

library(gplots)

#####
# Load in Tables and process

spec1 <- "Human"
spec2 <- "Mouse"
UCSCspec= "hg19"


rem.un <- "yes"

mb <- 1000000
bin.size = 500000

#which columns to keep 
#keep.NGenes = "yes"
#keep.NG4s = "no"
#keep.NCpGI = "yes"
#keep.CpGBP = "no"
#keep.GC = "yes"
#SCALE = "yes"
# trying out the unscaled method
keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "yes"

#create objects into whcih i will store the binsizes
s1name <- paste(spec1, "_AllBinCounts.txt", sep = "")
s2name <- paste(spec2, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
s2 <- read.table(s2name, header = TRUE)
slist <- list(s1,s2)


for(i in seq(along=slist)){
      count <- slist[[i]]
      count <- count[count$Known >= bin.size,]
      count[,5:(length(count)-1)] <- (count[,5:(length(count)-1)]/count$Known) * mb   
      if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}
      if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
      if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
      if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
      if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
      if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
      #count <- count[,!(colnames(count) == "Known")]
      colnames(count)[1:4] <- c("chr", "binID", "start", "end")
      count$binID <- 1:dim(count)[1]
      slist[[i]] <- count
}

KnownS1 <- data.frame(slist[[1]]$binID, slist[[1]]$Known)
KnownS2	<- data.frame(slist[[2]]$binID,	slist[[2]]$Known)

s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]
s2 <- slist[[2]][,!(colnames(slist[[2]]) == "Known")]

if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}
if(rem.un == "yes"){
	if(length(grep("U", s2$chr)) > 0){s2 <- s2[-(grep( "U", s2$chr)),]}
	if(length(grep("_", s2$chr)) > 0){s2 <- s2[-(grep("_", s2$chr)),]}
	if(length(grep("M", s2$chr)) > 0){s2 <- s2[-(grep("M", s2$chr)),]}
	}
	
# load in Sbin	
S.bin.name <- paste(spec1, "_aligning_", spec2, "_fix_count", sep = "")
S.bin_s1_s2 <- read.table(S.bin.name,header = TRUE)

# make sure each bin has the df has the same bins
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2$S1.bin %in% s1$binID ,]
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2$S2.bin %in% s2$binID ,]
s1 <- s1[s1$binID %in% S.bin_s1_s2$S1.bin,]



s2TE <- s2[,c(2, 5:length(s2))]
merg.s1_s2 <- merge(S.bin_s1_s2, s2TE, by.x=3, by.y=1)	

# we have all the s2 coordinates on the merg table and which S2 bins and S1 bins they correspond to
# my understanding is that we want to know what are our s1 coordinates in S2
# we do this by identifying which bins are the same as the S1 bins

S2_TE_name <- colnames(s2TE)[2:length(s2TE)] 
for( i in S2_TE_name){
	merg.s1_s2[,i] <- merg.s1_s2[,i]*merg.s1_s2$S2.Proportion
}
# now we sum them
# then we sum the S1.proportion and divide em
newS2.coord <- matrix(data=rep(0, (length(s2TE)-1)*length(unique(merg.s1_s2$S1.bin)) ), ncol = (length(s2TE)-1), nrow = length(unique(merg.s1_s2$S1.bin)) )
colnames(newS2.coord) <- S2_TE_name
replot_s1_s2 <- data.frame(S1_bin = unique(merg.s1_s2$S1.bin),newS2.coord)
for(i in seq(dim(replot_s1_s2)[1])){	
	divide <- sum(merg.s1_s2[replot_s1_s2[i,'S1_bin'] == merg.s1_s2[,'S1.bin'],'S1.Proportion'])
	replot_s1_s2[i,2:length(replot_s1_s2)] <- colSums(merg.s1_s2[replot_s1_s2[i,'S1_bin'] == merg.s1_s2[,'S1.bin'],5:length(merg.s1_s2)]) / divide
}
replot_s1_s2 <- replot_s1_s2[order(replot_s1_s2[,'S1_bin']),]
# therefore replot_s1_s2 should be the S2 coordinates of the S1 bins

all(replot_s1_s2$S1_bin == s1$binID)

rownames(replot_s1_s2) <- 1:dim(replot_s1_s2)[1]
rownames(s1) <- 1:dim(s1)[1]

#Calculate distance matrix
# and get the pairwise_dist cor

S1_dist <- dist(s1[,5:length(s1)])
S1_s2_dist <- dist(replot_s1_s2[,2:length(replot_s1_s2)])

S1.dist.m <- as.matrix(S1_dist)
S1_s2.dist.m <- as.matrix(S1_s2_dist)
cor.all <- cor(S1.dist.m,S1_s2.dist.m)
cor.diag <- rep(0, length(rownames(cor.all)))
for(i in 1:length(rownames(cor.all))){
	cor.diag[i] <- cor.all[i,i]
}
s1$archi.cor.scale <- scale(cor.diag)

# seperate genome into contigous sections

s1$group = 0
group = 0
for(i in seq(dim(s1)[1])){
	if(all(s1$start[i] != s1$end[i-1]+1)){
		group = group+1
	}
	s1[i,'group'] <- group
}

# Now classify bins based the local variance of the scaled_r score

Wiggle <- NULL
for(z in seq(group)){
	scale_r <- s1[s1$group == z,'archi.cor.scale'] 
	zero <- rep(0, length(scale_r))
	wiggle <- data.frame(scale_r_score = zero, mean = zero , bottom = zero , top = zero)
	for(i in seq(along=scale_r)){	
		if((i == length(scale_r)) & (i == 1)){
			a <- mean(scale_r[i])
			s <- sd(s1$archi.cor.scale)
			n <- 1
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}else if(i == 1){
			a <- mean(scale_r[i:(i+1)])
			s <- sd(scale_r[i:(i+1)])
			n <- 2
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}else if(i == length(scale_r)){
			a <- mean(scale_r[(i-1):i])
			s <- sd(scale_r[(i-1):i])
			n <- 2
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}else if(i == 2){
			a <- mean(scale_r[(i-1):(i+1)])
			s <- sd(scale_r[(i-1):(i+1)])
			n <- 3
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}else if(i == length(scale_r) - 1){
			a <- mean(scale_r[(i-1):(i+1)])
			s <- sd(scale_r[(i-1):(i+1)])
			n <- 3
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}else{
			a <- mean(scale_r[(i-2):(i+2)])
			s <- sd(scale_r[(i-2):(i+2)])
			n <- 5
			error <- qnorm(0.975)*s/sqrt(n)
			left <- a-error
			right <- a+error
		}
		wiggle$scale_r_score[i] <- scale_r[i]
		wiggle$mean[i] <- a 
		wiggle$bottom[i] <- left		
		wiggle$top[i] <- right
	}
	Wiggle <- rbind(Wiggle,wiggle )
}


n1 <- 700
n2 <- 800
plot(Wiggle$mean[n1:n2], type = "l")
lines(Wiggle$top[n1:n2], col = 2)
lines(Wiggle$bottom[n1:n2], col = 2)
lines(rep(0,length(n1:n2)), col = 3)

Wiggle$r_class <- rep("M", dim(Wiggle)[1])
Wiggle$r_class[Wiggle$bottom > 0] <- "H"
Wiggle$r_class[Wiggle$top < 0] <- "L"




# final output will need to be mantal test and partial mantal test results

# a table describing each bins genomic content as well as the r_class for each bin 





require(ggplot2)

s1.pca <- prcomp(s1[5:(length(s1) - 2)], scale.=TRUE)
s1.re.pca <- prcomp(replot_s1_s2[2:(length(replot_s1_s2))], scale.=TRUE)
s2.pca <- prcomp(s2[5:(length(s2))], scale.=TRUE)



#qplot(s1.pca$x[,1], s1.pca$x[,2], colour = s1$archi.cor.scale) + scale_colour_gradient(limits=c(min(s1$archi.cor.scale), max(s1$archi.cor.scale)), high = "red", low = "blue")
#qplot(s1.re.pca$x[,1], s1.re.pca$x[,2], colour = s1$archi.cor.scale) + scale_colour_gradient(limits=c(min(s1$archi.cor.scale), max(s1$archi.cor.scale)), high = "red", low = "blue")


# how can i get the colours to go from 1 to -1


colours <- rep("1", dim(Wiggle)[1])
colours[Wiggle$r_class == "L"] <- 2
colours[Wiggle$r_class == "H"] <- 3

pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/", spec1, "_in_", spec2,"_analysis","_scaled_" , SCALE, sep = ""), onefile = TRUE)
plot(s1.re.pca$x[,1], s1.re.pca$x[,2], col = colours, pch = 16, cex = .5, main = paste(spec1 ," in ", spec2," pairwise distance cor", sep = ""))
biplot(s1.re.pca, xlabs=rep(".", dim(s1.re.pca$x)[1]), cex = c(2,1), main = paste(spec1 ," in ", spec2," biplot", sep = ""))
biplot(s2.pca, xlabs=rep(".", dim(s2.pca$x)[1]), cex = c(2,1), main = paste(spec2 ," biplot", sep = ""))
plot(s1.pca$x[,1], s1.pca$x[,2], col = colours, pch = 16, cex = .5, main = paste(spec1 ," pairwise distance cor", sep = ""))
biplot(s1.pca, xlabs=rep(".", dim(s1.pca$x)[1]), cex = c(2,1), main = paste(spec1 ," biplot", sep = ""))
heatmap.2(cor(s1[,5:(length(s1) - 2)], replot_s1_s2[,2:length(replot_s1_s2)]), trace = "none",margins = c(8,8), xlab = spec2, ylab = spec1, main = "TE corelations between species", col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
heatmap.2(cor(s1[,5:(length(s1) - 2)]), trace = "none",margins = c(8,8), main = paste("TE corelations in",spec1,sep=" "), col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
heatmap.2(cor(s2[,5:(length(s2))]), trace = "none",margins = c(8,8), main = paste("TE corelations in",spec2,sep=" "), col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
dev.off()


plot(s1.pca$x[colours > 1,1], s1.pca$x[colours > 1,2], col = colours[colours > 1], pch = 16, cex = .5, main = paste(spec1 ," pairwise distance cor", sep = ""))
plot(s1.re.pca$x[colours > 1,1], s1.re.pca$x[colours > 1,2], col = colours[colours > 1], pch = 16, cex = .5, main = paste(spec1 ," in ", spec2," pairwise distance cor", sep = ""))

# there is the posibility to get this going with mouse aswell 
# the things I need are the mouse phastcons, mouse exons in gtf format like I have with human and the mouse chrom info


# lets look at a 3D plot that also has conservation of conserved elements 
# there are three kinds of conserved element types we want to look at
# 1. exonic
# 2. non exonix
# 3. all 



con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
sLengths <- SEQ$size
names(sLengths) <- SEQ[,'chr'] 


library(GenomicRanges)
# exon informatio
gene  <- read.delim("~/Desktop/Repeat evolution/bin_synteny/genome_anotation/hg19_refgene.gtf", header = FALSE)
list <- strsplit(as.character(gene[,9]), " ")
library("plyr")
df <- ldply(list)
gene[,9] <- df[,2]
gene[,9] <- gsub(";", "", gene[,9])
gene[,10] <- gene[,5] - gene[,4] +1
CDs <- gene[gene[,3] == "exon",]

#human region
# - 1 on the end of a range so they don't merge
H.gr <- GRanges(seqnames = Rle(s1$chr), ranges = IRanges(start = s1$start, end = s1$end -1) ,ID = s1$binID)

# conserved elements
CE <- read.table("~/Desktop/Repeat evolution/Genome_region_conservation/phastConsElements46wayPlacental.txt")
CE <- CE[-(grep("_", CE[,2])),]
# remove the extra chromosomes because they are not used in downstream analysis
CE.gr <- reduce(GRanges(seqnames = Rle(CE[,2]), ranges = IRanges(start = CE[,3], end = CE[,4]),seqlengths = sLengths ))
exon <- gene[gene[,3] == "exon",]
exon.gr <- reduce(GRanges(seqnames = Rle(exon[,1]), ranges = IRanges(start = exon[,4], end = exon[,5]), seqlengths = sLengths))
nonexon.gr <- gaps(exon.gr)
nonexon.gr <- nonexon.gr[strand(nonexon.gr) == "*"]





# proportion of exonic bases that are conserved elements per genomic bin
	# for each bin(the amount of CE exonic bases / the amount of exonic bases)
exon.bin.gr <- intersect(H.gr, exon.gr)
exon.bin.CE.gr <- intersect(exon.bin.gr, CE.gr)
	# we have exons and exonCEs, sort them into their bins
bin.sort.exon <- as.matrix(findOverlaps(H.gr, exon.bin.gr))
bin.sort.exon <- data.frame(bin.ID = H.gr$ID[bin.sort.exon[,1]], exon_length = width(exon.bin.gr[bin.sort.exon[,2]]))
bin.sort.exon.CE <- as.matrix(findOverlaps(H.gr, exon.bin.CE.gr))
bin.sort.exon.CE <- data.frame(bin.ID = H.gr$ID[bin.sort.exon.CE[,1]], CE_exon_length = width(exon.bin.CE.gr[bin.sort.exon.CE[,2]]))

bin.ID <- H.gr$ID
merged.exon.CE <- data.frame(bin.ID = bin.ID, exon_length = rep(0, length(bin.ID)), CE_exon_length = rep(0, length(bin.ID)))
for(i in 1:length(bin.ID)){
	merged.exon.CE$exon_length[i] <- sum(bin.sort.exon[bin.sort.exon$bin.ID == bin.ID[i] ,'exon_length'])
	merged.exon.CE$CE_exon_length[i] <- sum(bin.sort.exon.CE[bin.sort.exon.CE$bin.ID == bin.ID[i] ,'CE_exon_length'])
}
# we have to remove bins that have 0s
merged.exon.CE$P_exonic_CE <- rep(0, dim(merged.exon.CE)[1])
# if there are 0 exonic bases that bin gets a Pscore of 0
merged.exon.CE$P_exonic_CE[merged.exon.CE[,2] > 0] <- merged.exon.CE[merged.exon.CE[,2] > 0,]$CE_exon_length / merged.exon.CE[merged.exon.CE[,2] > 0,]$exon_length

#boxplot(merged.exon.CE$P_exonic_CE[colours == 2], merged.exon.CE$P_exonic_CE[colours == 1], merged.exon.CE$P_exonic_CE[colours == 3])


# proportion of non exonic bases that are conserved elements per genomic bin

nonexon.bin.gr <- intersect(H.gr, nonexon.gr)
nonexon.bin.CE.gr <- intersect(nonexon.bin.gr, CE.gr)
	# we have nonexons and nonexonCEs, sort them into their bins
bin.sort.nonexon <- as.matrix(findOverlaps(H.gr, nonexon.bin.gr))
bin.sort.nonexon <- data.frame(bin.ID = H.gr$ID[bin.sort.nonexon[,1]], nonexon_length = width(nonexon.bin.gr[bin.sort.nonexon[,2]]))
bin.sort.nonexon.CE <- as.matrix(findOverlaps(H.gr, nonexon.bin.CE.gr))
bin.sort.nonexon.CE <- data.frame(bin.ID = H.gr$ID[bin.sort.nonexon.CE[,1]], CE_nonexon_length = width(nonexon.bin.CE.gr[bin.sort.nonexon.CE[,2]]))

bin.ID <- H.gr$ID
merged.nonexon.CE <- data.frame(bin.ID = bin.ID, nonexon_length = rep(0, length(bin.ID)), CE_nonexon_length = rep(0, length(bin.ID)))
for(i in 1:length(bin.ID)){
	merged.nonexon.CE$nonexon_length[i] <- sum(bin.sort.nonexon[bin.sort.nonexon$bin.ID == bin.ID[i] ,'nonexon_length'])
	merged.nonexon.CE$CE_nonexon_length[i] <- sum(bin.sort.nonexon.CE[bin.sort.nonexon.CE$bin.ID == bin.ID[i] ,'CE_nonexon_length'])
}
# we have to remove bins that have 0s
merged.nonexon.CE$P_nonexonic_CE <- rep(0, dim(merged.nonexon.CE)[1])
# if there are 0 nonexonic bases that bin gets a Pscore of 0
merged.nonexon.CE$P_nonexonic_CE[merged.nonexon.CE[,2] > 0] <- merged.nonexon.CE[merged.nonexon.CE[,2] > 0,]$CE_nonexon_length / merged.nonexon.CE[merged.nonexon.CE[,2] > 0,]$nonexon_length

#boxplot(merged.nonexon.CE$P_nonexonic_CE[colours == 2], merged.nonexon.CE$P_nonexonic_CE[colours == 1], merged.nonexon.CE$P_nonexonic_CE[colours == 3])


# lets make it 3D

#library(rgl)
#plot3d(s1.pca$x[,1], s1.pca$x[,2], merged.nonexon.CE$P_nonexonic_CE, col = colours, size = 5)
#plot3d(s1.pca$x[,1], s1.pca$x[,2], merged.exon.CE$P_exonic_CE, col = colours, size = 5)



# proportion of bases that are conserved elements per genomic bin


# do an intersect then find overlaps with CE.gr ?
# the intersect will break anything that needs to be broken

CE.inter.gr <- intersect(H.gr,CE.gr )
CE.bin.sort <- as.matrix(findOverlaps(H.gr, CE.inter.gr))
CE.bin.sort <- data.frame(bin.ID = H.gr$ID[CE.bin.sort[,1]], CE_length = width(CE.inter.gr[CE.bin.sort[,2]]))
merged.CE <- data.frame(bin.ID = bin.ID, bin_length = rep(0, length(bin.ID)), CE_length = rep(0, length(bin.ID)))
for(i in 1:length(bin.ID)){
	merged.CE$bin_length[i] <- KnownS1[KnownS1[,1] == bin.ID[i], 2]
	merged.CE$CE_length[i] <- sum(CE.bin.sort[CE.bin.sort$bin.ID == bin.ID[i] ,'CE_length'])
	}

merged.CE$P_CE <- merged.CE$CE_length/merged.CE$bin_length

#boxplot(merged.CE$P_CE[colours == 2],merged.CE$P_CE[colours == 1], merged.CE$P_CE[colours == 3])


#library(rgl)
#plot3d(s1.pca$x[,1], s1.pca$x[,2], merged.nonexon.CE$P_nonexonic_CE, col = colours, size = 5)
#plot3d(s1.pca$x[,1], s1.pca$x[,2], merged.exon.CE$P_exonic_CE, col = colours, size = 5)
#plot3d(s1.pca$x[,1], s1.pca$x[,2], merged.CE$P_CE, col = colours, size = 5)




# introduce the S1 pca graph used to initally get the dimensions



# we have the known and the S1 pca with the same amount of bins, make sure they are in the correct order and then we can r

Hclass <- read.table("~/Desktop/Dist_matrix_TE_div/RARateClassifiedBins.5_1.5.txt", header = TRUE)
Hclass$binID <- KnownS1[,1]
Hclass <- Hclass[Hclass$binID %in% s1$binID,]



colours2 <- rep(1, dim(Hclass)[1])
colours2[Hclass$Classification == "L"] <- 2
colours2[Hclass$Classification == "H"] <- 3


# make a boxplot of each combination

HH.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "H"]
HM.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "M"]
HL.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "L"]
MH.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "H"]
MM.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "M"]
ML.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "L"]
LH.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "H"]
LM.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "M"]
LL.nonexon <- merged.nonexon.CE$P_nonexonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "L"]

# make a boxplot of each combination but exonic


HH.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "H"]
HM.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "M"]
HL.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "H" & Hclass$Classification == "L"]
MH.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "H"]
MM.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "M"]
ML.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "M" & Hclass$Classification == "L"]
LH.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "H"]
LM.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "M"]
LL.exon <- merged.exon.CE$P_exonic_CE[Wiggle$r_class == "L" & Hclass$Classification == "L"]




# boxplot for total proportion of CE

HH.CE <- merged.CE$P_CE[Wiggle$r_class == "H" & Hclass$Classification == "H"]
HM.CE <- merged.CE$P_CE[Wiggle$r_class == "H" & Hclass$Classification == "M"]
HL.CE <- merged.CE$P_CE[Wiggle$r_class == "H" & Hclass$Classification == "L"]
MH.CE<- merged.CE$P_CE[Wiggle$r_class == "M" & Hclass$Classification == "H"]
MM.CE <- merged.CE$P_CE[Wiggle$r_class == "M" & Hclass$Classification == "M"]
ML.CE <- merged.CE$P_CE[Wiggle$r_class == "M" & Hclass$Classification == "L"]
LH.CE <- merged.CE$P_CE[Wiggle$r_class == "L" & Hclass$Classification == "H"]
LM.CE <- merged.CE$P_CE[Wiggle$r_class == "L" & Hclass$Classification == "M"]
LL.CE <- merged.CE$P_CE[Wiggle$r_class == "L" & Hclass$Classification == "L"]



# look at how the boxes change when you have the PC2 values
# also look at organiztion scores
# look at the boxplots for the HLM of both organization group and PC2 score group
# axis taken care of need to think of good titles fro graphs



group.names1 = c("H", "M", "L")
group.names2 = c("HH","HM", "HL", "MH", "MM", "ML", "LH", "LM", "LL")
group.names3 = c("HH", "HL", "LH", "LL")
x.name2 <- "organization group"
x.name1 <- "PC2 group"
x.name3 <- "organization group / PC2 group"
y.nameCE <- "proportion of base pairs overlapped by conserved elements in each bin"
y.nameExon <- "proportion of exonic base pairs overlapped by conserved elements in each bin"
y.nameNonExon <- "proportion of nonexonic base pairs overlapped by conserved elements in each bin"
main2 = paste("Negative selection in", spec1, "in regions of different \norganizational conservation between", spec1, "and", spec2)
main1 = paste("Negative selection in", spec1, "in regions of different ancestry")
main3 = paste("Negative selection in", spec1, "in regions of different \norganizational conservation between", spec1, "and", spec2, "and \nregions of different ancestry in" ,spec1)

pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/", spec1, "_and_", spec2,"_boxplot_CE_analysis", sep = ""), onefile = TRUE)

boxplot(merged.CE$P_CE[colours2 == 3], merged.CE$P_CE[colours2 == 1], merged.CE$P_CE[colours2 == 2], names = group.names1, xlab = x.name1, ylab = y.nameCE, main = main1, notch = TRUE)
boxplot(merged.exon.CE$P_exonic_CE[colours2 == 3],merged.exon.CE$P_exonic_CE[colours2 == 1], merged.exon.CE$P_exonic_CE[colours2 == 2], names = group.names1, xlab = x.name1, ylab = y.nameExon, main = main1, notch = TRUE)
boxplot(merged.nonexon.CE$P_nonexonic_CE[colours2 == 3],merged.nonexon.CE$P_nonexonic_CE[colours2 == 1], merged.nonexon.CE$P_nonexonic_CE[colours2 == 2], names = group.names1, xlab = x.name1, ylab = y.nameNonExon, main = main1, notch = TRUE)

boxplot(merged.CE$P_CE[colours == 3],merged.CE$P_CE[colours == 1], merged.CE$P_CE[colours == 2], names = group.names1, xlab = x.name2, ylab = y.nameCE, main = main2, notch = TRUE)
boxplot(merged.exon.CE$P_exonic_CE[colours == 3], merged.exon.CE$P_exonic_CE[colours == 1], merged.exon.CE$P_exonic_CE[colours == 2], names = group.names1, xlab = x.name2, ylab = y.nameExon, main = main2, notch = TRUE)
boxplot(merged.nonexon.CE$P_nonexonic_CE[colours == 3], merged.nonexon.CE$P_nonexonic_CE[colours == 1], merged.nonexon.CE$P_nonexonic_CE[colours == 2], names = group.names1, xlab = x.name2, ylab = y.nameNonExon, main = main2, notch = TRUE)

boxplot(HH.CE,HM.CE,HL.CE,MH.CE,MM.CE,ML.CE,LH.CE,LM.CE,LL.CE, names = group.names2, xlab = x.name3, ylab = y.nameCE, main = main3, notch = TRUE)
boxplot(HH.exon,HM.exon,HL.exon,MH.exon,MM.exon,ML.exon,LH.exon,LM.exon,LL.exon, names = group.names2, xlab = x.name3, ylab = y.nameExon, main = main3, notch = TRUE)
boxplot(HH.nonexon,HM.nonexon,HL.nonexon,MH.nonexon,MM.nonexon,ML.nonexon,LH.nonexon,LM.nonexon,LL.nonexon, names = group.names2, xlab = x.name3, ylab = y.nameNonExon, main = main3, notch = TRUE)

boxplot(HH.CE,HL.CE,LH.CE,LL.CE, names = group.names3, xlab = x.name3, ylab = y.nameCE, main = main3, notch = TRUE)
boxplot(HH.exon,HL.exon,LH.exon,LL.exon, names = group.names3, xlab = x.name3, ylab = y.nameExon, main = main3, notch = TRUE)
boxplot(HH.nonexon,HL.nonexon,LH.nonexon,LL.nonexon, names = group.names3, xlab = x.name3, ylab = y.nameNonExon, main = main3, notch = TRUE)

dev.off()


# write the results into a table so they can later be compared with each other. 
# the way to do this is 


# overlap analysis in human



# we want to look at the intersect of each region and potentially the union. 


# this will mean writing tables to do so
#Basicly write out tables for further analysis

#merge each of the conservation tables that measures the per bin conservation

M.bin.CE.class <- merge(merged.CE, merged.exon.CE, by.x = 1, by.y = 1)
M.bin.CE.class <- merge(M.bin.CE.class,merged.nonexon.CE, by.x = 1, by.y = 1 )
M.bin.CE.class <- merge(M.bin.CE.class, Hclass, by.x= 1, by.y = 6)
M.bin.CE.class$r_class <- Wiggle$r_class

write.table(M.bin.CE.class, file=paste("/Users/labadmin/Desktop/Dist_matrix_TE_div/merged.tables/", spec1, "_",spec2,"merged.CE.class", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE )





# The next thing to start looking at is GO terms and things to do with expression
# what does it mean for their function  to be conserved
# with GO terms we are lookikng to see if anything associates with something that is known to be different between two species. 
# What does it mean for genes to be conserved and what does it mean for their expression patterns to be conserved
# do the different pairwise comparisons overlap

# For Go term analysis we identify genes and get their gene names




















library("vegan")


mant <- mantel(S1_dist, S1_s2_dist)
mant.result <- c(Species_1 = spec1, Species_2 = spec2, Method = mant$method, Statisitic = mant$statistic, Signif= mant$signif)




# these can be stuck together if i just rename the TE names according to their species.
# Then ill be able to run this script in single run
# and Ill get various pairwise comparisons

results_s1_s2_partial_mantal_S2_TE <- data.frame(TE.name = rep(0, length(replot_s1_s2)-1), statistic = rep(0, length(replot_s1_s2)-1), signif = rep(0, length(replot_s1_s2)-1))
for(i in seq(along = S2_TE_name)){
	m.part <- mantel.partial(S1_dist,S1_s2_dist,dist(replot_s1_s2[,TE_name[i]]))
	results_s1_s2_partial_mantal_S2_TE[i,'TE.name'] = as.character(TE_name[i])
	results_s1_s2_partial_mantal_S2_TE[i,'statistic'] = m.part$statistic
	results_s1_s2_partial_mantal_S2_TE[i,'signif'] = m.part$signif
}

S1TE_name <- colnames(s1)[5:length(s1)]
results_s1_s2_partial_mantal_S1_TE <- data.frame(S1TE.name = rep(0, length(S1TE_name)), statistic = rep(0, length(S1TE_name)), signif = rep(0, length(S1TE_name)))


for(i in seq(along = S1TE_name)){
	m.part <- mantel.partial(S1_dist,S1_s2_dist,dist(s1[,S1TE_name[i]]))
	results_s1_s2_partial_mantal_S1_TE[i,'TE.name'] = as.character(S1TE_name[i])
	results_s1_s2_partial_mantal_S1_TE[i,'statistic'] = m.part$statistic
	results_s1_s2_partial_mantal_S1_TE[i,'signif'] = m.part$signif
}





write.table(mant.result,file= "../mant.txt")






m.cor <- mantel.correlog(S1_dist, S1_s2_dist)	
plot(m.cor)	





