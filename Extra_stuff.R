

# mantel stuff and CE stuff 

library("vegan")




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



###########
#
#    Mantel TEST
#
###########

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














#Basiclcly CE stuff we need to do only one for the root species







#####
#
#   Chrom info
#
####

con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
sLengths <- SEQ$size
names(sLengths) <- SEQ[,'chr'] 


library(GenomicRanges)


#########
#
# exon informatio
#
###########

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

#####
#
#   Conserved elements
#
####

CE <- read.table("~/Desktop/Repeat evolution/Genome_region_conservation/phastConsElements46wayPlacental.txt")




CE <- CE[-(grep("_", CE[,2])),]
# remove the extra chromosomes because they are not used in downstream analysis
CE.gr <- reduce(GRanges(seqnames = Rle(CE[,2]), ranges = IRanges(start = CE[,3], end = CE[,4]),seqlengths = sLengths ))
exon <- gene[gene[,3] == "exon",]
exon.gr <- reduce(GRanges(seqnames = Rle(exon[,1]), ranges = IRanges(start = exon[,4], end = exon[,5]), seqlengths = sLengths))





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



# proportion of non exonic bases that are conserved elements per genomic bin
# to get the non-exon rates we will take the exon rates from the overall rate

# A = known bases - exon length per bin
# B = know bases length with CE - exon length with CE
# non exon rate = B/A





#merge each of the conservation tables that measures the per bin conservation

M.bin.CE.class <- merge(merged.CE, merged.exon.CE, by.x = 1, by.y = 1)
M.bin.CE.class <- merge(M.bin.CE.class,merged.nonexon.CE, by.x = 1, by.y = 1 )
M.bin.CE.class$r_class <- Wiggle$r_class


# IT would be nice to add the TE content of each bin too


write.table(M.bin.CE.class, file=paste("/Users/labadmin/Desktop/Dist_matrix_TE_div/merged.tables/", spec1, "_",spec2,"merged.CE.class", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE )





