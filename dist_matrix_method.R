# get the S bins and calculate the distance matricies

# all options need to be at the top

# this script just does S1 to S2


# how we get our pairwise information
require(rtracklayer)
require(gplots)
require(GenomicRanges)
rm(list= ls())
setwd("~/Desktop/Dist_matrix_TE_div/")



#####
# Load in Tables and process

spec1 <- "Human"
if(spec1 == "Mouse"){UCSCspec = "mm10"}
if(spec1 == "Human"){UCSCspec = "hg19"}
if(spec1 == "Opossum"){UCSCspec = "monDom5"}
all.species <- c("Human", "Mouse", "Rat", "Dog", "Bovine", "Horse", "Opossum", "Platypus")
all.species <- all.species[!all.species == spec1]
all_sbin <- list.files(path = paste("S_bin/", spec1, "_sbin/", sep = ""))
useful.species <- NULL
for(x in 1:length(all.species)){
	# do a check up here to update the list and then download the file that are needed
	update <- NULL
	if(length(all_sbin[grep(all.species[x], all_sbin)]) == 1){
		update <- all.species[x]
	}
	useful.species <- c(useful.species, update)
}


R_classes <- NULL
for(x in 1:length(useful.species)){
spec2 <- useful.species[x]


# load in Sbin	
S.bin.name <- paste("S_bin/", spec1, "_sbin/", spec1, "_aligning_", spec2, "_select", sep = "")
S.bin_s1_s2 <- read.table(S.bin.name,header = TRUE)

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
SCALE = "no"

#create objects into whcih i will store the binsizes
s1name <- paste("count_tables/",spec1, "_AllBinCounts.txt", sep = "")
s2name <- paste("count_tables/",spec2, "_AllBinCounts.txt", sep = "")
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
Wiggle$r_class[Wiggle$bottom > 0] <- "L"
Wiggle$r_class[Wiggle$top < 0] <- "H"




s1.pca <- prcomp(s1[5:(length(s1) - 2)], scale.=TRUE)
s1.re.pca <- prcomp(replot_s1_s2[2:(length(replot_s1_s2))], scale.=TRUE)
s2.pca <- prcomp(s2[5:(length(s2))], scale.=TRUE)



colours <- rep("1", dim(Wiggle)[1])
colours[Wiggle$r_class == "H"] <- 2
colours[Wiggle$r_class == "L"] <- 3

#pdf(file = paste("plot_results/", spec1, "_in_", spec2,"_NGenes.", keep.NGenes ,"_NG4s.", keep.NG4s, "_NCpGI.", keep.NCpGI, "_CpGBP.", keep.CpGBP, "_GC.", keep.GC, "_SCALE.", SCALE,".pdf" ,sep = ""), onefile = TRUE)
#plot(s1.re.pca$x[,1], s1.re.pca$x[,2], col = colours, pch = 16, cex = .5, main = paste(spec1 ," in ", spec2," pairwise distance cor", sep = ""))
#biplot(s1.re.pca, xlabs=rep(".", dim(s1.re.pca$x)[1]), cex = c(2,1), main = paste(spec1 ," in ", spec2," biplot", sep = ""))
#biplot(s2.pca, xlabs=rep(".", dim(s2.pca$x)[1]), cex = c(2,1), main = paste(spec2 ," biplot", sep = ""))
#plot(s1.pca$x[,1], s1.pca$x[,2], col = colours, pch = 16, cex = .5, main = paste(spec1 ," pairwise distance cor", sep = ""))
#biplot(s1.pca, xlabs=rep(".", dim(s1.pca$x)[1]), cex = c(2,1), main = paste(spec1 ," biplot", sep = ""))
#heatmap.2(cor(s1[,5:(length(s1) - 2)], replot_s1_s2[,2:length(replot_s1_s2)]), trace = "none",margins = c(8,8), xlab = spec2, ylab = spec1, main = "TE corelations between species", col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
#heatmap.2(cor(s1[,5:(length(s1) - 2)]), trace = "none",margins = c(8,8), main = paste("TE corelations in",spec1,sep=" "), col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
#heatmap.2(cor(s2[,5:(length(s2))]), trace = "none",margins = c(8,8), main = paste("TE corelations in",spec2,sep=" "), col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)
#dev.off()

R_info <- data.frame(s1$binID,Wiggle$r_class, cor.diag)
names(R_info) <- c("s1.binID", paste(useful.species[x], "r_class", sep ="_"), paste(useful.species[x], "dist_cor", sep ="_"))
R_classes <- c(R_classes, list(R_info))

# keep it as a list and save as an R opject

}
names(R_classes) <- useful.species


#########
#
#      new section
#
######



# build a dataframe that has NA and the values associated with each bin

# merge on the all.x = TRUE, thereby extra rows will be added with an NA 

#create objects into whcih i will store the binsizes
s1 <- read.table(s1name, header = TRUE)
slist <- list(s1)


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

s1 <- slist[[1]][,!(colnames(slist[[1]]) == "Known")]

if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}

all.spec.class <- R_classes[[1]]
for( i in 1:(length(R_classes)-1)){
	all.spec.class <- merge(all.spec.class, R_classes[i+1], by.x = 1, by.y =1, all= TRUE)
	colnames(all.spec.class) <- gsub(paste(names(R_classes)[i+1], ".",sep = ""), "", colnames(all.spec.class), fixed = TRUE)
}


#A <- all.spec.class[,grep("_dist_cor", colnames(all.spec.class))]
#B <- as.data.frame(A)
#boxplot(B, notch=TRUE)

#C <- B[complete.cases(B[1:4]),1:4]
#plot(hclust(dist((C))))
#heatmap((as.matrix(C)), scale="none")

# probably just merge it with s1 then i can wtite it up and do more analysis

done <- merge(s1, all.spec.class, by.x = 2, by.y = 1)

final_name <- paste(spec1, "_R.cor.cons", "_NGenes.", keep.NGenes ,"_NG4s.", keep.NG4s, "_NCpGI.", keep.NCpGI, "_CpGBP.", keep.CpGBP, "_GC.", keep.GC, "_SCALE.", SCALE, sep="")

write.table(done, file=paste("r_class_table/", final_name,"_no_conservation.txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE )




#########
#
# exon informatio
#
###########

mySession <- browserSession("UCSC")
genome(mySession) <- UCSCspec

track.gene <- track(ucscTableQuery(mySession,track = "refGene", table= "refGene"), asRangedData = FALSE)
GFF <- asGFF(blocks(track.gene), parentType = "mRNA", childType = "exon")
exon <- GFF[elementMetadata(GFF)$type == "exon"]
strand(exon) <- "*"
exon.gr <- reduce(exon)






#human region
# - 1 on the end of a range so they don't merge
S1.gr <- GRanges(seqnames = Rle(s1$chr), ranges = IRanges(start = s1$start, end = s1$end -1) ,ID = s1$binID)



# out here we can get all the exon lengths per bin 


# calculate the binwise exon length

exon.bin.gr <- intersect(S1.gr, exon.gr)
bin.sort.exon <- as.matrix(findOverlaps(S1.gr, exon.bin.gr))
bin.sort.exon <- data.frame(bin.ID = S1.gr$ID[bin.sort.exon[,1]], exon_length = width(exon.bin.gr[bin.sort.exon[,2]]))
bin.ID <- S1.gr$ID
exon_length <-  data.frame(bin.ID = bin.ID, exon_length = rep(0, length(bin.ID)))
for(i in 1:length(bin.ID)){
	exon_length$exon_length[i] <- sum(bin.sort.exon[bin.sort.exon$bin.ID == bin.ID[i] ,'exon_length'])
	}


# calculate the non_binwise exon length and total bin length

# estimate because we assume that none of the exon ranges span any unkown bases
non_exon_length_estimate <- merge(KnownS1, exon_length, by.x = 1, by.y = 1)
bin_length <- data.frame(bin.ID = non_exon_length_estimate[,1], bin_length = non_exon_length_estimate[,2])
non_exon_length_estimate <- data.frame(bin.ID = non_exon_length_estimate[,1], non_exon_length_estimate = (non_exon_length_estimate[,2] - non_exon_length_estimate[,3]))

lengths <- data.frame(bin.ID = bin.ID, exon_length = exon_length$exon_length, bin_length = bin_length$bin_length, non_exon_length_estimate = non_exon_length_estimate$non_exon_length_estimate)
# conserved elements



#####
#
#   Conserved elements
#
####



if(spec1 == "Human"){
	spec.track <- "cons46way"
	spec.table <- c("phastConsElements46wayPrimates", "phastConsElements46wayPlacental", "phastConsElements46way")
	spec.names <- c("Primate_", "Placental_", "Vertebrate_")
	}

if(spec1 == "Mouse"){
	spec.track <- "cons60way"
	spec.table <- c("phastConsElements60wayEuarchontoGlires", "phastConsElements60wayPlacental", "phastConsElements60way")
	spec.names <- c("EuarchontoGlires_", "Placental_", "Vertebrate_")
	}

if(spec1 == "Opossum"){
	spec.track <- "phastConsElements9way"
	spec.table <- "phastConsElements9way"
	spec.names <- "Vertebrate_"
	}


s1.con <- s1
for(x in 1:length(spec.table)){
	
	track.CE <- track(ucscTableQuery(mySession,track = spec.track, table= spec.table[x]), asRangedData = FALSE)
	CE.gr <- reduce(track.CE)
	exon.bin.CE.gr <- intersect(exon.bin.gr, CE.gr)
	bin.sort.exon.CE <- as.matrix(findOverlaps(S1.gr, exon.bin.CE.gr))
	bin.sort.exon.CE <- data.frame(bin.ID = S1.gr$ID[bin.sort.exon.CE[,1]], CE_exon_length = width(exon.bin.CE.gr[bin.sort.exon.CE[,2]]))
	total.bin.CE.gr <- intersect(S1.gr, CE.gr)
	CE.bin.sort <- as.matrix(findOverlaps(S1.gr, total.bin.CE.gr))
	CE.bin.sort <- data.frame(bin.ID = S1.gr$ID[CE.bin.sort[,1]], CE_length = width(total.bin.CE.gr[CE.bin.sort[,2]]))
	merged_CE_length <- data.frame(bin.ID = bin.ID, CE_exon_length = rep(0, length(bin.ID)), CE_length = rep(0, length(bin.ID)))
	for(i in 1:length(bin.ID)){
		merged_CE_length$CE_exon_length[i] <- sum(bin.sort.exon.CE[bin.sort.exon.CE$bin.ID == bin.ID[i] ,'CE_exon_length'])
		merged_CE_length$CE_length[i] <- sum(CE.bin.sort[CE.bin.sort$bin.ID == bin.ID[i] ,'CE_length'])
	}
	
	merged_CE_length$CE_nonexon_length <- merged_CE_length$CE_length - merged_CE_length$CE_exon_length
	
	merged_CE_rate <- data.frame(bin.ID = bin.ID, CE_exon_rate = rep(NA, length(bin.ID)), CE_rate = rep(NA, length(bin.ID)), CE_nonexon_rate = rep(NA, length(bin.ID)))
	merged_CE_rate$CE_exon_rate[lengths$exon_length > 0] <- merged_CE_length$CE_exon_length[lengths$exon_length > 0]/lengths$exon_length[lengths$exon_length > 0]
	merged_CE_rate$CE_nonexon_rate[lengths$non_exon_length_estimate > 0] <- merged_CE_length$CE_nonexon_length[lengths$non_exon_length_estimate > 0] / lengths$non_exon_length_estimate[lengths$non_exon_length_estimate > 0]
	merged_CE_rate$CE_rate[lengths$bin_length > 0] <- merged_CE_length$CE_length[lengths$bin_length > 0] / lengths$bin_length[lengths$bin_length > 0]
	
	
	if(all(lengths$bin.ID == merged_CE_length$bin.ID) & all(lengths$bin.ID == merged_CE_rate$bin.ID) & all(lengths$bin.ID == s1$binID)){
		all.CE <- cbind(lengths[2:length(lengths)], merged_CE_length[,2:length(merged_CE_length)], merged_CE_rate[,2:length(merged_CE_rate)])
		colnames(all.CE) <- paste(spec.names[x], colnames(all.CE), sep = "")
		s1.con <- cbind(s1.con, all.CE)
	}
}

final_data <- merge(s1.con, all.spec.class, by.x= 2, by.y= 1, all.x = TRUE)




# remember to merge the multispecies comparisons here with the rest of the stuff and wirte up the table for futher analysis


final_name <- paste(spec1, "_R.cor.cons", "_NGenes.", keep.NGenes ,"_NG4s.", keep.NG4s, "_NCpGI.", keep.NCpGI, "_CpGBP.", keep.CpGBP, "_GC.", keep.GC, "_SCALE.", SCALE, sep="")

write.table(final_data, file=paste("r_class_table/", final_name,".txt", sep = ""), quote = FALSE, row.names = FALSE, col.names = TRUE )
	
	
	
	
	
	
	
	
	



