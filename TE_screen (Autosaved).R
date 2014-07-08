

# proportion of variance at a PC

# corealtio nbetween sspecies at that PC


# what would a negative corelation mean 
# so we are using PC2 to divide species 


# waht is meant by us recognizing groups based on a metric and then being surprised when groups differ on that metric.


# Low PC2 and low org in diff species appear to align with each other to a greater degree than what could be expected by chance. 

# is this result real or the product of something else. 


rm(list = ls())
require(GenomicRanges)
require(ggbio)
require(gplots)

setwd("~/Desktop/Dist_matrix_TE_div")

spec1 <- "Mouse"
if(spec1 == "Mouse"){UCSCspec = "mm10"}
if(spec1 == "Human"){UCSCspec = "hg19"}
if(spec1 == "Opossum"){UCSCspec = "monDom5"}


# file selection 
keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "yes"

file_name <- paste(spec1, "_R.cor.cons", "_NGenes.", keep.NGenes ,"_NG4s.", keep.NG4s, "_NCpGI.", keep.NCpGI, "_CpGBP.", keep.CpGBP, "_GC.", keep.GC, "_SCALE.", SCALE,".txt" ,sep="")

# read in the table and seperate out conservation, TE dist, and species comparison
# make sure each one segregate with the binID

# maybe get a list of usefule species too 
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

useful.species1 <- useful.species[!(useful.species == "Platypus")]



all.species1 <- read.table(paste("r_class_table/", file_name, sep = ""), header = TRUE)

# mobility information
s1.mob <- all.species1[,c("binID",useful.species)] 

#conservation information
s1.con <- all.species1[,c(1,grep("rate", colnames(all.species1)) , grep("length", colnames(all.species1)))]

# the TE coordinates and bin information
s1 <- all.species1[,-(c(grep("rate", colnames(all.species1)) , grep("length", colnames(all.species1))))]
s1 <- s1[,!(colnames(s1) %in% useful.species)]


pca <- list(NULL)
pca$x <- s1[,5:length(s1)]	


s1$group = 0
group = 0
for(i in seq(dim(s1)[1])){
	if(all(s1$start[i] != s1$end[i-1]+1)){
		group = group+1
	}
	s1[i,'group'] <- group
}


all.wiggle <- NULL
for(PC in seq(dim(pca$x)[2])){
	pcom <- pca$x[,PC]
	Wiggle <- NULL
	for(z in seq(group)){
		scale_r <- pcom[s1$group == z] 
		zero <- rep(0, length(scale_r))
		wiggle <- data.frame(PC_score = zero, mean = zero , bottom = zero , top = zero)
		for(i in seq(along=scale_r)){	
			if((i == length(scale_r)) & (i == 1)){
				a <- mean(scale_r[i])
				s <- sd(pcom)
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
			wiggle$PC_score[i] <- scale_r[i]
			wiggle$mean[i] <- a 
			wiggle$bottom[i] <- left		
			wiggle$top[i] <- right
		}
		Wiggle <- rbind(Wiggle,wiggle )
	}
	all.wiggle <- c(all.wiggle,list(Wiggle))
}
names(all.wiggle) <- colnames(pca$x)


# now a loop to run through the entire list and annotate the position of each element


class <- NULL
for(i in seq(length(all.wiggle))){
	
	all.wiggle[[i]]$r_class <- rep("M", dim(all.wiggle[[i]])[1])
	all.wiggle[[i]]$r_class[all.wiggle[[i]]$bottom > 0] <- "H"
	all.wiggle[[i]]$r_class[all.wiggle[[i]]$top < 0] <- "L"
	class <- cbind(class,all.wiggle[[i]]$r_class)
}
class <- as.data.frame(class)
names(class)  <- names(all.wiggle)
class$binID <- s1$binID



#####
#
#   Get genome information
#
######




con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
sLengths <- SEQ$size
names(sLengths) <- SEQ[,'chr']
sLengths <- sLengths[names(sLengths) %in% as.character(unique(s1$chr))] 

s1.bin.gr <- GRanges(seqnames = Rle(s1$chr), 
			ranges = IRanges(start = s1$start, end = s1$end), 
			seqlengths = sLengths,
			binID = s1$binID) 
s1.chr <- GRanges(seqnames = Rle(names(sLengths)), 
			ranges = IRanges(start = 1, end = sLengths), 
			seqlengths = sLengths)





#####
#
#   Everything up until this point is fine
#
######


# up here we need to choose species we want
# then we can subset


# run the regression stuff from here

useful.species1 <- c("Horse","Dog", "Human", "Rat")

for(t in useful.species1){
analysis.spec <- useful.species1[!(useful.species1 == t)]

mob_class <- s1.mob[complete.cases(s1.mob[,c("binID",analysis.spec)]),]
TE_class <- class[complete.cases(s1.mob[,c("binID",analysis.spec)]),]

# subset the class table by the correct bin ID

overlap.pc <- NULL
overlap.sd <- NULL
overlap.bin <- NULL
for( i in seq(dim(TE_class)[2]-1)){
	
	
	group.sys <- c("HH", "HL", "LH" , "LL")
	group.sys1 <- c("H", "H", "L" , "L")
	group.sys2 <- c("H", "L", "H" , "L")
	for(o in analysis.spec){
		for(p in 1:4){
			assign(paste(o,".",group.sys[p], sep =""),  (TE_class[(mob_class[,o] == group.sys1[p]) & (TE_class[,i] == group.sys2[p]),"binID"]))
		}		
	}
	
	#pairwise comparisons to be made
	spec.no <- 1:length(analysis.spec)
	pairs <- NULL
	for(x in spec.no){
		pairs <- rbind(pairs, data.frame(rep(x,length(spec.no) - 1),spec.no[!(spec.no == x)]))
	}
	
	pairs <- data.frame(pair.1 = analysis.spec[pairs[,1]], pair.2 = analysis.spec[pairs[,2]])

	no.ol.g <- NULL
	no.ol.sd <- NULL
	no.ol.bin <- NULL
	for(z in seq(length(group.sys))){
		
		
		bins.no <- NULL
		for(x in 1:dim(pairs)[1]){
			bins <- length(get(paste(pairs[x,"pair.1"],group.sys[z], sep = ".")))
			bins.no <- c(bins.no, bins)
		}
		
		intersect.no <- NULL
		for(x in 1:dim(pairs)[1]){
			intersect <- length(intersect(get(paste(pairs[x,"pair.1"],group.sys[z], sep = ".")), get(paste(pairs[x,"pair.2"],group.sys[z], sep = "."))))
			intersect.no <- c(intersect.no, intersect)
		}
		
		no.OL <- mean(intersect.no / bins.no)
		no.ol.g <- c(no.ol.g, no.OL)
		no.ol.sd <- c(no.ol.sd, sd(intersect.no / bins.no))
		no.ol.bin <- c(no.ol.bin, mean(bins.no))
	}
	overlap.pc <- rbind(overlap.pc,no.ol.g)
	overlap.sd <- rbind(overlap.sd,no.ol.sd)
	overlap.bin <- rbind(overlap.bin,no.ol.bin)
	#plots	
#	for(x in 1:length(group.sys)){
#			q <- autoplot((s1.chr),layout = "karyogram",fill = "white", color = "black") 
#			for(y in 1:length(analysis.spec)){
#				range <- get(paste(analysis.spec[y], group.sys[x], sep="."))
#				q <- q + layout_karyogram(
#						reduce(s1.bin.gr[elementMetadata(s1.bin.gr)[,1] %in% range]), 
#						geom = "rect", 
#						fill = y, 
#						color = NA,
#						ylim = c(((y*(10/length(analysis.spec))) - (10/length(analysis.spec))) ,(y*(10/length(analysis.spec)))))
#			}
#			q <- q + labs(title = paste("Overlap between", group.sys2[x], names(TE_class)[i],"in", spec1,"regions and", group.sys1[x], "pairwsie mobility" ))				
#			assign(group.sys[x],q)
#		}
#		pdf(file = paste("plot_results/karyotype_plots/",spec1, "_", names(TE_class)[i], ".pdf", sep = ""), onefile=TRUE)
#		print(HH)
#		print(HL)
#		print(LH)
#		print(LL)
#	dev.off()

	
}

colnames(overlap.pc) <- group.sys  
rownames(overlap.pc) <- names(all.wiggle)
colnames(overlap.bin) <- group.sys  
rownames(overlap.bin) <- names(all.wiggle)
colnames(overlap.sd) <- group.sys  
rownames(overlap.sd) <- names(all.wiggle)




# probably do a heat map regression analysis 
# in each species we control for the effect of one of the species

pdf(file = paste("plot_results/TE_heatmap/",spec1, "_", "excluding_", t, ".pdf", sep = ""), onefile=TRUE)
print(heatmap.2((overlap.pc), 
		trace = "none", 
		density.info= "none", 
		margins = c(8,8), 
		main = paste(spec1, "TE regions\n overlaped with pairwise\n mobility, excluding",t ), 
		xlab = "TE mobility / TE level", 
		ylab = "TE", 
		breaks=seq(from=0, to=1, by=.01), 
		scale = "none",
		symbreaks=F,
		symkey=T,
		Colv=NA
		)
		)
		
		print(heatmap.2((overlap.bin), 
		trace = "none", 
		density.info= "none", 
		margins = c(8,8), 
		main = paste(spec1, "TE regions\n overlaped with pairwise\n mobility, excluding , bin no",t ), 
		xlab = "TE mobility / TE level", 
		ylab = "TE", 
		#breaks=seq(from=0, to=1, by=.01), 
		scale = "none",
		symbreaks=F,
		symkey=T,
		Colv=NA
		)
		)
		
		print(heatmap.2((overlap.sd), 
		trace = "none", 
		density.info= "none", 
		margins = c(8,8), 
		main = paste(spec1, "TE regions\n overlaped with pairwise\n mobility, excluding , bin no",t ), 
		xlab = "TE mobility / TE level", 
		ylab = "TE", 
		#breaks=seq(from=0, to=1, by=.01), 
		scale = "none",
		symbreaks=F,
		symkey=T,
		Colv=NA
		)
		)


dev.off()

#heatmap.2(cor(s1[,5:(length(s1) - 2)], replot_s1_s2[,2:length(replot_s1_s2)]), trace = "none",margins = c(8,8), xlab = spec2, ylab = spec1, main = "TE corelations between species", col = redgreen, scale = ("none"), symbreaks=TRUE, density.info= "none", breaks=seq(from=-1, to=1, by=.01), symkey=TRUE)



}














# find a way to do a significance test to compare the mean overlap

# wilcoxin signed rank test becausee we are sampling from a distribution and cant be sure if it is normmaly distributed or somthing like that

# paired will be equall to true to get a signed rank test



# we need to calculate the values for the total amount of overlap.   
# this approach should clear up the majority of questions about the amount of overlap
# I can use various species supposedly because im now comparing distributions rather than single numbers and the whole distribtibution will change according to the different species






































p.var <- names(class)

xlab <- "TE"
ylab <- "proportion of bins that overlap between two species in each TE level group"
main <- "proportion of overlap for each level of organisation for each level of each TE"

barplot(overlap.pc)

barplot(c(overlap.pc[,2],overlap.pc[,1],overlap.pc[,3],overlap.pc[,4]), col = c(2,1,3,4), xlab = xlab, ylab = ylab, main = main, ylim = c(0,1))
points(overlap.pc[,1], col = 1)
points(overlap.pc[,3], col = 3)
points(overlap.pc[,4], col = 4)


barplot(t(overlap.pc), col = c(1,2,3,4), cex.names = .3, )
legend("topright",group.sys, fill = c(1,2,3,4))

heatmap(cor(data.frame(p.var,overlap.pc)), scale = "none")

# also lets see whats happening with phastcon scores


# the correlation between species at a PC

barplot(t(overlap.pc[,1]), col = c(1), cex.names = .5, las = 2, beside = TRUE)
barplot(t(overlap.pc[,2]), col = c(2), cex.names = .5, las = 2, beside = TRUE)
barplot(t(overlap.pc[,3]), col = c(3), cex.names = .5, las = 2, beside = TRUE)
barplot(t(overlap.pc[,4]), col = c(4), cex.names = .5, las = 2, beside = TRUE)

barplot(t(overlap.pc), col = c(1,2,3,4), cex.names = .5, las = 2, beside = FALSE)


plot(hclust(dist(overlap.pc)))

# seperate things into three groups 
# ances   SL    erv


ances <- overlap.pc[c('SINE2_MIR', 'LINE_Other', 'LINE_L2', 'DNA_All', 'LINE_CR1'),]
SL <-  overlap.pc[c('SINE1_7SL', 'GC', 'NGenes', 'NCpGI', 'SINE_Other', 'ERV_ERV2', 'LINE_L1'),]
erv <- overlap.pc[c('LTR_All', 'ERV3_MaLR', 'ERV_Other', 'ERV_ERV1', 'ERV_ERV3'),]


dat_mean<- matrix(nrow = 4, ncol = 3)
colnames(dat_mean) = c("ances", "SL", "erv")
rownames(dat_mean) = c("HH","HL", "LH", "LL")

for(i in colnames(dat_mean)){
	for(z in rownames(dat_mean)){
		m_cal <- get(i)
		
		dat_mean[z,i] <- mean(m_cal[,z])
	}
}


dat_var<- matrix(nrow = 4, ncol = 3)
colnames(dat_var) = c("ances", "SL", "erv")
rownames(dat_var) = c("HH","HL", "LH", "LL")

for(i in colnames(dat_var)){
	for(z in rownames(dat_var)){
		v_cal <- get(i)
		
		dat_var[z,i] <- sd(v_cal[,z])
	}
}

dat_median<- matrix(nrow = 4, ncol = 3)
colnames(dat_median) = c("ances", "SL", "erv")
rownames(dat_median) = c("HH","HL", "LH", "LL")

for(i in colnames(dat_median)){
	for(z in rownames(dat_median)){
		m_cal <- get(i)
		
		dat_median[z,i] <- median(m_cal[,z])
	}
}




heatmap(dat_mean, scale= "none")

library(gplots)
heatmap.2((dat_mean), trace = "none", col = redgreen, density.info= "none",Rowv=NA, Colv=NA, margins = c(8,8))
heatmap.2((dat_median), trace = "none", col = redgreen, density.info= "none",Rowv=NA, Colv=NA, margins = c(8,8))


# what is the total level of conservation for various regions of organisation

Bov.HO <- as.integer(rownames(all.spec[(all.spec$Bovine_r_class == "H"),]))
Mou.HO <- as.integer(rownames(all.spec[(all.spec$Mouse_r_class == "H"),]))
Dog.HO <- as.integer(rownames(all.spec[(all.spec$Dog_r_class == "H"),]))

Bov.LO <- as.integer(rownames(all.spec[(all.spec$Bovine_r_class == "L"),]))
Mou.LO <- as.integer(rownames(all.spec[(all.spec$Mouse_r_class == "L"),]))
Dog.LO <- as.integer(rownames(all.spec[(all.spec$Dog_r_class == "L"),]))


	group.sys <- c("HO", "LO")
	class.2 <- c("H", "L", "H", "L")
	no.ol.g <- NULL
	for(z in seq(length(group.sys))){
		#ln.pc.cla <- length(class[class[,i] == class.2[z],i])
		B <- get(paste("Bov.",group.sys[z], sep = ""))
		M <- get(paste("Mou.",group.sys[z], sep = ""))
		D <- get(paste("Dog.",group.sys[z], sep = ""))
		BM <- B[B%in%M]
		MB <- M[M%in%B]
		DM <- D[D%in%M]
		MD <- M[M%in%D]
		BD <- B[B%in%D]
		DB <- D[D%in%B]
		no.OL <- sum(length(unique(c(BM,MB))),  length(unique(c(DM,MD))), length(unique(c(BD,DB)))) / sum(length(B), length(M), length(D))
		no.ol.g <- c(no.ol.g, no.OL)
	}
T.percent.ol <- no.ol.g
#names(T.percent.ol) <- c("HO","LO")


dat_mean2 <- as.data.frame(dat_mean)
dat_mean2$total <- c(rep(T.percent.ol[1],2),rep(T.percent.ol[2],2))
heatmap.2(as.matrix(dat_mean2[1:2,]), trace = "none", col = redgreen, density.info= "none",Rowv=NA, Colv=NA, margins = c(8,8))
heatmap.2(as.matrix(dat_mean2[3:4,]), trace = "none", col = redgreen, density.info= "none",Rowv=NA, Colv=NA, margins = c(8,8))


par(mfrow=c(2,2))
boxplot(ances[,'HH'], SL[,'HH'], erv[,'HH'], main = "HH", names = c("ances", "SL", "erv"))
abline(h = T.percent.ol[1], col = 2)
boxplot(ances[,'HL'], SL[,'HL'], erv[,'HL'], main = "HL", names = c("ances", "SL", "erv"))
abline(h = T.percent.ol[1], col = 2)
boxplot(ances[,'LH'], SL[,'LH'], erv[,'LH'], main = "LH", names = c("ances", "SL", "erv"))
abline(h = T.percent.ol[2], col = 3)
boxplot(ances[,'LL'], SL[,'LL'], erv[,'LL'], main = "LL", names = c("ances", "SL", "erv"))
abline(h = T.percent.ol[2], col = 3)



#probably do some venn diagrams next to see how exclusive each group is

# what are the parts that could be overlaping
# should i just get the numbers for everyhting and stick them in a big venn diagram 
library("VennDiagram")

# probably need to get the union of all the elements in the same groups



p1 <- qplot(x = mpg, y= cyl, data = mtcars, color = carb)
p2 <- qplot(x = mpg, y= cyl, data = mtcars, color = wt)
p3 <- qplot(x = mpg, y= cyl, data = mtcars, color = qsec)
p4 <- qplot(x = mpg, y= cyl, data = mtcars, color = gear)
arrangeGrobByParsingLegend(p1, p2, p3, p4)
arrangeGrobByParsingLegend(p1, p2, p3, p4, ncol = 1)
arrangeGrobByParsingLegend(p1, p2, p3, p4, legend.idx = 2)

