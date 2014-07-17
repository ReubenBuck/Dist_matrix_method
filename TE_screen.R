

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

spec1 <- "Human"
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
useful.species1 <- useful.species1[!(useful.species1 == "Opossum")]


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

#useful.species1 <- c("Horse","Dog", "Human", "Rat")



for(t in useful.species1){
analysis.spec <- useful.species1[!(useful.species1 == t)]

mob_class <- s1.mob[complete.cases(s1.mob[,c("binID",analysis.spec)]),c("binID",analysis.spec)]
TE_class <- class[complete.cases(s1.mob[,c("binID",analysis.spec)]),]

# subset the class table by the correct bin ID

overlap.pc <- NULL
overlap.sd <- NULL
overlap.bin <- NULL
overlap.wi <- NULL
overlap.bs <- NULL
overlap.di <- NULL
overlap.med <- NULL
for( i in seq(dim(TE_class)[2]-1)){
	
	
	group.sys <- c("HH", "HL", "LH" , "LL")
	group.sys1 <- c("H", "H", "L" , "L")
	group.sys2 <- c("H", "L", "H" , "L")
	
	# make total mob overlaps here
	for(o in analysis.spec){
		for(p in 1:4){
			# remember that the TE_class object is irrelevent, it is hust where the binIDs are
			assign(paste(o,".",group.sys1[p],"O", sep =""),  (TE_class[(mob_class[,o] == group.sys1[p]),"binID"]))
		}		
	}
	
	# we don't do total TE overlaps because we are looking through the window of one species
	# food for thought is what if we were to include the remodeled TE distributions from the dist analysis
		
	# make TE overlaps here
	for(o in analysis.spec){
		for(p in 1:4){
			assign(paste(o,".",group.sys[p], sep =""),  (TE_class[(mob_class[,o] == group.sys1[p]) & (TE_class[,i] == group.sys2[p]),"binID"]))
			assign(paste(o,".",group.sys[p],".",names(TE_class)[i], sep =""),  (TE_class[(mob_class[,o] == group.sys1[p]) & (TE_class[,i] == group.sys2[p]),"binID"]))
		}		
	}

	#pairwise comparisons to be made
	spec.no <- 1:length(analysis.spec)
	pairs <- NULL
	for(x in spec.no){
		pairs <- rbind(pairs, data.frame(rep(x,length(spec.no) - 1),spec.no[!(spec.no == x)]))
	}
	
	pairs <- data.frame(pair.1 = analysis.spec[pairs[,1]], pair.2 = analysis.spec[pairs[,2]])
	
	
	# once i have pairs data i can get the totals 
	# this thing has to split into high and low
	# it already is, I just need to have four coppies of the distribtuion
	
	
	
	# get the TE overlaps here
	no.ol.g <- NULL
	no.ol.sd <- NULL
	no.ol.bin <- NULL
	no.ol.w <- NULL
	no.ol.bin.sd <- NULL 
	no.ol.d <- NULL
	no.ol.med <- NULL
	for(z in seq(length(group.sys))){
		
		bins.no.o <- NULL
		for(x in 1:dim(pairs)[1]){
			bins <- length(get(paste(pairs[x,"pair.1"], ".", group.sys1[z],"O", sep = "")))
			bins.no.o <- c(bins.no.o, bins)
		}

		bins.no <- NULL
		for(x in 1:dim(pairs)[1]){
			bins <- length(get(paste(pairs[x,"pair.1"],group.sys[z], sep = ".")))
			bins.no <- c(bins.no, bins)
		}
		
		intersect.no.o <- NULL
		for(x in 1:dim(pairs)[1]){
			intersect <- length(intersect(get(paste(pairs[x,"pair.1"], ".", group.sys1[z],"O", sep = "")), get(paste(pairs[x,"pair.2"], ".", group.sys1[z],"O", sep = ""))))
			intersect.no.o <- c(intersect.no.o, intersect)
		}

		
		intersect.no <- NULL
		for(x in 1:dim(pairs)[1]){
			intersect <- length(intersect(get(paste(pairs[x,"pair.1"],group.sys[z], sep = ".")), get(paste(pairs[x,"pair.2"],group.sys[z], sep = "."))))
			intersect.no <- c(intersect.no, intersect)
		}
		
		# may need a conditional for when wilcox dont work
		if(any(bins.no ==0)){
			Wil <- NULL
			Wil$p.value <- 1
		}else{
		Wil <- wilcox.test(intersect.no.o / bins.no.o,intersect.no / bins.no, paired= TRUE)
		}
		
		no.ol.w <- c(no.ol.w, Wil$p.value)
		no.ol.d <- c(no.ol.d, mean(intersect.no / bins.no) - mean(intersect.no.o / bins.no.o))
		no.ol.g <- c(no.ol.g, mean(intersect.no / bins.no))
		no.ol.sd <- c(no.ol.sd, sd(intersect.no / bins.no))
		no.ol.bin <- c(no.ol.bin, mean(bins.no))
		no.ol.bin.sd <- c(no.ol.bin.sd, sd(bins.no))
		no.ol.med <- c(no.ol.med, median(intersect.no / bins.no) - median(intersect.no.o / bins.no.o))
	}
	overlap.wi <- rbind(overlap.wi, no.ol.w)
	overlap.pc <- rbind(overlap.pc,no.ol.g)
	overlap.sd <- rbind(overlap.sd,no.ol.sd)
	overlap.bin <- rbind(overlap.bin,no.ol.bin)
	overlap.bs <- rbind(overlap.bs,no.ol.bin.sd)
	overlap.di <- rbind(overlap.di,no.ol.d )
	overlap.med <- rbind(overlap.med,no.ol.med )
	#plots	
	
	#karyotype plots
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
#			q <- q + labs(title = paste("Overlap between", group.sys2[x], names(TE_class)[i],"in", #spec1,"regions and", group.sys1[x], "pairwsie mobility" ))				
#			assign(group.sys[x],q)
#		}
#		pdf(file = paste("plot_results/karyotype_plots/",spec1, "_", names(TE_class)[i], #"_no_scale.pdf", sep = ""), onefile=TRUE)
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
colnames(overlap.wi) <- group.sys  
rownames(overlap.wi) <- names(all.wiggle)
colnames(overlap.bs) <- group.sys  
rownames(overlap.bs) <- names(all.wiggle)
colnames(overlap.di) <- group.sys  
rownames(overlap.di) <- names(all.wiggle)
colnames(overlap.med) <- group.sys  
rownames(overlap.med) <- names(all.wiggle)

# probably do a heat map regression analysis 
# in each species we control for the effect of one of the species


makeRects <- function(tfMat,border){
   cAbove = expand.grid(1:dim(overlap.wi)[1],1:dim(overlap.wi)[2]) #[tfMat,]
   cAbove = cAbove[tfMat,]
   cAbove <- cAbove[,c("Var2", "Var1")]
   cAbove[,2] <- dim(overlap.wi)[1] - cAbove[,2] + 1
   xl=cAbove[,1]-0.49
   yb=cAbove[,2]-0.49
   xr=cAbove[,1]+0.49
   yt=cAbove[,2]+0.49
   rect(xl,yb,xr,yt,border=border,lwd=3)
 }

sig1 <- overlap.wi<.05 & overlap.wi >=.01
sig2 <- overlap.wi<.01 & overlap.wi >=.001
sig3 <- overlap.wi<.001

overlap.bin <- data.frame(overlap.bin)
overlap.bs <- data.frame(overlap.bs)
bin.numbers <- data.frame(HH_bin_mean = overlap.bin$HH,
						  HH_bin_sd = overlap.bs$HH,
						  HL_bin_mean = overlap.bin$HL,
						  HL_bin_sd = overlap.bs$HL,
						  LH_bin_mean = overlap.bin$LH,
						  LH_bin_sd = overlap.bs$LH,
						  LL_bin_mean = overlap.bin$LL,
						  LL_bin_sd = overlap.bs$LL)
rownames(bin.numbers) <- rownames(overlap.bs)


pdf(file = paste("plot_results/TE_heatmap/",spec1, "_", "excluding_", t, "_scale.pdf", sep = ""), onefile=TRUE)
par(mfrow=c(1,1))
print(heatmap.2((overlap.pc), 
		trace = "none", 
		density.info= "none", 
		margins = c(8,8), 
		xlab = "TE mobility / TE level", 
		ylab = "TE", 
		breaks=seq(from=0, to=1, by=.01), 
		scale = "none",
		symbreaks=T,
		symkey=T,
		Colv=NA
		)
		)
		legend(5,5,NA)
title(main = paste(spec1, "mean overlap of TE regions\n and pairwise mobility regions"))
text(.5,.95, label = paste( "excluding",t ))
 		
		
print(heatmap.2((overlap.di), 
				trace = "none",
				margins = c(8,8), 
				col = redgreen, 
				scale = "none", 
				symbreaks=TRUE, 
				density.info= "none", 
 				Colv=NA,
 				Rowv = NA,
				breaks=seq(from=-.3, to=.3, by=.01),
				symkey=TRUE,
				xlab = "bin mobility / TE level", 
				ylab = "TE", 
 				add.expr = {makeRects(sig3,"orange");add.expr = makeRects(sig2,"pink");add.expr = makeRects(sig1,"yellow")}
 					)
 					)
		legend("topright", c("p < .05", "p < .01", "p < .001"), fill = c("yellow", "pink", "orange"))
title(main = paste(spec1, "TE regions\n overlaped with pairwise mobility regions"))
text(.5,.95, label = paste("mean divergence from mean overall,\n excluding",t ))


print(heatmap.2((overlap.med), 
				trace = "none",
				margins = c(8,8), 
				col = redgreen, 
				scale = "none", 
				symbreaks=TRUE, 
				density.info= "none", 
 				Colv=NA,
 				Rowv = NA,
				breaks=seq(from=-.3, to=.3, by=.01),
				symkey=TRUE,
				xlab = "bin mobility / TE level", 
				ylab = "TE", 
 				add.expr = {makeRects(sig3,"orange");add.expr = makeRects(sig2,"pink");add.expr = makeRects(sig1,"yellow")}
 					)
 					)
		legend("topright", title = "wilcoxon test",c("p < .05", "p < .01", "p < .001"), fill = c("yellow", "pink", "orange"))
title(main = paste(spec1, "TE regions\n overlaped with pairwise mobility regions"))
text(.5,.95, label = paste("median divergence from median overall,\n excluding",t ))



par(mfrow=c(2,1))

textplot(signif(bin.numbers[,1:4], digits = 3), mar=c(0, 0, 3, 0), main = "hmm")
title(main = "bin numbers in each overlap")
textplot(signif(bin.numbers[,5:8], digits = 3), mar=c(3,0, 0, 0),)



dev.off()


}




# what about a really large venn diagram, every thing should be in object from floating around it wouldnt be too hard to actually do a venn diagram

# I could get everything and put it in a list
# what am i wanting to answear 

# we've found which TEs in human have significant overlap across species
# in a venn diagram situation we can get a sense of what are TEs these groups consist of. 

# 15 TEs, 4 species, 4 classification groups
# therefore there will be 240 catagories
# and heaps of intersections 
# could there be a better way of storing such data?

# maybe on a pairwise level and do multiple ones so we can get an understanding of whats happening



library(VennDiagram)

# so far I loop through TE names and overwirte them. 
# It wouldnt be to hard to store it into memmory 

# also venn diagram usage is limited
# what would be good to see is intersections across each TE 

HHHH <- list(Bovine.HH,Mouse.HH,Horse.HH, Dog.HH,Bovine.LH)
venn(HHHH) 
fill = c("darkred", "darkblue", "darkgreen","darkyellow"))


# calculate the length of the intersect in each case 
# of course this approach will work
# the amount that are similar between species one and two 
# focus on getting all the TE distributions for each species 
# could put it up above




# now i need a way to get all the intersects

# maybe the expand grid method will get us all the compairsons


grid <- expand.grid(names(all.wiggle),names(all.wiggle))

for(s in useful.species1){
	for(r in group.sys){
		Sect= NULL
		for(i in 1:dim(grid)[1]){
			sect <- length(intersect(get(paste(s,".",r,".",grid[i,1], sep = "")),get(paste(s,".",r,".",grid[i,2], sep = "")))) / length(get(paste(s,".",r,".",grid[i,1], sep = "")))
			Sect <- c(Sect,sect)
		}
		inters <- matrix(data = Sect,nrow = 14, ncol = 14)
		colnames(inters) <- rownames(inters) <- names(all.wiggle)
		assign(paste("inters",s,r, sep = "."), inters)
	}
}

# if we have to plot this we need to know wchich direction it is going in
# who are intersexctors and who are intersecties 



heatmap.2(inters.Dog.LH, scale = "none", density.info = "none", trace = "none", margins = c(8,8))








# make a bigger grid

nam <- names(all.wiggle)
nam2 <- c(paste(".HH.", nam, sep = ""), paste(".HL.", nam, sep = ""), paste(".LH.", nam, sep = ""),paste(".LL.", nam, sep = ""))



grid <- expand.grid(nam2,nam2)

for(s in useful.species1){
		Sect= NULL
		for(i in 1:dim(grid)[1]){
			sect <- length(intersect(get(paste(s,grid[i,1], sep = "")),get(paste(s,grid[i,2], sep = "")))) / length(get(paste(s,grid[i,1], sep = "")))
			Sect <- c(Sect,sect)
		}
		inters <- matrix(data = Sect,nrow = length(nam2))
		colnames(inters) <- rownames(inters) <- names(nam2)
		assign(paste("inters",s,"all", sep = "."), inters)
}

# if we have to plot this we need to know wchich direction it is going in
# who are intersexctors and who are intersecties 

heatmap.2(inters, scale = "none", density.info = "none", trace = "none", margins = c(10,10), Colv=NA,Rowv=NA)



nam <- names(all.wiggle)
nam2 <- c(paste(".HH.", nam, sep = ""), paste(".HL.", nam, sep = ""), paste(".LH.", nam, sep = ""),paste(".LL.", nam, sep = ""))
nam3 <- c(paste("Bovine", nam2, sep = ""), paste("Mouse", nam2, sep = ""), paste("Dog", nam2, sep = ""),paste("Horse", nam2, sep = ""))



grid <- expand.grid(nam3,nam3)


		Sect= NULL
		for(i in 1:dim(grid)[1]){
			sect <- length(intersect(get(paste(grid[i,1], sep = "")),get(paste(grid[i,2], sep = "")))) / length(get(paste(grid[i,1], sep = "")))
			Sect <- c(Sect,sect)
		}
		inters <- matrix(data = Sect,nrow = length(nam3))
		colnames(inters) <- rownames(inters) <- nam3
		assign(paste("inters","all.spec", sep = "."), inters)



