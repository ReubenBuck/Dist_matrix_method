

# proportion of variance at a PC

# corealtio nbetween sspecies at that PC


# what would a negative corelation mean 
# so we are using PC2 to divide species 


# waht is meant by us recognizing groups based on a metric and then being surprised when groups differ on that metric.


# Low PC2 and low org in diff species appear to align with each other to a greater degree than what could be expected by chance. 

# is this result real or the product of something else. 

# what would happen if we were to look at this through different PC window

# are we looking in too small a window 


# to work out how we go with subdividing it over each PCA will be a bit of a chalange

# proportion of overlap between any two species oro all three species what would be the best measure of conservation across the species 


# get the PCA meaure and run the window calculation across human

# for this we will need just the human PCA matrix and so for each human bin we will know whcih group it belongs to based on PCA 


setwd("~/Desktop/Human_and_Mouse_TE.div/")

spec1 <- "Human"
UCSCspec= "hg19"


rem.un <- "yes"

mb <- 1000000
bin.size = 500000

#which columns to keep 
keep.NGenes = "yes"
keep.NG4s = "no"
keep.NCpGI = "yes"
keep.CpGBP = "no"
keep.GC = "yes"
SCALE = "yes"

#create objects into whcih i will store the binsizes
s1name <- paste(spec1, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
slist <- list(s1)


for(i in seq(along=slist)){
      count <- slist[[i]]
      count <- count[count$Known >= bin.size,]
      count[,5:(length(count)-1)] <- (count[,5:(length(count)-1)]/count$Known)*mb
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

	
pca <- (prcomp(s1[,5:length(s1)]))
pca$x[,1] <- pca$x[,1]* -1
pca$x[,2] <- pca$x[,2]* -1
pca$rotation[,1] <- pca$rotation[,1]* -1
pca$rotation[,2] <- pca$rotation[,2]* -1
# maybe downloading joys scoring system


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
names(all.wiggle) <- paste("PC",seq(from= 1, to = dim(pca$x)[2]), sep ="")


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

# so now we have class identifiers for each pc for each bin

# next involves crossing this with coordinates



# read in the species

setwd("~/Desktop/Dist_matrix_TE_div/merged.tables/")


Bovine <- read.table("Human_Bovinemerged.CE.class", header = TRUE)
colnames(Bovine)[c(14, 16)] <- c("Bovine_Classification","Bovine_r_class")
Dog <- read.table("Human_Dogmerged.CE.class", header = TRUE)
colnames(Dog)[c(14, 16)] <- c("Dog_Classification","Dog_r_class")
Dog.class <- Dog[,c('bin.ID', 'Dog_Classification', 'Dog_r_class')]
Mouse <- read.table("Human_Mousemerged.CE.class", header = TRUE)
colnames(Mouse)[c(14, 16)] <- c("Mouse_Classification","Mouse_r_class")
Mouse.class <- Mouse[,c('bin.ID', 'Mouse_Classification', 'Mouse_r_class')]

all.spec <- merge(Bovine, Dog.class, by.x=1, by.y=1)
all.spec <- merge(all.spec, Mouse.class, by.x = 1, by.y = 1)


UCSCspec = "hg19"
con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
sLengths <- SEQ$size
names(sLengths) <- SEQ[,'chr']
sLengths <- sLengths[names(sLengths) %in% as.character(unique(all.spec$Chromosome))] 

hg.bin.gr <- GRanges(seqnames = Rle(all.spec$Chromosome), 
			ranges = IRanges(start = all.spec$Start, end = all.spec$End), 
			seqlengths = sLengths, CE_P = all.spec$P_CE) 
hg.chr <- GRanges(seqnames = Rle(names(sLengths)), 
			ranges = IRanges(start = 1, end = sLengths), 
			seqlengths = sLengths)

# subset the class table by the correct bin ID

class$binID <- s1$binID
class <- class[class$binID %in% all.spec$bin.ID,]
class <- class[order(class$binID),]
rownames(class) <- seq(dim(class)[1])

overlap.pc <- NULL
for( i in seq(dim(class)[2]-1)){
	
	Bov.HH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Bovine_r_class == "H"),]))
	Mou.HH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Mouse_r_class == "H"),]))
	Dog.HH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Dog_r_class == "H"),]))
	
	Bov.HL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Bovine_r_class == "H"),]))
	Mou.HL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Mouse_r_class == "H"),]))
	Dog.HL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Dog_r_class == "H"),]))
	
	Bov.LH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Bovine_r_class == "L"),]))
	Mou.LH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Mouse_r_class == "L"),]))
	Dog.LH <- as.integer(rownames(all.spec[(class[,i] == "H") & (all.spec$Dog_r_class == "L"),]))
	
	Bov.LL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Bovine_r_class == "L"),]))
	Mou.LL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Mouse_r_class == "L"),]))
	Dog.LL <- as.integer(rownames(all.spec[(class[,i] == "L") & (all.spec$Dog_r_class == "L"),]))
	
	# count the amount of bins that can be found in two species for each group at each pca class system
	group.sys <- c("HH", "HL", "LH" , "LL")
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
	
	overlap.pc <- rbind(overlap.pc,no.ol.g)
	
	#plots

#pdf(file = paste("../plot_results/circle_plots/", names(class)[i], ".pdf", sep = ""), onefile=TRUE)
#q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white", main = "high PC2 and high organization")
#q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7, main = "thing")
#q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
#q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
#q <- q + layout_circle(reduce(hg.bin.gr[Bov.HH]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Mou.HH]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Dog.HH]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
#q <- q + annotate("text",x = 0 , y = 0 ,label = paste("H.org H.", names(class)[i], sep = ""))
#print(q)


#q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
#q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
#q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
#q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
#q <- q + layout_circle(reduce(hg.bin.gr[Bov.HL]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Mou.HL]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Dog.HL]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
#q <- q + annotate("text",x = 0 , y = 0 ,label = paste("H.org L.", names(class)[i], sep = ""))
#print(q)

#q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
#q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
#q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
#q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
#q <- q + layout_circle(reduce(hg.bin.gr[Bov.LH]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Mou.LH]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Dog.LH]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
#q <- q + annotate("text",x = 0 , y = 0 ,label = paste("L.org H.", names(class)[i], sep = ""))
#print(q)


#q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
#q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
#q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
#q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
#q <- q + layout_circle(reduce(hg.bin.gr[Bov.LL]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Mou.LL]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
#q <- q + layout_circle(reduce(hg.bin.gr[Dog.LL]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
#q <- q + annotate("text",x = 0 , y = 0 ,label = paste("L.org L.", names(class)[i], sep = ""))
#print(q)

#dev.off()
	
}

colnames(overlap.pc) <- group.sys
rownames(overlap.pc) <- names(all.wiggle)

spca <- summary(pca)
spca <- as.data.frame(as.matrix(spca$importance))
p.var <- t(spca[2,])

xlab <- "proportion of variance at each PC"
ylab <- "proportion of bins that overlap between two species at each PC"
main <- "proportion of overlap for each level of organisation for each PC"

	
plot(p.var,overlap.pc[,1], type = "l", col = 1, xlab = xlab, ylab = ylab, main = main, ylim = c(0,1))
lines(p.var,overlap.pc[,2], col = 2)
lines(p.var,overlap.pc[,3], col = 3)
lines(p.var,overlap.pc[,4], col = 4)

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


