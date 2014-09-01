
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



## so we got all TEs of a certain class in H, M, L

# what we should really do is look at the distribtuion of edges within each class


all.dist.d <- NULL
all.hist.d <- NULL
g_plot.d <- NULL
for(i in 1:(length(class)-1)){
	
	dist_d <- dist(s1[class[,i] == "H",5:14])
	
	all.dist.d <- c(all.dist.d, list(dist_d))
	all.hist.d <- c(all.hist.d, list(hist(dist_d)))
	dist_d <- data.frame(as.numeric(dist_d),colnames(class)[i])
	g_plot.d <- rbind(g_plot.d, dist_d)
}
colnames(g_plot.d) <- c("edge_length", "TE_name")

plot(all.hist.d[[1]]$mids, all.hist.d[[1]]$density,type="l")

for(i in 1:(length(all.hist.d)-1)){
	lines(all.hist.d[[i+1]]$mids, all.hist.d[[i+1]]$density,col = i+1)
}

test <- g_plot.d[g_plot.d[,"TE_name"] %in% c("LINE_L2", "ERV_ERV1", "ERV_ERV2", "ERV_ERV3", "SINE1_7SL", "LINE_L1"),]
ggplot(test, aes(edge_length, fill = TE_name)) + geom_density(alpha = 0.2)



test <- g_plot.d[g_plot.d[,"TE_name"] %in% c("LINE_L2",  "ERV_ERV2"),]
ggplot(test, aes(edge_length, fill = TE_name)) + geom_density(alpha = 0.2)



A <- NULL
for(i in 1:length(class)){
	A <- c(A,print(dim(class[class[,i] == "H",])[1]))
	
}


pca <- prcomp(s1[,c("LINE_L2", "ERV_ERV1", "ERV_ERV2", "ERV_ERV3", "SINE1_7SL", "LINE_L1")])





# maybe start to compare across species so we can get a result
cols <- rep(1,dim(class)[1])
cols[class$ERV_ERV2 == "H"] <- 2
cols[class$SINE1_7SL == "H"] <- 3
cols[class$LINE_L2 == "H"] <- 4
cols[class$LINE_L1 == "H"] <- 5
cols[class$ERV_ERV1 == "H"] <- 6
cols[class$ERV_ERV3 == "H"] <- 7

plot(pca$x[cols > 1,c(1,2)], col = cols[cols>1], pch = 16, cex=.5)









