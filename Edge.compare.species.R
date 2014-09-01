

# read in each genome and ajust it

# so we need to have each genome and we pick the anceztral region in human and the low ancestral region
# for each speceis we take the difference of the mean edge lengths

library(ggplot2)
library(gplots)

rm(list= ls())
setwd("~/Desktop/Dist_matrix_TE_div/")


# constants
mb <- 1000000
bin.size = 500000


#option
# filter on s.bin file
# do not transfrom bins where less than x% of bin length in s2 is not aligned to a bin in s1
# eg. unique(S.bin_s1_s2[S.bin_s1_s2[,4] > filter,1])
filter = 0.1
# the min amount of bins a speies must have after filtration to be used in further analysis
bin.limit <- 1500


#####
# Load in Tables and process

spec1 <- "Human"
all.species <- c("Human", "Mouse", "Rat", "Dog", "Bovine", "Horse", "Opossum", "Platypus", "Elephant")
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


# load in the first species

keep.NGenes = "yes"
keep.NG4s = "no"
keep.NCpGI = "yes"
keep.CpGBP = "no"
keep.GC = "yes"
SCALE = "yes"

s1name <- paste("count_tables/",spec1, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
count <- s1
count <- count[count$Known >= bin.size,]
count$GC <- count$GC/count$Known*100
count[,5:(length(count)-2)] <- sqrt((count[,5:(length(count)-2)]/count$Known) * mb) 
if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}
colnames(count)[1:4] <- c("chr", "binID", "start", "end")
count$binID <- 1:dim(count)[1]
count <- count[,1:(length(count)-1)]
s1 <- count

# get each species adjusted count
useful.species1 <- NULL
for(x in 1:length(useful.species)){

#### Read in species2
spec2 <- useful.species[x]
s2name <- paste("count_tables/",spec2, "_AllBinCounts.txt", sep = "")
s2 <- read.table(s2name, header = TRUE)
count <- s2
count <- count[count$Known >= bin.size,]
count$GC <- count$GC/count$Known*100
count[,5:(length(count)-2)] <- sqrt((count[,5:(length(count)-2)]/count$Known) * mb) 
if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}
colnames(count)[1:4] <- c("chr", "binID", "start", "end")
count$binID <- 1:dim(count)[1]
count <- count[,1:(length(count)-1)]
s2 <- count

#### Transform species 2
S.bin.name <- paste("S_bin/", spec1, "_sbin/", spec1, "_aligning_", spec2, "_select", sep = "")
S.bin_s1_s2 <- read.table(S.bin.name,header = TRUE)
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2[,1] %in% unique(S.bin_s1_s2[S.bin_s1_s2[,4] > filter,1]),]
s2TE <- s2[,c(2, 5:length(s2))]
merg.s1_s2 <- merge(S.bin_s1_s2, s2TE, by.x=3, by.y=1)	
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
rownames(replot_s1_s2) <- 1:dim(replot_s1_s2)[1]
rownames(s1) <- 1:dim(s1)[1]
if(dim(replot_s1_s2)[1] >= bin.limit){
	assign(spec2, replot_s1_s2)
	useful.species1 <- c(useful.species1, spec2) 
}
}
useful.species <- useful.species1
# maybe need to remove opossum because any bit of filtering takes it out of the running

IDs <- get(useful.species[1])[,"S1_bin"]
for(i in 2:length(useful.species)){
	IDs <- intersect(IDs, get(useful.species[i])[,"S1_bin"])
}
for(i in 1:length(useful.species)){
	spec <- get(useful.species[i])
	spec <- spec[spec[,"S1_bin"] %in% IDs,]
	spec <- spec[order(spec$S1_bin), ]
	assign(useful.species[i],spec)	
}





# we work out the regions for human and then see how each species looks in human regions
# do it on the original s1
pca <- prcomp(s1[,5:length(s1)])
if(sqrt(pca$rotation["SINE2_MIR","PC1"]^2) > sqrt(pca$rotation["SINE2_MIR","PC2"]^2)){
	ancestral <- 1
	gene.no <- 2
}else{
	ancestral <- 2
	gene.no <- 1
	}
if(cor(pca$x[,ancestral], s1$SINE2_MIR) < 0){
	pca$x[,ancestral] <- pca$x[,ancestral] * -1
	pca$rotation[,ancestral] <- pca$rotation[,ancestral] * -1
}
if(cor(pca$x[,gene.no], s1$NGenes) < 0){
	pca$x[,gene.no] <- pca$x[,gene.no] * -1
	pca$rotation[,gene.no] <- pca$rotation[,gene.no] * -1
}
	
PCs <- data.frame(gene.no = pca$x[,gene.no], ancestral = pca$x[,ancestral])
assign(paste(spec1, "PCs", sep = "_"), PCs)



Group = 0
group = 0
for(i in seq(dim(s1)[1])){
	if(all(s1$start[i] != s1$end[i-1]+1)){
		group = group+1
	}
	Group[i] <- group
}

pca$x <- PCs
all.wiggle <- NULL
for(PC in seq(dim(pca$x)[2])){
	pcom <- pca$x[,PC]
	Wiggle <- NULL
	for(z in seq(group)){
		scale_r <- pcom[Group == z] 
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

class.comp <- class[class$binID %in% IDs,]
class.comp <- class.comp[order(class.comp$binID),]

# get the s1 to be one of the species now

assign("Human", s1[s1$binID %in% IDs,c(2,5:length(s1))])

# simply take the things that are H and L and compare the edge lengths

# Mouse.AH.dist


A.H <- class.comp[class.comp["ancestral"] == "H","binID"]
A.L <- class.comp[class.comp["ancestral"] == "L","binID"]
G.H <- class.comp[class.comp["gene.no"] == "H","binID"]
G.L <- class.comp[class.comp["gene.no"] == "L","binID"]

region <- c("A.H", "A.L", "G.H", "G.L")

for(r in region){
	reg <- get(r)
	edge.df <- data.frame(bin_no = reg)
	for(s in useful.species){
		spec <- get(s)
		Dist <- colMeans(as.matrix(dist(spec[spec$S1_bin %in% reg,2:length(spec)])))
		edge.df <- cbind(edge.df, Dist)
	}
	colnames(edge.df) <- c("bin_no", useful.species)
	assign(paste(r, "_edges", sep = ""), edge.df)
}




for(r in region){
	reg <- get(r)
	edge.df <- data.frame(bin_no = reg)
	spec <- s1
	Dist <- colMeans(as.matrix(dist(spec[spec$binID %in% reg,5:length(spec)])))
	edge.df <- cbind(edge.df, Dist)
	colnames(edge.df) <- c("bin_no", spec1)
	assign(paste(r, "Human_edges", sep = ""), edge.df)
}


for(r in region){
	edge <- get(paste(r, "_edges", sep=""))
	H.edge <- get(paste(r, "Human_edges", sep=""))
	for(i in 2:length(edge)){
		edge[i] <- edge[i] - H.edge[,2]
	}
	assign(paste(r,".diff", sep = ""), edge)
}





A.H.all <- merge(A.H_edges, A.HHuman_edges, by.x = 1 , by.y = 1)
A.L.all <- merge(A.L_edges, A.LHuman_edges, by.x = 1 , by.y = 1)
G.H.all <- merge(G.H_edges, G.HHuman_edges, by.x = 1 , by.y = 1)
G.L.all <- merge(G.L_edges, G.LHuman_edges, by.x = 1 , by.y = 1)




heatmap(as.matrix(A.H_edges[,2:length(A.H_edges)]), scale = "none")
heatmap(as.matrix(A.L_edges[,2:length(A.H_edges)]), scale = "none")
heatmap(as.matrix(G.H_edges[,2:length(A.H_edges)]), scale = "none")
heatmap(as.matrix(G.L_edges[,2:length(A.H_edges)]), scale = "none")



heatmap(as.matrix(A.H.all[,2:length(A.H.all)]), scale = "none")
heatmap(as.matrix(A.L.all[,2:length(A.L.all)]), scale = "none")
heatmap(as.matrix(G.H.all[,2:length(G.H.all)]), scale = "none")
heatmap(as.matrix(G.L.all[,2:length(G.L.all)]), scale = "none")


heatmap.2(as.matrix(A.H.diff[,2:length(A.H.diff)]),
		  scale = "none", 
		  col = "greenred" ,
		  density.info = "none",
		  symkey = T,
		  symbreaks = T,
		  trace = "none",
		  breaks = seq(from = -5, to = +5, by = 0.1),
		  margins = c(8,4),
		  dendrogram = "column",
		  main = "High PC2",
		  labRow = NA,
		  xlab = "species",
		  ylab = "bins"
		  )

heatmap.2(as.matrix(A.L.diff[,2:length(A.H.diff)]),
		  scale = "none", 
		  col = "greenred" ,
		  density.info = "none",
		  symkey = T,
		  symbreaks = T,
		  trace = "none",
		  breaks = seq(from = -5, to = +5, by = 0.001),
		  margins = c(8,4),
		  dendrogram = "column",
		  main = "Low PC2",
		  labRow = NA,
		  xlab = "species",
		  ylab = "bins"
		  )
		  
heatmap.2(as.matrix(G.H.diff[,2:length(A.H.diff)]),
		  scale = "none", 
		  col = "greenred" ,
		  density.info = "none",
		  symkey = T,
		  symbreaks = T,
		  trace = "none",
		  breaks = seq(from = -5, to = +5, by = 0.001),
		  margins = c(8,4),
		  dendrogram = "column",
		  main = "High PC1",
		  labRow = NA,
		  xlab = "species",
		  ylab = "bins"
		  )
		  
heatmap.2(as.matrix(G.L.diff[,2:length(A.H.diff)]),
		  scale = "none", 
		  col = "greenred" ,
		  density.info = "none",
		  symkey = T,
		  symbreaks = T,
		  trace = "none",
		  breaks = seq(from = -5, to = +5, by = 0.001),
		  margins = c(8,4),
		  dendrogram = "column",
		  main = "Low PC1",
		  labRow = NA,
		  xlab = "species",
		  ylab = "bins"
		  )		  		 
		  
		  
		  
# maybe also plot these using ggplots 2
# just give every other bin a zero


# if we want to see where these bins are moving 


heatmap.2(cor(Elephant[,-1],Bovine[,-1]), col = "redgreen", margins=c(10,10), trace = "none", xlab = "Bovine", ylab = "Elephant", scale = "none")
heatmap.2(cor(Dog[,-1],Bovine[,-1]), col = "redgreen", margins=c(10,10), trace = "none", xlab = "Bovine", ylab = "Dog", scale = "none")
heatmap.2(cor(Elephant[,-1],Dog[,-1]), col = "redgreen", margins=c(10,10), trace = "none", xlab = "Dog", ylab = "Elephant", scale = "none")

heatmap.2(cor(Elephant[,-1],Human[,-1]), col = "redgreen", margins=c(10,10), trace = "none", xlab = "Bovine", ylab = "Elephant", scale = "none")

# this part can be looped and the images can be saved so we get a feel for everything




# make all the pca files
Human.pca <- prcomp(Human[,2:length(Human)])
Elephant.pca <- prcomp(Elephant[,2:length(Elephant)])
Bovine.pca <- prcomp(Bovine[,2:length(Bovine)])
Mouse.pca <- prcomp(Mouse[,2:length(Mouse)])
Horse.pca <- prcomp(Horse[,2:length(Horse)])
Dog.pca <- prcomp(Dog[,2:length(Dog)])





# make all the every files
for(r in region){
	every <- data.frame(binID <- IDs, matrix(data = NA,ncol = length(useful.species), nrow = length(IDs)))
	diff <- get(paste(r,".diff", sep =""))
	every[every[,1] %in% diff[,1],2:dim(every)[2]] <- diff[,2:dim(every)[2]]
	colnames(every)[2:dim(every)[2]] <- useful.species
	assign(paste(r, "_every", sep = ""), every)
}

# run the graphing through all the everies and all the speices

#biplot(Human.pca, xlabs = rep(".", length(IDs)), main = "Human")
#biplot(Elephant.pca, xlabs = rep(".", length(IDs)), main = "Elephant")

for(s in useful.species){
	spec <- get(paste(s,".pca", sep =""))
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "biplot_", s,".pdf", sep = ""))
	biplot(spec, xlabs = rep(".", length(IDs)), main = s)
	dev.off()
}

pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "biplot_", "Human",".pdf", sep = ""))
biplot(Human.pca, xlabs = rep(".", length(IDs)), main = "Human")
dev.off()



for(s in useful.species){
	
	qplot(Human.pca$x[,1], Human.pca$x[,2], color = A.H_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste("Human_", s, "_A.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "Human.pca_spec_A.H_", s,".pdf", sep = ""))
	
	spec <- get(paste(s,".pca", sep =""))
	
	qplot(spec$x[,1], spec$x[,2], color =  A.H_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste(s, "_A.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "spec_A.H_", s,".pdf", sep = ""))

}

for(s in useful.species){
	
	qplot(Human.pca$x[,1], Human.pca$x[,2], color = A.L_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste("Human_", s, "_A.L"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "Human.pca_spec_A.L_", s,".pdf", sep = ""))
	
	spec <- get(paste(s,".pca", sep =""))
	
	qplot(spec$x[,1], spec$x[,2], color =  A.L_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste(s, "_A.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "spec_A.L_", s,".pdf", sep = ""))

}

for(s in useful.species){
	
	qplot(Human.pca$x[,1], Human.pca$x[,2], color = G.H_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste("Human_", s, "_G.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "Human.pca_spec_G.H_", s,".pdf", sep = ""))
	
	spec <- get(paste(s,".pca", sep =""))
	
	qplot(spec$x[,1], spec$x[,2], color =  G.H_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste(s, "_A.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "spec_G.H_", s,".pdf", sep = ""))

}

for(s in useful.species){
	
	qplot(Human.pca$x[,1], Human.pca$x[,2], color = G.L_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste("Human_", s, "_G.L"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "Human.pca_spec_G.L_", s,".pdf", sep = ""))
	
	spec <- get(paste(s,".pca", sep =""))
	
	qplot(spec$x[,1], spec$x[,2], color =  G.L_every[,s]) + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) + labs(title = paste(s, "_A.H"))
	ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/", "spec_G.L_", s,".pdf", sep = ""))

}

