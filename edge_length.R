# need to compare the edge lengths for different things
# need to sort out what I'm looking for 

# we can do it across both PCs

# take joys data apply the same transformations and maybe look at the first few PCs
# we see if the edges within them are conserved or unconserved. 

# this is because a change in the edge lengths reflects a change in an underlying distribution

# do this to the H and L regions of each PC from 1 to 3
# PC 2 would be easy becaus the data has been prepared 

# so read in each species perfomr the transformations and define the regions 
# we end up wiht several tables describing each species 



# For every species this can be done 


# code up a bit of script that changes the sign if needed. 
# essentailly the sign should corealte with something we know is conserved for PC2 there should be a positive corelation with MIR
# for PC1 we can use GC or gene.no
# not sure for PC3 yet, havnt had a good look at what assoicates with it.


# need to get all the bin counting stuff onto the console
# only use the stuff that is count data

# we will use only , NGenes, NCpGI, GC
# we already have a seperate thing going on with G4s and they don't appear to interesting in our context 
rm(list = ls())

library(cluster)
library(vioplot)
library(ggplot2)

setwd("~/Desktop/Dist_matrix_TE_div/")

keep.NGenes = "yes"
keep.NG4s = "no"
keep.NCpGI = "yes"
keep.CpGBP = "no"
keep.GC = "yes"
SCALE = "yes"

rem.un <- "yes"

mb <- 1000000
bin.size = 500000

#"Elephant" is sstrange for now, maybe bad scaffold

all.species <- c("Human", "Mouse", "Dog", "Bovine",  "Elephant", "Horse")

for(i in seq(along = all.species)){
	#create objects into whcih i will store the binsizes
	spec1 <- all.species[i]
	
	s1name <- paste("count_tables/",spec1, "_AllBinCounts.txt", sep = "")
	s1 <- read.table(s1name, header = TRUE)
	count <- s1	
    count <- count[count$Known >= bin.size,]
	count$GC <- count$GC/count$Known*100
    count[,5:(length(count)-2)] <- sqrt((count[,5:(length(count)-2)]/count$Known) * mb) 
#	count[,5:(length(count)-2)] <- ((count[,5:(length(count)-2)]/count$Known) * mb)
    if(keep.NGenes == "no"){count <- count[,!(colnames(count) == "NGenes")]}
    if(keep.NG4s ==	"no"){count <- count[,!(colnames(count) == "NG4s")]}
    if(keep.NCpGI ==	"no"){count <- count[,!(colnames(count) == "NCpGI")]}
    if(keep.CpGBP  ==	"no"){count <- count[,!(colnames(count) == "CpGBP")]}
    if(keep.GC ==	"no"){count <- count[,!(colnames(count) == "GC")]}
	if(SCALE == "yes"){count[,5:(length(count)-1)] <- scale(count[,5:(length(count)-1)])}

    #count <- count[,!(colnames(count) == "Known")]
    colnames(count)[1:4] <- c("chr", "binID", "start", "end")
    count$binID <- 1:dim(count)[1]
    assign(spec1, count )
	assign(paste(spec1, "data", sep="."), count[,c(2,5:(length(count)-1))])
	# Do the PCA stuff for these guys too
	# PCA them and change sign if i need to 
	
	pca <- prcomp(count[,5:(length(count)- 1)])
	if(sqrt(pca$rotation["SINE2_MIR","PC1"]^2) > sqrt(pca$rotation["SINE2_MIR","PC2"]^2)){
		ancestral <- 1
		gene.no <- 2
	}else{
		ancestral <- 2
		gene.no <- 1
		}
	if(cor(pca$x[,ancestral], count$SINE2_MIR) < 0){
		pca$x[,ancestral] <- pca$x[,ancestral] * -1
		pca$rotation[,ancestral] <- pca$rotation[,ancestral] * -1
	}
	if(cor(pca$x[,gene.no], count$NGenes) < 0){
		pca$x[,gene.no] <- pca$x[,gene.no] * -1
		pca$rotation[,gene.no] <- pca$rotation[,gene.no] * -1
	}
	
	assign(paste(spec1, ".pca", sep = ""),pca)
	
	PCAs <- data.frame(gene.no = pca$x[,gene.no], ancestral = pca$x[,ancestral])
	assign(paste(spec1, "PCAs", sep = "_"), PCAs)
	
	PCs <- count[,5:(length(count)-1)]
	assign(paste(spec1, "PCs", sep = "_"), PCs)


}


# What happens if I take bovB out of elephant and cow genomes


# now for the PCs i need to run them trough their sig regions
# lets grab the sig L2MIR regions, the intercept and union


CLASS = NULL
for(s in seq(along = all.species)){
s1 <- get(all.species[s])

s1$group = 0
group = 0
for(i in seq(dim(s1)[1])){
	if(all(s1$start[i] != s1$end[i-1]+1)){
		group = group+1
	}
	s1[i,'group'] <- group
}

pca$x <- get(paste(all.species[s], "PCs", sep = "_"))

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

CLASS <- c(CLASS, list(class))
}
names(CLASS) <- all.species

# once i have their sig regions i need to get the distribtuoin of edge lengths and print them off

# So now i go through the list and get distribtuions of the edge lengths for H and L regions fo the differnt species


# start by getting the four basic groups

region <- c("A.H.i", "A.L.i", "A.H.u", "A.L.u", "S.1.m", "S.2.m")
region.i <- region[1:2]
region.u <- region[3:4]
region.s <- region[5:6]

for(s in all.species){
	spec <- get(paste(s,"data", sep = "."))
	cl <- CLASS[[s]]
	for(r.i in region.i){
		A <- gsub("A.", "", gsub(".i", "", r.i))
		A <- intersect(cl$binID[cl$SINE2_MIR == A], cl$binID[cl$LINE_L2 == A])
		A <- A[order(A)]		
		assign(paste(r.i, s,sep="_"), spec[spec$binID %in% A,])
		spec1 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)])))
		spec2 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)], stand = TRUE)))
		assign(paste(r.i, s, "edge_means", sep="_"), spec1)
		assign(paste(r.i, s, "edge_means_stand", sep="_"), spec2)
	}
	
	for(r.u in region.u){
		A <- gsub("A.", "", gsub(".u", "", r.u))
		A <- union(cl$binID[cl$SINE2_MIR == A], cl$binID[cl$LINE_L2 == A])
		A <- A[order(A)]		
		assign(paste(r.u, s,sep="_"), spec[spec$binID %in% A,])	
		spec1 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)])))
		spec2 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)], stand = TRUE)))
		assign(paste(r.u, s, "edge_means", sep="_"), spec1)
		assign(paste(r.u, s, "edge_means_stand", sep="_"), spec2)
	}
	
	for(r.s in region.s){
		A <- sample(cl$binID, 400)
		A <- A[order(A)]		
		assign(paste(r.u, s,sep="_"), spec[spec$binID %in% A,])	
		spec1 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)])))
		spec2 <- colMeans(as.matrix(daisy(spec[spec$binID %in% A,2:length(spec)], stand = TRUE)))
		assign(paste(r.s, s, "edge_means", sep="_"), spec1)
		assign(paste(r.s, s, "edge_means_stand", sep="_"), spec2)
	}
}


# now that i have all the distribtuions 
# I can start to capture the edge lengths using my two different methods
# standardised and non standardised 

# now we can look in each species and see what the variance in edge lengths is like

s <- "Elephant"

thing.s <- data.frame( edges = 
						c(get(paste("A.H.i", s, "edge_means_stand", sep = "_")), 
						get(paste("A.H.u", s, "edge_means_stand", sep="_")),
						get(paste("A.L.i", s, "edge_means_stand", sep="_")),
						get(paste("A.L.u", s, "edge_means_stand", sep="_"))
						),
						
					   region = 
						c(rep("A.H.i", length(get(paste("A.H.i", s, "edge_means_stand", sep = "_")))),
						rep("A.H.u", length(get(paste("A.H.u", s, "edge_means_stand", sep = "_")))),
						rep("A.L.i", length(get(paste("A.L.i", s, "edge_means_stand", sep = "_")))),
						rep("A.L.u", length(get(paste("A.L.u", s, "edge_means_stand", sep = "_"))))
						)
						)
						


ggplot(thing.s, aes(log(edges), color = region)) + geom_density(alpha = 0.2)





s <- "Human"

thing <- data.frame( edges = 
						c(get(paste("A.H.i", s, "edge_means", sep = "_")), 
						get(paste("A.H.u", s, "edge_means", sep="_")),
						get(paste("A.L.i", s, "edge_means", sep="_")),
						get(paste("A.L.u", s, "edge_means", sep="_")),
						get(paste("S.1.m", s, "edge_means", sep="_")),
						get(paste("S.2.m", s, "edge_means", sep="_"))
						
						),
						
					   region = 
						c(rep("A.H.i", length(get(paste("A.H.i", s, "edge_means", sep = "_")))),
						rep("A.H.u", length(get(paste("A.H.u", s, "edge_means", sep = "_")))),
						rep("A.L.i", length(get(paste("A.L.i", s, "edge_means", sep = "_")))),
						rep("A.L.u", length(get(paste("A.L.u", s, "edge_means", sep = "_")))),
						rep("S.1.m", length(get(paste("S.1.m", s, "edge_means", sep = "_")))),
						rep("S.2.m", length(get(paste("S.2.m", s, "edge_means", sep = "_"))))
						)
						)
						


ggplot(thing, aes(edges, color = region)) + geom_density() + ggtitle(s)




# so far its all the same kind of distribtuion to varying degrees
# maybe we need to compare species by putting them on the same scale
# then we can have a better look at their centre of gravity 


# if we do it speices wise, we need to see that the good guys hang with each other and bad guys hang with each other
# we also need to see more variability in the blue and pirple lines
r <- "S.1.m"
spec.thing.s <- NULL
for(s in all.species){
	spec <- data.frame(edges = get(paste(r, s, "edge_means", sep = "_")), species = rep(s,length(get(paste(r, s, "edge_means", sep = "_")))))
	spec.thing.s <- rbind(spec.thing.s, spec)	
}


#ggplot(spec.thing.s, aes(edges, color = species)) + geom_density() + ggtitle(r) 


# somehow to compare the variance would be good





ggplot(spec.thing.s, aes( y = edges, x= species)) + geom_boxplot() + ggtitle(r) 





























A.H.i <- A.L.i <- A.H.u <- A.L.u <- NULL
A.H.i.s <- A.L.i.s <- A.H.u.s <- A.L.u.s <- NULL
A.H.i.g <- A.L.i.g <- A.H.u.g <- A.L.u.g <- NULL
A.H.i.gs <- A.L.i.gs <- A.H.u.gs <- A.L.u.gs <- NULL
for(s in seq(along = all.species)){
	s1 <- get(all.species[s])
	
	A.H <- 
	A.H.i <- c(A.H.i,list(as.numeric(daisy(s1[CLASS[[s]]$SINE2_MIR == "H" ,5:(length(s1)-1)], stand = F))))
	
	
	A.L.i <- c(A.L,list(as.numeric(daisy(s1[CLASS[[s]]$ancestral == "L" ,5:(length(s1)-1)], stand = F))))
	
	
	A.H.u <- c(G.H,list(as.numeric(daisy(s1[CLASS[[s]]$gene.no == "H" ,5:(length(s1)-1)], stand = F))))
	
	
	A.L.u <- c(G.L,list(as.numeric(daisy(s1[CLASS[[s]]$gene.no == "L" ,5:(length(s1)-1)], stand = F))))
	
	
	A.H.i.s <- c(A.H,list(as.numeric(daisy(s1[CLASS[[s]]$ancestral == "H" ,5:(length(s1)-1)], stand = T))))
	A.L.i.s <- c(A.L,list(as.numeric(daisy(s1[CLASS[[s]]$ancestral == "L" ,5:(length(s1)-1)], stand = T))))
	A.H.i.s <- c(G.H,list(as.numeric(daisy(s1[CLASS[[s]]$gene.no == "H" ,5:(length(s1)-1)], stand = T))))
	A.L.i.s <- c(G.L,list(as.numeric(daisy(s1[CLASS[[s]]$gene.no == "L" ,5:(length(s1)-1)], stand = T))))
	
	
	# this makes the data frame
	A.H.g <- rbind(A.H.g ,data.frame(edge_length = A.H[[s]], species = all.species[s]))
	A.L.g <- rbind(A.L.g ,data.frame(edge_length = A.L[[s]], species = all.species[s]))
	G.H.g <- rbind(G.H.g ,data.frame(edge_length = G.H[[s]], species = all.species[s]))
	G.L.g <- rbind(G.L.g ,data.frame(edge_length = G.L[[s]], species = all.species[s]))
	
	A.H.gs <- rbind(A.H.gs ,data.frame(edge_length = A.H.s[[s]], species = all.species[s]))
	A.L.gs <- rbind(A.L.gs ,data.frame(edge_length = A.L.s[[s]], species = all.species[s]))
	G.H.gs <- rbind(G.H.gs ,data.frame(edge_length = G.H.s[[s]], species = all.species[s]))
	G.L.gs <- rbind(G.L.gs ,data.frame(edge_length = G.L.s[[s]], species = all.species[s]))

	
}
names(A.H) <- all.species
names(A.L) <- all.species
names(G.H) <- all.species
names(G.L) <- all.species
















pdf(file = "plot_results/edge_length_compare/edge_compare_1_all_species-BovB_E_C.pdf", onefile = T)
ggplot(A.H.g, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("High ancestral content (PC2)")
ggplot(A.L.g, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("Low ancestral content (PC2)")
ggplot(G.H.g, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("High gene number (PC1)")
ggplot(G.L.g, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("Low gene number (PC1)")
dev.off()


pdf(file = "plot_results/edge_length_compare/edge_compare_stand_1_all_species-BovB_E_C.pdf", onefile = T)
ggplot(A.H.gs, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("High ancestral content (PC2)")
ggplot(A.L.gs, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("Low ancestral content (PC2)")
ggplot(G.H.gs, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("High gene number (PC1)")
ggplot(G.L.gs, aes(edge_length, color = species)) + geom_density(alpha = 0.2) + ggtitle("Low gene number (PC1)")
dev.off()


hist(A.H)

pdf(file = "plot_results/edge_length_compare/biplots.pdf", onefile = T)
for(s in all.species){
	p <- get(paste(s,".pca", sep =""))
	biplot(p, xlabs = rep(".", dim(p$x)[1]), main = s)
}
dev.off()


ggplot(data = A.H.g,aes(species, edge_length), geom = "violin")



p <- ggplot(G.H.g, aes(factor(species), edge_length))
p <- p + geom_violin(aes(colour = "#1268FF"), alpha = 0.3)
q <- p + geom_violin(aes(y = edge_length, colour = "#3268FF"), alpha = 0.3)

# the idea is that the distribution for the high ancestral has conserved edge lengths




















