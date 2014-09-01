
library(ggplot2)
library(gplots)
library(pdist)
library(cluster)

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
bin.limit <- 0


#####
# Load in Tables and process

spec1 <- "Human"
all.species <- c("Human", "Mouse", "Rat", "Dog", "Bovine", "Horse", "Platypus", "Elephant")
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
	colnames(spec)[colnames(spec) == "S1_bin"] <- "binID"
	assign(useful.species[i],spec)	
}


# lets run it on the sampled data set


assign(spec1, s1[s1$binID %in% IDs,c(2,5:length(s1))])
s1.samp <- s1[s1$binID %in% IDs,]



Group = 0
group = 0
for(i in seq(dim(s1.samp)[1])){
	if(all(s1.samp$start[i] != s1.samp$end[i-1]+1)){
		group = group+1
	}
	Group[i] <- group
}

pca <- NULL
pca$x <- Human[,2:length(Human)]
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
class$binID <- s1.samp$binID



# so we have H,M,L for each TE. For the conservation analysis we need to consider the variation in the base distribtuion. 
# Or we use two together somehow. 
# We have edge lengths on this graph and we compare them 


# a H in both or a H in at least one

# or just a high in MIR


A.H.u <- union(class[class[c("LINE_L2")] == "H","binID"], class[class[c("SINE2_MIR")] == "H","binID"])
A.H.u <- A.H.u[order(A.H.u)]
A.H.i <- intersect(class[class[c("LINE_L2")] == "H","binID"], class[class[c("SINE2_MIR")] == "H","binID"])
A.H.i <- A.H.i[order(A.H.i)]
A.L.u <- union(class[class[c("LINE_L2")] == "L","binID"], class[class[c("SINE2_MIR")] == "L","binID"])
A.L.u <- A.L.u[order(A.L.u)]
A.L.i <- intersect(class[class[c("LINE_L2")] == "L","binID"], class[class[c("SINE2_MIR")] == "L","binID"])
A.L.i <- A.L.i[order(A.L.i)]

R.S.1 <- sample(IDs,length(A.H.u))
R.S.1 <- R.S.1[order(R.S.1)]
R.S.2 <- sample(IDs,length(A.H.i))
R.S.2 <- R.S.2[order(R.S.2)]
R.S.3 <- sample(IDs,length(A.L.u))
R.S.3 <- R.S.3[order(R.S.3)]
R.S.4 <- sample(IDs,length(A.L.i))
R.S.4 <- R.S.4[order(R.S.4)]

region <- c("A.H.i", "A.L.i", "A.H.u", "A.L.u", "R.S.1", "R.S.2", "R.S.3", "R.S.4")


for(s in c(useful.species, spec1)){
	spec <- get(s)
	spec.pca <- prcomp(spec[,2:length(spec)])	
	assign(paste(s, ".pca", sep = ""), spec.pca)
}






# get the ancestral distribution
for( s in c(useful.species,spec1)){
	assign(paste("ancestral.", s, sep = ""), get(s)[,c("binID","SINE2_MIR", "LINE_L2")])
}

# lets look at distnce moved rather than change in node length
ancestral.movement <- data.frame(bin_no = IDs)
for(s in c(useful.species,spec1)){
	spec <- get(paste("ancestral.", s, sep = ""))
	Dist <- (spec$SINE2_MIR - ancestral.Human$SINE2_MIR) + (spec$LINE_L2 - ancestral.Human$LINE_L2)
	ancestral.movement <- cbind(ancestral.movement, Dist)
}
colnames(ancestral.movement) <- c("binID", c(useful.species,spec1))

plot(hclust(dist(t(as.matrix(ancestral.movement)[ancestral.movement$binID %in% A.L.i,1:(length(useful.species)+1)+1]))))


# with the ancestral movement lets make it easier and put it into regions so it can be better analysed

for(r in region){
	reg <- get(r)
	name <- paste(r, "_ancestral_movement", sep ="")
	AM <- ancestral.movement[ancestral.movement$binID %in% reg,]
	assign(name, AM)
}



pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/dendrogram/ancestral.pdf", sep = ""), onefile = T)
plot(hclust(dist(t(as.matrix(A.H.i_ancestral_movement)[,1:(length(useful.species)+1)+1]))), main = "ancestral",sub = "", xlab = "High MIR/L2 region (intercept)")

plot(hclust(dist(t(as.matrix(A.L.i_ancestral_movement)[,1:(length(useful.species)+1)+1]))),  main = "ancestral",sub = "", xlab = "Low MIR/L2 region (intercept)")

plot(hclust(dist(t(as.matrix(R.S.1_ancestral_movement)[,1:(length(useful.species)+1)+1]))),  main = "ancestral",sub = "", xlab = "Random sample region 1")

plot(hclust(dist(t(as.matrix(A.H.u_ancestral_movement)[,1:(length(useful.species)+1)+1]))), main = "ancestral",sub = "", xlab = "High MIR/L2 region (union)")

plot(hclust(dist(t(as.matrix(A.L.u_ancestral_movement)[,1:(length(useful.species)+1)+1]))), main = "ancestral",sub = "", xlab = "Low MIR/L2 region (union)")

plot(hclust(dist(t(as.matrix(R.S.2_ancestral_movement)[,1:(length(useful.species)+1)+1]))),  main = "ancestral",sub = "", xlab = "Random sample region 2")

plot(hclust(dist(t(as.matrix(R.S.3_ancestral_movement)[,1:(length(useful.species)+1)+1]))),  main = "ancestral",sub = "", xlab = "Random sample region 3")

plot(hclust(dist(t(as.matrix(R.S.4_ancestral_movement)[,1:(length(useful.species)+1)+1]))),  main = "ancestral",sub = "", xlab = "Random sample region 4")


dev.off()







# By captuing the actual movement we might be able to accuratly adjust the edge lengths to refelect changes in the unknown area

# In other words controlling for uninteresting variation. 



# so now that we have the actual movement we need to get regional movement for the rest of the TEs















# step 1 is to get the new genomes sorted



for( s in c(useful.species,spec1)){
	spec <- get(s)
	spec <- spec[,!(colnames(spec) %in% c("SINE2_MIR", "LINE_L2"))]
	assign(paste("new.", s, sep = ""),spec)
}



# so here we are calculating the distances

# So calculate all standardised distances species wise
# then we can subset the distance matricies as needed 



for(s in c(useful.species,spec1)){
	spec <- get(paste(s,".pca", sep = ""))
	# we choose 17 becasue 17 is the lowest dimensionality
	spec <- spec$x[,1:17]
	Dist <- as.matrix(daisy(spec, stand = FALSE))
	rownames(Dist) <- colnames(Dist) <- IDs
	assign(paste(s, "Dists", sep = "."), Dist)
}




for(r in region){
	reg <- get(r)
	edge.df <- data.frame(binID = reg)
	for(s in c(useful.species,spec1)){
		spec <- get(paste(s,"Dists", sep = "."))
		Dist <- colMeans(spec[rownames(spec) %in% reg,colnames(spec) %in% reg])
		edge.df <- cbind(edge.df, Dist)
	}
	colnames(edge.df) <- c("binID", c(useful.species,spec1))
#	edge.df[2:dim(edge.df)[2]] <- scale(edge.df[2:dim(edge.df)[2]])
	assign(paste(r, "_new_edges", sep = ""), edge.df)
}


# Here we get the movement form human of edge lengths

for(r in region){
	reg <- get(r)
	unknown.movement <- data.frame(binID = reg)
	Human.dist <- get(paste(r, "_new_edges", sep = ""))[,"Human"]
	unknown.Human <- ancestral.Human[ancestral.Human$binID %in% reg,]
	for(s in c(useful.species,spec1)){
		dist <- get(paste(r, "_new_edges", sep = ""))[,s]
		Dist <- (dist - Human.dist)
		unknown.movement <- cbind(unknown.movement, Dist)
	}
	colnames(unknown.movement) <- c("binID", c(useful.species,spec1))
	assign(paste(r, "_unknown_movement", sep = ""), unknown.movement)
}


pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/dendrogram/new_move.pdf", sep = ""), onefile = T)
plot(hclust(dist(t(as.matrix(A.H.i_unknown_movement)[,1:(length(useful.species)+1)+1]))), main = "edge_length",sub = "", xlab = "High MIR/L2 region (intercept)")
plot(hclust(dist(t(as.matrix(A.H.u_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "High MIR/L2 region (union)")
plot(hclust(dist(t(as.matrix(A.L.i_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Low MIR/L2 region (intercept)")
plot(hclust(dist(t(as.matrix(A.L.u_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Low MIR/L2 region (union)")
plot(hclust(dist(t(as.matrix(R.S.1_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Random sample region 1")
plot(hclust(dist(t(as.matrix(R.S.2_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Random sample region 2")
plot(hclust(dist(t(as.matrix(R.S.3_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Random sample region 3")
plot(hclust(dist(t(as.matrix(R.S.4_unknown_movement)[,1:(length(useful.species)+1)+1]))),main = "edge_length",sub = "", xlab = "Random sample region 4")
dev.off()


# The unknonw movement is capturing differences from human 
# basicly we are not capturing theright information because the differences from human can only be positive or negative and how far
# this is tough because things that have a positive movement by the same amount can actually be quite different  
# therefore we need to cluster on two values rather than one




# next step is to introduce some of the ploting techniques i used earlier, just to see what it actually looks like
# sort of what i had before? maybe human is positioned diferently and the mouse now has more red






# make all the pca files

# make all the every files
for(r in region){
	every <- data.frame(binID <- IDs, matrix(data = NA,ncol = length(c(useful.species,spec1)), nrow = length(IDs)))
	
	diff <- get(paste(r,"_unknown_movement", sep =""))
	every[every[,1] %in% diff[,1],2:dim(every)[2]] <- diff[,2:dim(every)[2]]
	colnames(every)[2:dim(every)[2]] <- c(useful.species,spec1)
	assign(paste(r, "_every1", sep = ""), every)
}

for(r in region){
	every <- data.frame(binID <- IDs, matrix(data = NA,ncol = length(c(useful.species,spec1)), nrow = length(IDs)))
	
	diff <- get(paste(r,"_ancestral_movement", sep =""))
	every[every[,1] %in% diff[,1],2:dim(every)[2]] <- diff[,2:dim(every)[2]]
	colnames(every)[2:dim(every)[2]] <- c(useful.species,spec1)
	assign(paste(r, "_every2", sep = ""), every)
}




for(s in c(useful.species,spec1)){
	spec <- get(paste(s,".pca", sep =""))
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/test2/biplot/", "biplot_", s,".pdf", sep = ""))
	biplot(spec, xlabs = rep(".", length(IDs)), main = s)
	dev.off()
}


for(r in region){
	for(s in c(useful.species,spec1)){
		cols <- get(paste(r,"_every1", sep = ""))
		q <-  qplot(Human.pca$x[,1], Human.pca$x[,2], color = cols[,s]) 
		q <- q + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) 
		q <- q + labs(title = paste("Human", s,r, sep = "_"))
		ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/test1/",r,"/", "Human.pca_spec_edge_change_", r,"_", s,".pdf", sep = ""), plot = q)
		spec <- get(paste(s,".pca", sep =""))
		q <- qplot(spec$x[,1], spec$x[,2], color =  cols[,s]) 
		q <- q + scale_colour_gradient2(low=("green"), high=("red"), limits = c(-5,5)) 
		q <- q + labs(title = paste(s, r, sep = "_"))
		ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/test1/",r,"/", "spec_edge_change_", r, "_", s,".pdf", sep = ""))
	}
}



for(r in region){
	for(s in c(useful.species,spec1)){
		cols <- get(paste(r,"_every2", sep = ""))
		q <-  qplot(Human.pca$x[,1], Human.pca$x[,2], color = cols[,s]) 
		q <- q + scale_colour_gradient2(low=("red"), high=("green"), limits = c(-5,5)) 
		q <- q + labs(title = paste("Human", s,r, sep = "_"))
		ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/test2/",r,"/", "Human.pca_spec_ancestral_change_", r,"_", s,".pdf", sep = ""), plot = q)
		spec <- get(paste(s,".pca", sep =""))
		q <- qplot(spec$x[,1], spec$x[,2], color =  cols[,s]) 
		q <- q + scale_colour_gradient2(low=("red"), high=("green"), limits = c(-5,5)) 
		q <- q + labs(title = paste(s, r, sep = "_"))
		ggsave(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/test2/",r,"/", "spec_ancestral_change_", r, "_", s,".pdf", sep = ""))
	}
}




# scale area wont work because you cant use it to plot divergence
















# keep in mind the ordering when using the %in% tool
# basicly we can't just col bind on it randomly


# do every combination and breate a dist matrix out of it so i can cluster on it
# do the MIR L2 stuff and combine it with the other stuff for each region

BIG.dist <- matrix(nrow = length(c(useful.species,spec1)), ncol = length(c(useful.species,spec1)))
colnames(BIG.dist) <- rownames(BIG.dist) <- c(useful.species,spec1)
for(r in region){
	for(s1 in c(useful.species,spec1)){
		for(s2 in c(useful.species,spec1)){
			spec1.a <- get(paste(r, "_ancestral_movement", sep = ""))[s1]
			spec1.n <- get(paste(r, "_unknown_movement", sep =""))[s1]
			spec1.m <- as.matrix(cbind(spec1.a, spec1.n))
			spec2.a <- get(paste(r, "_ancestral_movement", sep = ""))[s2]
			spec2.n <- get(paste(r, "_unknown_movement", sep =""))[s2]
			spec2.m <- as.matrix(cbind(spec2.a, spec2.n))

			BIG.dist[s1,s2] <- norm(spec1.m-spec2.m, type="F")
			
		}
	}
	assign(paste(r,"_distance", sep = ""), as.dist(BIG.dist))
}

pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/dendrogram/merged.pdf", sep = ""), onefile = T)
plot(hclust(A.H.i_distance), main = "merge",sub = "", xlab = "High MIR/L2 region (intercept)")
plot(hclust(A.L.i_distance), main = "merge",sub = "", xlab = "Low MIR/L2 region (intercept)")
plot(hclust(R.S.1_distance), main = "merge",sub = "", xlab = "Random sample region 1")
plot(hclust(A.H.u_distance), main = "merge",sub = "", xlab = "High MIR/L2 region (union)")
plot(hclust(A.L.u_distance), main = "merge",sub = "", xlab = "Low MIR/L2 region (union)")
plot(hclust(R.S.2_distance), main = "merge",sub = "", xlab = "Random sample region 2")
plot(hclust(R.S.3_distance), main = "merge",sub = "", xlab = "Random sample region 3")
plot(hclust(R.S.4_distance), main = "merge",sub = "", xlab = "Random sample region 4")


dev.off()




# maybe as a measure of variation

sum(A.H.i_distance)
sum(A.H.u_distance)
sum(A.L.i_distance)
sum(A.L.u_distance)





# now all i need is the individual score 
# a way of scoring them so i can see the movement 





# how does the distance from human change over each area


# this may be better plotted as a contour map

pdf(file = "~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/2d_dist/div_from_human.pdf", onefile = T)
for(s in useful.species){
	
	plot(
		 c(A.H.i_ancestral_movement[,s], A.L.i_ancestral_movement[,s], 
		 R.S.2_ancestral_movement[,s], R.S.4_ancestral_movement[,s]),
		 
		 c(A.H.i_unknown_movement[,s], A.L.i_unknown_movement[,s],
		 R.S.2_unknown_movement[,s], R.S.4_unknown_movement[,s]),
		 
		 col = c(rep(1,length(A.H.i)), rep(2, length(A.L.i)),
		 rep(3,length(R.S.2)), rep(4, length(R.S.4))),
		 
		 pch = 16, cex= .7,
		 
		 main = paste(s, "Divergence from Human"),
		 xlab = "ancestral TE content (MIR/L2)",
		 ylab = "mean edge length",
		 xlim = c(-4,4),
		 ylim = c(-4,4)
		 )
	legend("topright" , c("A.H.i", "A.L.i", "R.S.2", "R.S.4"), fill=1:4)
	}
dev.off()



for(s in useful.species){
	thing <- data.frame(ancestral = c(A.H.i_ancestral_movement[,s], A.L.u_ancestral_movement[,s],
		 R.S.2_ancestral_movement[,s], R.S.4_ancestral_movement[,s]), 
	     edge_length = c(A.H.i_unknown_movement[,s], A.L.u_unknown_movement[,s],
		 R.S.2_unknown_movement[,s], R.S.4_unknown_movement[,s]),
		 region = c(rep("A.H.i",length(A.H.i)), rep("A.L.u", length(A.L.u)),
		 rep("R.S.2",length(R.S.2)), rep("R.S.4", length(R.S.4)))
		 )
	d <- ggplot(thing, aes(ancestral, edge_length, color = region, linetype = region))+ geom_density2d()
	d
	d <- d + title(main = paste(s, "Divergence from Human"))
	ggsave(filename = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/2d_dist/", s, "_div_from_human.pdf", sep = ""))
	d <- NULL

}




# Are there any sort of assocaition statisitics we could use to determine what acould be causing the differences

# Once we have a target we can then find a way of confirming it in the real data


for(s in useful.species){
	
	thing <- data.frame( 
	     edge_length = c(A.H.i_unknown_movement[,s], A.L.i_unknown_movement[,s], A.H.u_unknown_movement[,s], A.L.u_unknown_movement[,s],
		 R.S.1_unknown_movement[,s], R.S.2_unknown_movement[,s],R.S.3_unknown_movement[,s],R.S.4_unknown_movement[,s]),
		 region = c(rep("High MIR/L2 intercept",length(A.H.i)), rep("Low MIR/L2 intercept", length(A.L.i)), rep("High MIR/L2 union",length(A.H.u)), rep("Low MIR/L2 union", length(A.L.u)),
		 rep("Random sample 1",length(R.S.1)), rep("Random sample 2",length(R.S.2)),rep("Random sample 3",length(R.S.3)),rep("Random sample 4",length(R.S.4)))
		 )
	ggplot(thing, aes( edge_length, color = region), xlab = "divergence in mean edge length from human") + geom_density() + ggtitle(s) + xlab("mean edge length difference from human")
	ggsave(filename = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/2d_dist/", s, ".pdf", sep=""))
}





# convert distances to positive lengths


for(s in useful.species){
	
	thing <- data.frame( 
	     edge_length = c(A.H.i_unknown_movement[,s], A.L.i_unknown_movement[,s], A.H.u_unknown_movement[,s], A.L.u_unknown_movement[,s],
		 R.S.1_unknown_movement[,s], R.S.2_unknown_movement[,s],R.S.3_unknown_movement[,s],R.S.4_unknown_movement[,s]),
		 region = c(rep("High MIR/L2 intercept",length(A.H.i)), rep("Low MIR/L2 intercept", length(A.L.i)), rep("High MIR/L2 union",length(A.H.u)), rep("Low MIR/L2 union", length(A.L.u)),
		 rep("Random sample 1",length(R.S.1)), rep("Random sample 2",length(R.S.2)),rep("Random sample 3",length(R.S.3)),rep("Random sample 4",length(R.S.4)))
		 )
		 
	thing[,1] <- sqrt(thing[,1]^2)
	
	ggplot(thing, aes((edge_length), color = region), xlab = "divergence in mean edge length from human") + geom_density() + ggtitle(s) + xlab("mean edge length difference from human")
	ggsave(filename = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/2d_dist/pos.dist.", s, ".pdf", sep=""))
}






for(s in useful.species){
	
	thing <- data.frame( 
	     edge_length = c(A.H.i_unknown_movement[,s], A.L.i_unknown_movement[,s], A.H.u_unknown_movement[,s], A.L.u_unknown_movement[,s])		 ,
		 region = c(rep("High MIR/L2 intercept",length(A.H.i)), rep("Low MIR/L2 intercept", length(A.L.i)), rep("High MIR/L2 union",length(A.H.u)), rep("Low MIR/L2 union", length(A.L.u)))
		 )
		 
	thing[,1] <- sqrt(thing[,1]^2)
	
	ggplot(thing, aes((edge_length), color = region), xlab = "divergence in mean edge length from human") + geom_density() + ggtitle(s) + xlab("mean edge length difference from human")
	ggsave(filename = paste("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/2d_dist/pos.dist.", s, ".pdf", sep=""))
}








A <- matrix(data=c(1,1))
B <- matrix(data=c(2,2))

D <- (A %*% t(A)) - (B %*% t(B))

 C <- sqrt(sum(((A %*% t(A))^2 - (B %*% t(B))^2)))

C <- sqrt(sum(((A - B)^2)))

D <- rbind(t(A), t(B))
dist(D)

C <- norm(A-B, type="F")
C
# a possible way to capture the distance between various matricies
# It will then be possible to cluster 






set.seed(111)
A = matrix(rnorm(9), nr=3)
set.seed(222)
B = matrix(rnorm(9), nr=3)
set.seed(333)
C = matrix(rnorm(9), nr=3)


m1 = diag(A[1, ])
m2 = diag(A[2, ])
b1 = B[1, ]
b2 = B[2, ]
c1 = C[1, ]
c2 = C[2, ]
distance = mean(m1 %*% ( (diag(b1)-diag(b2)) %*% (diag(c1)-diag(c2)) %*% m2))

ff<-function(i,j) mean(diag(A[i,]) %*% ( (diag(B[i,])-diag(B[j,])) %*% (diag(C[i,])-diag(C[j,])) %*% diag(A[j,])))

distMatrix<-outer(1:3,1:3,Vectorize(ff))





