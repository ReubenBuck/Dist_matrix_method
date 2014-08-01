
# load in genomes and selecty axis on whihc to work
# make sure they are in the sam eorder 

# merge on colnames will automate this step 
rm(list= ls())
setwd("~/Desktop/Dist_matrix_TE_div/")



#####
# Load in Tables and process
filter = .2

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

x = 3

#R_classes <- NULL
#for(x in 1:length(useful.species)){
	spec2 <- useful.species[x]
	

# load in Sbin	
S.bin.name <- paste("S_bin/", spec1, "_sbin/", spec1, "_aligning_", spec2, "_select", sep = "")
S.bin_s1_s2 <- read.table(S.bin.name,header = TRUE)

rem.un <- "yes"

mb <- 1000000
bin.size = 500000


keep.NGenes = "no"
keep.NG4s = "no"
keep.NCpGI = "no"
keep.CpGBP = "no"
keep.GC = "no"
SCALE = "yes"

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

# filtering on sbin , filter option decides on which s2 value the cutoff should be made
S.bin_s1_s2 <- S.bin_s1_s2[S.bin_s1_s2[,1] %in% unique(S.bin_s1_s2[S.bin_s1_s2[,4] > filter,1]),]


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

#  data is lodaed in now we can do the merge 
# read out both species in the correct format 
# rows are edges 
# node 1, node 2, edge wight (distance)
# s1
# replot_s1_s2
# expand grid takes evrything
# maybe combn will get the right coordinates


DistS1 <- dist(s1[5:length(s1)])
DistS2 <- dist(replot_s1_s2[2:length(replot_s1_s2)])
# expnad the gird and add the distance scores


e.grid <- expand.grid(1:length(rownames(as.matrix(DistS1))), 1:length(rownames(as.matrix(DistS1))))
# rearange the cols and get rid of any that are equal 

e.grid <- e.grid[,2:1]


com <- t(combn( dim(s1)[1], m = 2))
com1 <- data.frame(com, as.numeric(DistS1)) 
com2 <- data.frame(com, as.numeric(DistS2)) 

colnames(com2) <- colnames(com1) <- c(dim(s1)[1], "", "")
# line 1 has to equal the amount of nodes

path <- "~/Desktop/persistance/"
write.table(com1, file = paste(path, spec1, "_edges.txt", sep = ""), quote = F, row.names = F, col.names = T)
write.table(com2, file = paste(path, spec2, "_edges.txt", sep = ""), quote = F, row.names = F, col.names = T)



# write the table

# we have our list of nodes, that is good







