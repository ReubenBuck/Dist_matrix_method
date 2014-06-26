# Random TE distribution simulation
#
# Thereby we will model the idea of what it looks like when you gain something and loose something.

# how do we decide what to get rid of and how many 
# The same with adding
# The best way is to select numbers that are representative of what we are really seeing
# 
# maybe i need to look at that paper that models it as a random process and implement their same things

#
# 10 chromosomes of 200 bins

# insertion-elimination model
# a random process used to explain the distribtuion of TEs

# I guess it will be important to establish that our classification system means that TEs are arranged according to power laws. 

# we can sort of model it in our bin wise fashion
# we have adjacent bins of varing TE family density. The TE family density should be reflective of powerlaw distributions

# 

rm(list = ls())

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
SCALE = "no"

#create objects into whcih i will store the binsizes
s1name <- paste(spec1, "_AllBinCounts.txt", sep = "")
s1 <- read.table(s1name, header = TRUE)
slist <- list(s1)


for(i in seq(along=slist)){
      count <- slist[[i]]
      count <- count[count$Known >= bin.size,]
# for this we will just look at the number of elements
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

s1 <- slist[[1]]#[,!(colnames(slist[[1]]) == "Known")]

if(rem.un == "yes"){
	if(length(grep("U", s1$chr)) > 0){s1 <- s1[-(grep( "U", s1$chr)),]}
	if(length(grep("_", s1$chr)) > 0){s1 <- s1[-(grep("_", s1$chr)),]}
	if(length(grep("M", s1$chr)) > 0){s1 <- s1[-(grep("M", s1$chr)),]}
	}


g1 <- s1
plot(g1$LINE_L1, type = "l")
hist((g1$LINE_L1+g1$SINE1_7SL))
hist(g1$SINE2_MIR)
# so a power law would look like there were lots of small spaces and not many big spaces.
# the L1 density is proportional to the amount of small spaces near each other
# plot the number of elements per bin and the propotion of the bin taken up by elements
# maybe just look at the number of elements
 TE <- g1$LINE_L2
A <- TE * 250
B <- 500000 - A
B <- B/ TE

plot( (B),(TE))
plot( log(B),log(TE))

# a power relationship between the number of speces in a bin and the size of the spaces in a bin 

# if we have the rate per 500 thou bp we can known the rate of space per 500 thou bp

# what if we have a mixed TE group one that has the same distributin has some taken away and some new ones added, can do this probably later


# step 1 is to randomly plot the same distribution in two genomes
##########
#
#  Building the ancestral genome coordinates
#
##########
F1 <- data.frame(matrix(nrow = 2000, ncol = 4))
colnames(F1) <- c("chr", "start", "end", "old_TE")

chromosome <- paste("chr", 1:10, sep = "")
start <- rep(1,200)
for(i in 1:199){
	start[i+1] <- start[i+1] + (i * 1500000)
}
end <- rep(1500000,200)
for(i in 1:199){
	end[i+1] <- end[i+1] + (i * 1500000)
}
for(i in seq(length(chromosome))){
	F1[(200 * i - 199):(200 * i),1] <- chromosome[i]
	F1[(200 * i - 199):(200 * i),2] <- start
	F1[(200 * i - 199):(200 * i),3] <- end
}

# we have our F1 genome coordinates set up
# lets put in some TEs.

# so it is easy to set this up but how do we sample a random number 2000 times quickly rather than sampling 2000 random numbers

# bin the each smaple then get sample less then the lower amount for that bin and do as many are in the bin
# keep in mind there is also a normal distribution of the amount in the bin 
# i wonder if thats a product of sampling



##########
#
#  Building the new and old TE distributions
#
##########



# the ancestral part 
# a random uniform distribution
# remove a random number of TEs 
#     decided by sampling the amount of TEs in each bin over 
F1[,4] <- as.integer(runif(2000,min = min(g1$LINE_L1), max = max(g1$LINE_L1)))

newspec <- matrix(nrow = 2000, ncol = 3)
for(y in 1:3) {	
	remove.all = NULL
	for(z in 1:10){
		remove = rep(0,2000)
		for(i in 1:2000){
			remove[i] <- (sample(0:F1[i,4],1))
			}
		remove.all <- cbind(remove.all, remove)
		}
	rm <- rowMeans(remove.all) 
	newspec[,y] = F1[,4] - rm
}

s1 <- F1
s1$old_TE <- newspec[,1]
s2 <- F1
s2$old_TE <- newspec[,2]
s3 <- F1
s3$old_TE <- newspec[,3]

s1$new_TE <- as.integer(rnorm(2000,mean = mean(g1$LINE_L1), sd = sd(g1$LINE_L1)))
s2$new_TE <- as.integer(rnorm(2000,mean = mean(g1$LINE_L1), sd = sd(g1$LINE_L1)))
s3$new_TE <- as.integer(rnorm(2000,mean = mean(g1$LINE_L1), sd = sd(g1$LINE_L1)))

cor <- cor(s1[4:5], s2[4:5])
# if we look at stuff through the window os stuff that hasnt moved and through the window of stuff that has moved what do we see in regard to old and new


TE <- s1$old_TE
A <- TE * 250
B <- 500000 - A
B <- B/ TE

plot( (B),(TE))
plot( log(B),log(TE))

##########
#
#  dist matrix methods
#
##########



dist.s1 <- dist(s1[,4:5])
dist.s2 <- dist(s2[,4:5])
dist.s3 <- dist(s3[,4:5])
# have a distance matrix for each one now i have to comapre the distances across s1 to s2 to identify s1 bins that have moved.


# get bins that are diff and the same between s1 and s2 and diff between s1 and s3

# then we look to see what the overlap loos like

cor(dist.s1,dist.s2)
cor(dist.s1,dist.s3)
cor(dist.s2,dist.s3)





S1.dist.m <- as.matrix(dist.s1)
S2.dist.m <- as.matrix(dist.s2)
S3.dist.m <- as.matrix(dist.s3)
cor.s12 <- cor(S1.dist.m,S2.dist.m)
cor.diag.s12 <- rep(0, length(rownames(cor.all)))
for(i in 1:length(rownames(cor.all))){
	cor.diag.s12[i] <- cor.s12[i,i]
}
s1$s2.archi.cor.scale <- scale(cor.diag)

cor.s13 <- cor(S1.dist.m,S3.dist.m)
cor.diag.s13 <- rep(0, length(rownames(cor.all)))
for(i in 1:length(rownames(cor.all))){
	cor.diag.s13[i] <- cor.s13[i,i]
}
s1$s3.archi.cor.scale <- scale(cor.diag.s13)

# seperate genome into contigous sections





All <- data.frame(chr = s1$chr
				  start = s1$start
				  end = s1$end
				  s2.archi.cor.scale = s1$s2.archi.cor.scale, 
				  s3.archi.cor.scale = s1$s3.archi.cor.scale,
				  s1.new.TE = s1$new_TE,
				  s2.new.TE = s2$new_TE,
				  s3.new.TE = s3$new_TE,
				  s1.old.TE = s1$old_TE,
				  s2.old.TE = s2$old_TE,
				  s3.old.TE = s3$old_TE)





s1$group = 0
group = 0
for(i in seq(dim(s1)[1])){
	if(all(s1$start[i] != s1$end[i-1]+1)){
		group = group+1
	}
	s1[i,'group'] <- group
}







# can probably loop the wiggle tree for every single repeat group  


# Now classify bins based the local variance of the scaled_r score

# get a wiggle thing for each and produce a list of wiggles


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


rows <- 500:600
plot(Wiggle[rows,'mean'], type = "l")
lines(Wiggle[rows, 'bottom'], col = 2)
lines(Wiggle[rows, 'top'], col = 2)
abline(h = 0 , col = 3)

# we need to run this wiggle loop on a few things to get the highs nd lows we need



