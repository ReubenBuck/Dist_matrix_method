# Random TE distribution simulation
#

# old TE model 
# random insertion, random deletion, divergence, random deletion

# new TE model
# 								 divergence, random insertion, random deletion

# mix TE model
# random insertion, random deletion, divergence, random insertion, random deletion 

# min and max number will be 0:100 to keep it simple for the simulation

# select sampling numbers based on evolutionary time
# also run the loop a number of times




# would be interesting to try a genome with coreations 1, 0, -1




# step 1 is to randomly plot the same distribution in two genomes

##########
#
#  Building the ancestral genome coordinates
#
##########
rm(list = ls())

olist <- NULL
for(q in 1:100){

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




##########
#
#  Building the new and old TE distributions
#
##########



# the ancestral part 
# a random uniform distribution
# remove a random number of TEs 
#     decided by sampling the amount of TEs in each bin over 
#old <- sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE)
#mix <- sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE)

#s1 <- data.frame(old = old - sample(0:100, 2000, replace = TRUE),
#			   new = sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE), 
#			   mix = mix + sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE))
#s2 <- data.frame(old = old - sample(0:100, 2000, replace = TRUE),
#			   new = sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE), 
#			   mix = mix + sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE))
#s3 <- data.frame(old = old - sample(0:100, 2000, replace = TRUE),
#			   new = sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE),
#			   mix = mix + sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE))
#s4 <- data.frame(old = old - sample(0:100, 2000, replace = TRUE), 
#			   new = sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE), 
#			   mix = mix + sample(0:100, 2000, replace = TRUE) - sample(0:100, 2000, replace = TRUE))


# maybe a new approach to building my species genomes

#option two where we measure three correations 1,0,-1
# the 0 one is close to 0 but not ecatly 0
# the ancestral part 
# a random uniform distribution
# remove a random number of TEs 
#     decided by sampling the amount of TEs in each bin over 
old <- rnorm(mean=100, sd = 30, 2000)
#mix <- old * -1

s1 <- data.frame(old = old ,
			   new = rnorm(mean=100, sd = 30, 2000)
			   #mix = rnorm(mean=100, sd = 30, 2000) 
			   )
s2 <- data.frame(old = old ,
			   new = rnorm(mean=100, sd = 30, 2000) 
			   #mix = rnorm(mean=100, sd = 30, 2000)
			   )
s3 <- data.frame(old = old ,
			   new = rnorm(mean=100, sd = 30, 2000)
			   #mix = rnorm(mean=100, sd = 30, 2000)
			   )
s4 <- data.frame(old = old ,
			   new = rnorm(mean=100, sd = 30, 2000)
			   #mix = rnorm(mean=100, sd = 30, 2000)
			   )

# maybe a new approach to building my species genomes



			   
##########
#
#  dist matrix methods
#
##########

S1.dist.m <- as.matrix(dist(s1))
S2.dist.m <- as.matrix(dist(s2))
S3.dist.m <- as.matrix(dist(s3))
S4.dist.m <- as.matrix(dist(s4))

cor.s12 <- cor(S1.dist.m,S2.dist.m)
cor.diag.s12 <- rep(0, length(rownames(cor.s12)))
for(i in 1:length(rownames(cor.s12))){
	cor.diag.s12[i] <- cor.s12[i,i]
}
s2.archi.cor.scale <- scale(cor.diag.s12)

cor.s13 <- cor(S1.dist.m,S3.dist.m)
cor.diag.s13 <- rep(0, length(rownames(cor.s13)))
for(i in 1:length(rownames(cor.s13))){
	cor.diag.s13[i] <- cor.s13[i,i]
}
s3.archi.cor.scale <- scale(cor.diag.s13)

cor.s14 <- cor(S1.dist.m,S4.dist.m)
cor.diag.s14 <- rep(0, length(rownames(cor.s14)))
for(i in 1:length(rownames(cor.s14))){
	cor.diag.s14[i] <- cor.s14[i,i]
}
s4.archi.cor.scale <- scale(cor.diag.s14)



##########
#
#  sig region finder
#
##########



All <- data.frame(s2.archi.cor.scale = s2.archi.cor.scale, 
				s3.archi.cor.scale = s3.archi.cor.scale,
				s4.archi.cor.scale = s4.archi.cor.scale,
				s1.new.TE = scale(s1$new),
				s1.old.TE = scale(s1$old)
				#s1.mix.TE = scale(s1$mix)
				)


group = 0
for(i in seq(dim(All)[1])){
	if(all(F1$start[i] != F1$end[i-1]+1)){
		group = group+1
	}
}

jiggle.list <- (NULL)
for(y in seq(length(All))){
	jiggler <- All[,y]
	for(z in seq(group)){
		scale_r <- jiggler 
		zero <- rep(0, length(scale_r))
		wiggle <- data.frame(scale_r_score = zero, mean = zero , bottom = zero , top = zero)
		Jiggle <- NULL
		for(i in seq(along=scale_r)){	
			if((i == length(scale_r)) & (i == 1)){
				a <- mean(scale_r[i])
				s <- sd(jiggler)
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
		Jiggle <- rbind(Jiggle,wiggle )
	}
	jiggle.list <- c(jiggle.list, list(Jiggle))
}
names(jiggle.list) = names(All)

# we need to run this wiggle loop on a few things to get the highs nd lows we need

HL <- matrix(nrow = 2000, ncol= length(jiggle.list))
colnames(HL) <- names(jiggle.list)
for(i in seq(length(jiggle.list))){
	HL[,i] <- rep("M", 2000)
	HL[jiggle.list[[i]][,"top"] < 0  ,i] <- "L"
	HL[jiggle.list[[i]]$bottom > 0 ,i] <- "H"
}
HL <- as.data.frame(HL)
HL$binID <- 1:2000

##########
#
#  sig region compare
#
##########


first <- c("H", "H", "L", "L")
second <- c("H", "L", "H", "L")

# cycle through the TE based on S1
TE.g <- c("s1.old.TE", "s1.new.TE")


#speceis <- c("s2", "s3")

# the amount of overlap between both species 

overlap.matrix <- data.frame(matrix(nrow = length(first), ncol = length(TE.g)))
colnames(overlap.matrix) <- TE.g
rownames(overlap.matrix) <- paste(first,second,sep = "")

for(z in seq(TE.g)){
	# get the TE type
	for(i in seq(first)){
		s2ar <- HL[(HL[,"s2.archi.cor.scale"] == first[i]) & (HL[,TE.g[z]] == second[i]),"binID"]
		s3ar <- HL[(HL[,"s3.archi.cor.scale"] == first[i]) & (HL[,TE.g[z]] == second[i]),"binID"]
		s4ar <- HL[(HL[,"s4.archi.cor.scale"] == first[i]) & (HL[,TE.g[z]] == second[i]),"binID"]
		# then calculate the proportion that are the same
		overlap.matrix[i,z] <- sum(length(intersect(s2ar,s3ar)), length(intersect(s2ar,s4ar)), length(intersect(s3ar,s4ar)))/sum(length(s2ar), length(s3ar), length(s4ar))
	}
}

# also need to get the correct total proportion over lap
total.overlap <- rep(0, length(first))
	for(i in seq(first)){
		s2ar <- HL[(HL[,"s2.archi.cor.scale"] == first[i]),"binID"]
		s3ar <- HL[(HL[,"s3.archi.cor.scale"] == first[i]),"binID"]
		s4ar <- HL[(HL[,"s4.archi.cor.scale"] == first[i]),"binID"]
		total.overlap[i] <- sum(length(intersect(s2ar,s3ar)), length(intersect(s2ar,s4ar)), length(intersect(s3ar,s4ar)))/sum(length(s2ar), length(s3ar), length(s4ar))
}
overlap.matrix$total.overlap <- total.overlap

olist <- c(olist,list(overlap.matrix))
}

# save the olist Robject
save(olist, file = "~/Desktop/Dist_matrix_TE_div/overlap.sim.list/old_new_rnorm_TEdist")







#it has worked and i can start looking at the properties of these two TE families according to random normal distributions 

# if I do it again it would be good to get a second list that has cor matrix for each situation to help with the analysis

# in order to understand what is happening and what we are looking at we need to model the process under randomness




color <- rep(1,2000)
color[HL[,"s4.archi.cor.scale"] == "L"] <- 2
color[HL[,"s4.archi.cor.scale"] == "H"] <- 3

color2 <- rep(1,2000)
color2[HL[,"s1.new.TE"] == "L"] <- 2
color2[HL[,"s1.new.TE"] == "H"] <- 3

color3 <- rep(1,2000)
color3[HL[,"s1.old.TE"] == "L"] <- 2
color3[HL[,"s1.old.TE"] == "H"] <- 3


plot(s1[(color == 2 ) & (color2 == 2 ),], col = 1, main = "new_TE", ylim = c(min(s1$new), max(s1$new)), xlim = c(min(s1$old), max(s1$old)))
points(s1[(color == 2 ) & (color2 == 3 ),], col = 2 )
points(s1[(color == 3 ) & (color2 == 2 ),], col = 3)
points(s1[(color == 3 ) & (color2 == 3 ),], col = 4)

legend("topleft", c("LL", "LH", "HL", "HH"), fill = c(1,2,3,4))

abline(h = mean(s1$new), v = mean(s1$old))



plot(s1[(color == 2 ) & (color3 == 2 ),], col = 1, main = "old_TE", ylim = c(min(s1$new), max(s1$new)), xlim = c(min(s1$old), max(s1$old)))
points(s1[(color == 2 ) & (color3 == 3 ),], col = 2 )
points(s1[(color == 3 ) & (color3 == 2 ),], col = 3)
points(s1[(color == 3 ) & (color3 == 3 ),], col = 4)

legend("topleft", c("LL", "LH", "HL", "HH"), fill = c(1,2,3,4))

abline(h = mean(s1$new), v = mean(s1$old))

plot(s1[(color == 2 )| color == 3, ])
