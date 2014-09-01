# we have ranked lists 
# we go down until they are duplicated 
# maybe on sample data fist

library(cluster)

rm(list = ls())

# data
x <- 9
# preservation
pr <- 0

S1 <- data.frame(rnorm(x))
C <- rnorm(x)
S2 <- data.frame(S1 - runif(x, min = -0.5,max = .5) , C = sqrt(C^2))

# correlations
cor(S1,S2)

# calculate weights in random data 
d.s1 <- daisy(S1,stand=F)
d.s2 <- daisy(S2,stand=F)
com <- t(combn( dim(S1)[1], m = 2))
colnames(com) <- c("node_1", "node_2")
com1 <- data.frame(com, weights = as.numeric(d.s1)) 
com2 <- data.frame(com, weights = as.numeric(d.s2)) 

# workout how to draw the lines that we want
s1 <- com1[order(com1$weights),]
s2 <- com2[order(com1$weights),]
both <- data.frame(s1[,1:2], s1_weights = s1$weights, s2_weights = s2$weights)
rownames(both) <- 1:dim(both)[1]
both$d.score <- (both$s2_weights - both$s1_weights) 
pre.score <- data.frame(row = c(both$node_1, both$node_2), col = rep(1:dim(both)[1],2), score = both$d.score)
pre.score$length <- ((pre.score$col - 1) * x)  + pre.score$row

# get score matrix
score <- rep(NA, dim(both)[1] * x)
score[pre.score$length] <- pre.score$score
score <- matrix(data = score,  nrow = x, ncol = dim(both)[1], byrow =F)
score <- score - pr



# identify where two nodes become one
score[score < 0] <- 0
cols <- colSums(score, na.rm = T)
two2one <- which(cols == 0)
score.nm <- score

# merge score matrix
# this bit needs to run backwards
if(length(two2one[two2one<length(cols)]) > 0){
# the row.merge step has to be done on the original matrix and saved
row.merge <- NULL
for(i in two2one[two2one<length(cols)]){
	row.merge <- rbind(row.merge,matrix(which(score[,i] == 0),nrow = 1))
}

# this part still isnt working exactly the way i want it to but only in a minor way
# can be fixed but will increase runnung time
	for(o in 1:length(two2one)){
		for(i in rev(1:length(two2one[two2one<length(cols)]))){
			slice <- matrix(score[row.merge[i,],(two2one[i]+1):length(cols)], nrow=2)
			line1 <- slice[1,!(is.na(slice[1,]))]
			line2 <- slice[2,!(is.na(slice[2,]))]
			rep1 <- which(!(is.na(slice[1,])))
			rep2 <- which(!(is.na(slice[2,])))
			slice[2,rep1] <- line1
			slice[1,rep2] <- line2
			score[row.merge[i,],(two2one[i]+1):length(cols)] <- slice
		}
	}
}
# check for repeated zero cols and remove them
# which ones are the zeros and are they exactly the same
if(length(two2one) > 1){
	List <- apply(score[,two2one] == 0,2,which)
	dup <- which(duplicated(List))
	dup <- two2one[dup]
	if(length(dup) > 0){
		score <- score[,-dup]
		cols <- colSums(score, na.rm = T)
		two2one <- which(cols == 0)
	}			
}


# calculate node numbers
z <- c(1,two2one,length(cols))
no <- NULL
for(i in 1:(length(z)-1)){
	no <- c(no,rep(i,z[i+1]-z[i]))
}
if(cols[length(cols)]==0){
	no <- c(no,no[length(no)]+1) -1
}else{
	no <- c(no,no[length(no)]) -1
}


node.no <- rep(x, length(no)) - no


# calculate theoretical edge numbers
edge.cal <- function(arg){(arg^2 - arg)/2}
end.edge.co <- edge.cal(node.no)


# old fashioned edge count
edge <- 1:length(cols)
edge.no <- edge - no

# end the score matrix here and then get the scores by counting everything after zero

score.end <- score[,1:which(edge.no == end.edge.co)]
for(i in 1:x){
	if(any(score.end[i,] == 0, na.rm = T)){
	zero <- max(which(score.end[i,] == 0))
	score.end[i,1:(zero-1)] <- NA 
	}
}

A <- rowSums(score.end, na.rm=T)








plot(S2, pch = as.character(1:x), ylim=c(0,max(S2[,2])))
points(data.frame(S1,0), col = 2, pch = as.character(1:x))
plot(A,sqrt(S2[,2])^2)

plot(rowSums(score.nm, na.rm = T), S2[,2])

library(ggplot2)


qplot(S2[,1], S2[,2], color = A)
qplot(S2[,1], S2[,2], color = rowSums(score.nm, na.rm = T))

# now we can do this 
# for now we only care about the ones after the last zero

# write up something that can read the score matrix

# getting the last zero
L0 <- min(which(node.no == min(node.no)))
end <- (end.edge.co[L0] - sum(is.na(score[,L0]))) + L0

score.end <- score[,1:end]

# for each row get the last zero, how many last zeros are there? 



# everyhting before the last zero of its row will ne NA
# maybe go with the simple version , simple version concentrates on preservetion of smaller edges and builds up progresivly 









# how to calculate edge numbers so we know when to stop
# the problem I can see if the edge number and node number ratio is satisfied 
# and there is a zero further down we should keep going



# may hav efound away
# we count the amount of final zeros and calcualte the difference from the total number of nodes, the result is teh number of edges drawn at the last zero.
# all edges behind those 0s don't exist

#####
#
# bloddy hell nither way works 



# graph plotting

library("qgraph")

Edges <- data.frame(
    from = rep(1:5,each=5),
    to = rep(1:5,times=5),
    thickness = abs(rnorm(25)))

Edges <- subset(Edges,from!=to)
qgraph(Edges,esize=5,gray=TRUE)

qgraph(com2, edge.labels = T)

