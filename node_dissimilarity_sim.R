# we have ranked lists 
# we go down until they are duplicated 
# maybe on sample data fist

library(cluster)

# data
x <- 10
S1 <- data.frame(A = rnorm(x), B = rnorm(x))
S2 <- data.frame(S1, C = rnorm(x))

# correlations
cor(S1,S2)

# calculate weights in random data 

d.s1 <- daisy(S1, stand = TRUE)
d.s2 <- daisy(S2, stand = TRUE)
com <- t(combn( dim(S1)[1], m = 2))
colnames(com) <- c("node_1", "node_2")
com1 <- data.frame(com, weights = as.numeric(d.s1)) 
com2 <- data.frame(com, weights = as.numeric(d.s2)) 
# our expectation would be that dissimilrity correlates with C in S2

# workout how to draw the lines that we want


s1 <- com1[order(com1$weights),]

s2 <- com2[order(com1$weights),]
both <- data.frame(s1[,1:2], s1_weights = s1$weights, s2_weights = s2$weights)
rownames(both) <- 1:dim(both)[1]

# find a way of finding guys that are the same

# what if i check the code at its present state and decide what can go together 

# i then redo the check everytime there is a merge

# trying to find edges that can be draw at the same time

# the loop is shit 
# try something else that is more vectorisable 


node.set1 <- both[,c("node_1", "node_2")]
node.set2 <- data.frame(node_1 = c(node.set1[2:dim(node.set1)[1],1],0), node_2 = c(node.set1[2:dim(node.set1)[1],2],0))

node.set3 <- data.frame(node_1 = c(node.set1[3:dim(node.set1)[1],1],rep(0,2)), node_2 = c(node.set1[3:dim(node.set1)[1],2],rep(0,2)))

# I need to mark down when the neodes arent found in the group next to them 


draw1 <- data.frame(node.set1[,1] == node.set2[,1], node.set1[,1] == node.set2[,2], node.set1[,2] == node.set2[,1], node.set1[,2] == node.set2[,2])

# rows that are = to true

Da <- apply(draw1, 1, sum)


draw2 <- data.frame(node.set1[,1] == node.set3[,1], node.set1[,1] == node.set3[,2], node.set1[,2] == node.set3[,1], node.set1[,2] == node.set3[,2])


Db <- apply(draw2, 1, sum)

# where ever we get a true something is not equall to the guy in front of it


data.frame(node.set1, 0,Da,Db)




# easily vectorisable to get all this infomation

# can we from this table go furhter, are any more iterations required 









check <- NULL
draw <- NULL
end <- data.frame(edge = 1:dim(both)[1], previous = 0)
for(i in 1:dim(both)[1]){	
	check <- c(check, as.numeric(both[i + 1,1:2]))
	nod <- as.numeric(both[i,1:2])
	end[i,2] <- length(check)/2
	if(!(any(nod %in% check))){
		draw <- c(draw,i+1)	
	}else{
		check <- NULL
		draw = i
		}
}
end



for(i in 1:dim(both)[1]){
	work.space <- both[i,]
	colnames(work.space) <- colnames(both)
	number.check <- c(work.space$node1, work.space$node_2)
	work.space1 <- NULL
	
	while(!(any(c(work.space1$node_1, work.space1$node_2) %in% number.check))){
		if(dim(work.space)[2] == 1){
			work.space1 <- work.space
		}
		work.space1 <- rbind(work.space1, both[i+1,])
		number.check <- c(work.space$node1, work.space$node_2)
	}
	print(work.space1)
}



# the other hard bit will be keeping track of joins 


# maybe make three seperate data structures that can update and maintain information 


# where do we stop?


