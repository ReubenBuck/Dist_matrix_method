# simulation results analysis


load("~/Desktop/Dist_matrix_TE_div/overlap.sim.list/old_new_rnorm_TEdist")


HH.old <- sapply( olist, function(m) m[1,1] ) 
HL.old <- sapply( olist, function(m) m[2,1] )
LH.old <- sapply( olist, function(m) m[3,1] )
LL.old <- sapply( olist, function(m) m[4,1] )


HH.new <- sapply( olist, function(m) m[1,2] ) 
HL.new <- sapply( olist, function(m) m[2,2] ) 
LH.new <- sapply( olist, function(m) m[3,2] ) 
LL.new <- sapply( olist, function(m) m[4,2] ) 


HO <- sapply( olist, function(m) m[1,3] )
LO <- sapply( olist, function(m) m[3,3] )




boxplot(HH.old, HL.old, LH.old, LL.old, HH.new, HL.new, LH.new, LL.new, HO,LO, notch = TRUE)

boxplot(HH.old - HO, HH.new - HO, HL.old - HO, HL.new - HO)

plot(HH.old - HO, type = "b")
lines(HH.new - HO, col = 2, type = "b")
lines(HL.old - HO, col = 3, type = "b")
lines(HL.new - HO, col = 4, type = "b")
abline(h = 0)


Hstuff <- data.frame(HH.old - HO, HH.new - HO, HL.old - HO, HL.new - HO)
plot(hclust(dist(Hstuff)))
abline(h = .8, col = 2)
plot((HH.old) , (HL.new))


A <- hclust(dist(data.frame(HH.old - HO, HH.new - HO, HL.old - HO, HL.new - HO)))
plot(A)
groups <- cutree(A, h=.8) # cut tree into 5 clusters
rect.hclust(A, h=.8, border="red")




boxplot((HH.old - HO)[groups == 2], (HH.new - HO)[groups == 2], notch = TRUE)

boxplot((HH.old - HO)[A$order[1:81]], (HH.new - HO)[A$order[1:81]], notch = TRUE)








plot(data.frame(HH.old, HL.old, LH.old, LL.old, HH.new, HL.new, LH.new, LL.new, HO,LO))


# what has happened to lead to various patterns. SO we are interested in situations where the old stuff is below the overall in H org regions

# we are randomly purturbing a TE distribution and then looking at bins that have a similar TE content acorss speecies.

# we get a number of different results so what is it about the pertubation that would give similar results as to what we are seeing


# It is safe to say there is no one consistent pattern of conservation of pairwise similarity at different TE values. 

# therefore ther is some variation in the pertubation that has led to the distribtuoin we are seeing. 

old <- rnorm(2000, mean = 100, sd = 30)
s1 <- data.frame(old = old, new = rnorm(2000, mean = 100, sd = 30))
s2 <- data.frame(old = old, new = rnorm(2000, mean = 100, sd = 30))
s3<- data.frame(old = old, new = rnorm(2000, mean = 100, sd = 30))


cor(s1,s2)
cor(s3,s2)
cor(s3,s1)


# so what could we look for that will explain what we are seeing
# some sort of proporty of the new distribtuion in regard to the old distribtuoion

plot(s1$old/s1$new, type = "l")
lines(s2$old/s2$new, col = 2)

plot(s1$old-s1$new,s2$old-s2$new)
plot(s1$old-s1$new,s3$old-s3$new)
plot(s2$old-s2$new,s3$old-s3$new)

plot((s1$old-s1$new) - (s2$old-s2$new), (s1$old-s1$new) - (s3$old/s3$new))

 plot(dist(s1), dist(s2))
# how to look at the differnet ways that the TE contnet is being purturbed. 
# we are looking at places in which the TE content is different between two species

