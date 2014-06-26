# this will be the same thing except using random data
# the idea of this should be to identify regions where distributions differ
# What kind of patterns would we get by random and what would we get if there is two different distributions 
# what would we get if they had the same distribution


# TE 2 follows a different distribtuion in each species
# we will say bins that have the same row number align together 

# however we have TE5 in S2 that is not in S1 and doesnt correlate much with anything

number <- 300


#TE1 and 2 are correlated

# TE 3 and 4 are aswell, but they are not correlated with 1 and 2

TE1_S1 = rnorm(number)
TE1_S1 = TE1_S1[order(TE1_S1)]
TE2_S1 = rnorm(number)
TE2_S1 = TE2_S1[order(TE2_S1)]

TE3_S1 = rnorm(number, mean = 2, sd = 4)
TE4_S1 = TE3_S1 + runif(number)

S1 <- data.frame(TE1_S1, TE2_S1, TE3_S1, TE4_S1)

plot(S1)


TE1_S2 = rnorm(number, mean = 6)
TE1_S2 = TE1_S2[order(TE1_S2)]
TE2_S2 = rnorm(number, mean = 6)
TE2_S2 = TE2_S2[order(TE2_S2)]

TE3_S2 = TE3_S1 + runif(number, min = -1, max = 1)
TE4_S2 = TE3_S2 + runif(number, min = -2, max = 1)

TE5_S2 = TE4_S2 + runif(number, min = -20, max = 20)


S2 <- data.frame(TE1_S2, TE2_S2, TE3_S2, TE4_S2, TE5_S2)

plot(S2)




plot(S1, pch = 16, cex = .5, ylim = c(-5, 10), xlim = c(-5, 10))
points(S2, col =2, pch = 16, cex = .5)
cor(S1,S2)


# find points and there nearest neighbor


# to get something out of what we are doing we need to essentialy build a distance matrix of the first species in another species

# that means we need to get the distance of A to B in the first species for all of A to B, shouldn't be that hard
# using the s_bin table we need to turn that into the S1 dist matrix in S2. 
# I already know how to get a distance from S1 and see what it looks like in S2 now i just need to do that on a bigger scale

# I also need to read up on regression analysis


distS1 <- dist(S1)
distS2 <- dist(S2)

plot(distS1,distS2)
cor(distS1,distS2)


# the sample data has created two distributions of points that are highly correlated, we can see what our test tells us about these in regards to a distance matrix
library("ade4")
mantel.rtest(distS1, distS2, nrepet = 1000)

library("vegan")


m.cor <- mantel.correlog(distS1,distS2)
plot(m.cor)


biplot(prcomp(S1, scale.=T ))
biplot(prcomp(S2, scale.=T ))



# so we can test using partial mantal tests the varaible that is most associated with the varaition between the distance matricies

TE5_distS2 <- dist(TE5_S2)
TE1_distS2 <- dist(TE1_S2)
TE2_distS2 <- dist(TE2_S2)
TE3_distS2 <- dist(TE3_S2)
TE4_distS2 <- dist(TE4_S2)

TE1_distS1 <- dist(TE1_S1)
TE2_distS1 <- dist(TE2_S1)
TE3_distS1 <- dist(TE3_S1)
TE4_distS1 <- dist(TE4_S1)



m.part <- mantel.partial(distS1,distS2,TE5_distS2)

mantel.partial(distS1,distS2,TE5_distS2)


m.cor <- mantel.correlog(distS1,distS2)
plot(m.cor)

m.cor1 <- mantel.correlog(distS2,distS2)
plot(m.cor1)


