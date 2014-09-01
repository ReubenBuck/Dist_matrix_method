library(cluster)

rm(list = ls())

# data
x <- 1000
# preservation
pr <- 0

S1 <- data.frame(rnorm(x), rnorm(x))
C <- rnorm(x)
S2 <- data.frame(S1 - runif(x, min = -0.5,max = .5) , C = sqrt(C^2))
S2 <- data.frame(rnorm(x), rnorm(x))


# correlations
cor(S1,S2)

# calculate weights in random data 
d.s1 <- daisy(data.frame(S1[S1[,1]>0.5,]),stand=F)
d.s2 <- daisy(S2[S1[,1]>0.5,],stand=F)
com <- t(combn( dim(S1)[1], m = 2))
colnames(com) <- c("node_1", "node_2")
com1 <- data.frame(com, weights = as.numeric(d.s1)) 
com2 <- data.frame(com, weights = as.numeric(d.s2)) 

# workout how to draw the lines that we want
s1 <- com1[order(com1$weights),]
s2 <- com2[order(com1$weights),]
both <- data.frame(s1[,1:2], s1_weights = s1$weights, s2_weights = s2$weights)





plot(S2, pch = as.character(1:x))
points(data.frame(S1), col = 2, pch = as.character(1:x))




h1 <- hist(d.s1)
h2 <- hist(d.s2)
h1.all <- hist(dist(S1))
h2.all <- hist(dist(S2))


plot(h1.all$mids, h1.all$density, type = "l")
lines(h2.all$mids,h2.all$density, col = 2)
lines(h1$mids, h1$density, col = 3)
lines(h2$mids,h2$density, col = 4)

plot(h1$mids, h1$density, type = "l")
lines(h2$mids,h2$density, col = 2)
wilcox.test(d.s1,d.s2, paired = T)


hm1 <- hist(as.matrix(d.s1)[,2])
hm2 <- hist(as.matrix(d.s2)[,2])

plot(hm1$mids, hm1$density, type = "l")
lines(hm2$mids,hm2$density, col = 2)



samp <- sample(1:dim(S1)[1], 50)

hs1 <- hist(dist(S1[samp,]))
hs2 <- hist(dist(S2[samp,]))


plot(hs1$mids, hs1$count, type = "l")
lines(hs2$mids,hs2$count, col = 2)
wilcox.test(dist(S1[samp,]),dist(S2[samp,]), paired = T)

boxplot(dist(S1[samp,]),dist(S2[samp,]) )

# their edge lengths shouldnt differ too much
# so if they do there has been a change in TE content

