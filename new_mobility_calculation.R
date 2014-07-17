



library(matrixStats)

rm(list = ls())

s1 <- data.frame(TE_1 = runif(1000, min=0, max=1000), 
				 TE_2 = rnorm(1000, mean = 500, sd = 100), 
				 TE_3 = runif(1000, min=0, max=1000), 
				 TE_4 = runif(1000, min=0, max=1000), 
				 TE_5 = runif(1000, min=0, max=1000),
				 TE_6 = runif(1000, min=0, max=1000),
				 TE_7 = rnorm(1000, mean = 100, sd = 10),
				 TE_8 = rnorm(1000, mean = 500, sd = 100),
				 TE_9 = runif(1000, min=0, max=1000),
				 TE_10 = runif(1000, min=0, max=1000)
			)


# how to make the s2
# similar to s1 but it needs some speciec changes
move <- 100

s2 <- data.frame(TE_1 = s1[,1] + runif(1000, min=0-move,max=move),
				 TE_2 = s1[,2] + runif(1000, min=0-move,max=move),
				 TE_3 = rnorm(1000, mean = 500, sd = 100),
				 TE_4 = s1[,4] + runif(1000, min=0-move,max=move),
				 TE_5 = s1[,5] + runif(1000, min=0-move,max=move),
				 TE_6 = rnorm(1000, mean = 500, sd = 100),
				 TE_7 = s1[,7] + runif(1000, min=0-move,max=move),
				 TE_8 = s1[,8] + runif(1000, min=0-move,max=move),
				 TE_9 = s1[,9] + runif(1000, min=0-move,max=move),
				 TE_10 = rnorm(1000, mean = 500, sd = 100)
				 )
				 
				 
				 

rm(list = ls())


C.ball = "no"
resid = "no"


s1 <- data.frame(TE_1 = rnorm(1000, mean = 100, sd = 30), 
				 TE_2 = rnorm(1000, mean = 100, sd = 30),
				 TE_3 = rnorm(1000, mean = 100, sd = 30),
				 TE_4 = rnorm(1000, mean = 100, sd = 30),
				 TE_5 = rnorm(1000, mean = 100, sd = 30),
				 TE_6 = rnorm(1000, mean = 500, sd = 100)			
			)

s2 <- data.frame(TE_1 = rnorm(1000, mean = 100, sd = 30), 
				 TE_2 = rnorm(1000, mean = 100, sd = 30),
				 TE_3 = rnorm(1000, mean = 100, sd = 30),
				 TE_4 = rnorm(1000, mean = 100, sd = 30),
				 TE_5 = rnorm(1000, mean = 100, sd = 30), 
				 TE_6 = rnorm(1000, mean = 500, sd = 100)			
			)


s2.s <- as.data.frame(scale(s2))
s1.s <- as.data.frame(scale(s1))

if(C.ball == "yes"){
change <- sample(1:1000,20)
for(i in 1:20){
	s2[change[i],] <- runif(dim(s2)[2], min=0, max=1000)
}
}

D1 <- as.matrix(dist(s1))
D2 <- as.matrix(dist(s2))

D1.s <- as.matrix(dist(s1.s))
D2.s <- as.matrix(dist(s2.s))



# diference method
CD <- D1-D2
m.dist1 <- colSums(sqrt((CD)^2))

CD <- D1.s-D2.s
m.dist1.s <- colMedians(sqrt(CD^2))

# Distance and scaled disctance
actual.D <- NULL
for(i in 1:1000){
	D <- dist(rbind(s1[i,], s2[i,]))
	actual.D<-c(actual.D, as.numeric(D[1]))
}


actual.D.s <- NULL
for(i in 1:1000){
	D <- dist(rbind(s1.s[i,], s2.s[i,]))
	actual.D.s<-c(actual.D.s, as.numeric(D[1]))
}


# residual method
if( resid == "yes"){
CD <- NULL
for(i in 1:1000){
	LM <- lm(D1[i,] ~ D2[i,])
	CD <- rbind(CD,LM$residuals)
}
m.dist2 <- colMedians(sqrt(CD^2))

CD <- NULL
for(i in 1:1000){
	LM <- lm(D1.s[i,] ~ D2.s[i,])
	CD <- rbind(CD,LM$residuals)
}
m.dist2.s <- colMedians(sqrt(CD^2))
}






plot((m.dist1), (actual.D))
cor(actual.D,log(m.dist1))
cor(actual.D,m.dist2)
cor(actual.D.s,m.dist1.s)
cor(actual.D.s,m.dist2.s)


# currently around 60% of the variation in our distrib

# 



# another idea is to take a rank of similarity and show differences in rank as a predictor of actual distance. 

# then again how is that different from the nearest neighbor idea i had before


plot(m.dist2 + m.dist1, actual.D)









# 

diff <- sqrt((D1 - D2)^2)
COL <- colMedians(diff)
pick <- order(COL)[1:1000]

plot(D1[,pick])
plot(D2[,pick])


# so now we go through and get the actuall distances



actual.D <- NULL
for(i in 1:1000){
	D <- dist(rbind(s1[i,], s2[i,]))
	actual.D<-c(actual.D, as.numeric(D[1]))
}


pick.D <- NULL
for(i in 1:1000){
	D <- dist(rbind(D1[i,pick], D2[i,pick]))
	pick.D<-c(pick.D, as.numeric(D[1]))
}

plot(log(pick.D),(actual.D))
cor(log(pick.D), actual.D)

hist(log(pick.D))

# still only about 65 % of varience evplained

# if we actually define things that have moved and look at association statistics

