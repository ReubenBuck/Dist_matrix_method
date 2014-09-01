rm(list = ls())


 setwd("~/Desktop/Dist_matrix_TE_div/plot_results/edge_species/example")
TE_A <- rnorm(1000, mean = 3)
TE_B <- rnorm(1000, mean = 3)

TE_C <- TE_B + runif(1000,min = -.1, max = .1)
TE_C[750:1000] <- TE_C[750:1000] + 3


S1 <- as.matrix(data.frame(TE_A, TE_B))
S2 <- as.matrix(data.frame(TE_A, TE_C))


plot(S2, xlim = c(0,5), ylim = c(0,2))
plot(S1, xlim = c(0,5), ylim = c(0,2))



d.s1 <- dist(S1)
d.s2 <- dist(S2)


library(ggplot2)

df <- data.frame(distances = c(d.s1,d.s2), species = c(rep("S1", length(d.s1)), rep("S2", length(d.s2))))





s1.mean <- colMeans(as.matrix(d.s1))
s2.mean <- colMeans(as.matrix(d.s2))

edge_length_difference <- s2.mean - s1.mean

edge_means <- data.frame(edges = c(s1.mean, s2.mean), species = c(rep("S1", length(s1.mean)), rep("S2", length(s2.mean))))

ggplot(data = edge_means, aes(edges, fill = species)) + geom_density(alpha = .3) + xlab("mean edge lengths")
ggsave(filename = "distribtuion.pdf")

qplot(data = data.frame(S1), x = TE_A,y = TE_B,  xlim = c(0,6), ylim = c(0,8))
 ggsave(filename = "S1.pdf")
 
qplot(data = data.frame(S2), x = TE_A,y = TE_C,  xlim = c(0,6), ylim = c(0,8))
ggsave(filename = "S2.pdf")





qplot(data = data.frame(S1), x = TE_A, y = TE_B, color = edge_length_difference, xlim = c(0,6), ylim = c(0,8)) 
ggsave(filename = "S1.change.pdf")


qplot(data = data.frame(S2), x = TE_A,y = TE_C, color = edge_length_difference, xlim = c(0,6), ylim = c(0,8))
ggsave(filename = "S2.change.pdf")



