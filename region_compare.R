
# read in data and see if it merges well adn rename columns then we can start plotting everyhting accordingly

# we will have 4 circle plots for each catagory and three species to see the overlap of each catagory 

rm(list = ls())
setwd("~/Desktop/Dist_matrix_TE_div/merged.tables/")


Bovine <- read.table("Human_Bovinemerged.CE.class", header = TRUE)
colnames(Bovine)[c(14, 16)] <- c("Bovine_Classification","Bovine_r_class")
Dog <- read.table("Human_Dogmerged.CE.class", header = TRUE)
colnames(Dog)[c(14, 16)] <- c("Dog_Classification","Dog_r_class")
Dog.class <- Dog[,c('bin.ID', 'Dog_Classification', 'Dog_r_class')]
Mouse <- read.table("Human_Mousemerged.CE.class", header = TRUE)
colnames(Mouse)[c(14, 16)] <- c("Mouse_Classification","Mouse_r_class")
Mouse.class <- Mouse[,c('bin.ID', 'Mouse_Classification', 'Mouse_r_class')]

all.spec <- merge(Bovine, Dog.class, by.x=1, by.y=1)
all.spec <- merge(all.spec, Mouse.class, by.x = 1, by.y = 1)

# all the data is merged we can begin to plot things to see what all the species look like in human

require(ggbio)
require(GenomicRanges)
require(topGO)
library(GO.db)
library(org.Hs.eg.db)

#load data into genomic ranges for conversion into plots

UCSCspec = "hg19"
con <- gzcon(url(paste("http://hgdownload.soe.ucsc.edu/goldenPath/",UCSCspec,"/database/chromInfo.txt.gz", sep="")))
txt <- readLines(con)
dat <- read.table(textConnection(txt))
SEQ <- dat[,1:2]
colnames(SEQ) <- c("chr", "size")
sLengths <- SEQ$size
names(sLengths) <- SEQ[,'chr']
sLengths <- sLengths[names(sLengths) %in% as.character(unique(all.spec$Chromosome))] 

hg.bin.gr <- GRanges(seqnames = Rle(all.spec$Chromosome), 
			ranges = IRanges(start = all.spec$Start, end = all.spec$End), 
			seqlengths = sLengths, CE_P = all.spec$P_CE) 
hg.chr <- GRanges(seqnames = Rle(names(sLengths)), 
			ranges = IRanges(start = 1, end = sLengths), 
			seqlengths = sLengths)


#groups

Bov.HH <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "H") & (all.spec$Bovine_r_class == "H"),]))
Mou.HH <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "H") & (all.spec$Mouse_r_class == "H"),]))
Dog.HH <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "H") & (all.spec$Dog_r_class == "H"),]))

Bov.HL <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "L") & (all.spec$Bovine_r_class == "H"),]))
Mou.HL <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "L") & (all.spec$Mouse_r_class == "H"),]))
Dog.HL <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "L") & (all.spec$Dog_r_class == "H"),]))

Bov.LH <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "H") & (all.spec$Bovine_r_class == "L"),]))
Mou.LH <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "H") & (all.spec$Mouse_r_class == "L"),]))
Dog.LH <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "H") & (all.spec$Dog_r_class == "L"),]))

Bov.LL <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "L") & (all.spec$Bovine_r_class == "L"),]))
Mou.LL <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "L") & (all.spec$Mouse_r_class == "L"),]))
Dog.LL <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "L") & (all.spec$Dog_r_class == "L"),]))


Bov.HA <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "H"),]))
Mou.HA <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "H"),]))
Dog.HA <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "H"),]))

Bov.LA <- as.integer(rownames(all.spec[(all.spec$Bovine_Classification == "L"),]))
Mou.LA <- as.integer(rownames(all.spec[(all.spec$Mouse_Classification == "L"),]))
Dog.LA <- as.integer(rownames(all.spec[(all.spec$Dog_Classification == "L"),]))

Bov.HO <- as.integer(rownames(all.spec[(all.spec$Bovine_r_class == "H"),]))
Mou.HO <- as.integer(rownames(all.spec[(all.spec$Mouse_r_class == "H"),]))
Dog.HO <- as.integer(rownames(all.spec[(all.spec$Dog_r_class == "H"),]))

Bov.LO <- as.integer(rownames(all.spec[(all.spec$Bovine_r_class == "L"),]))
Mou.LO <- as.integer(rownames(all.spec[(all.spec$Mouse_r_class == "L"),]))
Dog.LO <- as.integer(rownames(all.spec[(all.spec$Dog_r_class == "L"),]))



#plots


q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white", main = "high PC2 and high organization")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7, main = "thing")
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.HH]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.HH]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.HH]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)


q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.HL]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.HL]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.HL]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)


q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.LH]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.LH]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.LH]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)


q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.LL]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.LL]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.LL]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)



q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.HO]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.HO]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.HO]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)


q <- autoplot((hg.chr),layout = "circle",fill = "white", color = "white")
q <- q + layout_circle(hg.chr, geom = "ideo", fill = "white", color= "black", trackWidth= 6.6, radius = 8.7)
q <- q + layout_circle(hg.chr, geom = "scale", size = 2, radius = 15, trackWidth = 1) 
q <- q + layout_circle(hg.chr, geom = "text", aes(label = seqnames), vjust = 0, size = 3, radius = 16) 
q <- q + layout_circle(reduce(hg.bin.gr[Bov.LO]), geom = "rect", color= "blue", fill = "blue", radius = 13, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Mou.LO]), geom = "rect", color = "red", fill = "red", radius = 11, trackWidth = 2)
q <- q + layout_circle(reduce(hg.bin.gr[Dog.LO]), geom = "rect", color = "green", fill = "green", radius = 9, trackWidth = 2)
print(q)

p <- q + layout_circle(hg.bin.gr, geom = "bar", aes(y = CE_P), grid =FALSE, radius = 15) + scale_size(range = c(1, 2))
p


# the next step is to start thinking about what genes are in these bins with different repeat content
# also how to think about these bins and how much they overlap and wheater or not to use the union or the intersect
# when getting measurments rember to think what aprt of each species we are measureing 

# 
# 
#
# 



# In this script after I do all the circle work I cna start looking for GO terms
# it is probably best I loop this some how to go through each term

# get the genes that overlap a region 
# and find what processes they are enriched for
# make note that IRanges and topGO score functions overwrite each other



# download all refseq data and subset approriatly 
con <- gzcon(url("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"))
txt <- readLines(con)
hg19 <- read.table(textConnection(txt))
hg19eg <- unique(as.character(mget(as.character(hg19[,2]), org.Hs.egREFSEQ2EG,ifnotfound = NA)))
hg19eg <- hg19eg[!(hg19eg == "NA")]
all.genes.eg <- rep(0,length(hg19eg))
names(all.genes.eg) <- hg19eg
hg.gene.gr <- GRanges(seqnames = Rle(hg19[,3]),
				ranges = IRanges(start = hg19[,5], end = hg19[,6]))


# I have a way of selecting terms now I need a way of saving terms

region <- c("Bov.HH","Bov.HL","Bov.LH","Bov.LL",
		    "Mou.HH","Mou.HL","Mou.LH","Mou.LL",
		    "Dog.HH","Dog.HL","Dog.LH","Dog.LL")

# over each loop I change G to be a different region


for(i in seq(region)){
	all.genes.eg[1:length(all.genes.eg)] <- 0
	region.loop <- get(region[i])
	OL <- as.matrix(findOverlaps(hg.gene.gr,hg.bin.gr[region.loop])) 
	G <- as.character(hg19[OL[,1],2])
	m <- unique(as.character(mget(G, org.Hs.egREFSEQ2EG,ifnotfound = NA)))
	m <- m[!(m == "NA")]
	all.genes.eg[names(all.genes.eg) %in% m] <- 1

	GOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP", 
				allGenes = as.factor(all.genes.eg),  
				annot =  annFUN.org, 
				mapping = "org.Hs.eg.db", 
				ID = "entrez",
				nodeSize = 10)
	resultFisher.pc <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher")
	resultFisher.wt <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	resultFisher.wt1 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
	resultFisher.cl <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	
	# tables
	
	allRes.cl <- GenTable(GOdata, classic_Fisher = resultFisher.cl, orderBy = "classic_Fisher", ranksOf = "classic_Fisher", 	topNodes = 20)	
					
	allRes.wt <- GenTable(GOdata, weight_Fisher = resultFisher.wt, orderBy = "weight_Fisher", ranksOf = "weight_Fisher", 	topNodes = 20)				
					
	allRes.wt1 <- GenTable(GOdata, weight_Fisher_01 = resultFisher.wt1, orderBy = "weight_Fisher_01", ranksOf = "weight_Fisher_01", 	topNodes = 20)
	
	allRes.pc <- GenTable(GOdata, parent_child_Fisher = resultFisher.pc, orderBy = "parent_child_Fisher", ranksOf = "parent_child_Fisher", 	topNodes = 20)

	#save tables				
	
	write.table(allRes.cl, file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_enrich/classic/", region[i], "_enrichment_classic_sig.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE )
	
	write.table(allRes.wt, file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_enrich/weight/", region[i], "_enrichment_weight_sig.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE )
	
	write.table(allRes.wt1, file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_enrich/weight01/", region[i], "_enrichment_weight01_sig.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE )

	write.table(allRes.pc, file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_enrich/parent_child/", region[i], "_enrichment_parent_child_sig.txt",sep = ""), quote = FALSE, sep = "\t", row.names = FALSE )
	
	# save DAG plots
	
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_DAG/classic/", region[i], "_DAG_classic_plot.pdf",sep = ""), onefile = TRUE)
	showSigOfNodes(GOdata, score(resultFisher.cl), firstSigNodes = 5, useInfo = 'all')
	dev.off()
	
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_DAG/weight/", region[i], "_DAG_weight_plot.pdf",sep = ""), onefile = TRUE)
	showSigOfNodes(GOdata, score(resultFisher.wt), firstSigNodes = 5, useInfo = 'all')
	dev.off()
	
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_DAG/weight01/", region[i], "_DAG_weight01_plot.pdf",sep = ""), onefile = TRUE)
	showSigOfNodes(GOdata, score(resultFisher.wt1), firstSigNodes = 5, useInfo = 'all')
	dev.off()
	
	pdf(file = paste("~/Desktop/Dist_matrix_TE_div/plot_results/GO_DAG/parent_child/", region[i], "_DAG_parent_cild_plot.pdf",sep = ""), onefile = TRUE)
	showSigOfNodes(GOdata, score(resultFisher.pc), firstSigNodes = 5, useInfo = 'all')
	dev.off()
	
	
	
}
 

inter <- intersect(hg.bin.gr[Dog.LH], intersect(hg.bin.gr[Mou.LH], hg.bin.gr[Bov.LH]))
inter <- intersect(hg.bin.gr[Dog.LH], hg.bin.gr[Bov.LH])

goID <- allRes[1, "GO.ID"]
print(showGroupDensity(GOdata,goID, ranks = TRUE))



showGroupDensity(GOdata,allRes[10, "GO.ID"], ranks = T, rm.one = F)



# it may be worth doing a sanity check with the repeat classes in these bins 
# however i next step will probably come from the line of thinking of what does it all mean. 

