# This will be where i can get to learn about various GO terms.
#
#using the various R tools used for GO analysis
#
#

library(GO.db)
# a data base that stores GO information in a species independant manner
# in otherwords you access the database with GO terms, form there you can look at
# parents and children of that term.

library(GOsummaries)
# creats a word cloud of GO terms, useful for data presentation

library(topGO)
# a package used to do various GO term enrichment analyses

library(goProfiles)
# get profiles by slicing the GO graph at various levels

library(org.Hs.eg.db)
#this last one seems promising

# maybe it would be better to download a table then access GO terms later using NCBI2R

# the kind of analysis I would want to do is to see if any of the regions I'm interested in 
# are enriched for certain mammalin speciefic or species specific processes. 
# processes or functions that have diverged considarbly since species divergence


# I have little idea of what species specific process could be
# if there was some sort of data base in which i could compare my results against


# download all refseq data and subset approriatly 
con <- gzcon(url("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz"))
txt <- readLines(con)
hg19 <- read.table(textConnection(txt))


# over each loop I change G to be a different region
G = as.character(sample(hg19[,2], 100))
 



#this samll section should find me enriched go terms for each region Im looking at

m <- unique(as.character(mget(G, org.Hs.egREFSEQ2EG,ifnotfound = NA)))
m <- m[!(m == "NA")]
hg19eg <- unique(as.character(mget(as.character(hg19[,2]), org.Hs.egREFSEQ2EG,ifnotfound = NA)))
hg19eg <- hg19eg[!(hg19eg == "NA")]
all.genes.eg <- rep(0,length(hg19eg))
names(all.genes.eg) <- hg19eg
all.genes.eg[names(all.genes.eg) %in% m] <- 1
GOdata <- new("topGOdata", description = "Simple session", ontology = "BP", allGenes = as.factor(all.genes.eg),  annot =  annFUN.org, mapping = "org.Hs.eg.db", ID = "entrez")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)		
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')














plotWordcloud(allRes$Term,1 - as.numeric(allRes$classicFisher), random.order=F, colors = c("blue", "purple","black", "red", "brown"), algorithm = "circle", rot.per = 0)		


# need to work out what tests would be interesting to run and the best way to view results form our comparisons
# at which point do we vieew the tree as well
# a go at using Goprofiles
firstP <- basicProfile(m, onto = "BP", orgPackage = "org.Hs.eg.db")

#plot the results
plotProfiles(firstP)

# also have a word cloud of terms
plotWordcloud(firstP$BP[,1],firstP$BP[,3], random.order=F, colors = c("blue", "purple","black", "red", "brown"), algorithm = "circle", rot.per = 0)



# could also do the same with the reactome

library("clusterProfiler")