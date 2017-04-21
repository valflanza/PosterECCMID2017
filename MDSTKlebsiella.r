library(tidyverse)
library(cluster)
library(dendextend)
library(rgl)
library(ggtree)
library(apcluster)
library(ape)
library(corrplot)
#library(gridExtra)
library(pca3d)
library(randomcoloR)
library(phangorn)

A2A = read.table("./klebsiella/all.dist")
colnames(A2A) = c("Source","Target","Dist","pvalue","sketch")
A2A = A2A %>% separate(sketch,c("sketch","maxSketch"), sep = "\\/")
A2A$sketch = as.numeric(A2A$sketch)
A2A$maxSketch = as.numeric(A2A$maxSketch)
A2A$D2 = round(100*((A2A$maxSketch - A2A$sketch) / A2A$maxSketch))

propcl = 6


#### Making distance matrix and removing outliers  #####
S = A2A %>% group_by(Source) %>% summarise(totalDist = sum(Dist))
Q = quantile(S$totalDist)

A2A.matrix = A2A %>% group_by(Source) %>% mutate(totalDist = sum(Dist)) %>% ungroup() %>% group_by(Target) %>% mutate(totalDist2 = sum(Dist)) %>% filter(totalDist < Q[4], totalDist2 < Q[4]) %>% select(Source,Target,D2) %>% spread(Source,D2)
tmp = A2A.matrix$Target
rownames(A2A.matrix) = tmp

A2A.affinity.matrix = as.matrix(1-A2A.matrix[,-1])

########### Looking for references without clustering approach and just with less correlated strains
A2Acor = cor(A2A.matrix[,-1])
A2Acor = abs(A2Acor)

tmp = which(A2Acor == min(A2Acor), arr.ind = TRUE)
references =rownames(A2Acor)[tmp[1,]]
tmp.cor = A2Acor[!rownames(A2Acor) %in% references, colnames(A2Acor %in%references)]

for(i in 1:(propcl-2))
{
  tmp = which.min(rowSums(tmp.cor))
  references = c(references,names(tmp[1]))
  tmp.cor = A2Acor[!rownames(A2Acor) %in% references, colnames(A2Acor %in%references)]
}

A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
tmp = A2A.coords$Target
rownames(A2A.coords) = tmp
A2A.coords$color1 = rgb(A2A.coords[,2],A2A.coords[,3],A2A.coords[,4],maxColorValue = 100)
A2A.coords$color2 = rgb(A2A.coords[,5],A2A.coords[,6],A2A.coords[,7],maxColorValue = 100)

pca.correlation = prcomp(A2A.coords[,2:(propcl+1)])
D.correlation = dist(A2A.coords[,2:(propcl+1)])
tree.correlation = agnes(D.correlation) %>% as.dendrogram()

#tree = as.data.frame(as.phylo(H))
#tree.correlation = full_join(tree,select(A2A.coords, label = Target,color1,color2))



########### Looking for references with clustering approach and correlation
A2A.clara = clara(A2A.matrix[,-1],k = propcl *3)
A2A.cluster = A2A.matrix[,-1]
A2A.cluster = A2A.cluster[,A2A.clara$i.med]
A2A.cluster.cor = abs(cor(A2A.cluster))

tmp = which(A2A.cluster.cor == min(A2A.cluster.cor), arr.ind = TRUE)
references =rownames(A2A.cluster.cor)[tmp[1,]]
tmp.cor = A2A.cluster.cor[!rownames(A2A.cluster.cor) %in% references, colnames(A2A.cluster.cor %in%references)]

for(i in 1:(propcl-2))
{
  tmp = which.min(rowSums(tmp.cor))
  references = c(references,names(tmp[1]))
  tmp.cor = A2A.cluster.cor[!rownames(A2A.cluster.cor) %in% references, colnames(A2A.cluster.cor %in%references)]
}

A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
A2A.coords$color1 = rgb(A2A.coords[,2],A2A.coords[,3],A2A.coords[,4],maxColorValue = 100)
A2A.coords$color2 = rgb(A2A.coords[,5],A2A.coords[,6],A2A.coords[,7],maxColorValue = 100)
pca.clara.correlation = prcomp(A2A.coords[,2:(propcl+1)])
tmp = A2A.coords$Target
rownames(A2A.coords) = tmp
D.correlation.clara = dist(A2A.coords[,2:(propcl+1)])
tree.correlation.clara = agnes(D.correlation.clara) %>% as.dendrogram()
#tree = as.data.frame(as.phylo(H))
#tree.correlation.clara = full_join(tree,select(A2A.coords, label = Target,color1,color2))



########### Looking for references with just clustering approach ("clara" clustering)
A2A.clara = clara(A2A.matrix[,-1], k = propcl)

references = rownames(A2A.matrix)[A2A.clara$i.med]
A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
A2A.coords$color1 = rgb(A2A.coords[,2],A2A.coords[,3],A2A.coords[,4],maxColorValue = 100)
A2A.coords$color2 = rgb(A2A.coords[,4],A2A.coords[,5],A2A.coords[,6],maxColorValue = 100)
#A2A.coords$color3 = rgb(A2A.coords[,8],A2A.coords[,9],A2A.coords[,10],maxColorValue = 100)
pca.clara = prcomp(A2A.coords[,2:(propcl+1)])
tmp = A2A.coords$Target
rownames(A2A.coords) = tmp
D.clara = dist(A2A.coords[,2:(propcl+1)])
tree.clara = agnes(D.clara)%>% as.dendrogram()
#tree.clara = agnes(D.clara, )


########### Looking for references with just clustering approach ("Affinity Propagation" clustering)

A2A.affinityK = apclusterK(A2A.affinity.matrix, K= propcl, bimaxit = 30)
references = names(A2A.affinityK@exemplars)
A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
A2A.coords$color1 = rgb(A2A.coords[,2],A2A.coords[,3],A2A.coords[,4],maxColorValue = 100)
A2A.coords$color2 = rgb(A2A.coords[,5],A2A.coords[,6],A2A.coords[,7],maxColorValue = 100)

pca.affinity = prcomp(A2A.coords[,2:(propcl+1)])
tmp = A2A.coords$Target
rownames(A2A.coords) = tmp
D.affinity = dist(A2A.coords[,2:(propcl+1)])
tree.affinity = agnes(D.affinity)%>% as.dendrogram()
#tree = as.data.frame(as.phylo(H))
#tree.affinity = full_join(tree,select(A2A.coords, label = Target,color1,color2))



########### Looking for references with clustering approach ("Affinity Propagation" clustering) + correlation 
A2A.affinity = apcluster(A2A.affinity.matrix)
references = names(A2A.affinity@exemplars)
A2A.cluster = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)

A2A.cluster.cor = abs(cor(A2A.cluster[,-1]))

tmp = which(A2A.cluster.cor == min(A2A.cluster.cor), arr.ind = TRUE)
references =rownames(A2A.cluster.cor)[tmp[1,]]
tmp.cor = A2A.cluster.cor[!rownames(A2A.cluster.cor) %in% references, colnames(A2A.cluster.cor %in%references)]

for(i in 1:(propcl-2))
{
  tmp = which.min(rowSums(tmp.cor))
  references = c(references,names(tmp[1]))
  tmp.cor = A2A.cluster.cor[!rownames(A2A.cluster.cor) %in% references, colnames(A2A.cluster.cor %in%references)]
}

A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
A2A.coords$color1 = rgb(A2A.coords[,2],A2A.coords[,3],A2A.coords[,4],maxColorValue = 100)
A2A.coords$color2 = rgb(A2A.coords[,5],A2A.coords[,6],A2A.coords[,7],maxColorValue = 100)


tmp = A2A.coords$Target
rownames(A2A.coords) = tmp

pca.affinity.correlation = prcomp(A2A.coords[,2:(propcl+1)])

D.affinity.cor = dist(A2A.coords[,2:(propcl+1)])
tree.affinity.cor = agnes(D.affinity.cor) %>% as.dendrogram()
#tree = as.data.frame(as.phylo(H))
#tree.affinity.cor = full_join(tree,select(A2A.coords, label = Target,color1,color2))
#ggtree(tree.affinity.cor, layout = "circular") + geom_tippoint(aes(color = color1))

tmp = A2A %>% select(Source,Target,D2) %>% spread(Source,D2)
tmp2 = tmp$Target
rownames(tmp) = tmp2
D.all = as.dist(tmp[,-1])
tree.all = agnes(as.dist(D.all)) %>% as.dendrogram()
pca.all = prcomp(tmp[,-1])

###### Correlation matrix ######




listOfDist = list(Reference = D.all, Affinity = D.affinity, AffinityCor = D.affinity.cor, Clara = D.clara, Correlation = D.correlation, ClaraCorrelation = D.correlation.clara)
Corrs = c()
for (i in 1:length(listOfDist))
{
  
  Corrs = c(Corrs,cor(D.all,listOfDist[[i]]))
}
Coors = data.frame(Type = names(listOfDist),Correlations = Corrs)





###### Correlation vs number of cluster (usin clara and distance matrix) ######


corV = c()
for (i in 2:20)
{
  A2A.clara = clara(A2A.matrix[,-1], k = i)
  
  references = rownames(A2A.matrix)[A2A.clara$i.med]
  A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
  tmp = A2A.coords$Target
  rownames(A2A.coords) = tmp
  D = dist(A2A.coords[,2:i+1])
  corV = c(corV,cor(D.all,D))
}
corV = as.data.frame(corV)
colnames(corV) = c("Y")
corV$X = 1:19

pV = ggplot(corV, aes(x = X+1, y = Y)) + geom_point() + geom_smooth() + theme_bw() + scale_x_continuous(breaks = c(2:20)) + scale_y_continuous(breaks = seq(0,1,0.05)) + geom_hline(aes(yintercept  =0.8), color = "darkred", linetype = "dotdash") + xlab("Number of Clusters") + ylab("Correlation")


######Correlation vs number of cluster (usin clara and tree distance) ######
corT = c()
for (i in 2:20)
{
  A2A.clara = clara(A2A.matrix[,-1], k = i)
  
  references = rownames(A2A.matrix)[A2A.clara$i.med]
  A2A.coords = A2A %>% filter(Source %in% references) %>% select(Source,Target,D2) %>% spread(Source,D2)
  tmp = A2A.coords$Target
  rownames(A2A.coords) = tmp
  D = dist(A2A.coords[,2:i+1])
  tree.tmp = agnes(D)
  corT = c(corT,cor_cophenetic(tree.all,tree.tmp))
}
corT = as.data.frame(corT)
colnames(corT) = c("Y")
corT$X = 1:19

pT = ggplot(corT, aes(x = X+1, y = Y)) + geom_point() + geom_line() + theme_bw() + scale_x_continuous(breaks = c(2:20)) + scale_y_continuous(breaks = seq(0,1,0.05)) + geom_hline(aes(yintercept  =0.8), color = "darkred", linetype = "dotdash") + labs(x ="Number of Clusters", y ="Correlation", title = "Cophenetic Correlation")

###### Ploting trees with color ######
tree.all.nj = midpoint(nj(D.all))
tree.clara.nj = midpoint(nj(D.clara))

tree.all.dt = as.data.frame(tree.all.nj)
tmp = tree.all.dt %>% filter(isTip) %>% arrange(y) %>% select(label)
tmp$color = rainbow(length(tmp$label))

tree.all.ggtree = ggtree(tree.all.nj) %<+% tmp
tree.clara.ggtree = ggtree(tree.clara.nj) %<+% tmp

gridExtra::grid.arrange(tree.all.ggtree + geom_tippoint(aes(color = color))+ labs(title = "Tree All"),tree.clara.ggtree + geom_tippoint(aes(color = color))+ labs(title = "Tree Clara"), ncol = 2)

###### Ploting PCAs with color ######

tmp = clara(as.matrix(D.all), 20)
colores = distinctColorPalette(max(tmp$clustering))[tmp$clustering]

pairs(pca.all$x[,1:6],col = colores, lower.panel = NULL, pch = 16)
pairs(pca.clara$x[,1:6], col = colores, upper.panel = NULL, pch = 16)


#### Correlation plot ####

colnames(corT) = c("TreeCorrelation","X")
colnames(corV) = c("DistCorrelation","X")

tmp = full_join(corT,corV)
tmp %>% gather("Correlation","Value", -X) %>% ggplot(aes(x = X+1, y = Value, color = Correlation)) + geom_point() + geom_line() + theme_bw() + scale_x_continuous(breaks = c(2:20)) + scale_y_continuous(breaks = seq(0,1,0.05), limits = c(0.7,1)) + geom_hline(aes(yintercept  =0.8), color = "darkred", linetype = "dotdash") + labs(x ="Number of Clusters", y ="Correlation", title = "Correlation")
