
### ======================================================================== ###
### ------------------------------- LIBRARIES ------------------------------ ###
### ======================================================================== ###

library(cluster)
library(factoextra)
library(ggplot2)
library(hopkins)
library(infotheo)
library(mclust)
library(mixOmics)
library(NbClust)

### ======================================================================== ###
### ------------------------------- DIRECTORY ------------------------------ ###
### ======================================================================== ###

wd = "C:/...P. cinnamomi ML Mapping/Analysis" # SET YOUR WORKING DIRECTORY 
setwd(wd)

### ======================================================================== ###
### -------------------------- READ IN INPUT DATA -------------------------- ###
### ======================================================================== ###

scaled.IV.DV.DF <- IV.DV.DF <- readRDS("Data/Labelled.Predictors.rda")

scaled.IV.DV.DF[, 6:24] <- scale(IV.DV.DF[, 6:24])
scaled.IV.DV.DF$disease[which(scaled.IV.DV.DF$disease==1)] <- "A"
scaled.IV.DV.DF$disease[which(scaled.IV.DV.DF$disease==2)] <- "P"
scaled.IV.DV.DF$disease[which(scaled.IV.DV.DF$disease==3)] <- "N"


scaled.DV.DF <- scaled.IV.DV.DF[, -c(1,2,3,4,5)]

### ----------------- Extract rows of top and bottom polygons -------------- ###

top.polygon.rows <- which(scaled.IV.DV.DF$y>3000)
bottom.polygon.rows <- which(scaled.IV.DV.DF$y<3000)

top.polygon.DF <- scaled.IV.DV.DF[top.polygon.rows,]
bottom.polygon.DF <- scaled.IV.DV.DF[bottom.polygon.rows,]

### ======================================================================== ###
### --------------------------- INTERCORRELATIONS -------------------------- ###
### ======================================================================== ###

# remove the x and y columns for intercorrelations
intercorrelations <- cor(scaled.DV.DF)
# show heatmap of intercors
cor.HM <- pheatmap(intercorrelations, fontsize = 20) ### FIGURE 4 ###
# show dendrogram of features based on intercors
par(cex=1.5, mar=c(5, 8, 4, 1))
plot(cor.HM$tree_col, xlab = "") ### FIGURE 5 ###

### ======================================================================== ###
### ------------------------------- CLUSTERING ----------------------------- ###
### ======================================================================== ###

### ======================================================================== ###
### -------------------------------- K-Means ------------------------------- ###

# yield models using 1 to 15 clusters
km.models <- list()
km.WSS <- rep(0, 15)
for (k in 1:15) {
  
  km <- kmeans(scaled.DV.DF, 
               centers = k, 
               iter.max = 300, 
               algorithm = "MacQueen")
  
  km.models <- append(km.models, list(km))
  km.WSS[k] <- km$tot.withinss
}

### --------------------------- Silhouette Method -------------------------- ###

set.seed(15)
kmean.sil <- fviz_nbclust(scaled.DV.DF, kmeans, method = "silhouette")
kmean.sil # k = 5

### ------------------------------- Gap Method ----------------------------- ###

set.seed(15)
kmean.gap <- fviz_nbclust(scaled.DV.DF, kmeans, method = "gap_stat")
kmean.gap # k = 9, but k = 5 also seems to be a fairly good option

### ---------------------------- Ensemble Method --------------------------- ###

kmeans.ensemble <- NbClust(scaled.DV.DF, distance="euclidean", method="kmeans")
#clusternum # k = 5

### ------------------------------------------------------------------------ ###

# using k = 3 for the three disease levels
base.km.model <- km.models[[3]]
base.km <- list(mi = mutinformation(base.km.model$cluster, IV.DV.DF$disease),
                 ari = adjustedRandIndex(base.km.model$cluster, IV.DV.DF$disease))

# using k = 5
opt.km.model <- km.models[[5]]
opt.km <- list(mi = mutinformation(opt.km.model$cluster, IV.DV.DF$disease),
                  ari = adjustedRandIndex(opt.km.model$cluster, IV.DV.DF$disease))


### ======================================================================== ###
### ---------------------------------- PAM --------------------------------- ###

# yield models using 1 to 15 clusters
pam.models <- list()
for (k in 1:15) {
  
  pm <- pam(scaled.DV.DF,
               k=k)
  
  pam.models <- append(pam.models, list(pm))
}

### --------------------------- Silhouette Method -------------------------- ###

set.seed(15)
pam.sil <- fviz_nbclust(scaled.DV.DF, pam, method = "silhouette")
pam.sil # k = 2

### ------------------------------- Gap Method ----------------------------- ###

set.seed(15)
pam.gap <- fviz_nbclust(scaled.DV.DF, pam, method = "gap_stat")
pam.gap # k = 1

### ------------------------------------------------------------------------ ###

# using k = 3 for the three disease levels
base.pam.model <- pam.models[[3]]
base.pam <- list(mi = mutinformation(base.pam.model$clustering, IV.DV.DF$disease),
                 ari = adjustedRandIndex(base.pam.model$clustering, IV.DV.DF$disease))

# using k = 2
opt.pam.model <- pam.models[[2]]
opt.pam <- list(mi = mutinformation(opt.pam.model$clustering, IV.DV.DF$disease),
                  ari = adjustedRandIndex(opt.pam.model$clustering, IV.DV.DF$disease))

### ======================================================================== ###
### ----------------------------- Hierarchical ----------------------------- ###

dissim.mat <- dist(scaled.DV.DF)
 
hierarchical.clustering <- hclust(dissim.mat)

hier.models <- list()
for (k in 1:15) {
  
  hier.models[[k]] <- cutree(hierarchical.clustering, k = k)
}


### --------------------------- Silhouette Method -------------------------- ###

set.seed(16)
hier.sil <- fviz_nbclust(scaled.DV.DF, hcut, method = "silhouette", k.max = 15)
hier.sil # k = 9

### ------------------------------- Gap Method ----------------------------- ###

set.seed(16)
hier.pam <- fviz_nbclust(scaled.DV.DF, hcut, method = "gap_stat", k.max = 15)
hier.pam # k = 15

### ------------------------------------------------------------------------ ###

# using k = 3 for the three disease levels
base.hier.model <- hier.models[[3]]
base.hier <- list(mi = mutinformation(base.hier.model, IV.DV.DF$disease),
                  ari = adjustedRandIndex(base.hier.model, IV.DV.DF$disease))

# using k = 9
opt.hier.model <- hier.models[[11]]
opt.hier <- list(mi = mutinformation(opt.hier.model, IV.DV.DF$disease),
                  ari = adjustedRandIndex(opt.hier.model, IV.DV.DF$disease))







base.hier
opt.hier

base.km
opt.km

base.pam
opt.pam

models <- c("Hierarchical", "Hierarchical", "Kmeans", "Kmeans", "PAM", "PAM")
opts <- c("baseline", "optimised", "baseline", "optimised", "baseline", "optimised")
mi <- c(base.hier$mi, opt.hier$mi, base.km$mi, opt.km$mi, base.pam$mi, opt.pam$mi)
ari <- c(base.hier$ari, opt.hier$ari, base.km$ari, opt.km$ari, base.pam$ari, opt.pam$ari)

df <- data.frame(model = models,
                 Configuration = opts,
                 mi = mi,
                 ari = ari)

### SUPP. FIGURE 3 ###
p <- ggplot(data=df, aes(x=model, y=mi, fill=Configuration)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  labs(title = "(a) Mutual Information of Clustering algoirthms") +
  ylab("Mutual Information (MI)") + 
  xlab("Model Type")
p

p <- ggplot(data=df, aes(x=model, y=ari, fill=Configuration)) + 
  geom_bar(stat='identity', position=position_dodge()) +
  labs(title = "(b) Adjusted Rand Index of Clustering algoirthms") +
  ylab("Adjusted Rand Index (ARI)") + 
  xlab("Model Type")
p

plot(IV.DV.DF$long, IV.DV.DF$lat, type = "p", col = opt.hier.model)

### ======================================================================== ###
### ------------------------------- DIM RED ----------------------------- ###
### ======================================================================== ###

### SUPP. FIGURE 4 ###
pca.obj <- pca(scaled.DV.DF, ncomp = 5)
plotIndiv(pca.obj, group = scaled.IV.DV.DF$disease,
          title = "(a) Samples on PCA Components",
          ellipse = T, legend = T, legend.title = "Disease")

plsda.obj <- mixOmics::plsda(scaled.DV.DF, scaled.IV.DV.DF$disease)
plotIndiv(plsda.obj, group = scaled.IV.DV.DF$disease,
          title = "(b) Samples on PLSDA Components",
          ellipse = T, legend = T, legend.title = "Disease")



