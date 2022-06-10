
### ======================================================================== ###
### ------------------------------- DIRECTORY ------------------------------ ###
### ======================================================================== ###

wd = "C:/...P. cinnamomi ML Mapping/Analysis" # SET YOUR WORKING DIRECTORY 
setwd(wd)
source("Internals/Internals.Supervised.R")

### ======================================================================== ###
### -------------------------- READ IN INPUT DATA -------------------------- ###
### ======================================================================== ###

DF <- readRDS("Data/Labelled.Predictors.rda")

DF$disease[which(DF$disease==1)] <- "A" # convert class back to character
DF$disease[which(DF$disease==2)] <- "P"
DF$disease[which(DF$disease==3)] <- "N"

scaled.DF <- DF

scaled.DF[, 6:24] <- scale(scaled.DF[, 6:24]) # scale all predictors



### ======================================================================== ###
### -------------------------- NORAMLITY ANALYSIS -------------------------- ###
### ======================================================================== ###

# for each feature, run S-W test and extract p.value
shapiro.wilk.p.values <- list()
for (c in 6:24) {
  IV <- names(DF)[c]
  shapiro.wilk.p.values[IV] <- shapiro.test(DF[, IV])$p.value
}
shapiro.wilk.p.values <- unlist(shapiro.wilk.p.values)

which(shapiro.wilk.p.values < 0.05) # all are less then 0.05, none have normal distribution

### ======================================================================== ###
### ------------------------- ASSOCIATION ANALYSIS ------------------------- ###
### ======================================================================== ###

### -------------------------- Kruskal-Wallis Test ------------------------- ###

kruskal.wallis.p.values <- list()
DV <- DF$disease

# for each feature, run K-W test and extract p.value
for (c in 6:24) {
  IV.name <- names(DF)[c]
  IV <- as.numeric(DF[,c])
  tmp.df <- data.frame(matrix(NA, nrow=446, ncol=0))
  tmp.df$IV <- IV
  tmp.df$DV <- DV
  kruskal.wallis.p.values[IV.name] <- kruskal.test(IV ~ DV, data=tmp.df)$p.value
}
kruskal.wallis.p.values <- unlist(kruskal.wallis.p.values)

p <- ggplot(data=data.frame(Feature = names(kruskal.wallis.p.values),
                            KW.Pvalues = log(unname(kruskal.wallis.p.values))), 
            aes(x=Feature, y=KW.Pvalues)) +
  geom_bar(stat='identity') +
  geom_hline(yintercept = log(0.05), linetype="dashed", 
                color = "red", size=2) + 
  theme(axis.text.x = element_text(size=22, angle=45),
        title = element_text(size=15)) +
  labs(title = "Kruskal-Wallis P-value for each feature") +
  ylab("log(Kruskal-Wallis P-value)") + 
  xlab("Feature")
p ### FIGURE 6 ###

which(kruskal.wallis.p.values < 0.05)
# some are less than 0.05, 14 have significant difference between groups


### ------------------------- Pairwise-Wilcoxon Test ------------------------- ###

bin.mat <- matrix(0, nrow = 3, ncol = 14)
colnames(bin.mat) <- names(which(kruskal.wallis.p.values < 0.05)) # names(which(kruskal.wallis.p.values < 0.05))
rownames(bin.mat) <- c("P vs A", "N vs A", "P vs N")

for (feat in names(which(kruskal.wallis.p.values < 0.05))) { # names(which(kruskal.wallis.p.values < 0.05)) 
  print(feat)

  # pairwise.wilcox.pvalues
  pw.w.pv <- pairwise.wilcox.test(DF[,feat], DF$disease,
                                  p.adjust.method = "BH")$p.value
  sig.vals <- which(pw.w.pv < 0.05, arr.ind = T)
  print(pw.w.pv)
  if (nrow(sig.vals) > 0) {
    for (r in 1:nrow(sig.vals)) {
      col <- sig.vals[r, "col"]
      row <- sig.vals[r, "row"]
      
      
      print(row)
      print(col)
      
      pair <- paste0(rownames(pw.w.pv)[row], " vs ", colnames(pw.w.pv)[col])
      print(pair)
      bin.mat[pair, feat] <- 1
    }
  }
}

df <- expand.grid(X=names(which(kruskal.wallis.p.values < 0.05)), Y = c("P vs A", "N vs A", "P vs N")) # names(which(kruskal.wallis.p.values < 0.05))
df$val <- expand.grid(bin.mat)

p <- ggplot(df, aes(unlist(X), unlist(Y), fill=unlist(val))) +
      geom_tile(color = "black") + 
      theme(axis.text = element_text(size=22),
        axis.text.x = element_text(angle=45),
        title = element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18)) +
      labs(title = "Heatmap of significant Pairwise-Wilcoxon Tests") +
      ylab("Class Pairs") + 
      xlab("Feature") + 
  guides(fill=guide_legend(title="Significant?"))
p

### ----------------------------- Visualiation ----------------------------- ###

# use these to visualise for ONE GIVEN FEATURE
ggboxplot(DF, x = "disease", y = "DEM",
          color = "disease", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          add = c("mean_se", "jitter"),
          order = c("N", "P", "A"),
          ylab = "DEM", xlab = "Disease State")

ggline(DF, x = "disease", y = "PVR",
          add = c("mean_se", "jitter"),
          order = c("N", "P", "A"),
          ylab = "Blue Reflectance", xlab = "Disease State")


### ======================================================================== ###
### ------------------- MULTINOMIAL LOGISTIC REGRESSION -------------------- ###
### ======================================================================== ###

mnlr.DF <- scaled.DF[, -c(1,2,3,4)]

# ENSURE THIS LINE RUNS BEFORE ANY USE OF PlotMultiplemodelMetrics()
mnlr.metrics <- readRDS("Metrics/MNLR/mnlr.metrics.rda") 

PlotOneModelMetrics(mnlr.metrics$Standard, "MNLR", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "(a) MNLR, Standard Model Metrics")

PlotOneModelMetrics(mnlr.metrics$Stable.Features, "MNLR", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "MNLR, Stable Features Model Metrics")

PlotOneModelMetrics(mnlr.metrics$Regularised, "MNLR", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "MNLR, Regularised Model Metrics")

PlotMultipleModelMetrics(mnlr.metrics,
                         title = "MNLR, All Model Metrics", 
                         plot.type="dot")



mnlr.metrics.df <- readRDS("Metrics/MNLR/mnlr.metrics.df.rda")

mnlr.stnd.imps <- mnlr.metrics$Standard$stability[, -1]
PlotModelFeatImp(mnlr.stnd.imps, mnlr.DF$disease, "MNLR", title = "Stabilities of features in MNLR, Standard model")

# ### ------------------------------- Standard ------------------------------- ###
# 
# basic.mnlr <- RepCV.MNLR(mnlr.DF, nrepeat = 100, nfolds=5)
# 
# ### ----------------------- MNLR on Stable Features ------------------------ ###
# 
# AV.feats <- names(basic.mnlr$stability)[which(basic.mnlr$stability["P", -1] > 0.8)+1]
# AN.feats <- names(basic.mnlr$stability)[which(basic.mnlr$stability["N", -1] > 0.8)+1]
# PN.feats <- names(basic.mnlr$stability)[which(basic.mnlr$stability["A", -1] > 0.8)+1]
# 
# stab.mnlr.DF <- mnlr.DF[, c("disease", unique(c(AV.feats, AN.feats, PN.feats)))]
# stbl.mnlr <- RepCV.MNLR(stab.mnlr.DF, nrepeat = 100, nfolds=5)
# 
# ### ----------------------------- Regularised ------------------------------ ###
# 
# reg.mnlr <- RepCV.MNLR(mnlr.DF, nrepeat = 100, nfolds=5, reg.alpha=1)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# mnlr.metrics <- list(basic.mnlr,
#                      stbl.mnlr,
#                      reg.mnlr)
# names(mnlr.metrics) <- c("Standard", "Stable.Features", "Regularised")
# 
# saveRDS(mnlr.metrics, file = "Metrics/mnlr.metrics.rda")
# 
# mnlr.metrics.df <- MakeMetricsDFfromList(mnlr.metrics)
# 
# saveRDS(mnlr.metrics.df, file = "Metrics/mnlr.metrics.df.rda")

### ======================================================================== ###
### ----------------------------- NAIVE BAYES ------------------------------ ###
### ======================================================================== ###

nb.DF <- scaled.DF[, -c(1,2,3,4)]

nb.metrics <- readRDS("Metrics/NB/nb.metrics.rda")

PlotOneModelMetrics(nb.metrics$Disc.Freq, "NB", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "NB, Disc.Freq Model Metrics")

PlotOneModelMetrics(nb.metrics$Disc.Interval, "NB", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "NB, Disc.Interval Model Metrics")

PlotOneModelMetrics(nb.metrics$Density.Gaussian, "NB", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "NB, Density.Gaussian Model Metrics")

PlotOneModelMetrics(nb.metrics$Density.Kernel, "NB", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "NB, Density.Kernel Model Metrics")

nb.metrics.2 <- nb.metrics
nb.metrics.2$MNLR <- mnlr.metrics$Standard
nb.metrics.2 <- nb.metrics.2[c(5,1,2,3,4)]
PlotMultipleModelMetrics(nb.metrics.2,
                         title = "NB, All Model Metrics", 
                         plot.type="dot")

# ### ------------------------- Discretised by Frequency --------------------- ###
# 
# nb.freq.basic <- RepCV.NB(nb.DF, "frequency", nrepeat = 100, nfolds = 5)
# 
# ### ------------------------- Discretised by Interval ---------------------- ###
# 
# nb.intv.basic <- RepCV.NB(nb.DF, "interval", nrepeat = 100, nfolds = 5)
# 
# ### ----------------------- Density with Gaussian Dist. -------------------- ###
# 
# nb.norm.basic <- RepCV.NB(nb.DF, kernel = F, nrepeat = 100, nfolds = 5)
# 
# ### ------------------------ Density with Kernel Dist. --------------------- ###
# 
# nb.kern.basic <- RepCV.NB(nb.DF, kernel = T, nrepeat = 100, nfolds = 5)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# nb.metrics <- list(nb.freq.basic, 
#                    nb.intv.basic, 
#                    nb.norm.basic, 
#                    nb.kern.basic)
# names(nb.metrics) <- c("Disc.Freq", "Disc.Interval", "Density.Gaussian", "Density.Kernel")
# 
# saveRDS(nb.metrics, file = "Metrics/NB/nb.metrics.rda")
# 
# nb.metrics.df <- MakeMetricsDFfromList(nb.metrics.2)
# saveRDS(nb.metrics.df, file = "Metrics/NB/nb.metrics.df.rda")

### ======================================================================== ###
### ------------------------ K-NEAREST NEIGHBOURS -------------------------- ###
### ======================================================================== ###

knn.DF <- DF[, -c(1,2,3,4)]

knn.metrics <- readRDS("Metrics/kNN/knn.metrics.rda")

PlotOneModelMetrics(knn.metrics$Basic, "KNN", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "KNN, Basic Model Metrics")

PlotOneModelMetrics(knn.metrics$Euclidean, "KNN", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "KNN, Euclidean Model Metrics")

PlotOneModelMetrics(knn.metrics$Manhattan, "KNN", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "KNN, Manhattan Model Metrics")

PlotOneModelMetrics(knn.metrics$Minkowski, "KNN", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "KNN, Minkowski Model Metrics")

knn.metrics.2 <- knn.metrics
knn.metrics.2$MNLR <- mnlr.metrics$Standard
knn.metrics.2 <- knn.metrics.2[c(5,1,2,3,4)]
PlotMultipleModelMetrics(knn.metrics.2,
                         title = "KNN, All Model Metrics", 
                         plot.type="dot")

# knn.basic <- RepCV.kNN(knn.DF, k = 3, nrepeat = 100, nfolds = 5)
# 
# ### --------------------- Tuning k for each distance ----------------------- ###
# 
# knn.tune.eucl <- PlotTune.k.kNN(knn.DF, kmax=10, dist = "euclidean", nrepeat = 100)
# knn.tune.mhtn <- PlotTune.k.kNN(knn.DF, kmax=10, dist = "manhattan", nrepeat = 100)
# knn.tune.mnks <- PlotTune.k.kNN(knn.DF, kmax=10, dist = "minkowski", nrepeat = 100)
# 
# p <- ggplot(knn.tune.eucl, aes(x = k, y = mean)) + 
#     geom_line(aes(color = variable), size=1) + 
#     geom_point(aes(color = variable), size=2) +
#     scale_color_manual(values = c("darkred", "steelblue", "magenta1", "seagreen", "tomato")) +  
#     scale_x_continuous(breaks=1:10) + 
#     ylab("Metrics") + 
#   labs(title = "kNN Tuning: Balanced Accruacy, MNKS distance")
# p
# 
# # extract optimal k values
# eucl.opt.k <- 4
# mhtn.opt.k <- 4
# mnks.opt.k <- 4
# 
# ### ----------------------- Optimised k, Euclidean ------------------------- ###
# 
# knn.eucl <- RepCV.kNN(knn.DF, k = eucl.opt.k, dist = "euclidean", 
#                       nrepeat = 100, nfolds = 5)
# 
# ### ----------------------- Optimised k, Manhattan ------------------------- ###
# 
# knn.mhtn <- RepCV.kNN(knn.DF, k = mhtn.opt.k, dist = "manhattan", 
#                       nrepeat = 100, nfolds = 5)
# 
# ### ----------------------- Optimised k, Minkowski ------------------------- ###
# 
# knn.mnks <- RepCV.kNN(knn.DF, k = mnks.opt.k, dist = "minkowski", 
#                       nrepeat = 100, nfolds = 5)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# knn.metrics <- list(knn.basic,
#                     knn.eucl, 
#                     knn.mhtn, 
#                     knn.mnks)
# names(knn.metrics) <- c("Basic", "Euclidean", "Manhattan", "Minkowski")
# 
# saveRDS(knn.metrics, file = "Metrics/kNN/knn.metrics.rda")
# 
# knn.metrics.df <- MakeMetricsDFfromList(knn.metrics.2)
# knn.metrics.df
# saveRDS(knn.metrics.df, file = "Metrics/kNN/knn.metrics.df.rda")

### ======================================================================== ###
### --------------------------- DECISION TREES ----------------------------- ###
### ======================================================================== ###

dt.DF <- scaled.DF[, -c(1,2,3,4)]

dt.metrics <- readRDS("Metrics/DT/dt.metrics.rda")

PlotOneModelMetrics(dt.metrics$Gini.Basic, "DT", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "DT, Gini.Basic Model Metrics")

PlotOneModelMetrics(dt.metrics$IG.Basic, "DT", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "DT, IG.Basic Model Metrics")

PlotOneModelMetrics(dt.metrics$Gini.Optimised, "DT", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "DT, Gini.Optimised Model Metrics")

PlotOneModelMetrics(dt.metrics$IG.Optimised, "DT", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "DT, IG.Optimised Model Metrics")

dt.metrics.2 <- dt.metrics
dt.metrics.2$MNLR <- mnlr.metrics$Standard
#dt.metrics.2$KNN <- knn.metrics$Minkowski
dt.metrics.2 <- dt.metrics.2[c(5,1,2,3,4)] #  c(5,6,1,2,3,4)
PlotMultipleModelMetrics(dt.metrics.2,
                         title = "DT, All Model Metrics", 
                         plot.type="dot")

dt.gini.opt.imps <- dt.metrics$Gini.Optimised$importance + ((100-sum(dt.metrics$Gini.Optimised$importance))/19)
dt.gini.opt.imps <- dt.gini.opt.imps/max(dt.gini.opt.imps) * 100
PlotModelFeatImp(dt.gini.opt.imps, type="DT", title="Importance of features in DT, Optimised Gini model")

# DT overfitting
dt.over.fitting <- data.frame(Overfitting = unlist(lapply(dt.metrics, function(x){mean(x$over.fitting, na.rm=T)})),
                              Accuracy = unlist(lapply(dt.metrics, function(x){x$metrics$weighted.avg.total$`Balanced Accuracy`})))

saveRDS(dt.over.fitting, file = "Metrics/DT/dt.over.fitting.rda")

# ### ------------------------- Basic Gini Index ----------------------------- ###
# 
# dt.gini.basic <- RepCV.DT(dt.DF, split.crit = "gini", nrepeat = 100)
# 
# ### ---------------------- Basic Information Gain -------------------------- ###
# 
# dt.ig.basic <- RepCV.DT(dt.DF, split.crit = "information", nrepeat = 100)
# 
# ### ----------------------- Optimised Gini Index --------------------------- ###
# 
# gini.tuning <- PlotTune.cp.minsplit.DT(dt.DF, nrepeat = 100, metric = "Balanced Accuracy",
#                                                  split.crit = "gini")
# 
# gini.tuning[which(gini.tuning[,3] == max(gini.tuning[,3])),]
# gini.opt.params <- list(cp = 0.04, minsplit = 30)
# 
# dt.gini.opt <- RepCV.DT(dt.DF, split.crit = "gini", nrepeat = 100,
#                                   cp = gini.opt.params$cp,
#                                   minsplit = gini.opt.params$minsplit)
# 
# ### -------------------- Optimised Information Gain ------------------------ ###
# 
# ig.tuning <- PlotTune.cp.minsplit.DT(dt.DF, nrepeat = 100, metric = "Balanced Accuracy",
#                                                  split.crit = "information")
# 
# ig.tuning[which(ig.tuning[,3] == max(ig.tuning[,3])),]
# ig.opt.params <- list(cp = 0.02, minsplit = 10)
# 
# dt.ig.opt <- RepCV.DT(dt.DF, split.crit = "information", nrepeat = 100,
#                                 cp = ig.opt.params$cp,
#                                 minsplit = ig.opt.params$minsplit)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# dt.metrics <- list(dt.gini.basic, 
#                    dt.ig.basic, 
#                    dt.gini.opt,
#                    dt.ig.opt)
# names(dt.metrics) <- c("Gini.Basic", "IG.Basic", "Gini.Optimised", "IG.Optimised")
# 
# saveRDS(dt.metrics, file = "Metrics/DT/dt.metrics.rda")
# 
# dt.metrics.df <- MakeMetricsDFfromList(dt.metrics.2)
# dt.metrics.df
# saveRDS(dt.metrics.df, file = "Metrics/DT/dt.metrics.df.rda")

### ======================================================================== ###
### --------------------------- RANDOM FORESTS ----------------------------- ###
### ======================================================================== ###

rf.DF <- scaled.DF[, -c(1,2,3,4)]

rf.metrics <- readRDS("Metrics/RF/rf.metrics.rda")

PlotOneModelMetrics(rf.metrics$Optimised, "RF", 
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "RF, Optimised Model Metrics")

rf.metrics.2 <- rf.metrics
rf.metrics.2$MNLR <- mnlr.metrics$Standard
rf.metrics.2 <- rf.metrics.2[c(3,1,2)]
PlotMultipleModelMetrics(rf.metrics.2,
                         title = "RF, All Model Metrics", 
                         plot.type="dot")

rf.gini.opt.imps <- WeightedAverage(rf.DF$disease, rf.metrics$Optimised$importance)
rf.gini.opt.imps <- rf.gini.opt.imps + ((100-sum(rf.gini.opt.imps))/19)
rf.gini.opt.imps <- rf.gini.opt.imps/max(rf.gini.opt.imps) * 100
PlotModelFeatImp(rf.gini.opt.imps, type="RF", title="Importance of features in RF, Optimised model")

# ### -------------------------------- Basic --------------------------------- ###
# 
# rf.basic <- RepCV.RF(rf.DF, nrepeat = 100)
# 
# ### ------------------------------ Optimised ------------------------------- ###
# 
# rf.tuning <- PlotTune.mtry.ntree.RF(rf.DF, nrepeat = 50)
# 
# rf.tuning[which(rf.tuning[,3] == max(rf.tuning[,3])),]
# rf.opt.params <- list(ntree = 800, mtry = 8)
# 
# rf.opt <- RepCV.RF(rf.DF, nrepeat = 100, 
#                              ntree = rf.opt.params$ntree,
#                              mtry = rf.opt.params$mtry)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# rf.metrics <- list(rf.basic, 
#                    rf.opt)
# names(rf.metrics) <- c("Basic", "Optimised")
# 
# saveRDS(rf.metrics, file = "Metrics/RF/rf.metrics.rda")
# 
# rf.metrics.df <- MakeMetricsDFfromList(rf.metrics.2)
# rf.metrics.df
# saveRDS(rf.metrics.df, file = "Metrics/RF/rf.metrics.df.rda")


### ======================================================================== ###
### --------------------------------- SVM ---------------------------------- ###
### ======================================================================== ###

svm.DF <- scaled.DF[, -c(1,2,3,4)]

svm.metrics <- readRDS("Metrics/SVM/svm.metrics.rda")

PlotOneModelMetrics(svm.metrics$opt.sigmoid, "SVM",
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "SVM, Optimised Sigmoid Model Metrics")

svm.metrics.2 <- svm.metrics
svm.metrics.2[["MNLR"]] <- mnlr.metrics$Standard
svm.metrics.2 <- svm.metrics.2[c(9,1,2,3,4,5,6,7,8)]
PlotMultipleModelMetrics(svm.metrics.2,
                         title = "SVM, All Model Metrics", 
                         plot.type="dot")

# ### ----------------------------- Basic Linear ----------------------------- ###
# 
# svm.linear.basic <- RepCV.SVM(svm.DF, kernel = "linear", nrepeat = 100)
# 
# ### ----------------------------- Basic Radial ----------------------------- ###
# 
# svm.radial.basic <- RepCV.SVM(svm.DF, kernel = "radial", nrepeat = 100)
# 
# ### --------------------------- Basic Polynomial --------------------------- ###
# 
# svm.polyn.basic <- RepCV.SVM(svm.DF, kernel = "polynomial", nrepeat = 100)
# 
# ### ---------------------------- Basic Sigmoid ----------------------------- ###
# 
# svm.sigmoid.basic <- RepCV.SVM(svm.DF, kernel = "sigmoid", nrepeat = 100)
# 
# ### --------------------------- Opitmised Linear --------------------------- ###
# 
# svm.linear.tuning <- PlotTune.tolerance.cost.SVM(svm.DF, "linear", nrepeat = 100)
# 
# svm.linear.tuning[which(svm.linear.tuning[,3] == max(svm.linear.tuning[,3])),]
# svm.linear.opt.params <- list(tolerance = 0.01, cost = 3)
# 
# svm.linear.opt <- RepCV.SVM(svm.DF, kernel = "linear", nrepeat = 100,
#                             tolerance = svm.linear.opt.params$tolerance,
#                             cost = svm.linear.opt.params$cost)
# 
# ### --------------------------- Optimised Radial --------------------------- ###
# 
# svm.radial.tuning <- PlotTune.tolerance.cost.SVM(svm.DF, "radial", nrepeat = 100)
# 
# svm.radial.tuning[which(svm.radial.tuning[,3] == max(svm.radial.tuning[,3])),]
# svm.radial.opt.params <- list(tolerance = 0.02, cost = 10)
# 
# svm.radial.opt <- RepCV.SVM(svm.DF, kernel = "radial", nrepeat = 100,
#                             tolerance = svm.radial.opt.params$tolerance,
#                             cost = svm.radial.opt.params$cost)
# 
# ### ------------------------- Optimised Polynomial ------------------------- ###
# 
# svm.polyn.tuning <- PlotTune.tolerance.cost.SVM(svm.DF, "polynomial", nrepeat = 100)
# 
# svm.polyn.tuning[which(svm.polyn.tuning[,3] == max(svm.polyn.tuning[,3])),]
# svm.polyn.opt.params <- list(tolerance = 0.08, cost = 10)
# 
# svm.polyn.opt <- RepCV.SVM(svm.DF, kernel = "polynomial", nrepeat = 100,
#                            tolerance = svm.polyn.opt.params$tolerance,
#                            cost = svm.polyn.opt.params$cost)
# 
# ### -------------------------- Opitmised Sigmoid --------------------------- ###
# 
# svm.sigmmoid.tuning <- PlotTune.tolerance.cost.SVM(svm.DF, "sigmoid", nrepeat = 100)
# 
# svm.sigmmoid.tuning[which(svm.sigmmoid.tuning[,3] == max(svm.sigmmoid.tuning[,3])),]
# svm.sigmmoid.opt.params <- list(tolerance = 0.08, cost = 1)
# 
# svm.sigmoid.opt <- RepCV.SVM(svm.DF, kernel = "sigmoid", nrepeat = 100,
#                              tolerance = svm.sigmmoid.opt.params$tolerance,
#                              cost = svm.sigmmoid.opt.params$cost)
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# svm.metrics <- list(svm.linear.basic, svm.radial.basic,
#                     svm.polyn.basic, svm.sigmoid.basic,
#                     svm.linear.opt, svm.radial.opt,
#                     svm.polyn.opt, svm.sigmoid.opt)
# 
# names(svm.metrics) <- c("basic.linear", "basic.radial", 
#                            "basic.polyn", "basic.sigmoid", 
#                            "opt.linear", "opt.radial", 
#                            "opt.polyn", "opt.sigmoid")
# 
# saveRDS(svm.metrics, file = "Metrics/SVM/svm.metrics.rda")
# 
# svm.metrics.df <- MakeMetricsDFfromList(svm.metrics.2)
# svm.metrics.df
# saveRDS(svm.metrics.df, file = "Metrics/SVM/svm.metrics.df.rda")

### ======================================================================== ###
### --------------------------- NEURAL NETWORKS ---------------------------- ###
### ======================================================================== ###

nn.DF <- scaled.DF[, -c(1,2,3,4)]

nn.metrics <- readRDS("Metrics/NN/nn.metrics.rda")

PlotOneModelMetrics(nn.metrics$Opt.Tanh, "NN",
                    c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                    title = "NN, Optimised Tanh Model Metrics")

nn.metrics.2 <- nn.metrics
nn.metrics.2[["MNLR"]] <- mnlr.metrics$Standard
nn.metrics.2 <- nn.metrics.2[c(5,1,2,3,4)]
PlotMultipleModelMetrics(nn.metrics.2,
                         title = "NN, All Model Metrics",
                         plot.type="dot")

nn.opt.imps <- nn.metrics$Opt.Logistic$importance
PlotModelFeatImp(nn.opt.imps, type="NN", title="Importance of features in NN, Optimised Logistic model")

# ### --------------------------- Basic Logistic ----------------------------- ###
# 
# nn.logistic.basic <- RepCV.NN(nn.DF, act.fct = "logistic", nrepeat = 100)
# 
# ### ----------------------------- Basic Tanh ------------------------------- ###
# 
# nn.tanh.basic <- RepCV.NN(nn.DF, act.fct = "tanh", nrepeat = 100)
# 
# ### ------------------------- Optimised Logistic --------------------------- ###
# 
# nn.logistic.tuning <- PlotTune.threshold.hidden.NN(nn.DF, act.fct = "logistic", 
#                                                    nrepeat = 25, 
#                                                    thresholds = seq(0.15, 0.27, 0.02),
#                                                    hiddens=seq(3, 9, 2))
# 
# nn.logistic.tuning[which(nn.logistic.tuning[,3] == max(nn.logistic.tuning[,3])),]
# nn.logistic.opt.params <- list(threshold = 0.17, hidden = 9)
# 
# nn.logistic.opt <- RepCV.NN(nn.DF, act.fct = "logistic", nrepeat = 100,
#                                    threshold = nn.logistic.opt.params$threshold,
#                                    hidden = nn.logistic.opt.params$hidden)
# 
# ### ------------------------- Optimised Tanh --------------------------- ###
# 
# nn.tanh.tuning <- PlotTune.threshold.hidden.NN(nn.DF, act.fct = "tanh", 
#                                                nrepeat = 25, 
#                                                thresholds = seq(0.15, 0.27, 0.02),
#                                                hiddens=seq(3, 9, 2))
# 
# nn.tanh.tuning[which(nn.tanh.tuning[,3] == max(nn.tanh.tuning[,3])),]
# nn.tanh.opt.params <- list(threshold = 0.27, hidden = 9)
# 
# nn.tanh.opt <- RepCV.NN(nn.DF, act.fct = "tanh", nrepeat = 100,
#                                threshold = nn.tanh.opt.params$threshold,
#                                hidden = nn.tanh.opt.params$hidden)   
# 
# ### ------------------------------- Overall -------------------------------- ###
# 
# nn.metrics <- list(nn.logistic.basic,
#                    nn.tanh.basic,
#                    nn.logistic.opt,
#                    nn.tanh.opt)
# names(nn.metrics) <- c("Basic.Logistic", "Basic.Tanh", 
#                        "Opt.Logistic", "Opt.Tanh")
# 
# saveRDS(nn.metrics, file = "Metrics/NN/nn.metrics.rda")
# 
# nn.metrics.df <- MakeMetricsDFfromList(nn.metrics.2)
# nn.metrics.df
# saveRDS(nn.metrics.df, file = "Metrics/NN/nn.metrics.df.rda")

### ======================================================================== ###
### ------------------------------ TOP MODELS ------------------------------ ###
### ======================================================================== ###

best.metrics <- list(MNLR = mnlr.metrics$Standard,
                     NB = nb.metrics$Density.Kernel,
                     KNN = knn.metrics$Minkowski,
                     DT = dt.metrics$Gini.Optimised,
                     RF = rf.metrics$Optimised,
                     SVM = svm.metrics$opt.radial,
                     NN = nn.metrics$Opt.Logistic)
best.metrics.df <- MakeMetricsDFfromList(best.metrics)

PlotMultipleModelMetrics(best.metrics, plot.type="dot")

### ======================================================================== ###
### ---------------------------- THE BEST MODEL ---------------------------- ###
### ======================================================================== ###

set.seed(52)

preds <- c()
for (f in YieldFolds(scaled.DF$disease, 5)) {
  test.rows <- f

  test.df <- scaled.DF[test.rows,-c(3,4)]
  train.df <- scaled.DF[(1:446)[-test.rows],-c(1,2,3,4)]
  
  rf.model <- randomForest(as.factor(disease)~., train.df)
  
  preds <- c(preds, predict(rf.model, test.df[,-c(1,2,3)]))
}
preds[as.character(setdiff(1:446, sort(as.numeric(names(preds)))))] <- 3
preds <- c("A", "N", "P")[preds[order(as.numeric(names(preds)))]]

sample.preds <- data.frame(long=scaled.DF$long,
                           lat=scaled.DF$lat,
                           preds = preds,
                           Ground.Truth=scaled.DF$disease,
                           Correct.Prediction=as.character(preds==scaled.DF$disease))

site1.df <- sample.preds[which(sample.preds$lat>2346000), ]
site2.df <- sample.preds[which(sample.preds$lat<2346000), ]

class.m <- confusionMatrix(table(site2.df$preds, site2.df$Ground.Truth))$byClass

PlotOneModelMetrics(class.m, "RF",
                     c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"),
                     title = "Final RF Model")

ggplot(site2.df, aes(x=long, y=lat)) + 
  geom_point(aes(colour=preds, shape=Correct.Prediction), size=6) +  
  geom_point(aes(colour=Ground.Truth, shape=Correct.Prediction), size=4) +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        title = element_text(size=20)) +
  labs(colour = "Disease Class", shape = "Correct Prediction", x = "Longitude", y = "Lattitude",
       title = "Site 2 samples; ground truth and RF prediction",
       )
