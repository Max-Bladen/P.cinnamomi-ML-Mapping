
### ------------------------------------------------------------------------ ###
### ------------------------------- LIBRARIES ------------------------------ ###
### ------------------------------------------------------------------------ ###

library(arules)
library(biganalytics)
library(caret) 
library(class)
library(DescTools)
library(dplyr)
library(e1071)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(kknn)
library(mlogit)
library(naivebayes)
library(neuralnet)
library(NeuralNetTools)
library(nnet)
library(pheatmap)
library(plotly)
library(randomForest)
library(raster)
library(readxl)
library(rpart)
library(stringi)
library(tidyverse)

### ======================================================================== ###
### ------------------------------- FUNCTIONS ------------------------------ ###
### ======================================================================== ###

### ============================ Initialisation ============================ ###

### ------------------------------------------------------------------------ ###
## total: total number of models to be generated. nfold*nrepeat
#
# Generates series of empty dataframes. This structure is used by all RepCV.* 
# methods. The columns in each correspond to performance metrics from
# 'confusionMatrix()'. DFs for each individual model (no averaging), per class
# averages and total, weighted (by class) averages
#
## return: list of three empty dataframes 
InitialiseMetricsObject <- function(total) {
  
  metric.names <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", 
                     "Precision", "Recall", "F1", "Prevalence", "Detection Rate", 
                     "Detection Prevalence", "Balanced Accuracy")
  metrics <- list()
  
  metrics[["weighted.avg.all"]] <- data.frame(matrix(0, ncol = length(metric.names), 
                                                   nrow = total))
  metrics[["per.class.total"]] <- data.frame(matrix(0, ncol = length(metric.names), 
                                                  nrow = 3))
  metrics[["weighted.avg.total"]] <- data.frame(matrix(0, ncol = length(metric.names), 
                                                     nrow = 1))
  
  colnames(metrics[["weighted.avg.all"]]) <- metric.names
  colnames(metrics[["per.class.total"]]) <- metric.names
  colnames(metrics[["weighted.avg.total"]]) <-  metric.names
    
  
  rownames(metrics[["per.class.total"]]) <- c("Class: A", "Class: N", "Class: P")
  rownames(metrics[["weighted.avg.total"]]) <- c("All Classes")
  
  return(metrics)
}


### ------------------------------------------------------------------------ ###
## df: df of stability values from 'RepCV.MNLR().'
#
## return: empty stability dataframe. cols correspond to features and rows
##         correspond to class
InitialiseStabilityDF <- function(df) {
  cols <- c("(Intercept)", colnames(df))
  stability <- data.frame(matrix(0, nrow = 3, ncol=length(cols)))
  colnames(stability) <- cols
  rownames(stability) <- c("A", "N", "P")
  
  return(stability)
}


### ------------------------------------------------------------------------ ###
## cols: list of feature names
## byClass: determines if output has 3 rows (for each class) or one row (for
##          weighted averages)
#
## return: empty stability dataframe. cols correspond to features and rows
##         correspond to class
InitialiseImportanceDF <- function(cols, 
                                   byClass=F) {
  importance <- data.frame(matrix(0, nrow = ifelse(byClass, 3, 1), ncol=length(cols)))
  colnames(importance) <- cols
  if (byClass) {rownames(importance) <- c("A", "N", "P")}
  else {rownames(importance) <- c("Importance")}
  
  return(importance)
}

### ======================================================================== ###

### ================================ Tuning ================================ ###
### ALL TUNING FUNCTION PARAMETERS:
## df: dataframe of disease class and predictors
## nrepeat: number of iterations at each HP grid position
## metrics: what metrics to plot. Single-HP tuning functions can use multiple metrics.
##          multi-HP tuning functions can only use one, defaults to BA.
#
# Iterates over grid of HP values. At each set of HPs, 5-folds nrepeat CV applied.
# The specified metrics are extracted for each of these models and averaged, yielding 
# metric scores for each set of HPs. Tuning curve or surface is then plotted
#
## return: the means and standard deviations of metric scores at each unique HP set
#
### ------------------------------------------------------------------------ ###
## kmax: maximum k considered. number of clusters
#
# Single-HP Tuning function
PlotTune.k.kNN <- function(df,
                           kmax = 10,
                           dist="euclidean",
                           nrepeat=50,
                           metrics = c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy"),
                           error.bars = F,
                           ylab = "metric") 
  {
  
  tuning.means <- data.frame(matrix(0, nrow = kmax, ncol = length(metrics)+1))
  tuning.sds <- data.frame(matrix(0, nrow = kmax, ncol = length(metrics)+1))

  for (k in 1:kmax) {
    cat("\rCalcuating metrics for k =", k,"...")
    
    model <- RepCV.kNN(df, k = k, dist = dist, nrepeat = nrepeat,
                                 nfolds = 5, silent=T)
    
    tuning.means[k, ] <- c(k, model$weighted.avg.total[metrics])
    
    tuning.sds[k, ] <- c(k, unlist(lapply(model$weighted.avg.all, sd))[metrics])
  }
  
  metrics <- c("k", metrics)
  
  cat("\n")
  
  colnames(tuning.means) <- metrics
  colnames(tuning.sds) <- metrics
  
  tuning.df <- tuning.means %>% 
    dplyr::select(all_of(metrics)) %>% 
    gather(key = "variable", value = "mean", -k)
  
  tuning.sds <- tuning.sds %>% 
    dplyr::select(all_of(metrics)) %>% 
    gather(key = "variable", value = "sd", -k)
  
  tuning.df[,"sd"] <- tuning.sds[,"sd"]
  
  p <- ggplot(tuning.df, aes(x = k, y = mean)) + 
    geom_line(aes(color = variable), size=1) + 
    geom_point(aes(color = variable), size=2) +
    scale_color_manual(values = c("darkred", "steelblue", "magenta1", "seagreen", "tomato")) + 
    scale_x_continuous(breaks=1:10) + 
    ylab(ylab)
  
  if (error.bars) { p <- p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = variable, width = 0.2))}
  
  print(p)
  
  return(tuning.df)
}


### ------------------------------------------------------------------------ ###
## cps: range of cp values to test. cost penalisation of complex tree
## minsplits: range of minsplit values to test. number of samples required for split
## split.crit: split criterion to use for DT building
#
# Multi-HP Tuning function
PlotTune.cp.minsplit.DT <- function(df,
                                    cps = seq(0, 0.2, 0.02),
                                    minsplits = seq(10, 30, 5),
                                    split.crit = "gini",
                                    nrepeat=50,
                                    metric = "Balanced Accuracy")
{
  
  out <- data.frame(matrix(0, nrow = length(cps)*length(minsplits), ncol = 3))
  value.map <- matrix(0, nrow = length(cps), ncol = length(minsplits))
  
  if (length(metric) > 1) {stop("Can only visualise one metric in 3D!")}
  
  i <- 1
  for (x in 1:length(cps)) {
    cp = cps[x]
    
    cat("\rCalcuating metrics for cp =", cp, "...")
    
    for (y in 1:length(minsplits)) {
      
      minsplit = minsplits[y]
      
      model <- RepCV.DecisionTree(df, cp=cp, minsplit = minsplit,
                                  split.crit = split.crit, 
                                  nrepeat = nrepeat,
                                  nfolds = 5, silent=T)
      
      out[i,] <- c(cp, minsplit, model$metrics$weighted.avg.total[[metric]])
      
      value.map[x, y] <- model$metrics$weighted.avg.total[[metric]]
      
      i <- i + 1
    }
    
  }
  
  colnames(out) <- c("cp", "minsplit", metric)
  
  axx <- list(
    ticketmode = 'array',
    ticktext = as.character(minsplits),
    tickvals = 1:length(minsplits)-1,
    range = 1:length(minsplits)-1,
    title = "minimum node size for split"
  )

  axy <- list(
    ticketmode = 'array',
    ticktext = as.character(cps),
    tickvals = 1:length(cps)-1,
    range = 1:length(cps)-1,
    title = "complexity penalisation"
  )
  
  axz <- list(
    title = metric
  )
  
  p <- plot_ly(z = value.map, type = "surface")
  
  p <- p %>%
    layout(scene= list(xaxis=axx, yaxis=axy, zaxis=axz))
  
  
  
  print(p)
  
  return(out)
}


### ------------------------------------------------------------------------ ###
## cps: range of cp (cost penalisation) values to test. cost penalisation of complex tree
## split.crit: split criterion to use for DT building
#
# Single-HP Tuning function
PlotTune.cp.DT <- function(df,
                           cps = seq(0, 0.2, 0.02),
                           split.crit = "gini",
                           nrepeat=50,
                           metrics = c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy"),
                           error.bars = F,
                           ylab = "metric")
{

  tuning.means <- data.frame(matrix(0, nrow = length(cps), ncol = length(metrics)+1))
  tuning.sds <- data.frame(matrix(0, nrow = length(cps), ncol = length(metrics)+1))

  for (i in 1:length(cps)) {
    cp = cps[i]
    
    cat("\rCalcuating metrics for cp =", cp,"...")
    model <- RepCV.DecisionTree(df, cp=cp, 
                                split.crit = split.crit, 
                                nrepeat = nrepeat,
                                nfolds = 5, silent=T)

    tuning.means[i, ] <- c(cp, model$metrics$weighted.avg.total[metrics])

    tuning.sds[i, ] <- c(cp, unlist(lapply(model$metrics$weighted.avg.all, sd))[metrics])
  }

  metrics <- c("cp", metrics)

  cat("\n")

  colnames(tuning.means) <- metrics
  colnames(tuning.sds) <- metrics

  tuning.df <- tuning.means %>%
    dplyr::select(all_of(metrics)) %>%
    gather(key = "variable", value = "mean", -cp)

  tuning.sds <- tuning.sds %>%
    dplyr::select(all_of(metrics)) %>%
    gather(key = "variable", value = "sd", -cp)

  tuning.df[,"sd"] <- tuning.sds[,"sd"]


  p <- ggplot(tuning.df, aes(x = cp, y = mean)) +
    geom_line(aes(color = variable), size=1) +
    geom_point(aes(color = variable), size=2) +
    scale_color_manual(values = c("darkred", "steelblue", "magenta1", "seagreen", "tomato")) +
    scale_x_continuous(n.breaks=length(cps)) +
    ylab(ylab)

  if (error.bars) { p <- p + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = variable, width = 0.2))}

  print(p)

  return(tuning.df)
}


### ------------------------------------------------------------------------ ###
## ntrees: range of ntree values to test. number of trees built
## mtrys: range of mtry values to test. number of features considered at each split
#
# Multi-HP Tuning function
PlotTune.mtry.ntree.RF <- function(df,
                                   ntrees = seq(100, 900, 100),
                                   mtrys = seq(2, 8, 2),
                                   nrepeat=50,
                                   metric = "Balanced Accuracy")
{
  
  out <- data.frame(matrix(0, nrow = length(ntrees)*length(mtrys), ncol = 3))
  value.map <- matrix(0, nrow = length(ntrees), ncol = length(mtrys))
  
  if (length(metric) > 1) {stop("Can only visualise one metric in 3D!")}
  
  i <- 1
  for (x in 1:length(ntrees)) {
    ntree = ntrees[x]
    
    for (y in 1:length(mtrys)) {
      
      mtry = mtrys[y]
      
      cat("\rCalcuating metrics for ntree =", ntree," and mtry =", mtry)
      
      model <- RepCV.RandomForest(df, ntree=ntree, mtry = mtry,
                                  nrepeat = nrepeat,
                                  nfolds = 5, silent=T)
      
      out[i,] <- c(ntree, mtry, model$metrics$weighted.avg.total[[metric]])
      
      value.map[x, y] <- model$metrics$weighted.avg.total[[metric]]
      
      i <- i + 1
    }
    
  }
  
  colnames(out) <- c("ntrees", "mtrys", metric)
  
  axx <- list(
    ticketmode = 'array',
    ticktext = as.character(mtrys),
    tickvals = 1:length(mtrys)-1,
    range = 1:length(mtrys)-1,
    title = "features allow at split"
  )
  
  axy <- list(
    ticketmode = 'array',
    ticktext = as.character(ntrees),
    tickvals = 1:length(ntrees)-1,
    range = 1:length(ntrees)-1,
    title = "num. trees"
  )
  
  axz <- list(
    title = metric
  )
  
  p <- plot_ly(z = value.map, type = "surface")
  
  p <- p %>%
    layout(scene= list(xaxis=axx, yaxis=axy, zaxis=axz))
  
  print(p)
  
  return(out)
}


### ------------------------------------------------------------------------ ###
## tolerances: range of tolerance values to test. allowance for misclassification
## costs: range of cost values to test. penalisation of samples within margin
## kernel: kernel transformation to use for SVM building
#
# Multi-HP Tuning function
PlotTune.tolerance.cost.SVM <- function(df, 
                                        kernel="radial",
                                        tolerances = seq(0.002, 0.01, 0.002),
                                        costs = seq(1, 10, 1),
                                        nrepeat=50,
                                        metric = "Balanced Accuracy")
{
  
  out <- data.frame(matrix(0, nrow = length(tolerances)*length(costs), ncol = 3))
  value.map <- matrix(0, nrow = length(tolerances), ncol = length(costs))
  
  if (length(metric) > 1) {stop("Can only visualise one metric in 3D!")}
  
  i <- 1
  for (x in 1:length(tolerances)) {
    tol = tolerances[x]
    
    cat("\rCalcuating metrics for tol =", tol,"...")
    
    for (y in 1:length(costs)) {
      
      cost = costs[y]
      
      model <- RepCV.SVM(df, tolerance = tol, cost=cost,
                         kernel=kernel,
                         nrepeat = nrepeat,
                         nfolds = 5, silent=T)
      
      out[i,] <- c(tol, cost, model$metrics$weighted.avg.total[[metric]])
      
      value.map[x, y] <- model$metrics$weighted.avg.total[[metric]]
      
      i <- i + 1
    }
    
  }
  
  colnames(out) <- c("tolerances", "costs", metric)
  
  axx <- list(
    ticketmode = 'array',
    ticktext = as.character(costs),
    tickvals = 1:length(costs)-1,
    range = 1:length(costs)-1,
    title = "cost"
  )
  
  axy <- list(
    ticketmode = 'array',
    ticktext = as.character(tolerances),
    tickvals = 1:length(tolerances)-1,
    range = 1:length(tolerances)-1,
    title = "tolerance"
  )
  
  axz <- list(
    title = metric
  )
  
  p <- plot_ly(z = value.map, type = "surface")
  
  p <- p %>%
    layout(scene= list(xaxis=axx, yaxis=axy, zaxis=axz))
  
  print(p)
  
  return(out)
}


### ------------------------------------------------------------------------ ###
## hiddens: range of hidden values to test. number of hidden neurons
## thresholds: range of threshold values to test. act.fct threshold
## act.fct: activation function for neurons in NN
#
# Multi-HP Tuning function
PlotTune.threshold.hidden.NN <- function(df,
                                         act.fct="logistic",
                                         hiddens=seq(5, 15, 2),
                                         thresholds=seq(0.1, 2, 0.2),
                                         nrepeat=50,
                                         metric = "Balanced Accuracy")
{
  
  out <- data.frame(matrix(0, nrow = length(thresholds)*length(hiddens), ncol = 3))
  value.map <- matrix(0, nrow = length(thresholds), ncol = length(hiddens))
  
  if (length(metric) > 1) {stop("Can only visualise one metric in 3D!")}
  
  i <- 1
  for (x in 1:length(thresholds)) {
    threshold = thresholds[x]
    
    
    for (y in 1:length(hiddens)) {
      
      hidden = hiddens[y]
      
      cat("\rCalcuating metrics for threshold =", threshold," and hidden =", hidden)
      
      model <- RepCV.NeuralNet(df, threshold = threshold, hidden=hidden,
                               act.fct=act.fct,
                               nrepeat = nrepeat,
                               nfolds = 5, silent=F)
      
      out[i,] <- c(threshold, hidden, model$metrics$weighted.avg.total[[metric]])
      
      value.map[x, y] <- model$metrics$weighted.avg.total[[metric]]
      
      i <- i + 1
    }
    
  }
  
  colnames(out) <- c("tolerances", "hiddens", metric)
  
  axx <- list(
    ticketmode = 'array',
    ticktext = as.character(hiddens),
    tickvals = 1:length(hiddens)-1,
    range = 1:length(hiddens)-1,
    title = "hidden neurons"
  )
  
  axy <- list(
    ticketmode = 'array',
    ticktext = as.character(thresholds),
    tickvals = 1:length(thresholds)-1,
    range = 1:length(thresholds)-1,
    title = "threshold"
  )
  
  axz <- list(
    title = metric
  )
  
  p <- plot_ly(z = value.map, type = "surface")
  
  p <- p %>%
    layout(scene= list(xaxis=axx, yaxis=axy, zaxis=axz))
  
  print(p)
  
  return(out)
}

### ======================================================================== ###

### ================================ Rep.CV ================================ ###
### ALL RepCV FUNCTION PARAMETERS
## df: dataframe of disease class and predictors
## nfolds: number of folds of data to use at each nrepeat
## nrepeat: number of iterations of model building
## seed: used to control for RNG in folds generation. reproducibility
## silent: if TRUE, then progress tracker won't be printed
#
# For each repeat, nfolds are generated. On each of these, nfold-1 used for training
# and the remaining fold used for test. Metrics (and importances) extracted and 
# appended/added to metrics dataframes. uses specified HPs
#
## return: list of metrics type object for all models. if stability/importance available, 
##         list of metrics and importance measures for all models
#
### ------------------------------------------------------------------------ ###
## reg.alpha: regularisation alpha. If NULL, no regularisation. Otherwise, employs
##            elastic net regularisation using alpha
#
RepCV.MNLR <- function(df,
                       reg.alpha=NULL,
                       nrepeat=1,
                       nfolds=5,
                       seed=NULL,
                       silent=F) {
  
  # set up the stability assessment df
  stability <- InitialiseStabilityDF(df[,-1])
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  m <- 1
  
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      for (ref in sort(unique(df$disease))) {
        tr <- train.data
        tr$disease <- ifelse(tr$disease == ref, 1, 0)
        
        te <- test.data
        te$disease <- ifelse(te$disease == ref, 1, 0)
        
        if (is.null(reg.alpha)) {
          model <- glm(disease~., data=tr, family="binomial")
          preds <- predict(model, newdata = te[-1], type="response")
        }
        else {
          cv.lasso <- cv.glmnet(model.matrix(disease~., tr)[,-1], tr[,1], alpha = reg.alpha, family = "binomial")
          model <- glmnet(model.matrix(disease~., tr)[,-1], tr[,1], family="binomial", alpha = reg.alpha, lambda = cv.lasso$lambda.1se)
          preds <- predict(model, newx = model.matrix(disease~., te)[,-1], type="response")
        }
        
        
        preds.class <- ifelse(preds > 0.5, 1, 0)
        test.pred.table <- table(te[,1], preds.class)
        
        if (!("0" %in% colnames(test.pred.table))) {
          test.pred.table <- cbind(c(0,0), test.pred.table)
          colnames(test.pred.table) <- c("0", "1")
        }
        if (!("1" %in% colnames(test.pred.table))) {
          test.pred.table <- cbind(test.pred.table, c(0,0))
          colnames(test.pred.table) <- c("0", "1")
        }
        
        conf.mat <- confusionMatrix(test.pred.table)$byClass
        conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
        
        metrics$per.class.total[paste0("Class: ", ref),] <- metrics$per.class.total[paste0("Class: ", ref),] + conf.mat/(nrepeat*nfolds)
        
        metrics$weighted.avg.all[m,] <- metrics$weighted.avg.all[m,] + conf.mat * (1/length(unique(df$disease)))
        
        if (is.null(reg.alpha)) {
          # get p-values for coefficients:
          p <- coef(summary(model))[,4]
          stability[ref, which(p<0.05)] <- stability[ref,which(p<0.05)] + 1/(nrepeat*nfolds)
        }
        else {
          sel.feats <- rownames(coef(model))[which(coef(model)[-1]!=0)+1]
          stability[ref, sel.feats] <- stability[ref, sel.feats] + 1/(nrepeat*nfolds)
        }
      }
    
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
    
  }
  metrics$weighted.avg.total <- WeightedAverage(df$disease, metrics$per.class.total)
  
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              stability = stability))
}


### ------------------------------------------------------------------------ ###
## discretise: if NULL, predictors left as continuous variables. If 'frequency'
##             or 'interval', will make equal frequency or width bins of each predictor
## kernel: if TRUE, kernel density function will be used. discretise must equal NULL
#
RepCV.NB <- function(df, 
                     discretise = NULL,
                     kernel = F,
                     nrepeat=1, 
                     nfolds=5,
                     seed=NULL,
                     silent=F) {
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!is.null(discretise)) {
    for (c in 2:ncol(df)) {
      df[,c] <- discretize(df[,c], breaks = 5, labels = F, method = discretise) # by frequency
    }
  }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- naive_bayes(disease ~ ., data = train.data, usekernel = kernel)
      
      preds <- predict(model, test.data[,-1])
      truth <- test.data[,1]
      test.pred.table <- table(truth, preds)
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
        
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              importance = NULL))
}


### ------------------------------------------------------------------------ ###
## k: number of clusters
## dist: determines what distance measure model will use. Can be "euclidean", 
##       "manhattan" or "minkowski"
#
RepCV.kNN <- function(df, 
                      k=5,
                      dist="euclidean",
                      nrepeat=1, 
                      nfolds=5,
                      seed=NULL,
                      silent=F) {
  
  distances <- c("manhattan", "euclidean", "minkowski")
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (dist %in% distances) {d <- which(dist == distances)}
  else {stop("Invalid Distance")}
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- kknn(as.factor(disease)~., train.data, test.data, k = k, distance = d)
      preds <- model$fitted.values
      truth <- test.data[,1]
      test.pred.table <- table(truth, preds)
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
      
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              importance = NULL))
}


### ------------------------------------------------------------------------ ###
## cp: cost penalisation of complex tree
## minsplit: minimum samples required for split at a node
## split.crit: split criterion used 
#
RepCV.DT <- function(df, 
                     cp = 0.01,
                     minsplit = 20,
                     split.crit = "gini",
                     nrepeat=1, 
                     nfolds=5,
                     seed=NULL,
                     silent=F) {
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  # set up the importance objects
  importance <- InitialiseImportanceDF(colnames(df)[-1])
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  over.fitting <- c()
  
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- rpart(as.factor(disease)~., train.data, 
                     control = rpart.control(cp=cp, minsplit=minsplit),
                     parms = list(split=split.crit))
      preds <- predict(model, test.data[,-1], "class")
      truth <- test.data[,1]
      train.pred.table <- table(train.data[,1], predict(model, train.data[,-1], "class"))
      test.pred.table <- table(truth, preds)
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
      
      train.acc <- WeightedAverage(df$disease, confusionMatrix(train.pred.table)$byClass)[[11]]
      test.acc <- WeightedAverage(df$disease, confusionMatrix(test.pred.table)$byClass)[[11]]
      
      over.fitting <- c(over.fitting, test.acc-train.acc)
      
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      model.var.imp <- model$variable.importance/sum(model$variable.importance) * 100
      
      for (feat in names(model.var.imp)) {
        importance[feat] <- importance[feat] + model.var.imp[feat]/(nrepeat*nfolds)
      }
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              importance = importance,
              over.fitting = over.fitting))
}


### ------------------------------------------------------------------------ ###
## ntree: number of trees to build
## mtry: how many features are to be considered at each split
#
RepCV.RF <- function(df, 
                     ntree=500,
                     mtry=4,
                     nrepeat=1, 
                     nfolds=5,
                     seed=NULL,
                     silent=F) {
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  # set up the importance objects
  importance <- InitialiseImportanceDF(colnames(df)[-1], T)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  all.preds <- data.frame(matrix(NA, ncol=1, nrow=nrow(df)))
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- randomForest(as.factor(disease)~., train.data, importance=T,
                            ntree = ntree, mtry = mtry)
      preds <- predict(model, test.data[,-1], "class")
      truth <- test.data[,1]
      test.pred.table <- table(truth, preds)
      
      all.preds[f, ] <- preds
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
      
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      model.var.imp <- t(model$importance[,1:3])
      model.var.imp <- model.var.imp/rowSums(model.var.imp)*100
      
      importance <- importance + model.var.imp/(nrepeat*nfolds)
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              importance = importance,
              preds=all.preds))
}


### ------------------------------------------------------------------------ ###
## kernel: kernel applied to data to transform it
## cost: penalisation of samples within the margin
## tolerance: allowance of misclassifications
#
RepCV.SVM <- function(df, 
                      kernel="radial",
                      cost=1,
                      tolerance=0.001,
                      nrepeat=1, 
                      nfolds=5,
                      seed=NULL,
                      silent=F) {
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- svm(as.factor(disease)~., train.data, kernel = kernel,
                   cost = cost, tolerance = tolerance)
      preds <- predict(model, test.data[,-1])
      truth <- test.data[,1]
      test.pred.table <- table(truth, preds)
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
      
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics))
}


### ------------------------------------------------------------------------ ###
## act.fct: activation function
## hidden: number of hidden neurons
## threshold: threshold of activation function
# 
RepCV.NN <- function(df,
                     act.fct="logistic",
                     hidden=11,
                     threshold=0.1,
                     nrepeat=1, 
                     nfolds=5,
                     seed=NULL,
                     silent=F) {
  
  # set up the classification metrics object
  metrics <- InitialiseMetricsObject(nfolds*nrepeat)
  
  importance <- InitialiseImportanceDF(colnames(df)[-1], F)
  
  if (!is.null(seed)) { set.seed(seed) }
  
  if (!silent) {cat("\rProgress: ", 0, "%")}
  
  m <- 1
  for (i in 1:nrepeat) {
    
    folds <- YieldFolds(df$disease, nfolds)
    
    for (f in folds) {
      train.data <- df[(1:nrow(DF))[-f],]
      test.data <- df[f,]
      
      model <- neuralnet(as.factor(disease)~., train.data, 
                         hidden = hidden, act.fct = act.fct,
                         threshold = threshold)
      preds <- predict(model, test.data[,-1])
      preds <- c("A", "N", "P")[apply(preds, 1, which.max)]
      truth <- test.data[,1]
      test.pred.table <- table(truth, preds)
      
      CompleteConfTable <- function(table) {
        for (c in c("A", "N", "P")) { 
          if (!(c %in% colnames(table))) {
            table <- cbind(table, c(0,0,0))
            colnames(table)[length(colnames(table))] <- c
          }
        }
        table <- table[, c("A", "N", "P")]
        return(table)
      }
      
      if (length(colnames(test.pred.table))!=3) {test.pred.table <- CompleteConfTable(test.pred.table)}
      
      conf.mat <- confusionMatrix(test.pred.table)$byClass
      conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
      
      metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
      
      avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
      metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
      metrics$weighted.avg.all[m,] <- avgd.fold.metrics
      
      imps <- t(suppressWarnings(olden(model, bar_plot=F)))
      imps <- imps/sum(abs(imps))*100
      
      importance <- importance + imps/(nrepeat*nfolds)
      
      m <- m+1
      if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
    }
  }
  if (!silent) {cat("\n")}
  
  return(list(metrics = metrics,
              importance = importance))
}

### ------------------------------------------------------------------------ ###
## params: HP values to use for each classifier
#
RepCV.EnsembleLearner <- function(df, 
                                  nfolds=5, 
                                  nrepeat=1,
                                  params,
                                  silent=F) {
        
        # set up the classification metrics object
        metrics <- InitialiseMetricsObject(nfolds*nrepeat)
        
        folds <- YieldFolds(enbl.DF$disease, 5)
        
        if (!silent) {cat("\rProgress: ", 0, "%")}
        
        m <- 1
        for (i in 1:nrepeat) {
          
          for (f in folds) {
            train.data <- df[(1:nrow(DF))[-f],]
            test.data <- df[f,]
            
            nb.model <- naive_bayes(disease ~ ., data = train.data,
                                    kernel = params$nb$kernel)
            
            if (params$knn$dist %in% c("manhattan", "euclidean", "minkowski")) 
              {d <- which(params$knn$dist == c("manhattan", "euclidean", "minkowski"))}
            knn.model <- kknn(as.factor(disease)~., train.data, test.data,
                              k = params$knn$k,
                              distance = d)
            
            dt.model <- rpart(as.factor(disease)~., train.data, 
                              control = rpart.control(cp=params$dt$cp, minsplit=params$dt$minsplit),
                              parms = list(split=params$dt$split.crit))
            
            rf.model <- randomForest(as.factor(disease)~., train.data,
                                     ntrees = params$rf$ntrees,
                                     mtry = params$rf$mtry)
            
            svm.model <- svm(as.factor(disease)~., train.data,
                             kernel = params$svm$kernel,
                             cost = params$svm$cost,
                             tolerance = params$svm$tolerance)
            
            nn.model <- neuralnet(as.factor(disease)~., train.data,
                                  act.fct = params$nn$act.fct,
                                  hidden = params$nn$hidden,
                                  threshold = params$nn$threshold)
            
            
            
            nb.preds <- predict(nb.model, test.data[,-1])
            knn.preds <- knn.model$fitted.values
            dt.preds <- unname(predict(dt.model, test.data[,-1], "class"))
            rf.preds <- unname(predict(rf.model, test.data[,-1], "class"))
            svm.preds <- unname(predict(svm.model, test.data[,-1]))
            nn.preds <- as.factor(c("A", "N", "P")[apply(predict(nn.model, test.data[,-1]), 1, which.max)])
            
            all.preds <- cbind(nb.preds,knn.preds,dt.preds,rf.preds,svm.preds,nn.preds) # 
            
            enbl.preds <- c()
            for (r in 1:nrow(all.preds)) {
              row.pred <- as.integer(names(which(table(all.preds[r,])==max(table(all.preds[r,])))))
              if (length(row.pred)!=1) {
                row.pred <- "P"
              }
              enbl.preds <- c(enbl.preds, c("A", "N", "P")[row.pred])
            }
            
            test.pred.table <- table(test.data[,1], enbl.preds)
            
            conf.mat <- confusionMatrix(test.pred.table)$byClass
            conf.mat[which(is.na(conf.mat), arr.ind=T)] <- 0
            
            metrics$per.class.total <- metrics$per.class.total + conf.mat/(nrepeat*nfolds)
            
            avgd.fold.metrics <- WeightedAverage(df$disease, conf.mat)
            metrics$weighted.avg.total <- metrics$weighted.avg.total + avgd.fold.metrics/(nrepeat*nfolds)
            metrics$weighted.avg.all[m,] <- avgd.fold.metrics
            
            m <- m+1
            if (!silent) {cat("\rProgress: ", (m-1)/(nrepeat*nfolds)*100, "%")}
          }
        }
        
        
        return(metrics)
      }

### ======================================================================== ###

### =============================== Plotting =============================== ###


### ------------------------------------------------------------------------ ###
## metrics: metrics type object for a single model
## model.type: type of model to be used. "MNLR", "NB", "kNN", "DT", "RF", "SVM" or "NN"
## d.m: desired metric(s) to be displayed
## title: figure title
#
# Creates a bar graph depicting the performance metrics of a model, separated by
# class.
#
## return: ggplot2 object
PlotOneModelMetrics <- function(metrics, 
                                model.type, 
                                d.m = c("Balanced Accuracy"), 
                                title = "placeholder"
                                ) {
  
  if (is.null(model.type)) {stop("need model type")}
  if (!(model.type %in% MODEL.TYPES)) {stop("invalid model type")}
  
  per.class <- metrics$metrics$per.class.total[c(1,3,2), ] # rearrange of rows to fit CLASSES
  rownames(per.class) <- CLASSES
  
  vals <- c()
  
  for (m in d.m) {
    vals <- c(vals, per.class[, m])
  }
  
  df <- data.frame(Class = rep(CLASSES, length(d.m)),
                   Metric = as.vector(sapply(d.m, function(x) {rep(x, length(CLASSES))})),
                   Value = vals)
  
  p <- ggplot(data=df, aes(x=Class, y=Value, fill=Metric)) +
      geom_bar(stat='identity', position=position_dodge()) +
      labs(title = title) +
      ylab("Metric Scores") +
      xlab("Metric & Class")

  p
}


### ------------------------------------------------------------------------ ###
## metrics.list: list of metrics type object
## d.m: desired metric(s) to be displayed
## title: figure title
## plot.type:  controls the type of plot to produce. can be "bar" or "dot"
# 
# Produces a figure which allows the visual and numerical comparison of models
# via depictions of metric means and standard deviations. 
#
## return: ggplot2 object
PlotMultipleModelMetrics <- function(metrics.list,
                                     d.m = c("Specificity", "Sensitivity", "Precision", "F1", "Balanced Accuracy"), # Desired Metrics to be shown,
                                     title = "placeholder",
                                     plot.type = "dot"
) {
  
  
  
  df <- data.frame(matrix(0, nrow = length(metrics.list)*length(d.m), ncol = 4))
  colnames(df) <- c("Model", "Metric", "Mean", "SD")
  i <- 1
  for (m in 1:length(metrics.list)) {
    
    m.name <- names(metrics.list)[m]
    
    
    all.runs <- metrics.list[[m]]$metrics$weighted.avg.all
    
    for (metric in d.m) {

      df[i, ]$Model <- m.name
      df[i, ]$Metric <- metric
      df[i, ]$Mean <- metrics.list[[m]]$metrics$weighted.avg.total[[metric]]#mean(all.runs[, metric])
      df[i, ]$SD <- sd(all.runs[, metric])
      
      i<-i+1
    }
  }
  
  if (plot.type=="bar") {
    p <- ggplot(data=df, aes(x=factor(Model, levels = names(metrics.list)),
                                       y=Mean, fill=Metric)) +
      geom_bar(stat='identity', position=position_dodge()) +
      labs(title = title) +
      ylab("Metric Scores") +
      xlab("Model") + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, width = 0.9),
                  position=position_dodge())
  }
  else if (plot.type=="dot") {
    mean.range <- c(min(df$Mean),max(df$Mean))
    p <- ggplot(df, aes(factor(Model, levels = names(metrics.list)), 
                        forcats::fct_rev(Metric), 
                        fill = Mean, 
                        size = SD)) +
        geom_point(shape = 21, stroke = 0) +
        geom_hline(yintercept = seq(.5, 4.5, 1), size = .2) +
        geom_text(label=as.character(round(df$Mean, 3)), size=7, nudge_y = -.35) +
        scale_x_discrete(position = "top") +
        scale_radius(range = c(10, 25)) +
        scale_fill_gradient(low = "orange", high = "blue", breaks = mean.range, 
                            labels = as.character(round(mean.range, 3)), limits = mean.range) + 
        theme_minimal() +
        theme(legend.position = "right", 
              panel.grid.major = element_blank(),
              legend.text = element_text(size = 18),
              legend.title = element_text(size = 20),
              axis.text = element_text(size = 22), 
              axis.text.x = element_text(angle=20)) +
        guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                                   label.position = "left",
                                   title.position = "top", 
                                   order = 1),
               fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
        labs(size = "Standard Deviation", fill = "Mean:", x = NULL, y = NULL)
  }
  
  

  p
}


### ------------------------------------------------------------------------ ###
## w.avg.imps: weighted average importances. importance object returned from RepCV.*()
## d: disease class vector
## type: type of model to be used. "MNLR", "NB", "kNN", "DT", "RF", "SVM" or "NN"
## title: figure title
#
# Produces a bar plot of all feature importance values
#
## return: ggplot2 object
PlotModelFeatImp <- function(w.avg.imps, 
                             d=NULL, 
                             type, 
                             title='placeholder') {
  
  if (type=="MNLR" && !is.null(d)) {
    
    w.avg.imps <- WeightedAverage(d, w.avg.imps)*100
  }
    
  w.avg.imps <- data.frame(Feature=unlist(names(w.avg.imps)),
                           W.Avg=unlist(unname(w.avg.imps)))
  w.avg.imps <- w.avg.imps[order(-w.avg.imps$W.Avg), ]
  
  if (type=="MNLR") {
    
    ylab = "Weighted Average of Stability (%)"
    bar.labs <- paste0(round(w.avg.imps$W.Avg, 1), "%")
    
  }
  else if (type %in% c("DT","RF")) {
    
    ylab = "Weighted Average of Importance (unitless)"
    bar.labs <- round(w.avg.imps$W.Avg, 2)
  }
  
  else if(type=="NN") {
    ylab = "Average of Olden Measure (unitless)"
    bar.labs <- round(w.avg.imps$W.Avg, 2)
  }
  
  p <- ggplot(data=w.avg.imps,
    aes(x=Feature,
                    y=W.Avg)) +
          geom_bar(stat='identity') +
      scale_x_discrete(limits=w.avg.imps$Feature) +
          labs(title = title) +
          ylab(ylab) +
      theme(axis.text.y = element_text(size=18),
            axis.text.x = element_text(size=22, angle=45),
            title = element_text(size=18),
            axis.title = element_text(size = 20))+
      geom_text(aes(label=bar.labs), vjust=-0.3, size=4.5)
    
  p
}

### ======================================================================== ###

### =============================== Helpers ================================ ###


### ------------------------------------------------------------------------ ###
## metrics.list: list of metrics type objects
#
## return: df  containing weighted average means and sds of each metric
MakeMetricsDFfromList <- function(metrics.list) {
  means <- data.frame(matrix(NA, nrow = length(metrics.list), ncol = 5))
  sds <- data.frame(matrix(NA, nrow = length(metrics.list), ncol = 5))
  for (i in 1:length(metrics.list)) {
    means[i,] <- metrics.list[[i]]$metrics$weighted.avg.total[c(1,2,5,7,11)]
    sds[i,] <- apply(metrics.list[[i]]$metrics$weighted.avg.all[c(1,2,5,7,11)],2,sd)
  }
  colnames(means) <- names(metrics.list[[1]]$metrics$weighted.avg.total)[c(1,2,5,7,11)]
  rownames(means) <- names(metrics.list)
  colnames(sds) <- colnames(means)
  rownames(sds) <- rownames(means)
  
  
  return(list(m=means,
              sd=sds))
}


### ------------------------------------------------------------------------ ###
## d: disease classes. Used to stratify folds
## folds: number of folds
#
# Undergoes stratified folding of the data in folds sub sets
#
## return: list of lists of indices to denote which samples are in which fold
YieldFolds <- function(d, 
                       folds = 5) {
  
  remain.d <- d
  n <- length(d)
  fold.size <- floor(n/folds)
  
  fold.list <- list()
  for (fold in 1:folds) {
    fold.samples <- c()
    for (class in sort(unique(d))) {
      class.samples <- which(remain.d==class)
      
      class.fold.samples <- sample(class.samples, floor(table(d)[class]/folds), replace=F)
      fold.samples <- c(fold.samples, class.fold.samples)
      remain.d[fold.samples] <- NA
    }
    
    fold.list[[fold]] <- fold.samples
  }
  
  return(fold.list)
}


### ------------------------------------------------------------------------ ###
## classes: vector of class labels of each sample
## m: per class metrics to be averaged
#
# returns: vector of weighted averages by class 
WeightedAverage <- function(classes, m) {
  
  avgs <- m[1,]
  
  
  class.weights <- unname(table(classes)/length(classes))
  
  for (c in names(avgs)) {
    avgs[c] <- sum(class.weights * m[, c])
  }
  
  return(avgs)
}
