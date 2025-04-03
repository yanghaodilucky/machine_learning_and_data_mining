# student number: 2398736 
# date: 2024.5.15 
#
#
# Part 1: Load Packages
#
# Part 2: Load Data and Processing
#    - Load Data
#    - Organize the Participants
#    - Restrict clinical data set to clinical variables of interest 
#    - Merge the data sets 
#    - Missing Data 
#    - Normalization 
#
# Part 3: Machine Learning Models
#    - XgBoost with clinic data 
#    - XgBoost with clinic and protein data 
#    - XgBoost with clinic and miRNA data 
#    - XgBoost with clinic and mRNA data 
#    - XgBoost with clinic and mutation data 
#    - XgBoost with clinic and all above omic data 
#    - Elastic Net with clinic and mutation data 
#    - Decision Tree with clinic and mutation data 
#    - Random Forest with clinic and mutation data 
#    - LASSO with clinic and mutation data 
#
#
# Part 4: Results and Visualization
#    - Table comparing different data sets
#    - Plot comparing different data sets 
#    - Table comparing different models 
#    - Plot comparing different models 



# ==============================================================================
# Part 1: Load Packages
# ==============================================================================

rm(list=ls()) # clear the environment

library(data.table) # load data
library(dplyr) 
library(VIM) # for missing data
library(mlr3)
library(mlr3pipelines)
library(mlr3verse)
library(Matrix)
library(glmnet) # for elastic net and LASSO
library(precrec)
library(ggplot2)
library(xgboost)
library(rpart) # for tree
library(gridExtra) # for combining the plots

# ==============================================================================
# Part 2: Load the Data and processing
# ==============================================================================


# ---- Load Data ----
data.dir <- "data/"

read.dataset=function(filename, ...) {       # define the read.dataset function  
  cat("reading", basename(filename), "... ")
  ## read in tab-delimited spreadsheet
  x=fread(
    filename,
    header=T,
    stringsAsFactors=F,
    sep="\t",
    check.names=F,
    ...)
  ## remove any duplicate rows (identified by the first column)
  x=x[match(unique(x[[1]]), x[[1]]),]
  ## make the first column the rownames of the data frame
  x=data.frame(x,row.names=1,stringsAsFactors=F,check.names=F)
  cat(nrow(x), "x", ncol(x), "\n")
  x
}

filename <- "clinical.txt"
full_path <- file.path(data.dir, filename)
clin <- read.dataset(full_path)


filename <- "mrna.txt"
full_path <- file.path(data.dir, filename)
mrna <- read.dataset(full_path)

filename <- "mirna.txt"
full_path <- file.path(data.dir, filename)
mirna <- read.dataset(full_path)

filename <- "protein.txt"
full_path <- file.path(data.dir, filename)
protein <- read.dataset(full_path)

filename <- "mutations.txt"
full_path <- file.path(data.dir, filename)
mutations <- read.dataset(full_path)


# ---- organize the participants ----

# define extract id function
extract.participant <- function(id)
  sub("TCGA-[^-]+-([^-]+)-.*", "\\1", id) 

# find the ids in every dataset
protein.ids=extract.participant(colnames(protein))
mrna.ids=extract.participant(colnames(mrna))
mirna.ids=colnames(mirna)
mutations.ids=colnames(mutations)
clinical.ids=rownames(clin)

# find these patients who have each kind of record
common.ids=intersect(clinical.ids, protein.ids)
common.ids=intersect(common.ids, mrna.ids)
common.ids=intersect(common.ids, mirna.ids)
common.ids=intersect(common.ids, mutations.ids)

# reserve those with all of the records
clin=clin[match(common.ids, clinical.ids),]
protein=protein[,match(common.ids, protein.ids)]
mrna=mrna[,match(common.ids, mrna.ids)]
mirna=mirna[,match(common.ids, mirna.ids)]
mutations=mutations[,match(common.ids, mutations.ids)]


# ---- restrict clinical dataset to clinical variables of interest ----
clinical.vars=c(
  "age.at.diagnosis","estrogen.receptor.status", "stage",
  "progesterone.receptor.status","her2.status",
  "lymphocyte.infiltration","necrosis.percent")
target.var="pfi" ## outcome variable

clin=clin[,c(target.var,clinical.vars)]

clin$estrogen.receptor.status=ifelse(
  clin$estrogen.receptor.status=="positive",1,0)

clin$her2.status=ifelse(
  clin$her2.status=="positive",1,0)

clin$progesterone.receptor.status=ifelse(
  clin$progesterone.receptor.status=="positive",1,0)

change_stage <- function(stage) { 
  case_when(
    stage == "stage i"    ~ 1,
    stage == "stage ia"   ~ 2,
    stage == "stage ib"   ~ 3,
    stage == "stage ii"   ~ 4,
    stage == "stage iia"  ~ 5,
    stage == "stage iib"  ~ 6,
    stage == "stage iii"  ~ 7,
    stage == "stage iiia" ~ 8,
    stage == "stage iiib" ~ 9,
    stage == "stage iiic" ~ 10,
    stage == "stage iv"   ~ 11,
    TRUE                  ~ NA_integer_  # 处理 NA 和其他未指定的情况
  )
}
clin$stage <- change_stage(clin$stage)

clin$pfi <- as.logical(clin$pfi)


# ---- merge the datasets ----


# in the training section, we will try different combination of the datasets

data1_clin <- clin

data2_clin_mrna = data.frame(clin, t(mrna)) 

data3_clin_prot = data.frame(clin, t(protein)) 

data4_clin_mut = data.frame(clin, t(mutations))

data5_clin_mirna = data.frame(clin, t(mirna))

data_mrna_mirna = data.frame(t(mrna), t(mirna))

data_mrna_mirna_mut = data.frame(data_mrna_mirna, t(mutations))

data6_all = data.frame(data3_clin_prot, data_mrna_mirna_mut)


# ---- Missing Data ----

# and we have two steps to deal with it: 

### first remove features with > 20% missing values
missing.pct=sapply(data1_clin, function(v) mean(is.na(v)))
data1_clin=data1_clin[,missing.pct < 0.2]

missing.pct=sapply(data2_clin_mrna, function(v) mean(is.na(v)))
data2_clin_mrna=data2_clin_mrna[,missing.pct < 0.2]

missing.pct=sapply(data3_clin_prot, function(v) mean(is.na(v)))
data3_clin_prot=data3_clin_prot[,missing.pct < 0.2]

missing.pct=sapply(data4_clin_mut, function(v) mean(is.na(v)))
data4_clin_mut=data4_clin_mut[,missing.pct < 0.2]

missing.pct=sapply(data5_clin_mirna, function(v) mean(is.na(v)))
data5_clin_mirna=data5_clin_mirna[,missing.pct < 0.2]

missing.pct=sapply(data6_all, function(v) mean(is.na(v)))
data6_all=data6_all[,missing.pct < 0.2]

### second replace any remaining missing values with the mean value
missing.pct=sapply(data1_clin, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data1_clin[[i]]))
    new.value=mean(data1_clin[[i]], na.rm=T)
    data1_clin[[i]][missing.idx]=new.value
  }

missing.pct=sapply(data2_clin_mrna, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data2_clin_mrna[[i]]))
    new.value=mean(data2_clin_mrna[[i]], na.rm=T)
    data2_clin_mrna[[i]][missing.idx]=new.value
  }

missing.pct=sapply(data3_clin_prot, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data3_clin_prot[[i]]))
    new.value=mean(data3_clin_prot[[i]], na.rm=T)
    data3_clin_prot[[i]][missing.idx]=new.value
  }

missing.pct=sapply(data4_clin_mut, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data4_clin_mut[[i]]))
    new.value=mean(data4_clin_mut[[i]], na.rm=T)
    data4_clin_mut[[i]][missing.idx]=new.value
  }

missing.pct=sapply(data5_clin_mirna, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data5_clin_mirna[[i]]))
    new.value=mean(data5_clin_mirna[[i]], na.rm=T)
    data5_clin_mirna[[i]][missing.idx]=new.value
  }


missing.pct=sapply(data6_all, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
  for (i in which(missing.pct > 0)) {
    missing.idx=which(is.na(data6_all[[i]]))
    new.value=mean(data6_all[[i]], na.rm=T)
    data6_all[[i]][missing.idx]=new.value
  }

# ---- normalization ----

data1_clin <- data1_clin %>%
  mutate_if(is.numeric, scale)

data2_clin_mrna <- data2_clin_mrna %>%
  mutate_if(is.numeric, scale)

data3_clin_prot <- data3_clin_prot %>%
  mutate_if(is.numeric, scale)

data4_clin_mut <- data4_clin_mut %>%
  mutate_if(is.numeric, scale)

data4_clin_mut <- data4_clin_mut[, colSums(is.na(data4_clin_mut)) == 0]

data5_clin_mirna <- data5_clin_mirna %>%
  mutate_if(is.numeric, scale)

data5_clin_mirna <- data5_clin_mirna[, colSums(is.na(data5_clin_mirna)) == 0]

data6_all <- data6_all %>%
  mutate_if(is.numeric, scale)




# ==============================================================================
# Part 3: Models
# ==============================================================================

# ---- XgBoost with clinic data ----
d1x.task = as_task_classif(
  x = data1_clin,
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data"
)

set.seed(43)

dat.parts = partition(d1x.task, ratio = 4/5)

d1x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 6,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d1x$train(d1x.task, dat.parts$train)

# test the model
preds <- d1x$predict(d1x.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd1x = resample(d1x.task, d1x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd1x$aggregate(metrics)

roc_d1x <- autoplot(rrd1x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d1x <- roc_d1x + ggtitle("Xgboost ROC Curve with Clinic Data")
auc_d1x <- autoplot(rrd1x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic data.png", 
       plot = roc_d1x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)


# ---- XgBoost with clinic and protein data ----

d3x.task = as_task_classif(
  x = data3_clin_prot,
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and protein abundance"
)

set.seed(43)

dat.parts = partition(d3x.task, ratio = 4/5)

d3x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 6,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d3x$train(d3x.task, dat.parts$train)

# test the model
preds <- d3x$predict(d3x.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd3x = resample(d3x.task, d3x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd3x$aggregate(metrics)

roc_d3x <- autoplot(rrd3x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d3x <- roc_d3x + ggtitle("Xgboost ROC Curve with Clinic and Protein Data")
auc_d3x <- autoplot(rrd3x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic and protein data.png", 
       plot = roc_d3x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)


# ---- XgBoost with clinic and miRNA data ----

d5x.task = as_task_classif(
  x = data5_clin_mirna,
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and mirna abundance"
)

set.seed(43)

dat.parts = partition(d5x.task, ratio = 4/5)

d5x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 6,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d5x$train(d5x.task, dat.parts$train)

# test the model
preds <- d5x$predict(d5x.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd5x = resample(d5x.task, d5x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd5x$aggregate(metrics)

roc_d5x <- autoplot(rrd5x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d5x <- roc_d5x + ggtitle("Xgboost ROC Curve with Clinic and miRNA Data")
auc_d5x <- autoplot(rrd5x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic and mirna data.png", 
       plot = roc_d5x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)

# ---- XgBoost with clinic and mRNA data ----

d2x.task = as_task_classif(
  x = data2_clin_mrna,
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and mrna data"
)

set.seed(43)

dat.parts = partition(d2x.task, ratio = 4/5)

d2x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 6,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d2x$train(d2x.task, dat.parts$train)

# test the model
preds <- d2x$predict(d2x.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd2x = resample(d2x.task, d2x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd2x$aggregate(metrics)

roc_d2x <- autoplot(rrd2x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d2x <- roc_d2x + ggtitle("Xgboost ROC Curve with Clinic and mRNA Data")
auc_d2x <- autoplot(rrd2x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic and mrna data.png", 
       plot = roc_d2x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)

# ---- XgBoost with clinic and mutation data ----

d4.task = as_task_classif(
  x = data4_clin_mut, 
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and mutation data"
)

set.seed(43)

dat.parts = partition(d4.task, ratio = 3/4)

d4x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 15,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d4x$train(d4.task, dat.parts$train)

# test the model
preds <- d4x$predict(d4.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd4x = resample(d4.task, d4x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd4x$aggregate(metrics)

roc_d4x <- autoplot(rrd4x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d4x <- roc_d4x + ggtitle("Xgboost ROC Curve with Clinic and Mutation Data")
auc_d4x <- autoplot(rrd4x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic and mutation data.png", 
       plot = roc_d4x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)

# ---- XgBoost with clinic and all above omic data ----

d6x.task = as_task_classif(
  x = data6_all, 
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and all omic data"
)

set.seed(43)

dat.parts = partition(d6x.task, ratio = 4/5)

d6x = lrn(
  "classif.xgboost",
  booster = "gbtree",
  eta = 0.1,
  max_depth = 6,
  nrounds = 100,
  predict_type = "prob"
)

# train the model
d6x$train(d6x.task, dat.parts$train)

# test the model 
preds <- d6x$predict(d6x.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd6x = resample(d6x.task, d6x, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd6x$aggregate(metrics)

roc_d6x <- autoplot(rrd6x, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d6x <- roc_d6x + ggtitle("Xgboost ROC Curve with Clinic and all omic Data")
auc_d6x <- autoplot(rrd4x, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of xgboost with clinic and all of the omic data.png", 
       plot = roc_d6x, 
       width = 10, height = 8, 
       units = "in", dpi = 300)


# ---- Elastic Net with clinic and mutation data ----

d4.task = as_task_classif(
  x = data4_clin_mut, 
  target = "pfi",
  id = "Breast cancer PFI predicted by clinic data and mutation data"
)

set.seed(43)
dat.parts = partition(d4.task, 3/4)


d4e = lrn(
  "classif.cv_glmnet",
  alpha=0.5,
  predict_type="prob")
#train the model
d4e$train(d4.task, dat.parts$train)


preds <- d4e$predict(d4.task, dat.parts$test)  # test the model in test data
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd4e = resample(d4.task, d4e, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd4e$aggregate(metrics)

roc_d4e <- autoplot(rrd4e, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d4e <- roc_d4e + ggtitle("Elastic Net ROC Curve with Clinic and Mutation Data")
auc_d4e <- autoplot(rrd4e, measure=msr("classif.auc"))

ggsave(filename = 
         "visulization/result of elastic net with clinic data and mutation data.png", 
       plot = roc_d4e, 
       width = 10, height = 8, 
       units = "in", dpi = 300)






# ---- Decision Tree with clinic and mutation data ----

set.seed(43)

dat.parts = partition(d4.task, ratio = 3/4)

d4dt = lrn(
  "classif.rpart",
  minsplit=20,
  predict_type = "prob"
)

# train the model
d4dt$train(d4.task, dat.parts$train)

# test the model
preds <- d4dt$predict(d4.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd4dt = resample(d4.task, d4dt, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd4dt$aggregate(metrics)

roc_d4dt <- autoplot(rrd4dt, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d4dt <- roc_d4dt + ggtitle("Decision Tree ROC Curve with Clinic and Mutation Data")

ggsave(filename = 
         "visulization/result of decision tree with clinic and mutation data.png", 
       plot = roc_d4dt, 
       width = 10, height = 8, 
       units = "in", dpi = 300)


# ---- Random Forest with clinic and mutation data ----

n = d4.task$nrow
p = d4.task$ncol
set.seed(43)

dat.parts = partition(d4.task, ratio = 4/5)

d4rf = lrn(
  "classif.ranger",
  oob.error=TRUE,
  importance="impurity",
  num.trees=100,
  mtry=floor(p/10),
  max.depth=floor(log2(p/16)),
  min.node.size=16,
  predict_type="prob"
)

# train the model
d4rf$train(d4.task, dat.parts$train)

# test the model
preds <- d4rf$predict(d4.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd4rf = resample(d4.task, d4rf, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd4rf$aggregate(metrics)

roc_d4rf <- autoplot(rrd4rf, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d4rf <- roc_d4rf + ggtitle("Random Forest ROC Curve with Clinic and Mutation Data")

ggsave(filename = 
         "visulization/result of random forest with clinic and mutation data.png", 
       plot = roc_d4rf, 
       width = 10, height = 8, 
       units = "in", dpi = 300)

# ---- LASSO with clinic and mutation data ----

set.seed(43)

dat.parts = partition(d4.task, ratio = 4/5)

d4ls = lrn(
    "classif.cv_glmnet",
    alpha=1,
    predict_type="prob"
)

# train the model
d4ls$train(d4.task, dat.parts$train)

# test the model
preds <- d4ls$predict(d4.task, dat.parts$test)
confusion_matrix <- preds$confusion
print(confusion_matrix)

# evaluate the performance of this learner using cross-validation

rrd4ls = resample(d4.task, d4ls, rsmp("cv",folds=3))

metrics = c("auc","acc","sensitivity","specificity",
            "precision","recall","fbeta")
metrics = paste("classif", metrics, sep=".")

metrics = lapply(metrics, msr)

rrd4ls$aggregate(metrics)

roc_d4ls <- autoplot(rrd4ls, type = "roc") +
  theme_minimal() +  
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),  
    plot.background = element_rect(fill = "white", colour = "white")  
  )
roc_d4ls <- roc_d4ls + ggtitle("LASSO ROC Curve with Clinic and Mutation Data")

ggsave(filename = 
         "visulization/result of LASSO with clinic and mutation data.png", 
       plot = roc_d4ls, 
       width = 10, height = 8, 
       units = "in", dpi = 300)


# ==============================================================================
# Part 4: output the result
# ==============================================================================

# ---- table comparing clinic with or without omic data ----

list_of_rrx <- list(rrd1x, rrd2x, rrd3x, rrd4x, rrd5x, rrd6x)
res_x <- data.frame()

for(rrd in list_of_rrx) {  # combine the result together
  
  aggregated_metrics <- rrd$aggregate(metrics)
 
  res_x <- rbind(res_x, aggregated_metrics)
}

print(res_x)

# output the csv in result folder

colnames(res_x) <- c(
  "auc", "acc", "sensitivity", "specificity", "precision", "recall", "fbeta")

rownames(res_x) <- c(
  "clinic", "clinic & mrna", 
  "clinic & protein", "clinic & mutation",
  "clinic & mirna", "clinic & mrna,mirna,protein, mutation")

write.csv(res_x, "result/table1.csv", row.names = TRUE)

# ---- plot comparing clinic with or without omic data ----

combined_plot <- grid.arrange(
  roc_d1x, roc_d2x, roc_d3x, roc_d4x, roc_d5x, roc_d6x, 
  ncol = 2)

ggsave("visulization/roc_clinic_omic.png", 
       plot = combined_plot, width = 10, height = 16, units = "in", dpi = 300)


# ---- table comparing different models ----

list_of_rr4 <- list(rrd4x, rrd4e, rrd4dt, rrd4rf, rrd4ls)
res_x2 <- data.frame()

for(rrd in list_of_rr4) {  # combine the result together
  
  aggregated_metrics <- rrd$aggregate(metrics)
  
  res_x2 <- rbind(res_x2, aggregated_metrics)
}

print(res_x2)

# output the csv in result folder

colnames(res_x2) <- c(
  "auc", "acc", "sensitivity", "specificity", "precision", "recall", "fbeta")

rownames(res_x2) <- c(
  "xgboost", "elasric net", "decision tree", "random forest", "LASSO")

write.csv(res_x2, "result/table2.csv", row.names = TRUE)

# ---- Plot comparing different models ----

combined_plot2 <- grid.arrange(
  roc_d4x, roc_d4e, roc_d4dt, roc_d4rf, roc_d4ls, 
  ncol = 2)

ggsave("visulization/roc_different_models.png", 
       plot = combined_plot2, width = 10, height = 16, units = "in", dpi = 300)


