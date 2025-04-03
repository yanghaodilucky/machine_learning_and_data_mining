data.dir <- "/user/work/ms13525/tcga-bc-dataset"

## Function for loading tab-delimited spreadsheets
library(data.table)
read.dataset=function(filename, ...) {
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

## The format of sample identifiers/barcodes is described here:
## https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
##
## Function extracts the participant identifier from a sample id/barcode.
extract.participant <- function(id)
    sub("TCGA-[^-]+-([^-]+)-.*", "\\1", id)


## harmonize datasets
protein.ids=extract.participant(colnames(protein))
clinical.ids=rownames(clin)

common.ids=intersect(clinical.ids, protein.ids)
clin=clin[match(common.ids, clinical.ids),]
protein=protein[,match(common.ids, protein.ids)]

## restrict clinical dataset to clinical variables of interest
clinical.vars=c(
    "age.at.diagnosis","estrogen.receptor.status",
    "progesterone.receptor.status",
    "lymphocyte.infiltration","necrosis.percent")
target.var="pfi" ## outcome variable

clinical.dat=clinical.dat[,c(target.var,clinical.vars)]
clinical.dat$estrogen.receptor.status=ifelse(clinical.dat$estrogen.receptor.status=="positive",1,0)
clinical.dat$progesterone.receptor.status=ifelse(clinical.dat$progesterone.receptor.status=="positive",1,0)

## merge datasets
bc.dat=data.frame(clinical.dat, t(protein.dat))

## remove features with > 20% missing values
missing.pct=sapply(bc.dat, function(v) mean(is.na(v)))
bc.dat=bc.dat[,missing.pct < 0.2]

## replace any remaining missing values with the mean value
## (note: there are better ways to handle missing values!)
missing.pct=sapply(bc.dat, function(v) mean(is.na(v)))
if (any(missing.pct > 0))
    for (i in which(missing.pct > 0)) {
        missing.idx=which(is.na(bc.dat[[i]]))
        new.value=mean(bc.dat[[i]], na.rm=T)
        bc.dat[[i]][missing.idx]=new.value
    }

## train a basic elastic net model
library(mlr3verse)

bc.task = as_task_classif(
    x=bc.dat,
    target="pfi",
    id="Breast cancer PFI predicted by protein abundance")

set.seed(43)
dat.parts = partition(bc.task, 3/4)
## (note: there may be data leakage due to
## the fact that we replaced missing values
## using the mean across training and testing subsets
## above)

elnet = lrn(
    "classif.cv_glmnet",
    alpha=0.5,
    predict_type="prob")

elnet$train(bc.task, dat.parts$train)

## test the elastic net model
preds=elnet$predict(bc.task, dat.parts$train)
preds$confusion

