## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide"---------------------
devtools::install_github("seanmcm/RDendrom")

## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide"---------------------
library(RDendrom)
data(INPUT_data)
head(INPUT.data)

## ----eval=FALSE, echo=TRUE, tidy=TRUE, results="hide"--------------------
## band.index <- as.numeric(table(INPUT.data$BAND_NUM))
## INPUT.data$BAND_NUM <- unlist(mapply(rep, seq(length(band.index)), length.out = band.index))
## 
## INPUT.data <- subset(INPUT.data, complete.cases(INPUT.data$GAP_WIDTH))

## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide"---------------------
get.optimized.dendro(INPUT.data,
  no.neg.growth = TRUE, cutoff = 9,
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "r.squared", "ts.sd"),
  OUTPUT.folder = ".",
  param.table.name = sprintf("Param_table.csv"),
  Dendro.data.name = sprintf("Dendro_data.Rdata"),
  Dendro.split.name = sprintf("Dendro_data_split.Rdata"))


## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide"---------------------
param.table.name = sprintf("Param_table.csv")
Dendro.data.name = sprintf("Dendro_data.Rdata")
Dendro.split.name = sprintf("Dendro_data_split.Rdata")
OUTPUT.folder <- c(".")
param.table <- read.csv(file = paste(OUTPUT.folder, param.table.name, sep = "/"))
load(file = paste(OUTPUT.folder, Dendro.data.name, sep = "/")) # loads Dendro.complete
load(file = paste(OUTPUT.folder, Dendro.split.name, sep = "/")) #loads Dendro.split
load(file = paste(OUTPUT.folder, "Dendro_Tree.Rdata", sep = "/")) #loads Dendro.tree

get.extra.metrics(
  param.table,
  Dendro.split,
  OUTPUT.folder      = ".",
  param.table.name = sprintf("Param_table_complete.csv"),
  Dendro.data.name = sprintf("Dendro_data_complete.Rdata"),
  Dendro.split.name = sprintf("Dendro_data_split_complete.Rdata"),
  Quantile.hull.name = "Quantile.hull.Rdata")


## ----eval=TRUE, echo=3:4, tidy=TRUE, results="markup"--------------------
param.table.name = sprintf("Param_table_complete.csv")
OUTPUT.folder <- c(".")
param.table.extended <- read.csv(file = paste(OUTPUT.folder, param.table.name, sep = "/"))
head(param.table.extended)


## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide", fig.width=5, fig.height=3, fig.cap="Individual time series, such as bands in years for individuals."----
make.dendro.plot.ts(ts.data = Dendro.split[[3]], params = param.table[3, ],
  day = seq(365))


## ----eval=TRUE, echo=TRUE, tidy=TRUE, results="hide", fig.width=7, fig.height=4, fig.cap="Single tree with multiple bands over years."----
make.dendro.plot.tree(Dendro.ind = Dendro.tree[[1]], param.tab = subset(param.table, TREE_ID == Dendro.tree[[1]]$TREE_ID[1]))



