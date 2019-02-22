
plot.dendro <- function(ts.data, params, r.square = NULL,
  main.title = "", day = seq(365)) {
  params <- as.numeric(params)

  plot(ts.data$DOY, ts.data$DBH_TRUE, pch = 19, cex = 0.7, main = main.title,
    cex.main = 0.9, ylab = "DBH (cm)", xlab = "Day of year", xlim = c(0, 365),
    col = "gray")
  lines(day, lg5.pred(params = params,
     doy = day), col = "red", lty = 3)

  # lines(day, lg5.pred(params = params, day), col = "red", lty = 2)
  segments(0, params[6], pred.doy(params, params[6]), params[6], col = "black")
  segments(pred.doy(params, params[7]), params[7], 365, params[7], col = "black")

  if(!is.null(r.square)) mtext(sprintf("R-2: %.2f",
    r.square), 3, line = -2, adj = 0.2, cex = 0.7)
}



# plot.dendro.fit <- function(ts.data, p.means, params.opt,
#   l_p = NULL, params.ab = NULL, main.title = "", day = seq(365)) {

#   n.p.means <- dim(p.means)[2]
#   cols <- c(brewer.pal((n.p.means - 1), "Dark2"), "red", "black")

#   plot(ts.data$DOY, ts.data$DBH_TRUE, pch = 19, cex = 0.7, main = main.title,
#     cex.main = 0.9, ylab = "DBH (cm)", xlab = "Day of year", xlim = c(0, 365),
#     col = "gray")
#   for(i in 1:(n.p.means - 1)) {
#     lines(day, lg5.pred(params = params,
#      doy = day), col = cols[2], lty = 3)
#     lp <- p.means[7, i]
#   }
#   lines(day, lg5.pred(params = params, day), col = "red", lty = 2)
#   abline(h = ab)
#     col = cols,
#     lty = c(rep(3, n.p.means - 1), 1, 1), cex = 0.6, lwd = 2)
# }

#' Plotting Dendrometer Band Time Series
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TURE} (numeric) and \emph{DOY} (integer).
#' @param params dataframe vector
#' @param day integer vector
#' @param outlier logical
#' @param par.names  A character vector of the column names of the params vector. This is used to pull out actual parameters from other information.
#' @return A plot of the dbh and doy of a single band in a year, with, optionally, a fitted line from the optimize output and outlier denotion in red.
#' @seealso \code{\link{get.param}} which creates Dendro.split, a list vector containing time-series dataframes for every tree, year, band.
#' @export
plot.dendro.ts <- function(ts.data, params = NULL, day = seq(365), outlier = TRUE,
                           par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "alt.a")) {

  plot(ts.data$DOY, ts.data$DBH_TRUE, pch = 19, cex = 0.7,
       cex.main = 0.9, ylab = "DBH (cm)", xlab = "Day of year", xlim = c(0, 365),
       col = "gray",
       cex.main = 0.8,
       cex.axis = 0.8,
       cex.lab = 0.8,
       main = sprintf("SITE: %s | TREE_ID: %s | YEAR: %s",
                      ts.data$SITE[1],
                      ts.data$TREE_ID[1],
                      ts.data$YEAR[1]))

  if(!is.null(params)) {

    param.col <- names(params) %in% par.names[1:7]
    para <- as.numeric(params[param.col])

    lines(day, lg5.pred(params = para, doy = day), col = "red", lty = 3)

    segments(0, para[6], pred.doy(para, para[6]), para[6], col = "black")
    segments(pred.doy(para, para[7]), para[7], 365, para[7], col = "black")
    if(outlier == TRUE) {
    }
    flag.0 <- subset(ts.data, REMOVE == 1)
    points(flag.0$DOY, flag.0$DBH_TRUE, col = "red")
  }
}

#' Plotting Multiple Dendrometer Band Time Series
#'
#' @param Dendro.tree  A dataframe of a time series of a single tree over multiple years and bands.
#' Must have column variables \emph{DBH_TRUE} (numeric), \emph{DOY} (integer), and \emph{YEAR} (integer).
#' @param params dataframe vector
#' @param day integer vector
#' @param outlier logical
#' @param par.names  A character vector of the column names of the params vector. This is used to pull out actual parameters from other information.
#' @return A plot of the dbh and doy of a single band in a year, with, optionally, a fitted line from the optimize output and outlier denotion in red.
#'
#' @description Plots an extended time series of dendrometer band measurements (translated into DBH) with
#' fits, outliers, band position movements, slippage corrections, etc. Can be used for presentation
#' of single trees, or for diagnostics and fit assessments.
#' @export
plot.dendro.tree <- function(Dendro.tree, params = NULL,
    par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "alt.a"),
    par.dim = c(3, 1), print.pdf = FALSE, out.file = "FIGURES/TEST_PLOTS.pdf",
    day = seq(365)) {
  param.col <- names(params) %in% par.names[1:7]

  if(print.pdf == TRUE) pdf(file = out.file)
  par(mfrow = par.dim)

  for(i in 1:length(Dendro.data.ls)) {
    ts.data <- Dendro.data.ls[[i]]
    plot.dendro(ts.data,
      params = as.numeric(params[i, param.col]),
      r.square = as.numeric(params$r.square[i]),
      main.title = sprintf("SITE: %s | TREE_ID: %s | YEAR: %s",
        ts.data$SITE[1],
        ts.data$TREE_ID[1],
        ts.data$YEAR[1]))

    flag.0 <- subset(ts.data, Outlier == "O")
    points(flag.0$DOY, flag.0$DBH_TRUE, col = "red")
  }
  if(print.pdf == TRUE) dev.off()
}


find.outliers <- function(ts.data, sd.lim = 3) {
  ol.id <- which(abs(ts.data$sd.resids) > sd.lim)
  ts.data$REMOVE[ol.id] <- 1
  return(ts.data)
}

#' Identify outliers in dendrometer band time series
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer), as well as other designations.
#' @param params dataframe vector
#' @param day integer vector
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#'
#' @return none
#' @export
#'
#' @description An interactive figure for the identification and designation of slipped bands, outliers, and problematic bands.
#' @examples none
id.outliers <- function(ts.data, params, day = seq(365),
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "alt.a")) {

  param.col <- names(params) %in% par.names[1:7]

  for(i in 1:dim(params)[1]) {
    ts.data <- Dendro.data.ls[[i]]
    plot.dendro(ts.data,
      params = as.numeric(params[i, param.col]),
      r.square = NULL,
      main.title = sprintf("SITE: %s | TREE_ID: %s | YEAR: %s",
        ts.data$SITE[1],
        ts.data$TREE_ID[1],
        ts.data$YEAR[1]))
    flag.0 <- subset(ts.data, Outlier == 1)
    points(flag.0$DOY, flag.0$DBH_TRUE, col = "red")

    y.pos <- (min(ts.data$DBH_TRUE) +
      ((max(ts.data$DBH_TRUE) - min(ts.data$DBH_TRUE)) * 0.8))
    y.pos2 <- (min(ts.data$DBH_TRUE) +
      ((max(ts.data$DBH_TRUE) - min(ts.data$DBH_TRUE)) * 0.9))
    y.pos3 <- (min(ts.data$DBH_TRUE) +
      ((max(ts.data$DBH_TRUE) - min(ts.data$DBH_TRUE)) * 0.7))

    text(60, y.pos, labels = "Skip", col = "purple", cex = 1, pos = 2, offset = 1)
    points(60, y.pos, col = "purple", cex = 3, pch = 1)

    SK.id <- identify(60, y.pos, labels = "Will skip")

    if(length(SK.id) != 0) next

    text(60, y.pos2, labels = "Slip band", col = "orange", pos = 2, offset = 1)
    points(60, y.pos2, col = "orange", cex = 3, pch = 1)
    SL.id <- identify(60, y.pos2, labels = "ID slipped values")
    if(length(SL.id) != 0)
    SL.vals <- identify(ts.data$DOY, ts.data$DBH_TRUE, labels = "SL value")

    text(60, y.pos3, labels = "Now ID Outliers", col = "purple", cex = 1, pos = 2, offset = 1)
    OL.id <- identify(ts.data$DOY, ts.data$DBH_TRUE, labels = ts.data$DOY)

    Dendro.data.ls[[i]]$Outlier[OL.id] <- 1
    }
    return(do.call(rbind, Dendro.data.ls))
}

integrate.outlier.id <- function(Dendro.complete, new.outlier.df) {
  TREE.ID.YR <- paste(as.character(Dendro.complete$SITE), as.character(Dendro.complete$TREE_ID),
    as.character(Dendro.complete$YEAR), sep = "_")
  TREE.ID.YR.OL <- paste(as.character(new.outlier.df$SITE), as.character(new.outlier.df$TREE_ID),
    as.character(new.outlier.df$YEAR), sep = "_")
  new.OL.id <- match(TREE.ID.YR.OL, TREE.ID.YR)
}

# make.WD.plot <- function(ts.data, params, QHull.ls) {
#
# }
# #pdf("FIGURES/All_CH.pdf")
# layout(matrix(c(1,1,1,1,2,2,3,3,4,4), ncol = 2, byrow = TRUE))
# for(t in 1) { #:dim(param.table)[1]) {
#   params <- param.table[i, ]
#   dbh <- as.numeric(ts.data[i, ])
#   doy.full <- get.doy(ts.data)
#   doy <- doy.full[complete.cases(dbh)]
#   dbh <- dbh[complete.cases(dbh)]
#   OH.list <- fit.outer.hull(dbh, doy, params, quant = 0.8)

#   doyP2 <- OH.list$doyP2
#   dbhP2 <- OH.list$dbhP2
#   doyP <- OH.list$doyP
#   dbhP <- OH.list$dbhP
#   OH.fit <- OH.list$OH.fit
# # TODO: figure out start and stop values ...
#   start.d <- start.diam(params = as.numeric(params), seq.l = seq.l,  doy = doy, dbh, deviation.val = deviation.val, figure = FALSE, resid.sd)
#   end.d <- end.diam(params = as.numeric(params), seq.l, doy, dbh, deviation.val, figure = FALSE, resid.sd)
#   asym <- c(start.d[1], end.d[1])

#   plot(doy, dbh, xlab = "", ylab = "DBH (cm)", pch = 19, col = "gray15", main = sprintf("Annual Growth for tree %i", t), cex = 1)
#   points(doyP2, dbhP2, pch = 19, col = "tomato")
#   days <- seq(365)
#   lines(days, lg5.pred(params = OH.fit$par, doy = days), col = cols[1], lty = 1, lwd = 1)
#   #lines(days, lg5.pred.a(asym, params = OH.fit$par, doy = days, asymptote = "upper"), col = cols[1], lty = 1, lwd = 1)
#   lines(days, lg5.pred.a(asym, params = as.numeric(params), doy = days), col = cols[2], lty = 2, lwd = 1)
#   legend("bottomright", lty = c(1,2), col = cols[1:2], legend = c("Quantile Hull", "ML fit"))
#   text(110, 55.2, labels = "a)")

#   plot(doyP, OH.list$Deficit, pch = 19, type = "b", xlim = range(doy), col = "steelblue", xlab = "", ylab = "Deficit")
#   abline(h = 0)
#   text(110, min(OH.list$Deficit), "b)", pos = 3)

#   plot(doyP, OH.list$Weighted.deficit, pch = 19, type = "b", xlim = range(doy), col = "steelblue", xlab = "", ylab = "Weighted deficit")
#   abline(h = 0)
#   text(110, min(OH.list$Weighted.deficit), "c)", pos = 3)
#   ##
# }
# precip.data <- read.csv("WaterBalance.csv", header = TRUE)
# plot(precip.data$doy, precip.data$cum.NET.3.1, col = "orange", type = "l",
#      xlim = range(doy), pch = 19, xlab = "Day of the Year", ylab = "Water balance (mm)")
#   text(110, min(precip.data$cum.NET.3.1), "d)", pos = 3)



#############################
