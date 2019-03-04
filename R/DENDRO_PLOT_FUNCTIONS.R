###################################################
###  functions for Intensive dendrometer plotting
###################################################

#' Plotting Dendrometer Band Time Series
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer).
#' @param params dataframe vector (i.e., there are names associated with the vector elements)
#' @param day integer vector
#' @param outlier logical
#' @param par.names  A character vector of the column names of the params vector. This is used to pull out actual parameters from other information.
#' @return A plot of the dbh and doy of a single band in a year, with, optionally, a fitted line from the optimize output and outlier denotion in red.
#' @seealso \code{\link{get.params}}, which creates Dendro.split, a list vector containing time-series dataframes for every tree, year, band.
#' @export
make.dendro.plot.ts <- function(ts.data, params = NULL, day = seq(365), outlier = TRUE,
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

.sum.doy <- function(x, doy.diff) {
  return(rep(doy.diff, each = length(x)))
}


#' Plotting dendrometer Band Time Series for a single tree over years
#'
#' @param Dendro.ind  A dataframe of a time series of a single tree over multiple years and bands.
#' Must have column variables \emph{DBH_TRUE} (numeric), \emph{DOY} (integer), and \emph{YEAR} (integer).
#' @param param.tab dataframe of fit parameters
#'
#' @description Plots an extended time series of dendrometer band measurements (translated into DBH) with
#' fits, outliers, band position movements, slippage corrections, etc. Can be used for presentation
#' of single trees, or for diagnostics and fit assessments.
#' @export
make.dendro.plot.tree <- function(Dendro.ind, param.tab) {

  get.years      <- unique(Dendro.ind$YEAR) #Just the years
  start.date     <- paste(get.years[1], "-01-01", sep = "")
  end.date       <- paste(get.years[length(get.years)], "-12-31", sep = "")
  date.seq       <- format(seq(as.Date(start.date), as.Date(end.date),
                    by = "1 day")) # vector of dates y-m-d
  doy.seq        <- seq_along(date.seq) # sequence length of dates (for x axis)
  year.index     <- droplevels(cut(as.Date(date.seq), breaks = "year",
        right = FALSE, drop = TRUE))
  ln.year.ind    <- length(levels(year.index))
  year.split.seq <- split(doy.seq, year.index) # splits dates into years
  year.dates     <- split(date.seq, year.index) # splits dates into years
  year.doy.seq   <- lapply(year.dates, function(x) seq(length(x)))

  doy.shift <- cumsum(sapply(year.doy.seq, length)) - length(year.doy.seq[[1]])

  if (length(get.years) != ln.year.ind) return(NULL)
  # SETTING UP PLOTTING

  doy.ls <- split(Dendro.ind$DOY, Dendro.ind$YEAR)
  doy.new.tmp <- vector("list", length(doy.ls))
  for(i in 1:length(doy.ls)) {
    doy.new.tmp[[i]] <- doy.ls[[i]] + .sum.doy(doy.ls[[i]], doy.shift[i])
  }
  doy.4.plotting <- as.integer(unlist(doy.new.tmp))

  # GETTING AXES AND MARKERS TOGETHER
  axis.date.at <- as.integer(which(unlist(year.doy.seq) == 180))
  year.at <- as.integer(which(unlist(year.doy.seq) == 1))
  axis.dates <- date.seq[axis.date.at]

  # Need to put any par() adjustments or meta-plot info
  cols <- ifelse(Dendro.ind$ADJUST != 0, "gold", "black")
  # Plot the data and establish the axes
  plot(doy.4.plotting, Dendro.ind$DBH_TRUE, col = cols, pch = 18, cex = 0.5,
    main = sprintf(" %s | %s | %s", Dendro.ind$SITE[1], Dendro.ind$TREE_ID[1],
      Dendro.ind$SP[1]),
    ylab = "DBH (cm)", xlab = "Date", axes = FALSE)
  axis(2)
  axis(1, labels = get.years, at = axis.date.at, tick = FALSE)
  box()
  abline(v = year.at, col = "gray", lty = 3)

  ##################################################
  # Plot functional fits to the data
  # with a and b from fits
  cols <- c("steelblue", "tomato")
  p.years <- param.tab$YEAR
  fit.ls <- vector("list", length(p.years))
  pred.dbh <- vector("list", dim(param.tab)[1])
  pred.dbh2 <- vector("list", dim(param.tab)[1])
  ##################################################
  # Plot new band data

  ##################################################
  # Plot extensions of last year's b to next year's a
  #DONE AS PART OF alt.
  par.base <- c("L", "K", "doyip", "r", "theta")
  par.ab <- c("a", "b")
  param.base.col <- names(param.tab) %in% par.base
  param.ab.col <- names(param.tab) %in% par.ab
  param.cc <- complete.cases(param.tab[, param.base.col])

  for(y in 1:dim(param.tab)[1]) {
    if(param.cc[y] == FALSE) next
    param.y <- as.numeric(param.tab[y, param.base.col])
    param.ab <- as.numeric(param.tab[y, param.ab.col])
    param.ab.alt <- c(param.tab$alt.a[y], param.ab[2]) #as.numeric(param.tab[y, c(11:12)])
    pred.dbh[[y]] <- lg5.pred.a(a = param.ab, params = param.y,
      doy = year.doy.seq[[y]])
    pred.dbh2[[y]] <- lg5.pred.a(a = param.ab.alt, params = param.y,
      doy = year.doy.seq[[y]])
    lines(year.split.seq[[y]] , pred.dbh2[[y]], col = cols[2], lty = 2)
    lines(year.split.seq[[y]] , pred.dbh[[y]], col = cols[1])
  }

  # param.tab.final <- param.tab

  # if (is.b.band) {
  #   b.band.yrs <- b.band.param.tab$Year
  #   slot.year <- match(b.band.yrs, get.years)
  #   for(y in 1:length(slot.year)) {
  #     param.y <- as.numeric(b.band.param.tab[y, c(6:10)])
  #     param.ab <- as.numeric(b.band.param.tab[y, 11:12]) #as.numeric(param.tab[y, c(11:12)])
  #     pred.dbh[[y]] <- lg5.pred.a(a = param.ab, params = param.y,
  #       doy = year.doy.seq[[slot.year[y]]])
  #     lines(year.split.seq[[slot.year[y]]] , pred.dbh[[y]], col = cols[1])
  #   }
  #   param.y <- as.numeric(b.band.param.tab[y, c(6:10)])
  #   param.ab <- param.y[1:2] #as.numeric(param.tab[y, c(11:12)])
  #   pred.dbh[[y]] <- lg5.pred.a(a = param.ab, params = param.y,
  #     doy = year.doy.seq[[y]])
  #   lines(year.split.seq[[y]] , pred.dbh[[y]], col = cols[3])
  #   b.band.param.tab$alt.a <- param.tab$alt.a[slot.year]
  #   param.tab.final <- rbind(param.tab, b.band.param.tab)
  # }
}


#' Automatically identify outliers in dendrometer band time series
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer),
#' as well as other designations.
#' @param params dataframe vector
#' @param day integer vector
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#'
#'
#' @description Identifies outliers that are beyond the determined standard
#' deviation of the residuals. Assumes residuals are from a fit model.
#' @export
find.outliers <- function(ts.data, sd.lim = 3) {
  ol.id <- which(abs(ts.data$sd.resids) > sd.lim)
  ts.data$REMOVE[ol.id] <- 1
  return(ts.data)
}

#' Identify outliers in dendrometer band time series
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer),
#' as well as other designations.
#' @param params dataframe vector
#' @param day integer vector, defaults to seq(365)
#' @param par.names  A character vector of the column names of the \emph{params}
#' vector argument.
#' This is used to pull out actual parameters from other information.
#'
#' @export
#'
#' @description An interactive figure for the identification and designation of
#' slipped bands, outliers, and problematic bands.
id.outliers <- function(ts.data, params, day = seq(365),
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "alt.a")) {

  param.col <- names(params) %in% par.names[1:7]

  for(i in 1:dim(params)[1]) {
    ts.data <- Dendro.data.ls[[i]]
    make.dendro.plot(ts.data,
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

.integrate.outlier.id <- function(Dendro.complete, new.outlier.df) {
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
