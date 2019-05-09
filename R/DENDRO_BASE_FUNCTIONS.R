###################################################
###  functions for Intensive dendrometer workflow
###################################################

#############################
#############################

#' Fit the LG5 function
#'
#' @param INPUT.data A data.frame of a specific form (see @details)
#' @param no.neg.growth logical denoting whether to skip any time series with negative trends (from a call to \emph{lm()})
#' @param cutoff Integer denoting the minimum sample size required to perform the optimization call.
#' @param units Character string denoting the units that the dbh is in. Defaults to "cm".
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#' @param OUTPUT.folder Character string for the target output folder.
#' @param param.table.name Character string for the parameter table output.
#' @param Dendro.data.name Character string for the data object, Dendro.complete, that is a complete data.frame table output.
#' @param Dendro.split.name Character string for the data object, Dendro.split, that is a list vector where
#' every entry is a separate time series of a band on a stem in a year.
#'
#' @description This is the wrapper for \emph{get.params()} which uses base \emph{optim()} (with two sequential methods) to
#' estimate the parameters of the logistic function. This follows McMahon and Parker 2014. This function organizes output, skips removed datasets, fits a linear model to detect negative or non-significant growth, eliminates time series with small sample sizes, and removes identified outliers (see @identify.outliers ).
#' @return Objects Saves four files to the OUTPUT folder as named above and collected stems.
#' @export
get.optimized.dendro <- function(INPUT.data,
  no.neg.growth = TRUE, cutoff = 9, units = "cm",
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "r.squared", "ts.sd"),
  OUTPUT.folder = ".",
  param.table.name = "Param_table.csv",
  Dendro.data.name = "Dendro_data.Rdata",
  Dendro.split.name = "Dendro_data_split.Rdata") {

  options(warn = -1)
  TREE.ID.YR <- paste(as.character(INPUT.data$SITE), as.character(INPUT.data$TREE_ID),
    as.character(INPUT.data$YEAR), sep = "_")
  Dendro.split <- vector("list", length = length(unique(TREE.ID.YR)))
  ind.dendro <- split(INPUT.data, f = INPUT.data$TREE_ID)

  n.obs <- length(ind.dendro)

  param.table <-c()
  pb <- txtProgressBar(style = 3)
  Dendro.tree <- vector("list", length(ind.dendro))
  ct <- 1
  for(i in 1:n.obs) {
    setTxtProgressBar(pb, i / n.obs, title = NULL, label = NULL)
    ind.data <- ind.dendro[[i]] # loads an individual (multiple years)
    band.index <- as.numeric(table(ind.data$BAND_NUM))
    ind.data$BAND_NUM <- unlist(mapply(rep, seq(length(band.index)), length.out = band.index))

    ind.data$DBH_TRUE[1] <- ind.data$ORG_DBH[1]
    for(v in 2:length(ind.data$DBH)) {
      ind.data$DBH_TRUE[v] <- gettruedbh(gw1 = ind.data$GAP_WIDTH[v - 1],
        gw2 = ind.data$GAP_WIDTH[v], dbh1 = ind.data$DBH_TRUE[v - 1], units = units)
    }

    if(max(ind.data$BAND_NUM) > 1) {

      # FIX NEW BAND ISSUE
      nb.index <- which(!duplicated(ind.data$BAND_NUM))[-1]

      ind.data.band <- split(ind.data, ind.data$BAND_NUM)

      # NOW GET ALL OF THE CORRECT "DBHs" for all bands
      for(b in 2:length(ind.data.band)) {
        year.nb <- ind.data.band[[b]]$YEAR[1]
        year.ob <- ind.data.band[[b - 1]]$YEAR
        old.band.data <- subset(ind.data.band[[b - 1]], YEAR == year.nb)
        if (dim(old.band.data)[1] != 0) {
         doy.band <- max(which(old.band.data$DOY <= ind.data.band[[b]]$DOY[1]))
         ind.data.band[[b]]$DBH_TRUE[1] <- old.band.data$DBH_TRUE[doy.band]
       } else {
         ind.data.band[[b]]$DBH_TRUE[1] <-
         ind.data.band[[b - 1]]$DBH_TRUE[nrow(ind.data.band[[b - 1]])]
       }

       if(nrow(ind.data.band[[b]]) < 2) next #for case with a band with one measurement

       for(v in 2:nrow(ind.data.band[[b]])){
        ind.data.band[[b]]$DBH_TRUE[v] <- gettruedbh(gw1 = ind.data.band[[b]]$GAP_WIDTH[v - 1],
          gw2 =  ind.data.band[[b]]$GAP_WIDTH[v], dbh1 = ind.data.band[[b]]$DBH_TRUE[v - 1],
          units = units)
      }

    }
    ind.data <- unsplit(ind.data.band, ind.data$BAND_NUM)
  }

  Dendro.tree[[i]] <- ind.data
  ind.year.band <- split(ind.data,
    f = list(YEAR = ind.data$YEAR, BAND_NUM = ind.data$BAND_NUM), drop = TRUE)

  params <- rep(NA, 7)
  r.squared <- -99
  param.mat <- matrix(NA, length(ind.year.band), length(par.names))
  years <- vector("integer", length(ind.year.band))
  band.no <- vector("integer", length(ind.year.band))

  for(t in 1:length(ind.year.band)) {
    ts.data.tmp <- ind.year.band[[t]]
    if(any(ts.data.tmp$ADJUST == 1)) {
      ts.data.tmp <- .make.adjust(ts.data.tmp)
    }
    ts.data <- subset(ts.data.tmp, REMOVE == 0)
    ts.sd <- sd(ts.data$DBH_TRUE, na.rm = TRUE)
    years[t] <- ts.data$YEAR[1]
    band.no[t] <- ts.data$BAND_NUM[1]

    if (sum(!is.na(ts.data$DBH_TRUE)) < cutoff | any(ts.data$SKIP == 1)) {
      param.mat[t, ] <- c(params, r.squared, ts.sd)
      ct <- ct + 1
      next
    }

    if (min(ts.data$DOY) > 140 | max(ts.data$DOY) < 250) {
      param.mat[t, ] <- c(params, r.squared, ts.sd)
      ct <- ct + 1
      next
    }

    if (no.neg.growth == TRUE) {
      lm.out <- coef(lm(ts.data$DBH_TRUE ~ ts.data$DOY ))
      if(as.numeric(lm.out[2]) < 0) {
        ts.data$SKIP[1] <- 1
        param.mat[t, ] <- c(params, r.squared, ts.sd)
        ct <- ct + 1
        next
      }
    }
    params <- get.params(ts.data)
    r.squared <- summary(lm(ts.data$DBH_TRUE ~ lg5.pred.a(params[6:7],
      params, doy = ts.data$DOY)))$r.squared
    param.mat[t, ] <- c(params, r.squared, ts.sd)
    Dendro.split[[ct]] <- ts.data
    ct <- ct + 1
  }

  param.tab.tmp <- data.frame(SITE = ts.data$SITE[1],
    YEAR = years, TREE_ID = ts.data$TREE_ID[1],
    BAND_NUM = band.no, UNIQUE_ID = ts.data$UNIQUE_ID[1],
    SP = ts.data$SP[1], param.mat)
  param.table <- rbind(param.table, param.tab.tmp)

}
close(pb)
Dendro.split <- Dendro.split[1:(ct - 1)]
par.ind <- grep("X", names(param.table))
names(param.table)[par.ind] <- par.names ## MAKE SURE THIS IS RIGHT ALWAYS
alt.a <- get.alt.a(param.table)
param.table$alt.a <- alt.a
Dendro.complete <- do.call(rbind, Dendro.tree)
write.csv(param.table, file = paste(OUTPUT.folder, param.table.name, sep = "/"),
  quote = FALSE, row.names = FALSE)
save(Dendro.complete, file = paste(OUTPUT.folder, Dendro.data.name, sep = "/"))
save(Dendro.split, file = paste(OUTPUT.folder, Dendro.split.name, sep = "/"))
save(Dendro.tree, file = paste(OUTPUT.folder, "Dendro_Tree.Rdata", sep = "/"))
}

.make.adjust <- function(ts.data) {
  which.adjust <- c(which(ts.data$ADJUST == 1), (nrow(ts.data) + 1))
  which.adjust.dbh <- cumsum(ts.data$DBH_TRUE[which.adjust] - ts.data$DBH_TRUE[which.adjust - 1])
  for(r in 1:(length(which.adjust) - 1)) {
      ts.data$DBH_TRUE[which.adjust[r] : (which.adjust[r + 1] - 1)] <-
        ts.data$DBH_TRUE[which.adjust[r] : (which.adjust[r + 1] - 1)] - which.adjust.dbh[r]
  }
  return(ts.data)
}


#' Estimates the maximum likelihood parameters for the 5-parameter logsitic function fit to dendrometer band data.
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer), as well as other designations.
#' @description This is the core fitting algorithm, which uses base \emph{optim()} (with two sequential methods) to
#' estimate the parameters of the logistic function. This follows McMahon and Parker 2014.
#' @return Returns a numeric vector of parameters
#' @export
get.params <- function(ts.data) {
  dbh      <- ts.data$DBH_TRUE
  doy      <- as.integer(ts.data$DOY)
    # date     <- Dates(ts.data$DATE, "%m/%d/%y") # need to fix this for plots
  doy[doy < 75] <- NA
  complete <- complete.cases(dbh, doy)
  dbh      <- as.numeric(dbh[complete])
  doy      <- as.numeric(doy[complete])

  doy.ip.hat <- doy[(which(dbh > mean(dbh)))[1]]

  optim.min <- c((min(dbh, na.rm = TRUE) * 0.99),
    quantile(dbh, 0.5, na.rm = TRUE), -50, 0, 0.01)

  optim.max <- c(min(dbh, na.rm = TRUE), max(dbh, na.rm = TRUE),
    350, 0.1, 15)
  resid.sd <- 0.02

  par.list <- list(L = min(dbh, na.rm = TRUE), K = max(dbh, na.rm = TRUE),
    doy.ip = doy.ip.hat, r = .08, theta = 1)
  params.start <- as.numeric(unlist(par.list))
  params <- params.start

  ##  THESE ARE THE CALLS TO OPTIM  ##
  # weighted values have false ML estimates, so the estimate is re-assessed based on the optimized parameters in an unweighted call

  lg5.output.LB <- optim(par = params, doy = doy, dbh = dbh,
    resid.sd = resid.sd, fn = .get.lg5.ML, method = "L-BFGS-B",
    lower = optim.min, upper = optim.max,
    hessian = FALSE, control = list(trace = 0))
  params <- lg5.output.LB$par

  lg5.output.LB.wt <- optim(par = params, doy = doy, dbh = dbh,
    resid.sd = resid.sd, fn = .get.lg5.ML.wt, method = "L-BFGS-B",
    lower = optim.min, upper = optim.max, hessian = TRUE,
    control = list(trace = 0))
  lg5.output.LB.wt$value <- .get.lg5.ML(params = lg5.output.LB.wt$par, doy, dbh,
    resid.sd = resid.sd)
  lg5.output.NM.wt <- optim(par = params, fn = .get.lg5.ML.wt, resid.sd = resid.sd,
    method = "Nelder-Mead", hessian = TRUE, control = list(trace = 0),
    doy = doy, dbh = dbh)
  lg5.output.NM.wt$value <- .get.lg5.ML(lg5.output.NM.wt$par, doy, dbh,
    resid.sd = resid.sd)

  ## CONSOLIDATE THE RESULTS  ##

  optim.output <- rbind(
    c(params.start, NA),
    # c(lg5.output.LB$par, lg5.output.LB$value),
    c(lg5.output.LB.wt$par, lg5.output.LB.wt$value),
    # c(lg5.output.NM$par, lg5.output.NM$value),
    c(lg5.output.NM.wt$par, lg5.output.NM.wt$value)
    # c(lg5.output.BFGS$par, lg5.output.BFGS$value)
    )

  winner.int <- 1 + match(min(optim.output[-1, 6], na.rm = T), optim.output[-1, 6])
  winner <- rep(".", length = dim(optim.output)[1])
  winner[1] <- NA
  winner[winner.int] <- "*"
  params <- optim.output[winner.int[1], c(1:5) ]
  ab <- optim(par = c(min(dbh), max(dbh)), fn = .lg5.ML.a, resid.sd = resid.sd,
    method = "Nelder-Mead", hessian = TRUE, control = list(trace = 0),
    doy = doy, dbh = dbh, params = params)$par
  ab[2] <- ifelse(ab[2] > params[2], params[2], ab[2])
  return(c(params, ab))

}


#' Extract extra growth metrics.
#' @param param.table The output data.frame from the optimization *seealso* get.optimized.dendro
#' @param Dendro.split A list vector where every entry is a separate time series of a band on a stem in a year.
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#' @param OUTPUT.folder Character string for the target output folder.
#' @param param.table.name Character string for the parameter table output.
#' @param Dendro.data.name Character string for the data object, Dendro.complete, that is a complete data.frame table output.
#' @param Quantile.hull.name Character string for the data object, QH.Rdata, that contains the Quantile Hull results.
#' @description Fits a Quantile Hull to the data and extracts extra growth and
#'    phenology metrics from dendrometer band time series.
#'
#' @return Saves four files to the OUTPUT folder as named above and collected stems.
#' @export
get.extra.metrics <- function(
  param.table,
  Dendro.split,
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "r.squared",
    "ts.sd", "alt.a"),
  OUTPUT.folder      = "OUTPUT",
  param.table.name   = "Param_table_complete.csv",
  Dendro.split.name  = "Dendro_split.Rdata",
  Dendro.data.name   = "Dendro_data.Rdata",
  Quantile.hull.name = "Quantile.hull.Rdata"
  ) {

  QH.list       <- vector("list", nrow(param.table))

  WD.sum          <- rep(NA, nrow(param.table))
  D.sum           <- rep(NA, nrow(param.table))
  RGR             <- rep(NA, nrow(param.table))
  GR              <- rep(NA, nrow(param.table))
  max.grow.day  <- rep(NA, nrow(param.table))
  max.grow.rate <- rep(NA, nrow(param.table))
  Size            <- rep(NA, nrow(param.table))
  Size.alt        <- rep(NA, nrow(param.table))
  fifty.doy       <- rep(NA, nrow(param.table))
  deficit         <- rep(NA, nrow(param.table))
  start.doy       <- rep(NA, nrow(param.table))
  start.doy.alt   <- rep(NA, nrow(param.table))
  stop.doy        <- rep(NA, nrow(param.table))
  doy.05          <- rep(NA, nrow(param.table))
  doy.95          <- rep(NA, nrow(param.table))
  doy.10          <- rep(NA, nrow(param.table))
  doy.90          <- rep(NA, nrow(param.table))
  DATA_SET        <- rep(NA, nrow(param.table))

  new.metrics.df <- data.frame() # this collects new parameters and will be joined
# with params.table at the bottom (in the form of Results.mat)

  seq.l <- 200
  resid.sd <- 0.1
  days <- seq(365)
  deviation.val <- 0.1
  pb <- txtProgressBar(style = 3)

  for(i in 1:nrow(param.table)) {

    setTxtProgressBar(pb, i / nrow(param.table),
      title = NULL, label = sprintf("%s :: %s", Site, Year))
    if(is.null(Dendro.split[[i]])) next
    ts.data <- Dendro.split[[i]]
    ts.data$resids.vec <- rep(NA, times = nrow(ts.data))
    param.data <- param.table[i, ]
    if(any(is.na(param.data[par.names[1:5]]))) {
      Dendro.split[[i]] <- ts.data
      next
    }

    par.col <-  names(param.table) %in% par.names
    params <- param.table[i, par.col]
    params.numeric <- as.numeric(params)

    Site <- as.character(param.data$SITE)
    Year <- as.integer(param.data$YEAR[1])

    site.year <- paste(Site, Year, sep = "-")

    dbh <- ts.data$DBH_TRUE
    doy <- ts.data$DOY


    start.doy[i]     <- pred.doy(params = params, a = params$a)
    stop.doy[i]      <- pred.doy(params = params, a = params$b)
    if(!is.na(params$alt.a) & params$alt.a > 0)
      start.doy.alt[i] <- pred.doy(params = params, a = params$alt.a)

    try.hull <- try(fit.quantile.hull(dbh, doy, params,
      quant = 0.8), TRUE)

    if(class(try.hull) != "try-error") {
      QH.list[[i]] <- try.hull
      D.sum[i] <- sum(try.hull$Deficit)
      WD.sum[i] <- sum(try.hull$Weighted.deficit)
    }

    # Other summary stats
    ts.data$resids.vec <- .get.lg5.resids(params.numeric, doy, dbh)
    ts.data$std.resids <- scale(ts.data$resids.vec)
    RGR[i] <- as.numeric(log(params$b) - log(params$a))
    GR[i] <- as.numeric(params$b - params$a)
    max.grow.day[i] <- max.growth.day(as.numeric(params))
    max.grow.rate[i] <- max.growth.rate(params)
    Size[i] <- params$a
    Size.alt[i] <- params$alt.a
    start.doy[i] <- round(pred.doy(params, params$a))
    stop.doy[i] <- round(pred.doy(params, params$b))
    .deriv <- lg5.deriv(paras = params, days, growth = params$b - params$a)
    fifty.doy[i] <- round(pred.doy(params, mean(c(params$a, params$b))))
    doy.05[i] <- round(pred.doy(params, params$a + 0.05 * GR[i]))
    doy.95[i] <- round(pred.doy(params, params$a + 0.95 * GR[i]))
    doy.10[i] <- round(pred.doy(params, params$a + 0.10 * GR[i]))
    doy.90[i] <- round(pred.doy(params, params$a + 0.90 * GR[i]))

    DATA_SET[i] <- as.character(ts.data$DATA_SET[1])
    Dendro.split[[i]] <- ts.data

  }
  close(pb)

  tmp.df <- data.frame(DATA_SET = DATA_SET, WD = WD.sum, RGR = RGR, GR = GR,
    Max.growth.day = max.grow.day, Max.growth.rate = max.grow.rate,
    Size.a = Size, Size.alt.a = Size.alt,
    Start.g = start.doy, Stop.g = stop.doy,
    Median.g = fifty.doy, DOY.05 = doy.05, DOY.95 = doy.95,
    DOY.10 = doy.10, DOY.90 = doy.90)

  new.metrics.df <- data.frame(cbind(param.table, tmp.df))

  new.metrics.df$GS_Length <- new.metrics.df$Stop.g - new.metrics.df$Start.g
  new.metrics.df$GS.90 <- new.metrics.df$DOY.95 - new.metrics.df$DOY.05
  new.metrics.df$GS.80 <- new.metrics.df$DOY.90 - new.metrics.df$DOY.10

  write.csv(new.metrics.df, file = paste(OUTPUT.folder, param.table.name, sep = "/"),
    quote = FALSE, row.names = FALSE)
  Dendro.complete <- do.call(rbind, Dendro.split)

  save(Dendro.complete, file = paste(OUTPUT.folder, Dendro.data.name, sep = "/"))
  save(Dendro.split, file = paste(OUTPUT.folder, Dendro.split.name, sep = "/"))
  save(QH.list, file = paste(OUTPUT.folder, Quantile.hull.name, sep = "/"))
}




#' Fit Quantile Hull to subset of time series data
#'
#' @param dbh Numeric vector of "TRUE_DBH" values
#' @param doy Integer vector of days of the year dbh values collected.
#' @param params Data.frame (vector with names) containing parameters.
#' @param quant Numeric scalar denoting the quantile threshold above which the hull will be estimated.
#' @param resid.sd Numeric scalar that constrains the optimization algorithm in \emph{optim()}.
#'
#' @return List with output from the optim function and other metrics.
#' @export
fit.quantile.hull <- function(dbh, doy, params, quant = 0.8, resid.sd = 0.02) {
  dbh <- as.numeric(dbh)
  complete <- complete.cases(dbh)
  dbh <- dbh[complete]
  doy <- doy[complete]
  paras <- as.numeric(params[1:5])
  a <- as.numeric(params[6])
  b <- as.numeric(params[7])
  growth <- b - a
  growth.quants <- quantile(seq(a, b, length = 1000),
    c(0.05, 0.1, 0.5, 0.9, 0.95))
  doy.quants <- pred.doy(paras, growth.quants)
  new.a <- growth.quants[1]
  new.a <- growth.quants[1]
  doy.start <- which(doy > doy.quants[1])
  doy.stop <- which(doy < doy.quants[5])

  curve.pure <- intersect(doy.start, doy.stop)

  doyP <- as.numeric(c(pred.doy(paras, growth.quants[1]), doy[curve.pure],
    pred.doy(paras, growth.quants[5])))
  dbhP <- as.numeric(c(growth.quants[1], dbh[curve.pure], growth.quants[5]))
  # residsP <- .get.lg5.resids(params = paras, doy = doy, dbh = dbh)
  residsP <- .get.lg5.resids(params = paras, doy = doyP, dbh = dbhP)
  ln.data <- length(residsP)
  top.resids <- unique(c(1, which(residsP >= quantile(residsP, quant)),
    length(residsP)))
  doyP2 <- doyP[top.resids]
  dbhP2 <- dbhP[top.resids]

  SSstart <- function(doy = doyP2, dbh = dbhP2) {
    lm.fit <- lm(dbhP2 ~ doyP2 + I(doyP2^2))
    new.doy <- seq(range(doyP2)[1], range(doyP2)[2])
    new.dbh <- predict(lm.fit, newdata = data.frame(doyP2 = new.doy),
      type = c("response"))
  }
  optim.min <- as.numeric(c((min(dbhP2, na.rm = TRUE) * 0.99),
    quantile(dbhP2, 0.5, na.rm = TRUE), 0, 0, 0.01))

  optim.max <- as.numeric(c(min(dbhP2, na.rm = TRUE), max(dbhP2, na.rm = TRUE),
    350, 0.1, 15))

  OH.fit <- optim(par = paras, fn = .get.lg5.ML, resid.sd = resid.sd,
    method = "L-BFGS-B", lower = optim.min, upper = optim.max,
    hessian = FALSE, control = list(trace = 0), doy = doyP2, dbh = dbhP2)
  deriv.list <- lg5.deriv(OH.fit$par, doyP, growth = (log(b) - log(a)),
    shift = 0.05)
  resids.hull <- scale(.get.lg5.resids(OH.fit$par, doyP, dbhP))
  weighted.deficit <- resids.hull * deriv.list
  OH.list <- list(doyP2 = doyP2, dbhP2 = dbhP2, doyP = doyP, dbhP = dbhP,
    OH.fit = OH.fit, Derivatives = deriv.list, Deficit = resids.hull,
    Weighted.deficit = weighted.deficit)
  return(OH.list)
}

##############################
## Auxiliary functions
##############################
#' Gets the day of the year for any dbh and paramater combination
#'
#' @param dbh Numeric scalar for diameter at breast height (translated from the GAP_WIDTH)
#' @param params Numeric vector of the parameter values of the LG5
#'
#' @return Returns day of the year for the input dbh and params.
#' @export
get_dcrit <- function(dbh,  params) {
  L <- params[1] # min(dbh, na.rm = T)
  K <- params[2]
  doyip <- params[3]
  r <- params[4]
  theta <- params[5]

  dcrit = doyip * r - theta * log((((K - L) /
    (dbh - L) - 1) * theta) ^ (1 / theta)) / r
  return(dcrit)
}

#' Gets the starting diameter from the year before
#'
#' @param param.tab data.frame of parameter values for every year for a particular TREE_ID
#'
#' @return Numeric vector of diameter values (or codes).
#' @export
get.alt.a <- function(param.tab) {
  TREE.ID.F <- factor(param.tab$UNIQUE_ID, levels = unique(param.tab$UNIQUE_ID))
  param.split <- split(param.tab, f = TREE.ID.F, drop = FALSE)
  new.alt.a <- c()
  dim1 <- c()
  ln.alt.a <- c()
  for(id in 1:length(param.split)) {
    alt.a <- rep(-99, length(param.split[[id]]$a))
    dim1[id] <- dim(param.split[[id]])[1]
    if(length(alt.a) > 1) {
      alt.a[2:length(alt.a)] <- param.split[[id]]$b[1:(length(alt.a) - 1)]
    }
    new.alt.a <- c(new.alt.a, alt.a)
    ln.alt.a[id] <- length(alt.a)
  }
  return(new.alt.a)
}



#' Corrects for the chord
#'
#' @param gw1 Numeric gap width at time 1
#' @param gw2 Numeric gap width at time 2
#' @param dbh1 DBH at time 1
#'
#' @return Numeric scalar for the corrected DBH
#' @export
gettruedbh <- function(gw1, gw2, dbh1, units = "cm") {
  dbh1 <- ifelse(units == "cm", dbh1 * 10, dbh1)
  rhs  <- dbh1 * (pi - asin(gw1 / dbh1))
  #rhs is the length of the dendrometer band at time 1
  dbh2 <- optimize(.difdendro, interval = c(0, dbh1 + 2 * gw2),
    gw2 = gw2, rhs = rhs)
  dbh2.min <- ifelse(units == "cm", dbh2$minimum / 10, dbh2$minimum)
  return(dbh2.min)
}

.difdendro <- function(dbh2, gw2, rhs) {
  lhs <- dbh2 * (pi - asin(gw2 / dbh2))
  return(abs(lhs - rhs))
}

##---------------------------------------------------------------
#' Gets dbh from gap width given org.dbh.
#'
#' @param gap.width Numeric vector of gap width measurements
#' @param org.dbh Numeric scalar of the original dbh of the tree (when bands were installed)
#' @param units Character string denoting the units that the dbh is in. Defaults to "cm".
#' @return Numeric vector of DBH_TRUE values
#' @export
gap2dbh <- function(gap.width, org.dbh, units = "cm") {
  if(units == "cm") {
    dbh.vec <- org.dbh + (((gap.width - gap.width[1]) / pi))
  } else {
    dbh.vec <- org.dbh + (((gap.width - gap.width[1]) / 10 / pi))
  }
  return(dbh.vec)
}

###################################################
### lg5.functions
###################################################
# These functions do many things related to the LG5 function (see vignette)

#' Predicts the dbh from doy and a parameter set
#'
#' @param params Numeric vector of parameter values
#' @param doy Integer vector of days of the year
#'
#' @return Numeric scalar or vector of dbh values
#' @export
lg5.pred <- function(params, doy) {
  paras <- names(params) %in% c("L", "K", "doyip", "r", "theta", "a", "b", "alt.a")
  L <- params[1] # min(dbh, na.rm = T)
  K <- params[2]
  doy.ip <- params[3]
  r <- params[4]
  theta <- params[5]
  dbh <- vector(length = length(doy))
  dbh <- L + ((K - L) / (1 + 1/theta * exp(-(r * (doy - doy.ip) / theta)) ^ theta))
  return(dbh)
}

.get.lg5.ML <- function(params, doy, dbh, resid.sd) {
	pred.dbh <- lg5.pred(params, doy)
	pred.ML <-  -sum(dnorm(dbh, pred.dbh, resid.sd, log = T))
	return(pred.ML)
}

.get.lg5.ML.wt <- function(params, doy, dbh, resid.sd) {
	wts <- 1 / dnorm(abs(seq(-2, 2, length = length(doy))), 0, 1)
	pred.dbh <- lg5.pred(params, doy)
	pred.ML <- -sum((wts * dnorm(dbh, pred.dbh, resid.sd, log = T)))
	return(pred.ML)
}

.get.lg5.resids <- function(params, doy, dbh) {
  para <- as.numeric(params)
  lg5.resid <- dbh - lg5.pred(para, doy)
  return(lg5.resid)
}

#' Predicts the day of the year given a diameter and parameter values
#'
#' @param params Numeric vector of parameter values.
#' @param a Numeric scalar of diameter for which a doy is wanted
#' @param diam.given Not sure
#'
#' @return returns day of year
#' @export
pred.doy <- function(params, a, diam.given = 0) {
	params <- as.numeric(params)
	L <- params[1] # min(dbh, na.rm = T)
	K <- params[2]
	doy.ip <- params[3]
	r <- params[4]
	theta <- params[5]
	a.par <- a
	.expr1 <- (K - L) / (a.par - L)
	.expr2 <- doy.ip * r - theta * log(((.expr1 - 1) * theta) ^ (1/theta))
	.expr3 <- .expr2 / r
	dcrit <- .expr3

	return(dcrit)
}

#' Returns starting and stopping diameter values for a year for a tree
#'
#' @param a Numeric scalar or vector of 2 elements with starting and stopping values
#' @param params Numeric vector of parameters
#' @param doy Integer (or numeric) vector of days of the year for which the curve was fit.
#' @param asymptote which side to get value for (should be removed).
#'
#' @return Numeric scalar (or vector of 2 scalars) for starting and stoping sizes
#' @export
lg5.pred.a <- function (a, params, doy, asymptote = "both") {
	asymptote <- ifelse(length(a) > 1, "both", asymptote)
	L <- params[1] # min(dbh, na.rm = T)
	K <- params[2]
	doy.ip <- params[3]
	r <- params[4]
	theta <- params[5]
	diam <- vector(length = length(doy))
	if(asymptote == "lower") {
		d.crit <- pred.doy(params, a)
		diam[which(doy <= d.crit)] <- a
		diam[which(doy > d.crit)] <- L + ((K - L) / (1 + 1/theta * exp(-(r * (doy[which(doy > d.crit)] - doy.ip) / theta)) ^ theta))
	}else{
		if(asymptote == "upper") {
			d.crit <- pred.doy(params, a)
			diam[which(doy >= d.crit)] <- a
			diam[which(doy < d.crit)] <- L + ((K - L) / (1 + 1/theta * exp(-(r * (doy[which(doy < d.crit)] - doy.ip) / theta)) ^ theta))
		}else{
			if(asymptote == "both") {
				d.crit <- pred.doy(params, a)
				is.nan.dc <- is.nan(d.crit)
				d.crit[is.nan.dc] <- 365
				diam[which(doy <= d.crit[1])] <- a[1]
				diam[which(doy >= d.crit[2])] <- a[2]
				doy.mid <- doy[which(doy > d.crit[1] & doy < d.crit[2])]
				diam[which(doy > d.crit[1] & doy < d.crit[2])] <-
				L + ((K - L) / (1 + 1/theta * exp(-(r * (doy.mid - doy.ip) / theta)) ^ theta))
			}
		}
	}
	return(diam)
}

.lg5.ML.a <- function(a, params, doy, dbh, resid.sd, asymptote = "both") {
	pred.dbh <- lg5.pred.a(a, params, doy)
	pred.ML <- -sum(dnorm(dbh, pred.dbh, resid.sd, log = T))
}


.lg5.QH <- function(paras, doyCP, dbhCP) {
	pred.dbh <- lg5.pred(paras, doyCP)
	pred.ND <-  sum(pred.dbh - dbhCP) # for "Negative Difference"
	return(pred.ND)
}

#' Get quantile hull 'residuals'
#'
#' @param params Numeric vector of parameters.
#' @param doy Day of the Year
#' @param dbh Sequence of DBH values
#' @param log.rate Logical whether rates are logged
#'
#' @return Numeric vector of residuals
#' @export
get.QH.resid <- function(params, doy, dbh, log.rate = FALSE) {
	lg5.pred <- lg5.pred(params, doy)
	rates <- lg5.deriv(params, doy)
	if(log.rate) {
		resids <- lg5.pred - dbh * log(rates + 0.000001)
	}else{
		resids <- lg5.pred - dbh * rates
	}
	return(resids)
}

#' Derivative of LG5 function
#'
#' @param paras Numeric vector of parameters for the lg5 function
#' @param doy Integer vector of days of the year dbh values collected.
#' @param growth Unclear.
#' @param shift Numeric scalar that gives half the size of the window the derivitive is taken over (in days).
#' @description This function takes the numerical derivative.
#' Default values return a single day of growth. Using the 'growth'
#' argument, the derivative can be returned, scaled by annual growth.
#' @return returns numeric vectore of derivatives
#' @export
lg5.deriv <- function(paras, doy, growth = 1, shift = 0.5) {
  paras = as.numeric(paras)
  .loVal <- lg5.pred(paras, (doy - shift))
  .hiVal <- lg5.pred(paras, (doy + shift))
  deriv.lg5 <- (.hiVal - .loVal) / (2 * shift)
  return(deriv.lg5 / growth)
}

#' Find day of the year of maximum growth
#'
#' @param params Numeric vector of parameters.
#'
#' @return Numeric scalar of the fastest day.
#' @export
max.growth.day <- function(params) {
  paras <- as.numeric(params)
  days <- seq(365)
  .deriv <- lg5.deriv(paras, days)
  fastest.day <- max(days[which( .deriv == max(.deriv))], na.rm = TRUE)
  return(fastest.day)
}

#' Find rate of maximum growth
#'
#' @param params Numeric vector of parameters.
#'
#' @return Numeric scalar of the fastest growth rate
#' @export
max.growth.rate <- function(params) {
  paras <- as.numeric(params)
  doy.a <- pred.doy(params, params$a)
  days <- try(seq(round(doy.a), 365), silent = TRUE)
  if(class(days) == "try-error") {
    return(NA)
  } else {
    .deriv <- lg5.deriv(paras, days)
    growth.rate <- max(.deriv, na.rm = TRUE)
    return(growth.rate)
  }
}

