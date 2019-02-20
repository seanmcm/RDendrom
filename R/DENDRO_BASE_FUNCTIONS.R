###################################################
###  functions for Intensive dendrometer workflow -- OPTIM
###################################################

#' Estimates the maximum likelihood parameters for the 5-parameter logsitic function fit to dendrometer band data.
#'
#' @param ts.data  A dataframe of a time series of a single tree in a year.
#' Must have column variables \emph{DBH_TRUE} (numeric) and \emph{DOY} (integer), as well as other designations.
#' @description This is the core fitting algorithm, which uses base \emph{optim()} (with a series of different methods) to
#' estimate the parameters of the logistic function. This follows McMahon and Parker 2015.
#' @return Returns a numeric vector of parameters
#' @export Nothing
#'
#' @examples none
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
    resid.sd = resid.sd, fn = get.lg5.ML, method = "L-BFGS-B",
    lower = optim.min, upper = optim.max,
    hessian = FALSE, control = list(trace = 0))
  params <- lg5.output.LB$par

  lg5.output.LB.wt <- optim(par = params, doy = doy, dbh = dbh,
    resid.sd = resid.sd, fn = get.lg5.ML.wt, method = "L-BFGS-B",
    lower = optim.min, upper = optim.max, hessian = TRUE,
    control = list(trace = 0))
  lg5.output.LB.wt$value <- get.lg5.ML(params = lg5.output.LB.wt$par, doy, dbh,
    resid.sd = resid.sd)
  lg5.output.NM.wt <- optim(par = params, fn = get.lg5.ML.wt, resid.sd = resid.sd,
    method = "Nelder-Mead", hessian = TRUE, control = list(trace = 0),
    doy = doy, dbh = dbh)
  lg5.output.NM.wt$value <- get.lg5.ML(lg5.output.NM.wt$par, doy, dbh,
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
  ab <- optim(par = c(min(dbh), max(dbh)), fn = lg5.ML.a, resid.sd = resid.sd,
    method = "Nelder-Mead", hessian = TRUE, control = list(trace = 0),
    doy = doy, dbh = dbh, params = params)$par
  ab[2] <- ifelse(ab[2] > params[2], params[2], ab[2])
  return(c(params, ab))

}

## Window differences are not differences in arcs but chords -----------
# Using KC's code as correction formulas -------------
# INPUT: initial window measurement (c1),
# final window measurement (c2), initial diameter measurement (d1)
.difdendro <- function(dbh2, gw2, rhs) {
  lhs <- dbh2 * (pi - asin(gw2 / dbh2))
  return(abs(lhs - rhs))
}

gettruedbh2 <- function(gw1, gw2, dbh1) {
	rhs  <- dbh1 * (pi - asin(gw1 / dbh1))
	#rhs is the length of the dendrometer band at time 1
	dbh2 <- optimize(.difdendro, interval = c(0, dbh1 + 2 * gw2),
		gw2 = gw2, rhs = rhs)
	return(dbh2$minimum)
}


##---------------------------------------------------------------
gap2dbh <- function(gap.width, org.dbh) {
  # computes a vector of dbh values given gap width
  # and starting dbh.
  #
  # Args:
  #   gap.width
  #   org.dbh
  #
  # Returns:
  #   a vector of dbh values
  # TODO: treat gap as a chord and not an arc and translate
  #       to dbh that way

  dbh.vec <- org.dbh + (((gap.width - gap.width[1]) / 10 / pi))
  return(dbh.vec)
}

###################################################
### lg5.functions
###################################################
# This function predicts a diameter given the day of the year and a vector of parameters for the lg5 model.
# It is called by the lg5.ss and lg5.plot functions.

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

get.lg5.ML <- function(params, doy, dbh, resid.sd) {
	pred.dbh <- lg5.pred(params, doy)
	pred.ML <-  -sum(dnorm(dbh, pred.dbh, resid.sd, log = T))
	return(pred.ML)
}

get.lg5.ML.wt <- function(params, doy, dbh, resid.sd) {
	wts <- 1 / dnorm(abs(seq(-2, 2, length = length(doy))), 0, 1)
	pred.dbh <- lg5.pred(params, doy)
	pred.ML <- -sum((wts * dnorm(dbh, pred.dbh, resid.sd, log = T)))
	return(pred.ML)
}

get.doy <- function(x) {
	names.data <- names(x)
	doy.1 <- as.numeric(unlist(strsplit(names.data[grep("X", names.data)], "X")))
	doy <- doy.1[!is.na(doy.1)]
	return(doy)
}

get.lg5.resids <- function(params, doy, dbh) {
  para <- as.numeric(params)
	lg5.resid <- dbh - lg5.pred(para, doy)
	return(lg5.resid)
}

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

lg5.ML.a <- function(a, params, doy, dbh, resid.sd, asymptote = "both") {
	pred.dbh <- lg5.pred.a(a, params, doy)
	pred.ML <- -sum(dnorm(dbh, pred.dbh, resid.sd, log = T))
}

make.seq <- function(param, params, deviation = 0.1, len.seq = 50, CI = c(0, 0), asymptote = "lower", min.val = NULL, max.val = NULL) {
	if(asymptote == "lower") {
		if(CI[1] > 0) {
			lower.lim <- max(min.val, CI[1] * (1 - deviation), na.rm = T)
			par.seq <- seq(lower.lim, (CI[2] * (1 + deviation)), length = len.seq)
		}else{
			par.seq <- seq(min.val, (param + deviation * param), length = len.seq)

		}
	}else{
		if(CI[1] > 0) {
			upper.lim <- min(max.val, CI[2] * (1 + deviation), na.rm = T)
			par.seq <- seq((CI[1] * (1 - deviation)), upper.lim, length = len.seq)
		}else{
			par.seq <- seq((param - deviation * param), max.val, length = len.seq)
		}
	}
	return(par.seq)
}

lg5.QH <- function(paras, doyCP, dbhCP) {
	pred.dbh <- lg5.pred(paras, doyCP)
	pred.ND <-  sum(pred.dbh - dbhCP) # for "Negative Difference"
	return(pred.ND)
}

get.QH.resid <- function(rate, params, doy, dbh, log.rate = FALSE) {
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
#' @param shift Numeric scaler that gives half the size of the window the derivitive is taken over (in days).
#' @description This function takes the numerical derivative.
#' Default values return a single day of growth. Using the 'growth'
#' argument, the derivative can be returned, scaled by annual growth.
#' @return returns numeric vectore of derivatives
#' @export Exports nothing
#'
#' @examples none
lg5.deriv <- function(paras, doy, growth = 1, shift = 0.5) {
  paras = as.numeric(paras)
	.loVal <- lg5.pred(paras, (doy - shift))
	.hiVal <- lg5.pred(paras, (doy + shift))
	deriv.lg5 <- (.hiVal - .loVal) / (2 * shift)
	return(deriv.lg5 / growth)
}

max.growth.day <- function(params) {
  paras <- as.numeric(params)
	days <- seq(365)
	.deriv <- lg5.deriv(paras, days)
	fastest.day <- max(days[which( .deriv == max(.deriv))],
    start.day, na.rm = TRUE)
	return(fastest.day)
}

max.growth.rate <- function(params) {
  paras <- as.numeric(params)
	days <- seq(round(pred.doy(params, params$a)), 365)
	.deriv <- lg5.deriv(paras, days)
	growth.rate <- max(.deriv, na.rm = TRUE)
	return(growth.rate)
}

#' Fit Quantile Hull to subset of time series data
#'
#' @param dbh Numeric vector of "TRUE_DBH" values
#' @param doy Integer vector of days of the year dbh values collected.
#' @param params Data.frame (vector with names) containing parameters.
#' @param quant Numeric scaler denoting the quantile threshold above which the hull will be estimated.
#' @param resid.sd Numeric scaler that constrains the optimization algorithm in \emph{optim()}.
#'
#' @return List with output from the optim function and other metrics.
#'
#' @examples none
outer.hull <- function(dbh, doy, params, quant = 0.8, resid.sd = 0.02) {
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
  residsP <- get.lg5.resids(params = paras, doy = doyP, dbh = dbhP)
	residsP <- get.lg5.resids(params = paras, doy = doy, dbh = dbh)
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

	OH.fit <- optim(par = paras, fn = get.get.lg5.ML, resid.sd = resid.sd,
		method = "L-BFGS-B", lower = optim.min, upper = optim.max,
		hessian = FALSE, control = list(trace = 0), doy = doyP2, dbh = dbhP2)
	deriv.list <- lg5.deriv(OH.fit$par, doyP, growth = (log(b) - log(a)),
		shift = 0.05)
	resids.hull <- get.lg5.resids(OH.fit$par, doyP, dbhP)
	weighted.deficit <- resids.hull * deriv.list
	OH.list <- list(doyP2 = doyP2, dbhP2 = dbhP2, doyP = doyP, dbhP = dbhP,
		OH.fit = OH.fit, Derivatives = deriv.list, Deficit = resids.hull,
		Weighted.deficit = weighted.deficit)
	return(OH.list)
}

fit.outer.hull <- function(dbh, doy.full, params, quant = 0.8) {
	dbh <- as.numeric(dbh)
	complete <- complete.cases(dbh)
	dbh <- dbh[complete]
	doy <- doy.full[complete]
	out.fit <- outer.hull(dbh, doy, params)
}


get_dcrit <- function( diam,  params) {
  L <- params[1] # min(dbh, na.rm = T)
  K <- params[2]
  doyip <- params[3]
  r <- params[4]
  theta <- params[5]

  dcrit = doyip * r - theta * log((((K - L) /
    (diam - L) - 1) * theta) ^ (1 / theta)) / r
  return(dcrit)
}

get.alt.a <- function(param.tab) {

  TREE.ID.F <- factor(param.tab$UNIQUE_ID, levels = unique(param.tab$UNIQUE_ID))
  param.split <- split(param.tab, f = TREE.ID.F, drop = TRUE)
  new.alt.a <- c()
  for(id in 1:length(param.split)) {
    alt.a <- rep(-99, length(param.split[[id]]$a))
    if(length(alt.a) > 1) {
      alt.a[2:length(alt.a)] <- param.split[[id]]$b[1:(length(alt.a) - 1)]
    }
    new.alt.a <- c(new.alt.a, alt.a)
  }
  return(new.alt.a)
}


#############################
#############################

#' Fit the LG5 function
#'
#' @param INPUT.dendro A data.frame of a specific form (see @details)
#' @param no.neg.growth logical denoting whether to skip any time series with negative trends (from a call to \emph{lm()})
#' @param cutoff Integer denoting the minimum sample size required to perform the optimization call.
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#' @param OUTPUT.folder Character string for the target output folder.
#' @param param.table.name Character string for the parameter table output.
#' @param Dendro.data.name Character string for the data object, Dendro.complete, that is a complete data.frame table output.
#' @param Dendro.split.name Character string for the data object, Dendro.split, that is a list vector where
#' every entry is a separate time series of a band on a stem in a year.
#'
#' @description Estimate LG5 fit, structure output, and get summary statistics for dendrometer band time series.
#'
#' @return Returns nothing.
#' @export Objects Saves four files to the OUTPUT folder as named above and collected stems.
#'
#' @examples none
get.optimized.dendro <- function(INPUT.dendro,
  no.neg.growth = TRUE, cutoff = 9,
  par.names = c("L", "K", "doyip", "r", "theta", "a", "b", "r.squared", "ts.sd"),
  OUTPUT.folder = "OUTPUT",
  param.table.name = "Param_table.csv",
  Dendro.data.name = "Dendro_data.Rdata",
  Dendro.split.name = "Dendro_data_split.Rdata") {

  options(warn = -1)
  TREE.ID <- paste(INPUT.dendro$SITE, INPUT.dendro$TREE_ID,
    INPUT.dendro$BAND_NUM, sep = "_")
  TREE.ID.YR <- paste(as.character(INPUT.dendro$SITE), as.character(INPUT.dendro$TREE_ID),
    as.character(INPUT.dendro$YEAR), sep = "_")
  Dendro.split <- vector("list", length = length(unique(TREE.ID.YR)))
  ind.dendro <- split(INPUT.dendro, f = TREE.ID)

  n.obs <- length(ind.dendro)

  param.table <-c()
  pb <- txtProgressBar(style = 3)
  Dendro.tree <- vector("list", length(ind.dendro))
  ct <- 1
  for(i in 1:n.obs) {
    setTxtProgressBar(pb, i / n.obs, title = NULL, label = NULL)
    ind.data <- ind.dendro[[i]] # loads an individual (multiple years)
    OD <- which(!is.na(ind.data$ORG_DBH))
    ind.data$ORG_DBH <- rep(ind.data$ORG_DBH[OD[1]], length(ind.data$ORG_DBH))
    ind.data$DBH <- gap2dbh(ind.data$GAP_WIDTH, ind.data$ORG_DBH[1])
    ind.data$DBH_TRUE[1] <- ind.data$ORG_DBH[1]
    for(v in 2:length(ind.data$DBH)){
      ind.data$DBH_TRUE[v] <- gettruedbh2(gw1 = 0.1 * ind.data$GAP_WIDTH[v - 1],
        gw2 =  0.1 * ind.data$GAP_WIDTH[v], dbh1 = ind.data$DBH_TRUE[v - 1])
    }

    # Here we fix the new band, slipped band problems (these will be separated)
    # if(any(ind.data$NEW_BAND == 1 | ind.data$ADJUST == 1, na.rm = TRUE)) {
    #   which.NB.tot <- which(ind.data$NEW_BAND == 1 | ind.data$ADJUST == 1)
    #   if(length(which.NB.tot) == 1) {
    #     which.NB <- which.NB.tot
    #     NB.ind <- c(which.NB:dim(ind.data)[1])

    #     NB.date <- as.Date(ind.data$DATE[which.NB], format = "%m/%d/%y")
    #     pre.NB.date <- as.Date(ind.data$DATE[1:(which.NB - 1)], format = "%m/%d/%y")
    #     start.date.dbh <- ind.data$DBH_TRUE[max(which(pre.NB.date <= NB.date))]
    #     NB.data <- ind.data$DBH_TRUE[NB.ind]
    #     NB.data.new <- NB.data - NB.data[1] + start.date.dbh
    #     ind.data$DBH_TRUE[NB.ind] <- NB.data.new
    #     } else {
    #       which.NB.tot <- c(which.NB.tot, length(ind.data$FLAG))
    #       for(j in 1:(length(which.NB.tot) - 1)) {
    #         which.NB <- which.NB.tot[j]
    #         NB.ind   <- c(which.NB:which.NB.tot[j + 1])

    #         NB.date        <- as.Date(ind.data$DATE[which.NB], format = "%m/%d/%y")
    #         pre.NB.date    <- as.Date(ind.data$DATE[1:(which.NB - 1)], format = "%m/%d/%y")
    #         start.date.dbh <- ind.data$DBH_TRUE[max(which(pre.NB.date <= NB.date))]
    #         NB.data        <- ind.data$DBH_TRUE[NB.ind]
    #         NB.data.new    <- NB.data - NB.data[1] + start.date.dbh
    #         ind.data$DBH_TRUE[NB.ind] <- NB.data.new

    #       }

    #     }

    #   }
      # if(any(ind.data$is.B.band != "0")) {
      #   org.B.band <- ind.data$DBH_TRUE[min(which(ind.data$is.B.band != 0)) - 1]
      #   band.split <- split(ind.data, ind.data$is.B.band)
      #   band.split[[2]]$DBH_TRUE <- gap2dbh(band.split[[2]]$GAP_WIDTH, org.B.band)
      #   ind.data <- unsplit(band.split, f = ind.data$is.B.band)
      # }

    Dendro.tree[[i]] <- ind.data
    ind.year <- split(ind.data,
      f = list(YEAR = ind.data$YEAR, BAND_NUM = ind.data$BAND_NUM), drop = TRUE)

    params <- rep(NA, 7)
    r.squared <- 0
    param.mat <- matrix(NA, length(ind.year), length(par.names))

    for(t in 1:length(ind.year)) {

      ts.data.tmp <- ind.year[[t]]
      ts.data <- subset(ts.data.tmp, REMOVE == 0)
      ts.sd <- sd(ts.data$DBH_TRUE, na.rm = TRUE)
      if (sum(!is.na(ts.data$DBH_TRUE)) < cutoff) {
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
      YEAR = ts.data$YEAR[1], TREE_ID = ts.data$TREE_ID[1],
      BAND_NUM = ts.data$BAND_NUM[1], UNIQUE_ID = ts.data$UNIQUE_ID[1],
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



#' Extract extra growth metrics.
#' @param param.table The output data.frame from the optimization @seealso get.optimized.dendro
#' @param Dendro.split A list vector where every entry is a separate time series of a band on a stem in a year.
#' @param par.names  A character vector of the column names of the \emph{params} vector argument.
#' This is used to pull out actual parameters from other information.
#' @param OUTPUT.folder Character string for the target output folder.
#' @param param.table.name Character string for the parameter table output.
#' @param Dendro.data.name Character string for the data object, Dendro.complete, that is a complete data.frame table output.
#' @param Quantile.hull.name Character string for the data object, QH.Rdata, that contains the Quantile Hull results.
#'
#' @return Returns nothing.
#' @export Objects Saves four files to the OUTPUT folder as named above and collected stems.

#'
#' @return nothing
#' @export Objects Saves four files to the OUTPUT folder as named above and collected stems.
#' @description Fits a Quantile Hull to the data and extracts extra growth and phenology metrics from dendrometer band time series.
#' @examples none
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
  max.growth.day  <- rep(NA, nrow(param.table))
  max.growth.rate <- rep(NA, nrow(param.table))
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

    resids.vec    <- vector("numeric", length = nrow(ts.data))

    Site <- as.character(param.data$SITE)
    Year <- as.integer(param.data$YEAR[1])

    site.year <- paste(Site, Year, sep = "-")

    dbh <- ts.data$DBH_TRUE
    doy <- ts.data$DOY


    start.doy[i]     <- pred.doy(params = params, a = params$a)
    stop.doy[i]      <- pred.doy(params = params, a = params$b)
    if(!is.na(params$alt.a) & params$alt.a > 0)
      start.doy.alt[i] <- pred.doy(params = params, a = params$alt.a)

    try.hull <- try(fit.outer.hull(dbh, doy, params,
      quant = 0.8), TRUE)

    if(class(try.hull) != "try-error") {
      QH.list[[i]] <- try.hull
      D.sum[i] <- sum(QH.list$Deficit)
      WD.sum[i] <- sum(QH.list$Weighted.deficit)
    }

    # Other summary stats
    ts.data$resids.vec <- get.lg5.resids(params.numeric, doy, dbh)
    ts.data$sd.resids <- scale(ts.data$resids.vec)
    RGR[i] <- as.numeric(log(params$b) - log(params$a))
    GR[i] <- as.numeric(params$b - params$a)
    max.growth.day[i] <- as.numeric(params)
    max.growth.rate[i] <- as.numeric(params)
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

    Dendro.split[[i]] <- ts.data

  }
  close(pb)

  tmp.df <- data.frame(WD = WD.sum, RGR = RGR, GR = GR,
    Max.growth.day = max.growth.day, Max.growth.rate = max.growth.rate,
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

