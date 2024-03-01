# Post-processing functions to be sourced

move.mean <- function(x) {
  x.1 <- rbind(x[-length(x)], x[seq(2, length(x))])
  y <- apply(x.1, MAR = 2, mean)
  return(y)
}

run.tot.plot <- function(main.title = "", suppress.lines = TRUE) {
  plot(range(dff$doy, na.mr = TRUE), c(0, 1), xlab = "Day of the year",
    ylab = "Proportion of growth", type = "b",
      pch = 19, col = "white", lwd = 0.2,
      main = main.title, cex = 0.1)

  for(t in 1:dim(res.df)[1]) {
    if(t %% 100 == 0) print(t)
    obj.tree <- res.df[t,]
    params <- as.numeric(res.df[t, c("L", "K", "ip", "r", "theta")])

    dbh.delt <- lg5.pred.a(a = as.numeric(obj.tree[c("Start", "Stop")]),
      params = params,
      doy, asymptote = "both")

      dbh.norm <- diff(dbh.delt) / max(diff(dbh.delt))
      curve.mat[t, ] <- dbh.norm
      if(!suppress.lines) {
      lines(doy.mm, dbh.norm, type = "l",
        pch = 19, col = "gray75", lwd = 0.2, cex = 0.1)
      }

  }
  return(curve.mat)
}


make.line <- function(subset.t, cols = "tomato", lty.s = 1, lwd.s = 2) {
  lines(doy.mm, apply(subset.t, 2, median, na.rm = TRUE), lwd = lwd.s,
    col = cols, lty = lty.s)
}


make.polygon <- function(doy.mm, curve.mat, hi, lo, cols,
  alpha = 50, bord = NULL) {
  curve.50 <- apply(curve.mat, 2, median, na.rm = TRUE)
  curve.hi <- apply(curve.mat, 2, quantile, prob = hi, na.rm = TRUE)
  curve.lo <- apply(curve.mat, 2, quantile, prob = lo, na.rm = TRUE)

  polygon(c(doy.mm, rev(doy.mm)), c(curve.hi, rev(curve.lo)),
    col = paste(cols, alpha, sep = ""), border = bord)
  lines(doy.mm, curve.50, col = cols, lwd = 2)
}

make.ba <- function(d) {
  y <- pi * (d * 0.5) ^ 2
  return(y)
}

stand <- function(x) {
  y <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  return(y)
}


#  plotting function for any grouped variable.

make.line.plot <- function(dat, ind = list(res.df$Year, res.df$Site),
  main.title = "", make.date = "DOY", probs.vec = c(0.05, 0.25, 0.5, 0.75, 0.95) ) {

  dat.quant <- sapply(split(dat, ind, drop = TRUE), quantile,
    probs = prob.vec)
  gp.names <- matrix(unlist(strsplit(colnames(dat.quant), "\\.")),
    nrow = length(ind), byrow = FALSE)
  n.gps <- sapply(sapply(ind, unique), length)

  if(n.gps[2] > 8) {
    cols <- rep(brewer.pal(10, "Spectral"),
      length = n.gps[2])[as.factor(gp.names[2,])]
  } else {
    if(n.gps[2] == 2) {
      cols <- brewer.pal(3, "Dark2")[2:3][as.factor(gp.names[2,])]
    } else {
      cols <- brewer.pal(n.gps[2], "Dark2")[as.factor(gp.names[2,])]
    }
  }
  cols <-  brewer.pal(3, "Dark2")[3]
  gps <- dim(dat.quant)[2]
  y.ax <- rev(seq(gps))
  plot(dat.quant["50%", ], y.ax, xlim = range(dat.quant), pch = 19,
    xlab = "Day of the Year", axes = 0, col = cols, ylab = "",
      main = main.title)

  if(make.date == "DATE") {
    doy.seq <- seq(min(dat.quant), max(dat.quant), length = 5)
    mdy <- month.day.year(doy.seq, origin = c(1, 1, 2011))
    x.ax <- paste(mdy[[1]], mdy[[2]], sep = "/")
    axis(1, x.ax, at = doy.seq)
  } else {
    axis(1)
  }


  axis(2, at = y.ax, labels = gp.names[1,], las = 2, cex.axis = 0.5)
  box()
  segments(dat.quant[2, ], y.ax, dat.quant[4, ], y.ax, lwd = 3, col = cols)
  segments(dat.quant[1, ], y.ax, dat.quant[5, ], y.ax, lwd = 1,
    adj = 1, col = cols)
}

make.dot.plot <- function(dff1) {
  doy.mult <- dff1$doy + 365 * (1 + (dff1$YEAR - min(dff1$YEAR)))

  plot(doy.mult, dff1$DBH, xlab = "Day of the year", ylab = "DBH (cm)",
      pch = 19, col = cols[as.factor(dff1$YEAR)],
      main = sprintf("Growth of %s at %s", dff1$TREE_ID[1],
          dff1$SITE[1]),
      cex = 1)

  # days <- seq(365)
  # lines(days, lg5.pred(params = params.best, doy = days),
  #     col = "tomato", lty = 1, lwd = 1)
}


make.quant.reg <- function(res.df, group.var = c("Site", "Year"), var.1, var.2,
  prob.vec = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  require(plyr)
  eval(parse(text = paste("quant.dat <- ddply(.data = res.df, .variables =
    group.var, .fun = summarize, quantile.1 = quantile(", var.1, ", probs =
    prob.vec),
    quantile.2 = quantile(", var.2, ", probs = prob.vec))", sep = "")))

  quant.dat$quants <- rep(as.character(prob.vec), times = 8)
  quant.sp <- split(quant.dat, quant.dat$quants)
  return(quant.sp)
}

make.quant.pic <- function(quant.sp, var.1, var.2, buffx = 2, buffy = 2,
  x.lab = var.2, y.lab = var.1, maintitle = "", pdf.out = FALSE,
  ab.line = TRUE) {

  med.lo.x <- quant.sp[[3]]$quantile.1 - buffx
  med.hi.x <- quant.sp[[3]]$quantile.1 + buffx
  med.lo.y <- quant.sp[[3]]$quantile.2 - buffy
  med.hi.y <- quant.sp[[3]]$quantile.2 + buffy
  if(pdf.out == TRUE) {
    pdf(paste("FIGURES/", var.1, var.2, ".pdf", sep = "_"),
      height = 5, width = 6)
  }

  pchs <- seq(21, 25)
  cols <- brewer.pal(5, "Dark2")

  plot(quant.sp[[3]]$quantile.1, quant.sp[[3]]$quantile.2,
    col = cols[quant.sp[[3]]$Site], pch = pchs[quant.sp[[3]]$Year],
    ylab = y.lab, xlab = x.lab,
    xlim = c(80, 160),  #range(c(quant.sp[[1]]$quantile.1, quant.sp[[5]]$quantile.1)),
    ylim = range(c(quant.sp[[1]]$quantile.2, quant.sp[[5]]$quantile.2)),
    log = "")

  legend("topright", legend = levels(quant.sp[[1]]$Site), col = cols, lwd = 2,
    cex = 1)
  legend("right", legend = unique(quant.sp[[1]]$Year), col = 1, pch = pchs,
    cex = 1)

  segments(quant.sp[[1]]$quantile.1, quant.sp[[3]]$quantile.2, med.lo.x,
    quant.sp[[3]]$quantile.2, col = cols[quant.sp[[3]]$Site], lty = 3,
    lwd = 0.5)

  segments(quant.sp[[5]]$quantile.1, quant.sp[[3]]$quantile.2, med.hi.x,
    quant.sp[[3]]$quantile.2, col = cols[quant.sp[[3]]$Site], lty = 3,
    lwd = 0.5)

  segments(quant.sp[[3]]$quantile.1, quant.sp[[1]]$quantile.2,
    quant.sp[[3]]$quantile.1, med.lo.y, col = cols[quant.sp[[3]]$Site],
    lty = 3, lwd = 0.5)

  segments(quant.sp[[3]]$quantile.1, quant.sp[[5]]$quantile.2,
    quant.sp[[3]]$quantile.1, med.hi.y, col = cols[quant.sp[[3]]$Site],
    lty = 3, lwd = 0.5)


  segments(quant.sp[[2]]$quantile.1, quant.sp[[3]]$quantile.2, med.lo.x,
    quant.sp[[3]]$quantile.2, col = cols[quant.sp[[3]]$Site], lty = 1,
    lwd = 1)

  segments(quant.sp[[4]]$quantile.1, quant.sp[[3]]$quantile.2, med.hi.x,
    quant.sp[[3]]$quantile.2, col = cols[quant.sp[[3]]$Site], lty = 1,
    lwd = 1)

  segments(quant.sp[[3]]$quantile.1, quant.sp[[2]]$quantile.2,
    quant.sp[[3]]$quantile.1, med.lo.y, col = cols[quant.sp[[3]]$Site],
    lwd = 1, lty = 1)

  segments(quant.sp[[3]]$quantile.1, quant.sp[[4]]$quantile.2,
    quant.sp[[3]]$quantile.1, med.hi.y, col = cols[quant.sp[[3]]$Site],
    lwd = 1, lty = 1)

  med.mod <- lm(quant.sp[[3]]$quantile.2 ~ quant.sp[[3]]$quantile.1)
  if(ab.line == TRUE)  abline(med.mod, col = "red")

  points(quant.sp[[3]]$quantile.1, quant.sp[[3]]$quantile.2,
    col = cols[quant.sp[[3]]$Site], pch = pchs[quant.sp[[3]]$Year], cex = 1.5)

  if(pdf.out == TRUE) dev.off()

  return(med.mod)
}


#' Shift the day of the year that defines annual growth for a
#'  time series of an individual
#'
#' @param dff Data frame for shift
#' @param doy.shift Integer value of the day of the year that becomes the first.
#'
#' @return Data frame with a new column of the days of the year of the shift
#' @export

shift.year <- function(dff, doy.shift) {
  new.doy <- ifelse(dff$doy < doy.shift, 365 - doy.shift + dff.doy, 
                dff$doy - doy.shift)
}
