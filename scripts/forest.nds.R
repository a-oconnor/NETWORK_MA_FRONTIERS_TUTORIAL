forest.nds<-function (x, vi, sei, ci.lb, ci.ub, annotate = TRUE, showweight = FALSE, 
    xlim = NULL, alim = NULL, ylim = NULL, at = NULL, steps = 5, 
    level = 95, digits = 2, refline = 0, xlab = NULL, slab = NULL, 
    ilab = NULL, ilab.xpos = NULL, ilab.pos = NULL, subset = NULL, 
    transf = FALSE, atransf = FALSE, targs, rows = NULL, efac = 1, 
    pch = 15, psize = NULL, cex = NULL, cex.lab = NULL, cex.axis = NULL, 
    ...) 
{
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs)) 
        targs <- NULL
    alpha <- (100 - level)/100
    yi <- x
    if (is.null(attributes(yi)$measure)) {
        measure <- "GEN"
    }
    else {
        measure <- attributes(yi)$measure
    }
    if (hasArg(ci.lb) && hasArg(ci.ub)) {
        if (length(ci.lb) != length(ci.ub)) 
            stop("Length of ci.lb and ci.ub do not match.")
        if (missing(vi) && missing(sei)) {
            vi <- ((ci.ub - ci.lb)/(2 * qnorm(alpha/2, lower.tail = FALSE)))^2
        }
        else {
            if (missing(vi)) 
                vi <- sei^2
        }
        if (length(ci.lb) != length(vi)) 
            stop("Length of vi (or sei) does not match length of (ci.lb, ci.ub) pairs.")
    }
    else {
        if (missing(vi)) {
            if (missing(sei)) {
                stop("Must specify either vi, sei, or (ci.lb, ci.ub) pairs.")
            }
            else {
                vi <- sei^2
                ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
                ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                  sei
            }
        }
        else {
            ci.lb <- yi - qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
            ci.ub <- yi + qnorm(alpha/2, lower.tail = FALSE) * 
                sqrt(vi)
        }
    }
    if (length(yi) != length(vi)) 
        stop("Length of yi does not match the length of vi, sei, or the (ci.lb, ci.ub) pairs.")
    k <- length(yi)
    if (is.null(slab)) 
        slab <- paste("Study ", 1:k)
    if (length(yi) != length(slab)) 
        stop("Number of outcomes does not correspond to the length of the slab argument.")
    if (is.vector(ilab)) 
        ilab <- cbind(ilab)
    if (length(pch) == 1L) 
        pch <- rep(pch, k)
    if (length(pch) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the pch argument.")
    if (!is.null(psize)) {
        if (length(psize) == 1L) 
            psize <- rep(psize, k)
        if (length(psize) != length(yi)) 
            stop("Number of outcomes does not correspond to the length of the psize argument.")
    }
    if (!is.null(subset)) {
        yi <- yi[subset]
        vi <- vi[subset]
        ci.lb <- ci.lb[subset]
        ci.ub <- ci.ub[subset]
        slab <- slab[subset]
        ilab <- ilab[subset, , drop = FALSE]
        pch <- pch[subset]
        psize <- psize[subset]
    }
    k <- length(yi)
    if (is.null(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) 
            rows <- rows:(rows - k + 1)
    }
    if (length(rows) != length(yi)) 
        stop("Number of outcomes does not correspond to the length of the rows argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    ci.lb <- ci.lb[k:1]
    ci.ub <- ci.ub[k:1]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yivi.na <- is.na(cbind(yi, vi))
    if (any(yivi.na)) {
        not.na <- apply(yivi.na, MARGIN = 1, sum) == 0L
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            ci.lb <- ci.lb[not.na]
            ci.ub <- ci.ub[not.na]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in 1:length(rows.na)) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >= 
                  rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail") 
            stop("Missing values in results.")
    }
    k <- length(yi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
        }
    }
    if (showweight) {
        weights <- 1/vi
        weights <- 100 * weights/sum(weights, na.rm = TRUE)
    }
    if (is.null(psize)) {
        if (any(vi <= 0, na.rm = TRUE)) {
            psize <- rep(1, k)
        }
        else {
            wi <- 1/sqrt(vi)
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize, 
                na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        if (showweight) {
            plot.multp.l <- 2
            plot.multp.r <- 2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
            axis.multp.l <- 0.2
            axis.multp.r <- 0.2
        }
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
        axis.multp.l <- 0.2
        axis.multp.r <- 0.2
    }
    if (is.null(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l, 
            max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits)
    }
    if (is.null(alim)) {
        if (is.null(at)) {
            alim <- c(min(ci.lb, na.rm = TRUE) - rng * axis.multp.l, 
                max(ci.ub, na.rm = TRUE) + rng * axis.multp.r)
            alim <- round(alim, digits)
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] > min(yi, na.rm = TRUE)) {
        alim[1] <- min(yi, na.rm = TRUE)
    }
    if (alim[2] < max(yi, na.rm = TRUE)) {
        alim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (is.null(ylim)) {
        ylim <- c(0.5, k + 3)
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        at <- seq(alim[1], alim[2], length = steps)
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits, 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits, format = "f")
        }
    }
    else {
        at.lab <- round(at.lab, digits)
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 1] <- 1
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", 
        yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    lheight <- strheight("O")
    cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 * 
        k * lheight), 1)
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab)) 
            cex.lab <- cex
        if (is.null(cex.axis)) 
            cex.axis <- cex
    }
    if (is.null(cex.lab)) 
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis)) 
        cex.axis <- par("cex") * cex.adj
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis, 
        ...)
    if (is.null(xlab)) 
        xlab <- .setxlab(measure, transf.char, atransf.char, 
            gentype = 1)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, 
            line = 2.5, cex = cex.lab, ...)
    for (i in 1:k) {
        if (is.na(yi[i]) || is.na(ci.lb)[i] || is.na(ci.ub)[i]) 
            next
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i], 
            alim[2]), rows[i], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex * 
                efac, ci.lb[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[1]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex * 
                efac, ci.ub[i], rows[i] + (height/150) * cex * 
                efac, ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex * 
                (xlim[2] - xlim[1]), alim[2]), y = c(rows[i], 
                rows[i] + (height/150) * cex * efac, rows[i] - 
                  (height/150) * cex * efac, rows[i]), col = "black", 
                ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos)) 
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != dim(ilab)[2]) 
            stop("Number of 'ilab' columns does not match length of 'ilab.xpos' argument.")
        for (l in 1:NCOL(ilab)) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l], 
                cex = cex, ...)
        }
    }
    if (is.numeric(refline)) 
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2, 
            lty = "dotted", ...)
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                annotext <- round(cbind(sapply(yi, atransf), 
                  sapply(ci.lb, atransf), sapply(ci.ub, atransf)), 
                  digits)
            }
            else {
                annotext <- round(cbind(sapply(yi,atransf,targs), 
                  sapply(ci.lb, atransf,targs), sapply(ci.ub, 
                    atransf,targs)),digits)
            }
        }
        else {
            annotext <- round(cbind(yi,ci.lb,ci.ub),digits)
        }
        if (showweight) {
            annotext <- cbind(round(weights,digits),annotext)
            annotext <- matrix(apply(annotext,2,format,nsmall=digits), 
                ncol = 4)
            annotext <- cbind(annotext[,1],annotext[,2],'(',annotext[,3],", ",annotext[,4],')')
        }
        else {
            annotext <- matrix(apply(annotext,2,format, nsmall = digits), ncol = 3)
            annotext <- cbind(annotext[,1],'(',annotext[,2],", ",annotext[,3],")")
        }
	 #annotext[-1,3]<-substr(annotext[-1,3],2,5)
	
        annotext <- apply(annotext,1, paste,collapse="")
        text(x = xlim[2], rows, labels = annotext, pos = 2, cex = cex, 
            ...)
    }
    points(yi,rows, pch = pch, cex = cex * psize, ...)
    invisible()
}

