## This code is a modification of TukeyHSD.R from R stats package
TukeyHSD.aov <-
  function(x, which = seq_along(tabs), ordered = FALSE,
           conf.level = 0.95, ...)
  {
    mm <- model.tables(x, "means")
    if(is.null(mm$n))
      stop("no factors in the fitted model")
    tabs <- mm$tables[-1L]
    tabs <- tabs[which]
    ## mm$n need not be complete -- factors only -- so index by names
    nn <- mm$n[names(tabs)]
    nn_na <- is.na(nn)
    if(all(nn_na))
      stop("'which' specified no factors")
    if(any(nn_na)) {
      warning("'which' specified some non-factors which will be dropped")
      tabs <- tabs[!nn_na]
      nn <- nn[!nn_na]
    }
    out <- vector("list", length(tabs))
    names(out) <- names(tabs)
    MSE <- sum(resid(x)^2)/x$df.residual
    for (nm in names(tabs)) {
      tab <- tabs[[nm]]
      means <- as.vector(tab)
      nms <- if(length(d <- dim(tab)) > 1L) {
        dn <- dimnames(tab)
        apply(do.call("expand.grid", dn), 1L, paste, collapse=":")
      } else names(tab)
      n <- nn[[nm]]
      ## expand n to the correct length if necessary
      if (length(n) < length(means)) n <- rep.int(n, length(means))
      if (as.logical(ordered)) {
        ord <- order(means)
        means <- means[ord]
        n <- n[ord]
        if (!is.null(nms)) nms <- nms[ord]
      }
      center <- outer(means, means, "-")
      ratios <- outer(means, means, "/")
      keep <- lower.tri(center)
      center <- center[keep]
      ratios <- ratios[keep]
      width <- qtukey(conf.level, length(means), x$df.residual) *
        sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep]
      est <- center/(sqrt((MSE/2) * outer(1/n, 1/n, "+"))[keep])
      pvals <- ptukey(abs(est),length(means),x$df.residual,lower.tail=FALSE)
      dnames <- list(NULL, c("diff", "lwr", "upr","p adj", "ratio"))

      if (!is.null(nms)) dnames[[1L]] <- c(outer(nms, nms, paste, sep = "/")[keep])
      out[[nm]] <- array(c(center, center - width, center + width, pvals, ratios),
                         c(length(width), 5), dnames)
      names(means) <- paste("mean (", nms, ")", sep="")
      out[["means"]] <- means
    }
    class(out) <- c("multicomp", "TukeyHSD")
    attr(out, "orig.call") <- x$call
    attr(out, "conf.level") <- conf.level
    attr(out, "ordered") <- ordered
    out
  }
