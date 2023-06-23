eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- q.mod <- rep(NA,n)
  ids <- which(!is.na(p.ord))
  k <- 0
  q1 <- qvalue(p.ord[!is.na(p.ord)])$q
  q2 <- qvalue(p.mod[!is.na(p.mod)])$q
  for(i in ids){ 
    k <- k+1
    q.ord[i] <- q1[k] 
    q.mod[i] <- q2[k]
  }
  results.eb.mult <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.0, df.r, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}