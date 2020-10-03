RankCorrelationPlots = function ( K, Klist, X, corr.type ) {
  par(mfrow = c(K/2, K/2))
  for (k in 1:K) {
    U = Klist[[k]]$U
    V = Klist[[k]]$V
    X.hat = U %*% t(V)
    
    rank.corr = rep(x = NA, times = M)
    for (m in 1:M) {
      rank.corr[m] = cor(X[m, ], X.hat[m, ], method = corr.type)
    }
    summary(rank.corr)
    hist(rank.corr, breaks = 25, xlim = c(-1, 1), col = "darkorange", 
         xlab = paste(corr.type, " correlation", sep = ""), 
         ylab = "Frequency", 
         main = paste("k = ", k, sep = ""))
  }
  par(mfrow = c(1, 1))
}

RankBoxPlots = function ( K, Klist, X ) {
  par(mfrow = c(K/2, K/2))
  for (k in 1:K) {
    U = Klist[[k]]$U
    V = Klist[[k]]$V
    X.hat = U %*% t(V)
    
    v = matrix(data = NA, nrow = M, ncol = N)
    v.hat = v
    for (m in 1:M) {
      # v[m, ] = order(X[m, ])
      # v.hat[m, ] = order(X.hat[m, ])
      v[m, ] = rank(X[m, ])
      v.hat[m, ] = rank(X.hat[m, ])
    }
    cv = c(v)
    cv.hat = c(v.hat)
    boxplot(cv.hat ~ cv, 
            col = "darkorange", 
            outline = FALSE, 
            xlab = "True ranks", 
            ylab = "Fitted ranks", 
            main = paste("k = ", k, sep = ""))
  }
  par(mfrow = c(1, 1))
}

RankStackedPlots = function ( K, Klist, X, nq ) {
  par(mfrow = c(K/2, K/2))
  for (k in 1:K) {
    U = Klist[[k]]$U
    V = Klist[[k]]$V
    X.hat = U %*% t(V)
    
    v = matrix(data = NA, nrow = M, ncol = N)
    v.hat = v
    for (m in 1:M) {
      # v[m, ] = order(X[m, ])
      # v.hat[m, ] = order(X.hat[m, ])
      v[m, ] = rank(X[m, ])
      v.hat[m, ] = rank(X.hat[m, ])
    }
    cv = c(v)
    cv.hat = c(v.hat)
    
    q = matrix(data = NA, nrow = N, ncol = nq-1)
    thres = seq(from = 1/nq, to = 1-1/nq, by = 1/nq)
    for (n in 1:N) {
      q[n, ] = quantile(x = cv.hat[cv == n], probs = thres)
    }
    cols = rev(sequential_hcl(n = nq, "OrRd", rev = TRUE))
    
    plot(1:N, q[, 1], type = 'l', col = cols[1], ylim = c(1, N), xlim = c(1, N),
         xlab = "True ranks", ylab = "Fitted ranks",
         main = paste("k = ", k, sep = ""))
    polygon(c(1:N, rev(1:N)), c(q[, 1], rev(rep(x = 0, times = N))), col = cols[1])
    for (iq in 2:(nq-1)) {
      lines(1:N, q[, iq], type = 'l', col = cols[iq])
      polygon(c(1:N, rev(1:N)), c(q[, iq], rev(q[, iq-1])), col = cols[iq])
    }
    polygon(c(1:N, rev(1:N)), c(rep(x = N, times = N), rev(q[, nq-1])), col = cols[nq])
  }
  par(mfrow = c(1, 1))
}

RMSEPlot = function ( K, X.RS, X ) {
  rmse = rep(x = NA, times = K)
  names(rmse) = as.character(1:K)
  for (k in 1:K) {
    X.hat = X.RS[[k]]$U %*% t(X.RS[[k]]$V)
    X.err = X-X.hat
    rmse[k] = sqrt(sum(X.err^2, na.rm = TRUE)/sum(!is.na(X)))
  }

  barplot(rmse, col = "darkorange",
          xlab = "Number of latent factors used by RS",
          ylab = "RMSE", 
          main = "RMSE vs. number of factors in RS")
}

FitRS = function( K, X ) {
  X.RS = list()
  for (k in 1:K) {
    cat("k = ", k, "\n")
    X.RS[[k]] = funkSVD(x = X, k = k)
  }
  return(X.RS)
}