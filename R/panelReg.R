# Various algorithms for semiparameteric panel count regression
# List of internal functions :
#   doPanelFit.AEE
#   doPanelFit.AEEX
#   doPanelFit.HWZ
#   doPanelFit.SW
#   doPanelFit.Engine.Bootstrap
#   doPanelFit.AEE.Impute


##############################################################################
# Augmented Estimating Equations (AEE)
##############################################################################
doPanelFit.AEE <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    eStep <- function(lambda) {
        e <- matrix(0, N, K)

        for (i in 1:N) {
            end <- which(!is.na(panelMatrix[i, ]))
            start <- c(1, head(end, -1) + 1)

            for (j in which(panelMatrix[i, end] > 0)) {
                sq <- seq(start[j], end[j])
                e[i, sq] <- panelMatrix[i, end[j]] * lambda[sq] / sum(lambda[sq])
            }
        }
        e
    }

    # ncol(X) dimensional nonlinear equation
    f <- function(beta, e) {
        lambda <- c(colSums(e)) / c(t(r) %*% exp(X %*% beta))
        c(t(X) %*% (rowSums(e) - c(exp(X %*% beta)) * c(r %*% lambda)))
    }

    sStep <- function(f, beta, e) {
        if (ncol(X) == 1) {
            beta <- uniroot(f, engine@interval, e=e)$root
        } else {
            beta <- nleqslv(beta, function(x) f(x, e))$x
        }

        lambda <- colSums(e) / c(t(r) %*% exp(X %*% beta))
        list(beta=beta,
             lambda=lambda)
    }

    ##############################
    # Initialize e and r matrix
    e <- r <- matrix(0, N, K)
    for (i in 1:N) {
        set <- which(!is.na(panelMatrix[i, ]))
        mi <- tail(set, 1)
        dset <- diff(c(0, set))

        e[i, 1:mi] <- rep(panelMatrix[i, set] / dset, dset)
        r[i, 1:mi] <- 1
    }

    convergence <- 1
    sRes <- sStep(f, engine@betaInit, e)
    for (i in 2:engine@maxIter) {
        e <- eStep(sRes$lambda)

        betaPre <- sRes$beta
        sRes <- sStep(f, sRes$beta, e)
        s <- sRes$beta - betaPre

        if (max(abs(s)) < engine@absTol | max(abs(s / betaPre)) < engine@relTol) {
            convergence <- 0
            break
        }
    }
    iter <- i

    list(beta=sRes$beta,
         baseline=stepfun(timeGrid, cumsum(c(0, sRes$lambda))),
         lambda=sRes$lambda,
         convergence=convergence,
         iter=iter)
}

##############################################################################
# Extension of Augmented Estimating Equations (AEEX)
##############################################################################
doPanelFit.AEEX <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    eStep <- function(lambda) {
        e <- matrix(0, N, K)

        for (i in 1:N) {
            end <- which(!is.na(panelMatrix[i, ]))
            start <- c(1, head(end, -1) + 1)

            for (j in which(panelMatrix[i, end] > 0)) {
                sq <- seq(start[j], end[j])
                e[i, sq] <- panelMatrix[i, end[j]] * lambda[sq] / sum(lambda[sq])
            }

            if (tail(end, 1) < K) {
                sq <- seq(tail(end, 1) + 1, K)
                e[i, sq] <- (sum(panelMatrix[i, end]) + engine@a) * lambda[sq] /
                    (sum(lambda[-sq]) + engine@a)
            }
        }
        e
    }

    # ncol(X) dimensional nonlinear equation
    f <- function(beta, e) {
        lambda <- c(colSums(e)) / sum(exp(X %*% beta))
        c(t(X) %*% (rowSums(e) - c(exp(X %*% beta)) * sum(lambda)))
    }

    sStep <- function(f, beta, e) {
        if (ncol(X) == 1) {
            beta <- uniroot(f, engine@interval, e=e)$root
        } else {
            beta <- nleqslv(beta, function(x) f(x, e))$x
        }

        lambda <- colSums(e) / sum(exp(X %*% beta))
        list(beta=beta,
             lambda=lambda)
    }

    ##############################
    # Initialize e matrix
    e <- matrix(0, N, K)
    for (i in 1:N) {
        sq <- which(!is.na(panelMatrix[i, ]))
        mi <- tail(sq, 1)
        dsq <- diff(c(0, sq))

        e[i, 1:mi] <- rep(panelMatrix[i, sq] / dsq, dsq)

        if (mi < K) {
            e[i, (mi + 1):K] <- sum(panelMatrix[i, sq]) / mi
        }
    }

    convergence <- 1
    sRes <- sStep(f, engine@betaInit, e)
    for (i in 2:engine@maxIter) {
        e <- eStep(sRes$lambda)

        betaPre <- sRes$beta
        sRes <- sStep(f, sRes$beta, e)
        s <- sRes$beta - betaPre

        if (max(abs(s)) < engine@absTol | max(abs(s / betaPre)) < engine@relTol) {
            convergence <- 0
            break
        }
    }
    iter <- i

    list(beta=sRes$beta,
         baseline=stepfun(timeGrid, cumsum(c(0, sRes$lambda))),
         lambda=sRes$lambda,
         convergence=convergence,
         iter=iter)
}

##############################################################################
# Huang, Wang and Zhang's method, Biometrika 2006 (HWZ)
##############################################################################
doPanelFit.HWZ <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    NPEM <- function(p) {
        e <- matrix(0, N, K)

        for (i in 1:N) {
            end <- which(!is.na(panelMatrix[i, ]))
            start <- c(1, head(end, -1) + 1)

            for (j in which(panelMatrix[i, end] > 0)) {
                sq <- seq(start[j], end[j])
                e[i, sq] <- panelMatrix[i, end[j]] * p[sq] / sum(p[sq])
            }

            yi <- tail(end, 1)
            if (yi < K) {
                sq1 <- seq(1, yi)
                sq2 <- seq(yi + 1, K)
                e[i, sq2] <- sum(panelMatrix[i, end]) / sum(p[sq1]) * p[sq2]
            }
        }

        d <- colSums(e)
        p <- d / sum(d)
    }

    ##############################
    convergence <- 1
    p <- rep(1/K, K)
    for (i in 1:engine@maxIter) {
        pPre <- p
        p <- NPEM(p)
        if (max(abs(pPre - p)) < engine@absTol) {
            convergence <- 0
            break
        }
    }
    iter <- i

    m <- rowSums(panelMatrix, na.rm=TRUE)  # N*1
    y <- apply(panelMatrix, 1, function(x) tail(which(!is.na(x)), 1))
    F <- cumsum(p)
    n <- m / F[y]

    f <- function(beta, w) {
        totalBase <- sum(w * n) / sum(w * c(exp(X %*% beta)))
        c((w * t(X)) %*% (n - totalBase * c(exp(X %*% beta))))
    }

    sStep <- function(f, beta, w) {
        if (ncol(X) == 1) {
            beta <- uniroot(f, engine@interval, w=w)$root
        } else {
            beta <- nleqslv(beta, function(x) f(x, w))$x
        }

        totalBase <- sum(w * n) / sum(w * c(exp(X %*% beta)))
        list(beta=beta,
             totalBase=totalBase)
    }

    w <- rep(1, N)
    if (engine@adjust == "W")
        w[which(F[y] < 0.005)] <- 0
    if (engine@adjust == "H")
        n[which(F[y] < 1e-4)] <- 0

    sRes <- sStep(f, engine@betaInit, w)

    if (!engine@unitWeight) {
        w <- sRes$totalBase * c(exp(X %*% sRes$beta)) /
            mean((n - sRes$totalBase * c(exp(X %*% sRes$beta)))^2)

        if (engine@adjust == "W")
            w[which(F[y] < 0.005)] <- 0
        if (engine@adjust == "H")
            n[which(F[y] < 1e-4)] <- 0

        w <- w / sum(w)
        sRes <- sStep(f, sRes$beta, w)
    }

    list(beta=sRes$beta,
         baseline=stepfun(timeGrid, c(0, F * sRes$totalBase)),
         convergence=convergence,
         iter=iter)
}

##############################################################################
# Zhang's pseudo likelihood method, Biometrika (2002)
##############################################################################
doPanelFit.MPL <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    obsMatrix <- !is.na(panelMatrix)
    NMatrix <- matrix(0, N, K)
    for (i in 1:N) {
        sq <- !is.na(panelMatrix[i, ])
        NMatrix[i, sq] <- cumsum(panelMatrix[i, sq])
    }

    wNbar <- colSums(NMatrix)

    # Cumulative sum upper triangle matrix, [i, j] element = sum_{s=i}^j(x_s)
    cumSumMatrix <- function(x) {
        m <- outer(rep(1, length(x)), x)
        m[row(m) > col(m)] <- 0
        t(apply(m, 1, cumsum))
    }

    cumSumMatrix.wNbar <- cumSumMatrix(wNbar)

    # Isotonic regression estimates of Lambda, use max-min formula
    # Same as pava(wNbar / wAbar, wAbar) in 'Iso' package
    iso <- function(beta) {
        wAbar <- colSums(obsMatrix * c(exp(X %*% beta)))
        m <- cumSumMatrix.wNbar / cumSumMatrix(wAbar)
        m[row(m) > col(m)] <- 0

        for (i in 1:(K - 1))
            m[, K - i] <- pmin(m[, K - i], m[, K - i + 1])

        apply(m, 2, max)
    }

    f <- function(beta, Lambda) {
        c(t(X) %*% rowSums(NMatrix - outer(c(exp(X %*% beta)), Lambda) * obsMatrix))
    }

    beta <- engine@betaInit
    convergence <- 1
    for (i in 1:engine@maxIter) {
        Lambda <- iso(beta)
        betaPre <- beta

        if (ncol(X) == 1) {
            beta <- uniroot(f, engine@interval, Lambda)$root
        } else {
            beta <- nleqslv(beta, function(x) f(x, Lambda))$x
        }

        s <- beta - betaPre
        if (max(abs(s)) < engine@absTol | max(abs(s / betaPre)) < engine@relTol) {
            convergence <- 0
            break
        }
    }
    iter <- i

    list(beta=beta,
         baseline=stepfun(timeGrid, c(0, Lambda)),
         convergence=convergence,
         iter=iter)
}

##############################################################################
# Sun and Wei's method (SW), independent observation and censoring times
##############################################################################
doPanelFit.SW <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)
    X <- X - outer(rep(1, nrow(X)), colMeans(X))

    Nbar <- apply(panelMatrix, 1, function(x) sum(cumsum(x[!is.na(x)])))
    f <- function(beta) {
        c(t(X) %*% c(exp(- X %*% beta) * Nbar))
    }

    if (ncol(X) == 1) {
        beta <- uniroot(f, engine@interval)$root
    } else {
        beta <- nleqslv(engine@betaInit, f)$x
    }

    list(beta=beta,
         baseline=function(x) rep(NA, length(x)),
         convergence=0)
}

##############################################################################
# Bootstrap variance estimation for any engine
##############################################################################
doPanelFit.Engine.Bootstrap <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)
    res <- doPanelFit(panelMatrix, timeGrid, X, engine, NULL)
    engine@betaInit <- res$beta

    R <- stdErr@R
    betaMatrix <- matrix(0, R, length(res$beta))
    baselineMatrix <- matrix(0, R, K)
    convergence <- rep(0, R)

    for (i in 1:R) {
        subjects <- sort(sample(1:N, size=N, replace=TRUE))
        panelMatrix2 <- panelMatrix[subjects, ]
        X2 <- as.matrix(X[subjects, ])
        subCol <- which(colSums(!is.na(panelMatrix2)) > 0)
        panelMatrix2 <- panelMatrix2[, subCol]
        timeGrid2 <- timeGrid[subCol]

        res2 <- doPanelFit(panelMatrix2, timeGrid2, X2, engine, NULL)
        betaMatrix[i, ] <- res2$beta
        baselineMatrix[i, ] <- res2$baseline(timeGrid)
        convergence[i] <- res2$convergence
    }

    converged <- which(convergence == 0)
    betaVar <- var(betaMatrix[converged, ])
    betaSE <- sqrt(diag(as.matrix(betaVar)))
    baselineSE <- stepfun(timeGrid, c(0, sd(baselineMatrix[converged, ])))

    c(res, list(betaSE=betaSE, betaVar=betaVar,
                baselineSE=baselineSE, R=length(converged)))
}

##############################################################################
# Imputation based variance estimation for AEE only
##############################################################################
doPanelFit.AEE.Impute <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    eTime <- timeGrid
    sTime <- c(0, head(eTime, -1))

    # Impute the exact recurrent event time, use multinomial assumption
    imputeSurvData <- function(lambda) {
        survData <- NULL

        for (i in 1:N) {
            eIndex <- as.numeric(which(!is.na(panelMatrix[i,])))
            sIndex <- c(0, head(eIndex, -1)) + 1
            impTime <- event <- NULL

            for (j in 1:length(eIndex)) {
                if (panelMatrix[i, eIndex[j]] > 0) {
                    indexSeq <- seq(sIndex[j], eIndex[j])
                    numEvent <- c(rmultinom(1,
                                            size=panelMatrix[i, eIndex[j]],
                                            prob=lambda[indexSeq]))

                    for (k in which(numEvent > 0)) {
                        impTime <- c(impTime, sort(runif(numEvent[k],
                                                         sTime[sIndex[j] - 1 + k],
                                                         eTime[sIndex[j] - 1 + k])))
                    }
                    event <- c(event, rep(1, panelMatrix[i, eIndex[j]]))
                } else {
                    impTime <- c(impTime, eTime[eIndex[j]])
                    event <- c(event, 0)
                }
            }
            survData <- rbind(survData, data.frame(ID=i,
                                                   start=c(0, head(impTime, -1)),
                                                   end=impTime,
                                                   event=event))
        }
        survData
    }

    ##############################
    R <- stdErr@R
    betaMatrix <- matrix(0, R, ncol(X))
    betaVarMatrix <- matrix(0, ncol(X), ncol(X))

    baselineMatrix <- matrix(0, R, K)
    baselineVarMatrix <- matrix(0, R, K)

    res <- doPanelFit(panelMatrix, timeGrid, X, engine, NULL)

    for (i in 1:R) {
        survData <- imputeSurvData(res$lambda)
        freq <- as.numeric(table(survData$ID))

        survData <- cbind(survData, apply(X, 2, function(y) rep(y, freq)))
        attr(survData, "names")[-(1:4)] <- xname <- paste("X", 1:ncol(X), sep="")

        fml <- as.formula(paste("Surv(start, end, event) ~ ",
                                paste(xname, collapse="+"),
                                "+ cluster(ID)"))

        resCoxph <- coxph(fml, survData)

        # beta
        betaMatrix[i, ] <- as.numeric(resCoxph$coefficients)
        betaVarMatrix <- betaVarMatrix + resCoxph$var

        # baseline
        newdf <- data.frame(t(rep(0, ncol(as.matrix(X)))))
        attr(newdf, "names") <- xname
        sfit <- survfit(resCoxph, newdata=newdf)

        baselineMatrix[i, ] <- stepfun(sfit$time, c(0, -log(sfit$surv)))(timeGrid)
        baselineVarMatrix[i, ] <- (stepfun(sfit$time, c(0, sfit$std.err))(timeGrid))^2
    }

    betaVar <- betaVarMatrix / R   + (1 + 1/R) * var(betaMatrix)
    betaSE <- sqrt(diag(as.matrix(betaVar)))
    baselineSE <- stepfun(timeGrid,
                          c(0, sqrt(colMeans(baselineVarMatrix) +
                                    (1 + 1/R) * sd(baselineMatrix)^2)))

    c(res, list(betaSE=betaSE, betaVar=betaVar, baselineSE=baselineSE))
}

##############################################################################
# Imputation based variance estimation for AEEX
##############################################################################
doPanelFit.AEEX.Impute <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)
    res <- doPanelFit(panelMatrix, timeGrid, X, engine, NULL)

    # Impute the number of count between the actuall censoring time and tau
    for (i in 1:N) {
        end <- which(!is.na(panelMatrix[i, ]))
        if (tail(end, 1) < K) {
            sq <- seq(tail(end, 1) + 1, K)
            panelMatrix[i, K] <- (sum(panelMatrix[i, end]) + engine@a) *
                sum(res$lambda[sq]) / (sum(res$lambda[-sq]) + engine@a)
        }
    }

    # Use doPanelFit.AEE.Impute
    seRes <- doPanelFit(panelMatrix, timeGrid, X,
                        new("AEE", betaInit=engine@betaInit), stdErr)

    c(res, list(betaSE=seRes$betaSE, betaVar=seRes$betaVar, baselineSE=seRes$baselineSE))
}

##############################################################################
# Observed information matrix based variance estimation for AEE
##############################################################################
doPanelFit.AEE.Sandwich <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    res <- doPanelFit(panelMatrix, timeGrid, X, engine, NULL)
    beta <- res$beta
    lambda <- res$lambda

    atRiskMatrix <- matrix(0, N, K)
    lastObs <- apply(panelMatrix, 1, function(x) tail(which(!is.na(x)), 1))
    atRiskMatrix[col(atRiskMatrix) <= lastObs] <- 1

    # A is the complete information matrix for all subjects
    A11 <- diag(c(t(exp(X %*% beta)) %*% atRiskMatrix))
    A21 <- t(c(exp(X %*% beta)) * X) %*% atRiskMatrix
    A22 <- t(X) %*% (X * c(exp(X %*% beta)) * c(atRiskMatrix %*% lambda))
    A <- rbind(cbind(A11, t(A21) * lambda),
               cbind(A21, A22))

    # B is the missing information matrix for all subjects
    B <- matrix(0, K + ncol(X), K + ncol(X))
    for (i in 1:N) {
        sq <- which(!is.na(panelMatrix[i, ]))
        mi <- panelMatrix[i, sq]

        if (is.na(panelMatrix[i, K])) {
            sq <- c(sq, K)
            mi <- c(mi, 0)
        }

        dsq <- diff(c(0, sq))
        ndsq <- length(dsq)

        # normalize lambda, multinomial
        p <- lambda / rep(diff(c(0, cumsum(lambda)[sq])), dsq)
        p[which(p == Inf)] <- 1
        blkp <- p * diag(1, ndsq)[rep(1:ndsq, dsq), ]

        # atRisk (rij) is taken care of
        Xi <- X[i, ]
        B11 <- rep(mi, dsq) * (diag(p) - blkp %*% t(blkp))
        B12 <- rowSums(B11) %*% t(Xi)
        B22 <- sum(B11) * Xi %*% t(Xi)
        B <- B + rbind(cbind(B11, B12),
                       cbind(t(B12), B22))
        # matrix.plot(B)
    }

    # Inverse of observed information matrix
    V <- solve(A - B)

    # Regularization
    dgV <- diag(V)
    dgV[which(dgV < 0)] <- 0
    diag(V) <- dgV

    # Sandwich estimator
    e <- matrix(0, N, K)
    for (i in 1:N) {
        end <- which(!is.na(panelMatrix[i, ]))
        start <- c(1, head(end, -1) + 1)

        for (j in which(panelMatrix[i, end] > 0)) {
            sq <- seq(start[j], end[j])
            e[i, sq] <- panelMatrix[i, end[j]] * lambda[sq] / sum(lambda[sq])
        }
    }
    U1 <- (t(e) - outer(lambda, c(exp(X %*% beta)))) * t(atRiskMatrix)
    U2 <- t(colSums(U1) * X)
    U <- rbind(U1, U2)
    V <- V %*% (U %*% t(U)) %*% t(V)

    ##
    betaVar <- V[-c(1:K), -c(1:K)]
    betaSE <- sqrt(diag(as.matrix(betaVar)))

    lowOne <- matrix(0, K, K)
    lowOne[row(lowOne) >= col(lowOne)] <- 1
    vLambda <- diag(lowOne %*% V[1:K, 1:K] %*% t(lowOne))
    baselineSE <- stepfun(timeGrid, c(0, sqrt(vLambda)))

    c(res, list(betaSE=betaSE, betaVar=betaVar, baselineSE=baselineSE))
}

##############################################################################
# Observed information matrix based variance estimation for AEEX
##############################################################################
doPanelFit.AEEX.Sandwich <- function(panelMatrix, timeGrid, X, engine, stdErr) {
    N <- nrow(panelMatrix)
    K <- ncol(panelMatrix)

    res <- doPanelFit(panelMatrix, timeGrid, X, engine, NULL)
    beta <- res$beta
    lambda <- res$lambda

    atRiskMatrix <- matrix(1, N, K)

    # A is the complete information matrix for all subjects
    A11 <- diag(c(t(exp(X %*% beta)) %*% atRiskMatrix))
    A21 <- t(c(exp(X %*% beta)) * X) %*% atRiskMatrix
    A22 <- t(X) %*% (X * c(exp(X %*% beta)) * c(atRiskMatrix %*% lambda))
    A <- rbind(cbind(A11, t(A21) * lambda),
               cbind(A21, A22))

    # B is the missing information matrix for all subjects
    B <- matrix(0, K + ncol(X), K + ncol(X))
    for (i in 1:N) {
        sq <- which(!is.na(panelMatrix[i, ]))
        mi <- panelMatrix[i, sq]

        if (is.na(panelMatrix[i, K])) {
            y1 <- sum(mi)
            mu1 <- sum(lambda[seq(1, tail(mi, 1))])
            mu2 <- sum(lambda[seq(tail(mi, 1) + 1, K)])

            sq <- c(sq, K)
            mi <- c(mi, (y1 + engine@a) * mu2 / (mu1 + engine@a))
        }

        dsq <- diff(c(0, sq))
        ndsq <- length(dsq)

        # normalize lambda, multinomial
        p <- lambda / rep(diff(c(0, cumsum(lambda)[sq])), dsq)
        p[which(p == Inf)] <- 1
        blkp <- p * diag(1, ndsq)[rep(1:ndsq, dsq), ]

        Xi <- X[i, ]
        B11 <- rep(mi, dsq) * (diag(p) - blkp %*% t(blkp))

        if (is.na(panelMatrix[i, K])) {
            p[seq(1, K - tail(dsq, 1))] <- 0
            B11 <- B11 + outer(p, p) * (y1 + engine@a) *
                mu2 / (mu2 + engine@a) * (1 + mu1 / (mu2 + engine@a))
        }

        B12 <- rowSums(B11) %*% t(Xi)
        B22 <- sum(B11) * Xi %*% t(Xi)
        B <- B + rbind(cbind(B11, B12),
                       cbind(t(B12), B22))
        # matrix.plot(B)
    }

    # Inverse of observed information matrix
    V <- solve(A - B)

    # Regularization
    dgV <- diag(V)
    dgV[which(dgV < 0)] <- 0
    diag(V) <- dgV

    # Sandwich estimator
    e <- matrix(0, N, K)
    for (i in 1:N) {
        end <- which(!is.na(panelMatrix[i, ]))
        start <- c(1, head(end, -1) + 1)

        for (j in which(panelMatrix[i, end] > 0)) {
            sq <- seq(start[j], end[j])
            e[i, sq] <- panelMatrix[i, end[j]] * lambda[sq] / sum(lambda[sq])
        }
    }
    U1 <- (t(e) - outer(lambda, c(exp(X %*% beta)))) * t(atRiskMatrix)
    U2 <- t(colSums(U1) * X)
    U <- rbind(U1, U2)
    V <- V %*% (U %*% t(U)) %*% t(V)

    ##
    betaVar <- V[-c(1:K), -c(1:K)]
    betaSE <- sqrt(diag(as.matrix(betaVar)))

    lowOne <- matrix(0, K, K)
    lowOne[row(lowOne) >= col(lowOne)] <- 1
    vLambda <- diag(lowOne %*% V[1:K, 1:K] %*% t(lowOne))
    baselineSE <- stepfun(timeGrid, c(0, sqrt(vLambda)))

    c(res, list(betaSE=betaSE, betaVar=betaVar, baselineSE=baselineSE))
}

##############################################################################
# Class Definition
##############################################################################
setClass("Engine",
         representation(betaInit="numeric", interval="numeric",
                        maxIter="numeric", absTol="numeric", relTol="numeric"),
         prototype(betaInit=0, interval=c(-5, 5),
                   maxIter=150, absTol=1e-6, relTol=1e-6),
         contains="VIRTUAL")

setClass("AEE",
         representation(),
         prototype(),
         contains="Engine")

setClass("AEEX",
         representation(a="numeric"),
         prototype(maxIter=500, a=0.1),
         contains="Engine")

setClass("HWZ",
         representation(unitWeight="logical", adjust="character"),
         prototype(maxIter=500, absTol=1e-4,
                   unitWeight=TRUE, adjust="W"),
         validity=function(object) {
             object@adjust <- match.arg(object@adjust,
                                        choices=c("W", "H"))
             return(TRUE)
         },
         contains="Engine")

setClass("MPL",
         contains="Engine")

setClass("SW",
         contains="Engine")

setClass("StdErr")

setClass("Bootstrap",
         representation(R="numeric"),
         prototype(R=30),
         contains="StdErr")

setClass("Impute",
         representation(R="numeric"),
         prototype(R=30),
         contains="StdErr")

setClass("Sandwich",
         representation(),
         prototype(),
         contains="StdErr")

##############################################################################
# Method Dispatch
##############################################################################
setGeneric("doPanelFit",
           function(panelMatrix, timeGrid, X, engine, stdErr) {
               standardGeneric("doPanelFit")
           })

setMethod("doPanelFit",
          signature(engine="AEE", stdErr="NULL"),
          doPanelFit.AEE)

setMethod("doPanelFit",
          signature(engine="AEEX", stdErr="NULL"),
          doPanelFit.AEEX)

setMethod("doPanelFit",
          signature(engine="HWZ", stdErr="NULL"),
          doPanelFit.HWZ)

setMethod("doPanelFit",
          signature(engine="MPL", stdErr="NULL"),
          doPanelFit.MPL)

setMethod("doPanelFit",
          signature(engine="SW", stdErr="NULL"),
          doPanelFit.SW)

setMethod("doPanelFit",
          signature(engine="Engine", stdErr="Bootstrap"),
          doPanelFit.Engine.Bootstrap)

setMethod("doPanelFit",
          signature(engine="AEE", stdErr="Impute"),
          doPanelFit.AEE.Impute)

setMethod("doPanelFit",
          signature(engine="AEE", stdErr="Sandwich"),
          doPanelFit.AEE.Sandwich)

setMethod("doPanelFit",
          signature(engine="AEEX", stdErr="Impute"),
          doPanelFit.AEEX.Impute)

setMethod("doPanelFit",
          signature(engine="AEEX", stdErr="Sandwich"),
          doPanelFit.AEEX.Sandwich)

##############################################################################
# Print a panelReg object
##############################################################################
print.panelReg <- function(x, digits=max(options()$digits - 4, 3), ...) {
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$beta
    if (is.null(x$betaSE))
        se <- rep(NA, length(coef))
    else
        se <- x$betaSE

    tmp <- data.frame(coef, exp(coef), se, coef/se,
                      signif(1 - pchisq((coef/ se)^2, 1), digits -1))

    dimnames(tmp) <- list(names(coef),
                          c("coef", "exp(coef)", "se(coef)", "z", "p"))

    # Print results
    cat("\n")
    cat("Call:\n")
    dput(x$call)
    cat("\n")
    prmatrix(tmp)
    cat("\n")
    invisible()
}

##############################################################################
# Plot a PanelReg object
##############################################################################
plot.panelReg <- function(x, se=TRUE, ...) {
    time <- knots(x$baseline)
    y <- x$baseline(time)
    yse <- x$baselineSE(time)
    low <- stepfun(time, c(0, y - 1.96 * yse))
    high <- stepfun(time, c(0, y + 1.96 * yse))
    if (se) {
        plot(high, do.points=FALSE, lty=2, xlab="Time",
             ylab=expression(hat(Lambda[0])(t)), main="Cumulative Baseline Mean")
        plot(x$baseline, do.points=FALSE, add=TRUE)
        plot(low, do.points=FALSE, lty=2, add=TRUE)
    }
    else
        plot(fit$baseline, do.points=FALSE, xlab="Time",
             ylab=expression(hat(Lambda[0])(t)), main="Cumulative Baseline Mean")
}

##############################################################################
# User's Main Function
##############################################################################
## panelReg <- function(formula, data, engine, stdErr=NULL) {

##     Call <- match.call()

##     # A PanelSurv object
##     response <- eval(formula[[2]], data)
##     if (!is.PanelSurv(response))
##         stop("Response must be a PanelSurv object")

##     # Get design matrix, remove intercept column and abundant rows
##     id <- eval(formula[[2]][[2]], data)
##     formula[[2]] <- NULL
##     X <- as.matrix(ddply(cbind(id, model.matrix(formula, data)[, -1]),
##                          "id", head, n=1)[, -1])

##     if (length(engine@betaInit) == 1 & ncol(X) > 1)
##         engine@betaInit <- rep(engine@betaInit, ncol(X))
##     if (length(engine@betaInit) > 1 & length(engine@betaInit) != ncol(X))
##         stop("Invalid length of initial beta values!")

##     fit <- doPanelFit(panelMatrix=response$panelMatrix, timeGrid=response$timeGrid,
##                    X=X, engine=engine, stdErr=stdErr)

##     names(fit$beta) <- attr(terms(formula), "term.labels")
##     fit$call <- Call

##     class(fit) <- "panelReg"
##     fit
## }

panelReg <- function(formula, data, method=c("AEE", "AEEX", "HWZ", "MPL", "SW"),
                     se=c("NULL", "Bootstrap", "Impute", "Sandwich"), control=list()) {

    method <- match.arg(method)
    se <- match.arg(se)
    Call <- match.call()

    # A PanelSurv object
    response <- eval(formula[[2]], data)
    if (!is.PanelSurv(response))
        stop("Response must be a PanelSurv object")

    # Get design matrix, remove intercept column and abundant rows
    id <- eval(formula[[2]][[2]], data)
    formula[[2]] <- NULL
    X <- as.matrix(ddply(cbind(id, model.matrix(formula, data)[, -1]),
                         "id", head, n=1)[, -1])

    # Create an Engine object
    engine.control <- control[names(control) %in% names(attr(getClass(method), "slots"))]
    engine <- do.call("new", c(list(Class=method), engine.control))

    if (length(engine@betaInit) == 1 & ncol(X) > 1)
        engine@betaInit <- rep(engine@betaInit, ncol(X))
    if (length(engine@betaInit) > 1 & length(engine@betaInit) != ncol(X))
        stop("Invalid length of initial beta values!")

    # Create a StdErr object
    if (se == "NULL")
        stdErr <- NULL
    else {
        stdErr.control <- control[names(control) %in% names(attr(getClass(se), "slots"))]
        stdErr <- do.call("new", c(list(Class=se), stdErr.control))
    }

    # Dispatch
    fit <- doPanelFit(panelMatrix=response$panelMatrix, timeGrid=response$timeGrid,
                      X=X, engine=engine, stdErr=stdErr)

    names(fit$beta) <- attr(terms(formula), "term.labels")
    fit$call <- Call

    class(fit) <- "panelReg"
    fit
}

