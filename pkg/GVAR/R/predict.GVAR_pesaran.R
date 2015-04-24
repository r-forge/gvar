predict.GVAR2 <- function (x, steps, start = NULL, exo.var = NULL, exo.type = "exo"){
    bigT <- x$bigT
    max.lag <- length(x$H)
    X <- x$X
    d.dat <- t(x$Data[[length(x$Data)]])
    d <- x$arguments$d
    G_1 <- solve(x$G)
    H.tilde <- list()
    for (i in 1:max.lag) {
        H.tilde[[i]] <- G_1 %*% x$H[[i]]
    }
    Upsilon.0.tilde <- NULL
    Upsilon.tilde <- NULL
#    Upsilon.0.tilde <- G_1 %*% x$Upsilon.0
#    Upsilon.tilde <- list()
#    for (i in 1:length(x$Upsilon)) {
#        Upsilon.tilde[[i]] <- G_1 %*% x$Upsilon[[i]]
#    }
    c.0.tilde <- G_1 %*% x$c.0
    c.1.tilde <- G_1 %*% x$c.1
    if (exo.type == "from.sub") {
        region <- x$subsys[exo.var[1]]
        region.length <- nchar(region)
        regions.cut <- substr(rownames(X), 1, region.length)
        exo.index <- ((1:dim(X)[1])[is.element(regions.cut, region)])[exo.var[-1]]
    }
    if (is.null(start)) {
        loc <- dim(X)[2] + 1
    } else {
        loc <- start
    }
    for (i in 1:steps) {
#        cat(i, ": ", loc, "\n", X[exo.index, loc - 1], "\n")
        new.data <- matrix(0, nrow = dim(X)[1], ncol = 1)
        if (!is.null(c.0.tilde)) {
            new.data <- new.data + c.0.tilde
        }
        t <- loc + 1
        if (!is.null(c.1.tilde)) {
            new.data <- new.data + c.1.tilde * t
        }
        if (exo.type == "from.sub") {
#            new.d <- new.data[exo.index, ]
#            for (j in 1:max.lag) {
#                new.d <- new.d + H.tilde[[j]][exo.index, ] %*%
#                  X[, (loc - j)]
#            }
#            if (is.null(start)) {
#                d.dat <- cbind(d.dat, new.d)
#            }
#            else {
#                if (loc <= dim(X)[2]) {
#                  d.dat[, loc] <- new.d
#                }
#                else {
#                  d.dat <- cbind(d.dat, new.d)
#                }
#            }
        } else {
            d.dat <- cbind(d.dat, exo.var[, i])
        }
        for (j in 1:max.lag) {
            new.data <- new.data + H.tilde[[j]] %*% X[, (loc -
                j)]
        }
        if (!is.null(Upsilon.0.tilde)) {
            new.data <- new.data + Upsilon.0.tilde %*% d.dat[,
                loc]
        }
        if (!is.null(Upsilon.tilde)) {
            for (j in 1:length(x$Upsilon)) {
                new.data <- new.data + Upsilon.tilde[[j]] %*%
                  d.dat[, (loc - j)]
            }
        }
        if (is.null(start)) {
            X <- cbind(X, new.data)
        }
        else {
            if (loc <= dim(X)[2]) {
                X[, loc] <- new.data
            }
            else {
                X <- cbind(X, new.data)
            }
        }
        loc <- loc + 1
    }
    if (is.null(start)) {
        res1 <- X[, (bigT + 1):(bigT + steps)]
        colnames(res1) <- paste("T+", 1:steps, sep = "")
        res2 <- d.dat[, (bigT + 1):(bigT + steps)]
    } else {
        res1 <- X[, start:(start + steps - 1)]
        colnames(res1) <- paste(start - 1, "+", 1:steps, sep = "")
#        res2 <- d.dat[, start:(start + steps - 1)]
        res2 <- NULL
    }
    rownames(res1) <- rownames(X)
    res <- list(X = res1, d = res2)
    return(res)
}
