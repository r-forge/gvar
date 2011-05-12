GVAR=function (Data, tw = NULL, p, q = p, r = NULL, weight, Case,
    exo.var = FALSE, d = NULL, lex = NULL, endo = NULL, ord = NULL,
    we = NULL, method = "max.eigen", caseTest = FALSE, weTest = FALSE)
{
    check.temp <- 0
    if (exo.var)
        check.temp <- 1
    if (dim(weight)[1] != (length(Data) - check.temp)) {
        stop("Data and weight matrix not of the same dimension")
    } else {
        if (any(colnames(weight) != names(Data[1:(length(Data) -
            check.temp)])))
            stop("Data and weight matrix are not identically ordered")
    }
    cmodel <- list()
    N <- length(Data) - 1
    dims <- vector()
    for (i in 1:(N + 1)) {
        if (!is.null(dim(Data[[i]]))) {
            dims[i] <- dim(Data[[i]])[1]
        } else {
            dims[i] <- length(Data[[i]])
        }
    }
    max.dim <- max(dims)
    tsi <- tsp(Data[[((1:length(dims))[dims == max(dims)])[1]]])
    if (is.null(tw)) {
        start.ts <- tsi[1]
        end.ts <- tsi[2]
    } else {
        start.ts <- tw[1]
        end.ts <- tw[2]
    }
    freq <- tsi[3]
    dt <- 1/freq
    n.exo <- 0
    ex <- 0
    n.ex <- rep(0, N + 1)
    if (exo.var) {
        N <- N - 1
        ex <- 1
        n.exo <- dim(Data[[length(Data)]])[2]
        if (is.null(n.exo)) {
            n.exo <- 1
        }
        if (is.null(d)) {
            d <- vector("list", N + 1)
            for (i in 1:(N + 1)) {
                d[[i]] <- 1:n.exo
                n.ex[i] <- n.exo
            }
        } else {
            for (i in 1:(N + 1)) {
                n.ex[i] <- length(d[[i]])
            }
        }
    }
    if ((length(p) == 1) && (N >= 1)) {
        p <- rep(p, (N + 1))
    }
    if ((length(q) == 1) && (N >= 1)) {
        q <- rep(q, (N + 1))
    }
    if ((length(Case) == 1) && (N >= 1)) {
        Case <- rep(Case, (N + 1))
    }
    if (length(lex) == 0) {
        lex <- rep(0, (N + 1))
    } else if ((length(lex) == 1) && (N >= 1)) {
        lex <- rep(lex, (N + 1))
    }
    if (is.null(endo)) {
        endo <- list()
        for (i in 1:(N + 1)) {
            endo[[i]] <- 1:dim(Data[[i]])[2]
        }
    }
    n <- vector()
    for (i in 1:(N + 1)) {
        n[i] <- length(endo[[i]])
    }
    if (is.null(ord)) {
        ord <- list()
        for (i in 1:(N + 1)) {
            ord[[i]] <- 1:dim(Data[[i]])[2]
        }
    }
    if (is.null(we)) {
        we <- list()
        fr <- table(ord)
        for (i in 1:(N + 1)) {
            if (max(fr) == (N + 1)) {
                we[[i]] <- as.numeric(names(fr[fr == max(fr)]))
            } else {
                we[[i]] <- vector()
            }
        }
    }
    vars <- list()
    index <- vector()
    for (i in 1:(N + 1)) {
        vars[[i]] <- vector()
        for (j in 1:(N + 1)) {
            if (j != i) {
                vars[[i]] <- c(vars[[i]], we[[j]])
            }
        }
        ordered <- order(ord[[i]])[is.element(ord[[i]], sort(intersect(ord[[i]],
            unique(c(endo[[i]], vars[[i]])))))]
        ordered <- ordered[is.element(sort(intersect(ord[[i]],
            unique(c(endo[[i]], vars[[i]])))), ord[[i]])]
        ord[[i]] <- sort(intersect(ord[[i]], unique(c(endo[[i]],
            vars[[i]]))))
        Data[[i]] <- Data[[i]][, ordered]
        index <- c(index, ord[[i]])
    }
    m <- vector()
    for (i in 1:(N + 1)) {
        m[i] <- n[i] + length(we[[i]])
    }
    X <- matrix(nrow = 0, ncol = max.dim)
    for (i in 1:(length(Data) - ex)) {
        X <- rbind(X, t(Data[[i]]))
    }
    empty <- (1:dim(X)[1])[apply(X, 1, sum) == 0]
    if (length(empty) > 0)
        X <- X[-empty, ]
    colnames(X) <- (1 - max(p, q)):(max.dim - max(p, q))
    if (!is.null(rownames(X))) {
        if (is.null(names(Data))) {
            names(Data) <- paste("R", 0:N, sep = "")
        }
        nam <- vector()
        idx <- vector()
        for (i in 1:(N + 1)) {
            idx[i] <- length(ord[[i]])
            nam <- c(nam, paste(names(Data)[i], rownames(X)[(sum(idx[1:i]) -
                idx[i] + 1):sum(idx[1:i])], sep = "."))
        }
        rownames(X) <- nam
    } else {
        nam <- vector()
        for (i in 0:N) {
            for (j in 1:dim(Data[[i + 1]])[2]) {
                nam <- c(nam, paste("x", i, ".", j, sep = ""))
            }
        }
        rownames(X) <- nam
    }
    W <- list()
    idx <- vector()
    for (i in 1:(N + 1)) {
        W[[i]] <- matrix(0, nrow = m[i], ncol = dim(X)[1])
        temp <- intersect(ord[[i]], sort(unique(c(endo[[i]],
            vars[[i]]))))
        idx[i] <- length(temp)
        counter <- 0
        pos <- (1:length(temp))[is.element(temp, endo[[i]])]
        for (j in pos) {
            counter <- counter + 1
            if (i == 1) {
                l <- 0
            } else {
                l <- sum(idx[1:i - 1])
            }
            W[[i]][counter, l + j] <- 1
        }
        if (m[i]>n[i])
        { 
            for (j in (n[i] + 1):m[i]) {
              tf <- vector()
              check.where <- vector()
              for (l in 1:(N + 1)) {
                  tf <- c(tf, is.element(we[[i]][j - n[i]], ord[[l]]))
              }
              W[[i]][j, ][index == we[[i]][j - n[i]]] <- as.numeric((weight[i,
                tf]))
            }
        }
        W[[i]] <- W[[i]]/apply(W[[i]], 1, sum)
    }
    names.we <- vector(length = length(unique(unlist(endo))))
    for (i in 1:length(names.we)) {
        for (j in 1:(N + 1)) {
            if (is.element(i, endo[[j]])) {
                temp <- colnames(Data[[j]])[is.element(endo[[j]],i)]
                break
            }
        }
        names.we[i] <- temp
    }
    
    
    # calc models
    
    
    models <- caseList <- list()
    if (is.null(r)) {
        r <- vector()
    }
    rr <- vector()
    for (i in 1:(N + 1)) {
        if (exo.var) {
            if (n.exo == 1 && !is.null(d[[i]])) {
                z <- ts(cbind(t(W[[i]] %*% X), Data[[N + 2]]),
                  start = start.ts, freq = freq)
            } else {
                z <- ts(cbind(t(W[[i]] %*% X), Data[[N + 2]][,
                  d[[i]]]), start = start.ts, freq = freq)
            }
        } else {
            z <- ts(t(W[[i]] %*% X), start = start.ts, freq = freq)
        }
        cols <- c(colnames(Data[[i]])[is.element(ord[[i]], endo[[i]])])
        if (m[i]-n[i]>0) {
            cols <- c(cols, paste(names.we[we[[i]]], "*", sep = ""))
        }
        if (exo.var) {
            if (!is.null(d[[i]])) {
                cols <- c(cols, paste(colnames(Data[[length(Data)]])[d[[i]]],
                  "**", sep = ""))
            }
        }
        if ((n[i] + length(we[[i]]) + length(d[[i]])) == length(cols)) {
            colnames(z) <- cols
        }
        etw <- list(start = start.ts + max(p, q) * dt, end = end.ts,
            freq = freq)
        if (m[i]-n[i]>0) {    
            ranks <- rank.test.we(z.ts = z, etw = etw, p = p[i],
            q = q[i], n = n[i], case = Case[i], ex = n.ex[i],
            lex = lex[i])
        } else {
            ranks <- rank.test.vecm(Y.ts = z,etw,p = p[i],case = Case[i], 
            ex = n.ex[i], lex = lex[i], season=NULL,season.start.time=NULL)
        }
        if (is.na(r[i])) {
            if (method == "trace") {
                for (r_ in 0:(n[i] - 1)) {
                  CV <- ranks$CV.trace[[paste("rank", r_, "vs.",
                    n[i])]][1, 1]
                  if ((ranks$LR.trace[[paste("rank", r_, "vs.",
                    n[i])]][1] < CV) && is.na(r[i])) {
                    r[i] <- r_
                  }
                }
            } else {
                for (r_ in 0:(n[i] - 1)) {
                  CV <- ranks$CV.maxeigen[[paste("rank", r_,
                    "vs.", r_ + 1)]][1, 1]
                  if ((is.na(r[i]) && (ranks$LR.maxeigen[[paste("rank",
                    r_, "vs.", r_ + 1)]][1] < CV))) {
                    r[i] <- r_
                  }
                }
            }
            if (is.na(r[i])) {
                r[i] <- n[i]
            }
        }
        rr[i] <- r[i]
        if (r[i] == 0) {
            rr[i] <- 1
        }
        if (caseTest) {
            cat(names(Data)[i], ":\n")
            if (m[i]-n[i]>0) {
                caseList[[i]] = case.test(z.ts = z, etw = etw, p = p[i],
                  q = q[i], n = n[i], r = rr[i], type = "we.vecm")
            } else {
                caseList[[i]] = case.test(z.ts = z, etw = etw, p = p[i],
                n = n[i], r = rr[i], type = "vecm")
            }    
            cat("\n")
        }
        exo <- NULL
        if (!is.null(d[[i]])) {exo = length(d[[i]])}
        if(n.ex[i]==0) {lex.temp <- NULL} else {lex.temp <- lex[i]}
        if (m[i]-n[i]>0) {
            mdls <- est.we.mdls(z.ts = z, etw = etw, p = p[i], q = q[i], r = rr[i], n = n[i], case = Case[i], ex = n.ex[i], lex = lex.temp, we.test = weTest)
        } else {
            mdls <- est.vecm.mdls(Y.ts = z, etw = etw, p = p[i], r = rr[i], ex = n.ex[i], lex = lex.temp, case = Case[i])
        }
        models[[i]] <- mdls
        names(models)[i] <- names(Data)[i]
        cmodel[[i]] <- set.mdl(mdls, exo = exo)
    }
    G <- NULL
    for (i in 1:(N + 1)) {
        Ac <- cbind(diag(n[i]), -cmodel[[i]]$VAR$B_0)
        Wc <- W[[i]]
        G <- rbind(G, Ac %*% Wc)
    }
    H <- vector("list", max(p, q))
    for (i in 1:max(p, q)) {
        for (j in 1:(N + 1)) {
            if (p[j] < i) {
                AA <- matrix(0, dim(cmodel[[j]]$VAR$A[[1]])[1],
                  dim(cmodel[[j]]$VAR$A[[1]])[2])
            }
            else {
                AA <- cmodel[[j]]$VAR$A[[i]]
            }
            if (q[j] < i) {
                BB <- matrix(0, dim(cmodel[[j]]$VAR$B[[1]])[1],
                  dim(cmodel[[j]]$VAR$B[[1]])[2])
            }
            else {
                BB <- cmodel[[j]]$VAR$B[[i]]
            }
            Ac <- cbind(AA, BB)
            H[[i]] <- rbind(H[[i]], Ac %*% W[[j]])
        }
    }
    c.0 <- NULL
    c.1 <- NULL
    if (exo.var) {
        Upsilon.0 <- NULL
        Upsilon <- vector("list", max(lex))
    }
    for (i in 1:(N + 1)) {
        if (!is.null(cmodel[[i]]$VAR$c.0)) {
            c.0 <- rbind(c.0, as.matrix(cmodel[[i]]$VAR$c.0))
        }
        else c.0 <- rbind(c.0, matrix(0, ncol = 1, nrow = cmodel[[i]]$VAR$n))
        if (!is.null(cmodel[[i]]$VAR$c.1)) {
            c.1 <- rbind(c.1, as.matrix(cmodel[[i]]$VAR$c.1))
        }
        else c.1 <- rbind(c.1, matrix(0, ncol = 1, nrow = cmodel[[i]]$VAR$n))
        if (exo.var) {
            if (is.null(dim(Data[[N + 2]])[2])) {
                Y.0 <- matrix(0, nrow = n[i], ncol = 1)
            }
            else {
                Y.0 <- matrix(0, nrow = n[i], ncol = dim(Data[[N +
                  2]])[2])
            }
            if (!is.null(d[[i]])) {
                Y.0[, d[[i]]] <- as.matrix(cmodel[[i]]$VAR$Upsilon_0)
            }
            Upsilon.0 <- rbind(Upsilon.0, Y.0)
            for (j in 1:max(lex)) {
                if (is.null(dim(Data[[N + 2]])[2])) {
                  Y.x <- matrix(0, nrow = n[i], ncol = 1)
                }
                else {
                  Y.x <- matrix(0, nrow = n[i], ncol = dim(Data[[N +
                    2]])[2])
                }
                if (!is.null(d[[i]])) {
                  if (lex[i] >= j) {
                    Y.x[, 1:d[[i]]] <- as.matrix(cmodel[[i]]$VAR$Upsilon[[j]])[,
                      1:d[[i]]]
                  }
                }
                Upsilon[[j]] <- rbind(Upsilon[[j]], Y.x)
            }
        }
    }
    bigT <- dim(X)[2]
    U <- G %*% X[, (max(p, q) + 1):bigT] - matrix(c.0, nrow = sum(n),
        ncol = bigT - max(p, q)) - c.1 %*% t(1:(bigT - max(p, q)))
    for (k in 1:max(p, q)) {
        U <- U - H[[k]] %*% X[, (max(p, q) + 1 - k):(bigT - k)]
    }
    if (exo.var) {
        if (is.null(dim(Data[[length(Data)]]))) {
            U <- U - Upsilon.0 %*% t(Data[[length(Data)]][(max(p,
                q) + 1):bigT])
            for (k in 1:max(lex)) {
                U <- U - Upsilon[[k]] %*% t(Data[[length(Data)]][(max(p,
                  q) + 1 - k):(bigT - k)])
            }
        }
        else {
            U <- U - Upsilon.0 %*% t(Data[[length(Data)]][(max(p,
                q) + 1):bigT, ])
            for (k in 1:max(lex)) {
                U <- U - Upsilon[[k]] %*% t(Data[[length(Data)]][(max(p,
                  q) + 1 - k):(bigT - k), ])
            }
        }
    }
    U.mean <- apply(U, 1, mean)
    U.cov <- tcrossprod(U)/(bigT - max(p, q))
    w.m <- weight/apply(weight, 1, sum)
    rownames(w.m) <- colnames(w.m) <- names(Data)[1:(N + 1)]
    if (exo.var == FALSE) {
        Upsilon.0 <- Upsilon <- NULL
    }



    if (caseTest & length(models) == length(caseList)) {
        names(caseList) = names(models)
    }

    # this function constructs a list object of the arguments submitted to (in this case) the GVAR function
  constrarglist=function (funobj, framepos = -1)
  {
    namedlist = formals(funobj)
    argnames = names(namedlist)

    for (argn in 1:length(namedlist)) {
       #if (exists(argnames[argn], where = sys.frame(framepos))) {
           hihi = (get(argnames[argn], envir = sys.frame(framepos)))
          if (is.null(hihi)) namedlist[[argn]]="NULL" else namedlist[[argn]]=hihi
       #}
    }
    return(as.vector(namedlist,mode="list"))
  }

    res <- list(subsys = names(Data), Data = Data,we.vecms = models, X = X,bigT=bigT, r = rr,
                Case = Case, W = W, G = G, H = H, Upsilon.0 = Upsilon.0, Upsilon = Upsilon,
                c.0 = c.0, c.1 = c.1, caseTest = caseList,
                weight = w.m, U = U, U.cov = U.cov, arguments=constrarglist(GVAR),
                GVAR.call = match.call(GVAR, sys.call(0)))
    class(res) <- "GVAR"
    return(invisible(res))
}