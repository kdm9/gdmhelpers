library(gdm)
library(foreach)
library(fossil)
library(dplyr)


#' Perform Gower's covariance to euclidean distance transform
#'
#' @param C Square, symmetrical covariance matrix
#'
#' @return D, normalised euclidean distance
#' @export
gower_cov2dist = function(C) {
    D = matrix(0, nrow(C), nrow(C))
    colnames(D) = rownames(D) = rownames(C)
    for (i in seq_len(nrow(C))) {
        for (j in seq_len(nrow(C))) {
            D[i,j] = sqrt(C[i, i] + C[j, j] - 2 * C[i, j])
        }
    }
    D
}


make_sitepair = function(samples, respDist, longLatPreds, minGeoDist=1) {
    resp = cbind(ID=samples, as.data.frame(respDist))
    preds = cbind(ID=samples, longLatPreds)
    sp = formatsitepair(
        resp, 3, predData=preds, siteColumn="ID", XColumn="Longitude",
        YColumn="Latitude"
    )
    geodist = fossil::deg.dist(sp$s1.xCoord, sp$s1.yCoord, sp$s2.xCoord, sp$s2.yCoord)
    sp.filt = sp[geodist > minGeoDist,]
}

run_gdm_with_subset = function(samples, keep, respDist, longLatPreds, minGeoDist, splines=NULL, knots=NULL, ret.sp=FALSE) {
    subsamp = samples[keep]

    sp.filt = make_sitepair(subsamp, respDist[keep, keep], longLatPreds[keep,], minGeoDist)
    mdl = gdm::gdm(sp.filt,geo=T, splines=splines, knots=knots)
    if (ret.sp) {
        return(list(gdm=mdl, sp=sp.filt))
    } else {
        return(mdl)
    }
}

#' Extract splines from a GDM model object into a tidy dataframe
#'
#' @param gdm.mdl GDM model object
#'
#' @return A tidy data frame
#' @export
tidysplines = function(gdm.mdl) {
    spl = gdm::isplineExtract(gdm.mdl)
    xs = tidyr::gather(as.data.frame(spl$x), "predictor", "x")
    ys = tidyr::gather(as.data.frame(spl$y), "predictor", "y")
    data.frame(predictor=xs$predictor, x=xs$x, y=ys$y)
}


#' Runs a site-level jack-knifed GDM analysis
#'
#' @param samples Sample names
#' @param respDist Response distance
#' @param longLatPreds Predictors. Cols must be Longitude, Latitude then other preds
#' @param minGeoDist Minimum geographic distance to allow between sites (filters site-pair table), in kilometers
#' @param nperm Number of permutations
#' @param sampleFrac Keep data at `sampleFrac*nSites` sites
#'
#' @return Data frame of per-predictor splines
#' @export
gdm_jackknife_filtered = function (samples, respDist, longLatPreds, minGeoDist=1, nperm=100, sampleFrac=0.8, splines=NULL, knots=NULL) {
    saf = options(stringsAsFactors=FALSE)
    respDist = as.matrix(as.dist(respDist))
    Nsamp = length(samples)
    longLat=longLatPreds[,1:2]

    # site detection
    geoDist = fossil::earth.dist(longLat)
    sites = stats::cutree(stats::hclust(geoDist), h=minGeoDist)
    uniquesites = unique(sites)

    # Run full gdm
    full = tidysplines(run_gdm_with_subset(samples, 1:Nsamp, respDist, longLatPreds, minGeoDist, splines=splines, knots=knots))
    full = cbind(full, perm="true", is.perm=F)

    # Do jacknives
    withperms = foreach(i=seq_len(nperm), .combine="rbind", .init=full, .errorhandling="remove") %dopar% {
        keep = sites %in% sample(uniquesites, length(uniquesites)*sampleFrac, replace = F)
        spl = tidysplines(run_gdm_with_subset(samples, keep, respDist, longLatPreds, minGeoDist, splines=splines, knots=knots))
        if (is.null(spl)) {
            return(NULL)
        }
        cbind(spl, perm=i, is.perm=T)
    }

    # Here we fit a linear model to predict great-circle distance from our points
    lldist = data.frame(eucl=as.numeric(dist(longLat)),
                        gcd=as.numeric(fossil::earth.dist(longLat)))
    lldmdl = lm(gcd~0+eucl, data=lldist)
    isgeogr = withperms$predictor == "Geographic"
    withperms$x[isgeogr] = predict(lldmdl, data.frame(eucl=withperms$x[isgeogr]))
    suppressWarnings(options(saf))
    withperms
}

#' Cross-validates a GDM model at site level
#'
#' @param samples Sample names
#' @param respDist Response distance
#' @param longLatPreds Predictors. Cols must be Longitude, Latitude then other preds
#' @param minGeoDist Minimum geographic distance to allow between sites (filters site-pair table), in kilometers
#' @param nperm Number of permutations
#' @param sampleFrac Train on data at `sampleFrac*nSites` sites
#'
#' @return Data frame of per-predictor splines
#' @export
gdm_crossval = function (samples, respDist, longLatPreds, minGeoDist=1, nperm=100, sampleFrac=0.8, splines=NULL, knots=NULL) {
    saf = options(stringsAsFactors=FALSE)
    respDist = as.matrix(as.dist(respDist))
    Nsamp = length(samples)
    longLat=longLatPreds[,1:2]

    # site detection
    geoDist = fossil::earth.dist(longLat)
    sites = stats::cutree(stats::hclust(geoDist), h=minGeoDist)
    uniquesites = unique(sites)

    # Do jacknives
    folds = foreach(i=seq_len(nperm), .combine="rbind", .errorhandling="pass") %dopar% {
        keep = sites %in% sample(uniquesites, length(uniquesites)*sampleFrac, replace = F)
        gdmsp = run_gdm_with_subset(samples, keep, respDist, longLatPreds, minGeoDist, splines=splines, knots=knots, ret.sp=T)
        pred = predict(gdmsp$gdm, data=gdmsp$sp)
        real = gdmsp$sp$distance
        r = cor(pred, real)
        p = cor(pred, real, method="spearman")
        data.frame(perm=i, r=r, rho=p)
    }

    suppressWarnings(options(saf))
    invisible(folds)
}


#' Perform forwards selection on GDM model context
#'
#' @param samples Sample names
#' @param respDist Response distance
#' @param longLatPreds Predictors. Cols must be Longitude, Latitude then other preds
#' @param minGeoDist Minimum geographic distance to allow between sites (filters site-pair table), in kilometers
#' @param nperm Number of permutations
#' @param sampleFrac Train on data at `sampleFrac*nSites` sites
#'
#' @return Data frame of per-predictor splines
#' @export
gdm_fwdsel = function (sitepair, explained.tolerance=0.1, max.vars=NULL) {
    vars.in = NULL
    nvar = (ncol(sitepair)-6)/2
    if (nvar < 1) stop("No environmental varaibles")
    if (is.null(max.vars)) max.vars = nvar

    sort.sp = function (sp) {
        if (ncol(sp) <= 6) return(sp)
        geo = colnames(sp)[1:6]
        envs = sort(colnames(sp)[7:ncol(sp)])
        sp[,c(geo, envs)]
    }
    sp = sitepair[,1:6]
    last.model = gdm(sp, geo=T)
    progress = data.frame(i=0, varname="geography", deviance.explained = last.model$explained)
    for (i in seq_len(max.vars)) {
        to.test = setdiff(seq_len(nvar), vars.in)
        model.devs = NULL
        for (v in to.test) {
            j = 6 + v # variable offset
            this.sp = sort.sp(cbind(sp, sitepair[,c(j, j+nvar)]))
            class(this.sp) = c("gdmData", class(this.sp))
            this.model = gdm(this.sp, geo=T)
            added.dev = this.model$explained - last.model$explained
            model.devs = c(model.devs, added.dev)
        }
        best.var = to.test[which.max(model.devs)]
        best.var.col = 6 + best.var
        best.var.name = sub("s1.", "", colnames(sitepair)[best.var.col])
        vars.in = c(vars.in, best.var)
        progress = rbind(progress, data.frame(i=i, varname=best.var.name, deviance.explained= max(model.devs)))
        sp = sort.sp(cbind(sp, sitepair[,c(best.var.col, best.var.col + nvar)]))
        class(sp) = c("gdmData", class(sp))
        last.model = gdm(sp, geo=T)
        cat(paste("Iteration", i, ", best var", best.var.name, "added", round(max(model.devs), 2), "\n"))
        if (max(model.devs) < explained.tolerance) break
    }
    invisible(progress)
}

fixed_varimp = function (spTable, geo, splines = NULL, knots = NULL, fullModelOnly = FALSE,
          nPerm = 50, parallel = FALSE, cores = 2, sampleSites = 1,
          sampleSitePairs = 1, outFile = NULL)
{
    k <- NULL
    if (class(spTable)[1] != "gdmData") {
        warning("spTable class does not include type 'gdmData'. Make sure your data is in site-pair format or the gdm model will not fit.")
    }
    if (!(class(spTable)[1] == "gdmData" | class(spTable)[1] ==
          "matrix" | class(spTable)[1] == "data.frame")) {
        stop("spTable argument needs to be gdmData, a matrix, or a data frame")
    }
    if (ncol(spTable) < 6) {
        stop("Not enough columns in data. (Minimum need: Observed, weights, X0, Y0, X1, Y1)")
    }
    if (nrow(spTable) < 1) {
        stop("Not enough rows in data")
    }
    if (!(geo == TRUE | geo == FALSE)) {
        stop("geo argument must be either TRUE or FALSE")
    }
    if (is.null(splines) == FALSE & class(splines) != "numeric") {
        stop("argument splines needs to be a numeric data type")
    }
    if (is.null(knots) == FALSE & class(knots) != "numeric") {
        stop("argument knots needs to be a numeric data type")
    }
    if (!(fullModelOnly == TRUE | fullModelOnly == FALSE)) {
        stop("fullModelOnly argument must be either TRUE or FALSE")
    }
    if ((is.null(nPerm) == FALSE & is.numeric(nPerm) == FALSE) |
        nPerm < 1) {
        stop("argument nPerm needs to be a positive integer")
    }
    if (!(parallel == TRUE | parallel == FALSE)) {
        stop("parallel argument must be either TRUE or FALSE")
    }
    if (parallel == TRUE & is.null(cores) == TRUE) {
        stop("If parallel==TRUE, the number of cores must be specified")
    }
    if ((is.null(cores) == FALSE & is.numeric(cores) == FALSE) |
        cores < 1) {
        stop("argument cores needs to be a positive integer")
    }
    if (is.numeric(sampleSites) == FALSE | sampleSites < 0 |
        sampleSites > 1) {
        stop("argument sampleSites needs to be a positive number between 0 and 1")
    }
    if (is.numeric(sampleSitePairs) == FALSE | sampleSitePairs <
        0 | sampleSitePairs > 1) {
        stop("argument sampleSitePairs needs to be a positive number between 0 and 1")
    }
    if (sampleSites == 0) {
        stop("a sampleSites value of 0 will remove all sites from the analysis")
    }
    if (sampleSitePairs == 0) {
        stop("a sampleSitePairs value of 0 will remove all sites from the analysis")
    }
    if (is.null(outFile) == FALSE) {
        if (is.character(outFile) == FALSE) {
            stop("argument outFile needs to be a character string of the directory and file name you wish the tables to be written to")
        }
        outFileChar <- nchar(outFile)
        if (substr(outFile, outFileChar - 5, outFileChar) !=
            ".RData") {
            outFile <- paste(outFile, ".RData", sep = "")
        }
        if (length(strsplit(outFile, "/")[[1]]) > 1) {
            splitOutFile <- strsplit(outFile, "/")[[1]][-length(strsplit(outFile,
                                                                         "/")[[1]])]
            dir.create(paste(splitOutFile, collapse = "/"))
        }
        else {
            outFile <- paste("./", outFile, sep = "")
        }
    }
    nPerm <- as.integer(nPerm)
    cores <- as.integer(cores)
    if (sampleSites < 1) {
        spTable <- removeSitesFromSitePair(spTable, sampleSites = sampleSites)
        if (sampleSitePairs < 1) {
            warning("You have selected to randomly remove sites and/or site-pairs.")
        }
    }
    if (sampleSitePairs < 1) {
        numRm <- sample(1:nrow(spTable), round(nrow(spTable) *
                                                   (1 - sampleSitePairs)))
        spTable <- spTable[-c(numRm), ]
    }
    rtmp <- spTable[, 1]
    if (length(rtmp[rtmp < 0]) > 0) {
        stop("Response spTable has negative values. Must be between 0 - 1.")
    }
    if (length(rtmp[rtmp > 1]) > 0) {
        stop("Response spTable has values greater than 1. Must be between 0 - 1.")
    }
    nVars <- (ncol(spTable) - 6)/2
    varNames <- colnames(spTable[c(7:(6 + nVars))])
    varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
    if (geo == TRUE) {
        nVars <- nVars + 1
        varNames <- c("Geographic", varNames)
    }
    sortMatX <- sapply(1:nrow(spTable), function(i, spTab) {
        c(spTab[i, 3], spTab[i, 5])
    }, spTab = spTable)
    sortMatY <- sapply(1:nrow(spTable), function(i, spTab) {
        c(spTab[i, 4], spTab[i, 6])
    }, spTab = spTable)
    sortMatNum <- sapply(1:nrow(spTable), function(i) {
        c(1, 2)
    })
    sortMatRow <- sapply(1:nrow(spTable), function(i) {
        c(i, i)
    })
    fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY),
                         as.vector(sortMatNum), as.vector(sortMatRow), rep(NA,
                                                                           length(sortMatX)))
    siteByCoords <- as.data.frame(unique(fullSortMat[, 1:2]))
    numSites <- nrow(siteByCoords)
    for (i in 1:numSites) {
        fullSortMat[which(fullSortMat[, 1] == siteByCoords[i,
                                                           1] & fullSortMat[, 2] == siteByCoords[i, 2]), 5] <- i
    }
    indexTab <- matrix(NA, nrow(spTable), 2)
    for (iRow in 1:nrow(fullSortMat)) {
        indexTab[fullSortMat[iRow, 4], fullSortMat[iRow, 3]] <- fullSortMat[iRow,
                                                                            5]
    }
    rm(fullSortMat)
    rm(sortMatX)
    rm(sortMatY)
    rm(sortMatNum)
    rm(sortMatRow)
    rm(siteByCoords)
    exBySite <- lapply(1:numSites, function(i, index, tab) {
        rowSites <- which(index[, 1] %in% i)
        if (length(rowSites) < 1) {
            rowSites <- which(index[, 2] %in% i)
        }
        exSiteData <- tab[rowSites[1], ]
        return(exSiteData)
    }, index = indexTab, tab = spTable)
    outSite <- which(!(1:numSites %in% indexTab[, 1]))
    for (i in 1:length(exBySite)) {
        siteRow <- exBySite[[i]]
        if (i %in% outSite) {
            siteRow <- siteRow[grep("s2.", colnames(siteRow))]
            colnames(siteRow) <- sapply(strsplit(colnames(siteRow),
                                                 "s2."), "[[", 2)
        }
        else {
            siteRow <- siteRow[grep("s1.", colnames(siteRow))]
            colnames(siteRow) <- sapply(strsplit(colnames(siteRow),
                                                 "s1."), "[[", 2)
        }
        exBySite[[i]] <- siteRow
    }
    siteData <- do.call("rbind", exBySite)
    modelTestValues <- matrix(NA, 4, nVars, dimnames = list(c("Model deviance",
                                                              "Percent deviance explained", "Model p-value", "Fitted permutations"),
                                                            c("fullModel", paste("fullModel-", seq(1, nVars - 1),
                                                                                 sep = ""))))
    devReductVars <- matrix(NA, nVars, nVars - 1)
    rownames(devReductVars) <- varNames
    colnames(devReductVars) <- c("fullModel", paste("fullModel-",
                                                    seq(1, nVars - 2), sep = ""))
    pValues <- numPermsFit <- devReductVars
    currSitePair <- spTable
    nullGDMFullFit <- 0
    for (v in 1:nVars) {
        fullGDM <- gdm(currSitePair, geo = geo, splines = splines,
                       knots = knots)
        if (is.null(fullGDM) == TRUE) {
            warning(paste("The model did not fit when testing variable number ",
                          v, ". Terminating analysis. Returning output objects filled up to the point of termination.",
                          sep = ""))
            break
        }
        if (parallel == TRUE) {
            permSitePairs <- foreach(k = 1:nPerm, .verbose = F, .packages = c("gdm"),
				     .export = c("permutateSitePair", "currSitePair", "siteData", "indexTab", "varNames")) %dopar%
                permutateSitePair(currSitePair, siteData, indexTab,
                                  varNames)
            permGDM <- try(foreach(k = 1:length(permSitePairs),
                                   .verbose = F, .packages = c("gdm")) %dopar%
                               gdm(permSitePairs[[k]], geo = geo, splines = NULL,
                                   knots = NULL))
        }
        else {
            permSitePairs <- lapply(1:nPerm, function(i, csp,
                                                      sd, it, vn) {
                permutateSitePair(csp, sd, it, vn)
            }, csp = currSitePair, sd = siteData, it = indexTab,
            vn = varNames)
            permGDM <- lapply(permSitePairs, gdm, geo = geo,
                              splines = NULL, knots = NULL)
        }
        permModelDev <- sapply(permGDM, function(mod) {
            mod$gdmdeviance
        })
        modPerms <- length(which(sapply(permModelDev, is.null) ==
                                     TRUE))
        if (modPerms > 0) {
            permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,
                                                               is.null) == T))])
        }
        modelTestValues[1, v] <- fullGDM$gdmdeviance
        modelTestValues[2, v] <- fullGDM$explained
        modelTestValues[3, v] <- sum(permModelDev <= fullGDM$gdmdeviance)/(nPerm -
                                                                               modPerms)
        modelTestValues[4, v] <- nPerm - modPerms
        if (length(varNames) < 2) {
            break
        }
        if (geo == TRUE) {
            noGeoGDM <- gdm(currSitePair, geo = FALSE, splines = NULL,
                            knots = NULL)
            if (parallel == TRUE) {
                permSitePairs <- foreach(k = 1:nPerm, .verbose = F, .packages = c("gdm"),
				         .export = c("permutateSitePair", "currSitePair", "siteData", "indexTab", "varNames")) %dopar%
                    permutateSitePair(currSitePair, siteData,
                                      indexTab, varNames)
                permGDM <- try(foreach(k = 1:length(permSitePairs),
                                       .verbose = F, .packages = c("gdm")) %dopar%
                                   gdm(permSitePairs[[k]], geo = geo, splines = NULL,
                                       knots = NULL))
            }
            else {
                permSitePairs <- lapply(1:nPerm, function(i,
                                                          csp, sd, it, vn) {
                    permutateSitePair(csp, sd, it, vn)
                }, csp = currSitePair, sd = siteData, it = indexTab,
                vn = varNames)
                permGDM <- lapply(permSitePairs, gdm, geo = geo,
                                  splines = NULL, knots = NULL)
            }
            permModelDev <- sapply(permGDM, function(mod) {
                mod$gdmdeviance
            })
            modPerms <- length(which(sapply(permModelDev, is.null) ==
                                         TRUE))
            if (modPerms > 0) {
                permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,
                                                                   is.null) == T))])
            }
            if (is.null(noGeoGDM$gdmdeviance) == TRUE) {
                permDevReduct <- -9999
                devReductVars[1, v] <- -9999
                pValues[1, v] <- -9999
            }
            else {
                permDevReduct <- noGeoGDM$gdmdeviance - permModelDev
                devReductVars[1, v] <- 100 * abs((noGeoGDM$explained -
                                                      fullGDM$explained)/fullGDM$explained)
                pValues[1, v] <- sum(permDevReduct >= (noGeoGDM$gdmdeviance -
                                                           fullGDM$gdmdeviance))/(nPerm - modPerms)
            }
            numPermsFit[1, v] <- nPerm - modPerms
        }
        for (varChar in varNames) {
            if (varChar != "Geographic") {
                testVarCols1 <- grep(paste("^s1.", varChar,
                                           "$", sep = ""), colnames(currSitePair))
                testVarCols2 <- grep(paste("^s2.", varChar,
                                           "$", sep = ""), colnames(currSitePair))
                testSitePair <- currSitePair[, -c(testVarCols1,
                                                  testVarCols2)]
                noVarGDM <- gdm(testSitePair, geo = geo, splines = NULL,
                                knots = NULL)
                if (parallel == TRUE) {
                    noVarSitePairs <- foreach(k = 1:nPerm, .verbose = F, .packages = c("gdm"),
				              .export = c("permutateSitePair", "currSitePair",
							  "siteData", "indexTab", "varChar")) %dopar%
                        permutateVarSitePair(currSitePair, siteData,
                                             indexTab, varChar)
                    permGDM <- try(foreach(k = 1:length(noVarSitePairs),
                                           .verbose = F, .packages = c("gdm")) %dopar%
                                       gdm(noVarSitePairs[[k]], geo = geo, splines = NULL,
                                           knots = NULL))
                }
                else {
                    noVarSitePairs <- lapply(1:nPerm, function(i,
                                                               csp, sd, it, vn) {
                        permutateVarSitePair(csp, sd, it, vn)
                    }, csp = currSitePair, sd = siteData, it = indexTab,
                    vn = varChar)
                    permGDM <- lapply(noVarSitePairs, gdm, geo = geo,
                                      splines = NULL, knots = NULL)
                }
                permModelDev <- sapply(permGDM, function(mod) {
                    mod$gdmdeviance
                })
                modPerms <- length(which(sapply(permModelDev,
                                                is.null) == TRUE))
                if (modPerms > 0) {
                    permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,
                                                                       is.null) == T))])
                }
                if (is.null(noVarGDM$gdmdeviance) == TRUE) {
                    permDevReduct <- -9999
                    ggg <- which(rownames(devReductVars) %in%
                                     varChar)
                    devReductVars[ggg, v] <- rep(-9999, times = length(ggg))
                    pValues[ggg, v] <- rep(-9999, times = length(ggg))
                }
                else {
                    permDevReduct <- noVarGDM$gdmdeviance - permModelDev
                    devReductVars[which(rownames(devReductVars) %in%
                                            varChar), v] <- 100 * abs((noVarGDM$explained -
                                                                           fullGDM$explained)/fullGDM$explained)
                    pValues[which(rownames(pValues) %in% varChar),
                            v] <- sum(permDevReduct >= (noVarGDM$gdmdeviance -
                                                            fullGDM$gdmdeviance))/(nPerm - modPerms)
                }
                numPermsFit[which(rownames(numPermsFit) %in%
                                      varChar), v] <- nPerm - modPerms
            }
        }
        if (fullModelOnly == TRUE) {
            break
        }
        tempPVals <- as.numeric(pValues[c(1:nVars), v])
        tempDevs <- as.numeric(devReductVars[c(1:nVars), v])
        tempPVals <- tempPVals[!is.na(tempPVals)]
        tempDevs <- tempDevs[!is.na(tempDevs)]
        varToOmit <- which.max(tempPVals)
        for (iCheck in 1:length(varNames)) {
            if (tempPVals[iCheck] == tempPVals[varToOmit]) {
                if (tempDevs[iCheck] < tempDevs[varToOmit]) {
                    varToOmit <- iCheck
                }
            }
        }
        if (varToOmit == 1 & geo == TRUE) {
            geo <- FALSE
            varNames <- varNames[-1]
        }
        else {
            nameToRemove <- varNames[varToOmit]
            varNames <- varNames[-varToOmit]
            removeFromSitePs1 <- grep(paste("^s1.", nameToRemove,
                                            "$", sep = ""), colnames(currSitePair))
            removeFromSitePs2 <- grep(paste("^s2.", nameToRemove,
                                            "$", sep = ""), colnames(currSitePair))
            currSitePair <- currSitePair[, -c(removeFromSitePs1,
                                              removeFromSitePs2)]
        }
    }
    outObject <- list(modelTestValues, devReductVars, pValues,
                      numPermsFit)
    if (is.null(outFile) == FALSE) {
        save(outObject, file = outFile)
    }
    return(outObject)
}
