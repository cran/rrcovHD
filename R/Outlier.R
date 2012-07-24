setMethod("getClassLabels", "Outlier", function(obj, i=1)
{
        return(which(obj@grp == levels(obj@grp)[i]))
})

setMethod("getDistance", "Outlier", function(obj)
{
   return(NULL)
})

setMethod("getFlag", "Outlier", function(obj, prob=0.975)
{
   return(obj@flag)
})

setMethod("getWeight", "Outlier", function(obj)
{
   return(obj@wt)
})

setMethod("getOutliers", "Outlier", function(obj)
{
   return(which(obj@flag == 0))
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("show", "Outlier", function(object){
    cat("\nCall:\n")
    print(object@call)
    cat("-> Method: ", object@method, "\n")
    if(is.list(object@singularity))
        cat(strwrap(robustbase:::singularityMsg(object@singularity, object@n.obs)), sep ="\n")

    fl <- getFlag(object)
    nout <- length(which(fl == 0))
    digits = max(3, getOption("digits") - 3)
    cat("\nNumber of outliers detected:", nout, "\n")

    print(which(fl==0))

    invisible(object)
})

##
## Follow the standard methods: show, summary, plot
##
setMethod("plot", signature(x="Outlier", y="missing"), function(x, y="missing",
                                class=1,
                                id.n=3,
                                ...){

    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 2))
    on.exit(par(op))

    ind <- getClassLabels(x, class)
    dist <- getDistance(x)[ind]
    const <- getCutoff(x)[class]
    flag <- getFlag(x)[ind]
    plot(dist, xlab = "Index", ylab = "Distance")
    abline(h = const)
    plot(flag, xlab = "Index", ylab = "0/1 weights", ylim = c(0, 1), ...)

})

.getGrouping <- function(grouping=NULL, n)
{
    stopifnot(!missing(n))

    if(missing(grouping) || is.null(grouping))
        grouping <- rep(0, n)

    if(length(grouping) == 1) {
        # this is the number of groups and the groups are of equal size
        ng = grouping
        ni = n/ng
        if(ng*ni < n)
            stop("nrow(x) is not divisible by the number of groups")
        grouping <- rep(0,0)
        for(i in 1:ng)
            grouping <- c(grouping, rep(i,ni))
    }else if(length(grouping) > 1 && length(grouping) < n) {
        # grouping contains a vector with the group sizes
        ng <- length(grouping)
        if(sum(grouping) != n)
            stop("nrow(x) is not equal to n1+n2+...+nn")

        gx <- rep(0,0)
        for(i in 1:ng)
            gx <- c(gx, rep(i,grouping[i]))
        grouping <- gx
    }

    if(n != length(grouping))
        stop("nrow(x) and length(grouping) are different")

    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))

    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),"are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    names(g) <- NULL

    list(grouping=g, counts=counts, lev=lev1, ng=ng, proportions=proportions)
}

## from mvoutlier
.pbb <- function(map = "bss.background", add.plot = FALSE, ...)
{
    all = get(map)
    xrange = c(min(all[, 1], na.rm = TRUE), max(all[, 1], na.rm = TRUE))
    yrange = c(min(all[, 2], na.rm = TRUE), max(all[, 2], na.rm = TRUE))
    if (!add.plot) {
        plot(1, 1, xlim = xrange, ylim = yrange, xlab = "", ylab = "",  ...)
    }
    lines(all, ...)
}

.test.outlier <- function()
{

    bsstop <- NULL
    load(system.file("data","bsstop.rda", package="mvoutlier"))
    load(system.file("data","bss.background.rda", package="mvoutlier"))

    x=bsstop[,5:14]

    ## visualize multivariate outliers in the map
    op <- par(mfrow=c(2,2))

    ## identify multivariate outliers
    x.out=OutlierPCOut(x)
    cat("\nPCOUT: ", length(which(!getFlag(x.out))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.out) + 2)
    title("PCOUT")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))


    ## compare with outlier detection based on MCD:
    x.mcd=OutlierMahdist(x, control="mcd")
    cat("\nMCD: ", length(which(!getFlag(x.mcd))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO,bsstop$YCOO,pch=16,col=getFlag(x.mcd) + 2)
    title("MCD")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    ## compare with outlier detection based on PCDIST:
    x.pcdist <- OutlierPCDist(x, explvar=0.99)
    cat("\nPCDIST: ", length(which(!getFlag(x.pcdist))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.pcdist)+2)
    title("PCDIST")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    ## compare with outlier detection based on SIGN1:
    x.sign1 <- OutlierSign1(x)
    cat("SIGN1: ", length(which(!getFlag(x.sign1))), "\n")
    .pbb(asp=1)
    points(bsstop$XCOO, bsstop$YCOO, pch=16, col=getFlag(x.sign1)+2)
    title("SIGN1")
    legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))

    par(op)
}
