library(chiptsne2)
ct2 = exampleChIPtsne2.with_meta()
ct2 = dimReduceTSNE(ct2)
ct2.by_cell = split(ct2, "cell")
rowRanges(ct2)
ct2.by_peak = split(ct2, "peak_MCF10CA1_CTCF")

#### cbind ####

# rowToRowMat = rowToRowMat,
# colToRowMatCols = colToRowMatCols,
# name_VAR = name_VAR,
# position_VAR = position_VAR,
# value_VAR = value_VAR,
# region_VAR = region_VAR,
# fetch_config = fetch_config

setMethod("cbind", "ChIPtsne2", function(..., deparse.level=1) {
    args <- list(...)
    all.cv <- lapply(args, colVec, withDimnames=FALSE)
    all.ccm <- lapply(args, colToColMat, withDimnames=FALSE)
    all.rcm <- lapply(args, rowToColMat, withDimnames=FALSE)

    all.cv <- do.call(c, all.cv)
    all.ccm <- do.call(cbind, all.ccm)
    all.rcm <- do.call(rbind, all.rcm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.rv <- rowVec(ref, withDimnames=FALSE)
    ref.rrm <- rowToRowMat(ref, withDimnames=FALSE)
    ref.crm <- colToRowMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(ref.rv, rowVec(x, withDimnames=FALSE))
            || !identical(ref.rrm, rowToRowMat(x, withDimnames=FALSE))
            || !identical(ref.crm, colToRowMat(x, withDimnames=FALSE)))
        {
            stop("per-row values are not compatible")
        }
    }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, colVec=all.cv,
                                colToColMat=all.ccm, rowToColMat=all.rcm,
                                check=FALSE)
})

#### rbind ####



setMethod("rbind", "ChIPtsne2", function(..., deparse.level=1) {
    args <- list(...)
    all.rv <- lapply(args, rowVec, withDimnames=FALSE)
    all.rrm <- lapply(args, rowToRowMat, withDimnames=FALSE)
    all.crm <- lapply(args, colToRowMat, withDimnames=FALSE)

    all.rv <- do.call(c, all.rv)
    all.rrm <- do.call(rbind, all.rrm)
    all.crm <- do.call(cbind, all.crm)

    # Checks for identical column state.
    ref <- args[[1]]
    ref.cv <- colVec(ref, withDimnames=FALSE)
    ref.ccm <- colToColMat(ref, withDimnames=FALSE)
    ref.rcm <- rowToColMat(ref, withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(ref.cv, colVec(x, withDimnames=FALSE))
            || !identical(ref.ccm, colToColMat(x, withDimnames=FALSE))
            || !identical(ref.rcm, rowToColMat(x, withDimnames=FALSE)))
        {
            stop("per-column values are not compatible")
        }
    }

    old.validity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(old.validity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out, rowVec=all.rv,
                                rowToRowMat=all.rrm, colToRowMat=all.crm,
                                check=FALSE)
})

