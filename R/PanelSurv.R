##############################################################################
# Create a PanelSurv object
##############################################################################
PanelSurv <- function(ID, time, count) {
    if (sum(time <= 0) > 0)
        stop("Observation time must be positive.")

    index <- which(!duplicated(ID))
    N <- length(index)
    uniqueID <- ID[index]

    timeGrid <- sort(unique(time))

    panelMatrix <- matrix(NA, N, length(timeGrid))
    for (i in 1:N) {
        rowSet <- which(ID == uniqueID[i])
        panelMatrix[i, which(timeGrid %in% time[rowSet])] <- count[rowSet]
    }

    ps <- list(psDF=data.frame(ID=ID, time=time, count=count),
               timeGrid=timeGrid, panelMatrix=panelMatrix)
    class(ps) <- "PanelSurv"
    ps
}

is.PanelSurv <- function(x) inherits(x, "PanelSurv")

plot.PanelSurv <- function(x, ...) {
    ## time <- x$df$time
    ## ID <- x$df$ID
    ## count <- x$df$count
    ggplot(x$psDF, aes(time, ID)) + geom_tile(aes(fill=count))
}
