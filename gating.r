library(flowCore)
library(openCyto)

singlet_gate <- function(data, nr, nc, channel_names, gating_channels, ...) {
    library(openCyto)
    library(flowCore)
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs=data)
    gating_channels <- unlist(gating_channels)
    g <- openCyto:::.singletGate(fr, channels = gating_channels)
    res <- Subset(fr, g)
    return(exprs(res))
}