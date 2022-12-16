library(flowCore)
library(openCyto)

singlet_gate <- function(data, nr, nc, channel_names, gating_channels, ...) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs=data)
    gating_channels <- unlist(gating_channels)
    g <- openCyto:::.singletGate(fr, channels = gating_channels)
    res <- Subset(fr, g)
    return(exprs(res))
}

clust_2d_gate <- function(data, nr, nc, channel_names, gating_channels, target, K=2, quantile=0.95) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs=data)
    gating_channels <- unlist(gating_channels)
    target <- unlist(target)
    g <- openCyto:::.flowClust.2d(fr, channels = gating_channels, K=K, target=target, quantile=quantile)
    res <- Subset(fr, g)
    return(exprs(res))
}

transitional_gate <- function(data, nr, nc, channel_names, gating_channels, target, K=2, quantile=0.95, translation = 0.15) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs=data)
    gating_channels <- unlist(gating_channels)
    target <- unlist(target)
    g <- openCyto:::.flowClust.2d(fr, channels = gating_channels, K=K, transitional=TRUE, target=target, quantile=quantile, translation=translation, pp_res = NULL)
    res <- Subset(fr, g)
    return(exprs(res))
}