import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

from . import r_gating_general

r = ro.r

r('''

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

clust_2d_gate <- function(data, nr, nc, channel_names,
                        gating_channels, target, K = 2, quantile = 0.95) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs = data)
    gating_channels <- unlist(gating_channels)
    target <- unlist(target)
    g <- openCyto:::.flowClust.2d(fr, channels = gating_channels,
                                K = K, target = target, quantile = quantile)
    res <- Subset(fr, g)
    return(exprs(res))
}

transitional_gate <- function(data, nr, nc, channel_names,
                            gating_channels, target, K = 2,
                            quantile = 0.95, translation = 0.15) {
    data <- matrix(unlist(data), ncol = nc, nrow = nr)
    channel_names <- unlist(channel_names)
    colnames(data) <- channel_names
    fr <- flowFrame(exprs = data)
    gating_channels <- unlist(gating_channels)
    target <- unlist(target)
    g <- openCyto:::.flowClust.2d(fr, channels = gating_channels,
                                K = K, transitional = TRUE, target = target,
                                quantile = quantile, translation = translation,
                                pp_res = NULL)
    res <- Subset(fr, g)
    return(exprs(res))
}

''')

def singlet_gate(
        data_to_gate: dict[str, np.ndarray],
        gating_channels: list[str],
        r_ready: bool = False,
        **kwargs
    ) -> dict[str, np.ndarray]:
    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")
    
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        fcs_data = data_copy[key]
        fcs_data_nparr = np.asarray(fcs_data)

        r_func = ro.globalenv['singlet_gate']

        num_rows,num_cols = fcs_data_nparr.shape

        r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = num_rows, nc = num_cols, channel_names = list(fcs_data.channels), gating_channels = gating_channels)

        num_rows, num_cols = r_result.shape
        fcs_data = fcs_data[0:num_rows, :]
        fcs_data[:, :] = r_result
    return data_copy

def clust_2d_gate(
        data_to_gate:dict[str, np.ndarray],
        gating_channels: list[str],
        target: list[float, float],
        quantile: float = 0.95,
        K: bool = 2,
        r_ready: bool = False,
        **kwargs
    ) -> dict[str, np.ndarray]:
    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")

    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        fcs_data = data_copy[key]
        fcs_data_nparr = np.asarray(fcs_data)

        r_func = ro.globalenv['clust_2d_gate']

        num_rows,num_cols = fcs_data_nparr.shape

        r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = num_rows, nc = num_cols, channel_names = list(fcs_data.channels), gating_channels = gating_channels, target = target, K = K, quantile = quantile)

        num_rows, num_cols = r_result.shape
        fcs_data = fcs_data[0:num_rows, :]
        fcs_data[:, :] = r_result
    return data_copy

def transitional_gate(
        data_to_gate: dict[str, np.ndarray],
        gating_channels: list[str],
        target: list[float],
        quantile: float = 0.95,
        K: int = 2,
        translation: float = 0.15,
        r_ready: bool = False,
        **kwargs
    ) -> dict[str, np.ndarray]:
    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")

    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        fcs_data = data_copy[key]
        fcs_data_nparr = np.asarray(fcs_data)

        r_func = ro.globalenv['transitional_gate']

        num_rows,num_cols = fcs_data_nparr.shape

        r_result = r_func(
            data = fcs_data_nparr.T.tolist(), nr = num_rows,
            nc = num_cols, channel_names = list(fcs_data.channels),
            gating_channels = gating_channels, target = target,
            K = K, quantile = quantile, translation = translation)

        num_rows, num_cols = r_result.shape
        fcs_data = fcs_data[0:num_rows, :]
        fcs_data[:, :] = r_result
    return data_copy
    
def __test__():
    import FlowCal as fc
    ungated_data = {'test_fcs': fc.io.FCSData("./PreDox_B_002.fcs")}
    gated_data = singlet_gate(ungated_data, ["FSC-A", "FSC-H"], r_ready=True)
    print(gated_data.shape)
    gated_data = clust_2d_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4], r_ready=True)
    print(gated_data.shape)
    gated_data = transitional_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4], r_ready=True)
    print(gated_data.shape)