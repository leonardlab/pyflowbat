import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r
r['source']('gating.r')

def singlet_gate(data, gating_channels, r_ready, **kwargs):

    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['singlet_gate']

    num_rows,num_cols = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = num_rows, nc = num_cols, channel_names = list(fcs_data.channels), gating_channels = gating_channels)

    num_rows, num_cols = r_result.shape
    fcs_data = fcs_data[0:num_rows, :]
    fcs_data[:, :] = r_result
    return fcs_data

def clust_2d_gate(data, gating_channels, target, quantile=0.95, K=2, r_ready = False, **kwargs):

    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['clust_2d_gate']

    num_rows,num_cols = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = num_rows, nc = num_cols, channel_names = list(fcs_data.channels), gating_channels = gating_channels, target = target, K = K, quantile = quantile)

    num_rows, num_cols = r_result.shape
    fcs_data = fcs_data[0:num_rows, :]
    fcs_data[:, :] = r_result
    return fcs_data

def transitional_gate(data, gating_channels, target, quantile = 0.95, K = 2, translation = 0.15, r_ready = False, **kwargs):

    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['transitional_gate']

    num_rows,num_cols = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = num_rows, nc = num_cols, channel_names = list(fcs_data.channels), gating_channels = gating_channels, target = target, K = K, quantile = quantile, translation = translation)

    num_rows, num_cols = r_result.shape
    fcs_data = fcs_data[0:num_rows, :]
    fcs_data[:, :] = r_result
    return fcs_data
    
def __test__():
    import FlowCal as fc
    ungated_data = fc.io.FCSData("./test-file-no-upload/PreDox_B_002.fcs")
    gated_data = singlet_gate(ungated_data, ["FSC-A", "FSC-H"])
    print(gated_data.shape)
    gated_data = clust_2d_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4])
    print(gated_data.shape)
    gated_data = transitional_gate(ungated_data, ["FSC-A", "SSC-A"], target = [7.5*10**4, 5*10**4])
    print(gated_data.shape)