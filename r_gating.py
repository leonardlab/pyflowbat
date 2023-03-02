import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r
r['source']('gating.r')

def singlet_gate(data, gating_channels, **kwargs):

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['singlet_gate']

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
    fcs_data[:, :] = r_result
    return fcs_data

def clust_2d_gate(data, gating_channels, target, quantile=0.95, K=2, **kwargs):

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['clust_2d_gate']

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels, target = target, K = K, quantile = quantile)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
    fcs_data[:, :] = r_result
    return fcs_data

def transitional_gate(data, gating_channels, target, quantile = 0.95, K = 2, translation = 0.15, **kwargs):

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['transitional_gate']

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels, target = target, K = K, quantile = quantile, translation = translation)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
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