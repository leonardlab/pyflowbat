import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r
r['source']('gating.r')


def singlet_gate(data, gating_channels):

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv['singlet_gate']

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
    fcs_data[:, :] = r_result
    return fcs_data
    
def __test__():
    import FlowCal as fc
    ungated_data = fc.io.FCSData("./test-file-no-upload/PreDox_B_002.fcs")
    gated_data = singlet_gate(ungated_data, ["FSC-A", "FSC-H"])
    print(gated_data.shape)