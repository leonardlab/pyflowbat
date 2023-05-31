import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r

def general_r_gate(data, gating_channels, r_file, r_function, arguments, r_ready = False, **kwargs):

    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")
    
    r['source'](r_file)

    fcs_data = data.copy()
    fcs_data_nparr = np.asarray(data)

    r_func = ro.globalenv[r_function]

    nr,nc = fcs_data_nparr.shape

    r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels, arguments = arguments)

    nr, nc = r_result.shape
    fcs_data = fcs_data[0:nr, :]
    fcs_data[:, :] = r_result
    return fcs_data