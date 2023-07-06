import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
rpy2.robjects.numpy2ri.activate()

r = ro.r

def general_r_gate(
        data_to_gate: dict[str, np.ndarray],
        gating_channels: list[str],
        r_file: str,
        r_function: str,
        arguments: dict,
        r_ready: bool = False,
        **kwargs
    ):
    if not r_ready:
        raise RuntimeError("R functionality has not been initialized")
    
    r['source'](r_file)
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        fcs_data = data_copy[key]
        fcs_data_nparr = np.asarray(fcs_data)

        r_func = ro.globalenv[r_function]

        nr,nc = fcs_data_nparr.shape

        r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channels, arguments = arguments)

        nr, nc = r_result.shape
        fcs_data = fcs_data[0:nr, :]
        fcs_data[:, :] = r_result
    return data_copy