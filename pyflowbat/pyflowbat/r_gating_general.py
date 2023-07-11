import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import numpy as np
from .pyflowbat import SampleCollection
rpy2.robjects.numpy2ri.activate()

r = ro.r

def general_r_gate(
        data_to_gate: SampleCollection,
        gating_channel_names: list[str],
        r_file: str,
        r_function: str,
        arguments: dict,
        _r_ready: bool = False,
        **kwargs
    ) -> SampleCollection:
    """
    PyFlowBAT wrapper for an arbitrary R function.

    :param data_to_gate: the PyFlowBAT sample collection to gate;
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user
    :type data_to_gate: pyflowbat.pyflowbat.SampleCollection
    :param gating_channel_names: the channels to use for gating
    :type gating_channel_names: list[str]
    :param r_file: the path to the R file with the R function to use
    :type r_file: str
    :param r_function: the function in the R file to use
    :type r_function: str
    :param arguments: the arguments to pass to the R function
    :type arguments: dict
    :param _r_ready: a boolean describing whether or not R functionality
        has been initialized forthis Workspace;
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user
    :type _r_ready: bool
    :returns: the gated PyFlowBAT sample collection
    :rtype: pyflowbat.pyflowbat.SampleCollection
    """
    if not _r_ready:
        raise RuntimeError("R functionality has not been initialized")
    
    r['source'](r_file)
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        fcs_data = data_copy[key]
        fcs_data_nparr = np.asarray(fcs_data)

        r_func = ro.globalenv[r_function]

        nr,nc = fcs_data_nparr.shape

        r_result = r_func(data = fcs_data_nparr.T.tolist(), nr = nr, nc = nc, channel_names = list(fcs_data.channels), gating_channels = gating_channel_names, arguments = arguments)

        nr, nc = r_result.shape
        fcs_data = fcs_data[0:nr, :]
        fcs_data[:, :] = r_result
        data_copy[key] = fcs_data
    return data_copy