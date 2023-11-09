import numpy as np
import statsmodels.api as sm
from scipy.signal import find_peaks, peak_prominences
from .pyflowbat import Workspace, SampleCollection
from typing import Optional, Union

#####################
# PERCENTILE GATING #
#####################

def find_percentile(
        workspace: Workspace,
        sample_collection_name: str,
        sample_name: str,
        channel_name: str,
        percentile: float,
        **kwargs
    ) -> float:
    """
    Finds the specified percentile of the data in a given channel
    a PyFlowBAT sample.

    :param workspace: the PyFlowBAT Workspace containing the sample
        whose percentile is desired.
    :type workspace: pyflowbat.pyflowbat.Workspace
    :param sample_collection_name: the name of the sample collection
        containing the sample whose percentile is desired
    :type sample_collection_name: str
    :param sample_name: the name of the sample whose percentile
        is desired
    :type sample_name: str
    :param channel_name: the name of the channel whose percentile
        is desired
    :type channel_name: str
    :param percentile: the desired percentile
    :type percentile: float
    :returns: the desired percentile of the data in the specified channel
        of the PyFlowBAT sample
    :rtype: float
    """
    flow_sample = workspace.sample_collections[sample_collection_name][sample_name]
    return np.percentile(flow_sample[:, channel_name], percentile)

def gate_high_low(
        data_to_gate: SampleCollection,
        gating_channel_name: str,
        high: Optional[float] = None,
        low: Optional[float] = None,
        **kwargs
    ) -> SampleCollection:
    """
    Gates all samples in a collection within an upper and lower
    bound in a specified channel.

    :param data_to_gate: the PyFlowBAT sample collection to gate;
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user
    :type data_to_gate: pyflowbat.pyflowbat.SampleCollection
    :param gating_channel_name: the name of the channel to gate
        the samples
    :type gating_channel_name: str
    :param high: the upper bound
    :type high: Optional[float]
    :param low: the lower bound
    :type low: Optional[float]
    :returns: the gated PyFlowBAT sample collection
    :rtype: pyflowbat.pyflowbat.SampleCollection
    """
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        flow_sample = data_copy[key]
        if high is None:
            high = np.max(flow_sample[:, gating_channel_name])
        if low is None:
            low = np.min(flow_sample[:, gating_channel_name])
        flow_sample = flow_sample[flow_sample[:, gating_channel_name] > low]
        flow_sample = flow_sample[flow_sample[:, gating_channel_name] < high]
        data_copy[key] = flow_sample
    return data_copy

##################
# SINGLET GATING #
##################

def gate_singlets(
        data_to_gate: SampleCollection,
        a: float = 10**10,
        b: float= 2*10**4,
        **kwargs
    ) -> SampleCollection:
    """
    Gates all samples in a collection for singlets by drawing a tilted
    ellipse.

    :param data_to_gate: the PyFlowBAT sample collection to gat
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user
    :type data_to_gate: pyflowbat.pyflowbat.SampleCollection
    :param a: the long axis of the ellipse surrounding the singlets,
        defaults to 1*10**10
    :type a: float
    :param b: the short axis of the ellipse surrounding the singlets,
        defaults to 2*10**4
    :type b: float
    :returns: the gated PyFlowBAT sample collection
    :rtype: pyflowbat.pyflowbat.SampleCollection
    """
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        flow_sample = data_copy[key]
        copied_sample = np.copy(np.asarray([flow_sample[:, 'FSC-A'], flow_sample[:, 'FSC-H']])).T
        data_ch = np.asarray([copied_sample[:, 0], copied_sample[:, 1]]).T
        lin_reg_res = sm.OLS([0, np.median(copied_sample[:, 1])], [0, np.median(copied_sample[:, 0])]).fit()
        center = np.array(np.asarray([np.median(copied_sample[:, 0]), np.median(copied_sample[:, 1])]))
        theta = np.arctan(lin_reg_res.params[0])
        data_centered = data_ch - center
        R = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
        data_rotated = np.dot(data_centered, R.T)
        mask = ((data_rotated[:,0]/a)**2 + (data_rotated[:,1]/b)**2 <= 1)
        flow_sample = flow_sample[mask]
        data_copy[key] = flow_sample
    return data_copy

##############
# HEK GATING #
##############

def _find_HEK_gate(data_to_gate):
    """
    Computes cutoffs for a HEK gate in a single sample of
    HEK 293FT cells by finding the most extreme valleys
    in the FSC-A and SSC-A channels.

    Called by `gate_heks`.
    """
    data_copy = data_to_gate.copy()
    H_fsc, edges_fsc = np.histogram(data_copy[:, 'FSC-A'], bins=1024)
    inv_H_fsc = np.max(H_fsc)-H_fsc
    peaks_fsc, _ = find_peaks(inv_H_fsc)
    peaks_fsc, _ = find_peaks(inv_H_fsc, prominence = np.max(peak_prominences(inv_H_fsc, peaks_fsc)[0])-1)
    H_ssc, edges_ssc = np.histogram(data_copy[:, 'SSC-A'], bins=1024)
    inv_H_ssc = np.max(H_ssc)-H_ssc
    peaks_ssc, _ = find_peaks(inv_H_ssc)
    peaks_ssc, _ = find_peaks(inv_H_ssc, prominence = np.max(peak_prominences(inv_H_ssc, peaks_ssc)[0])-1)
    return edges_fsc[peaks_fsc[0]], edges_ssc[peaks_ssc[0]]

def _gate_heks_one_sample(data_to_gate, manual_cutoffs, limits):
    """
    Gates a single sample of HEK 293FT cells by finding the
    most extreme valleys in the FSC-A and SSC-A channels.

    Called by `gate_heks` on each sample in a SampleCollection.
    """
    max_FSC = limits[0]
    max_SSC = limits[1]
    if manual_cutoffs is not None:
        data_copy = data_to_gate.copy()
        gated_data = data_copy[np.logical_and(data_copy[:, 'FSC-A'] > manual_cutoffs['FSC-A'], data_copy[:, 'SSC-A'] > manual_cutoffs['SSC-A'])]
    else:
        data_copy = data_to_gate.copy()
        H_fsc, edges_fsc = np.histogram(data_copy[:, 'FSC-A'], bins=1024)
        inv_H_fsc = np.max(H_fsc)-H_fsc
        peaks_fsc, _ = find_peaks(inv_H_fsc)
        peaks_fsc, _ = find_peaks(inv_H_fsc, prominence = np.max(peak_prominences(inv_H_fsc, peaks_fsc)[0])-1)
        H_ssc, edges_ssc = np.histogram(data_copy[:, 'SSC-A'], bins=1024)
        inv_H_ssc = np.max(H_ssc)-H_ssc
        peaks_ssc, _ = find_peaks(inv_H_ssc)
        peaks_ssc, _ = find_peaks(inv_H_ssc, prominence = np.max(peak_prominences(inv_H_ssc, peaks_ssc)[0])-1)
        gated_data = data_copy[np.logical_and(data_copy[:, 'FSC-A'] > edges_fsc[peaks_fsc[0]], data_copy[:, 'SSC-A'] > edges_ssc[peaks_ssc[0]])]
    gated_data = gated_data[np.logical_and(gated_data[:, 'FSC-A'] < max_FSC, gated_data[:, 'SSC-A'] < max_SSC)]
    return gated_data

def gate_heks(
        data_to_gate: SampleCollection,
        limits: dict[str, list[float]],
        method: str = "same",
        samples: Optional[list[str]] = None,
        **kwargs
    ) -> SampleCollection:
    """
    Gates all samples in a collection for HEK 293FT by finding the
    most extreme valleys in the FSC-A and SSC-A channels.

    :param data_to_gate: the PyFlowBAT sample collection to gate;
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user
    :type data_to_gate: pyflowbat.pyflowbat.SampleCollection
    :param limits: the limits of the FSC-A and SSC-A channels;
        NOTE: this parameter is provided to by the
        `pyflowbat.pyflowbat.Workspace.apply_gate` method and should
        NOT be specified by the user, to specify this parameter's
        value, change the value of `limits` within the Workspace in
        which this gate is being applied
    :type limits: dict[str, list[float]]
    :param method: whether or not to use a unique gate for each
        sample in the collection or the same gate for all samples,
        options are \"unique\" or \"same\",
        defaults to \"same\"
    :type method: str
    :param samples: the samples to use to define the gate if the
        \"same\" method is used,
        defaults to None
    :type samples: Optional[list[str]]
    :returns: the gated PyFlowBAT sample collection
    :rtype: pyflowbat.pyflowbat.SampleCollection
    """
    limits = [limits["FSC-A"][1], limits["SSC-A"][1]]
    data_copy = data_to_gate.copy()
    if method == 'unique':
        for key in data_to_gate.keys():
            data_copy[key] = _gate_heks_one_sample(data_copy[key], None, limits)
    elif method == 'same':
        avg_FSC = []
        avg_SSC = []
        for i in samples:
            tmp1, tmp2 = _find_HEK_gate(data_copy[i])
            avg_FSC.append(tmp1)
            avg_SSC.append(tmp2)
        avg_FSC = np.mean(avg_FSC)
        avg_SSC = np.mean(avg_SSC)
        for key in data_to_gate.keys():
            data_copy[key] = _gate_heks_one_sample(data_copy[key], {'FSC-A': avg_FSC, 'SSC-A': avg_SSC}, limits)
    else:
        raise ValueError("method must be either \"unique\" or \"same\"")
    return data_copy