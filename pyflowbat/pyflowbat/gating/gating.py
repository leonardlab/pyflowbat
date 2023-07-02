import numpy as np
import statsmodels.api as sm
from scipy.signal import find_peaks, peak_prominences

#####################
# PERCENTILE GATING #
#####################

def find_percentile(workspace, sample_collection, sample_name, channel, percentile, **kwargs):
    flow_sample = workspace.sample_collections[sample_collection][sample_name]
    return np.percentile(flow_sample[:, channel], percentile)

def gate_high_low(data_to_gate, channel, high = None, low = None, **kwargs):
    data_copy = data_to_gate.copy()
    for key in data_copy.keys():
        flow_sample = data_to_gate[key]
        if high is None:
            high = np.max(flow_sample[:, channel])
        if low is None:
            low = np.min(flow_sample[:, channel])
        flow_sample = flow_sample[flow_sample[:, channel] > low]
        flow_sample = flow_sample[flow_sample[:, channel] < high]
    return data_copy

##################
# SINGLET GATING #
##################

def gate_singlets(data_to_gate, a = 10**10, b = 2*10**4, **kwargs):
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
    return data_copy

##############
# HEK GATING #
##############

def find_HEK_gate(data_to_gate):
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

def gate_heks_helper(data_to_gate, manual_cutoffs, limits):
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

def gate_heks(data_to_gate, method, samples, limits, **kwargs):
    limits = [limits["FSC-A"][1], limits["SSC-A"][1]]
    data_copy = data_to_gate.copy()
    if method == 'unique':
        for key in data_to_gate.keys():
            data_copy[key] = gate_heks_helper(data_copy[key], None, limits)
    else:
        avg_FSC = []
        avg_SSC = []
        for i in samples:
            tmp1, tmp2 = find_HEK_gate(data_copy[i])
            avg_FSC.append(tmp1)
            avg_SSC.append(tmp2)
        avg_FSC = np.mean(avg_FSC)
        avg_SSC = np.mean(avg_SSC)
        for key in data_to_gate.keys():
            data_copy[key] = gate_heks_helper(data_copy[key], {'FSC-A': avg_FSC, 'SSC-A': avg_SSC}, limits)
    return data_copy