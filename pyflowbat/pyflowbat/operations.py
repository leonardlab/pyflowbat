import numpy as np

####################
# STATS OPERATIONS #
####################

def non_negative(x, **kwargs):
    return max(x, 0)

def apply_conversion_factor(x, factor, **kwargs):
    return x * factor

def compute_conversion_factor_stdErr(x, factor, factor_err, **kwargs):
    return x[0] * factor * np.sqrt((x[1]/x[0])**2+(factor_err/factor)**2)

####################
# STATS EXTRACTION #
####################

def split_sample_name(name, by, index, **kwargs):
    return name.split(by)[index]

def channel_mean(data, channel, **kwargs):
    return np.mean(data[:, channel])
