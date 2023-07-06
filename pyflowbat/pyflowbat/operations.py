import numpy as np

####################
# STATS OPERATIONS #
####################

def non_negative(
        x: float,
        **kwargs
    ) -> float:
    """
    Given a number x, returns the x if x is non-negative and 0 otherwise

    :param x: a number
    :returns: x if x is non-negative, else 0
    """
    return max(x, 0)

def apply_conversion_factor(
        x: float,
        factor: float,
        **kwargs
    ) -> float:
    """
    Given a number x and a conversion factor, returns x multiplied by the factor

    :param x: the value to convert
    :param factor: the conversion factor to apply to x
    :returns: x multiplied by the conversion factor
    """
    return x * factor

def compute_conversion_factor_stdErr(
        x: float,
        factor: float,
        factor_err: float,
        **kwargs
    ) -> float:
    """
    Computes the standard error of the conversion factor multiplication 

    :param x: the value being converted
    :param factor: the conversion factor being applied
    :param factor_err: the standard error of the conversion factor
    :returns: the standard error of the conversion
    """
    return x[0] * factor * np.sqrt((x[1]/x[0])**2+(factor_err/factor)**2)

####################
# STATS EXTRACTION #
####################

def split_sample_name(
        name: str,
        by: str,
        index: int,
        **kwargs
    ) -> str:
    """
    Splits a string by a character 

    :param str: the string to split
    :param by: the character to split by
    :param index: the index of the split string to return
    :returns: the string at the specified index of the provided string split by the specified character
    """
    return name.split(by)[index]

def channel_mean(
        data: np.ndarray,
        channel: str,
        **kwargs
    ) -> float:
    """
    Computes the mean value of a given channel of FCS data

    :param data: the FCS data with the desired mean
    :param channel: the channel of the desired mean
    :returns: the mean value of the specified channel of the specified FCS sample
    """
    return np.mean(data[:, channel])
