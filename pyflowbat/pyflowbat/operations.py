import numpy as np
from FlowCal.io import FCSData

####################
# STATS OPERATIONS #
####################

def non_negative(
        x: float,
        **kwargs
    ) -> float:
    """
    Given a number x, returns the x if x is non-negative and 0 otherwise.

    :param x: a number
    :type x: float
    :returns: x if x is non-negative and 0 otherwise
    :rtype: float
    """
    return max(x, 0)

def apply_conversion_factor(
        x: float,
        factor: float,
        **kwargs
    ) -> float:
    """
    Given a number x and a conversion factor, returns x multiplied by the factor.

    :param x: the value to convert
    :type x: float
    :param factor: the conversion factor to apply to x
    :type factor: float
    :returns: x multiplied by the conversion factor
    :rtype: float
    """
    return x * factor

def compute_conversion_factor_stdErr(
        x: float,
        factor: float,
        factor_err: float,
        **kwargs
    ) -> float:
    """
    Computes the standard error of the conversion factor multiplication.

    :param x: the value being converted
    :type x: float
    :param factor: the conversion factor being applied
    :type factor: float
    :param factor_err: the standard error of the conversion factor
    :type factor_err: float
    :returns: the standard error of the conversion
    :rtype: float
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
    Splits a string by a character.

    :param name: the string to split
    :type name: str
    :param by: the character to split by
    :type by: str
    :param index: the index of the split string to return
    :type index: int
    :returns: the string at the specified index of the provided string split by the specified character
    :rtype: str
    """
    return name.split(by)[index]

def channel_mean(
        data: FCSData,
        channel_name: str,
        **kwargs
    ) -> float:
    """
    Computes the mean value of a given channel of FCS data.

    :param data: the FCS data with the desired mean
    :type data: FlowCal.io.FCSData
    :param channel_name: the name of the channel of the desired mean
    :type channel_name: str
    :returns: the mean value of the specified channel of the specified FCS sample
    :rtype: float
    """
    return np.mean(data[:, channel_name])
