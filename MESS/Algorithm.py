"""
Author: Haotian Hong
Date: 2025/4/22
Description: This module provides utility functions for calculating parameters for arbitrary Steady-State Free Precession (SSFP) sequences, including support for spoiler gradients.
"""


def arbitarySSFP(start_echo: int, end_echo: int, ascend=None, balance=False):
    """
    Calculate parameters for an arbitrary Steady-State Free Precession (SSFP) sequence.

    Args:
        start_echo (int): Starting echo index.
        end_echo (int): Ending echo index.
        ascend (bool, optional): Direction of the echo sequence. If None, it is determined automatically.
        balance (bool, optional): Whether the sequence is balanced.

    Returns:
        tuple: Parameters (a, b, c) for the SSFP sequence.

    Raises:
        ValueError: If invalid arguments are provided.
    """
    # Validate input arguments
    if start_echo == end_echo and ascend is None and not balance:
        raise ValueError(
            "when only single echo, the orientation must be specified")
    if balance and start_echo != 0 and end_echo != 0:
        raise ValueError("when balance, start_echo and end_echo must be 0")

    # Calculate the number of echoes and determine the direction
    num_echo = abs(start_echo - end_echo) + 1
    if num_echo != 1:
        ascend = (start_echo < end_echo)

    # Initialize parameters based on direction
    if ascend:
        idx_echo0 = -start_echo
        sum = -2
    else:
        idx_echo0 = start_echo
        sum = 2

    # Adjust parameters for balanced sequences
    if balance:
        sum = 0

    # Compute the parameters a, b, and c
    a = -(idx_echo0*2+1)
    b = num_echo*2
    c = sum-a-b

    return a, b, c


def arbitarySSFP_Spoiler(start_echo: int, end_echo: int, ascend=None, balance=False, spoiler_portion=0.0):
    """
    Calculate parameters for an arbitrary SSFP sequence with a spoiler gradient.

    Args:
        start_echo (int): Starting echo index.
        end_echo (int): Ending echo index.
        ascend (bool, optional): Direction of the echo sequence. If None, it is determined automatically.
        balance (bool, optional): Whether the sequence is balanced.
        spoiler_portion (float, optional): Portion of the spoiler gradient.

    Returns:
        tuple: Parameters (a, b, c) for the SSFP sequence with spoiler.
    """
    # Get base SSFP parameters
    a, b, c = arbitarySSFP(start_echo=start_echo,
                           end_echo=end_echo, ascend=ascend, balance=balance)

    # Calculate the number of echoes
    num_echo = abs(start_echo - end_echo) + 1

    # Adjust parameters for the spoiler gradient
    a = a*(1+spoiler_portion)+spoiler_portion
    b = ([2, spoiler_portion*2]*num_echo)[:-1]
    c = c*(1+spoiler_portion)+spoiler_portion

    return a, b, c
