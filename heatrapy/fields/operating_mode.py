# -*- coding: utf-8 -*-
from __future__ import unicode_literals
"""Contains the function operating_mode.

Used to sweep applied and removed fields

"""

import numpy as np


def operating_mode(mode, time_before, time_after, field_steps, freq, j):

    period = 1./freq

    # mode 1
    if mode == 'constant_right':
        time_interval = ((1 - time_after - time_before) * period / 2.)
        delta_t = (time_interval / (field_steps))

        return delta_t

    # mode 2
    if mode == 'accelerated_right':
        delta_t = ((1 / (2 * (1. / (1. - time_after - time_before)) *
                    freq * np.sqrt(field_steps))) *
                   (np.sqrt(j + 1) - np.sqrt(j)))

        return delta_t

    # mode 3
    if mode == 'decelerated_right':
        delta_t = ((1 / (2 * (1. / (1. - time_after - time_before)) *
                    freq * np.sqrt(field_steps))) *
                   (np.sqrt(field_steps - j) - np.sqrt(field_steps - j - 1)))

        return delta_t

    # mode 4
    if mode == 'constant_left':
        time_interval = ((1 - time_after - time_before) * period / 2.)
        delta_t = (time_interval / (field_steps))

        return delta_t

    # mode 5
    if mode == 'accelerated_left':
        delta_t = ((1 / (2 * (1. / (1. - time_after - time_before)) *
                    freq * np.sqrt(field_steps))) *
                   (np.sqrt(field_steps - j) - np.sqrt(field_steps - j - 1)))

        return delta_t

    # mode 6
    if mode == 'decelerated_left':
        delta_t = ((1 / (2 * (1. / (1. - time_after - time_before)) *
                    freq * np.sqrt(field_steps))) *
                   (np.sqrt(j + 1) - np.sqrt(j)))

        return delta_t

    else:
        print 'invalid operating mode!'
