#! /usr/bin/env python

import math, sys, os
from subprocess import Popen, PIPE

from .logger import get_logger
logger = get_logger(__name__)


def soft_max(x):
    x_max = max(x)
    return(x_max + math.log(sum([math.exp(y - x_max) for y in x]))) 


def median(numbers):
    if len(numbers) == 0:
        logger.warning("Vector of zero length was put to median function. Return 0")
        return 0
    # return (sorted(numbers)[int(round((len(numbers) - 1) / 2.0))] + sorted(numbers)[int(round((len(numbers) - 1) // 2.0))]) / 2.0
    return sorted(numbers)[int(round((len(numbers) - 1) / 2.0))]


def quartile(numbers, prob):
    if len(numbers) == 0:
        logger.warning("Vector of zero length was put to quartile function. Return 0")
        return 0
    return sorted(numbers)[int((len(numbers) - 1) * prob)]


def generate_configurations(dim):

    conf = [[0], [1]]

    for k in range(dim - 1):

        new_conf = []
        for elm in conf:
            new_conf.append(elm + [0])
            new_conf.append(elm + [1])
            conf = new_conf

    return conf


def is_tool(executable):

    proc = Popen(["which", executable], stdout = PIPE, stderr = PIPE)
    out, err = proc.communicate()
    if proc.returncode == 1:
        logger.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True


