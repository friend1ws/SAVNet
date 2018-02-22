#! /usr/bin/env python

import math, sys, os, logging
from subprocess import Popen, PIPE

logging.basicConfig(format='%(asctime)s %(message)s', datefmt="%Y-%m-%d %I:%M:%S", level=logging.INFO)

def soft_max(x):
    x_max = max(x)
    return(x_max + math.log(sum([math.exp(y - x_max) for y in x]))) 


def median(numbers):
    if len(numbers) == 0:
        logging.error("Vector of zero length was put to median function. Return 0")
        return 0
    return (sorted(numbers)[int(round((len(numbers) - 1) / 2.0))] + sorted(numbers)[int(round((len(numbers) - 1) // 2.0))]) / 2.0


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
        logging.error("Executable does not exist: " + executable)
        sys.exit(1) 

    return True


