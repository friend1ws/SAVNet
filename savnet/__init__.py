#! /usr/bin/env python

from parser import create_parser
from run import *

def main():

    parser = create_parser()
    args = parser.parse_args()
    expression_main(args)

