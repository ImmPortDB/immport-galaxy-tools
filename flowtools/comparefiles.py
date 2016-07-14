#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import pandas as pd

from argparse import ArgumentParser

def compare_file(file_un, file_dx):
    df1 = pd.read_table(file_un)
    df2 = pd.read_table(file_dx)
    if df1.equals(df2):
        print("yay")
    else:
        print("ohno")

if __name__ == "__main__":
    parser = ArgumentParser(
             prog="CompareFile",
             description="compares 2 files")

    parser.add_argument(
            '-u',
            dest="input1",
            required=True,
            help="File location for the first file.")

    parser.add_argument(
            '-d',
            dest="input2",
            required=True,
            help="File location for the second file.")

    args = parser.parse_args()
    
    compare_file(args.input1, args.input2)
    sys.exit(0)
