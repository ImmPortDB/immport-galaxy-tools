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
    return
        
def compare_text(file1, file2):
    with open(file1, "r") as un, open(file2, "r") as dx:
        i = 0
        for line in un:
            ln = dx.readline()
            if not line == ln:
                print ("ohno")
                i += 1
                with open("compare_report.txt", "a") as rt:
                    rt.write(file1 + "\n" +  line + "\n" + file2 + "\n" + ln + "\n\n")
                if i==10:
                    sys.exit(2)
        if i == 0:
            print ("yay")
    return

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

    parser.add_argument(
            '-t',
            dest="as_text",
            help="compare files as text files")

    args = parser.parse_args()
    
    if args.as_text:
        compare_text(args.input1, args.input2)
    else:
        compare_file(args.input1, args.input2)
    sys.exit(0)
