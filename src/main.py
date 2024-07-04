import argparse
import re
from parameter import Parameter
from spaced import SpacedQmer, MultiSpacedQmer
from input import Input, Sequence, FileScan

def main():
    parser = argparse.ArgumentParser(description="MISSH application")
    parser.add_argument('-si', type=str, help="Single-end input file")
    parser.add_argument('-pi', nargs=2, help="Paired-end input files")
    parser.add_argument('-q', type=str, help="Path to spaced seeds file")
    parser.add_argument('-dirO', type=str, default="../output/", help="Output directory")

    args = parser.parse_args()

    dir_output = args.dirO
    param = Parameter()

    if args.si:
        if not param.init(args.si, ""):
            print("Please enter an input filename single-end: -si <AbsPathFile>")
            return
        sequence = True
    elif args.pi:
        if not param.init(args.pi[0], args.pi[1]):
            print("Please enter input filenames paired-end: -pi <AbsPathFile1> <AbsPathFile2>")
            return
        sequence = True
    elif args.q:
        lines = get_lines(args.q)
        correct_qmer = [bool(re.match("^1(0|1)*1$", line.strip())) for line in lines]
