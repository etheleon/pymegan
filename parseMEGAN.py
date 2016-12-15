#!/usr/bin/env python

from MEGAN.process import Parser
import argparse

parser = argparse.ArgumentParser(description='Command line tool for processing blast2lca outputs')

parser.add_argument('rootDir', metavar='root', type=str,
                    help="the root directory")

parser.add_argument('sampleDir', metavar='sampledir', type=str,
                    help="path to sample directory")

parser.add_argument('sampleName', metavar='sample', type=str,
                    help="sample name")

parser.add_argument('tax', metavar='taxonomy', type=str,
                    help='blast2lca tax output - has to be in taxIDs d__2')

parser.add_argument('ko', metavar='kegg', type=str,
                    help='blast2lca ko output filename')

parser.add_argument('--verbose', dest='verbose', action='store_true',
                    help='to switch on verbose mode')

args = parser.parse_args()

#print(args)
#lcaparser = Parser(
    #"/export2/home/uesu/mouseData/",#rootPath
    #"data/trimmed/NUSM01AD00_M01_1_Day0/", #sampledirectory
    #"NUSM01AD00_M01_1_Day0", #sampleName
    #"KOoutput",#kofile
    #"taxoutput3"#taxfile
#)

lcaparser = Parser(args.rootDir, args.sampleDir, args.sampleName, args.ko, args.tax, args.verbose)

lcaparser.singleComparison()
lcaparser.combined()
