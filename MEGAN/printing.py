#!/usr/bin/env python

import datetime
from collections import defaultdict
import sys

class io:
    def __init__ (self, totalReads, rootDir, sampleName, sampleDir, readInfo, outputFile, verbose):
        self.verbose = verbose
        self.root = rootDir
        self.outfile = outputFile
        self.readInfo = readInfo
        self.heading = [
                ('@Creator'         , 'blast2lca2megansummary'),
                ('@CreationDate'    , datetime.datetime.utcnow().strftime("%a %b %d %H:%M:%S UTC %Y")),
                ('@ContentType'     , 'Summary4'),
                ('@Names'           , sampleName),
                ('@BlastMode'       , 'BlastX'),
                ('@Uids'            , 1480988728000), #honestly what is this file for?
                ('@Sizes'           , totalReads),
                ('@TotalReads'      , totalReads),
                ('@AdditionalReads' , 0),
                ('@Collapse'        , 'Taxonomy        -1'),
                ('@Algorithm'       , 'Taxonomy        LCA'),
                ('@Parameters'      , "minScore=50.0 maxExpected='0.01' minPercentIdentity='0.0' topPercent=10.0 minSupportPercent=0.01 minSupport=393 weightedLCAPercent=80.0 minComplexity=0.0 fNames= { Taxonomy KEGG }"),
                ('@NodeStyle'       , 'Taxonomy        Circle'),
                ('@ColorTable'      , 'Fews8   White-Green'),
                ('@ColorEdits'      , ''),
            ]
        self.singleEnding = "END_OF_DATA_TABLE\n #SampleID       @Source\n %s.daa      %s/%s.daa\n" % (sampleName, sampleDir, sampleName)

    def printMeganSummary(self):
        print("Writing .megan summary file to output: %s" % self.outfile)
        with open("%s%s" % (self.outfile, ".megan"), 'w') as outfile:
            print("Printing header...")
            for elem in self.heading:
                line  = "%s\t%s" % (elem[0], elem[1])
                outfile.write(line + "\n")
            print("Printing taxa summary...")
            self.__summariseTaxa(outfile)
            print("Printing KO summary...")
            self.__summariseKO(outfile)
            outfile.write(self.singleEnding + "\n")

    def printCombinedAnalysis(self):
        translate = {
            'p': 'phylum',
            'c': 'class',
            'o': 'order',
            'f': 'family',
            'g': 'genus',
            's': 'species',
        }
        with open("%s%s" % (self.outfile, "-combined.txt"),'w') as outfile:
            for rank in ['p', 'c', 'o', 'f', 'g', 's']:
                rankDict = defaultdict(lambda: defaultdict(int))
                for indiv in self.readInfo:
                    try:
                        taxa    = self.readInfo[indiv]['taxa'][rank]
                        ko      = self.readInfo[indiv]['ko']
                        rankDict[taxa][ko] +=1
                    except KeyError:
                        if self.verbose:
                            sys.stderr.write("%s rank is not available for read: %s\n" % (rank, indiv))
                for taxon in rankDict:
                    for ko in rankDict[taxon]:
                        outfile.write("%s\t%s\tK%s\t%s\n" % (translate[rank], taxon, ko, rankDict[taxon][ko]))

    def __summariseTaxa(self, fh):
        taxahash = defaultdict(lambda:1)
        translate = {
            'p': 'phylum',
            'c': 'class',
            'o': 'order',
            'f': 'family',
            'g': 'genus',
            's': 'species',
        }
        for rank in ['p', 'c', 'o', 'f', 'g', 's']:
            print("Taxonomy Processing: %s" % translate[rank])
            for indiv in self.readInfo:
                try:
                    taxa = self.readInfo[indiv]['taxa'][rank]
                    taxahash[taxa] +=1
                except KeyError:
                    if self.verbose:
                        sys.stderr.write("%s rank is not available for read: %s\n" % (rank, indiv))
        for taxa in taxahash:
            line = "tax\t%s\t%s" % (taxa,taxahash[taxa])
            fh.write(line + "\n")
        print("%s taxa sumamrised"% len(taxahash))

    def __summariseKO(self, fh):
        kohash = defaultdict(lambda:1)
        print("Function: KEGG - Processing...")
        for read in self.readInfo:
            ko = self.readInfo[read]['ko']
            kohash[ko] += 1
        for ko in kohash:
            count = kohash[ko]
            line = "KEGG\t%s\t%s" % (ko,count)
            fh.write(line + "\n")
        print("%s ko sumamrised"% len(kohash))
