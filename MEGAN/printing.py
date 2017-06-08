#!/usr/bin/env python

import re
import datetime
from collections import defaultdict
import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.INFO)

class io:
    def __init__ (self, totalReads, rootDir, sampleName, sampleDir, readInfo, outputFile, verbose):
        self.verbose    = verbose
        self.root       = rootDir
        self.outfile    = outputFile
        self.readInfo   = readInfo
        self.sampleName = sampleName
        self.sampleDir  = sampleDir
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

    def Contig2LCAedgelist(self, altOutputfile = None):
        logging.info("Writing Edgelist for LCA assignments: %s" % self.outfile)

        readID_LCA_Map= {}
        for readID in self.readInfo:
            readID_LCA_Map[readID] = self.readInfo[readID]['minTaxa']
        df = pd.DataFrame.from_dict(readID_LCA_Map, orient='index')
        df.index.name = "readID"
        df.reset_index(inplace=True)
        df['readID'] = df['readID'].str.replace("\|", ":")
        filename = altOutputfile if altOutputfile is not None else self.outfile
        df.to_csv(filename, header=['contigid', 'taxid'], index=False)
        return df

    def printMeganSummary(self):
        logging.info("Writing .megan summary file to output: %s" % self.outfile)
        with open("%s%s" % (self.outfile, ".megan"), 'w') as outfile:
            logging.info("Printing header...")
            for elem in self.heading:
                line  = "%s\t%s" % (elem[0], elem[1])
                outfile.write(line + "\n")
            logging.info("Printing taxa summary...")
            self.__summariseTaxa(outfile)
            logging.info("Printing KO summary...")
            self.__summariseKO(outfile)
            outfile.write("END_OF_DATA_TABLE\n")
            outfile.write("#SampleID       @Source\n")
            outfile.write("%s.daa      %s/%s.daa\n" % (self.sampleName, self.sampleDir, self.sampleName))

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
            line = "TAX\t%s\t%s" % (taxa,taxahash[taxa])
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
