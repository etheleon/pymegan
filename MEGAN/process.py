#!/usr/bin/env python

import re
import bz2
from collections import deque, defaultdict

from printing import io

class Parser:
    """
    For processing BLAST2LCA outputs

    Allows one to parse via command line MEGAN tools blast2lca outputs:
    # Taxonomy
        'HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;'
    # KO - KEGG
        HISEQ:327:HN35KBCXX:2:1101:16910:2803/1; ;
        HISEQ:327:HN35KBCXX:2:1101:16768:6954/1; ; [1] K04079: 100 # 1
    """

    def __init__ (self, rootPath, sampleDir, sampleName, kofile, taxfile, verbose):
        #sample details
        self.verbose = verbose
        self.rootPath   = rootPath
        self.sampleName = sampleName
        self.sampleDir = sampleDir
        self.outputFile = "%s/%s/%s" % (rootPath, sampleDir, sampleName)
        print("outputFile: %s" % self.outputFile)
        self.kofile     = "%s/%s/%s" % (rootPath, sampleDir, kofile)
        print("kofile: %s" % self.kofile)
        self.taxfile    = "%s/%s/%s" % (rootPath, sampleDir, taxfile)
        print("taxfile: %s" % self.taxfile)
        #data
        self.totalReads   = 0
        self.reads        = defaultdict(dict)
        self.kohash       = defaultdict(lambda: 1)
        self.taxonomyhash = defaultdict(lambda:1)

    def singleComparison(self):
        koreads     = self.__parseKO()
        taxareads   = self.__parseTAXA()
        if koreads == taxareads:
            self.totalReads = koreads
        else:
            raise ValueError("total number of reads dont tally")
        summary = io(
            self.totalReads,
            self.rootPath,
            self.sampleName,
            self.sampleDir,
            self.reads,
            self.outputFile,
            self.verbose
        )
        summary.printMeganSummary()

    def combined(self):
        koreads     = self.__parseKO()
        taxareads   = self.__parseTAXA()
        if koreads == taxareads:
            self.totalReads = koreads
        else:
            raise ValueError("total number of reads dont tally")
        summary = io(
            self.totalReads,
            self.rootPath,
            self.sampleName,
            self.sampleDir,
            self.reads,
            self.outputFile,
            self.io
        )
        summary.printCombinedAnalysis()

    def __parseKO(self):
        """
        HISEQ:327:HN35KBCXX:2:1101:16910:2803/1; ;
        HISEQ:327:HN35KBCXX:2:1101:16768:6954/1; ; [1] K04079: 100 # 1
        """
        print("Processing KOs....")
        totalReads = 0
        isBZ = bool(re.search(".bz2$", self.kofile))
        if isBZ:
            fh = bz2.BZ2File(self.kofile)
        else:
            fh = open(self.kofile, 'r')
        with fh as koFile:
            for line in koFile:
                totalReads += 1
                elements = deque(line.split(";"))
                readID = elements.popleft()
                elements.popleft()
                data = elements.popleft()
                found = bool(re.search("K(\d{5})", data))
                if found:
                    ko = re.search("K(\d{5})", data).groups()[0]
                    self.kohash[ko] += 1
                    self.reads[readID]['ko'] = ko
                else:
                    self.reads[readID]['ko'] = '00000'
        print("Total number of reads processed (ko): %s" % len(self.reads))
        return totalReads

    def __parseTAXA(self):
        """
        HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;
        """
        print("Processing TAXA....")
        totalReads = 0
        isBZ = bool(re.search(".bz2$", self.kofile))
        if isBZ:
            fh = bz2.BZ2File(self.taxfile)
        else:
            fh = open(self.taxfile, 'r')
        with fh as taxFile:
            for line in taxFile:
                phylahash = {}
                totalReads +=1
                elements    = deque(line.split(";"))
                readID      = elements.popleft()
                elements.popleft()
                while elements:
                    try:
                        if len(elements) == 1 :
                            break
                        #what is the case?
                        else:
                            try:
                                assignment = elements.popleft()
                                rank, taxa = assignment.split("__")
                                score = int(elements.popleft())
                                if score < 50 :
                                    break
                                else:
                                    phylahash[rank] = taxa
                                    self.taxonomyhash[taxa] += 1
                            except TypeError:
                                print("this is the assignment : %s" % assignment)
                    except IndexError:
                        break
                self.reads[readID]['taxa'] = phylahash
        print("Total number of reads processed (taxonomy): %s" % len(self.reads))
        return totalReads
