#!/usr/bin/env python

import re
import bz2
import logging
from collections import deque, defaultdict
import numpy as np
#import fire

from MEGAN.printing import io

logging.basicConfig(level=logging.INFO, filename="logfile", filemode="a+",
                    format="%(asctime)-15s %(levelname)-8s %(message)s")
class Parser:
    """
    For processing BLAST2LCA ko and taxonomy annotation outputs

    Allows one to parse via command line MEGAN tools blast2lca outputs:
    # Taxonomy
        'HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;'
    # KO - KEGG
        HISEQ:327:HN35KBCXX:2:1101:16910:2803/1; ;
        HISEQ:327:HN35KBCXX:2:1101:16768:6954/1; ; [1] K04079: 100 # 1
    Parameters:
    -----------

        rootPath: str
            path to the root directory of the project
        sampleDir: str
            sample directory which should be below the root directory
        sampleName: str
            the name of the sample so that the output files will be named as such
        taxfile: str
            name of the taxonomy file should be sitting below the root
        kofile: str
    """

    def __init__ (self, rootPath, sampleDir, sampleName, taxfile,kofile, verbose = True):
        #sample details
        self.verbose = verbose
        self.rootPath   = rootPath
        self.sampleName = sampleName
        self.sampleDir = sampleDir
        self.outputFile = "%s/%s/%s" % (rootPath, sampleDir, sampleName)
        logging.info("outputFile: %s" % self.outputFile)
        self.kofile     = "%s/%s/%s" % (rootPath, sampleDir, kofile)
        logging.info("kofile: %s" % self.kofile)
        self.taxfile    = "%s/%s/%s" % (rootPath, sampleDir, taxfile)
        logging.info("taxfile: %s" % self.taxfile)
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
            self.verbose
        )
        summary.printCombinedAnalysis()

    def LCA2neo4j(self, outfile = None):
        """
        For neo4j

        >>> rootDir = "/w/simulation_fr_the_beginning/reAssemble/everybodyelse/out/diamond/"
        >>> sampleDir = ""
        >>> tax = "newb29.allKOs.blast2lca.taxOutput"
        >>> ko = "newb29.allKOs.blast2lca.koOutput"
        >>> lcaparser = Parser(rootDir, sampleDir, "", tax, ko, False)
        >>> lcaparser.LCA2neo4j("neo4j_contig_taxa_lca_full")
        """
        self.__justTaxa()
        out = outfile if outfile is not None else self.outputFile
        summary = io(
            self.totalReads,
            self.rootPath,
            self.sampleName,
            self.sampleDir,
            self.reads,
            out,
            self.verbose
        )
        df = summary.Contig2LCAedgelist()
        return df

    def __parseKO(self):
        """
        HISEQ:327:HN35KBCXX:2:1101:16910:2803/1; ;
        HISEQ:327:HN35KBCXX:2:1101:16768:6954/1; ; [1] K04079: 100 # 1
        careful not all asignments are 1 to 1; maybe 1 query to multiple KOs.
        Have not incorporated that.
        Maybe some user flag to switch between taking the 1st KO (default) and double counting for the subsequent KOs, the default of MEGAN records up to 4 KOs
        K00001|contig00030; ; [1] K00001: 50.0 [2] K00100: 50.0 # 2
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
                #print(data)
                found = bool(re.search("K(\d{5})", data))
                if found:
                    ko = re.search("K(\d{5})", data).groups()[0]
                    self.kohash[ko] += 1
                    self.reads[readID]['ko'] = ko
                    #print("readID: %s ko: %s" % (readID, ko))
                else:
                    self.reads[readID]['ko'] = '00000'
                    #print("readID: %s ko: K00000" % readID)
        print("Total number of queries processed (ko): %s" % len(self.reads))
        return totalReads

    def __parseTAXA(self):
        """
        __parseTaxa()

        Parses the blast2lca outputs

        HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;

        """
        print("Processing TAXA....")
        print("Only storing the following taxonomic ranks: Domain, Phylum, Class, Order, Family, Genus, Species")
        ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        totalReads = 0
        isBZ = bool(re.search(".bz2$", self.taxfile))
        if isBZ:
            fh = bz2.BZ2File(self.taxfile)
        else:
            fh = open(self.taxfile, 'r')
        with fh as taxFile:
            for line in taxFile:
                phylahash = {}
                totalReads += 1
                elements    = deque(line.split(";"))
                readID      = elements.popleft()
                elements.popleft()
                for rank in ranks:
                    phylahash[rank] = 0 #default as unclassified
                while elements:
                    try:
                        hasAssignment = len(elements) != 1
                        if hasAssignment:
                            try:
                                assignment = elements.popleft()
                                rank, taxa = assignment.split("__")
                                score = int(elements.popleft())
                                if score < 50 :
                                    break
                                else:
                                    if rank in ranks:
                                        phylahash[rank] = taxa
                                        self.taxonomyhash[taxa] += 1
                                        #print("taxa: %s", taxa)
                            except (TypeError, ValueError) as err:
                                print("%s: this is the errorneous entry: %s" % (err, assignment))
                        else:
                            break

                    except IndexError:
                        break
                self.reads[readID]['taxa'] = phylahash
                #print("readID %s", (readID))
        print("Total number of queries processed (taxonomy): %s" % len(self.reads))
        return totalReads

    def __justTaxa(self, line = None):
        """
        justTaxa(line)

        >>> lines = ["HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;", "K00001|contig00060; ;d__2; 100;p__1224; 100;c__28216; 100;g__327159; 100;s__327160; 100;"]
        >>> for line in lines:
        >>>     justTaxa(line)
        """
        logging.info("Processing TAXA....")
        logging.info("Only storing the following taxonomic ranks: Domain, Phylum, Class, Order, Family, Genus, Species")
        ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        totalReads = 0
        isBZ = bool(re.search(".bz2$", self.taxfile))
        if isBZ:
            fh = bz2.BZ2File(self.taxfile)
        else:
            fh = open(self.taxfile, 'r')
        with fh as taxFile:
            queries = 0
            for line in taxFile:
                queries = queries + 1
                phylahash = {}
                elements    = deque(line.split(";"))
                readID      = elements.popleft()
                elements.popleft()
                for rank in ranks:
                    phylahash[rank] = 0 #default as unclassified
                while elements:
                    try:
                        hasAssignment = len(elements) != 1
                        if hasAssignment:
                            try:
                                assignment = elements.popleft()
                                rank, taxa = assignment.split("__")
                                score = int(elements.popleft())
                                if score < 50: #any lower I'll be uncertain
                                    break
                                else:
                                    if rank in ranks:
                                        phylahash[rank] = taxa
                                        self.taxonomyhash[taxa] += 1
                                        logging.debug("taxa: " +taxa)
                            except (TypeError, ValueError) as err:
                                print("%s: this is the errorneous entry: %s" % (err, assignment))
                        else:
                            break
                    except IndexError:
                        break
                noAssignment = True if np.sum([int(phylahash[ranks[i]]) for i in reversed(range(len(ranks)))]) == 0 else False
                if noAssignment:
                    logging.info("no assignment for this query: {} assigning to taxid:0 (unclassified)".format(readID))
                    self.reads[readID]['minTaxa'] = 0
                for index, taxid in enumerate([phylahash[ranks[i]] for i in reversed(range(len(ranks)))]):
                    if taxid != 0:
                        self.reads[readID]['minTaxa'] = taxid
                        logging.debug("this is the min taxa {}".format(taxid))
                        break
                logging.debug("readID %s", (readID))
        difference = queries - len(self.reads)
        logging.info("#queries missed (taxonomy): {} check log file".format(difference))


#if __name__ == '__main__':
    #fire.Fire(Parser)
