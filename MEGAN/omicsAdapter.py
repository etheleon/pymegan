#!/usr/bin/env python

import logging
import re
#import os
from collections import deque
#import numpy
import pandas as pd
#import csv
import argparse
from Bio import SeqIO


#This is experimental as this point in time for contigs based analyses only
class NodeEdge:
    def __init__(self, root, sampleDir, fasta, taxOutput, koOutput):
        self.root      = root
        self.sampleDir = sampleDir
        self.taxOutput = taxOutput
        self.kofile = koOutput
        self.fasta = "%s/%s/%s" % (root, sampleDir, fasta)
        #self.df = pd.read_table("allContigs-combined.txt", names=["rank", "taxon", "ko", "count"])
        self.contigs = {}
        self.kohash = {}

    def __parseFA(self):
        '''
        the fasta for the concatenated
        '''
        with open(self.fasta,'r') as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                keggBIN = record.id.split("|")[0]
                newRecordID = re.sub(r"\|", ":", record.id)
                #print(newRecordID)
                self.contigs[newRecordID] = "ko:%s" % keggBIN
        contigDF = pd.DataFrame(self.contigs.items(), columns=['contig:ID','bin'])
        contigDF['l:label'] = 'contigs'
        self.contigDF = contigDF

    def __parseAnnotation(self, type='contig'):
        """
        HISEQ:327:HN35KBCXX:2:1101:16768:6954/1; ; [1] K04079: 100 # 1
        K00001|contig00030; ; [1] K00001: 50.0 [2] K00100: 50.0 # 2
        #only takes first assignment POC only
        """
        logging.info("Processing KOs....")
        totalReads = 0
        fh = open(self.kofile, 'r')
        with fh as koFile:
            for line in koFile:
                totalReads += 1
                elements = deque(line.split(";"))
                readID = elements.popleft()
                readID = re.sub(r"\|", ":", readID) #because neo4j probably doesnt welcome this
                elements.popleft()
                data = elements.popleft()
                #might have more KOs attached
                #print(data)
                found = bool(re.search("K(\d{5})", data))
                if found:
                    ko = re.search("K(\d{5})", data).groups()[0]
                    #self.kohash[ko] += 1
                    self.kohash[readID] = "ko:K%s" % ko.zfill(5)
                    #print("readID: %s ko: %s" % (readID, ko))
                else:
                    self.kohash[readID] = 'ko:K00000'
                    #print("readID: %s ko: K00000" % readID)
        #print("Total number of queries processed (ko): %s" % len(self.reads))
        return totalReads

    def outputNodes(self, outputFile = "contig.node"):
        self.__parseFA()
        self.contigDF.to_csv(outputFile, sep="\t", index=False)

    def outputEdges(self, outputFile="contig.rels"):
        self.__parseAnnotation()
        logging.info("Writing edges to %s" % outputFile)
        relsDF = pd.DataFrame(self.kohash.items(), columns=['contig:START_ID', 'ko:END_ID'])
        relsDF['relationship:TYPE'] = 'assignment'
        relsDF.to_csv(outputFile, sep="\t", index=False)
#outputs

## rels
#contig:ID  ko  type:string
#K00001:contig00002      ko:K00001       assignment
#K00001:contig00004      ko:K00001       assignment
#K00001:contig00005      ko:K00001       assignment
#K00001:contig00007      ko:K00001       assignment
#K00001:contig00008      ko:K00001       assignment
#K00001:contig00010      ko:K00001       assignment
#K00001:contig00011      ko:K00001       assignment
#K00001:contig00014      ko:K00001       assignment
#K00001:contig00015      ko:K00001       assignment

## Nodes
#contig:ID bin l:label
#K00001:contig00002     ko:K00001      contigs
#K00001:contig00004     ko:K00001      contigs
#K00001:contig00005     ko:K00001      contigs
#K00001:contig00007     ko:K00001      contigs
#K00001:contig00008     ko:K00001      contigs
#K00001:contig00010     ko:K00001      contigs
#K00001:contig00011     ko:K00001      contigs
#K00001:contig00014     ko:K00001      contigs
#K00001:contig00015     ko:K00001      contigs

if __name__ =="__main__":
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="output blast2lca results into node and edge")
    parser.add_argument('root', help="root directory")
    parser.add_argument('sampleDir', help="sampleDir directory")
    parser.add_argument('fasta', help="fasta file")
    parser.add_argument('taxOutput', help="taxOutput file")
    parser.add_argument('koOutput', help="koOutput file")

    args = parser.parse_args()

    ne = NodeEdge(args.root, args.sampleDir, args.fasta, args.taxOutput, args.koOutput)
    ne.outputNodes()
    ne.outputEdges()
