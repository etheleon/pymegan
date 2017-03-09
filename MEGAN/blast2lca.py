#!/usr/bin/env python

import subprocess
import re

class Blast2lca:

    """
    input .daa (Diamond output file) should be MEGANISED first
    """

    def __init__ (self,
                  jar = '/export2/home/uesu/local/megan6/tools/blast2lca',
                  #gi2kegg   = '/export2/home/uesu/local/megan6/tools/kegg/gi2kegg.map',
                  ref2kegg   = '/export2/home/uesu/github/tools/ncbi/ref2ko.map',
                  gi2taxid  = '/export2/home/uesu/simulation_fr_the_beginning/data/classifier/gi2taxid.refseq.map'):
        self.jar = jar
        self.gi2taxid = gi2taxid
        #self.gi2kegg = gi2kegg
        self.ref2kegg = ref2kegg


    def run(self, root, sampleDir, inputFile, taxOutput='taxOutput', koOutput='koOutput', debug=False):
        #timeoutlimit = 7200
        def testType(inputFile):
            if bool(re.search("\.daa", inputFile)):
                print("input format: DIAMOND")
                return "DAA"
            elif bool(re.search("\.m8", inputFile)):
                print("input format BlastTab")
                return "BlastTab"
            else:
                return False
        inputtype = testType(inputFile)

        if bool(inputtype):
            daaPath = "%s/%s/%s" % (root, sampleDir, inputFile)
            koPath = "%s/%s/%s" % (root, sampleDir, koOutput)
            taxPath = "%s/%s/%s" % (root, sampleDir, taxOutput)
            #cmd = "%s -i %s -f %s -k -g2kegg %s -g2t %s -o %s -ko %s -tid" % (self.jar, daaPath, inputtype, self.gi2kegg, self.gi2taxid, taxPath, koPath)
            cmd = "%s -i %s -f %s -k -a2kegg %s -g2t %s -o %s -ko %s -tid" % (self.jar, daaPath, inputtype, self.ref2kegg, self.gi2taxid, taxPath, koPath)
            #try:
            if debug:
                print(cmd)
            else:
                print(cmd)
                subprocess.call(cmd, shell=True)
            #except subprocess.CalledProcessError as err:
                #print("blast2lca eorr:\n", err.output)
