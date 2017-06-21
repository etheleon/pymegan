#!/usr/bin/env python

import subprocess
import re

class Blast2lca:

    """
    input .daa (Diamond output file should be MEGANISED first) or m8 tabbed blast format.
    """

    def __init__ (self,
            blast2lca = '/export2/home/uesu/local/megan6ce/tools/blast2lca',
            gi2kegg   = '/export2/home/uesu/local/megan6ce/tools/kegg/gi2kegg.map',
            gi2taxid  = '/export2/home/uesu/simulation_fr_the_beginning/data/classifier/gi2taxid.refseq.map',
            new=False):
        self.script= blast2lca
        self.gi2taxid = gi2taxid
        self.gi2kegg = gi2kegg
        self.new = new

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
            inputPath = "%s/%s/%s" % (root, sampleDir, inputFile)
            koPath = "%s/%s/%s" % (root, sampleDir, koOutput)
            taxPath = "%s/%s/%s" % (root, sampleDir, taxOutput)
            if(self.new==False):
                cmd = "%s -i %s -f %s -k -g2kegg %s -g2t %s -o %s -ko %s -tid" % (self.script, inputPath, inputtype, self.gi2kegg, self.gi2taxid, taxPath, koPath)
            else:
                cmd = "%s -i %s -f %s -k -a2kegg %s -a2t %s -o %s -ko %s -tid" % (self.script, inputPath, inputtype, self.gi2kegg, self.gi2taxid, taxPath, koPath)
            if debug:
                print(cmd)
            else:
                print(cmd)
                subprocess.call(cmd, shell=True)
