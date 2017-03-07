#!/usr/bin/env python

import subprocess

class Blast2lca:

    """
    input .daa (Diamond output file) should be MEGANISED first
    """

    def __init__ (self,
                  jar = '/export2/home/uesu/local/megan6/tools/blast2lca',
                  gi2kegg   = '/export2/home/uesu/local/megan6/tools/kegg/gi2kegg.map',
                  gi2taxid  = '/export2/home/uesu/simulation_fr_the_beginning/data/classifier/gi2taxid.refseq.map'):
        self.jar = jar
        self.gi2taxid = gi2taxid
        self.gi2kegg = gi2kegg

    def run(self, root, sampleDir, daaFileName, taxOutput='taxOutput', koOutput='koOutput', debug=False):
        #timeoutlimit = 7200
        daaPath = "%s/%s/%s" % (root, sampleDir, daaFileName)
        koPath = "%s/%s/%s" % (root, sampleDir, koOutput)
        taxPath = "%s/%s/%s" % (root, sampleDir, taxOutput)
        cmd = "%s -i %s -f DAA -k -g2kegg %s -g2t %s -o %s -ko %s -tid" % (self.jar, daaPath, self.gi2kegg, self.gi2taxid, taxPath, koPath)
        #try:
        if debug:
            print(cmd)
        else:
            subprocess.call(cmd, shell=True)
        #except subprocess.CalledProcessError as err:
            #print("blast2lca eorr:\n", err.output)
