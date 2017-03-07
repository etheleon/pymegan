#!/usr/bin/env bash

#tests
#├── fullTest
#│   └── Day6_ANOX_Sample21
#│       └── Day6_ANOX_Sample21.daa
#└── trimmed
#    └── NUSM01AD00_M01_1_Day0
#        ├── KOoutput.bz2
#        ├── LSP008_ORF-tax.txt
#        ├── LSP008_ORF-tax.txt.bak
#        ├── NUSM01AD00_M01_1_Day0-combined.txt
#        ├── NUSM01AD00_M01_1_Day0.megan
#        └── taxoutput.bz2


python fullPipeline.py $PWD/tests/fullTest Day6_ANOX_Sample21 Day6_ANOX_Sample21 Day6_ANOX_Sample21.daa taxOutput koOutput
