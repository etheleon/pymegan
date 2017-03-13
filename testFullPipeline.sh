#!/usr/bin/env bash

echo "/export2/home/uesu/github/MEGAN/fullPipeline --rootDir /export2/home/uesu/simulation_fr_the_beginning/reAssemble/everybodyelse/out/diamond --sampleDir . --sampleName allContigs --inputFile newb29.allKOs.m8 --taxOutput meganWrapper.taxOutput --koOutput meganWrapper.koOutput"
echo "/export2/home/uesu/github/MEGAN/parseMEGAN /export2/home/uesu/simulation_fr_the_beginning/reAssemble/everybodyelse/out/diamond . allContigs meganWrapper.taxOutput meganWrapper.koOutput"
