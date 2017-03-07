# Python Module for parsing blast2lca outputs

this is a collection of python classes for dealing with MEGAN6 UE’s blast2lca outputs

## Usage

```
usage: parseMEGAN.py [-h] [--verbose] root sampledir sample taxonomy kegg

Command line tool for processing blast2lca outputs

positional arguments:
  root        the root directory
  sampledir   relative path from root directory to sample directory
  sample      sample name, could be same name as sample directory
  taxonomy    blast2lca taxonomy output filename - has to be in taxIDs d__2
  kegg        blast2lca ko output filename

optional arguments:
  -h, --help  show this help message and exit
  --verbose   to switch on verbose mode
```

Accepts `.bz2` files

```
#example
./parseMEGAN.py $PWD \ 
                tests/trimmed/NUSM01AD00_M01_1_Day0 \
                NUSM01AD00_M01_1_Day0 \
                taxoutput.bz2 \
                KOoutput.bz2
```

## Incorporating python class into your python script

### Initialize

```python
from MEGAN.process import Parser
lcaparser = Parser(
    "/export2/home/uesu/mouseData/",#rootPath
    "data/trimmed/NUSM01AD00_M01_1_Day0/", #sampledirectory
    "NUSM01AD00_M01_1_Day0", #sampleName
    "KOoutput",#kofile
    "taxoutput3"#taxfile
)
```

### Directory structure

```
/root/directory
└── sampleDirectory/
     ├── KOoutput
     ├── taxoutput
     └── sample.daa
```

### Input files

__NOTE__: The readID itself should not include any semi-colon character, if there is please remove it before running the parser.

#### KEGG
```
HISEQ:327:HN35KBCXX:2:2216:17590:101139/2; ; [1] K16363: 100 # 1
HISEQ:327:HN35KBCXX:2:2216:17505:101147/2; ;
HISEQ:327:HN35KBCXX:2:2216:18101:101081/2; ; [1] K16053: 100 # 1
HISEQ:327:HN35KBCXX:2:2216:18469:101047/2; ;
HISEQ:327:HN35KBCXX:2:2216:18275:101110/2; ; [1] K02621: 100 # 1
HISEQ:327:HN35KBCXX:2:2216:19012:101012/2; ;
HISEQ:327:HN35KBCXX:2:2216:19208:101184/2; ;
HISEQ:327:HN35KBCXX:2:2216:20173:101054/2; ;
HISEQ:327:HN35KBCXX:2:2216:20417:101000/2; ; [1] K02837: 100 # 1
HISEQ:327:HN35KBCXX:2:2216:20656:101185/2; ;
```

#### Taxonomy

```
HISEQ:327:HN35KBCXX:2:1101:17405:2046/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 80;g__838; 80;s__1262917; 20;
HISEQ:327:HN35KBCXX:2:1101:20056:2050/1; ;d__2; 100;p__1239; 100;c__186801; 100;o__186802; 100;f__186803; 76;g__189330; 40;s__1263073; 4;
HISEQ:327:HN35KBCXX:2:1101:19475:2499/1; ;d__2; 100;p__1239; 100;c__91061; 92;o__186826; 92;f__81852; 92;g__1350; 92;s__1351; 92;
HISEQ:327:HN35KBCXX:2:1101:16910:2803/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__815; 33;g__816; 33;s__1410607; 33;
HISEQ:327:HN35KBCXX:2:1101:20670:2813/1; ;d__2; 100;p__1239; 100;c__186801; 100;o__186802; 100;f__186803; 75;g__572511; 25;s__1226324; 13;
HISEQ:327:HN35KBCXX:2:1101:5294:3036/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 100;g__838; 100;s__1263102; 50;
HISEQ:327:HN35KBCXX:2:1101:10938:3042/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171552; 100;g__838; 67;s__52227; 67;
HISEQ:327:HN35KBCXX:2:1101:7569:5698/1; ;d__2; 59;p__1239; 33;c__186801; 33;o__186802; 33;f__541000; 22;g__1263; 19;s__40519; 7;
HISEQ:327:HN35KBCXX:2:1101:11823:5595/1; ;d__2; 100;p__1239; 92;c__186801; 68;o__186802; 68;f__31979; 20;g__1485; 20;s__1262806; 4;
HISEQ:327:HN35KBCXX:2:1101:7035:5916/1; ;d__2; 100;p__976; 100;c__200643; 100;o__171549; 100;f__171551; 100;g__283168; 100;s__1263090; 50;
```


### `.megan` summary format



```python
lcaparser.singleComparison()
```

### Combined format

We count reads which have the same taxons and KOs annotations

```python
lcaparser.combined()
```

```
phylum  67820   K00000  4
phylum  1224    K06937  1
phylum  1224    K00656  6
phylum  1224    K04564  2
phylum  1224    K06934  1
phylum  1224    K12524  24
phylum  1224    K00558  7
phylum  1224    K02674  1
phylum  1224    K06694  1
phylum  1224    K01785  12
...

...
species 1262910 K00033  1
species 1262911 K00000  525
species 35760   K00000  14
species 7462    K15421  1
species 7462    K00000  1
species 1262918 K00000  365
species 1262919 K02429  5
species 1262919 K07133  1
species 1262919 K03800  1
species 1262919 K00000  529
```

