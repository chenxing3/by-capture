by-capture is a tool to assembly mitochondrial geneome or rDNA or other small gene regions. 

## Citation


## Dependencies
depends on 
+ [Python 64-bit]
It wroten using python 3 and related packages:

(1) subprocess
(2) argparse
(3) tensorflow
(4) biopython
(5) multiprocessing



Besides python packages, the following software also required:

(1) mrsfast
(2) eqtk
(3) spades
(4) mreps(optinal, only for rDNA)
(5) mira(optinal, alternative assembler)


## This software support fastq or compressed files as inputs. Reference sequence also required

# A quick run for mitochondrial genome:

    $ python begion.py 

# A quick run for rDNA:

    $ python begion.py --organelle rDNA
