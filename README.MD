by-capture is a tool to assembly mitochondrial geneome or rDNA or other small gene regions, which depend on SPAdes assember.

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

https://github.com/sfu-compbio/mrsfast/

(2) seqtk

https://github.com/lh3/seqtk

(3) spades

https://github.com/ablab/spades

(4) mreps(optinal, only for rDNA)

http://mreps.univ-mlv.fr/index.html

(5) mira(optinal, alternative assembler)

http://www.chevreux.org/projects_mira.html

Please specify the these software in the script "begin.py" from row 25 to 30.

## Getting the source code

To obtain the source code from github, let us assume you want to clone this repo into a
directory named `Waveman`:

    git clone https://github.com/chenxing3/by-capture
    cd ./by-capture
    
## This software support fastq or compressed files as inputs. Reference sequence also required

please specify the sample name, reference file, and two reads files in the "sample.list" file.

Here, we summarize how to setup this software package and run the scripts on a small dataset from firefly for mitogenome assembly.


# A quick run for mitochondrial genome:

    $ python begion.py 

# A quick run for rDNA:

    $ python begion.py --organelle rDNA

The result file is contigs.fasta.filter.fasta in each iteration folder (such as /by-catpure/begin_assembly/tmp/iter_16/spades). 

## License
Code and audio data are available for research purposes only. For any other use of the software or data, please contact the us.

Our email: chenxing3753@qq.com
