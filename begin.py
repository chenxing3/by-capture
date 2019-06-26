#-*- coding: UTF-8 -*- 

##########################################
##                                      ##
##   Version 2 by Chx at 06/17/2019     ##
##                                      ##
##########################################

import os
import argparse
import time
import tensorflow as tf
from Main_process import PROCESS

def configure():
    flags = tf.app.flags
    # data
    flags.DEFINE_string('Workdir', './begin_assembly/', 'work dir')
    flags.DEFINE_string('SampleList', './sample.list', 'input paired sample')
    flags.DEFINE_string('TempDir', 'tmp/', "tempory file for analysis")
    # flags.DEFINE_string('ResultDir', 'Result/','Result dir')
    flags.DEFINE_string('assembler', 'spades','aseember software')  

    # software path
    flags.DEFINE_string('mrsfast', '/home/chx/NGS/Softwares/mrsfast/mrsfast', 'mapping software')
    flags.DEFINE_string('seqtk', '/home/chx/NGS/Softwares/seqtk/seqtk', 'extract sequence according to title software')
    flags.DEFINE_string('spades', '/home/chx/NGS/Softwares/SPAdes-3.13.1/spades.py', 'spades assembler')
    flags.DEFINE_string('cut', '/usr/bin/cut', 'cut software')
    flags.DEFINE_string('mreps', '/home/chx/NGS/Softwares/mreps/mreps', 'Repeat detect software')
    flags.DEFINE_string('mira', '/usr/bin/mira', 'mira assembler')

    ## parameters
    flags.DEFINE_integer('IterNumber', 16, 'Number of iteration')
    flags.DEFINE_integer('Filterlength', 150, 'initial length to filter short contigs')
    flags.DEFINE_integer('Add_length', 50, 'add a specific basepair for each iterative')
    flags.DEFINE_integer('SeqDepth', 3, 'initial sequence depth contigs')
    flags.DEFINE_integer('thread', 4, 'Number of thread')
    flags.DEFINE_integer('crop', 50, 'crop min length for mrsfast')
    flags.DEFINE_integer('disimilarity', 30, 'the percent of not similar between ref and reads')    
    flags.DEFINE_integer('AddDepth', 1, 'add sequence depth in iterative')
    
    ## debug
    flags.FLAGS.__dict__['__parsed'] = False
    return flags.FLAGS

def main(_):
    parser = argparse.ArgumentParser()
    parser.add_argument('--organelle', dest='organelle', type=str, default='mitochondria',
                        help='organelle: mitochondria, rDNA')
                    
    args = parser.parse_args()
    
    if args.organelle not in ['mitochondria', 'rDNA']:
        print('Invalid input: ', args.organelle)
        print('Please input a organelle: mitochondria, rDNA')
    else:
        start = time.time()
        print("\nby-capture begin to process reads...... \n####################################### ")
        model = PROCESS(configure())
        getattr(model, args.organelle)()
        end = time.time()
        print('Time cost: ', end-start)

if __name__ == '__main__':
    tf.app.run()
