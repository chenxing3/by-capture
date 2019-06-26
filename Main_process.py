#-*- coding: UTF-8 -*- 

##########################################
##                                      ##
##   Version 2 by Chx at 06/17/2019     ##
##                                      ##
##########################################

import os
import re
import sys
import shutil
from Bio import SeqIO
from multiprocessing.dummy import Pool as ThreadPool
import subprocess
import time

class PROCESS(object):
    
    def __init__(self, conf):
        self.conf = conf
        if not os.path.exists(conf.Workdir):
            os.makedirs(conf.Workdir)
        self.tmpDir = self.conf.Workdir + conf.TempDir
        if not os.path.exists(self.tmpDir):
            os.makedirs(self.tmpDir)
        self.database = self.tmpDir + "sum.fastq"

    def seq_reformat(self, input, output):
        command_4seq = self.conf.seqtk + ' seq -l0 ' + input + ' > ' + output
        # print("command_4seq: ", command_4seq)
        os.system(command_4seq)


    def Read_input(self):
        for record in open(self.conf.SampleList, "r"):
            record = record.strip()
            records = record.split("\t")
            self.name = records[0]
            self.ref = records[1]
            
            # pass these step for mapping directly
            flag = 1
            for eachfile in records[2:]:
                input_file = eachfile
                check_file = self.tmpDir + os.path.basename(input_file)
                self.seq_reformat(input_file, check_file)
                command_rename  = 'cat ' + check_file + ' | awk ' + ''''{if(NR%4==1) $0=sprintf("@%012d/"'''\
                            + str(flag) +''',(1+i++)); print;}' >> ''' + self.database
                os.system(command_rename)
                os.remove(check_file)
                flag += 1

    def title_list(self, input):
        list_1 = input + '.1'
        list_2 = input + '.2'
        list_tmp1 = input+'.id.list'
        list_tmp2 = input+'.sort'
        # extract id for file 1 and file 2
        command_extract = '''awk -F '/' '{if($2 == 1 || $2 ==2){print $1}}' ''' + input + " > " + list_tmp1
        # filter repeat
        command_repeat = 'cut -f 1 ' + list_tmp1 + ' |sort -u  > ' + list_tmp2
        # generate list 1 and 2
        command_file1 = '''awk '{print $1"/1"}' ''' + list_tmp2 + ' > ' + list_1
        command_file2 = '''awk '{print $1"/2"}' ''' + list_tmp2 + ' > ' + list_2

        os.system(command_extract)
        os.system(command_repeat)
        os.system(command_file1)
        os.system(command_file2)
        return (list_1, list_2)
                     
    def mapping(self, mapdir, ref):        
        mrs_res = mapdir + "/mrs_res.out"
        log_file = mapdir + "/log"        
        ref_command = self.conf.mrsfast + ' --index ' + ref + \
                        ' --threads ' + str(self.conf.thread) + ' > ' + log_file
        os.system(ref_command) 
        mrsfast_command = self.conf.mrsfast + ' --search ' + ref + ' --seq ' + self.database + \
                            ' --crop ' + str(self.conf.crop) + ' --threads ' + str(self.conf.thread) + \
                            ' -e ' + str(self.conf.disimilarity) + ' -o ' + mrs_res + \
                            ' --disable-nohits >> ' + log_file
        os.system(mrsfast_command)
        if os.path.exists(mrs_res):
            result_list = mapdir + "/index.list"
            command_res = self.conf.cut + ' -f 1 ' + mrs_res + '| sort -u > ' + result_list
            os.system(command_res)
            result_lists  = self.title_list(result_list)
            iter_fastq1 = mapdir + '/iter_selected_1.fastq'
            iter_fastq2 = mapdir + '/iter_selected_2.fastq'
            command_list1_extract = self.conf.seqtk + ' subseq ' + self.database + ' ' + \
                                    result_lists[0] + ' > ' + iter_fastq1
            command_list2_extract = self.conf.seqtk + ' subseq ' + self.database + ' ' + \
                                    result_lists[1] + ' > ' + iter_fastq2
            os.system(command_list1_extract)
            os.system(command_list2_extract)
            return (iter_fastq1, iter_fastq2)        

    def assembly(self, inputs, assembly_dir, ref):
        if self.conf.assembler == 'mira':
            ## assembly in mira
            mira_dir = assembly_dir + '/mira/'
            if not os.path.exists(mira_dir):
                os.makedirs(mira_dir)
            os.chdir(mira_dir)
            mira_config = mira_dir + 'mira.conf'
            command_contents = '''project = mira
job = genome,denovo,accurate
parameters = -NW:cmrnl=warn -NW:cac=warn -NW:cmpm=warn -AS:sd=yes -GE:not=8 -SK:mmhr=2 SOLEXA_SETTINGS -AS:mrpc=5

readgroup = Reads

data = ''' + inputs[0] + ' ' + inputs[1] + '''
technology = solexa
strain = Seq
autopairing
template_size = 50 500 autorefine
segment_placement = ---> <---'''

            mira_conf_handle = open(mira_config, "w")
            mira_conf_handle.writelines(command_contents)
            mira_conf_handle.close()
            command_mira = self.conf.mira + ' ' + mira_config + ' > log'
            os.system(command_mira)
            mira_res_temp = mira_dir + 'mira_assembly/mira_d_results/mira_out.unpadded.fasta'
            ## title change
            mira_res = mira_dir + 'mira_assembly/mira_d_results/mira_out.unpadded.fasta.title'
            mira_res_handle = open(mira_res, 'w')
            for mira_line in open(mira_res_temp, 'r'):
                if re.search(r'^>', mira_line):
                    p=re.compile('\s+')
                    mira_line = re.sub(p,'|', mira_line) + '\n'
                    mira_res_handle.writelines(mira_line)
                else:
                    mira_res_handle.writelines(mira_line)
            mira_res_handle.close()
            if os.path.exists(mira_res):
                return mira_res
            else:
                print("mira failed!!")
                
        else:
            spades_dir = assembly_dir + '/spades/'
            if not os.path.exists(spades_dir):
                os.makedirs(spades_dir)
            # os.chdir(spades_dir)
            command_spades = self.conf.spades + ' --pe1-1 ' + inputs[0] + ' --pe1-2 ' + inputs[1] + \
                                ' --untrusted-contigs ' + ref + ' -o ' + spades_dir + ' > ' + \
                                spades_dir + 'assembly.log'
            # print("command_spades: ", command_spades)            
            os.system(command_spades)
            spades_res = spades_dir + 'contigs.fasta'
            if os.path.exists(spades_res):
                return spades_res
            else:
                print("spades failed!!")

    def repeative_fasta(self, _input, _output):
        _rep_res = open(_output, 'w')
        _rep_temp_log = _input + '.log'
        _clf_handle = open(_input, "rU")
        for _clf_record in SeqIO.parse(_clf_handle, 'fasta'):
            # print(_clf_record.id, len(_clf_record.seq))
            _flag = 0
            _starts, _ends = [], []
            _rep_temp_fasta = _input + '.temp.fasta'
            _rep_temp_log_handle = open(_rep_temp_fasta, "w")
            _rep_temp_content = '>' + _clf_record.id + '\n' + _clf_record.seq + '\n'
            _rep_temp_log_handle.writelines(_rep_temp_content)
            _rep_temp_log_handle.close()
            _command_mreps = self.conf.mreps + ' -res 5 -exp 3.0 -minsize 220 -from 1 -to ' + str(len(_clf_record.seq)-1) + ' -fasta ' + _rep_temp_fasta + ' > ' + _rep_temp_log
            os.system(_command_mreps)
    
            _rep_chk_log = open(_rep_temp_log, "r")
            for _rep_chk_line in _rep_chk_log:
                if re.search('^RESULTS: There are no', _rep_chk_line):
                    _flag = 1
                elif re.search('^RESULTS: There \w+ \d+', _rep_chk_line):
                   _flag = 2
            _rep_chk_log.close()
    
            if _flag == 1:
                _rep_res.writelines(_rep_temp_content)
            elif _flag == 2:
                print(_clf_record.id, "@@@@@")
                _rep_chk_log = open(_rep_temp_log, "r")
                for _rep_chk_line in _rep_chk_log:
                    if re.search(r'^\s+(\d+)\s+\-\>\s+(\d+)', _rep_chk_line):
                        _starts.append(re.search(r'^\s+(\d+)\s+\-\>\s+(\d+)', _rep_chk_line).group(1))
                        _ends.append(re.search(r'^\s+(\d+)\s+\-\>\s+(\d+)', _rep_chk_line).group(2))
                _rep_chk_log.close()
            if _starts != []:
                _start = int(min(_starts))
                _end = int(max(_ends))
                print(_start,_end)
                if _start < len(_clf_record.seq)/2 and _end < len(_clf_record.seq)/2:
                    _cutted_seq = '>' + _clf_record.id + '\n' + _clf_record.seq[_end:] + '\n'
                    _rep_res.writelines(_cutted_seq)
                elif _start >= len(_clf_record.seq)/2 and _end >= len(_clf_record.seq)/2:
                    _cutted_seq = '>' + _clf_record.id + '\n' + _clf_record.seq[:_start] + '\n'
                    _rep_res.writelines(_cutted_seq)
                else:
                    print("The seqence have two ends of tandem, we will exclude the sequence", _clf_record.id)
        _rep_res.close()

    def filter(self, contigs_file, length, seqdepth, type):
        filter_result = contigs_file + '.filter.fasta'
        filter_handle = open(filter_result, 'w')
        contig_handle = open(contigs_file, "rU")
        for contig_record in SeqIO.parse(contig_handle, 'fasta'):
            if self.conf.assembler == 'mira':
                contig_temp = re.findall(r'.*\|?cov=(\d+\.*\d*)\|?len=(\d+)\|?', contig_record.id)
                temp_length = float(contig_temp[0][1])
                temp_coverage = float(contig_temp[0][0])            
            else:
                contig_temp = re.findall(r'^NODE_\d+_length_(\d+)_cov_(\d+\.*\d*)', contig_record.id)
                temp_length = float(contig_temp[0][0])
                temp_coverage = float(contig_temp[0][1])               
            if temp_length >= length and temp_coverage >= seqdepth:
                ## for the rep filter to convert S,Y,R to N
                command_notATGC = '''/bin/echo "''' + contig_record.seq +  '''" | /bin/sed -e '/^[^>]/s/[^ATCGatcg]/N/g' '''
                temp_seq = subprocess.Popen([str(command_notATGC)], stdout=subprocess.PIPE, shell = True)
                temp_out = re.findall(r"\'(.*)\\n", str(temp_seq.stdout.read()))
                temp_fasta = str('>' + contig_record.id + '\n' + temp_out[0] + '\n')
                filter_handle.write(temp_fasta)
        filter_handle.close(), contig_handle.close()
        
        if type == 'mito':
            return filter_result
        elif type == 'rDNA':
            repeat_final = contigs_file + '.final.fasta'
            self.repeative_fasta(filter_result , repeat_final)
            return repeat_final
    
    def mitochondria(self):
        self.Read_input()
        ref = self.ref
        for iter in range(self.conf.IterNumber):
            print("\nStart iteration " + str(iter) + '......')
            iter_dir = self.tmpDir + 'iter_' + str(iter+1)
            if not os.path.exists(iter_dir):
                os.makedirs(iter_dir)
            # print("iter_dir: ", ref)
            fastq_files = self.mapping(iter_dir, ref)
            assembled_result = self.assembly(fastq_files, iter_dir, ref)
            result = self.filter(assembled_result, self.conf.Filterlength, self.conf.SeqDepth, 'mito')
            self.conf.Filterlength += self.conf.Add_length
            self.conf.SeqDepth += self.conf.AddDepth
            
            ref = result
            print("Done!')


    def rDNA(self):        
        self.Read_input()
        ref = self.ref
        for iter in range(self.conf.IterNumber):
            print("\nStarting iteration " + str(iter + 1) + '......')
            iter_dir = self.tmpDir + 'iter_' + str(iter+1)
            if not os.path.exists(iter_dir):
                os.makedirs(iter_dir)
            fastq_files = self.mapping(iter_dir, ref)
            assembled_result = self.assembly(fastq_files, iter_dir, ref)
            result = self.filter(assembled_result, self.conf.Filterlength, self.conf.SeqDepth, 'rDNA')
            self.conf.Filterlength += self.conf.Add_length
            self.conf.SeqDepth += self.conf.AddDepth
            
            ref = result
            print("Done!')
