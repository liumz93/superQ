#!/usr/bin/python
#encoding:UTF-8
#author: Mengzhu Liu
#Date: 2017.6.29 7:25PM
"""

"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#standard library package
import os
import sys
import optparse
from time import time
import shlex
from subprocess import Popen, PIPE
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastMultx(object):
    """
    
    
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "FastMultx 0.1"
    USAGE = "Usage: %prog -i INDEX -1 FASTQ_R1 -2 FASTQ_R2 [...]"
    
    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """
        
        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        # optparser.add_option('-p','--index1', dest="index1",
        #      help= "Path to index1 file")
        optparser.add_option('-i','--index', dest="index",
             help= "Path to index file to demultx different samples!")
        optparser.add_option('-1','--fastq_R1', dest="fastq_R1",
             help= "Path to the fastq_R1 file")
        optparser.add_option('-2', '--fastq_R2', dest="fastq_R2",
             help= "Path to the fastq_R2 file") 
        optparser.add_option('-c', '--cutadapt_opt', dest="cutadapt_opt", default= "-m 25 -q 30,30 --trim-n",
             help= "[facultative] cutadapt options for the qc step (facultative and quoted) (default: -m 25 -q 30,30 --trim-n)")
        optparser.add_option('-d', '--deML_opt', dest="deML_opt", default= "--mm 1",
             help= "[facultative] deML options for demultiplex index1 step (facultative and quoted) (default: -mm 1)")
        optparser.add_option('-m', '--multx_opt', dest="multx_opt", default= "-m 0 -d 2 -x -b",
             help= "[facultative] deML options for the demultx1 step (facultative and quoted) (default: -mm 1)")
        optparser.add_option('-t', '--thread', dest="thread", default=1,
             help= "[facultative] Number of thread to use (default: 1)")
        optparser.add_option('--skip_demultx1', dest="skip_demultx1", action="store_true", default=False,
             help= "[facultative] Skip the demultiplex index1 step (default: False)")
        optparser.add_option('-r', '--run', dest="run", action="store_true", default=False,
            help= "[facultative] Run command lines else the sotware run in demo mode (default: False)")  
        
        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return FastMultx (
            # index1=opt.index1,
            index2=opt.index,
            fastq_R1=opt.fastq_R1,
            fastq_R2=opt.fastq_R2,
            thread=int(opt.thread),
            deML_opt=opt.deML_opt,
            multx_opt=opt.multx_opt,
            cutadapt_opt=opt.cutadapt_opt,
            run=opt.run,
            skip_demultx1=opt.skip_demultx1)
        
        #~~~~~~~FONDAMENTAL METHODS~~~~~~~#
        
    def __init__(self,
        # index1=None,
        index2=None,
        fastq_R1=None,
        fastq_R2=None,
        forward_adapter="ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        # forward_adapter="GCCTAGACATTTGGGAAGGACTGACTCTCTG",
        reverse_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC",
        thread=1,
        deML_opt="--mm 1",
        multx_opt="-m 0 -d 2 -x -b",
        cutadapt_opt="-m 25 -q 30,30 --trim-n",
        run=False,
        skip_demultx1=False):
        """
        General initialization function for import and command line
        """
        
        ### Verifications
        # assert index1, "A path to the index1 file is mandatory"
        assert index2, "A path to the index file is mandatory"
        assert fastq_R1, "A path to the fastq_R1 and fastq_R2 files is mandatory"
        assert fastq_R2, "A path to the fastq_R1 and fastq_R2 files is mandatory"
    
    
        print("\nInitialize FastqMultiplex")
    
        ### Storing Variables
        # self.index1 = index1
        self.index2 = index2
        self.fastq_R1 = fastq_R1
        self.fastq_R2 = fastq_R2
        self.thread = int(thread)
        self.cutadapt_opt = cutadapt_opt
        self.multx_opt = multx_opt
        self.deML_opt = deML_opt
        self.forward_adapter = forward_adapter
        self.reverse_adapter = reverse_adapter
        self.run = run
        self.skip_demultx1=skip_demultx1
        self.basename = fastq_R1.rpartition('/')[2].partition('.')[0]
        
    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        General function launching 
        """
        
        start_time = time()

        if self.run:
            print ("\nRunning in real mode")
        else:
            print ("\nRunning in demo mode")
        
        count = self.process_data_preparing()
        
        # Generate Reports
        if self.run:
            print ("\nGenerate a Report")
            with open (self.basename+"_FastqMultx_report.csv", "w") as report:
                report.write ("Program {}\tDate {}\n".format(self.VERSION,str(datetime.today())))
                report.write ("\nRUN PARAMETERS\n")
                report.write("  Index basename:\t{}\n".format(self.index))
                report.write("  Fastq R1 path:\t{}\n".format(self.fastq_R1))
                report.write("  Fastq R2 path:\t{}\n".format(self.fastq_R2 if self.fastq_R2 else "None"))
                report.write("  Number of thread:\t{}\n".format(self.thread))
                report.write("  Cutadapt options:\t{}\n".format(self.cutadapt_opt))
                report.write("  Forward Adapter:\t{}\n".format(self.forward_adapter if self.foward_adapter else "None"))
                report.write("  Forward Adapter:\t{}\n".format(self.reverse_adapter if self.reverse_adapter else "None"))
                report.write("\nREAD COUNT PER CATEGORY\n")
                for key, value in count.items():
                    report.write("  {}:\t{}\n".format(key, value))

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        
        return(0)
            
    
    def process_data_preparing (self):
        
        os.system("mkdir trim")
        count = OrderedDict()
        
        # #### deML demultiplex index1 ####
        #
        # if self.skip_demultx1:
        #     print ("\nSkiping demultx1 step!")
        #     demultx_R1 = self.fastq_R1
        #     demultx_R2 = self.fastq_R2
        #
        # else:
        #     print "\nStarting demultiplexing with index1!"
        #
        #     deML_report = self.basename+"_deML_report.txt"
        #
        #     cmd = "deML {0} -i {1} -f {2} -r {3} -if1 {4} -o demultx1 ".format(self.deML_opt, \
        #     self.index1, self.fastq_R1, self.fastq_R2, "todemultiplex.i1.gz")
        #     print cmd
        #
        #     demultx_R1 = "demultx_*r1.fq.gz"
        #     demultx_R2 = "demultx_*r2.fq.gz"
        #
        #     with open (deML_report, "w") as fout:
        #         for line in self.yield_cmd(cmd):
        #             fout.write(line)
        #
        
        demultx_R1 = self.fastq_R1 
        demultx_R2 = self.fastq_R2  
         
        #### cut i5,i7 adapters ####

        print "\nStarting trimming adapters!"

        trimmed_R1 = self.basename+"_trim_1.fastq.gz"
        trimmed_R2 = self.basename+"_trim_2.fastq.gz"
        cutadapt_report = self.basename+"_trim_report.txt"

        cmd = "cutadapt {0} -a {1} -A {2} -o {3} -p {4} {5} {6}".format (self.cutadapt_opt,\
        self.forward_adapter, self.reverse_adapter, "trim/"+trimmed_R1, "trim/"+trimmed_R2, demultx_R1, demultx_R2)
        print cmd

        with open (cutadapt_report, "w") as fout:
            for line in self.yield_cmd(cmd):
                fout.write(line)

        # Extract information from cutadapt_report
        with open (cutadapt_report, "r") as fin:
            for line in fin:
                if line.startswith("Total reads processed:"):
                    count["Total reads before trimming"] = int(line.split()[-1].replace(",",""))
                if line.startswith("Reads with adapters:"):
                    count["Reads with adapters"] = int(line.split()[-2].replace(",",""))
                if line.startswith("Reads that were too short:"):
                    count["Reads that were too short"] = int(line.split()[-2].replace(",",""))
                if line.startswith("Reads written (passing filters):"):
                    count["Reads after trimming"] = int(line.split()[-2].replace(",",""))


        #### demultiplex index2 ####

        print "\nStart demultiplexing index2!"

        multx_report = self.basename+"_multx_report.txt"

        cmd = "fastq-multx {0} -B {1} {2} {3} -o {4} {5}".format (self.multx_opt,\
        self.index2, "trim/"+trimmed_R1, "trim/"+trimmed_R2, "%_R1.fq.gz", "%_R2.fq.gz")
        print cmd

        with open (multx_report, "w") as fout:
            for line in self.yield_cmd(cmd):
                fout.write(line)


                            
        return count


    def yield_cmd (self, cmd):

        # Split the commands
        split_cmd = shlex.split(cmd)

        # Prepare the popen object
        proc = Popen(split_cmd, stdout=PIPE)

        # yield results line by line
        for line in proc.stdout:

             yield line
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    fastq_demultiplex = FastMultx.class_init()
    fastq_demultiplex()