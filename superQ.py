#!/usr/bin/python
#encoding:UTF-8
#author: Mengzhu Liu
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
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class runPipeline(object):
    """


    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "superQ V1.0 by Mengzhu"
    USAGE = "Usage: %prog -i inputdir -m metadatafile [...]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-i','--inputdir', dest="inputdir",
             help= "input directory containing fastq files (usually ./preprocess) ")
        optparser.add_option('-m','--metadata', dest="metadata",
             help= "Path to metadata file")
        optparser.add_option('-w','--workdir', dest="workdir", default= "./pipeline_results",
             help= "[facultative] Path to run the pipeline (default: ./results)")
        optparser.add_option('-p', '--pipeline_opt', dest="pipeline_opt", default= "--random-barcode 17 --find-adapter 1",
             help= "[facultative] pipeline options for running TranslocPipeline.pl (default: --random-barcode 17 --find-adapter 1) use TranslocPipeline.pl -h to check")
        optparser.add_option('-t', '--thread', dest="thread", default=1,
             help= "[facultative] Number of thread to use (default: 1)")
        optparser.add_option('-n', dest = 'BarcodeLength', type = str, default=17,
             help = "Random Molecular barcode length (defualt:17)")


        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return runPipeline (
            meta=opt.metadata,
            wd=opt.workdir,
            inputdir=opt.inputdir,
            thread=int(opt.thread),
            pipeline_opt=opt.pipeline_opt,
            bl=opt.BarcodeLength)

        #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
        meta=None,
        wd="./pipeline_results",
        inputdir=None,
        thread=1,
        pipeline_opt="--random-barcode 17 --find-adapter 1",
        bl = "17"):
        """
        General initialization function for import and command line
        """

        ### Verifications
        assert meta, "A path to metadata file is mandatory"
        assert inputdir, "A input dir(containing fastq files) is mandatory"


        print("\nInitialize runPipeline\n")

        ### Storing Variables
        self.meta = meta
        self.wd = wd
        self.inputdir = inputdir
        self.thread = int(thread)
        self.pipeline_opt = pipeline_opt
        self.bl = bl

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        General function launching
        """

        start_time = time()

        count = self.process_experiment()

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))

        return(0)


    def process_experiment (self):
        
        count = OrderedDict()
        
        print "Welcome use superQ, by Mengzhu :)"
        
        # read in metadata file

        metadata = pd.read_table(self.meta,sep='\t')
        sample = metadata['Library']
        primer = metadata['Primer']
        adapter = metadata['Adapter']
        index = metadata['MID']
        assembly = metadata['Assembly']
        chromosome = metadata['Chr']
        start = metadata['Start']
        end = metadata['End']
        strand = metadata['Strand']
        # check metadata file completeness

        if (sample.count() != primer.count()) or \
           (sample.count() != adapter.count()) or \
           (sample.count() != index.count()) or \
           (sample.count() != assembly.count()) or \
           (sample.count() != chromosome.count()) or \
           (sample.count() != start.count()) or \
           (sample.count() != end.count()) or \
           (sample.count() != strand.count()):
            print sample.count()
            print primer.count()
            print adapter.count()
            print index.count()
            print assembly.count()
            print chromosome.count()
            print start.count()
            print end.count()
            print "Error:  metadata file is malformed! Please check!"
            exit()
    

        for i in range(0,sample.count()):
            
            # mkdir sample dir

            print "\nStart Processing Sample:"+sample[i]+"\n"
            os.system("mkdir " + self.wd)
            os.system("mkdir " + self.wd + "/" + sample[i])
            os.system("mkdir " + self.wd + "/" + sample[i] + "/" + sample[i] +"_sequences")

            # generate fasta file

            Seq = primer[i]
            with open(self.wd + "/" + sample[i] + "/" + sample[i] +"_sequences/primer.fa", "w") as file_output:
                file_output.write(">Primer"+"\n")
                file_output.write(Seq+"\n")
            Seq = adapter[i]
            with open(self.wd + "/" + sample[i] + "/" + sample[i] +"_sequences/adapter.fa", "w") as file_output:
                file_output.write(">Adapter"+"\n")
                file_output.write(Seq+"\n")
            Seq = index[i]
            with open(self.wd + "/" + sample[i] + "/" + sample[i] +"_sequences/mid.fa", "w") as file_output:
                file_output.write(">MID"+"\n")
                file_output.write(Seq+"\n")

            # call TranslocPipeline.pl

            cmd = "TranslocPipeline.pl --workdir {0} --read1 {1} --read2 {2} --assembly {3} --chr {4} --start {5} --end {6} --strand {7} --primer {8} --adapter {9} {10}".format(\
                  self.wd+"/"+sample[i], self.inputdir+"/"+sample[i]+"_R1.fq.gz", self.inputdir+"/"+sample[i]+"_R2.fq.gz", assembly[i], chromosome[i], start[i], end[i], strand[i], self.wd+"/"+sample[i]+"/"+sample[i]+"_sequences/primer.fa",\
                  self.wd+"/"+sample[i]+"/"+sample[i]+"_sequences/adapter.fa", self.pipeline_opt)
            print "\n"+cmd+"\n"

            os.system(cmd)
            

            # call super MB filter method

            cmd = "SuperFilter.py -i {0} -o {1} -n {2} > {3}".format(self.wd+"/"+sample[i]+"/"+sample[i]+".tlx", self.wd+"/"+sample[i], self.bl, self.wd+"/"+sample[i]+"/"+sample[i]+"_transloc_filter.log")
            print "\n"+cmd+"\n"
            os.system(cmd)

            # call uncut filter method
                          
            cmd = "{0} -i {1} -o {2} -n {3} > {4}".format("GermlineAnalysis_superQ.py", self.wd+'/'+sample[i]+'/'+sample[i]+'_uncut_raw.tlx', self.wd+"/"+sample[i], self.bl, self.wd+"/"+sample[i]+"/"+sample[i]+"_uncut_filter.log")          
            print "\n"+cmd+"\n"
            os.system(cmd)
                        
            # find small indels in uncut
            
            #set how many bases range from cutsite
            n = 5
            
            if strand[i] == '+':
                cutsite = end[i]
            else:
                cutsite = start[i]
                
            cmd = "FindSmallIndel.py -i {0} -o {1} -n {2} -c {3} -s {4} > {5}".format(self.wd+"/"+sample[i]+"/"+sample[i]+'_uncut_raw.tlx',self.wd+"/"+sample[i]+"/"+sample[i]+'_SmallIndel_raw.tlx',n,cutsite,strand[i], self.wd+"/"+sample[i]+"/"+sample[i]+"_rawIndel.log")
            print "\n"+cmd+"\n"
            os.system(cmd)
            
            cmd = "FindSmallIndel.py -i {0} -o {1} -n {2} -c {3} -s {4} > {5}".format(self.wd+"/"+sample[i]+"/"+sample[i]+'_uncut_raw_result.tlx',self.wd+"/"+sample[i]+"/"+sample[i]+'_SmallIndel_result.tlx',n,cutsite,strand[i], self.wd+"/"+sample[i]+"/"+sample[i]+"_resultIndel.log")
            print "\n"+cmd+"\n"
            os.system(cmd)
            
            # clean file
            # os.system("rm {0} {1} {2} {3} {4} {5}".format(self.wd+"/"+sample[i]+"/*txt",
                                                      # self.wd+"/"+sample[i]+"/*sam",
                                                      # self.wd+"/"+sample[i]+"/*bam",
                                                      # self.wd+"/"+sample[i]+"/*raw*",
                                                      # self.wd+"/*",
                                                      # "*txt"))
            # os.system("rm {0} {1}".format(self.wd+"/"+sample[i]+"/*sam", self.wd+"/"+sample[i]+"/*bam"))
            
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

    run_pipeline = runPipeline.class_init()
    run_pipeline()