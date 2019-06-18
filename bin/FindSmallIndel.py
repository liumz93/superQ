#!/usr/bin/python
#encoding:UTF-8
#author: Mengzhu Liu
#date: 2017.10.23

import pandas as pd
import numpy as np
import sys, os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(version='1.0', description='This is a filter program to find small indels in sequence. by Mengzhu')
    parser.add_argument("-i", dest = "tlxfile", type = str, required = True,
                              help = "tlx file" )
    parser.add_argument("-o", dest = "outputfile", type = str, required = True,
                              help = "output file name" )
    parser.add_argument("-n", dest = "range", type = int, required = True,
                              help = "how many bp +- cutsite you wish to define small indels" )
    parser.add_argument("-c", dest = "cutsite", type = int, required = True,
                              help = "provide cutsite of your bait" )
    parser.add_argument("-s", dest = "strand", type = str, required = True,
                              help = "provide your primer strand" )  

    args = parser.parse_args()
    return args
    

def CigarParser(args,cutsite,strand,start,n,cigar):
    
    n = args.range
    SmIn = 0
    letter = re.findall('\D', cigar)
    number = re.findall('\d+', cigar)
    
    if strand ==  '-':
        letter = list(reversed(letter))
        number = list(reversed(number))
        
    count = 0
    for i in range(0,len(letter)):
        if letter[i] in ['M','D','N','X']:
            count = count + int(number[i])
            if count+start-1 >= cutsite-n:
                if letter[i] in ['D','X']:
                    SmIn = 1
                    break
                if (i+1 > len(letter)-1):
                    break
                if (count+start-1 < cutsite+n-1) and (letter[i+1] in ['D','I','X']):
                    SmIn = 1
                    break
                else:
                    break
    return SmIn

def FindSmallIndel(args):
    
    #read values from args
    tlxfile = args.tlxfile
    output = args.outputfile
    basename = tlxfile.rpartition('/')[2].partition('.')[0]
    n = args.range
    cutsite = args.cutsite
    strand = args.strand
    
    #read in tlx file
    data = pd.read_table(tlxfile, sep='\t')
    
    #define small indels according to B_Cigar value

    data['SmallIndel'] = data.apply(lambda row: 0 if row['uncut'] < 10 else CigarParser(args,cutsite,strand,row['B_Rstart'],n,row['B_Cigar']),axis=1)
    
    #output tlx file
    # data.to_csv(output,header=True,sep='\t',index=False)
    SmallIndel = data[data['SmallIndel'] == 1]
    SmallIndel.to_csv(output,header=True,sep='\t',index=False)
    print "{0} SmallIndel {1}".format (output,SmallIndel['Qname'].count())
    
def main():
    args = parse_args()
    FindSmallIndel(args)

main()