#!/usr/bin/python
#author: Mengzhu Liu

import pandas as pd
import numpy as np
import sys, os
import argparse
import Levenshtein
from time import time

def parse_args():
    parser = argparse.ArgumentParser(version='1.0', description='This is a filter program to analyze and filter reads with molecular barcode.')
    parser.add_argument("-i", dest = "tlxfile", type = str, required = True,
                              help = "tlx file" )
    parser.add_argument("-o", dest = "outputdir", type = str, required = True,
                              help = "output directory" )
    parser.add_argument("-n", dest = "BarcodeLength", type = str, required = True,
                              help = "random molecular barcode length" )
    
    args = parser.parse_args()
    return args
    
def rmDup(args):
    
    pd.options.mode.chained_assignment = None
        
    # initiated args  
    
    tlxfile = args.tlxfile
    wd = args.outputdir
    bl = int(args.BarcodeLength)
    
    # file basename
    
    basename = tlxfile.rpartition('/')[2].partition('.')[0]
    
    # pandas read in tlx file
    
    data = pd.read_table(args.tlxfile,sep = '\t')
    data.reset_index(level=data.index.names, inplace=True)
    
    # add prey length into the original dataframe
    
    #Pl <- data$Rend-data$Rstart+1
    data['P_len'] = data['Rend'] - data['Rstart'] + 1
    
    rawdata = data
    
    # get reads without Barcode as null
    
    null = data[data['Barcode'].isnull()]
    print 'Barcode Null Number {0}\n'.format(null['Qname'].count())
    data = data[data['Barcode'].notnull()]
    print 'Barcode Notnull Number {0}\n'.format(data['Qname'].count())
          
    #**************** Molecular Barcode filter method *****************#
    
    # #change barcode length
    # data['Barcode'] = data['Barcode'].str[:8]
    # print data['Barcode'].head()
    
    print "Starting Molecular Barcode filter method!\n"
    
    #~~~ drop truncated barcode(l<(bl-2)) ~~~#
    
    # get barcode length >=(bl-2) reads
    length = data['Barcode'].astype('str')
    length = data['Barcode'].str.len()
    data = data[length >=(bl-2)]
    print 'Barcode >=(bl-2) Number {0}\n'.format(data['Qname'].count())
    
    #~~~ analysis of barcode frequency ~~~#
    
    # get barcode frequency
    counts = pd.value_counts(data[u'Barcode'])
    counts = counts.to_frame()
    counts.reset_index(level=counts.index.names, inplace=True)
    counts.columns = ['barcode','counts']

    # add freq into the original dataframe
    data['freq'] = data.groupby('Barcode')['Barcode'].transform('count')
    
    #freq max and min value
    print ">=(bl-2) length Barcode freq range from {0} to {1}\n".format(counts['counts'].min(),counts['counts'].max())
    num_max = counts['counts'].max()
    num_min = counts['counts'].min()

    ###################################################################################
    #                Main Method to analyze molecular barcode                         #
    ###################################################################################
    # For barcode >=(bl-2) reads                                                      #
    # Find Rname, strand, Offset B_junction duplicates and rm dup barcode             #
    ###################################################################################

    # define offset, B_junction

    null.loc[:,'Offset'] = np.where(null['Strand'] == 1, null.Junction - null.Qstart, null.Junction + null.Qstart)
    null.loc[:,'B_junction'] = np.where(null['B_Strand'] == 1, null.B_Rend, null.B_Rstart)
    
    data.loc[:,'Offset'] = np.where(data['Strand'] == 1, data.Junction - data.Qstart, data.Junction + data.Qstart)
    data.loc[:,'B_junction'] = np.where(data['B_Strand'] == 1, data.B_Rend, data.B_Rstart)

    #~~~ step 2 ~~~#

    # get Rname, strand, Offset B_junction and p_len duplicates

    dup = data[data.duplicated(['Rname', 'Strand', 'Offset', 'B_junction'],keep=False)]
    other = data[~data['Qname'].isin(dup['Qname'])]
    other = other.drop_duplicates(['Rname', 'Strand', 'Offset', 'B_junction', 'Barcode'])
    dup = dup.drop_duplicates(['Rname', 'Strand', 'Offset', 'B_junction', 'Barcode'])
    dup = dup.sort_values(by=['Rname', 'Strand', 'B_junction', 'Offset'])
    dup = dup.reset_index(drop=True)

    # find barcodes should be rm
    
        
    print "Processing Method1 to compare barcode:\n"
    
    start_time = time()
    
    Qname_rm = pd.DataFrame({'Qname':['begin']})
    q = 0
    i = 0
    
    while i < dup['Qname'].count():
            
        n=0
        while dup.ix[i,'Rname'] == dup.ix[(i+n),'Rname'] and dup.ix[i,'Strand'] == dup.ix[(i+n),'Strand'] and \
        dup.ix[i,'B_junction'] == dup.ix[(i+n),'B_junction'] and dup.ix[i,'Offset'] == dup.ix[(i+n),'Offset']:
            n = n+1
            if i+n > dup['Qname'].count()-1:
                print 'out of range'
                break
        n=n-1
        # print "n {0}".format(n)
            
        for m in range(i,i+n):
            for j in range(m+1,i+n+1):
                if Levenshtein.distance(dup['Barcode'].iloc[m],dup['Barcode'].iloc[j]) <=2 :
                    Qname_rm.loc[q] = dup['Qname'].iloc[m]
                    q = q+1
                    # print "m {0}".format(m)
        i = i+n+1
        # print "i {0}".format(i)
    
    print ("\nBarcode align done in {}s".format(round(time()-start_time, 3)))
        
    #write to output
    
    dup = dup[~dup['Qname'].isin(Qname_rm['Qname'])]
    dup = dup.drop_duplicates(['Rname', 'Strand', 'Offset', 'B_junction', 'Barcode'])
    dup.to_csv(wd+'/'+basename+'_uni.tlx', header = True, sep = '\t', columns = [u'Qname', u'Rname', u'Strand', u'B_junction', u'Offset', u'Barcode', u'freq'])

    # clean reads (other + dup)
    data = other.append([dup])
    data = data.drop_duplicates(['Qname'])

    # reads removed
    duplicates = rawdata[~rawdata['Qname'].isin(data['Qname'])]

    # write result
    data.to_csv(wd+'/'+basename+'_dedup.tlx', header = True, sep = '\t', index=False,
       columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
       u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
       u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
       u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
       u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
       u'repeatseq', u'duplicate', u'B_junction', u'Offset', u'freq'])

    duplicates.to_csv(wd+'/'+basename+'_dup.tlx', header = True, sep = '\t', index=False,
       columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
       u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
       u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
       u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
       u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
       u'repeatseq', u'duplicate', u'B_junction', u'Offset', u'freq'])


def main():
    args = parse_args()
    rmDup(args)

main()