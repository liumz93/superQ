#!/usr/bin/python
#encoding:UTF-8
#author: Mengzhu Liu

import pandas as pd
import numpy as np
import sys, os
import argparse
import matplotlib
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Processing grab tlx that filtered!')
    parser.add_argument("-i", dest = "tlxfile", type = str, required = True,
                              help = "tlx file" )
    parser.add_argument("-o", dest = "outputdir", type = str, required = True,
                              help = "output directory" )
    parser.add_argument("-n", dest = "BarcodeLength", type = str, required = True,
                              help = "random molecular barcode length" )
    # parser.add_argument("-d", dest = "dedup", type = str, required = True,
    #                           help = "choose a dedup method: 1 or 2 ; 1: ld=1" )
    args = parser.parse_args()
    return args

def SuperFilter(args):
    
    tlxfile = args.tlxfile
    wd = args.outputdir
    basename = tlxfile.rpartition('/')[2].partition('.')[0]
    bl = int(args.BarcodeLength)
    
    # dedup = args.dedup
    # if dedup == "1":
    #     dedupmethod = 'DeDup_MB_v1.3.1.py'
    # else:
    #     dedupmethod = 'DeDup_MB_v1.3.2.py'

    #read in tlx
    data = pd.read_table(tlxfile, sep='\t')
    #,error_bad_lines=False, index_col=False, dtype='unicode')
                                                  # dtype={'Qname': object, 'JuncID': object, 'Rname': object, 'Junction': int,
                                                  # 'Strand': int, 'Rstart': int, 'Rend': int, 'B_Rname': object,
                                                  # 'B_Rstart': int, 'B_Rend': int, 'B_Strand': int,'B_Qstart': int,
                                                  # 'B_Qend': int,'Qstart': int,'Qend': int,'Qlen': int,'B_Cigar': object,
                                                  # 'Cigar': object,'Seq': object,'J_Seq': object,'Barcode': object,
                                                  # 'unaligned': int,'baitonly': int,'uncut': int,'misprimed': int,
                                                  # 'freqcut': int,'largegap': int,'mapqual': int,'breaksite': int,
                                                  # 'sequential': int,'repeatseq': int,'duplicate': int})
    print "Total {0}".format (data['Qname'].count())

    #unique Qname
    data = data.drop_duplicates(['Qname'])
    print "Unique_Qname {0}".format (data['Qname'].count())

    #************************** Process Filter *************************#

# ###test###
#     print "old pipeline:"
#     uncut = data[data['uncut'] == 1]
#     print "uncut    {0}".format (uncut['Qname'].count())
#     unaligned = data[data['unaligned'] == 1]
#     print "unaligned    {0}".format (unaligned['Qname'].count())
#     misprimed = data[data['misprimed'] < 10]
#     print "misprimed    {0}".format (misprimed['Qname'].count())
#     mapqual = data[data['mapqual'] == 1]
#     print "mapqual  {0}".format (mapqual['Qname'].count())
#     largegap = data[data['largegap'] > 30]
#     print "largegap {0}".format (largegap['Qname'].count())
#     sequential = data[data['sequential'] == 1]
#     print "sequential   {0}".format (sequential['Qname'].count())
#     baitonly = data[data['baitonly'] == 1]
#     print "baitonly {0}".format (baitonly['Qname'].count())
#     freqcut = data[data['freqcut'] == 1]
#     print "freqcut  {0}".format (freqcut['Qname'].count())
#     breaksite = data[data['breaksite'] == 1]
#     print "breaksite    {0}".format (breaksite['Qname'].count())
#     repeatseq = data[data['repeatseq'] == 1]
#     print "repeatseq    {0}".format (repeatseq['Qname'].count())
#     print "New Method:"
# ###test###

    # #~~~ First of ALL: filter junction is null and Adapters ~~~#
    #
    # data = data[data['Junction'].notnull()]
    # print 'None-Junction {0}'.format(data['Qname'].count())
    # data = data[data.Rname != 'Adapter']
    # print 'Adapter {0}'.format(data['Qname'].count())

    #~~~~ 1.unaligned ~~~~#

    unaligned = data[data['unaligned'] == 1]
    data = data[data['unaligned'] == 0]
    print "unaligned {0}".format (unaligned['Qname'].count())

    # # get rows without Barcode in unaligned
    # null = unaligned[unaligned['Barcode'].isnull()]
    # print "noBarcode_uaigned {0}".format (null['Qname'].count())


    #~~~~ 2.mapqual ~~~~#

    mapqual = data[data['mapqual'] == 1]
    data = data[data['mapqual'] == 0]
    print "mapqual  {0}".format (mapqual['Qname'].count())

    #~~~~ 3.misprimed ~~~~#

    misprimed = data[data['misprimed'] < 10]
    data = data[data['misprimed'] >= 10]
    print "misprimed    {0}".format (misprimed['Qname'].count())

    #~~~~ 4.uncut ~~~~#

    uncut = data[data['uncut'] >= 10]
    data = data[data['uncut'] < 10]
    print "germline {0}".format (uncut['Qname'].count())

    #~~~~ 5.baitonly ~~~~#

    baitonly = data[data['baitonly'] == 1]
    data = data[data['baitonly'] == 0]
    print "baitonly {0}".format (baitonly['Qname'].count())


    #output uncut and transloc raw data
    uncut.to_csv(wd+'/'+basename+'_uncut_raw.tlx',header=True,sep='\t',index=False,
    columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
    u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
    u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
    u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
    u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
    u'repeatseq', u'duplicate', u'B_junction', u'Offset', u'freq', u'P_len'])

    data.to_csv(wd+'/'+basename+'_transloc_raw.tlx',header=True,sep='\t',index=False,
    columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
    u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
    u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
    u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
    u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
    u'repeatseq', u'duplicate', u'B_junction', u'Offset', u'freq', u'P_len'])

    #~~~~ 6.Translocations Filter ~~~~#

    # print "{0} -i {1} -o {2} > {3}".format(dedupmethod,wd+'/'+basename+'_transloc_raw.tlx',wd,wd+'/'+basename+'_transloc_dup.log')
    os.system("{0} -i {1} -o {2} -n {3} > {4}".format("DeDup_MB_superQ.py",wd+'/'+basename+'_transloc_raw.tlx',wd,bl,wd+'/'+basename+'_transloc_dup.log'))

    #################################################################################################################
    # parallel: uncut analysis                                                                                      #
    #                                                                                                               #
    # "analyze uncut use command: DeDup_MB.py -i {0} -o {1} > transloc_dup.log".format(basename+'_uncut_raw.tlx',wd)#
    #                                                                                                               #
    #################################################################################################################

    
    #read in tlx
    data = pd.read_table(wd+'/'+basename+'_transloc_raw_dup.tlx', sep='\t')
    print "dup_MB {0}".format (data['Qname'].count())
    data = pd.read_table(wd+'/'+basename+'_transloc_raw_dedup.tlx', sep='\t')
    
    #~~~~ 7.sequential ~~~~#

    sequential = data[data['sequential'] == 1]
    data = data[data['sequential'] == 0]
    print "sequential {0}".format (sequential['Qname'].count())

    #~~~~ 8.largegap ~~~~#

    largegap = data[data['largegap'] > 30]
    data = data[data['largegap'] <= 30]
    print "largegap {0}".format (largegap['Qname'].count())

    #~~~~ 9.freqcut ~~~~#

    freqcut = data[data['freqcut'] == 1]
    data = data[data['freqcut'] == 0]
    print "freqcut {0}".format (freqcut['Qname'].count())

    #~~~~ 10.breaksite ~~~~#

    breaksite = data[data['breaksite'] == 1]
    data = data[data['breaksite'] == 0]
    print "breaksite  {0}".format (breaksite['Qname'].count())

    #~~~~ 11.repeatseq ~~~~#

    repeatseq = data[data['repeatseq'] == 1]
    data = data[data['repeatseq'] == 0]
    print "repeatseq  {0}".format (repeatseq['Qname'].count())

    #****************Transloc output *********************#
    
    data = data.sort_values(by=['Junction'])
    data.to_csv(wd+'/'+basename+'_transloc_result.tlx',header=True,sep='\t',index=False,
    columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
    u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
    u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
    u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
    u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
    u'repeatseq', u'duplicate', u'B_junction', u'Offset', u'freq', u'P_len'])
    print "final_hits {0}".format (data['Qname'].count())
    
    
def main():
    args = parse_args()
    SuperFilter(args)
main()