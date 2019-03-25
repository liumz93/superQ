#!/usr/bin/python
# Author: Mengzhu
# Date:
# Note: v1.2 + split barcode(Test)

import pandas as pd
import numpy as np
import sys, os
import argparse
import Levenshtein
import matplotlib
from time import time
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(version='1.0', description='This is a filter program to analyze and filter germline reads with molecular barcode.')
    parser.add_argument("-i", dest = "tlxfile", type = str, required = True,
                              help = "tlx file" )
    parser.add_argument("-o", dest = "outputdir", type = str, required = True,
                              help = "output directory" )
    parser.add_argument("-n", dest = "BarcodeLength", type = str, required = True,
                              help = "random molecular barcode length" )
    
    args = parser.parse_args()
    return args

def Uncut(args):
    
    pd.options.mode.chained_assignment = None

    # initiated args

    tlxfile = args.tlxfile
    wd = args.outputdir
    bl = int(args.BarcodeLength)
    ld = 2
    # defined_ambg = args.defined_ambg

    # file basename

    basename = tlxfile.rpartition('/')[2].partition('.')[0]

    # ~~ pandas read in tlx file ~~ #

    data = pd.read_table(args.tlxfile,sep = '\t')

    # #change barcode length
    # data['Barcode'] = 'AAAAAAAAAAAAAA'
    # print data['Barcode'].head()
    
    # ~~ define B_junc ~~ #
    data.loc[:,'B_junc'] = np.where(data['B_Strand'] == 1, data.B_Rend, data.B_Rstart)
    # rawdata = data
    
    print "Drop SEVERLY DAMAGED barcode!\n"
    
    # ~~ get reads without Barcode as null ~~ #
    null = data[data['Barcode'].isnull()]
    print 'Barcode_Null_Number {0}\n'.format(null['Qname'].count())
    data = data[data['Barcode'].notnull()]
    print 'Barcode_Number {0}\n'.format(data['Qname'].count())

    # ~~~~~~~~ catch barcode length >=(ld-2) as rawdata ~~~~~~~~ #
    
    # get barcode length >=(ld-2) reads ~~ #
    length = data['Barcode'].astype('str')
    length = data['Barcode'].str.len()
    data = data[length >=(bl-2)]
    print 'Barcode >=(n-2) Number {0}\n'.format(data['Qname'].count())
    
    # ~~ analysis of barcode frequency ~~ #

    # get barcode frequency
    counts = pd.value_counts(data[u'Barcode'])
    counts = counts.to_frame()
    counts.reset_index(level=counts.index.names, inplace=True)
    counts.columns = ['barcode','counts']
    counts.to_csv(wd+'/'+basename+'_barcode_counts.txt',header=True,sep='\t',index=False)

    # add freq into the original dataframe
    # data['freq'] = data.groupby('Barcode')['Barcode'].transform('count')
    data['freq'] = data.groupby(['B_junc','Barcode']).Barcode.transform('count')
    # data.to_csv(basename+'_freq.tlx', header = True, sep = '\t')
    
    ############# TEST Barcode Length #############
    
    # # add a barcode column to test barcode length
    # data['Barcode'] = data['Barcode'].map(lambda x: str(x)[0:9])
    
     ##############################################
    
    print data['Barcode'].head()
    
    counts = pd.value_counts(data[u'freq'])
    counts = counts.to_frame()
    counts.reset_index(level=counts.index.names, inplace=True)
    counts.columns = ['freq','counts']
    counts.to_csv(basename+'_freq_counts.txt',header=True,sep='\t',index=False)


    # freq max and min value
    print ">=(n-2) length Barcode freq range from {0} to {1}\n".format(counts['freq'].min(),counts['freq'].max())
    num_max = counts['freq'].max()
    num_min = counts['freq'].min()
    
    rawdata = data
    data = data.sort_values(by=['B_junc','Barcode'])
    # data.to_csv(wd+'/'+basename+'_l101112.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode',u'freq'])
    
    ################################################################
    # ~~~~~~~~ keep freq > num_max*0.05 as real barcode ~~~~~~~~~~ #
    ################################################################
    
    cutoff = num_max*0.05
    if cutoff < 3:
        cutoff = 3
    print "cutoff={0}".format(cutoff)
    Real = data[data['freq'] >= cutoff]
    data = data[data['freq'] < cutoff]
    Real = Real.drop_duplicates(['Qname'])
    Real = Real.sort_values(by=['B_junc','Barcode'])
    Real.reset_index(level=Real.index.names, inplace=True, drop=True)

    #############################################
    # ~~~~~~~~ Analysis l=12 Barocde ~~~~~~~~~~ #
    #############################################
        
    print "Begin analyze full length barcode!\n"
    
    # ~~ extract barcode l = 12 ~~ #
    
    Full = data.loc[length == bl]
    FullRaw = Full  
    barcode_list = map(str, range(bl))

    # ~~ split barcode in to single base within l = bl ~~ #
    
    x = range(0,bl)
    y = range(0,bl)
    couple = zip(x,y)
    for i,j in couple:
        Full.loc[:,barcode_list[i]] = Full['Barcode'].str[j]
    
    Full = Full.sort_values(by=['B_junc','Barcode']) 
     
    # ~~ drop only two base different barcode within l=bl ~~ #

    list_len = bl
    for i in range(1,bl):
        list_len = list_len - 1
        for j in range(1,list_len):
            column_list = ['B_junc']+ map(str, range(bl))
            del column_list[i]
            del column_list[j]
            # print "i {0} j {1}".format(i,j)
            Full = Full.drop_duplicates(column_list)

    Full = Full.sort_values(by=['B_junc','Barcode'])
    Full.reset_index(level=Full.index.names, inplace=True)
    
    
    
    #############################################
    # ~~~~~~~~ Analysis l=(bl-1) Barocde ~~~~~~~~~~ #
    #############################################

    print "Begin analyze l=(n-1) length barcode!\n"

    # ~~  extract barcode l = (bl-1) ~~ #

    Trunc1 = data.loc[length == (bl-1)]
    barcode_list = map(str, range(bl-1))

    # ~~  split barcode in to single base ~~ #

    x = range(0,(bl-1))
    y = range(0,(bl-1))
    couple = zip(x,y)
    for i,j in couple:
        Trunc1.loc[:,barcode_list[i]] = Trunc1['Barcode'].str[j]

    Trunc1 = Trunc1.sort_values(by=['B_junc','Barcode'])

    # ~~ drop only two base different barcode within l =(bl-1) ~~ #

    list_len = bl-1
    for i in range(1,(bl-1)):
        list_len = list_len - 1
        for j in range(1,list_len):
            column_list = ['B_junc']+map(str, range(bl-1))
            del column_list[i]
            del column_list[j]
            # print "i {0} j {1}".format(i,j)
            Trunc1 = Trunc1.drop_duplicates(column_list)

    Trunc1 = Trunc1.sort_values(by=['B_junc','Barcode'])
    Trunc1.reset_index(level=Trunc1.index.names, inplace=True)
    # Trunc1.to_csv(wd+'/'+basename+'_l12Uni.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])
    Trunc1Raw = Trunc1

    # ~~ drop l=(bl-1) and l=bl with the same B_junc and Levenshtein.distance = 2 ~~ #

    # get junctions in both l=11 and l=12
    Trunc1 = Trunc1[Trunc1['B_junc'].isin(Full['B_junc'])]
    Trunc1 = Trunc1.sort_values(by=['B_junc','Barcode'])
    Trunc1 = Trunc1.drop_duplicates(['Qname'])
    Trunc1.reset_index(level=Trunc1.index.names, inplace=True)
    # Trunc1.to_csv(wd+'/'+basename+'_l11tofilter.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])

    print "Processing Levenshtein.distance Filtering between l=(n-1) and l=bl\n"

    start_time = time()

    Qname_rm = pd.DataFrame({'Qname':['begin']})
    q = 0
    total1 = Trunc1['Qname'].count()
    total2 = Full['Qname'].count()

    # count l = 12 with the same junction
    i = 0
    j = 0
    while i < total1:

        while  Trunc1.ix[i, 'B_junc'] != Full.ix[(j) ,'B_junc']:
            j = j + 1
            # print 'j {0}'.format(j)

        m = 0
        while Trunc1.ix[i, 'B_junc'] == Full.ix[(j+m) ,'B_junc']:
            m = m + 1
            if j+m > total2-1:
                # print 'j+m out of range'
                break

        # print "l=12: j {0} j+m-1 {1}".format(j,j+m-1)

        # count l = 11 with the same junction
        n = 0
        while Trunc1.ix[i, 'B_junc'] == Trunc1.ix[(i+n) ,'B_junc']:
            n = n + 1
            if i+n > total1-1:
                # print 'i+n out of range'
                break

        # print "l=11: i {0} i+n-1 {1}".format(i,i+n-1)

        end1=i+n
        end2=j+m
        for i in range(i,end1):
            for j in range(j,end2):
                if Levenshtein.distance(Trunc1['Barcode'].iloc[i],Full['Barcode'].iloc[j]) <=3 :
                    Qname_rm.loc[q] = Trunc1['Qname'].iloc[i]
                    # print "i {0} Qname {1} Barcode {2}".format(i,Trunc1['Qname'].iloc[i],Trunc1['Barcode'].iloc[i])
                    q=q+1
                    break
        i = i+1
        j = j+1
        # print "i {0} j {1}".format(i,j)

    print ("l = n-1 Barcode Levenshtein.distance Align Done in {}s!\n".format(round(time()-start_time, 3)))

    Trunc1 = Trunc1Raw[~Trunc1Raw['Qname'].isin(Qname_rm['Qname'])]
    Trunc1 = Trunc1.drop_duplicates(['Qname'])
    Trunc1 = Trunc1.sort_values(by=['B_junc','Barcode'])
    Trunc1.reset_index(level=Trunc1.index.names, inplace=True)
    # Trunc1.to_csv(wd+'/'+basename+'_l11Uni_Result.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])

    #############################################
    # ~~~~~~~~ Analysis l=bl-2 Barocde ~~~~~~~~~~ #
    #############################################

    print "Begin analyze l=n-2 length barcode!\n"

    # ~~ extract barcode l = bl-2 ~~ #

    Trunc2 = data.loc[length == (bl-2)]
    Trunc2Raw = Trunc2
    barcode_list = map(str, range(bl-2))

    # ~~ split barcode in to single base ~~ #

    x = range(0,bl-2)
    y = range(0,bl-2)
    couple = zip(x,y)
    for i,j in couple:
        Trunc2.loc[:,barcode_list[i]] = Trunc2['Barcode'].str[j]

    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    # Trunc2.to_csv(wd+'/'+basename+'_l11.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode',u'Barcode1',u'Barcode2',u'Barcode3',
    #                                                                                u'Barcode4',u'Barcode5',u'Barcode6',u'Barcode7',u'Barcode8'])
    # ~~ drop only two base different barcode within l =10 ~~ #

    list_len = bl-2
    for i in range(1,(bl-2)):
        list_len = list_len - 1
        for j in range(1,list_len):
            column_list = ['B_junc']+map(str, range(bl-2))
            del column_list[i]
            del column_list[j]
            # print "i {0} j {1}".format(i,j)
            Trunc2 = Trunc2.drop_duplicates(column_list)

    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    Trunc2.reset_index(level=Trunc2.index.names, inplace=True)
    # Trunc2.to_csv(wd+'/'+basename+'_l11Uni.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])
    Trunc2Raw = Trunc2

    # ~~ drop l=bl-2 and l=bl-1 with Levenshtein.distance = 2 ~~ #

    # get junctions in both l=11 and l=12
    Trunc2 = Trunc2[Trunc2['B_junc'].isin(Trunc1['B_junc'])]
    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    Trunc2 = Trunc2.drop_duplicates(['Qname'])
    Trunc2.reset_index(level=Trunc2.index.names, inplace=True)
    # Trunc2.to_csv(wd+'/'+basename+'_l11tofilter.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])

    print "Processing Levenshtein.distance<=3 Filtering between l=(n-2) and l=(n-1)\n"

    start_time = time()

    Qname_rm = pd.DataFrame({'Qname':['begin']})
    q = 0
    total1 = Trunc2['Qname'].count()
    total2 = Trunc1['Qname'].count()

    # count l = 12 with the same junction
    i = 0
    j = 0
    while i < total1:

        while  Trunc2.ix[i, 'B_junc'] != Trunc1.ix[(j) ,'B_junc']:
            j = j + 1
            # print 'j {0}'.format(j)

        m = 0
        while Trunc2.ix[i, 'B_junc'] == Trunc1.ix[(j+m) ,'B_junc']:
            m = m + 1
            if j+m > total2-1:
                # print 'j+m out of range'
                break

        # print "l=12: j {0} j+m-1 {1}".format(j,j+m-1)

        # count l = 11 with the same junction
        n = 0
        while Trunc2.ix[i, 'B_junc'] == Trunc2.ix[(i+n) ,'B_junc']:
            n = n + 1
            if i+n > total1-1:
                # print 'i+n out of range'
                break

        # print "l=11: i {0} i+n-1 {1}".format(i,i+n-1)

        end1=i+n
        end2=j+m
        for i in range(i,end1):
            for j in range(j,end2):
                if Levenshtein.distance(Trunc2['Barcode'].iloc[i],Trunc1['Barcode'].iloc[j]) <=2 :
                    Qname_rm.loc[q] = Trunc2['Qname'].iloc[i]
                    # print "i {0} Qname {1} Barcode {2}".format(i,Trunc2['Qname'].iloc[i],Trunc2['Barcode'].iloc[i])
                    q=q+1
                    break
        i = i+1
        j = j+1
        # print "i {0} j {1}".format(i,j)

    print ("l = n-2 Barcode Levenshtein.distance<=3 Align Done in {}s!\n".format(round(time()-start_time, 3)))

    Trunc2 = Trunc2Raw[~Trunc2Raw['Qname'].isin(Qname_rm['Qname'])]
    Trunc2 = Trunc2.drop_duplicates(['Qname'])
    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    Trunc2.reset_index(level=Trunc2.index.names, inplace=True)
    # Trunc2.to_csv(wd+'/'+basename+'_l10Uni_Result1.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])
    Trunc2Raw = Trunc2

    # ~~ drop l=10 and l=12 with Levenshtein.distance = 3 ~~ #

    # get junctions in both l=11 and l=12
    Trunc2 = Trunc2[Trunc2['B_junc'].isin(Full['B_junc'])]
    Trunc2 = Trunc2.drop_duplicates(['Qname'])
    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    Trunc2.reset_index(level=Trunc2.index.names, inplace=True, drop=True)
    # Trunc2.to_csv(wd+'/'+basename+'_l10tofilter.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])

    print "Processing Levenshtein.distance<=4 Filtering between l=n-2 and l=n-1\n"

    start_time = time()

    Qname_rm = pd.DataFrame({'Qname':['begin']})
    q = 0
    total1 = Trunc2['Qname'].count()
    total2 = Full['Qname'].count()

    # count l = 12 with the same junction
    i = 0
    j = 0
    while i < total1:

        while  Trunc2.ix[i, 'B_junc'] != Full.ix[(j) ,'B_junc']:
            j = j + 1
            # print 'j {0}'.format(j)

        m = 0
        while Trunc2.ix[i, 'B_junc'] == Full.ix[(j+m) ,'B_junc']:
            m = m + 1
            if j+m > total2-1:
                # print 'j+m out of range'
                break

        # print "l=12: j {0} j+m-1 {1}".format(j,j+m-1)

        # count l = 11 with the same junction
        n = 0
        while Trunc2.ix[i, 'B_junc'] == Trunc2.ix[(i+n) ,'B_junc']:
            n = n + 1
            if i+n > total1-1:
                # print 'i+n out of range'
                break

        # print "l=11: i {0} i+n-1 {1}".format(i,i+n-1)

        end1=i+n
        end2=j+m
        for i in range(i,end1):
            for j in range(j,end2):
                if Levenshtein.distance(Trunc2['Barcode'].iloc[i],Full['Barcode'].iloc[j]) <=4 :
                    Qname_rm.loc[q] = Trunc2['Qname'].iloc[i]
                    # print "i {0} Qname {1} Barcode {2}".format(i,Trunc2['Qname'].iloc[i],Trunc2['Barcode'].iloc[i])
                    q=q+1
                    break
        i = i+1
        j = j+1
        # print "i {0} j {1}".format(i,j)

    print ("l = n-2 and l =n Barcode Levenshtein.distance<=4 Align Done in {}s!\n".format(round(time()-start_time, 3)))

    Trunc2 = Trunc2Raw[~Trunc2Raw['Qname'].isin(Qname_rm['Qname'])]
    Trunc2 = Trunc2.drop_duplicates(['Qname'])
    Trunc2 = Trunc2.sort_values(by=['B_junc','Barcode'])
    Trunc2.reset_index(level=Trunc2.index.names, inplace=True, drop=True)
    # Trunc2.to_csv(wd+'/'+basename+'_l10Uni_Result2.tlx', header = True, sep = '\t', columns = [u'Qname', u'B_junc', u'Barcode'])

    # get Filter step unique

    unique = Full.append([Trunc1])
    unique = unique.append([Trunc2])
    unique = unique.drop_duplicates(['Qname'])
    unique = unique.sort_values(by=['B_junc', 'Barcode'])
    unique.reset_index(level=Trunc2.index.names, inplace=True, drop=True)
    uniqueRaw = unique
    
    ###########################################################################
    # ~~ compare unique and Real barcode and drop Levenshtein.distance <=3 ~~ #
    ###########################################################################
    
    
    unique = unique[unique['B_junc'].isin(Real['B_junc'])]
    unique = unique.drop_duplicates(['Qname'])
    unique = unique.sort_values(by=['B_junc','Barcode'])
    unique.reset_index(level=unique.index.names, inplace=True, drop=True)
    
    print "Processing Levenshtein.distance<=3 Filtering between unique and Real\n"

    start_time = time()

    Qname_rm = pd.DataFrame({'Qname':['begin']})
    q = 0
    total1 = unique['Qname'].count()
    total2 = Real['Qname'].count() 
        
    # count l = 12 with the same junction 
    i = 0
    j = 0 
    while i < total1:
        
        while  unique.ix[i, 'B_junc'] != Real.ix[(j) ,'B_junc']:
            j = j + 1
            # print 'j {0}'.format(j)
            
        m = 0    
        while unique.ix[i, 'B_junc'] == Real.ix[(j+m) ,'B_junc']:
            m = m + 1           
            if j+m > total2-1:
                # print 'j+m out of range'
                break
                               
        # print "l=12: j {0} j+m-1 {1}".format(j,j+m-1)
        
        # count l = 11 with the same junction
        n = 0
        while unique.ix[i, 'B_junc'] == unique.ix[(i+n) ,'B_junc']:
            n = n + 1
            if i+n > total1-1:
                # print 'i+n out of range'
                break
        
        # print "l=11: i {0} i+n-1 {1}".format(i,i+n-1)
        
        end1=i+n
        end2=j+m
        for i in range(i,end1):
            for j in range(j,end2):
                if Levenshtein.distance(unique['Barcode'].iloc[i],Real['Barcode'].iloc[j]) <=3 :
                    Qname_rm.loc[q] = unique['Qname'].iloc[i]
                    # print "i {0} Qname {1} Barcode {2}".format(i,unique['Qname'].iloc[i],unique['Barcode'].iloc[i])
                    q=q+1
                    break
        i = i+1
        j = j+1
        # print "i {0} j {1}".format(i,j)

    print ("unique and Real Barcode Levenshtein.distance<=3 Align Done in {}s!\n".format(round(time()-start_time, 3)))
        
    unique = uniqueRaw[~uniqueRaw['Qname'].isin(Qname_rm['Qname'])]
    unique = unique.append([Real])
    unique = unique.drop_duplicates(['Qname'])
    unique = unique.drop_duplicates(['B_junc','Barcode'])
    unique = unique.sort_values(by=['B_junc','Barcode'])
    unique.reset_index(level=unique.index.names, inplace=True, drop=True)  
    duplicates = rawdata[~rawdata['Qname'].isin(unique['Qname'])]
    duplicates = duplicates.drop_duplicates(['Qname'])
    duplicates = duplicates.sort_values(by=['B_junc', 'Barcode'])
    duplicates.reset_index(level=duplicates.index.names, inplace=True)
    print 'Unique {0}\n'.format(unique['Qname'].count())
    print 'duplicates {0}\n'.format(duplicates['Qname'].count())
    
    # write result

    unique.to_csv(wd+'/'+basename+'_result.tlx', header = True, sep = '\t', index=False,
       columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
    u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
    u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
    u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
    u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
    u'repeatseq', u'duplicate', u'B_junc',u'P_len'])

    duplicates.to_csv(wd+'/'+basename+'_dup.tlx', header = True, sep = '\t', index=False,
       columns = [u'Qname', u'JuncID', u'Rname', u'Junction', u'Strand', u'Rstart',
    u'Rend', u'B_Rname', u'B_Rstart', u'B_Rend', u'B_Strand', u'B_Qstart',
    u'B_Qend', u'Qstart', u'Qend', u'Qlen', u'B_Cigar', u'Cigar', u'Seq',
    u'J_Seq', u'Barcode', u'unaligned', u'baitonly', u'uncut', u'misprimed',
    u'freqcut', u'largegap', u'mapqual', u'breaksite', u'sequential',
    u'repeatseq', u'duplicate', u'B_junc',u'P_len'])

    # B_junc frequency
    unique = pd.read_table(wd+'/'+basename+'_result.tlx',sep = '\t')
    duplicates = pd.read_table(wd+'/'+basename+'_dup.tlx',sep = '\t')

    uni = pd.value_counts(unique[u'B_junc']) 
    uni = uni.to_frame()
    uni.reset_index(level=uni.index.names, inplace=True)
    uni.columns = ['B_junc','uni']
    uni_fre = uni
    # uni_fre.to_csv(wd+'/'+basename+'_uni_fre.txt',header=True,sep='\t',index=False)

    dup = pd.value_counts(duplicates[u'B_junc'])
    dup = dup.to_frame()
    dup.reset_index(level=dup.index.names, inplace=True)
    dup.columns = ['B_junc','dup']
    dup_fre = dup
    # dup_fre.to_csv(wd+'/'+basename+'_dup_fre.txt',header=True,sep='\t',index=False)

    merge = pd.merge(uni_fre, dup_fre, on='B_junc', how='outer')
    merge = merge.sort_values(by=['B_junc'])
    merge.reset_index(level=merge.index.names, inplace=True)
    merge.to_csv(wd+'/'+basename+'_germfre.txt',header=True,sep='\t',index=False,columns = [u'B_junc',u'uni',u'dup'])

    # # plot
    #
    # print 'Begin plotting!'
    #
    # #rereadin merge file
    # merge = pd.read_table(wd+'/'+basename+'_germfre.txt',sep = '\t')
    # # Remove the plot frame lines. They are unnecessary chartjunk.
    # ax = plt.subplot(111)
    # ax.spines["top"].set_visible(False)
    # ax.spines["bottom"].set_visible(False)
    # ax.spines["right"].set_visible(False)
    # ax.spines["left"].set_visible(False)
    # # Ticks on the right and top of the plot are generally unnecessary chartjunk.
    # ax.get_xaxis().tick_bottom()
    # ax.get_yaxis().tick_left()
    # ax.set_yscale('log')
    # # Remove the tick marks; they are unnecessary with the tick lines we just plotted.
    # plt.tick_params(axis="both", which="both", bottom="off", top="off",
    #                 labelbottom="on", left="off", right="off", labelleft="on")
    #
    #
    #
    # merge.set_index(['B_junc'],inplace=True)
    # germfre = merge[['uni', 'dup']].plot(kind = 'bar', edgecolor = "none", title ="germline analysis", figsize=(15, 10), legend=True, fontsize=5, rot=45,
    #                                      color=[(44/255., 160/255., 44/255.), (214/255., 39/255., 40/255.)])
    #                                      #kind = 'bar', edgecolor = "none",
    #                                      #lw=1.5
    # germfre.set_xlabel('B_junc', fontsize=12)
    # germfre.set_ylabel("counts", fontsize=12)
    # germfre.get_figure().savefig(wd+'/'+basename+'_germfre.pdf', format='pdf')
    #
     
       
def main():
    args = parse_args()
    Uncut(args)

main() 
    