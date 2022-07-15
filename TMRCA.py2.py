# -*- coding: utf-8 -*-
## python2.7

import argparse
import pandas
import numpy as np
import sys
import socket
import os
import gzip
import time
import subprocess
import itertools

def readregion(samplefile, regionfile, poplist, pairsfile):
    regions = pandas.read_csv(regionfile,sep='\s+',header=None,usecols=[0,1,2],names=['chr','start','end'])
    regions.drop_duplicates(inplace=True)

    if samplefile == 'all':
        pairs = pandas.DataFrame({'pop1':['pop'],'pop2':['pop']})
    else:
        if pairsfile == 'all':
            pairs = pandas.DataFrame([[p1,p2] for p1 in poplist for p2 in poplist])
            pairs.columns = ['pop1','pop2']
        else:
            pairs = pandas.read_csv(pairsfile,sep='\s+',header=None,usecols=[0,1],names=['pop1','pop2'])
            pairs = pairs[(pairs['pop1'].isin(poplist)) & (pairs['pop2'].isin(poplist))]
    
    regiondf = []
    for i in list(pairs.index):
        pop1, pop2 = list(pairs.loc[i, ['pop1','pop2']])
        tmp = regions.copy()
        tmp['pop1'] = pop1; tmp['pop2'] = pop2
        regiondf += [tmp]

    regiondf = pandas.concat(regiondf, ignore_index=True)
    regiondf['chr'] = regiondf['chr'].astype(str)

    return regiondf

def readsample(vcffile, samplefile):
    with gzip.open(vcffile) as f:
        line = f.readline()
        while line[:2] == '##':
            line = f.readline()
    samplelist = line.strip().split('\t')[9:]

    if samplefile == 'all':
        sampledf = pandas.DataFrame({0:[s+'_'+x for s in samplelist for x in ['1','2']]})
        sampledf[1] = 'pop'
    else:
        sampledf = pandas.read_csv(samplefile,sep='\s+',header=None,usecols=[0,1])
        samplelist = list(set(sampledf[0].apply(lambda x: x[:-2])) & set(samplelist))
        sampledf = sampledf[sampledf[0].isin([s+'_'+x for s in samplelist for x in ['1','2']])]
    
    poplist = list(sampledf[1].unique())

    return samplelist, sampledf, poplist

def readvcf(vcffile, regiondf):
    ## read vcf
    headerline = 0
    with gzip.open(vcffile) as f:
        line = f.readline()
        while line[:2] == '##':
            headerline += 1
            line = f.readline()

    samplelist = line.strip().split('\t')[9:]
    datatype = dict(zip(['POS']+samplelist, ['int32']+['category']*len(samplelist)))
    data = pandas.read_csv(vcffile,sep='\t',skiprows=range(headerline),usecols=['#CHROM','POS']+samplelist, dtype=datatype)
    data['#CHROM'] = data['#CHROM'].astype(str)

    chrlist = list(set(data['#CHROM']) & set(regiondf['chr']))
    data = data[data['#CHROM'].isin(chrlist)]
    region = regiondf[regiondf['chr'].isin(chrlist)][['chr','start','end']].copy()
    region.sort_values(['chr','start','end'],ascending=True,inplace=True)
    region.drop_duplicates(inplace=True)
    region.index = range(region.shape[0])

    ## extract vcf
    # merge region
    if region.shape[0] == 1:
        pass
    else:
        for index in list(region.index)[:-1]:
            chrom1, start1, end1 = list(region.loc[index])[:3]
            chrom2, start2, end2 = list(region.loc[index+1])[:3]

            if ((chrom2 == chrom1) & (start2 <= end1+1)):
                new_end = max(end1,end2)
                region.loc[index+1,'start'] = start1
                region.loc[index+1,'end'] = new_end
                region.drop(index,inplace=True)
            else:
                continue

    # extract vcf
    data = pandas.concat(list(region.apply(lambda x: data[(data['#CHROM']==x['chr']) & (data['POS']>=x['start']) & (data['POS']<=x['end'])].copy(), axis=1)),ignore_index=True)
    data.drop_duplicates(['#CHROM','POS'],inplace=True)

    ## convert format
    pos = data[['#CHROM','POS']].copy()

    #geno = data[samplelist].applymap(lambda x: x[:3].replace('|','').replace('/',''))    
    geno = pandas.DataFrame(data[samplelist].apply(lambda x: list('|'.join(x[:]).replace('|','')),axis=1).tolist(),columns=[s+'_'+x for s in samplelist for x in ['1','2']]).astype('category')

    data = pandas.concat([pos,geno],axis=1)

    return data

def readape(apefile,chrlist):
    apedf = {}
    for chrom in chrlist:
        apedf[chrom] = pandas.read_csv(apefile.replace('@',str(chrom)), sep='\s+',header=None)

    return apedf

## tmrca = m/(2uL), m=pi, N=1, L stand for the region length
## u=d/(L*T*2)
## d: divergence between reference genomes of human and ape
## T = 13 million years, time from commone ancestor of human and ape
def tmrca(data, divT):
    data['T'] = divT  # 13e6
    data['tmrca'] = data['pi'].values * divT / data['d'].values

    return data

def pop_pi(hapdata, sampledf, pop1, pop2):
    haplist1 = list(sampledf[sampledf[1]==pop1][0])
    df1 = hapdata[haplist1].copy()
    
    ## allele count
    df1count1 = (df1=='1').apply(sum,axis=1); df1count0 = (df1=='0').apply(sum,axis=1)
    hap1count = pandas.concat([df1count0, df1count1],axis=1).astype('int64')
    hap1count.rename(columns=lambda x: str(x),inplace=True)
    # ref allele count, alt allele count
    rac1 = hap1count['0'].values; aac1 = hap1count['1'].values
    count1 = rac1 + aac1

    ## pairwise difference
    if pop1 == pop2:
        pi = np.nansum(rac1 * aac1 * 2.0 / (count1*(count1-1.0)))
    else:
        haplist2 = list(sampledf[sampledf[1]==pop2][0])
        df2 = hapdata[haplist2].copy()

        df2count1 = (df2=='1').apply(sum,axis=1); df2count0 = (df2=='0').apply(sum,axis=1)
        hap2count = pandas.concat([df2count0, df2count1],axis=1).astype('int64')
        hap2count.rename(columns=lambda x: str(x),inplace=True)
        
        rac2 = hap2count['0'].values; aac2 = hap2count['1'].values
        count2 = rac2 + aac2

        ## pairwise difference
        pi = np.nansum((rac1*aac2 + rac2*aac1) * 1.0/(count1 * count2))
    
    return pi

def ind_pi_tmp(df):
    df.columns = [0,1]

    diff = (((df[0].values == '1') & (df[1].values == '0')) | ((df[0].values == '0') & (df[1].values == '1'))).sum()

    return diff

def ind_pi(hapdata, sampledf, chrom, start, end, pop1, pop2):
    sampleinfo = dict(zip(list(sampledf[0]),list(sampledf[1])))

    haplist1 = list(sampledf[sampledf[1]==pop1][0])
    haplist2 = list(sampledf[sampledf[1]==pop2][0])
    haplisttot = list(set(haplist1) | set(haplist2))

    annotdf = pandas.DataFrame([[h1,h2] for h1 in haplist1 for h2 in haplist2 if h1!=h2])
    if annotdf.empty:
        return annotdf
    else:
        pass
    annotdf.columns = ['hap1','hap2']

    annotdf['chr'] = chrom; annotdf['start'] = start; annotdf['end'] = end
    annotdf['pop1'] = annotdf['hap1'].apply(lambda x: sampleinfo[x])
    annotdf['pop2'] = annotdf['hap2'].apply(lambda x: sampleinfo[x])
    annotdf = annotdf[['chr','start','end','pop1','pop2','hap1','hap2']]

    ## pairwise difference
    hd = hapdata[haplisttot].copy()
    annotdf['pi'] = annotdf.apply(lambda x: ind_pi_tmp(hd[[x['hap1'], x['hap2']]].copy()), axis=1)

    return annotdf

def calculate_pi_d(vcfdata, apedf, samplelist, sampledf, poplist, regions, Tind, outprefix):    
    tmp = regions.copy()  # <chr> <start> <end> <pop1> <pop2>
    tmp['pi'] = np.nan; tmp['d'] = np.nan
    tmp['index'] = tmp.apply(lambda x: '{}:{}:{}'.format(x['chr'],x['start'],x['end']),axis=1)
    indexlist = list(tmp['index'].unique())

    ## dict, region:d
    r = tmp[['chr','start','end']].drop_duplicates()
    r.index = r.apply(lambda x: '{}:{}:{}'.format(x['chr'],x['start'],x['end']),axis=1)
    r['d'] = r.apply(lambda x: apedf[x['chr']][(apedf[x['chr']][1]>=x['start']) & (apedf[x['chr']][1]<=x['end'])].shape[0], axis=1)
    r = dict(zip(list(r.index),list(r['d'])))
    
    ## pi & d, pop
    annotpops = []
    for index in indexlist:
        annot = tmp[tmp['index']==index].copy()

        chrom, start, end = list(annot.iloc[0,:3])
        hapdata = vcfdata[(vcfdata['#CHROM']==chrom) & (vcfdata['POS']>=start) & (vcfdata['POS']<=end)].copy()

        annot['pi'] = annot.apply(lambda x: pop_pi(hapdata, sampledf, x['pop1'], x['pop2']), axis=1)
        annot['d'] = r[index]

        annotpops += [annot]
    
    annotpops = pandas.concat(annotpops, ignore_index=True)
    annotpops.drop(['index'],axis=1,inplace=True)
    annotpops.sort_values(by=['chr','start','end','pop1','pop2'],ascending=True,inplace=True)
    
    ## pi & d, ind
    if Tind == 'F':
        annotinds = pandas.DataFrame()
    else:
        annotinds = []
        for index in indexlist:
            annot = tmp[tmp['index']==index].copy()
            chrom, start, end = list(annot.iloc[0,:3])
            hapdata = vcfdata[(vcfdata['#CHROM']==chrom) & (vcfdata['POS']>=start) & (vcfdata['POS']<=end)].copy()
            
            annot = pandas.concat(list(annot.apply(lambda x: ind_pi(hapdata, sampledf, x['chr'], x['start'], x['end'], x['pop1'], x['pop2']), axis=1)))
            annot['d'] = r[index]

            annotinds += [annot]

        annotinds = pandas.concat(annotinds, ignore_index=True)
        annotinds.sort_values(by=['chr','start','end','pop1','pop2','hap1','hap2'],ascending=True,inplace=True)
    
    return annotpops, annotinds

def main():
    ## input
    parser = argparse.ArgumentParser()
    parser.add_argument("--gzvcf", type=str, required = True, \
                        help="/path/to/phased.vcf.gz file.")
    parser.add_argument("--samples", type=str, required=False, default='all', \
                        help="<sample_1/2> <pop>, which haplotype to be used, e.g., sampleX_2 popX, stands for the 2nd haplotype from sampleX. additional columns will be ignored, no header. if not given, all samples will be used, and considered as one pop.")
    parser.add_argument("--Tind", type=str, required=False, default='F',choices=['T','F'], \
                        help="whether to estimate pairwise TMRCA between individuals/haplotypes")
    parser.add_argument("--region", type=str, required = True, \
                        help="<chr> <start> <end>, no header line, tab or space sperated, additional columns will be ignored")
    parser.add_argument("--pairs", type=str, required = False, default='all', \
                        help="<pop1> <pop2>, pop name(s) should be in the --samples file. if not given, all possible pairs of populations will be considered. no header line, tab or space sperated, additional columns will be ignored.")
    parser.add_argument("--ape",type=str, required = True, \
                        help="humna vs. ape diff")
    parser.add_argument("--divT",type=float, required = False, default=13e6, \
                        help="divergence time between humna and ape (year)")
    parser.add_argument("--out", type=str, required = False, default='out', \
                        help="/out/file/prefix")
    args = parser.parse_args()

    ## log
    with open(args.out+'.logfile','w') as log:
        log.write('python {}\n'.format(sys.argv[0]))
        log.write('{}--gzvcf        {}\n'.format(' '*8, args.gzvcf))
        log.write('{}--samples      {}\n'.format(' '*8, args.samples))
        log.write('{}--Tind         {}\n'.format(' '*8, args.Tind))
        log.write('{}--region       {}\n'.format(' '*8, args.region))
        log.write('{}--pairs        {}\n'.format(' '*8, args.pairs))
        log.write('{}--ape          {}\n'.format(' '*8, args.ape))
        log.write('{}--divT         {}\n'.format(' '*8, args.divT))
        log.write('{}--out          {}\n\n'.format(' '*8, args.out))
        
        log.write('Hostname: '+socket.gethostname()+'\n')
        log.write('Working directory: '+os.getcwd()+'\n')
        log.write('Start time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    ## sample, region, and frequency
    # samplelist from overlaping samples in vcf and sample file
    # sampledf: <hap> <pop>
    # "pop" if no pop name is provided 
    samplelist, sampledf, poplist = readsample(args.gzvcf, args.samples)
    
    # <chr> <start> <end> <pop1> <pop2>, all possible pairs of populations
    regions = readregion(args.samples, args.region, poplist, args.pairs)
    if regions.empty:
        with open(args.out+'.logfile','w') as log:
            log.write('NO region remained.\n')
        exit()
    else:
        pass

    ## vcf data, <#CHROM> <POS> <hap...>
    vcfdata = readvcf(args.gzvcf, regions)

    ## prune vcf data and region
    chrlist = list(set(vcfdata['#CHROM']) & set(regions['chr']))
    regions = regions[regions['chr'].isin(chrlist)]

    ## human - ape divergence
    apedata = readape(args.ape, chrlist)

    ## annotdf: regions + <pi> <d>
    # if Tind == T, output inddf, regions + <ind1> <ind2> <pi> <d>
    annotdf, inddf = calculate_pi_d(vcfdata, apedata, samplelist, sampledf, poplist, regions, args.Tind, args.out)

    ## TMRCA
    res = tmrca(annotdf, args.divT) # 13e6
    res.to_csv(args.out+'.tmrca.txt',sep='\t',index=None)

    if args.Tind == 'T':
        resind = tmrca(inddf, args.divT)
        resind.to_csv(args.out+'.tmrca.ind.txt',sep='\t',index=None)
    else:
        pass

    with open(args.out+'.logfile','a') as log:
        log.write("Done TMRCA calculation.\n")
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    print('Have a Nice Day!')

if __name__ == '__main__':
    main()
