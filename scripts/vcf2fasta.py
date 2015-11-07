from __future__ import print_function
from __future__ import division

import numpy as np
import sys
import argparse
import time
import re
import gzip
import os
import logging
from collections import defaultdict

from operator import itemgetter

__version__ = "1.0"

def main():

    parser=argparse.ArgumentParser(description='vcf2fasta (diploid)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-r', '--ref', dest='ref_fasta', type=str, required=True, help='input reference file (fasta)')
    parser.add_argument('-v', '--vcf', dest='vcf_file', type=str, required=True, help='input vcf file (vcf)')
    parser.add_argument('-n', '--name', dest='name', type=str, required=True, help='sample name (column header)')
    parser.add_argument('--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    ref_fasta=args.ref_fasta
    vcf_file=args.vcf_file
    name=args.name
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    # process VCF file
    
    verboseprint("processing VCF file")
    
    snps=defaultdict(list)
    snp_count=0
    last_chrom=None
    chrs=dict()
    
    # open vcf file
    if vcf_file.endswith('.gz'):
        vcf_fh=gzip.open(vcf_file,'r')
    else:
        vcf_fh=open(vcf_file,'r')
    
    # iterate over vcf file
    for linenum,line in enumerate(vcf_fh):
        if line.startswith("##"): # skip headers
            continue
        
        line=line.rstrip("\n")
        
        # get header line
        if line.startswith("#"):
            header=line.lstrip("#").split("\t");
            header2index=dict([(h,i) for i,h in enumerate(header)])
        
            # ensure FORMAT column is included
            if "FORMAT" not in header2index:
                print("FORMAT field not specified in VCF file!")
                print(header2index)
                sys.exit('error')
        
            # ensure user-specified sample name column is included
            if name not in header2index:
                print(name,"field not specified in VCF file!")
                print(header2index)
                sys.exit('error')
            
            continue
        
        tmp=line.split("\t")
        genotype=tmp[header2index[name]]
        
        format=tmp[header2index["FORMAT"]].split(":")
        index2field=dict([(f,i) for i,f in enumerate(format)])            
        # ensure GT field id included in FORMAT column
        if "GT" not in index2field:
            print("GT field not specified in FORMAT!")
            print(index2field)
            sys.exit('error')
                
        genotype_list=genotype.split(":")
        gt=genotype_list[index2field["GT"]]
        
        pattern = re.compile('[\|\/]')
        (a,b) = pattern.split(gt)
        
        if(a != b):
            sys.exit('error: non-homo SNP found @ line# '+linenum+'\n')
           
        c=a=b
        c=int(c)
        
        chrom=tmp[header2index["CHROM"]]
        pos=int(tmp[header2index["POS"]])
        ref=tmp[header2index["REF"]]
        alt=tmp[header2index["ALT"]].split(",")
            
        if(c == 0):
            snps["chr"+chrom].append((pos,ref))
            snp_count+=1
        elif(c >= 1):
            snps["chr"+chrom].append((pos,alt[c-1]))
            snp_count+=1
        
        if(chrom != last_chrom):
            if(chrom not in chrs):
                verboseprint("\tchr",chrom)
                chrs[chrom]=1
                
        last_chrom=chrom
        
    vcf_fh.close()
    
    verboseprint("found",snp_count,"snps")
    
    # process VCF file
    
    verboseprint("")
    
    # ensure all snps are sorted by position (sorted seperatley for each chromosome)
    verboseprint("sorting by position")
    for chr in snps: # ensure sorted by pos
        snp_positions=snps[chr]
        verboseprint("\t",chr," ... ",len(snp_positions)," snps",sep="")
        sorted_snp_positions=sorted(snp_positions, key=itemgetter(0)) 
        snps[chr]=sorted_snp_positions
    verboseprint("")
    
    # process REFERENCE file
  
    verboseprint("processing REF file")
    
    # get output name
    ref_fasta_name=os.path.basename(ref_fasta)
    ref_fasta_name=re.sub(".gz", "", ref_fasta_name)
    ref_fasta_name=re.sub(".fasta", "", ref_fasta_name)
    ref_fasta_name=re.sub(".fa", "", ref_fasta_name)
    out_fh=open(ref_fasta_name+'__'+name+'.fa',"w")
    
    placed_snps=0
    total_placed_snps=0
    current_snp=(None,None)
    pos=1
    last_chrom=None
    tmp_pos_list=[(None,None)]
    
    # open reference fasta file
    if ref_fasta.endswith('.gz'):
        ref_fh=gzip.open(ref_fasta,'r')
    else:
        ref_fh=open(ref_fasta,'r')
    
    # iterate over fasta file
    for linenum,line in enumerate(ref_fh):
        line=line.rstrip("\n")
        
        # search for > (contig name)
        regexp = re.compile('>')
        if regexp.search(line):
            if line.startswith(">"):
                chrom=line.lstrip(">")
                pos=1
                print(line,"-",name,file=out_fh,sep="")
                continue
            else: # random > found in line - issue with cat ?
                sys.exit('error with fasta file'+'\n'+str(line))
        
        if(chrom != last_chrom):
            
            tmp_pos_list=[]
               
            if(last_chrom != None):
                verboseprint(" ... ",placed_snps," / ",possible_snps,sep="")
                
            tmp_pos_list=[(None,None)]            
            possible_snps=0
            if(chrom in snps):
                tmp_pos_list=snps[chrom]
                possible_snps=len(tmp_pos_list)
                
            verboseprint("\t",chrom,sep="",end="")
            current_snp=tmp_pos_list.pop(0)
            
            total_placed_snps += placed_snps
            placed_snps=0
           
        tmp_len=len(line)
        start=pos
        end=pos+tmp_len-1
        
        while((current_snp[0] < start) and (len(tmp_pos_list) > 0)):
            print("ERROR: missed snp!",current_snp,"\t",start,"-",end,">",current_snp[0])
            current_snp=tmp_pos_list.pop(0)
        
        if((current_snp[0] == None) or (current_snp[0] > end)):
            print(line,file=out_fh)
        else:
            char_list=list(line)

            snp_offset=current_snp[0]-start
            if((snp_offset < 0) or (snp_offset > len(char_list))): # check to ensure SNP overlaps interval
                sys.exit('error '+str(current_snp)+' '+str(snp_offset)+' '+str(start)+'-'+str(end))
                
            # replace snp in char arr
            char_list[snp_offset]=current_snp[1]
            placed_snps+=1
            
            if(len(tmp_pos_list) == 0):
                current_snp=(None,None)
                    
            # handle multiple SNPs per FASTA line (normally 50 chars /buffer/)
            if(len(tmp_pos_list) > 0):
                current_snp=tmp_pos_list.pop(0)
                while((current_snp[0] <= end) and (len(tmp_pos_list) > 0)):
                    snp_offset=current_snp[0]-start
                    
                    # replace snp in char arr
                    char_list[snp_offset]=current_snp[1]
                    placed_snps+=1
                    current_snp=tmp_pos_list.pop(0)
                
                if((current_snp[0] <= end) and (len(tmp_pos_list) == 0)):
                    snp_offset=current_snp[0]-start
                    char_list[snp_offset]=current_snp[1]
                    placed_snps+=1
                    current_snp=(None,None)
            
            # char list to string, and print
            print(''.join(char_list),file=out_fh)
            
        pos += tmp_len
        last_chrom=chrom
        
    ref_fh.close()
    
    # handle last line
    verboseprint(" ... ",last_chrom," ",placed_snps," / ",len(snps[last_chrom]),sep="")
           
    # process REFERENCE file
    
if __name__=="__main__":
    main()

    