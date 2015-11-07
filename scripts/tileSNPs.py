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
import math
from collections import defaultdict

from operator import itemgetter

__version__ = "1.0"

def main():

    parser=argparse.ArgumentParser(description='vcf2fasta (diploid)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-r', '--ref', dest='ref_fasta', type=str, required=True, help='input reference file (fasta)')
    parser.add_argument('-v', '--vcf', dest='vcf_file', type=str, required=True, help='input vcf file (vcf)')
    parser.add_argument('-n', '--name', dest='name', type=str, required=True, help='sample name (column header)')
    parser.add_argument('-rl', '--readlength', dest='read_length', type=int, default=50, help='length of reads, tiles all through SNP')
    parser.add_argument('--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)

    args=parser.parse_args()

    ref_fasta=args.ref_fasta
    vcf_file=args.vcf_file
    name=args.name
    read_length=args.read_length
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    # process VCF file
    
    verboseprint("")
    
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
    
    verboseprint("")
    
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
    out_fastq=output_wrapper(ref_fasta_name+'__'+name+'.fastq.gz')
    
    start=1
    end=0
    n_processed_snps=0
    pc_resolution=10000
    
    last_chrom=None
    
    buf=[]
    line_width=50 
    offset=0
    max_buf=math.ceil(((read_length*2)+1)/line_width)
    
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
                offset=0
                continue
            else: # random > found in line - issue with cat ?
                sys.exit('error with fasta file'+'\n'+str(line))
                
        if(chrom != last_chrom):
            
            if(last_chrom != None):
                verboseprint(last_chrom,n_processed_snps)
                verboseprint("")
            
            tmp_snps=[]
            tmp_snps=[(None,None)]            
            
            n_possible_snps=0
            if(chrom in snps):
                tmp_snps=snps[chrom]
                n_possible_snps=len(tmp_snps)
                pc_resolution=int(math.ceil(n_possible_snps/1000))
                verboseprint("\t",chrom,sep="")
        
        line=line.upper()
        line_start=end+1
        line_end=end+len(line)
        
        # change ref, embed SNPs
        while offset < len(tmp_snps) and tmp_snps[offset][0] != None and ( ( int(tmp_snps[offset][0] >= line_start) and int(tmp_snps[offset][0] <= line_end) ) or int(tmp_snps[offset][0] < line_start) ):
            snp_offset=int(tmp_snps[offset][0])-line_start
            
            if snp_offset < 0:
                offset += 1
                continue
            
            line = line[:snp_offset] + tmp_snps[offset][1].upper() + line[snp_offset+1:]
            offset += 1

        # build buffer to tile arounds SNPs
        if len(buf) < max_buf:
        
            buf.append(line)
        
        else:
                       
            # find all SNP positions
            
            while( len(tmp_snps) != 0 and int(tmp_snps[0][0] >= start+read_length) and int(tmp_snps[0][0] <= end-read_length) ):
                if n_processed_snps % pc_resolution == 0:
                    pc=(n_processed_snps/(n_possible_snps))*100
                    verboseprint("\r",""*50,"\r\t\tprocessing "+str(n_processed_snps)+" / "+str(n_possible_snps)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                    if verbose: sys.stdout.flush()
                
                snp_offset=int(tmp_snps[0][0])-start
                seq=''.join(buf)
                
                # ensure SNP is inserted!
                if seq[snp_offset] != tmp_snps[0][1].upper():
                    sys.exit('\n\nerror with SNP pos('+str(tmp_snps[0][0])+') snp['+str(tmp_snps[0][1])+'] vs ref['+seq[snp_offset]+'] in file'+'\n'+str(line))                

                snp_offset_start=snp_offset-read_length+1
                snp_offset_end=snp_offset+1
                for i in xrange(snp_offset_start,snp_offset_end):
                    
                    fastq_seq=seq[i:i+read_length]
                    seq_start=i+start
                    seq_end=seq_start+len(fastq_seq)-1
                    tmp_snp_offset=snp_offset-i
                    
                    fastq_header="@tileSNP:"+str(n_processed_snps)+":"+name+":"+str(tmp_snps[0][0])+":"+str(tmp_snps[0][1])+":"+str(tmp_snp_offset)+":"+last_chrom+":"+str(seq_start)+":"+str(seq_end)
                    fastq_qual='J'*len(fastq_seq)
                    
                    print(fastq_header,"\n",fastq_seq,"\n","+","\n",fastq_qual,sep="",file=out_fastq)
                
                n_processed_snps += 1                
                tmp_snps.pop(0)
                offset=0
               
            tmp_buf=buf[1:]
            tmp_buf.append(line)
            buf=tmp_buf
            start += len(line)
             
        end += len(line)
        last_chrom=chrom
                
    verboseprint(last_chrom,n_processed_snps)
        
    ref_fh.close()
    
def input_wrapper(infile):
    if infile.endswith('.gz'):
        fh=gzip.open(infile,'r')
    else:
        fh=open(infile,'r')
        
    return fh
    
def output_wrapper(outfile,append=False,suppress_comments=False):
    
    if outfile.endswith('.gz'):
        if append:
            fh=gzip.open(outfile,'a')
        else:
            fh=gzip.open(outfile,'w')   
    else:
        if append:
            fh=open(outfile,'a')
        else:
            fh=open(outfile,'w')
    
    # disable comment(s)if (UCSC format file)
    if outfile.endswith('.bed'):
        suppress_comments = True
    if outfile.endswith('.bed.gz'):
        suppress_comments = True
    if outfile.endswith('.bedGraph'):
        suppress_comments = True
    if outfile.endswith('.bedGraph.gz'):
        suppress_comments = True
    if outfile.endswith('.wig'):
        suppress_comments = True
    if outfile.endswith('.wig.gz'):
        suppress_comments = True
    if outfile.endswith('.sam'):
        suppress_comments = True
    if outfile.endswith('.sam.gz'):
        suppress_comments = True
    if outfile.endswith('.bam'):
        suppress_comments = True
    if outfile.endswith('.fastq.gz'):
        suppress_comments = True

    if not suppress_comments:
        print("## ",os.path.basename(__file__),sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Dekker Lab",sep="",file=fh)
        print("## Contact: Bryan R. Lajoie",sep="",file=fh)
        print("## https://github.com/blajoie",sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Version:\t",__version__,sep="",file=fh)
        print("## Date:\t",get_date(),sep="",file=fh)
        print("## Host:\t",get_compute_resource(),sep="",file=fh)
    
    return(fh)

def get_date():
    time=datetime.now()
    date=time.strftime('%I:%M:%S %p, %m/%d/%Y')
    
    return date
    
if __name__=="__main__":
    main()
