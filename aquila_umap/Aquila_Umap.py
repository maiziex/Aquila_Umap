#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pysam
import time
import pickle
import os
from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing import Pool,cpu_count,active_children,Manager


# In[3]:


parser = ArgumentParser(description="Author: xzhou15@cs.stanford.edu liuyichen@std.uestc.edu.cn\n",usage='use "python3 %(prog)s --help" for more information')
parser.add_argument('--fa_folder','-fa',help="Required parameter; The folder path where fasta files are saved eg: /path/to/fasta/",required=True)
parser.add_argument('--fa_name','-fan',help="Required parameter; The file name of fasta file eg: sample.fa",required=True)
parser.add_argument('--out_dir','-o',help="Required parameter; The output path eg: /path/to/result/",required=True)
parser.add_argument('--chr_start','-start',type=int,help="chromosome start by,default = 1", default=1)
parser.add_argument('--chr_end','-end',type=int,help="chromosome end by,default = 23", default=23)
parser.add_argument('--kmer_len','-k', type=int,help="The length of generated kmers,default = 100", default=100)
parser.add_argument('--mapq_thres','-mapq', type=int,help="The MAPQ threshold to filter bowtie2 map result, default = 20", default=20)
parser.add_argument('--chr_thread','-cthread',type=int,help="number of threads for processing chromosome, default = 8 (recommended)", default=8)
parser.add_argument('--bowtie_thread','-bthread',type=int,help="number of threads for bowtie2 mapping, default = 20", default=20)

args = parser.parse_args()


# In[28]:


def QualStr (kmer_len):
    #Generate a fake quality string, default k = 100
    qual = ''
    for i in range(kmer_len):
        qual = qual+'~'
    return qual


# In[29]:


def GetGenome (fa_folder,fa_name):
    genome ={}# key = chrnum ; item = chrseq
    chro = []
    with open(fa_folder+fa_name) as f:
        for line in f:
            if line.startswith('>'):
                if chro:
                    genome[chrnum] = ''.join(chro)
                    chro = []
                    fw.close()
                chrnum = line[1:].split(' ')[0].strip('\n')
                fw = open(fa_folder+chrnum+'.fa','w')
                fw.write(line)
            else:
                chro.append(line.strip('\n'))
                fw.write(line)
        genome[chrnum] = ''.join(chro)
        chro = []
    return genome # key : chrnum (chr21) value : chrseq


# In[30]:


def GenerateFqFile (chrseq,chrnum,qual,kmer_len):
    with open(chrnum+'/'+chrnum+'.fq','w') as fw:
        for n in range(len(chrseq)+1-kmer_len):
            seq = chrseq[n:n+kmer_len]
            if 'N' not in seq:
                fw.write('@%s\n%s\n+\n%s\n'%(n+1,seq,qual))


# In[ ]:


def Bowtie2Map (chrnum,fa_folder,bowtie_thread):
    thread = str(bowtie_thread)
    os.chdir(chrnum+'/')
    os.mkdir('./ref'+chrnum)
    os.chdir('./ref'+chrnum)
    #index folder is refchrnum eg:refchr21
        
    index_build = os.popen('bowtie2-build --threads '+thread +' '+fa_folder+chrnum+'.fa'+' ref'+chrnum,'r')
    print(index_build.read())
        
    map_result = os.popen('bowtie2 -p '+thread+' -x '+'ref'+chrnum+' -U '+'../'+chrnum+'.fq '+'-S ../'+chrnum+'.sam','r')
    print(map_result.read())
    #output is chrnum.sam eg:chr21.sam
        
    os.chdir('..')
    out1 = os.popen('samtools view -bS '+chrnum+'.sam'+' > '+chrnum+'.bam')
    print(out1.read())
    out2 = os.popen('samtools sort '+chrnum+'.bam -o '+chrnum+'.sorted.bam')
    print(out2.read())
    out3 = os.popen('samtools view -H '+chrnum+'.sorted.bam > header.sam')
    print(out3.read())
    out4 = os.popen('samtools view -F 4 '+chrnum+'.sorted.bam | grep -v "XS:" | cat header.sam - | samtools view -b - > unique'+chrnum+'.bam')
    print(out4.read())
    out5 = os.popen('rm header.sam')
    print(out5.read())
    out6 = os.popen('samtools index unique'+chrnum+'.bam')
    print(out6.read())
    os.chdir('..')


# In[32]:


def Filter(chrnum,mapq_thres,kmer_len):
    filtered = []
    bamfile = pysam.AlignmentFile(chrnum+'/'+'unique'+chrnum+'.bam',"rb")
    for read in bamfile:
        if int(read.mapq) >= mapq_thres and not (read.is_reverse):
            filtered.append([int(read.pos),int(read.pos)+kmer_len-1])
    return filtered


# In[33]:


def Merge(chrnum,filtered):
    start = 0
    with open(chrnum+'/'+'merged.bed',"w") as fm:
        for line in filtered:
            if start == 0:
                start,end = line[0],line[1]
            elif line[0] > end+1:
                fm.write("%s\t%s\t%s\n"%(chrnum,start,end))
                start,end = line[0],line[1]
            else:
                end = line[1] 
        fm.write("%s\t%s\t%s\n"%(chrnum,start,end))


# In[34]:


def Get_uniqness(chrnum):
    uniq_map = defaultdict(int)
    with open(chrnum+'/'+"merged.bed","r") as f:
        with open(chrnum+'/'+"500merged.bed","w") as fw:
            for line in f:
                data = line.rsplit()
                _start = int(data[1])
                _end = int(data[2])
                block_len = _end - _start
                if block_len >= 500:
                    use_start = _start + 10
                    use_end = _end - 10
                    fw.write('%s\t%s\t%s\n'%(chrnum,use_start, use_end))
                    for step in range(use_start, use_end+1):
                        uniq_map[step] = 1
    pickle.dump(uniq_map,open("Uniqness_map/"+"uniq_map_"+chrnum+".p","wb"))


# In[35]:


def run(fa_folder,out_dir,chrseq,chrnum,kmer_len,qual,mapq_thres,bowtie_thread,xin):
    print("Starting to process "+chrnum)
    t = time.time()
    os.mkdir(chrnum)
    
    GenerateFqFile (chrseq,chrnum,qual,kmer_len)
    print(chrnum,":Generate .fq DONE")
    Bowtie2Map (chrnum,fa_folder,bowtie_thread)
    print(chrnum,":Bowtie2 mapping DONE")
    filtered = Filter(chrnum,mapq_thres,kmer_len)
    print(chrnum,":MAPQ filter DONE")
    Merge(chrnum,filtered)
    print(chrnum,":Merge DONE")
    Get_uniqness(chrnum)
    print(chrnum,":Get uniqness DONE")
    #-----------------------------------------------------------------------------------------------
    out7 = os.popen('rm -R '+chrnum+'/')#These two lines delete the intermediate results. 
    print(out7.read())                   #If you want to keep those results, comment out these two.
    #-----------------------------------------------------------------------------------------------
    print(chrnum,"DONE! Time used:", time.time()-t)


# In[1]:


def main()::
    fa_folder = args.fa_folder+"/"
    fa_name = args.fa_name
    out_dir = args.out_dir+"/"
    chr_start = args.chr_start - 1
    chr_end = args.chr_end
    kmer_len = args.kmer_len
    mapq_thres = args.mapq_thres
    chr_thread = args.chr_thread
    bowtie_thread = args.bowtie_thread
    
    #==============================================================================================
    
    qual = QualStr(kmer_len)
    Genome = GetGenome(fa_folder,fa_name)
    chr_list = list(Genome.keys())
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    os.chdir(out_dir)
    os.mkdir('Uniqness_map')
    
    count = 1
    pool = Pool(chr_thread)
    for chrnum in chr_list[chr_start:chr_end]:
        chrseq = Genome[chrnum]
        count += 1
        pool.apply_async(run,(fa_folder,out_dir,chrseq,chrnum,kmer_len,qual,mapq_thres,bowtie_thread,"xin"))
        if (count-1)%chr_thread == 0 or count-1 == (chr_end - chr_start):
            pool.close()
            while len(active_children()) > 1:
                time.sleep(0.5)
            pool.join()
            if (count - 1) == (chr_end - chr_start):
                print("finished!")
            else:
                pool = Pool(chr_thread)

    for chrnum in chr_list:
        os.popen('rm '+fa_folder+chrnum+'.fa')

    print("ALL DONE")


if __name__ == "__main__":
	main()
    

