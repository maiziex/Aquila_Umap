# Aquila_Umap
```
aquila_umap --fa_folder /path/to/fasta/ --fa_name sample.fa  --out_dir /path/to/result/ --chr_start 1 --chr_end 1 --chr_thread 4
```
#### *Required parameters
##### --fa_folder: the directory which you save the reference fa files. 
##### --fa_name: the name of the reference fa file. 
#####  --out_dir: the directory to save the output results.

#### *Optional parameters
##### --kmer_len , default = 100, the length of generated kmers.
##### --mapq_thres: default = 50000 (50kb), the MAPQ threshold to keep unique mapping kmers. 
##### ---chr_thread , default = 2, number of threads for processing chromosome,
##### --chr_start, --chr_end: if you only process some chromosomes or only one chromosome. For example: use "--chr_start 1 --chr_end 5"  will process chromsomes 1,2,3,4,5. Use "--chr_start 2 --chr_end 2" will only assemlby chromosome 2. 


#### A practical example

Download Rhesus macaque (Macaca mulatta, ftp://ftp.ensembl.org/pub/release-97/fasta/macaca_mulatta/dna/ ) genome fasta file from  <a href="http://xinzhouneuroscience.org/wp-content/uploads/2019/08/macaca_mulatta.fa">HERE</a>.

Your folder structure should be as follows :
```
Rhesus_macaque
    |-fasta
        |-macaca_mulatta.fa
```
"macaca_mulatta.fa" is the fasta file you have just downloaded.

Then edit the second block and run the whole notebook.
```
fa_folder = "path/to/Rhesus_macaque/fasta/"
fa_name = "macaca_mulatta.fa"
out_dir = "path/to/Rhesus_macaque/output/"
start = (your start chromosome)
end = (your end chromosome)
kmer_len = 100 
mapq_thres = 20
bowtie_thread = 20 
```

When finished, you will get: 
```
Rhesus_macaque
|-fasta
|   |-macaca_mulatta.fa
|
|-output
    |-Uniqness_map
        |-uniq_map_chrxxx.p
        |-uniq_map_chrxxx.p
        |-uniq_map_chrxxx.p
        ...
```
The final results are stored in Uniqness_map folder.
