# <center>centromFind</center>
## Introducion
centromFind is a python program to finding centromeres from chromosomes based on Hi-C matrix. The result will be return and print in dictionary format. This program requires scipy, pandas and numpy.  

## Input
This program have 4 required parameters which can be obtained from [Hi-Cpro](https://github.com/nservant/HiC-Pro.git).  

```
usage: findCentre.py [-h] [--MINGAP MINGAP] [--SIGNAL_THRESHOLD SIGNAL_THRESHOLD] [-o O] matrix1 matrix2 bed1 bed2

A program for find centremere position based on Hi-C.

positional arguments:
  matrix1               Path to matrix of 100000 file
  matrix2               Path to matrix of 20000 file
  bed1                  Path to bed file of matrix1
  bed2                  Path to bed file of matrix2

optional arguments:
  -h, --help            show this help message and exit
  --MINGAP MINGAP       Minimum gap value n*100000 (default: 2)
  --SIGNAL_THRESHOLD SIGNAL_THRESHOLD
                        Signal threshold value (default: 0.7)
  -o O                  Output file
```  
We used 100k and 20k windows for chromosome finding, respectively, to locate candidate mitotic regions through the 100k window, and then used the 20k window to more accurately locate the coordinates of these regions. So user must modify the parameters of HiCpro to generate the 100k matrix, bed file and 20k matrix, bed file. The usage of [Hi-Cpro](https://github.com/nservant/HiC-Pro.git) can be visited in [there](https://github.com/nservant/HiC-Pro.git). We provide a test sample and config of HiC-pro in "sample" path for the user's information.  
If you hava other methods to obtained a HiC matrix, such as HiCexplorer, juicebox or your custom script,  you can still enter files(matrix, bed) conforming to input specifications into this programme specification. The format of input file as follow:  

```  
# matrix (this line excludes in the file)
# id1   id2 weight(count)       (this line excludes in the file)
1   1   34
1   2   15
2   2   30
2   3   5

# bed (this line excludes in the file)
# Chromosone    start   end id      (this line excludes in the file)
chr1    1   100000  1
chr1    100000  200000  2
chr2    1   100000  3
chr2    100000  200000  4
``` 

## Example  
We give user a complete example in this section. We suppose that there is a chromosome level reference sequence and a Hi-C reads. The annotation file  must be enter config file with absolute path.  

```  
# Preparing annotation files and index file
# build index for reference. This index will be used in HiCPro reference
bowtie2-build ref.fna ref  

# Obtain the restriction endonuclease sites
HiC-Pro_PATH/bin/utils/digest_genome.py -r ^GATC -r ^AATT ref.fna -o workpath/enzyme.bed 

# Obtain the chromosome list and their length. We supply a script for this step.
seq_list.py ref.fna -o ref.list

```

After prepared annotation file, We can run hicpro with the config as follow. In the config, The ***"MIN_MAPQ"*** must be set 30 or bigger. ***BOWTIE2_IDX_PATH*** is path of bowtie-build index exclude index file itself. ***BIN_SIZE*** must be set as 20000 100000 like sample. 
```
# config-hicpro.txt
#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 30

BOWTIE2_IDX_PATH = /data/yangjinbao/projects/centromereFind/index/  # This is path of index of bowtie-build. Only path, but not index.
BOWTIE2_GLOBAL_OPTIONS = -p 32 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS = -p 32 --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = ref  # This is name of index of bowtie-build 
GENOME_SIZE = ALL_PATH/ref.list  # This only use absolute path of chromosome and length list

#######################################################################
## Digestion Hi-C                                                                                                       #######################################################################
GENOME_FRAGMENT = ALL_PATH/enzyme.bed #  The absolute path of bed the enzyme restriction site
LIGATION_SITE = GATCGATC  

#######################################################################                                                 
## Contact Maps                                                                                                         #######################################################################
BIN_SIZE = 20000 100000
MATRIX_FORMAT = upper
```

Then run HiC-Pro. HiC-Pro must be use absolute path to run and the path of hic 
reads must be in the secondary directory of the input folder. The sample is as follow:  

```
ALL_PATH_HIC-Pro/bin/HiC-Pro -i hic_path -c config-hicpro.txt -o result_path 
# the  hic_path diractory tree structure is as follow:
./hic_path
--hic
----read_hic_R1.fastq
----read_hic_R2.fastq
```

At last we get the matrix and bed file after running HiC-Pro.  
matrix1 : result_path/hic_results/matrix/hic/raw/100000/hic_100000.matrix  
bed1 : result_path/hic_results/matrix/hic/raw/100000/hic_100000_abs.bed  
matrix2 : result_path/hic_results/matrix/hic/raw/20000/hic_20000.matrix  
bed2 : result_path/hic_results/matrix/hic/raw/20000/hic_20000_abs.bed  

```
findcentre.py matrix1 matrix2 bed1 bed2
```

Thanks your reading

## Cite
If you use this program, please cite the atricle of [https://www.sciencedirect.com/science/article/pii/S2590346224003870](https://www.sciencedirect.com/science/article/pii/S2590346224003870)
