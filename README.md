# StrainScan
One efficient and accurate strain-level composition analysis tool based on reference genomes and k-mers.


### Contributor: Liao Herui and Ji Yongxin (Ph.D of City University of Hong Kong, EE)
### E-mail: heruiliao2-c@my.cityu.edu.hk / yxjijms@gmail.com
### Version: V1.0

---------------------------------------------------------------------------
### Dependencies:
* Python ==3.7.x
* R
* Sibeliaz ==1.2.2 (https://github.com/medvedevgroup/SibeliaZ)
* Required python package: numpy==1.17.3, pandas==1.0.1, biopython==1.74, scipy==1.3.1, sklearn==0.23.1, bidict==0.21.3, treelib==1.6.1

Make sure these programs have been installed before using StrainScan. 

## Install (Linux or ubuntu only)

Currently, yon can install StrainScan via [Anaconda](https://anaconda.org/) using the commands below:<BR/>
####
`git clone https://github.com/liaoherui/StrainScan.git`<BR/>
`cd StrainScan`<BR/>

`conda env create -f environment.yaml`<BR/>
`conda activate strainscan`<BR/>

`chmod 755 library/jellyfish-linux`<BR/>
`chmod 755 library/dashing_s128`<BR/>
####

Or, you can install all dependencies of StrainScan mannually and then run the commands below.

`git clone https://github.com/liaoherui/StrainScan.git`<BR/>
`cd StrainScan`<BR/>

`chmod 755 library/jellyfish-linux`<BR/>
`chmod 755 library/dashing_s128`<BR/>

## Pre-built databases download
The table below offers information about the pre-built databases of 6 bacterial species used in the paper. Users can download these databases and use them to identify strain directly.

Species   |	Source  | Number of Strains |	Number of Clusters |	Download link
------------ | -------------| ------------- | ------------- | ------------- 
A. muciniphila |  NCBI | 157 | 53  | [Google drive](https://drive.google.com/file/d/1Io096TYo_9vyvqOZi_fFTMKQi5sjA9vQ/view?usp=sharing)
C. acnes |  NCBI | 275 | 28  | [Google drive](https://drive.google.com/file/d/1Tvoz7qsfwBlj7gbRoW9ejIGYhszRua6b/view?usp=sharing)
P. copri |  NCBI | 112 | 51  | [Google drive](https://drive.google.com/file/d/1_ep-rSX0-bzvECpyuo-viuO7nXb1dkyz/view?usp=sharing)
E. coli |  NCBI | 1433 | 823  | [Google drive](https://drive.google.com/file/d/14c0rqVJNlhYb4T2JdhuVeeMsG_E0fsS8/view?usp=sharing)
M. tuberculosis |  NCBI | 792 | 25  | [Google drive](https://drive.google.com/file/d/1529TjKmPrG_kqTnE2NJcXjm0Bm04c-ww/view?usp=sharing)
S. epidermidis |  NCBI | 995 | 378  | [Google drive](https://drive.google.com/file/d/1fwGTe5WuzmJOYn1bsc8UzsIGqsmMyZb2/view?usp=sharing)




## Usage
One example about database construction and identification commands can be found in "<b>test_run.sh</b>".

### Use StrainScan to build your own custom database.<BR/>
  `python StrainScan_build.py -i <Input_Genomes> -d <Database_Dir>`<BR/>

### Use StrainScan to identify bacterial strains in short reads.
  `python StrainScan.py -i <Input_reads> -d <Database_Dir> -o <Output_Dir>`<BR/>
 
### Full command-line options
<!---(Note: The initial idea of development of StrainScan is "Simpler is better". We do not want to burden users due to complicated usage of StrainScan. So the default parameters (some are inside the program) are simple but have good performance in our test, however, more useful parameters will be added for users who need them.)-->

Identification - StrainScan.py (Default k-mer size: 31)
```
StrainScan - A kmer-based strain-level identification tool.

Example: python StrainScan.py -i Sim_Data/GCF_003812785.fq -d DB_Small -o Test_Sim/GCF_003812785

required arguments:
    -i, --input_fastq             Input fastq data.
    -d, --database_dir            Path of StrainScan database.

optional arguments:
    -h, --help                    Show help message and exit.
    -o, --output_dir              The output directory. (Default: ./StrainScan_Result)
    -k, --kmer_size               The size of k-mer, should be odd number. (Default: k=31)
    -l, --low_dep                 This parameter can be set to "1" if the sequencing depth of input data is very low (e.g. < 5x). For super low depth ( < 1x ), you can use "-l 2" (default: -l 0)
    -s, --minimum_snv_num         The minimum number of SNV at Layer-2 identification. (Default: s=40)
```
Build database - StrainScan_build.py (Default k-mer size: 31)
```
StrainScan - A kmer-based strain-level identification tool.

Example:  python StrainScan_build.py -i <Input_genomes> -o <Database_Dir>

required arguments:
     -i, --input_fasta             The path of input genomes. ("fasta" format)
     
optional arguments:
     -o, --output_dir              The output directory of constructed database. (Default: ./StrainScan_DB)
     -k, --kmer_size               The size of k-mer, should be odd number. (Default: k=31)
     -u, --uk_num                  The maximum number of unique k-mers in each genome to extract. (Default: 100000)
     -g, --gk_ratio                The ratio of group-specific k-mers to extract. (Default: g=1.0)        
     -m, --strainest_sample        If this parameter is 1, then the program will search joint kmer sets from msa generated by Strainest. To use this parameter, you have to make sure Strainest can run normally. (Default: 0)
     -n, --mink_cutoff             Minimum k-mer number cutoff in a node of the cluster search tree (CST). (Default: n=1000)
     -x, --maxk_cutoff             Maximum k-mer number cutoff in a node of the cluster search tree (CST). (Default: x=30000)
     -r, --maxn_cutoff             Maximum cluster number for node reconstruction of the cluster search tree(CST). (Default: r=3000)
```

## Output Format
The output of StrainScan contains two parts. The first part is the final identification report file in text format. This file contains all identified strains and their predicted depth and relative abundance, etc. The second part is the strain identification report files inside each cluster.

For your reference, two output files are given as example in the folder "Output_Example" in this repository. These files contain identification results of one single-strain and one two-strain (depth: 5X and 5X) simulated datasets, respectively.

Explaination about the headers in the final identification report file (E.g. "Output_Example/GCA_000144385_5X_GCF_008868325_5X/final_report.txt") of StrainScan.
Header    |	Description
------------ | ------------- 
Strain_ID | The numerical id of identified strains in the ascending order.
Strain_Name | The name of identified strains. (In the example output, the name refers to the NCBI RefSeq accession id)
Cluster_ID  | The cluster id of identified strains. (For cluster information, users can check "<Database_Dir>/Cluster_Result/hclsMap_95_recls.txt")
Relative_Abundance_Inside_Cluster | The predicted relative abundance of identified strains inside the cluster.
Predicted_Depth (Enet) | The predicted sequencing depth of identified strains inside the cluster using elastic net model.
Predicted_Depth (Ab\*cls_depth) | The final predicted sequencing depth of identified strains.
Coverage  | The estimated k-mer-based coverage of identified strains.
Coverd/Total_kmr  | The number of "covered" and "total" k-mers of identified strains.
Valid_kmr | The valid k-mer refers to the k-mer belonging to the identified strain during the iterative matrix multiplication. More valid k-mers there are, more likely this strain exist.
Remain_Coverage | The coverage calculated by "covered" / "total" k-mers during the iterative matrix multiplication.
CV  | The number of "covered" and "valid" k-mers of identified strains.
Exist evidence  | By default, identified strains with "relative abundance > 0.02 and coverage >0.7" will be marked as "\*". Strains with "\*" are more likely to exist. However, for low-depth strains, this parameter is not useful.
