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
