# StrainScan
One efficient and accurate strain-level composition analysis tool based on reference genomes and k-mers.


### E-mail: heruiliao2-c@my.cityu.edu.hk
### Version: V1.0

---------------------------------------------------------------------------
### Dependencies:
* Python ==3.7.x
* R
* Sibeliaz ==1.2.2 (https://github.com/medvedevgroup/SibeliaZ)
* Required python package: numpy==1.17.3, pandas==1.0.1, biopython==1.74, scipy==1.3.1, sklearn==0.23.1, bidict==0.21.3, psutil==5.7.0

Make sure these programs have been installed before using StrainScan. 

## Install (Linux or ubuntu only)
Currently, yon can install StrainScan via Anaconda using the command below:<BR/>
`conda env create -f environment.yaml`<BR/>
`conda activate strainscan`<BR/>

Or, you can install StrainScan mannually.

####
`git clone https://github.com/liaoherui/StrainScan.git`<BR/>
`cd StrainScan`<BR/>
`chmod 755 library/jellyfish-linux`<BR/>
`chmod 755 library/dashing_s128`<BR/>
####

## Usage
One example about database construction and identification command can be found in "test_run.sh".

### Use StrainScan to build your own custom database.<BR/>
  `python StrainScan_build.py -i <Input_Genomes> -d <Database_Dir>`<BR/>

### Use StrainScan to identify bacterial strains in short reads.
  `python StrainScan.py -i <Input_reads> -d <Database_Dir> -o <Output_Dir>`<BR/>
