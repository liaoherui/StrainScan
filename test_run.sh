python StrainScan_build.py -i Test_genomes -o DB_Small

python StrainScan.py -i Sim_Data/GCF_003812785.fq -d DB_Small -o Test_Sim/GCF_003812785

python StrainScan.py -i Sim_Data_mul/GCA_000144385_5X_GCF_008868325_5X.fq -d  DB_Small -o Test_Sim/GCA_000144385_5X_GCF_008868325_5X 
