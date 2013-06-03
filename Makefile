data/genomes/NC_00913.fna : 
	curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna > data/genomes/NC_00913.fna

data/genomes : data/genomes/NC_00913.fna
	bowtie-build e_coli_k12_mg1655 data/genomes/NC_00913.fna

results/AraC/maps/AraCa.map : data/genomes/e_coli_k12_mg1655* data/Sample_AMS-ara-gal-Ec_AraCa/AMS-ara-gal-Ec_AraCa_CAGATC_L006_R1_001.fastq
	bowtie -m 1 data/genomes/e_coli_k12_mg1655 data/Sample_AMS-ara-gal-Ec_AraCa/AMS-ara-gal-Ec_AraCa_CAGATC_L006_R1_001.fastq > results/AraC/maps/AraCa.map

results/AraC/maps/AraCb.map : data/genomes/e_coli_k12_mg1655* data/Sample_AMS-ara-gal-Ec_AraCb/AMS-ara-gal-Ec_AraCb_CTTGTA_L006_R1_001.fastq
	bowtie -m 1 data/genomes/e_coli_k12_mg1655 data/Sample_AMS-ara-gal-Ec_AraCb/AMS-ara-gal-Ec_AraCb_CTTGTA_L006_R1_001.fastq > results/AraC/maps/AraCb.map 

results/AraC/AraCa.csv : src/map2csv.py data/Sample_AMS-ara-gal-Ec_AraCa/AraCa.map
	python src/map2csv.py results/AraC/ results/AraCa.csv
results/AraC/AraCb.csv : src/map2csv.py data/Sample_AMS-ara-gal-Ec_AraCb/AraCb.map
	python src/map2csv.py data/Sample_AMS-ara-gal-Ec_AraCb/AraCb.map results/AraCb.csv