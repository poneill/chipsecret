# Get the E. coli genome
data/genomes/NC_00913.fna : 
	curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna > data/genomes/NC_00913.fna

# Build a Bowtie index for it
data/genomes : data/genomes/NC_00913.fna
#	bowtie-build e_coli_k12_mg1655 data/genomes/NC_00913.fna
	bowtie-build data/genomes/NC_00913.fna e_coli_k12_mg1655 
# Map reads from AraCa replicate
results/AraC/maps/AraCa.map : data/genomes/e_coli_k12_mg1655* data/Sample_AMS-ara-gal-Ec_AraCa/AMS-ara-gal-Ec_AraCa_CAGATC_L006_R1_001.fastq
	bowtie -m 1 data/genomes/e_coli_k12_mg1655 data/Sample_AMS-ara-gal-Ec_AraCa/AMS-ara-gal-Ec_AraCa_CAGATC_L006_R1_001.fastq 1> results/AraC/maps/AraCa.map 2> results/AraC/maps/AraCa.out

# Map reads from AraCb replicate
results/AraC/maps/AraCb.map : data/genomes/e_coli_k12_mg1655* data/Sample_AMS-ara-gal-Ec_AraCb/AMS-ara-gal-Ec_AraCb_CTTGTA_L006_R1_001.fastq
	bowtie -m 1 data/genomes/e_coli_k12_mg1655 data/Sample_AMS-ara-gal-Ec_AraCb/AMS-ara-gal-Ec_AraCb_CTTGTA_L006_R1_001.fastq 1> results/AraC/maps/AraCb.map 2> results/AraC/maps/AraCb.out

# Compile map file into read density csv for AraCa replicate
results/AraC/AraCa.csv : src/map2csv.py results/AraC/maps/AraCa.map
	python src/map2csv.py results/AraC/maps/AraCa.map results/AraC/AraCa.csv

# Compile map file into read density csv for AraCa replicate
results/AraC/AraCb.csv : src/map2csv.py results/AraC/maps/AraCb.map
	python src/map2csv.py results/AraC/maps/AraCb.map results/AraC/AraCb.csv

all: results/AraC/AraCa.csv results/AraC/AraCb.csv
