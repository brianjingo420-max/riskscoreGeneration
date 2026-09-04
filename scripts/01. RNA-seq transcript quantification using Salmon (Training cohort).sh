# Prepare Salmon indexing
# INDEX_PATH="/path/to/salmon_index"
# OUTPUT_DIR="/path/to/salmon_output/sample01"

# Ensembl release 114
FASTA=Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa

# Session
# Salmon v1.10.2

# Salmon indexing
salmon index -t "${FASTA}" -i "${INDEX_PATH}" -k 31

# Salmon quantification
salmon quant \
  -i "${INDEX_PATH}" \
  -l A \
  -1 "sample01_R1.fastq" \
  -2 "sample01_R2.fastq" \
  -p 8 \
  --validateMappings \
  -o "${OUTPUT_DIR}"

# Repeat the quantification command for each sample.
