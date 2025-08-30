CPU=20
. ~/enable_conda.rc
conda activate kfri_genome_pipeline
T=/tmp/flye_${1##*/}
C=/tmp/${1##*/}_contigs.fasta

flye --meta --nano-raw /home/caefs/microbiome/projects/metagenome_analysis/KMAGs/assembly_long/reads/${1}.nohuman.fastq --out-dir $T --threads $CPU
mkdir -p $1 && mv -f /tmp/flye_${1##*/} $1

conda deactivate
