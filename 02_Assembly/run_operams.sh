CPU=40
. ~/enable_conda.rc
conda activate operams
T=/tmp/op_${4##*/}
operams --short-read-assembler spades --long-read $1 --short-read1 $2 --short-read2 $3 --out-dir $T --num-processors $CPU &> ${4##*/}.log
mkdir -p ${4%/*} && mv -f /tmp/op_${4##*/} $4
conda deactivate
