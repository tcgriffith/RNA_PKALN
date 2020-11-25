

casedir=$1


caseid=$(basename $casedir)


cd $casedir

## prep input files

mkdir -p input

esl-reformat a2m $caseid.sto > input/$caseid.a2m

esl-reformat afa $caseid.sto > input/$caseid.afa

esl-reformat fasta $caseid.sto > input/$caseid.unaligned.fasta


## run cmalign

mkdir -p cmalign

cmalign -o cmalign/$caseid.cmalign.a2m --outformat A2M $caseid.cm  input/$caseid.unaligned.fasta 
# cmalign $caseid.cm $caseid.cm


## run gremlin

# gremlin_cpp -alphabet rna -i input/$caseid.afa -o input/$caseid.afa.gremlincpp -gap_cutoff 0.5 -mrf_o input/$caseid.afa.mrf
