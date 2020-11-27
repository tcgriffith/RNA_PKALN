


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)


cd $casedir

## prep input files

if [ ! -f input/$caseid.unaligned.fasta ]; then

    mkdir -p input

    esl-reformat  a2m $caseid.sto > input/$caseid.a2m

    esl-reformat  afa $caseid.sto > input/$caseid.afa

    esl-reformat fasta $caseid.sto > input/$caseid.unaligned.fasta


else
    echo "# input ok"

fi

## run cmalign

if [ ! -f cmalign/$caseid.cmalign.a2m ]; then

mkdir -p cmalign

cmalign -o cmalign/$caseid.cmalign.a2m --outformat A2M $caseid.cm  input/$caseid.unaligned.fasta 

else
    echo "# cmalign ok"

fi


if [ ! -f input/$caseid.afa.mrf ]; then
  gremlin_cpp -alphabet rna -i input/$caseid.afa -o input/$caseid.afa.gremlincpp  -mrf_o input/$caseid.afa.mrf

else
    echo "# gremlin ok"

fi

## run RNAmrfalign

if [ ! -f rnamrf/$caseid.rnamrf.a2m ]; then

    mkdir -p rnamrf

   time $DIR_PROJROOT/R/RUN_mrfaln.R input/$caseid.afa.mrf input/$caseid.unaligned.fasta  rnamrf/$caseid.rnamrf.a2m
else
    echo "# mrf ok"
fi