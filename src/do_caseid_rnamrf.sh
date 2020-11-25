

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)


cd $casedir

## run cmalign

mkdir -p cmalign

cmalign -o cmalign/$caseid.cmalign.a2m --outformat A2M $caseid.cm  input/$caseid.test.fasta 

## run mrfalign

if [ ! -f input/$caseid.afa.mrf ]; then
    gremlin_cpp -alphabet rna -i input/$caseid.afa -o input/$caseid.afa.gremlincpp -gap_cutoff .8 -mrf_o input/$caseid.afa.mrf
fi


mkdir -p rnamrf

$DIR_PROJROOT/R/RUN_mrfaln.R input/$caseid.afa.mrf input/$caseid.test.fasta  rnamrf/$caseid.rnamrf.a2m
