

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)


cd $casedir



## run mrfalign

gremlin_cpp -alphabet rna -i input/$caseid.afa -o input/$caseid.afa.gremlincpp -gap_cutoff 0.5 -mrf_o input/$caseid.afa.mrf


mkdir -p rnamrf

$DIR_PROJROOT/R/RUN_mrfaln.R input/$caseid.afa.mrf input/$caseid.unaligned.fasta rnamrf/$caseid.rnamrf.a2m