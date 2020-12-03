


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)


cd $casedir

## prep input files


## For simulation
if [ -f seq1000_gaps.a2m ]; then

    seqss=$(cat seq1000_gaps.ss)

    esl-reformat stockholm seq1000_gaps.a2m > seq1000_gaps.sto.tmp

    head -n -1 seq1000_gaps.sto.tmp > seq1000_gaps.sto
    echo "#=GC SS_cons $seqss" >> seq1000_gaps.sto
    echo "//" >> seq1000_gaps.sto

    cp seq1000_gaps.sto $caseid.sto
    cmbuild --hand $caseid.cm $caseid.sto

fi


if [ ! -f input/$caseid.unaligned.fasta ]; then

    mkdir -p input

    esl-reformat  a2m $caseid.sto > input/$caseid.a2m

    esl-reformat  afa $caseid.sto > input/$caseid.afa
    
    # esl-reformat --replace acgturyswkmbdhvn:................ a2m $caseid.sto > input/$caseid.rmgap.a2m


    esl-reformat fasta $caseid.sto > input/$caseid.unaligned.fasta


else
    echo "# input ok"

fi

## run cmalign

if [ ! -f cmalign/$caseid.cmalign.a2m ]; then

mkdir -p cmalign

cmalign -o cmalign/$caseid.cmalign.a2m --outformat A2M $caseid.cm  input/$caseid.unaligned.fasta 



esl-reformat  afa cmalign/$caseid.cmalign.a2m > cmalign/$caseid.cmalign.afa


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


    esl-reformat  afa rnamrf/$caseid.rnamrf.a2m > rnamrf/$caseid.rnamrf.afa


else
    echo "# mrf ok"
fi
