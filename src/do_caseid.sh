


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1

## a case dir requires:
## $caseid.sto
## $caseid.ss
## $caseid.cm


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

    rm seq1000_gaps.sto.tmp
else
    echo "# cm ok"
fi

if [ ! -f $caseid.ss ]; then
   cat $caseid.sto |grep "^#.*SS_cons" | awk '{print $3}' > $caseid.ss
fi

if [ ! -f input/$caseid.cm ]; then 
    cmbuild --hand input/$caseid.cm $caseid.sto
fi

## input seqs
if [ ! -f input/$caseid.unaligned.fasta ]; then

    mkdir -p input



    esl-reformat  a2m $caseid.sto > input/$caseid.a2m # for benchmark 
    esl-reformat  afa $caseid.sto > input/$caseid.afa # for training GREMLIN; GREMLIN only knows fasta
    esl-reformat fasta $caseid.sto > input/$caseid.unaligned.fasta # for alignment, gaps removed


else
    echo "# input ok"

fi

## run cmalign

if [ ! -f cmalign/$caseid.cmalign.a2m ]; then

    mkdir -p cmalign
    cmalign -o cmalign/$caseid.cmalign.sto input/$caseid.cm  input/$caseid.unaligned.fasta 
    esl-reformat a2m cmalign/$caseid.cmalign.sto > cmalign/$caseid.cmalign.a2m
else
    echo "# cmalign ok"
fi



## mrf from GREMLIN
mkdir -p input/mrfs/

if [ ! -f input/mrfs/rnamrf.mrf ]; then
  gremlin_cpp -alphabet rna -i input/$caseid.afa -o input/mrfs/$caseid.afa.gremlincpp  -mrf_o input/mrfs/rnamrf.mrf 

else
    echo "# gremlin mrf ok"

fi

## mrf from SS
if [ ! -f input/mrfs/rnamrf_ssmrf.mrf ]; then
   $DIR_PROJROOT/R/ss2mrf.R $caseid.ss > input/mrfs/rnamrf_ssmrf.mrf 

else
    echo "# SSmrf ok"

fi


if [ -d input/mrfs/ ]; then

    for filemrf in input/mrfs/*.mrf; do

        amrf=$(basename $filemrf) # ssmrf.mrf
        wd=${amrf%.mrf}    # ssmrf, working directory

        if [ ! -f $wd/$caseid.$wd.a2m ]; then

            mkdir -p $wd
            echo "# Running $wd ok"

            cp $filemrf $wd/$caseid.$wd.mrf

            time $DIR_PROJROOT/R/RUN_mrfaln.R $wd/$caseid.$wd.mrf input/$caseid.unaligned.fasta  $wd/$caseid.$wd.a2m
            esl-reformat  afa $wd/$caseid.$wd.a2m > $wd/$caseid.$wd.afa
        else
            echo "# $wd ok"
        fi

    done

fi