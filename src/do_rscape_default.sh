

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)



cd $casedir

## extract SS
# cat $caseid.sto |grep "SS_cons" |sed 's/\.//g' > ss.dbn

cat $caseid.sto |grep "^#" > tmp.sto
echo "//" >> tmp.sto

## cmalign ss uses the ref
esl-reformat --keeprf --mingap stockholm tmp.sto |grep "#=GC SS" > ss.cmalign

## mrf ss may be slightly different
$DIR_PROJROOT/R/ss2mrfidx.R $caseid.sto input/$caseid.afa.mrf ss.rnamrf

# esl-reformat --replace acgturyswkmbdhvn:................ afa cmalign/$caseid.cmalign.a2m > cmalign/$caseid.cmalign.afa

for subdir in cmalign rnamrf;do
    ## keep ref
    ~/bin/esl-reformat --replace acgturyswkmbdhvn:................ a2m $subdir/$caseid.$subdir.a2m > $subdir/$caseid.$subdir.rmgap
    ## add ss
    ~/bin/esl-reformat stockholm $subdir/$caseid.$subdir.rmgap > $subdir/$caseid.$subdir.rmgap.tmp
    head -n -1 $subdir/$caseid.$subdir.rmgap.tmp > $subdir/$caseid.$subdir.rmgap.tmp2
    cat $subdir/$caseid.$subdir.rmgap.tmp2 ss.$subdir > $subdir/$caseid.$subdir.rmgap.sto
    echo "//" >> $subdir/$caseid.$subdir.rmgap.sto

    ## benchmark by R-scape
    R-scape --nofigures --outdir $subdir --outname default.$subdir $subdir/$caseid.$subdir.rmgap.sto

done

R-scape  --nofigures --outdir input --outname default.input $caseid.sto 

# if [ ! -f cmalign/$caseid.cmalign. ]; then
#     esl-reformat afa cmalign/$caseid.cmalign.a2m > cmalign/$caseid.cmalign.afa
# fi
