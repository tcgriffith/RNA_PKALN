

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

DIR_PROJROOT=$(dirname $DIR)

casedir=$1


caseid=$(basename $casedir)


cd $casedir


## cmalign ss uses the ref

cat $caseid.sto |grep "^#" > tmp.sto
echo "//" >> tmp.sto
esl-reformat --keeprf --mingap stockholm tmp.sto |grep "#=GC SS" > ss.cmalign

## mrf ss may be slightly different
$DIR_PROJROOT/R/ss2mrfidx.R $caseid.sto input/$caseid.afa.mrf ss.rnamrf

for subdir in cmalign rnamrf;do

    if [ ! -f $subdir/$caseid.$subdir.rmgap.sto ]; then

	    ## keep ref
	    ~/bin/esl-reformat --replace acgturyswkmbdhvn:................ a2m $subdir/$caseid.$subdir.a2m > $subdir/$caseid.$subdir.rmgap
	    ## add ss
	    ~/bin/esl-reformat stockholm $subdir/$caseid.$subdir.rmgap > $subdir/$caseid.$subdir.rmgap.tmp
	    head -n -1 $subdir/$caseid.$subdir.rmgap.tmp > $subdir/$caseid.$subdir.rmgap.tmp2
	    cat $subdir/$caseid.$subdir.rmgap.tmp2 ss.$subdir > $subdir/$caseid.$subdir.rmgap.sto
	    echo "//" >> $subdir/$caseid.$subdir.rmgap.sto
    fi     

    scorescif=$(RNAalifold -q --sci $subdir/$caseid.$subdir.rmgap.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')


	echo "$subdir $caseid  $scorescif"
done

score3scif=$(RNAalifold -q --sci $caseid.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')


echo "input $caseid  $score3scif"

# ref=input/$caseid.afa

# test1=cmalign/$caseid.cmalign.a2m

# test2=rnamrf/$caseid.rnamrf.a2m

# ~/bin/esl-reformat stockholm $test1 >${test1%.a2m}.sto

# ~/bin/esl-reformat stockholm $test2 >${test2%.a2m}.sto

# # esl-reformat afa $test1 >${test1%.a2m}.afa

# # esl-reformat afa $test2 >${test2%.a2m}.afa

# # score1=$(compalignp -r $ref -t ${test1%.a2m}.afa)
# # RNAalifold 2.4.13
# score1scif=$(RNAalifold -q --sci ${test1%.a2m}.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')

# # score2=$(compalignp -r $ref -t ${test2%.a2m}.afa)

# score2scif=$(RNAalifold -q --sci ${test2%.a2m}.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')

# score3scif=$(RNAalifold -q --sci $caseid.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')


# # echo "method caseid SPS SCIF"
# echo "cmalign $caseid  $score1scif"

# echo "rnamrf $caseid  $score2scif"

# echo "input $caseid  $score3scif"
