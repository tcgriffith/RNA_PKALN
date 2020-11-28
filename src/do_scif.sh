casedir=$1


caseid=$(basename $casedir)


cd $casedir

ref=input/$caseid.afa

test1=cmalign/$caseid.cmalign.a2m

test2=rnamrf/$caseid.rnamrf.a2m

~/bin/esl-reformat stockholm $test1 >${test1%.a2m}.sto

~/bin/esl-reformat stockholm $test2 >${test2%.a2m}.sto

# esl-reformat afa $test1 >${test1%.a2m}.afa

# esl-reformat afa $test2 >${test2%.a2m}.afa

# score1=$(compalignp -r $ref -t ${test1%.a2m}.afa)
# RNAalifold 2.4.13
score1scif=$(RNAalifold -q --sci ${test1%.a2m}.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')

# score2=$(compalignp -r $ref -t ${test2%.a2m}.afa)

score2scif=$(RNAalifold -q --sci ${test2%.a2m}.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')

score3scif=$(RNAalifold -q --sci $caseid.sto |grep "sci = "|sed -r 's/.*sci = (.*)]/\1/')


# echo "method caseid SPS SCIF"
echo "cmalign $caseid  $score1scif"

echo "rnamrf $caseid  $score2scif"

echo "input $caseid  $score3scif"
