

casedir=$1


caseid=$(basename $casedir)


cd $casedir

ref=input/$caseid.afa

test1=cmalign/$caseid.cmalign.afa
test2=rnamrf/$caseid.rnamrf.afa

score1=$(compalignp -r $ref -t $test1)

score2=$(compalignp -r $ref -t $test2)

echo "cmalign $caseid $score1"

echo "rnamrf $caseid $score2"
