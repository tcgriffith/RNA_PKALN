

casedir=$1


caseid=$(basename $casedir)


cd $casedir



# esl-reformat --replace acgturyswkmbdhvn:................ afa cmalign/$caseid.cmalign.a2m > cmalign/$caseid.cmalign.afa


# esl-reformat  afa cmalign/$caseid.cmalign.a2m > cmalign/$caseid.cmalign.afa


# gremlin_cpp -alphabet rna -i cmalign/$caseid.cmalign.afa -o cmalign/$caseid.cmalign.afa.gremlincpp -gap_cutoff 0.5 -mrf_o cmalign/$caseid.cmalign.afa.mrf


# # esl-reformat --replace acgturyswkmbdhvn:................ afa rnamrf/$caseid.rnamrf.a2m > rnamrf/$caseid.rnamrf.afa

# esl-reformat  afa rnamrf/$caseid.rnamrf.a2m > rnamrf/$caseid.rnamrf.afa


# gremlin_cpp -alphabet rna -i rnamrf/$caseid.rnamrf.afa -o rnamrf/$caseid.rnamrf.afa.gremlincpp -gap_cutoff 0.5 -mrf_o rnamrf/$caseid.rnamrf.afa.mrf
