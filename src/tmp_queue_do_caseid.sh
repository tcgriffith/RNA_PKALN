exe=/export/home/s5029518/GIT/RNA_PKALN/src/do_caseid.sh

dirlist=$1

if [ -f $dirlist ];then
    while read inputdir;do
        cd $inputdir
        queue_job.py -c1 --mem 6 -q aspen -- $exe $inputdir
    done <$dirlist
fi
