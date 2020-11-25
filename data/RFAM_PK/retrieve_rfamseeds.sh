while read rfamid; do
    #statements
   echo $rfamid
   mkdir -p RFAM/$rfamid

   wget http://rfam.org/family/$rfamid/cm -O RFAM/$rfamid/$rfamid.cm

   # wget http://rfam.org/family/$rfamid/alignment/stockholm -O RFAM/$rfamid/$rfamid.sto

done < $1