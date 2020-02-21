

Align an RNA to a family

- blastn: not good, sequence only

- profile HMM: CM without SS. no pseudoknots due to HMM
 
- covariance model: ordered-tree based with pair info. generalized HMM but still no pseudoknots.

- mapalign: align two matrices. with contacts


###

$ g++ -O3 -std=c++0x -o MRFaln MRFaln.cpp