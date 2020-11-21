# What's the issue

The best RNA homology search program: Infernal, doesn't include PK

# What's new

Applied MRF model to allow for RNA seq-to-family alignment

# Why is it important

- Pseudoknots is prevalent in ncRNA, MRFalign is suitable for aligning RNAs

- Also, as MRF models preserve pairwise information, homologous RNAs can also be well aligned




Align an RNA to a family

- blastn: not good, sequence only

- profile HMM: CM without SS. no pseudoknots due to HMM
 
- covariance model: ordered-tree based with pair info. generalized HMM but still no pseudoknots.

- mapalign: align two matrices. with contacts


###

$ g++ -O3 -std=c++0x -o MRFaln MRFaln.cpp