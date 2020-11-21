#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double str_sim(CharacterVector v1,CharacterVector v2){
  double len =v1.length();
  double sim=sum(v1==v2)/len;
  return(sim);
}

// [[Rcpp::export]]
NumericVector calc_wt(CharacterMatrix seqmat){
  int nseq=seqmat.nrow();
  NumericVector wt(nseq);
  
  for (int i=0;i < seqmat.nrow();i++){
    CharacterVector vi = seqmat(i,_);
    for(int j=i+1;j<seqmat.nrow();j++){
      CharacterVector vj = seqmat(j,_);
      double sim = str_sim(vi,vj);
      if (sim > 0.8){
        wt[i]++;
        wt[j]++;
      }
    }
  }
  
  return(wt+1);
}