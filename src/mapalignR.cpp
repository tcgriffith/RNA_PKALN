#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

NumericVector test(NumericVector x){
	int n=x.size();
	NumericVector out(n);

	out[0]= x[0];
	for (int i = 0; i < n; ++i)
	{
		out[i]=out[i-1]+x[i];
	}

	return out;
}


// [[Rcpp::export]]

double Falign(NumericMatrix sco_mtx){

	int rows=sco_mtx.nrow();
	int cols=sco_mtx.ncol();



    double max_sco = 0;
    double sco[rows+1][cols+1]; memset(sco, 0, sizeof(sco));
    for (int i = 1; i <= rows; i++){
        for (int j = 1; j <= cols; j++){
            double A = sco[i-1][j-1] + sco_mtx(i-1,j-1);
            double D = sco[i-1][j];
            double R = sco[i][j-1];
            
            if(A >= R){if(A >= D){sco[i][j] = A;}else{sco[i][j] = D;}}
            else{if(R >= D){sco[i][j] = R;}else{sco[i][j] = D;}}
            
            if(sco[i][j] > max_sco){max_sco = sco[i][j];}
        }
    }
    return(max_sco);
}