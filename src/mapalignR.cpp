#include <Rcpp.h>
// #include <math.h>
// #include <algorithm>
using namespace Rcpp;


double exp_fast(double x){
  // WARNING fails if |x| > 1024
  //https://codingforspeed.com/using-faster-exponential-approximation/
  x = 1 + x/1024;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}

double gaussian(double mean, double stdev, double x){return exp_fast(-pow((x - mean),2)/(2*(pow(stdev,2))));}


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

int id2_to_id1(int i, int j, int dim_j){
    int id1=(j-1) *dim_j + i;
    return(id1);
}


IntegerVector id1_to_id2(int id1, int dim_i){
    IntegerVector id2(2);
    
    int mod=id1 % dim_i;
    
    if (mod == 0){
      id2(0) = dim_i;
      id2(1) = id1 / dim_i;
    } else{
      id2(0) = mod;
      id2(1) = id1 / dim_i +1;
    }

    return(id2);

}



// [[Rcpp::export]]

double retrieve_matj(int i, int a, int j, int b, NumericMatrix mat_j, int len_a){
    return(mat_j( (i)*len_a + a, (j)*len_a + b));
}

double retrieve_mat_h(int i, int a, NumericMatrix mat_h){
  return(mat_h(i,a));
}


double sepw(double sep){if(sep <= 4){return 0.50;}else if(sep == 5){return 0.75;}else{return 1.00;}}


// [[Rcpp::export]]

NumericMatrix ini_SCO(IntegerVector seq, NumericMatrix mrf_mat, int mrf_len, double sep_x, double sep_y) {


  // SCO = matrix(0, nrow = length(seq), ncol = mrf_len)
  NumericMatrix SCO(seq.length(), mrf_len);
  

  
  
  for (int ai = 0; ai < seq.length(); ++ai){
    for (int bi = 0; bi < mrf_len; ++bi){
      NumericMatrix M(seq.length(), mrf_len);
        for (int aj = 0; aj < seq.length(); aj++){
            for (int bj=0;bj < mrf_len;bj++){

                int nt_ai=seq(ai);
                int nt_aj=seq(aj);
                double score_a2b;

                if (bi > bj){
                    score_a2b=retrieve_matj(bj,nt_aj, bi,nt_ai,mrf_mat, 5);
                }
                else{
                    score_a2b=retrieve_matj(bi,nt_ai,bj,nt_aj,mrf_mat,5);
                }

                int sep_a = abs(ai - aj);
                int sep_b = abs(bi - bj);
                int sep_D = abs(sep_a - sep_b);
                int sep_M = std::min(sep_a, sep_b);
                double sep_std = sep_y * (1 + pow(sep_M - 2,sep_x));

                if (sep_D/sep_std <6){
                    M(aj,bj) =  score_a2b *  sepw(sep_M) * gaussian(0,sep_std,sep_D);
                }else{
                    M(aj,bj) = 0;
                }
            }
        }

        SCO(ai,bi)=Falign(M);
    }
  }
  
  // NumericMatrix SCO_h(seq.length(),mrf_len);
  // 
  // for (int ai = 0; ai < seq.length(); ++ai){
  //   for (int bi = 0; bi < mrf_len; ++bi){
  //     int nt_ai=seq(ai);
  //     SCO_h(ai,bi) = mrf_h(bi, nt_ai) +SCO(ai,bi);
  //   }
  // }
  

  return(SCO);
  
}




// 
// vec_int mod_SCO(double do_it, vec_double &gap_a, vec_double &gap_b, double gap_e, mtx_double &SCO, mtx_double &P_SCO,
//                 vec_int &vec_a_div,vec_int &vec_b_div,vec_int &vec_a, vec_int &vec_b,mtx_int &vec_a_i, mtx_int &vec_b_i,mtx_double &mtx_a, mtx_double &mtx_b){
//     // iterate
//     vec_int a2b_tmp;
//     for(int it=0; it < do_it; it++)
//     {
//         // align
//         a2b_tmp = align(gap_a,gap_b,gap_e,SCO,P_SCO);
//         
//         // update similarity matrix
//         double IT = (double)it + 1;
//         double s1 = (IT/(IT+1)); double s2 = (1/(IT+1));
//         for(int a=0; a < vec_a.size(); a++){ // go through columns (vec_a) in map_a that has contacts
//             int ai = vec_a[a];
//             for(int b=0; b < vec_b.size(); b++){ // go through columns (vec_b) in map_b that has contacts
//                 int bi = vec_b[b];
//                 double sco_contact = 0;
//                 for(int n=0; n < vec_a_i[ai].size(); n++){ // go through contacts in vec_a
//                     int aj = vec_a_i[ai][n];
//                     int bj = a2b_tmp[aj]; // get mapping
//                     if(bj != -1){ // if mapping exists
//                         if((ai > aj and bi > bj) or (ai < aj and bi < bj)){ // if ai-aj in same direction as bi-bj
//                             double sep_M = min(abs(ai-aj),abs(bi-bj));
//                             sco_contact += mtx_a[ai][aj] * mtx_b[bi][bj] * sepw(sep_M);
//                         }
//                     }
//                 }
//                 SCO[ai][bi] = s1*SCO[ai][bi] + s2*sco_contact;
//             }
//         }
//     }
//     return(a2b_tmp);
// }
// 
// 
// 
// 
// IntegerVector align(vec_double &gap_a, vec_double &gap_b, double &gap_e, mtx_double &sco_mtx){
//     // LOCAL_ALIGN
//     // Start    0
//     // [A]lign  1
//     // [D]own   2
//     // [R]ight  3
//     
//     double max_sco = 0;
//     int rows = sco_mtx.size();
//     int cols = sco_mtx[0].size();
//     
//     // vec_int a2b(rows,-1);
//     IntegerVector a2b(rows);
//     
//     mtx_double sco(rows+1,vector<double>(cols+1,0));
//     mtx_int label(rows+1,vector<int>(cols+1,0));
//     
//     int max_i = 0;int max_j = 0;
//     for (int i = 1; i <= rows; i++){
//         
//         for (int j = 1; j <= cols; j++){
//             double A = sco[i-1][j-1] + sco_mtx[i-1][j-1]; if(add_prf == true){A += p_sco_mtx[i-1][j-1];}
//             double D = sco[i-1][j];
//             double R = sco[i][j-1];
//             
//             if(label[i-1][j] == 1){D += gap_b[j-1];}else{D += gap_b[j-1] * gap_e;}
//             if(label[i][j-1] == 1){R += gap_a[i-1];}else{R += gap_a[i-1] * gap_e;}
//             
//             if(A <= 0 and D <= 0 and R <= 0){label[i][j] = 0;sco[i][j] = 0;}
//             else{
//                 if(A >= R){if(A >= D){label[i][j] = 1;sco[i][j] = A;}else{label[i][j] = 2;sco[i][j] = D;}}
//                 else{if(R >= D){label[i][j] = 3;sco[i][j] = R;}else{label[i][j] = 2;sco[i][j] = D;}}
//                 if(sco[i][j] > max_sco){max_i = i;max_j = j;max_sco = sco[i][j];}
//             }
//         }
//     }
//     int i = max_i;int j = max_j;
//     while(1){
//         if(label[i][j] == 0){break;}
//         else if(label[i][j] == 1){a2b[i-1] = j-1;i--;j--;}
//         else if(label[i][j] == 2){i--;}
//         else if(label[i][j] == 3){j--;}
//     }
//     return(a2b);
// }