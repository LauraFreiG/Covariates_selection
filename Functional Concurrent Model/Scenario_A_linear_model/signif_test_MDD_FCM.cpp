
#include <RcppArmadillo.h>
#include <math.h>
#include <vector>
using namespace arma;
using namespace Rcpp;
using namespace std;


//declarations
vec UCenter_X(int n, int p, vec X, vec XX);
vec UCenter_Y(int n, int q, vec Y, vec YY);



/*
// MDD_n calculation
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double MDD_n_fun(mat X, vec Y, int n, int p, int combi){
	vec MDD_vec(p), X_jj(n);
	mat A(n,n), B(n,n);
	int c_n;
	double MDD=0, suma_S=0.0, S=0.0, T_n=0.0, A_j[n][n][p];
	for(int jj=0; jj<p; jj++){
		X_jj=X.col(jj);
		double MDD=0.0;
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
			double A_l=0.0, B_l=0.0, A_k=0.0, B_k=0.0, A_kl=0.0, B_kl=0.0;
			for(int l=0; l<n; l++){ 
				A_l+= abs(X_jj[i]-X_jj[l]);
				B_l+= pow(Y[i]-Y[l],2)/2;
				for(int k=0; k<n; k++){	A_kl += abs(X_jj[k]-X_jj[l]);  B_kl += pow(Y[k]-Y[l],2)/2; }
			}
			for(int k=0; k<n; k++){ A_k+= abs(X_jj[k]-X_jj[j]);  B_k += pow(Y[k]-Y[j],2)/2; }
			A(i,j)=	(X_jj[i]-X_jj[j]) -(1/(n-2))*A_l -(1/(n-2))*A_k +(1/(n-1)*(n-2))*A_kl;
			A_j[i][j][jj]=A(i,j);
			B(i,j)= pow(Y[i]-Y[j],2)/2 -(1/(n-2))*B_l -(1/(n-2))*B_k +(1/(n-1)*(n-2))*B_kl;
			if(j!=i){MDD += A(i,j)*B(i,j);}
			}
		}
		MDD_vec[jj] = (1/(n*(n-3)))*MDD ;
	}
	
	MDD=sum(MDD_vec);
	
	c_n= pow(n-3,4)/pow(n-1,4) + (2*pow(n-3,4))/(pow(n-1,4)*pow(n-2,3)) + (2*(n-3))/(pow(n-1,4)*pow(n-2,3));
	
	for(int k=0; k<n; k++){
		for(int l=k+1; l<n; l++){
			for(int jj=0; jj<p; jj++){
				for(int jj2=0; jj2<p; jj2++){
					suma_S += A_j[k][l][jj]*A_j[k][l][jj2]*pow(B(k,l),2);
				}	
			}			
		}	
	}
	
	
	S= (2/(n*(n-1)*c_n))* suma_S;
	T_n=sqrt(combi)*(MDD/sqrt(S));
	
	return T_n;
} 
*/



/*
// MDD_n calculation
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void MDD_UCenter(int *N, int *P, int *Q, double *X, double *Y, double *V) {
  int n = N[0];
  int p = P[0];
  int q = Q[0];

  double XX[n * n];
  memset(XX, 0, n * n * sizeof(double));

  double YY[n * n];
  memset(YY, 0, n * n * sizeof(double));

  UCenter_X(n, p, X, XX);
  UCenter_Y(n, q, Y, YY);

  //V[0] = inner_UCenter(n, XX, YY);
  //List L = List::create( Named("XX")=XX, Named("YY")=YY);
  
  //return L;
  
}
*/



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec UCenter_X(int n, int p, vec X, vec XX) {
  double row_sum[n];
  memset(row_sum, 0, n * sizeof(double));

  double col_sum[n];
  memset(col_sum, 0, n * sizeof(double));

  double total_sum = 0;
  double a;
  double temp;

  int i;
  int j;
  int k;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      if (i != j) {
        a = 0;

        // XX[i, j] = |X[i, ] - X[j, ]|
        for (k = 0; k < p; ++k) {
          temp = X[i + k * n] - X[j + k * n];
          a += temp * temp;
        }

        a = sqrt(a);

        XX[i + j * n] = a;
        row_sum[i] += a;
        col_sum[j] += a;
        total_sum += a;
      }
    }
  }

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // diag(XX) = 0
      if (i != j) {
        XX[i + j * n] = XX[i + j * n] - row_sum[i] / (n - 2) - col_sum[j] / (n - 2) + total_sum / (n - 1) / (n - 2);
      }
    }
  } 
  
  return XX;

}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec UCenter_Y(int n, int q, vec Y, vec YY) {
  double row_sum[n];
  memset(row_sum, 0, n * sizeof(double));

  double col_sum[n];
  memset(col_sum, 0, n * sizeof(double));

  double total_sum = 0;
  double b;
  double temp;

  int i;
  int j;
  int k;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      if (i != j) {
        b = 0;

        // YY[i, j] = 0.5 * |Y[i, ] - Y[j, ]|^2 or YY[i, j] = -Y[i, ]^T * Y[j, ]
        for (k = 0; k < q; ++k) {
          temp = Y[i + k * n] - Y[j + k * n];
          b += temp * temp;
        }

        b = 0.5 * b;

        YY[i + j * n] = b;
        row_sum[i] += b;
        col_sum[j] += b;
        total_sum += b;
      }
    }
  }

  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      // diag(YY) = 0
      if (i != j) {
        YY[i + j * n] = YY[i + j * n] - row_sum[i] / (n - 2) - col_sum[j] / (n - 2) + total_sum / (n - 1) / (n - 2);
      }
    }
  }
  
    return YY;

}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double suma_S_boot(vec e, int n, int p, List A, mat B) {
	double suma=0.0;
	mat A1(n,n), A2(n,n);
	
		
	for (int k = 0; k < n; k++) {
			for (int l = k+1; l < n; l++) {
				for (int j1 = 0; j1 < p; j1++){
					A1= as<arma::mat>(A[j1]);
					for (int j2 = 0; j2 < p; j2++){
						A2= as<arma::mat>(A[j2]);
						suma += (A1(k,l)*A2(k,l)*pow(B(k,l),2)*pow(e[k],2)*pow(e[l],2) );						
					}
				}		
		}
		
	}
	
	
	return suma;	

}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double suma_S_boot_parcial(vec e, int n, int p, mat A, mat B) {
	double suma=0.0;
	
		
	for (int k = 0; k < n; k++) {
			for (int l = k+1; l < n; l++) {
				suma += (pow(A(k,l),2)*pow(B(k,l),2)*pow(e[k],2)*pow(e[l],2) );						
				}		
		}
		
	return suma;	

}








