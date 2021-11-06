// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Individual_Score_Test_cond(arma::mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::mat X_adj, arma::vec residuals)
{
	int i;

	// number of markers
	int p = G.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;
	// log(p-value)
	arma::vec pvalue_log;
	pvalue_log.zeros(p);

	// SE of Uscore
	arma::vec Uscore_se;
	Uscore_se.zeros(p);
	// Effect size estimation
	arma::vec Est;
	Est.zeros(p);
	// SE of Effect size estimation
	arma::vec Est_se;
	Est_se.zeros(p);

	double test_stat = 0;

	int q = Sigma_iX.n_cols;

	// arma::mat tSigma_iX_G;
	// tSigma_iX_G.zeros(q,p);

	arma::mat Cov;
	Cov.zeros(p,p);

	int PX_adj_col = X_adj.n_cols;
	int PX_adj_row = Sigma_i.n_rows;

	arma::mat PX_adj;
	PX_adj.zeros(PX_adj_row,PX_adj_col);
	PX_adj = Sigma_i*X_adj - Sigma_iX*cov*(trans(Sigma_iX)*X_adj);
	Cov = trans(Sigma_i*G)*G - trans(trans(Sigma_iX)*G)*cov*trans(Sigma_iX)*G - trans(trans(X_adj)*G)*inv(trans(X_adj)*X_adj)*(trans(PX_adj)*G) - trans(trans(PX_adj)*G)*inv(trans(X_adj)*X_adj)*(trans(X_adj)*G) + trans(trans(X_adj)*G)*inv(trans(X_adj)*X_adj)*trans(PX_adj)*X_adj*inv(trans(X_adj)*X_adj)*(trans(X_adj)*G);

	for(i = 0; i < p; i++)
	{
		Uscore_se(i) = sqrt(Cov(i , i));

		if (Cov(i , i) == 0)
		{
			pvalue_log(i) = 0;
			Uscore_se(i) = 0;
			Est(i) = 0;
			Est_se(i) = 0;
		}
		else
		{
			test_stat = pow(Uscore(i),2)/Cov(i,i);
			pvalue_log(i) = -R::pchisq(test_stat,1,false,true);
			
			Uscore_se(i) = sqrt(Cov(i,i));
			Est(i) = Uscore(i)/Cov(i,i);
			Est_se(i) = 1/Uscore_se(i);
		}

	}

	return List::create(Named("Score") = trans(Uscore), Named("Score_se") = Uscore_se, Named("pvalue_log") = pvalue_log, Named("Est") = Est, Named("Est_se") = Est_se);
}

