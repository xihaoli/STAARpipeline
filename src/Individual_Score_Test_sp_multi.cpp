// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Individual_Score_Test_sp_multi(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, int n_pheno=1)
{
	int i,k;
	int p = G.n_cols;

	// number of markers
	int pp = p/n_pheno;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;
	// log(p-value)
	arma::vec pvalue_log;
	pvalue_log.zeros(pp);

	arma::mat Uscore_cov;
	Uscore_cov.zeros(n_pheno,n_pheno);

	arma::uvec id_single;
	id_single.zeros(n_pheno);

	double test_stat = 0;

	int q = Sigma_iX.n_cols;

	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(q,p);

	arma::mat Cov;
	Cov.zeros(p,p);

	tSigma_iX_G = trans(Sigma_iX)*G;
	Cov = trans(trans(Sigma_i*G)*G) - trans(tSigma_iX_G)*cov*tSigma_iX_G;

	arma::mat quad;
	quad.zeros(1,1);

	for(i = 0; i < pp; i++)
	{
		for(k = 0; k < n_pheno; k++)
		{
			id_single(k) = k*pp+i;
		}

		Uscore_cov = Cov(id_single,id_single);

		if (arma::det(Uscore_cov) == 0)
		{
			pvalue_log(i) = 0;
		}
		else
		{
			quad = trans(Uscore(id_single))*inv(Uscore_cov)*Uscore(id_single);
			test_stat = quad(0,0);
			pvalue_log(i) = -R::pchisq(test_stat,n_pheno,false,true);
		}

	}

	return List::create(Named("Score") = trans(Uscore), Named("pvalue_log") = pvalue_log);
}

