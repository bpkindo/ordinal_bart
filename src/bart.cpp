#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"



using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List ordinal_bart_main(arma::ivec const& y, 
					arma::mat const& x, 
					arma::mat const& xtest, 
					int burn, 
					int nd, 
					int m,
					int nc, 
					double pbd, 
					double pb, 
					double alpha, 
					double betap, 
					double kappa,
					arma::mat const& Ad, 
					double s, 
					arma::mat const& inc_root, 
					arma::vec const& dstarbar) {


  arma::ivec distincty = arma::unique(y);
  size_t K = distincty.size();
  

   //--------------------------------------------------
   //random number generation
  RNGScope scope;         // Initialize Random number generator
   

   //--------------------------------------------------
   //read in data
   //read y NOTE: we assume y is already centered at 0.
  
   sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.

	allys.sy = arma::as_scalar(arma::sum(y)); // sum of y
    allys.sy2 = arma::as_scalar(arma::trans(y)*y); // sum of y^2
   
   size_t n = y.size();
   if(n<1) {
      Rcpp::Rcout << "error n<1\n";
      return 1;
   }
   allys.n = n;
   Rcpp::Rcout << "\ny read in:\n";
   Rcpp::Rcout << "n: " << n << std::endl;

   //read x
   //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
   size_t p = x.size()/n;
   if(x.size() != n*p) {
      Rcpp::Rcout << "error: input x file has wrong number of values\n"; 
      return 1;
   }
   Rcpp::Rcout << "\nx read in:\n";
   Rcpp::Rcout << "p: " << p << std::endl;
   
   Rcpp::Rcout << "first row: " <<  x(0,0) << " ...  " << x(0,p-1) << std::endl;
   Rcpp::Rcout << "last row: " << x(n-1,0) << " ...  " << x(n-1,p-1) << std::endl;
   
   std::vector<double> xvec;
   
   for(size_t i=0; i<n;i++){
	 for(size_t j=0; j<p; j++){
		xvec.push_back(x(i,j));
	 }
   }
   

   //--------------------------------------------------
   //dinfo
   arma::vec allfit(n);
   for(size_t i=0;i<n;i++) allfit(i)=0;
   double* r = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
   double* ftemp = new double[n]; //fit of current tree
   dinfo di;
   di.n=n; di.p=p; di.x = &xvec[0]; di.y=r; //the y for each draw will be the residual 
   
   
   //x for predictions

   dinfo dip; //data information for prediction
   dip.n=0;
   std::vector<double> xpvec;     //prediction observations, stored like x
   
    size_t np = xtest.n_rows;
   for(size_t i=0; i < np;i++){
	 for(size_t j=0; j < p; j++){
		xpvec.push_back(xtest(i,j));
	 }
   }

      if(xtest.size() != np*p) {
         Rcpp::Rcout << "error, wrong number of elements in prediction data set\n";
         return 1;
      }
	 dip.n=np; dip.p=p; dip.x = &xpvec[0]; dip.y=0; //there are no y's!
   
   Rcpp::Rcout << "\nx for prediction read in:\n";
   Rcpp::Rcout << "n: " << dip.n << std::endl;
   
   if(dip.n) {
      Rcpp::Rcout << "first row: " <<  dip.x[0] << " ...  " << dip.x[p-1] << std::endl;
      Rcpp::Rcout << "last row: " << dip.x[(dip.n-1)*p] << " ...  " << dip.x[dip.n*p-1] << std::endl;
   }

   //--------------------------------------------------
   double lambda = 1.0; //this one really needs to be set
   double nu = 3.0;
   double kfac=kappa;
   
  
   Rcpp::Rcout <<"\nburn,nd,number of trees: " << burn << ", " << nd << ", " << m << std::endl;
   Rcpp::Rcout <<"\nlambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << std::endl;

   //--------------------------------------------------
   //x cutpoints
   xinfo xi;
   makexinfo(p,n,&xvec[0],xi,(size_t)nc);

   //--------------------------------------------------
   //trees
   std::vector<tree> t(m);
   for(size_t i=0;i<(size_t)m;i++) t[i].setm(0.0); //if you sum the fit over the trees you get the fit.
   //--------------------------------------------------
  //prior and mcmc
   pinfo pi;
   pi.pbd=pbd; //prob of birth/death move
   pi.pb=pb; //prob of birth given  birth/death

   pi.alpha=alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi.beta=betap; //2 for bart means it is harder to build big trees.
   pi.tau=(3.0)/(kfac*sqrt((double)m));
   pi.sigma=1;


   //ordinal probit 
   arma::vec z(n);
   arma::vec ztest(dip.n);
   Rcpp::List metropout;
   int ncuts = K+1;
   int ncut = ncuts-3;
   int ndstar = K-2;
   arma::vec cutoff1(n);
   arma::vec cutoff2(n);
   arma::mat cutdraw(nd, ncuts);

   arma::vec olddstar(ndstar); 
   olddstar.zeros();
   arma::vec cutoffs = dstartoc(olddstar); 

   arma::mat ucholinv = arma::solve(arma::trimatu(arma::chol(Ad)), arma::eye(ndstar,ndstar));
  //double oldll = lldstar(olddstar, y, X*betahat);
   arma::mat Adinv = ucholinv*arma::trans(ucholinv);
   arma::mat rootdi = chol(Adinv);
   arma::vec sigma(n); sigma.ones();

   double oldll = lldstar(olddstar, y, allfit);
   //--------------------------------------------------
   arma::imat vec_class_prob_train(n,nd);
   arma::imat vec_class_prob_test(n,nd);
   size_t pclass_tmp = 0;

   //out of sample fit
   arma::vec allfit_test; //posterior mean for prediction
   double* fpredtemp=0; //temporary fit vector to compute prediction
   if(dip.n) {
      allfit_test.resize(dip.n);
      fpredtemp = new double[dip.n];
      allfit_test.fill(0.0);
   }
   //for sigma draw
   //--------------------------------------------------
   //mcmc

   Rcpp::Rcout << "\nMCMC:\n";
   clock_t tp;
   tp = clock();
   for(size_t i=0;i<(nd+burn);i++) {
   
      if(i%100==0) Rcpp::Rcout << "i: " << i << std::endl;
	  
	      //draw z given allfit, sigma, y, cut-offs
    for (size_t k=0; k<di.n; k++){
      cutoff1(k) = cutoffs(y(k)-1);
      cutoff2(k) = cutoffs(y(k));
    }
    z = rtrunVec(allfit, sigma, cutoff1, cutoff2);

	
   //draw trees
      for(size_t j=0;j<(size_t)m;j++) {
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) {
            allfit(k) = allfit(k)-ftemp[k];
            r[k] = z(k)-allfit(k);
         }
         bd(t[j],xi,di,pi);
         drmu(t[j],xi,di,pi);
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) allfit(k) += ftemp[k];
      }
	  
	      //draw gamma given z - these are the cutpoints on the latent rv
    metropout = dstarRwMetrop(y,allfit,olddstar,s,inc_root,dstarbar,oldll,rootdi, ncut);   
    olddstar = Rcpp::as<arma::vec>(metropout["dstardraw"]); //conversion from Rcpp to Armadillo requires explict declaration of variable 
	oldll =  Rcpp::as<double>(metropout["oldll"]);
    cutoffs = dstartoc(olddstar);
	

  if(i>=burn) {
      for(size_t k=0; k<di.n;k++) z(k) = allfit(k) + pi.sigma*Rcpp::rnorm(1)[0];
  
		 for(size_t k=0;k<di.n;k++) {
			for(size_t l=0; l < K; l++){
				if(z(k) <= cutoffs(l+1) && z(k) > cutoffs(l) ) pclass_tmp = l+1;
			}
		vec_class_prob_train(k,i-burn) = pclass_tmp;
	 }
	 
         if(dip.n) {
			allfit_test.fill(0.0);
            for(size_t j=0;j<(size_t)m;j++) {
               fit(t[j],xi,dip,fpredtemp);
               for(size_t k=0;k<dip.n;k++) allfit_test(k) += fpredtemp[k];
            }
			
      for(size_t k=0; k<dip.n;k++) ztest(k) = allfit_test(k) + pi.sigma*Rcpp::rnorm(1)[0];
				
				
			for(size_t k=0;k<dip.n;k++){
				for(size_t l=0; l < K; l++){
					if(ztest(k) <= cutoffs(l+1) && ztest(k) > cutoffs(l)) pclass_tmp = l+1;
				}
						vec_class_prob_test(k,i-burn) = pclass_tmp;
			}
			
         }
		 
		cutdraw.row(i-burn) = arma::trans(cutoffs);
      }
	  
	  
   }
      

   tp=clock()-tp;
   double thetime = (double)(tp)/(double)(CLOCKS_PER_SEC);
   Rcpp::Rcout << "time for loop: " << thetime << std::endl;

   arma::mat class_prob_train(di.n,K);
   arma::mat class_prob_test(dip.n,K);
   
for(size_t k=0; k<di.n;k++){
  for(size_t comp = 0; comp < K; comp++){
    class_prob_train(k, comp)=arma::sum(vec_class_prob_train.row(k) == (comp+1));
  }
 }

 
 for(size_t k=0; k<dip.n;k++){
  for(size_t comp = 0; comp < K; comp++){
    class_prob_test(k, comp)=arma::sum(vec_class_prob_test.row(k) == (comp+1));
  }
 }
 
  class_prob_train /=(double)nd;
  class_prob_test /=(double)nd;
return Rcpp::List::create(
            Rcpp::Named("class_prob_train") = class_prob_train,
			Rcpp::Named("class_prob_test") = class_prob_test,
			Rcpp::Named("cutdraw") = cutdraw
			);
}
