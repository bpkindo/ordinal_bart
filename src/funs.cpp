#include <cmath>
#include "funs.h"
#include <map>



arma::ivec numcomp(arma::ivec const& indic, int k){
  
  //find the number of times each of k integers is in the vector indic
  arma::ivec ncomp(k);
  
  for(int comp = 0; comp < k; comp++){
    ncomp[comp]=sum(indic == (comp+1));
  }
  
  return(ncomp);
}


arma::vec rtrunVec(arma::vec const& mu,arma::vec const& sigma, arma::vec const& a, arma::vec const& b){
  

//function to draw from univariate truncated norm
//a is vector of lower bounds for truncation
//b is vector of upper bounds for truncation

  int n = mu.size();
  arma::vec FA(n);
  arma::vec FB(n);
  arma::vec out(n);
  for (int i=0; i<n; i++) {
    FA[i] = R::pnorm((a[i]-mu[i])/sigma[i],0,1,1,0);
    FB[i] = R::pnorm((b[i]-mu[i])/sigma[i],0,1,1,0);
    out[i] = mu[i]+sigma[i]*R::qnorm(R::runif(0,1)*(FB[i]-FA[i])+FA[i],0,1,1,0);
  }

  return(out);
}

double lndMvn(arma::vec const& x, arma::vec const& mu, arma::mat const& rooti){


// function to evaluate log of MV Normal density with  mean mu, var Sigma
// Sigma=t(root)%*%root   (root is upper tri cholesky root)
// Sigma^-1=rooti%*%t(rooti)   
// rooti is in the inverse of upper triangular chol root of sigma
//          note: this is the UL decomp of sigmai not LU!
//                Sigma=root'root   root=inv(rooti)

  arma::vec z = arma::vectorise(arma::trans(rooti)*(x-mu));
  
  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}


Rcpp::List dstarRwMetrop(arma::ivec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root, arma::vec const& dstarbar, double oldll, arma::mat const& rootdi, int ncut){ 

// function to execute rw metropolis for the dstar
// y is n vector with element = 1,...,j 
// X is n x k matrix of x values 
// RW increments are N(0,s^2*t(inc.root)%*%inc.root)
// prior on dstar is N(dstarbar,Sigma)  Sigma^-1=rootdi*t(rootdi)
//  inc.root, rootdi are upper triangular
//  this means that we are using the UL decomp of Sigma^-1 for prior 
// olddstar is the current
//
  int stay = 0;
  double unif;
  arma::vec dstardraw;

  arma::vec dstarc = olddstar + s*arma::trans(inc_root)*arma::vec(Rcpp::rnorm(ncut));
  double cll = lldstar(dstarc, y, mu);
  double clpost = cll + lndMvn(dstarc, dstarbar, rootdi);
  double ldiff = clpost - oldll - lndMvn(olddstar, dstarbar, rootdi);
  double alpha = exp(ldiff);
  
  if (alpha>1){
    alpha = 1.0;
  } 

  if (alpha<1){
    unif = unif_rand(); 
  }
  else{
    unif = 0;
  }
  
  if (unif<=alpha){
    dstardraw = dstarc; 
    oldll = cll;
  }
  else{
    dstardraw = olddstar;
    stay = 1;
  }
  
  return Rcpp::List::create(
      Rcpp::Named("dstardraw") = dstardraw,
      Rcpp::Named("oldll") = oldll,
      Rcpp::Named("stay") = stay
  );
}   


// compute conditional likelihood of data given cut-offs
double lldstar(arma::vec const& dstar, arma::ivec const& y, arma::vec const& mu){
  arma::vec gamma = dstartoc(dstar);
  
  int ny = y.size();
  Rcpp::NumericVector gamma1(ny);
  Rcpp::NumericVector gamma2(ny);
  for (int i=0; i<ny; i++){
    gamma1[i] = gamma(y[i]);
    gamma2[i] = gamma(y[i]-1);
  }
  Rcpp::NumericVector temp = Rcpp::pnorm(gamma1-Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu)))-Rcpp::pnorm(gamma2-Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu))); //pnorm takes Rcpp type NumericVector, NOT arma objects of type vec
  arma::vec arg = Rcpp::as<arma::vec>(temp);
  double epsilon = 1.0/(10^-50);
  for (int j=0; j<ny; j++){
    if (arg[j]<epsilon){
      arg[j] = epsilon;
    }
  }
  return (sum(log(arg)));
}

//dstartoc is a fuction to transfer dstar to its cut-off value    
arma::vec dstartoc(arma::vec const& dstar){
  int ndstar = dstar.size();
  arma::vec c(ndstar+3);
  c[0] = -100;
  c[1] = 0;
  c(arma::span(2,ndstar+1)) = arma::cumsum(arma::exp(dstar));
  c[ndstar+2] = 100;
  
  return (c);
}
//--------------------------------------------------
// normal density N(x, mean, variance)
double pn(double x, double m, double v)
{
   double dif = x-m;
   return exp(-.5*dif*dif/v)/sqrt(2*PI*v);
}
//--------------------------------------------------
// draw from discrete distributin given by p, return index
int rdisc(double *p)
{

   double sum;
   double u = unif_rand();

    int i=0;
    sum=p[0];
    while(sum<u) {
       i += 1;
       sum += p[i];
    }
    return i;
}
//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os) 
{
   size_t p = xi.size();
   if(p!=2) {
      Rcpp::Rcout << "error in grm, p !=2\n";
      return;
   }
   size_t n1 = xi[0].size();
   size_t n2 = xi[1].size();
   tree::tree_cp bp; //pointer to bottom node
   double *x = new double[2];
   for(size_t i=0;i!=n1;i++) {
      for(size_t j=0;j!=n2;j++) {
         x[0] = xi[0][i]; 
         x[1] = xi[1][j]; 
         bp = tr.bn(x,xi);
         os << x[0] << " " << x[1] << " " << bp->getm() << " " << bp->nid() << std::endl;
      }
   }
   delete[] x;
}
//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots)
{
   double pb;  //prob of birth to be returned
   tree::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++) 
      if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else { 
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pi.pb;
   }
   return pb;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi)
{
   if(cansplit(n,xi)) {
      return pi.alpha/pow(1.0+n->depth(),pi.beta);
   } else {
      return 0.0;
   }
}
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x
   double y;          //current y

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   sv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) bnmap[bnv[i]]=i;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      y=di.y[i];

      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(sv[ni].n);
      sv[ni].sy += y;
      sv[ni].sy2 += y*y;
   }
}
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         y = di.y[i];
         if(xx[v] < xi[v][c]) {
               sl.n++;
               sl.sy += y;
               sl.sy2 += y*y;
          } else {
               sr.n++;
               sr.sy += y;
               sr.sy2 += y*y;
          }
      }
   }
}
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr)
{
   double *xx;//current x
   double y;  //current y
   sl.n=0;sl.sy=0.0;sl.sy2=0.0;
   sr.n=0;sr.sy=0.0;sr.sy2=0.0;

   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==nl) {
         y = di.y[i];
         sl.n++;
         sl.sy += y;
         sl.sy2 += y*y;
      }
      if(bn==nr) {
         y = di.y[i];
         sr.n++;
         sr.sy += y;
         sr.sy2 += y*y;
      }
   }
}
//--------------------------------------------------
//log of the integrated likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau)
{
   double yb,yb2,S,sig2,d;
   double sum, rv;

   yb = sy/n;
   yb2 = yb*yb;
   S = sy2 - (n*yb2);
   sig2 = sigma*sigma;
   d = n*tau*tau + sig2;
   sum = S/sig2 + (n*yb2)/d;
   rv = -(n*LTPI/2.0) - (n-1)*log(sigma) -log(d)/2.0;
   rv = rv -sum/2.0;
   return rv;
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv)
{
   double *xx;
   tree::tree_cp bn;
   fv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv)
{
   double *xx;
   tree::tree_cp bn;
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      fv[i] = bn->getm();
   }
}
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv)
{
   double *xx;
   tree::tree_cp bn;
   pv.resize(di.n);
   for(size_t i=0;i<di.n;i++) {
      xx = di.x + i*di.p;
      bn = t.bn(xx,xi);
      pv[i] = bn->nid();
   }
}
//--------------------------------------------------
// draw all the bottom node mu's

void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi)
{
   tree::npv bnv;
   std::vector<sinfo> sv;
   allsuff(t,xi,di,bnv,sv);

   double a = 1.0/(pi.tau * pi.tau);
   double sig2 = pi.sigma * pi.sigma;
   double b,ybar;

   for(tree::npv::size_type i=0;i!=bnv.size();i++) {
      b = sv[i].n/sig2;
      ybar = sv[i].sy/sv[i].n;
      bnv[i]->setm(b*ybar/(a+b) + Rcpp::rnorm(1)[0]/sqrt(a+b));
   }
}


//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
   Rcpp::Rcout << "xinfo: \n";
   for(size_t v=0;v!=xi.size();v++) {
      Rcpp::Rcout << "v: " << v << std::endl;
      for(size_t j=0;j!=xi[v].size();j++) Rcpp::Rcout << "j,xi[v][j]: " << j << ", " << xi[v][j] << std::endl;
   }
   Rcpp::Rcout << "\n\n";
}
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
   double xinc;

   //compute min and max for each x
   std::vector<double> minx(p,INFINITY);
   std::vector<double> maxx(p,-INFINITY);
   double xx;
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}
// get min/max needed to make cutpoints
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx)
{
   double xx;

   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
}
//make xinfo = cutpoints give the minx and maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx)
{
	double xinc;
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}


// construct E[f] over sample of draws from the posterior, overwrites exisiting values in postp
// if there are any.
void makepostpred(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp)
{
	double* fpredtemp=0; //temporary fit vector to compute prediction
	double* ppredmean=0; //temporary fit vector for mean from 1 tree
	size_t m,ndraws;

	fpredtemp = new double[dip.n];
	ppredmean = new double[dip.n];
	ndraws=t.size();
	m=t[0].size();
	for(size_t i=0;i<dip.n;i++) { 
		ppredmean[i]=0.0;
		postp[i]=0.0;
	}

	for(size_t i=0;i<ndraws;i++) {
		for(size_t j=0;j<m;j++) {
			fit(t[i][j],xi,dip,fpredtemp);
			for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
		}
		if(i>0)
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] *= (int)i;
				postp[k] += ppredmean[k];
				postp[k] /= (int)(i+1);
				ppredmean[k] = 0.0;
			}
		else
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] = ppredmean[k];
				ppredmean[k] = 0.0;
			}
	}
}

// construct E[f] and E[f^2] over sample of draws from the posterior, overwrites exisiting values in postp,
// postp2 if there are any.
void makepostpred2(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp, double *postp2)
{
	double* fpredtemp=0; //temporary fit vector to compute prediction
	double* ppredmean=0; //temporary fit vector for mean from 1 tree
	size_t m,ndraws;

	fpredtemp = new double[dip.n];
	ppredmean = new double[dip.n];
	ndraws=t.size();
	m=t[0].size();
	for(size_t i=0;i<dip.n;i++) { 
		ppredmean[i]=0.0;
		postp[i]=0.0;
		postp2[i]=0.0;
	}

	for(size_t i=0;i<ndraws;i++) {
		for(size_t j=0;j<m;j++) {
			fit(t[i][j],xi,dip,fpredtemp);
			for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
		}
		if(i>0)
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] *= (int)i;
				postp[k] += ppredmean[k];
				postp[k] /= (int)(i+1);
				postp2[k] *= (int)i;
				postp2[k] += ppredmean[k]*ppredmean[k];
				postp2[k] /= (int)(i+1);
				ppredmean[k] = 0.0;
			}
		else
			for(size_t k=0;k<dip.n;k++)
			{
				postp[k] = ppredmean[k];
				postp2[k] = ppredmean[k]*ppredmean[k];
				ppredmean[k] = 0.0;
			}
	}
}



