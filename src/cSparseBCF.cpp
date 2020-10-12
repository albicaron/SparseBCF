#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using namespace Rcpp;


//' C++ function called by R wrapper SparseBCF.R
//' @export



// y = m(x) + b(x)z + e, e~N(0, sigma^2_y

//x_con is the design matrix for m. It should have n = rows
//x_mod is the design matrix for b. It should have n = rows
//data should come in sorted with all trt first, then control cases

// [[Rcpp::export]]
List cSparseBCF(NumericVector y_, NumericVector z_,
                NumericVector x_con_, NumericVector x_mod_, NumericVector x_mod_est_,
                List x_con_info_list, List x_mod_info_list,
                arma::mat random_des, //needs to come in with n rows no matter what(?)
                arma::mat random_var, arma::mat random_var_ix, //random_var_ix*random_var = diag(Var(random effects))
                double random_var_df,
                int burn, int nd, int thin, //Draw nd*thin + burn samples, saving nd draws after burn-in
                int ntree_mod, int ntree_con,
                double lambda, double nu, //prior pars for sigma^2_y
                double con_sd, // Var(m(x)) = con_sd^2 marginally a priori (approx)
                double mod_sd, // Var(b(x)) = mod_sd^2 marginally a priori (approx)
                double con_alpha, double con_beta,
                double mod_alpha, double mod_beta,
				double a,
				double b,
				double rho,
				double rho_tau,
				double theta,
				double omega,
				NumericVector weights_mu,
				NumericVector weights_tau,
                CharacterVector treef_name_,
				bool dart,
				bool aug,
                int status_interval=100,
                bool RJ= false, bool use_mscale=true, bool use_bscale=true, bool b_half_normal=true,
                double trt_init = 1.0)
{



  bool randeff = true;
  if(random_var_ix.n_elem == 1) {
    randeff = false;
  }

  if(randeff) Rcout << "Using random effects." << std::endl;

  std::string treef_name = as<std::string>(treef_name_);
  std::ofstream treef(treef_name.c_str());

  RNGScope scope;
  RNG gen; //this one random number generator is used in all draws




  bool const_theta = false;

  if(theta==0.){
	const_theta=false;
    theta=1.;
  }









  /*****************************************************************************
  /* Read, format y
  *****************************************************************************/
  std::vector<double> y; //storage for y
  double miny = INFINITY, maxy = -INFINITY;
  sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.

  for(NumericVector::iterator it=y_.begin(); it!=y_.end(); ++it) {
    y.push_back(*it);
    if(*it<miny) miny=*it;
    if(*it>maxy) maxy=*it;
    allys.sy += *it; // sum of y
    allys.sy2 += (*it)*(*it); // sum of y^2
  }
  size_t n = y.size();
  allys.n = n;

  double ybar = allys.sy/n; //sample mean
  double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation





  /*****************************************************************************
  /* Read, format X_con
  *****************************************************************************/
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> x_con;
  for(NumericVector::iterator it=x_con_.begin(); it!= x_con_.end(); ++it) {
    x_con.push_back(*it);
  }
  size_t p_con = x_con.size()/n;

  Rcout << "Using " << p_con << " control variables." << std::endl;

  //x cutpoints
  xinfo xi_con;

  xi_con.resize(p_con);
  for(int i=0; i<p_con; ++i) {
    NumericVector tmp = x_con_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_con[i] = tmp2;
  }





  /*****************************************************************************
  /* Read, format X_mod
  *****************************************************************************/
  int ntrt = 0;
  for(size_t i=0; i<n; ++i) {
    if(z_[i]>0) ntrt += 1;
  }
  std::vector<double> x_mod;
  for(NumericVector::iterator it=x_mod_.begin(); it!= x_mod_.end(); ++it) {
    x_mod.push_back(*it);
  }
  size_t p_mod = x_mod.size()/n;

  Rcout << "Using " << p_mod << " potential effect moderators." << std::endl;

  //x cutpoints
  xinfo xi_mod;

  xi_mod.resize(p_mod);
  for(int i=0; i<p_mod; ++i) {
    NumericVector tmp = x_mod_info_list[i];
    std::vector<double> tmp2;
    for(size_t j=0; j<tmp.size(); ++j) {
      tmp2.push_back(tmp[j]);
    }
    xi_mod[i] = tmp2;
  }




  /*****************************************************************************
  /* Read, format col_ind_mu_
  *****************************************************************************/
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> col_ind_mu_;
  for(NumericVector::iterator it=weights_mu.begin(); it!= weights_mu.end(); ++it) {
    col_ind_mu_.push_back(*it);
  }




  /*****************************************************************************
  /* Read, format col_ind_tau_
  *****************************************************************************/
  //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
  std::vector<double> col_ind_tau_;
  for(NumericVector::iterator it=weights_tau.begin(); it!= weights_tau.end(); ++it) {
    col_ind_tau_.push_back(*it);
  }






  /*****************************************************************************
  /* Setup the model
  *****************************************************************************/

  //--------------------------------------------------
  //trees
  //make trt_init a parameter later
  std::vector<tree> t_mod(ntree_mod);
  for(size_t i=0;i<ntree_mod;i++) t_mod[i].setm(trt_init/(double)ntree_mod);

  std::vector<tree> t_con(ntree_con);
  for(size_t i=0;i<ntree_con;i++) t_con[i].setm(ybar/(double)ntree_con);



  //--------------------------------------------------
  //prior parameters
  // PX scale parameter for b: b(x) = bscale*b_0(x), bscale ~ N(0,1/bscale_prec)
  double bscale_prec = 2;
  double bscale0 = -0.5;
  double bscale1 = 0.5;

  double mscale_prec = 1.0;
  double mscale = 1.0;
  double delta_con = 1.0;
  double delta_mod = 1.0;

  pinfo pi_mod;
  pi_mod.pbd = 1.0; //prob of birth/death move
  pi_mod.pb = .5; //prob of birth given  birth/death

  pi_mod.alpha = mod_alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
  pi_mod.beta  = mod_beta;  //2 for bart means it is harder to build big trees.
  pi_mod.tau   = con_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod)); //sigma_mu, variance on leaf parameters
  pi_mod.sigma = shat; //resid variance is \sigma^2_y/bscale^2 in the backfitting update

  pinfo pi_con;
  pi_con.pbd = 1.0; //prob of birth/death move
  pi_con.pb = .5; //prob of birth given  birth/death

  pi_con.alpha = con_alpha;
  pi_con.beta  = con_beta;
  pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con)); //sigma_mu, variance on leaf parameters

  pi_con.sigma = shat/fabs(mscale); //resid variance in backfitting is \sigma^2_y/mscale^2

  double sigma = shat;




  //--------------------------------------------------
  //dinfo for control function m(x)
  double* allfit_con = new double[n]; //sum of fit of all trees
  for(size_t i=0;i<n;i++) allfit_con[i] = ybar;
  double* r_con = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_con;
  di_con.n=n; di_con.p=p_con; di_con.x = &x_con[0]; di_con.y=r_con; //the y for each draw will be the residual



  // Initialize Splitting probabilities and counts for m(x)
   std::vector<size_t> nv_con;
   std::vector<double> pv_con;

   for(size_t j=0;j<p_con;j++){
     nv_con.push_back(0);
     pv_con.push_back(1/(double)p_con);
   }

   std::vector<double> lpv_con;




  //--------------------------------------------------
  //dinfo for trt effect function b(x)
  double* allfit_mod = new double[n]; //sum of fit of all trees
  for(size_t i=0;i<n;i++) allfit_mod[i] = (z_[i]*bscale1 + (1-z_[i])*bscale0)*trt_init;
  double* r_mod = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
  dinfo di_mod;
  di_mod.n=n; di_mod.p=p_mod; di_mod.x = &x_mod[0]; di_mod.y=r_mod; //the y for each draw will be the residual





  // Initialize Splitting probabilities and counts for b(x)
   std::vector<size_t> nv_mod;
   std::vector<double> pv_mod;

   for(size_t j=0;j<p_mod;j++){
     nv_mod.push_back(0);
     pv_mod.push_back(1/(double)p_mod);
   }


   std::vector<double> lpv_mod;







  //--------------------------------------------------
  //dinfo and design for trt effect function out of sample
  //x for predictions
  dinfo di_mod_est; //data information for prediction
  std::vector<double> x_mod_est;     //stored like x
  size_t n_mod_est;
  //  if(x_mod_est_.size()) {
  for(NumericVector::iterator it=x_mod_est_.begin(); it!=x_mod_est_.end(); ++it) {
    x_mod_est.push_back(*it);
  }
  n_mod_est = x_mod_est.size()/p_mod;
  if(x_mod_est.size() != n_mod_est*p_mod) stop("error, wrong number of elements in effect estimate data set\n");
  //if(n_mod_est)
  di_mod_est.n=n_mod_est; di_mod_est.p=p_mod; di_mod_est.x = &x_mod_est[0]; di_mod_est.y=0; //there are no y's!


 //  }
  //--------------------------------------------------
  //storage for ouput
  //in sample fit
  /*
  //out of sample fit
  double* ppredmean=0; //posterior mean for prediction
  double* fpredtemp=0; //temporary fit vector to compute prediction
  if(dip.n) {
  ppredmean = new double[dip.n];
  fpredtemp = new double[dip.n];
  for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
  }
  //for sigma draw
  double rss, restemp;
  */






  //--------------------------------------------------
  //setup for random effects
  size_t random_dim = random_des.n_cols;
  int nr=1;
  if(randeff) nr = n;

  arma::vec r(nr); //working residuals
  arma::vec Wtr(random_dim); // W'r

  arma::mat WtW = random_des.t()*random_des; //W'W
  arma::mat Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

  // PX parameters
  arma::vec eta(random_var_ix.n_cols); //random_var_ix is num random effects by num variance components
  eta.fill(1.0);

  for(size_t k=0; k<nr; ++k) {
    r(k) = y[k] - allfit_con[k] - allfit_mod[k];
  }

  Wtr = random_des.t()*r;
  arma::vec gamma = solve(WtW/(sigma*sigma)+Sigma_inv_random, Wtr/(sigma*sigma));
  arma::vec allfit_random = random_des*gamma;
  if(!randeff) allfit_random.fill(0);

  //--------------------------------------------------
  //storage for the fits
  double* allfit = new double[n]; //yhat
  for(size_t i=0;i<n;i++) {
    allfit[i] = allfit_mod[i] + allfit_con[i];
    if(randeff) allfit[i] += allfit_random[i];
  }
  double* ftemp  = new double[n]; //fit of current tree

  double* precs = new double[n]; // temp storage for conditional ''precisions'' in heteroskedastic updates

  NumericVector sigma_post(nd);
  NumericVector msd_post(nd);
  NumericVector bsd_post(nd);
  NumericVector mybscale1(nd);
  NumericVector mybscale0(nd);
  NumericMatrix m_post(nd,n);
  NumericMatrix yhat_post(nd,n);
  NumericMatrix b_post(nd,n);
  NumericMatrix b_est_post(nd,n_mod_est);
  arma::mat gamma_post(nd,gamma.n_elem);
  arma::mat random_var_post(nd,random_var.n_elem);




	// Save DART posterior inclusion prob and count
	// for m(x)
	IntegerMatrix varcnt_con(nd,p_con);
	NumericMatrix varprb_con(nd,p_con);

	// for b(x)
	IntegerMatrix varcnt_mod(nd,p_mod);
	NumericMatrix varprb_mod(nd,p_mod);




  // Write tree to file
  treef << xi_mod << std::endl; //cutpoints
  treef << ntree_mod << std::endl;  //number of trees
  treef << p_mod << std::endl;  //dimension of x's
  treef << (int)(nd / thin) << std::endl;








  //*****************************************************************************
  /* MCMC
   * note: the allfit objects are all carrying the appropriate scales
   */
  //*****************************************************************************
  Rcout << "\nBeginning MCMC:\n";

  time_t tp;
  int time1 = time(&tp);

  size_t save_ctr = 0;


  for(size_t i=0;i<(nd*thin+burn);i++) {

    if(i%status_interval==0) {
      Rcout << "iteration: " << i << " sigma: "<< sigma << endl;
    }



    //draw trees for m(x)
    for(size_t j=0;j<ntree_con;j++) {

      fit(t_con[j],xi_con,di_con,ftemp);

      for(size_t k=0;k<n;k++) {

        if(ftemp[k] != ftemp[k]) {

          Rcout << "control tree " << j <<" obs "<< k<<" "<< endl;
          Rcout << t_con[j] << endl;
          stop("nan in ftemp");

        }

        allfit[k] = allfit[k]-mscale*ftemp[k];
        allfit_con[k] = allfit_con[k]-mscale*ftemp[k];
        r_con[k] = (y[k]-allfit[k])/mscale;

        if(r_con[k] != r_con[k]) {

          Rcout << (y[k]-allfit[k]) << endl;
          Rcout << mscale << endl;
          Rcout << r_con[k] << endl;
          stop("NaN in resid");

        }
      }


      bd(t_con[j],xi_con,di_con,pi_con,nv_con,pv_con,gen);
      drmu(t_con[j],xi_con,di_con,pi_con,gen);

      fit(t_con[j],xi_con,di_con,ftemp);

      for(size_t k=0;k<n;k++) {

        allfit[k] += mscale*ftemp[k];
        allfit_con[k] += mscale*ftemp[k];

      }
    }



	// Update splitting probabilities from Dirichlet for m(x)
	if(i>=(burn/2)&&dart) {
		draw_s(nv_con,lpv_con,theta,col_ind_mu_,gen);
		draw_theta0(const_theta,theta,lpv_con,a,b,rho,gen);
		for(size_t j=0;j<p_con;j++) pv_con[j]=::exp(lpv_con[j]);
	}






    //draw trees for b(x)
    for(size_t k=0;k<ntrt;k++) {
      precs[k] = bscale1*bscale1/(sigma*sigma);
    }
    for(size_t k=ntrt;k<n;k++) {
      precs[k] = bscale0*bscale0/(sigma*sigma);
    }


    for(size_t j=0;j<ntree_mod;j++) {

      fit(t_mod[j],xi_mod,di_mod,ftemp);

      for(size_t k=0;k<n;k++) {

        if(ftemp[k] != ftemp[k]) {

          Rcout << "moderator tree " << j <<" obs "<< k<<" "<< endl;
          Rcout << t_mod[j] << endl;
          stop("nan in ftemp");

        }

        double bscale = (k<ntrt) ? bscale1 : bscale0;
        allfit[k] = allfit[k]-bscale*ftemp[k];
        allfit_mod[k] = allfit_mod[k]-bscale*ftemp[k];
        r_mod[k] = (y[k]-allfit[k])/bscale;

      }


      bdhet(t_mod[j],xi_mod,di_mod,precs,pi_mod,nv_mod,pv_mod,gen);
      drmuhet(t_mod[j],xi_mod,di_mod,precs,pi_mod,gen);

      fit(t_mod[j],xi_mod,di_mod,ftemp);

      for(size_t k=0;k<ntrt;k++) {
        allfit[k] += bscale1*ftemp[k];
        allfit_mod[k] += bscale1*ftemp[k];
      }

      for(size_t k=ntrt;k<n;k++) {
        allfit[k] += bscale0*ftemp[k];
        allfit_mod[k] += bscale0*ftemp[k];
      }

    }



	// Update splitting probabilities from Dirichlet for b(x)
	if(i>=(burn/2)&&dart) {
		draw_s(nv_mod,lpv_mod,theta,col_ind_tau_,gen);
		draw_theta0(const_theta,theta,lpv_mod,a,b,rho_tau, gen);
		for(size_t j=0;j<p_mod;j++) pv_mod[j]=::exp(lpv_mod[j]);
	}





    //update bscale0, bscale1, the PX parameters
    // bscale ~ N( [sum_i 1/w_i + bscale_prec]^{-1} \sum_i r_i/w_i, [sum_i 1/w_i + bscale_prec]^{-1})
    // where w_i = \sigma^2_y/b(x_i)^2
    //       r_i = (y_i-m(x_i))/b(x)
    // sums run over i:z_i=0 and i:z_i=1

    if(use_bscale) {
      double ww0 = 0.0, ww1 = 0.;
      double rw0 = 0.0, rw1 = 0.;
      double s2 = sigma*sigma;
      for(size_t k=0; k<n; ++k) {
        double bscale = (k<ntrt) ? bscale1 : bscale0;
        double w = s2*bscale*bscale/(allfit_mod[k]*allfit_mod[k]);

        if(w!=w) {
          Rcout << " w " << w << endl;
          stop("");
        }

        double randeff_contrib = randeff ? allfit_random[k] : 0.0;

        double r = (y[k] - allfit_con[k] - randeff_contrib)*bscale/allfit_mod[k];

        if(r!=r) {
          Rcout << "bscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
          stop("");
        }
        if(k<ntrt) {
          ww1 += 1/w;
          rw1 += r/w;
        } else {
          ww0 += 1/w;
          rw0 += r/w;
        }
      }

      double bscale1_old = bscale1;
      double bscale_fc_var = 1/(ww1 + bscale_prec);
      bscale1 = bscale_fc_var*rw1 + gen.normal(0., 1.)*sqrt(bscale_fc_var);

      double bscale0_old = bscale0;
      bscale_fc_var = 1/(ww0 + bscale_prec);
      bscale0 = bscale_fc_var*rw0 + gen.normal(0., 1.)*sqrt(bscale_fc_var);


//      Rcout << "bscale pars " << " " << 1/ww0 << " " << 1/ww1 << endl;

      for(size_t k=0; k<ntrt; ++k) {
        allfit_mod[k] = allfit_mod[k]*bscale1/bscale1_old;
      }
      for(size_t k=ntrt; k<n; ++k) {
        allfit_mod[k] = allfit_mod[k]*bscale0/bscale0_old;
      }


      // update delta_mod
      if(b_half_normal) {
        double ssq = 0.0;
        tree::npv bnv;
        typedef tree::npv::size_type bvsz;
        double endnode_count = 0.0;

        for(size_t j=0;j<ntree_mod;j++) {
          bnv.clear();
          t_mod[j].getbots(bnv);
          bvsz nb = bnv.size();
          for(bvsz ii = 0; ii<nb; ++ii) {
            double mm = bnv[ii]->getm(); //node parameter
            ssq += mm*mm/(pi_mod.tau*pi_mod.tau);
            endnode_count += 1.0;
          }
        }
        delta_mod = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));
      }
      pi_mod.tau   = con_sd/(sqrt(delta_mod)*sqrt((double) ntree_mod));



    } else {
      bscale0 = -0.5;
      bscale1 =  0.5;
    }
    pi_mod.sigma = sigma;


    //update mscale, the PX parameter
    // bscale ~ N( [sum_i 1/w_i + bscale_prec]^{-1} \sum_i r_i/w_i, [sum_i 1/w_i + bscale_prec]^{-1})
    // where w_i = \sigma^2_y/b(x_i)^2
    //       r_i = (y_i-m(x_i))/b(x)

    if(use_mscale) {
      double ww = 0.;
      double rw = 0.;
      double s2 = sigma*sigma;
      for(size_t k=0; k<n; ++k) {
        double w = s2*mscale*mscale/(allfit_con[k]*allfit_con[k]);
        if(w!=w) {
          Rcout << " w " << w << endl;
          stop("");
        }

        double randeff_contrib = randeff ? allfit_random[k] : 0.0;

        double r = (y[k] - allfit_mod[k]- randeff_contrib)*mscale/allfit_con[k];
        if(r!=r) {
          Rcout << "mscale " << k << " r " << r << " mscale " <<mscale<< " b*z " << allfit_mod[k]*z_[k] << " bscale " << bscale0 << " " <<bscale1 << endl;
          stop("");
        }
        ww += 1/w;
        rw += r/w;
      }

      double mscale_old = mscale;
      double mscale_fc_var = 1/(ww + mscale_prec);
      mscale = mscale_fc_var*rw + gen.normal(0., 1.)*sqrt(mscale_fc_var);

      //Rcout<< mscale_fc_var << " " << rw <<" " << mscale << endl;

      for(size_t k=0; k<n; ++k) {
        allfit_con[k] = allfit_con[k]*mscale/mscale_old;
      }

      // update delta_con

      double ssq = 0.0;
      tree::npv bnv;
      typedef tree::npv::size_type bvsz;
      double endnode_count = 0.0;

      for(size_t j=0;j<ntree_con;j++) {
        bnv.clear();
        t_con[j].getbots(bnv);
        bvsz nb = bnv.size();
        for(bvsz ii = 0; ii<nb; ++ii) {
          double mm = bnv[ii]->getm(); //node parameter
          ssq += mm*mm/(pi_con.tau*pi_con.tau);
          endnode_count += 1.0;
        }
      }

      delta_con = gen.gamma(0.5*(1. + endnode_count), 1.0)/(0.5*(1 + ssq));

      pi_con.tau   = con_sd/(sqrt(delta_con)*sqrt((double) ntree_con));

    } else {
      mscale = 1.0;
    }
    pi_con.sigma = sigma/fabs(mscale); //should be sigma/abs(mscale) for backfitting

    //sync allfits after scale updates, if necessary. Could do smarter backfitting updates inline
    if(use_mscale || use_bscale) {
      for(size_t k=0; k<n; ++k) {
        double randeff_contrib = randeff ? allfit_random[k] : 0.0;
        allfit[k] = allfit_con[k] + allfit_mod[k] + randeff_contrib;
      }
    }


    if(randeff) {
      //update random effects
      for(size_t k=0; k<n; ++k) {
        r(k) = y[k] - allfit_con[k] - allfit_mod[k];
        allfit[k] -= allfit_random[k];
      }

      Wtr = random_des.t()*r;

      arma::mat adj = diagmat(random_var_ix*eta);
      //    Rcout << adj << endl << endl;
      arma::mat Phi = adj*WtW*adj/(sigma*sigma) + Sigma_inv_random;
      arma::vec m = adj*Wtr/(sigma*sigma);
      //Rcout << m << Phi << endl << Sigma_inv_random;
      gamma = rmvnorm_post(m, Phi);

      //Rcout << "updated gamma";

      // Update px parameters eta

      arma::mat adj2 = diagmat(gamma)*random_var_ix;
      arma::mat Phi2 = adj2.t()*WtW*adj2/(sigma*sigma) + arma::eye(eta.size(), eta.size());
      arma::vec m2 = adj2.t()*Wtr/(sigma*sigma);
      //Rcout << m << Phi << endl << Sigma_inv_random;
      eta = rmvnorm_post(m2, Phi2);

      //Rcout << "updated eta";

      // Update variance parameters

      arma::vec ssqs   = random_var_ix.t()*(gamma % gamma);
      //Rcout << "A";
      arma::rowvec counts = sum(random_var_ix, 0);
      //Rcout << "B";
      for(size_t ii=0; ii<random_var_ix.n_cols; ++ii) {
        random_var(ii) = 1.0/gen.gamma(0.5*(random_var_df + counts(ii)), 1.0)*2.0/(random_var_df + ssqs(ii));
      }
      //Rcout << "updated vars" << endl;
      Sigma_inv_random = diagmat(1/(random_var_ix*random_var));

      allfit_random = random_des*diagmat(random_var_ix*eta)*gamma;

      //Rcout << "recom allfit vars" << endl;

      for(size_t k=0; k<n; ++k) {
        allfit[k] = allfit_con[k] + allfit_mod[k] + allfit_random(k); //+= allfit_random[k];
      }
    }




    //draw sigma
    double rss = 0.0;
    double restemp = 0.0;

    for(size_t k=0;k<n;k++) {
      restemp = y[k]-allfit[k];
      rss += restemp*restemp;
    }

    sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
    pi_con.sigma = sigma/fabs(mscale);
    pi_mod.sigma = sigma;




	// Save MCMC output
    if( ((i>=burn) & (i % thin==0)) )  {
		for (size_t j = 0; j < ntree_mod; j++) treef << t_mod[j] << std::endl;

      msd_post(save_ctr) = fabs(mscale)*con_sd;
      bsd_post(save_ctr) = fabs(bscale1-bscale0)*mod_sd;

	  mybscale1(save_ctr) = bscale1;
	  mybscale0(save_ctr) = bscale0;

      gamma_post.row(save_ctr) = (diagmat(random_var_ix*eta)*gamma).t();
      random_var_post.row(save_ctr) = (sqrt( eta % eta % random_var)).t();

      sigma_post(save_ctr) = sigma;
      for(size_t k=0;k<n;k++) {
        m_post(save_ctr, k) = allfit_con[k];
        yhat_post(save_ctr, k) = allfit[k];
      }
      for(size_t k=0;k<n;k++) {
        double bscale = (k<ntrt) ? bscale1 : bscale0;
        b_post(save_ctr, k) = (bscale1-bscale0)*allfit_mod[k]/bscale;
      }
      for(size_t k=0;k<di_mod_est.n;k++) {
        b_est_post(save_ctr, k) = (bscale1-bscale0)*fit_i(k, t_mod, xi_mod, di_mod_est);
      }


	  // Saving Splitting Dirichlet Probabilities
	  // For m(x)
	  for(size_t k=0;k<p_con;k++) {
		  varcnt_con(save_ctr, k) = nv_con[k];
		  varprb_con(save_ctr, k) = pv_con[k];
	  }

	  // For b(x)
	  for(size_t k=0;k<p_mod;k++) {
		  varcnt_mod(save_ctr,k) = nv_mod[k];
		  varprb_mod(save_ctr, k) = pv_mod[k];
	  }


      save_ctr += 1;


    }


  }


  int time2 = time(&tp);
  Rcout << "time for loop: " << time2 - time1 << endl;

  t_mod.clear(); t_con.clear();
  delete[] allfit;
  delete[] allfit_mod;
  delete[] allfit_con;
  delete[] r_mod;
  delete[] r_con;
  delete[] ftemp;

  treef.close();

  return(List::create(_["yhat_post"] = yhat_post, _["b_post"] = b_post, _["b_est_post"] = b_est_post,
                      _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post,
                      _["gamma"] = gamma_post, _["random_var_post"] = random_var_post, _["bscale1"] = mybscale1, _["bscale0"] = mybscale0,
					  _["varcnt_con"] = varcnt_con, _["varcnt_mod"] = varcnt_mod, _["varprb_con"] = varprb_con, _["varprb_mod"] = varprb_mod
  ));

}
