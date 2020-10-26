Rcpp::loadModule(module = "TreeSamples", TRUE)


.ident <- function(...){
  # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
  args <- c(...)
  if( length( args ) > 2L ){
    #  recursively call ident()
    out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
  }else{
    out <- identical( args[1] , args[2] )
  }
  return( all( out ) )
}

.cp_quantile = function(x, num=10000, cat_levels=8){
  nobs = length(x)
  nuniq = length(unique(x))

  if(nuniq==1) {
    ret = x[1]
    warning("A supplied covariate contains a single distinct value.")
  } else if(nuniq < cat_levels) {
    xx = sort(unique(x))
    ret = xx[-length(xx)] + diff(xx)/2
  } else {
    q = approxfun(sort(x),quantile(x,p = 0:(nobs-1)/nobs))
    ind = seq(min(x),max(x),length.out=num)
    ret = q(ind)
  }

  return(ret)
}

#' Fit Sparse Bayesian Causal Forests
#'
#' @references
#' Caron et. al. (2020) (TO BE FILLED)
#' \cr
#' \cr
#' Hahn, Murray, and Carvalho(2017). Bayesian regression tree models for causal inference: regularization, confounding, and heterogeneous effects.
#'  https://arxiv.org/abs/1706.09523. (Call citation("bcf") from the
#' command line for citation information in Bibtex format.)
#'
#' @details Fits Sparse version of Bayesian Causal Forest model (Hahn et. al. 2018) as in Caron et al. (2020): For a response
#' variable y, binary treatment z, and covariates x,
#' \deqn{y_i = \mu(x_i, \pi_i) + \tau(x_i, \pi_i)z_i + \epsilon_i}
#' where \eqn{\pi_i} is an (optional) estimate of the propensity score \eqn{\Pr(Z_i=1 | X_i=x_i)} and
#' \eqn{\epsilon_i \sim N(0,\sigma^2)}
#'
#' Some notes:
#' \itemize{
#'    \item x_control and x_moderate must be numeric matrices. See e.g. the makeModelMatrix function in the
#'    dbarts package for appropriately constructing a design matrix from a data.frame
#'    \item The set of arguments "inform_mu", "weights_mu", "inform_tau", "weights_tau" refers to the informative prior
#'    version of Sparse BCF. If not specified, the non-informative version of Sparse BCF is implemented
#'    \item sd_control and sd_moderate are the prior SD(mu(x)) and SD(tau(x)) at a given value of x (respectively). If
#'    use_muscale = FALSE, then this is the parameter \eqn{\sigma_\mu} from the original BART paper, where the leaf parameters
#'    have prior distribution \eqn{N(0, \sigma_\mu/m)}, where m is the number of trees.
#'    If use_muscale=TRUE then sd_control is the prior median of a half Cauchy prior for SD(mu(x)). If use_tauscale = TRUE,
#'    then sd_moderate is the prior median of a half Normal prior for SD(tau(x)).
#'    \item By default the prior on \eqn{\sigma^2} is calibrated as in Chipman, George and McCulloch (2008).
#'
#'
#' }
#' @param y Response variable
#' @param z Treatment variable
#' @param x_control Design matrix for the "prognostic" function mu(x)
#' @param x_moderate Design matrix for the covariate-dependent treatment effects tau(x)
#' @param pihat Length n estimates of
#' @param OOB Boolean for Out-of-Bag prediction of CATE (default is FALSE)
#' @param x_pred Design test-set matrix X for OOB prediction
#' @param z_pred Binary test-set treatment vector for OOB prediction
#' @param sparse Boolean for implementing Sparse BCF (default is TRUE). sparse=FALSE implements default BCF
#' @param a First parameter of the hyperprior Beta(a,b) distribution (default a=0.5)
#' @param b Second parameter of the Beta hyperprior Beta(a,b) distribution (default b=1)
#' @param inform_mu Boolean for informative prior weights in the Dirichlet on mu(x) (default is FALSE)
#' @param weights_mu Vector of informative weights for mu(x, pi(x)), last entry is the weight corresponding to pi(x). Length must be (P+1)
#' @param inform_tau Boolean for informative prior weights in the Dirichlet on tau(x)
#' @param weights_tau Vector of informative weights for tau(x). Length must be P
#' @param save_trees_dir File where trees info are saved for OOB prediction. If unspecified a temporary file is used, then deleted
#' @param nburn Number of burn-in MCMC iterations
#' @param nsim Number of MCMC iterations to save after burn-in
#' @param nthin Save every nthin'th MCMC iterate. The total number of MCMC iterations will be nsim*nthin + nburn.
#' @param update_interval Print status every update_interval MCMC iterations
#' @param ntree_control Number of trees in mu(x)
#' @param sd_control SD(mu(x)) marginally at any covariate value (or its prior median if use_muscale=TRUE)
#' @param base_control Base for tree prior on mu(x) trees (see details)
#' @param power_control Power for the tree prior on mu(x) trees
#' @param ntree_moderate Number of trees in tau(x)
#' @param sd_moderate SD(tau(x)) marginally at any covariate value (or its prior median if use_tauscale=TRUE)
#' @param base_moderate Base for tree prior on tau(x) trees (see details)
#' @param power_moderate Power for the tree prior on tau(x) trees (see details)
#' @param nu Degrees of freedom in the chisq prior on \eqn{sigma^2}
#' @param lambda Scale parameter in the chisq prior on \eqn{sigma^2}
#' @param sigq Calibration quantile for the chisq prior on \eqn{sigma^2}
#' @param sighat Calibration estimate for the chisq prior on \eqn{sigma^2}
#' @param theta Set theta sparsity parameter; zero means random
#' @param omega Set omega sparsity parameter; zero means random
#' @param rho_mu Sparse parameter for mu(x) alpha draws. Default is rho_mu=P+1, lower values increase sparsity
#' @param rho_tau Sparse parameter for tau(x) alpha draws. Default is rho_mu=ceiling(P/2), lower values increase sparsity
#' @param augment whether data augmentation is to be performed in sparse variable selection
#' @param include_pi Takes values "control", "moderate", "both" or "none". Whether to
#' include pihat in mu(x) ("control"), tau(x) ("moderate"), both or none. Values of "control"
#' or "both" are HIGHLY recommended with observational data.
#' @param use_muscale Use a half-Cauchy hyperprior on the scale of mu.
#' @param use_tauscale Use a half-Normal prior on the scale of tau.
#' @return A list with elements
#' \item{tau}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{mu}{\code{nsim} by \code{n} matrix of posterior samples of individual treatment effects}
#' \item{sigma}{Length \code{nsim} vector of posterior samples of sigma}
#' @examples
#'\donttest{
#'
#' # data generating process
#' p = 10    # control and moderating variables
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' # Higher number of MCMC iterations is recommended; here it is just 4000 for time complexity
#' SparseBCF_fit = SparseBCF(y, z, x, x, pihat, nburn=2000, nsim=2000)
#'
#' # Get posterior of treatment effects
#' tau_post = SparseBCF_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#' # Get Sparse BCF posterior splitting proabilities (variable importance) on mu(x) and tau(x)
#' SplitProb_mu = colMeans(SparseBCF_fit$varprb_mu)
#' SplitProb_tau = colMeans(SparseBCF_fit$varprb_tau)
#'
#' barplot(SplitProb_mu)
#' barplot(SplitProb_tau)
#'
#'}
#'\dontshow{
#'
#' # data generating process
#' p = 10    # control and moderating variables
#' n = 250
#' #
#' set.seed(1)
#'
#' x = matrix(rnorm(n*p), nrow=n)
#'
#' # create targeted selection
#' q = -1*(x[,1]>(x[,2])) + 1*(x[,1]<(x[,2]))
#'
#' # generate treatment variable
#' pi = pnorm(q)
#' z = rbinom(n,1,pi)
#'
#' # tau is the true (homogeneous) treatment effect
#' tau = (0.5*(x[,3] > -3/4) + 0.25*(x[,3] > 0) + 0.25*(x[,3]>3/4))
#'
#' # generate the response using q, tau and z
#' mu = (q + tau*z)
#'
#' # set the noise level relative to the expected mean function of Y
#' sigma = diff(range(q + tau*pi))/8
#'
#' # draw the response variable with additive error
#' y = mu + sigma*rnorm(n)
#'
#' # If you didn't know pi, you would estimate it here
#' pihat = pnorm(q)
#'
#' # Higher number of MCMC iterations is recommended; here it is just 4000 for time complexity
#' SparseBCF_fit = SparseBCF(y, z, x, x, pihat, nburn=10, nsim=10)
#'
#' # Get posterior of treatment effects
#' tau_post = SparseBCF_fit$tau
#' tauhat = colMeans(tau_post)
#' plot(tau, tauhat); abline(0,1)
#'
#' # Get Sparse BCF posterior splitting proabilities (variable importance) on mu(x) and tau(x)
#' SplitProb_mu = colMeans(SparseBCF_fit$varprb_mu)
#' SplitProb_tau = colMeans(SparseBCF_fit$varprb_tau)
#'
#' barplot(SplitProb_mu)
#' barplot(SplitProb_tau)
#'
#'}
#'#'
#' @useDynLib SparseBCF
#' @export
SparseBCF <-
  function(y, z, x_control, x_moderate=x_control, pihat,
           OOB = FALSE, x_pred, z_pred,
           sparse=TRUE,
           a=0.5,
           b=1,
           inform_mu = FALSE, weights_mu = NULL,
           inform_tau = FALSE, weights_tau = NULL,
           save_trees_dir,
           nburn, nsim, nthin = 1, update_interval = 1000,
           ntree_control = 200,
           sd_control = 2*sd(y),
           base_control = 0.95,
           power_control = 2,
           ntree_moderate = 50,
           sd_moderate = sd(y),
           base_moderate = 0.25,
           power_moderate = 3,
           nu = 3, lambda = NULL, sigq = .9, sighat = NULL,
           theta=0,
           omega=1,
           rho_mu=NULL,
           rho_tau = NULL,
           augment=FALSE,
           include_pi = "control", use_muscale=TRUE, use_tauscale=TRUE
) {

  pihat = as.matrix(pihat)
  if( !.ident(length(y),
              length(z),
              nrow(x_control),
              nrow(x_moderate),
              nrow(pihat)
  )
  ) {

    stop("Data size mismatch. The following should all be equal:
         length(y): ", length(y), "\n",
         "length(z): ", length(z), "\n",
         "nrow(x_control): ", nrow(x_control), "\n",
         "nrow(x_moderate): ", nrow(x_moderate), "\n",
         "nrow(pihat): ", nrow(pihat),"\n"
    )
  }

  if(any(is.na(y))) stop("Missing values in y")
  if(any(is.na(z))) stop("Missing values in z")
  if(any(is.na(x_control))) stop("Missing values in x_control")
  if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
  if(any(is.na(pihat))) stop("Missing values in pihat")

  if(any(!is.finite(y))) stop("Non-numeric values in y")
  if(any(!is.finite(z))) stop("Non-numeric values in z")
  if(any(!is.finite(x_control))) stop("Non-numeric values in x_control")
  if(any(!is.finite(x_moderate))) stop("Non-numeric values in x_moderate")
  if(any(!is.finite(pihat))) stop("Non-numeric values in pihat")

  if(!all(sort(unique(z)) == c(0,1))) stop("z must be a vector of 0's and 1's, with at least one of each")

  if(length(unique(y))<5) warning("y appears to be discrete")

  if(nburn<0) stop("nburn must be positive")
  if(nsim<0) stop("nsim must be positive")
  if(nthin<0) stop("nthin must be positive")
  if(nthin>nsim+1) stop("nthin must be < nsim")
  if(nburn<100) warning("A low (<100) value for nburn was supplied")


  p = ncol(x_control)
  if(length(rho_mu)==0) rho_mu=p
  if(length(rho_tau)==0) rho_tau=ceiling(p/2)  #  More sparsity on TAU



  # Check "informed" for Mu
  if (inform_mu == TRUE & is.null(weights_mu)) stop("informative prior for mu(x) option activated but no weights_mu provided")
  if (class(inform_mu) != "logical") stop("inform_mu must be boolean")
  if (inform_mu == TRUE & length(weights_mu) != (ncol(x_control)+1)) stop("weights_mu must be of length  ncol(x_control)+1")

  if (inform_mu == FALSE) weights_mu = rep(1, ncol(x_control)+1)



  # Check "informed" for Tau
  if (inform_tau == TRUE & is.null(weights_tau)) stop("informative prior for tau(x) option activated but no weights_tau provided")
  if (class(inform_tau) != "logical") stop("inform_tau must be boolean")
  if (inform_tau == TRUE & length(weights_tau) != ncol(x_moderate)) stop("weights_tau must be of same length as ncol(x_moderate)")

  if (inform_tau == FALSE) weights_tau = rep(1, ncol(x_moderate))




  ### TODO range check on parameters

  ###
  x_c = matrix(x_control, ncol=ncol(x_control))
  x_m = matrix(x_moderate, ncol=ncol(x_moderate))
  if(include_pi=="both" | include_pi=="control") {
    x_c = cbind(x_control, pihat)
  }
  if(include_pi=="both" | include_pi=="moderate") {
    x_m = cbind(x_moderate, pihat)
  }
  cutpoint_list_c = lapply(1:ncol(x_c), function(i) .cp_quantile(x_c[,i]))
  cutpoint_list_m = lapply(1:ncol(x_m), function(i) .cp_quantile(x_m[,i]))

  yscale = scale(y)
  sdy = sd(y)
  muy = mean(y)

  if(is.null(lambda)) {
    if(is.null(sighat)) {
      lmf = lm(yscale~z+as.matrix(x_c))
      sighat = summary(lmf)$sigma #sd(y) #summary(lmf)$sigma
    }
    qchi = qchisq(1.0-sigq,nu)
    lambda = (sighat*sighat*qchi)/nu
  }


  # Saving forest for OOB prediction
  save_trees_dir = tempfile(pattern = "forest", fileext = ".txt")


  perm = order(z, decreasing=TRUE)

  fitbcf = cSparseBCF(yscale[perm], z[perm], t(x_c[perm,]), t(x_m[perm,,drop=FALSE]), t(x_m[1,,drop=FALSE]),
                      cutpoint_list_c, cutpoint_list_m,
                      random_des = matrix(1),
                      random_var = matrix(1),
                      random_var_ix = matrix(1),
                      random_var_df = 3,
                      nburn, nsim, nthin,
                      ntree_moderate, ntree_control, lambda, nu,
                      con_sd = ifelse(abs(2*sdy - sd_control)<1e-6, 2, sd_control/sdy),
                      mod_sd = ifelse(abs(sdy - sd_moderate)<1e-6, 1, sd_moderate/sdy)/ifelse(use_tauscale,0.674,1), # if HN make sd_moderate the prior median
                      base_moderate, power_moderate, base_control, power_control,
                      a,
                      b,
                      rho_mu,
                      rho_tau,
                      theta,
                      omega,
                      weights_mu,
                      weights_tau,
                      save_trees_dir,
                      dart=sparse,
                      aug=augment,
                      status_interval = update_interval,
                      use_mscale = use_muscale,
                      use_bscale = use_tauscale,
                      b_half_normal = TRUE)

  #B = drop(fit$post_B)
  #B0 = fit$b0
  #EYs = fit$post_yhat

  #return(List::create(_["m_post"] = m_post, _["b_post"] = b_post, _["b_est_post"] = b_est_post,
  #                     _["sigma"] = sigma_post, _["msd"] = msd_post, _["bsd"] = bsd_post,
  #                     _["gamma"] = gamma_post, _["random_var_post"] = random_var_post

  m_post = muy + sdy*fitbcf$m_post[,order(perm)]
  tau_post = sdy*fitbcf$b_post[,order(perm)]
  #yhat_post = muy + sdy*fitbcf$m_post
  #yhat_post[,z==1] = yhat_post[,z==1] + fitbcf$b_post
  #yhat_post = (muy + sdy*fitbcf$m_post)
  #yhat_post[,z[perm]==1] = yhat_post[,z[perm]==1] + sdy*fitbcf$b_post
  #yhat_post = yhat_post[,order(perm)]


  if (OOB == F) {

    file.remove(save_trees_dir)

    # Returns
    return(
      list(sigma = sdy*fitbcf$sigma,
           yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
           mu  = m_post,
           tau = tau_post,
           mu_scale = fitbcf$msd*sdy,
           tau_scale = fitbcf$bsd*sdy,
           perm = perm,
           bscale1 = fitbcf$bscale1,
           bscale0 = fitbcf$bscale0,
           varcnt_mu = fitbcf$varcnt_con,
           varprb_mu = fitbcf$varprb_con,
           varcnt_tau = fitbcf$varcnt_mod,
           varprb_tau = fitbcf$varprb_mod
      )
    )
  }


  if (OOB == T) {

    perm_pred = order(z_pred, decreasing=TRUE)

    cat("\n\n###########################\n")
    cat("OOB loading forest\n")

    ts = TreeSamples$new()
    ts$load(save_trees_dir)

    insam = ts$predict(t(x_pred[perm_pred, ]))
    tau_pred = sdy*(fitbcf$bscale1 - fitbcf$bscale0)*insam
    tau_pred = tau_pred[, order(perm_pred)]

    # Returns
    return(
      list(sigma = sdy*fitbcf$sigma,
           yhat = muy + sdy*fitbcf$yhat_post[,order(perm)],
           mu  = m_post,
           tau = tau_post,
           tau_pred = tau_pred,
           mu_scale = fitbcf$msd*sdy,
           tau_scale = fitbcf$bsd*sdy,
           perm = perm,
           bscale1 = fitbcf$bscale1,
           bscale0 = fitbcf$bscale0,
           varcnt_mu = fitbcf$varcnt_con,
           varprb_mu = fitbcf$varprb_con,
           varcnt_tau = fitbcf$varcnt_mod,
           varprb_tau = fitbcf$varprb_mod
      )
    )

    rm(ts)
    file.remove(save_trees_dir)

  }


}
