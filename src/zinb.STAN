data {
  // number of observations
  int<lower=0> N;
  // number of predictors
  int<lower=1> K;
  // number of samples
  int<lower=2> S;
  // outcome counts
  int<lower=0> y[N];
  // predictors
  matrix[N, K] X;
  // for random effects
  int<lower=1,upper=S> s[N];
  // depth info for ZINB
  real d[S];
  // mean for zi prior
  real zi_mean;
  // sd for zi prior
  real<lower=0> zi_sd;
  // sd for coefficient prior
  real <lower=0> b_sd;
}

parameters {
  // intercept
  real beta0;
  // coefficients
  vector[K] beta;
  // dispersion
  real<lower=0> rphi;
  // random effects
  vector[S] alpha;
  // standard deviation for RE
  real sigma;
  // intercept for zi
  real beta0zi;
  // coefficient for zi
  real betazi;
}

transformed parameters  {
  real mu[N];
  real phi;
  real theta[S];
  for (n in 1:N) {
    real lpred = 0;
    for (k in 1:K) {
      lpred += beta[k] * X[n,k];
      }
    mu[n] = exp(lpred + beta0 + alpha[s[n]]);
  }
  for (j in 1:S) {
    theta[j] = inv_logit(beta0zi + betazi * d[j]);
  }
  phi = 1/rphi;
}

model {
  // priors on the regression parameters
  beta0 ~ normal(0, 3);
  beta ~ normal(0, b_sd);
  rphi ~ normal(0, 25);
  
  // RE part
  alpha ~ normal(0, sigma);
  sigma ~ normal(0,5);
  
  // ZI priors
  beta0zi ~ normal(0,3);
  betazi ~ normal(zi_mean,zi_sd);
  
  // error
  for (n in 1:N) {
  	if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta[s[n]]), bernoulli_lpmf(0 | theta[s[n]]) + neg_binomial_2_lpmf(y[n] | mu[n], phi));  	
    else
      target += bernoulli_lpmf(0 | theta[s[n]]) + neg_binomial_2_lpmf(y[n] | mu[n], phi);
      }
}

generated quantities {
  // log likelihood for loo
  real log_lik[N];
  
  for (n in 1:N) {
  	if (y[n] == 0)
  		log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[s[n]]), bernoulli_lpmf(0 | theta[s[n]]) + neg_binomial_2_lpmf(y[n] | mu[n], phi));
  else
	log_lik[n] = bernoulli_lpmf(0 | theta[s[n]]) + neg_binomial_2_lpmf(y[n] | mu[n], phi);
    }

}

