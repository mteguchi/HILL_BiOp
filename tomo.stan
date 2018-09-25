data {
int<lower=1> N;
int<lower=1> n_notna;
real y[n_notna];
int<lower=1> notna[n_notna];
int<lower=1> z[N];
}
parameters {
  real<lower=0> prosigma[2];
  real<lower=0> obssigma[1];
  real pro_devs[N-1];
  real X0;
}
transformed parameters {
  vector[N] x;
  x[1] = X0;
  for(i in 1:(N-1)) {
  	x[i+1] = x[i] + pro_devs[i];
  }
 }
model {
	prosigma ~ cauchy(0,3);
	obssigma ~ cauchy(0,3);
	X0 ~ normal(0, 1);
	for(i in 1:(N-1)) {
		pro_devs[i] ~ normal(0, prosigma[z[i]]);
	}
	
	for(i in 1:n_notna) {
		y[i] ~ normal(x[notna[i]], obssigma);
	}
}
generated quantities {
  // log_lik is for use with the loo package
  vector[n_notna] log_lik;
  for(i in 1:n_notna) {
  	log_lik[i] = normal_lpdf(y[i] | x[notna[i]], obssigma);
  }
}