data {
	int<lower=0> N; // Number of SNPs
	vector[N] beta_hat_1; //effect size estimates
	vector<lower=0>[N] seb1; 	//Standard errors
	vector[N] beta_hat_2; //effect size estimates
	vector<lower=0>[N] seb2; 	//Standard errors
	real<lower=0> sigma_1;
	real<lower=0> sigma_2;
	real<lower=0,upper=1> pb1_0;
	real<lower=0,upper=1> pb2_0;
}
parameters {
    //pb2_0_b1_n0 = 0
	real z; // z = arctanh(r)
  real<lower=0> k; // k = sqrt(1-rho^2)*tilde_sigma_2/sigma_2
}
transformed parameters {
	real rho;
    real tilde_sigma_2;
    vector[4] alpha;
    real pb2_0_b1_0;
    real tau;

	rho = tanh(z);
    tilde_sigma_2 = k*sigma_2/sqrt(1-rho^2);
    alpha[2] = 0;
    alpha[1] = pb2_0 - alpha[2];
    alpha[3] = pb1_0 - alpha[1];
    alpha[4] = 1-alpha[1]-alpha[2]-alpha[3];
    pb2_0_b1_0 = alpha[1]/pb1_0;
    tau = tilde_sigma_2*sqrt(1-rho^2);
}

model {
	z ~ normal(0, 0.5);
    k ~ chi_square(1);
	for(n in 1:N)
		target += log_mix(pb1_0,
			log_mix(pb2_0_b1_0,
			    normal_lpdf(beta_hat_1[n]| 0, seb1[n]) + //alpha1
                    normal_lpdf(beta_hat_2[n]| 0, seb2[n]),
				normal_lpdf(beta_hat_1[n]| 0,  seb1[n]) + //alpha3
                    normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 + seb2[n]^2))),
			normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 + seb1[n]^2)) + //alpha4
                  normal_lpdf(beta_hat_2[n]|
                     (beta_hat_1[n]*rho*sigma_1*tilde_sigma_2)/(sigma_1^2+ seb1[n]^2),
                     sqrt((tilde_sigma_2^2 + seb2[n]^2)-(((rho*sigma_1*tilde_sigma_2)^2)/(sigma_1^2 + seb1[n]^2))))
			);

}


generated quantities {
    vector[N] log_lik;
    for(n in 1:N)
        log_lik[n] = log_mix(pb1_0,
			log_mix(pb2_0_b1_0,
			    normal_lpdf(beta_hat_1[n]| 0, seb1[n]) + //alpha1
                    normal_lpdf(beta_hat_2[n]| 0, seb2[n]),
				normal_lpdf(beta_hat_1[n]| 0,  seb1[n]) + //alpha3
                    normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 + seb2[n]^2))),
			normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 + seb1[n]^2)) + //alpha4
                  normal_lpdf(beta_hat_2[n]|
                     (beta_hat_1[n]*rho*sigma_1*tilde_sigma_2)/(sigma_1^2+ seb1[n]^2),
                     sqrt((tilde_sigma_2^2 + seb2[n]^2)-(((rho*sigma_1*tilde_sigma_2)^2)/(sigma_1^2 + seb1[n]^2))))
			);

}

