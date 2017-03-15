data {
	int<lower=0> N; // Number of SNPs
	vector[N] beta_hat_1; //effect size estimates
	vector<lower=0>[N] b1; 	//Standard errors
	vector[N] beta_hat_2; //effect size estimates
	vector<lower=0>[N] b2; 	//Standard errors
	real<lower=0> sigma_1;
	real<lower=0> sigma_2;
	real<lower=0,upper=1> pb1_0;
	real<lower=0,upper=1> pb2_0;
}


model {
	for(n in 1:N)
		target += log_mix(pb1_0, 
                    normal_lpdf(beta_hat_1[n]| 0, b1[n]), 
                    normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 +b1[n]^2))) + 
                log_mix(pb2_0, 
                    normal_lpdf(beta_hat_2[n]| 0, b2[n]), 
                    normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 +b2[n]^2)));

}

generated quantities {
    vector[N] log_lik;
    for(n in 1:N)
        log_lik[n] = log_mix(pb1_0, 
                normal_lpdf(beta_hat_1[n]| 0, b1[n]), 
                normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 +b1[n]^2))) + 
             log_mix(pb2_0,
                normal_lpdf(beta_hat_2[n]| 0, b2[n]),
                normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 +b2[n]^2)));
}

