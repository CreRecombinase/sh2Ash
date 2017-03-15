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
	vector<lower=0>[2] ptilde; //beta hyper parameters
}
parameters {
	//real z; // z = arctanh(r)
    //real<lower=0> k; // k = sqrt(1-rho^2)tilde_sigma_2/sigma_2
    //real logit_pb2_0_b1_n0;
    real<lower=0,upper=1> pb2_0_b1_n0;
}
transformed parameters {
    vector[4] alpha;
    real pb2_0_b1_0;
    //real pb2_0_b1_n0;
    //pb2_0_b1_n0 = exp(logit_pb2_0_b1_n0)/(1 + exp(logit_pb2_0_b1_n0));
    alpha[2] = pb2_0_b1_n0 * (1-pb1_0);
    alpha[1] = pb2_0 - alpha[2];
    alpha[3] = pb1_0 - alpha[1];
    alpha[4] = 1-alpha[1]-alpha[2]-alpha[3];
    pb2_0_b1_0 = alpha[1]/pb1_0;
}

model {
	//pb2_0_b1_n0 ~ beta(ptilde[1], ptilde[2]);
	target += beta_lpdf(pb2_0_b1_n0 | ptilde[1], ptilde[2]);
	for(n in 1:N)
		target += log_mix(pb1_0,
			log_mix(pb2_0_b1_0,
			    normal_lpdf(beta_hat_1[n]| 0, seb1[n]) + //alpha1
                    normal_lpdf(beta_hat_2[n]| 0, seb2[n]),
				normal_lpdf(beta_hat_1[n]| 0,  seb1[n]) + //alpha3
                    normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 + seb2[n]^2))),
			log_mix(pb2_0_b1_n0,
				normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 +seb1[n]^2))+  //alpha2
				  normal_lpdf(beta_hat_2[n]| 0, seb2[n]),
				normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 + seb1[n]^2)) + //alpha4
                  normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 + seb2[n]^2))
			));

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
            log_mix(pb2_0_b1_n0,
                normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 +seb1[n]^2))+  //alpha2
                  normal_lpdf(beta_hat_2[n]| 0, seb2[n]),
                normal_lpdf(beta_hat_1[n]| 0, sqrt(sigma_1^2 + seb1[n]^2)) + //alpha4
                  normal_lpdf(beta_hat_2[n]| 0, sqrt(sigma_2^2 + seb2[n]^2))
            ));
}

