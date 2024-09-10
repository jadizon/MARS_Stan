// Code for Meta - Analysis of Reconstructed Survival (MARS) (Fu 2022, https: // arxiv.org / pdf / 2201.00281.pdf)
// This model pools:
// 1. Estimates from a piecewise exponential survival model with a Poisson likelihood using (reconstructed) survival data (T1), 
// 2. HR estimates from studies without survival curves (T2)
// 3. Survival probabilities (T3)

// The model contains parameters for:
// log - HRs (b) in each interval, up to J intervals; 
// log - baseline hazards in each interval (log_lambda), up to J intervals;
// study - specific intercepts (alpha), up to K studies;
// study - specific log - HR slopes (beta), up to K studies;
// 
// The log - HR parameters are modelled with a multivariate normal prior for the time - series correlation between these parameters
// The log - baseline hazard parameters are modelled with a multivariate normal prior for the time - series correlation between these parameters

// For the T1 component; the log - likelihood is specified as:
// poisson_log_lupmf(death[i]|logoffset[i] + log_lambda[int1[i]] + alpha[study1[i]] + (b[int1[i]] + beta[study1[i]]) .* z1[i]) * count[i];
// logoffset[i] denotes the porportion of the time a given observation contributes to an interval. 
// e.g., for interval length of 5 years; a study with 3 years followup would have log(3 / 5) as the offset for interval 1.
// a study with 8 years follow - up would have log(5 / 5) as the offset for interval 1, and log(3 / 5) as the offset for interval 2.

// For the T2 component; the log - likelihood is
// normal_lupdf(y2|m2, s2); i.e. the log - HR estimates from T2 studies are given a normal prior

// For the T3 component; the log - likelihood is
// normal_lupdf(y3|m3, s3); i.e. the logit - transformed survival probabilities from from T3 studies are given a normal prior
// 

// Priors and hyperpriors:
// The log HRs (b) is given a multivariate - normal distribution, with the vector of location hyperparameters (mu_b). The prior for the hyperparameters mu_b is normal with mean = 0 and SD = 3.
// The covariance of the MVN is generated from a cholesky correlation matrix (L_corr_b) with an LKJ prior, and a vector of error hyperparameters (b_raw) from a normal distribution with mean = 0 and SD = 1

// The baseline hazards (log_lambda) is given a multivariate - normal distribution, with the vector of location hyperparameters (mu_log_lambda). The prior for the hyperparameters mu_log_lambda is normal with mean = 0 and SD = 3.
// Its covariance of the MVN is generated from a cholesky correlation matrix (L_corr_log_lambda) with an LKJ prior, and a vector of error hyperparameters (log_lambda_raw) from a normal distribution with mean = 0 and SD = 1

// For efficiency and numerical stability, the correlation matrices are Cholesky factored.

// For sampling efficiency, non - centred parameterisation has been used for b and log_lambda, 

// For computational efficiency in calculating log - likelihoods, the T1 data set is reduced to unique rows of  [study], [treatment], [interval], [event], and [offset] variable combinations.
// the log - likelihood is computed once for each of unique combination of variables, then multiplied by the frequency of their occurences
data {
  int<lower = 1> J;         // Number if discrete intervals to model 
  int<lower = 0> J_step;    // Step sizes between interval e.g. J_step = 12 means there are 12 units (month, year, etc.) between intervals
  int<lower = 1> N1;        // Number of obervations in T1
  int<lower = 0> study_n1;  // Number of studies constituting T1
  vector[N1] x1;            // Treatment or exposure group
  array[N1] int int1;       // The interval in which the study ends e.g. if J_step = 12 and the study ended at 100 months then the interval will be set at 9 (ceiling(100 / 12))
  array[N1] int study1;     // Mapping from 1:N to 1:n_Study for reconstructed data
  array[N1] int<lower = 0> death; // outcome variable in T1 data set
  array[N1] int<lower = 0> count; // Frequency of each unique row in the T1 data set
  array[N1] real<lower = 0> time_in_interval;     // time in interval (as proportion). Persons surviving the entire duration of an interval will have a value of 1. 
                                                  // Persons who experience death within an interval will have value between 0 and 1 and is proportional to the time alive during the interval and the total duration of that interval 
                                                  //(e.g., if an interval is 12 months in duration and a person survived up to the 6th month of that interval, their time_in_interval will be 0.5 (6/12))
  
  int<lower = 0> N2;
  int<lower = 0> study_n2; 
  vector[N2] y2;
  vector[N2] s2;
  array[N2] int study2;
  array[N2] int<lower = 0> int2;
  array[N2] real<lower = 0> t2;

  int<lower = 0> N3;
  int<lower = 0> study_n3; 
  vector[N3] y3;
  vector[N3] s3;
  array[N3] int study3;
  vector[N3]  x3; 
  array[N3] int<lower = 0> int3;
  array[N3] real<lower = 0> t3;
  array[N3] real<lower = 0, upper = 1> event_prob;
  array[N3] real<lower = 0, upper = 1> event_se;
  
}


transformed data { 
  int<lower = 1> n_Study;
  n_Study = study_n1 + study_n2 + study_n3;
  vector[N1] logoffset = log(to_vector(time_in_interval));

  matrix[N2, J] t2_matrix = rep_matrix(0, N2, J);
  for (i in 1:N2){
      real temp2 = t2[i];

    for (j in 1:int2[i]){ 
      t2_matrix[i, j] = min({J_step * 1.0, temp2}) / t2[i];
      temp2 += -J_step; 
    }  
  }

}


parameters {

  real mu_alpha; //location hyper-parameteter for study-specific random intercepts
  real<lower=0> sigma_alpha; //scale hyper-parameter for random intercepts
  vector[study_n1 + study_n3] sigma_alpha_raw; 

  real mu_beta; // location hyper-parameter for study-specific random slopes
  real<lower=0> sigma_beta; // scale hyper-parameter for study-specific random slopes
  vector[n_Study] sigma_beta_raw;

  real mu_b; // location hyperparameter of interval-specific log HR
  real<lower=0> sigma_b; // scale hyper-parameter for interval-specific log HR
  vector[J] b_raw; 

  real  mu_log_lambda; // location hyperparameter of interval-specific (log) baseline hazards
  real<lower=0> sigma_log_lambda; // scale hyper-parameter for interval-specific (log) baseline hazards
  vector[J] log_lambda_raw; 

  real<lower = 0, upper=1> rho_b; //correlation parameter used in AR1 correlation for the log-hazards
  real<lower = 0, upper=1> rho_log_lambda; //correlation parameter used in AR1 correlation for the log-baseline hazards

}


transformed parameters { 

    matrix[J, J] ar1_corr_chol_b;
    matrix[J, J] ar1_corr_chol_log_lambda;

    //Direct computation of Cholesky root of an AR(1) matrix
    ar1_corr_chol_b = rep_matrix(0, J,J);
    ar1_corr_chol_log_lambda = rep_matrix(0, J,J);

    real scaling_factor_b = sqrt(1 - pow(rho_b, 2));
    real scaling_factor_log_lambda = sqrt(1 - pow(rho_log_lambda, 2));


    for(j in 1:J){
        ar1_corr_chol_b[j,1] = pow(rho_b, j-1);
        ar1_corr_chol_log_lambda[j,1] = pow(rho_log_lambda, j-1);
    }

    vector[J] v_b = scaling_factor_b * ar1_corr_chol_b[,1];
    vector[J] v_log_lambda = scaling_factor_log_lambda * ar1_corr_chol_log_lambda[,1];

    for (j in 2:J){
        ar1_corr_chol_b[j:J, j] = v_b[1:J-j+1];
        ar1_corr_chol_log_lambda[j:J, j] = v_log_lambda[1:J-j+1];
    }

    


  // non-centred reparameterisation of b and log_lambda
  vector[J] b = 5*mu_b + sigma_b * (ar1_corr_chol_b * b_raw); // interval-specific log HR parameters 
                                        // implies b~MVN(mu_b, sigma_b * [L_corr_b' * L_corr_b]); 
                                        // mu_b~normal(0, 5);
                                        // sigma_b~half-normal(0, 1)

  vector[J] log_lambda = 5*mu_log_lambda + sigma_log_lambda * (ar1_corr_chol_log_lambda * log_lambda_raw);   // interval-specific (log) baseline hazard parameters
                                                                              // implies log_lambda~MVN(mu_log_lambda, sigma_log_lambda*[L_corr_log_lambda' * L_corr_log_lambda]); 
                                                                              // mu_log_lambda~ normal(0, 5); 
                                                                              // sigma_log_lambda~half-normal(0, 1)

  //non-centred reparameterisation of alpha and beta
  vector[study_n1 + study_n3] alpha = 2*mu_alpha + sigma_alpha*sigma_alpha_raw; // study - specific intercept
                                                                    //implies alpha~N(mu_alpha, sigma_alpha)
                                                                    //mu_alpha ~ N(0,2)
                                                                    //sigma_alpha ~ half-Normal(0, 1)

  vector[n_Study] beta = 2*mu_beta + sigma_beta*sigma_beta_raw; // study - specific slopes
                                                                    //implies beta~N(mu_beta, sigma_beta)
                                                                    //mu_beta ~ N(0,2)
                                                                    //sigma_beta ~ half-Normal(0, 1)


}

model {
  mu_alpha~std_normal();
  sigma_alpha~std_normal();
  sigma_alpha_raw~std_normal();

  mu_beta~std_normal();
  sigma_beta~std_normal();
  sigma_beta_raw~std_normal();


  mu_b~std_normal();
  sigma_b~std_normal();
  b_raw~std_normal();

  mu_log_lambda~std_normal();
  sigma_log_lambda~std_normal();
  log_lambda_raw~std_normal();
  
    vector[N3] u3;
    vector[N3] m3;
    vector[N3] m3_t;
  for(i in 1:N3){
    real temp1 = 0;
    real temp2 = t3[i];
    vector[int3[i]] temp_vec;
    for (j in 1:int3[i]){
      temp_vec[j] = log_lambda[j] + alpha[study3[i]] + (beta[study3[i]] + b[j]) * x3[i];
    }
    temp1 = exp(log_sum_exp(temp_vec));
    m3[i] = -temp1 - log1m(exp(-temp1));

  }

 
  for(i in 1:N1){
    target += poisson_log_lupmf(death[i]|logoffset[i] + log_lambda[int1[i]] + alpha[study1[i]] + (b[int1[i]] + beta[study1[i]]) .* x1[i]) * count[i]; //multiplies the log-likelihood of unique rows by their frequencies (count[i])
  }
  target += normal_lupdf(y2|beta[study2] + t2_matrix * b, s2);
  target += normal_lupdf(y3|m3, s3); 
  



  
}

generated quantities { 

  vector[J] HR_x1;
    HR_x1 = exp(b + mean(beta));

  real mean_HR_x1;
    mean_HR_x1 = mean(HR_x1);

  vector[J] lambda_adj_x0;
  vector[J] lambda_adj_x1;
    lambda_adj_x0 = exp(log_lambda + mean(alpha));
    lambda_adj_x1 = exp(log_lambda + mean(alpha) + b);

  vector[J] L_x0;
  vector[J] L_x1;
  for (j in 1:J)
  {
    L_x0[j] = sum(lambda_adj_x0[1:j]);
    L_x1[j] = sum(lambda_adj_x1[1:j]);
  }

  vector[J] loglog_surv_rate_x0;

  vector[J] loglog_surv_rate_x1;
    

  for (j in 1:J){
    loglog_surv_rate_x0[j] = (log_lambda[j] + mean(alpha));
    loglog_surv_rate_x1[j] = (log_lambda[j] + mean(alpha) + b[j]  + mean(beta)) ; 
  }

  vector[J] surv_rate_x0;
  vector[J] surv_rate_x1;

  for(j in 1:J){
    surv_rate_x0[j] = exp(-exp(log_sum_exp(loglog_surv_rate_x0[1:j])));
    surv_rate_x1[j] = exp(-exp(log_sum_exp(loglog_surv_rate_x1[1:j])));
  }

  vector[J] RMST_x0;
    RMST_x0[1] = (surv_rate_x0[1] - 1) / (-exp(log_lambda[1] + mean(alpha)));

  vector[J] RMST_x1;
    RMST_x1[1] = (surv_rate_x1[1] - 1) / (-exp(log_lambda[1] + mean(alpha) + b[1] + mean(beta)));

  for (j in 2:J){
    RMST_x0[j] = (surv_rate_x0[j] - surv_rate_x0[j - 1]) / (-exp(log_lambda[j] + mean(alpha))) + RMST_x0[j - 1];
    RMST_x1[j] = (surv_rate_x1[j] - surv_rate_x1[j - 1]) / (-exp(log_lambda[j] + mean(alpha) + b[j] + mean(beta))) + RMST_x1[j - 1];
  }

  vector[J] RMST_delta;
    RMST_delta = RMST_x1 - RMST_x0;

  vector[J] RMST_ratio;
    RMST_ratio = RMST_x1 ./ RMST_x0;
}
