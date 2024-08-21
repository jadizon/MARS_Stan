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
// The log HRs (b) is given a multivariate - normal distribution, with the vector of location hyperparameters (m_b). The prior for the hyperparameters m_b is normal with mean = 0 and SD = 3.
// The covariance of the MVN is generated from a cholesky correlation matrix (L_corr_b) with an LKJ prior, and a vector of error hyperparameters (b_raw) from a normal distribution with mean = 0 and SD = 1

// The baseline hazards (log_lambda) is given a multivariate - normal distribution, with the vector of location hyperparameters (m_log_lambda). The prior for the hyperparameters m_log_lambda is normal with mean = 0 and SD = 3.
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
  vector[study_n1 + study_n3] alpha; // study - specific intercept
  vector[n_Study] beta; // study - specific log HR
  vector[J] log_lambda_raw; // standard deviation hyperparameters of interval-specific (log) baseline hazards
  vector[J] b_raw; // standard deviation hyperparameters of interval-specific log HR
  vector[J] m_b; // location hyperparameter of interval-specific log HR
  vector[J] m_log_lambda; // location hyperparameter of interval-specific (log) baseline hazards
  cholesky_factor_corr[J] L_corr_log_lambda; // cholesky factor of correlation matrix of interval-specific (log) baseline hazards
  cholesky_factor_corr[J] L_corr_b; // cholesky factor of correlation matrix of interval-specific log HR
}


transformed parameters { 
  // non-centred reparameterisation of b and log_lambda
  vector[J] b = m_b + L_corr_b * b_raw; // interval-specific log HR parameters 
                                        // implies b~MVN(m_b, [L_corr_b' * L_corr_b] * b_raw); 
                                        // m_b~normal(0, 3);
                                        // b_raw~normal(0, 1)

  vector[J] log_lambda = m_log_lambda + L_corr_log_lambda * log_lambda_raw;   // interval-specific (log) baseline hazard parameters
                                                                              // implies log_lambda~MVN(m_log_lambda, [L_corr_log_lambda' * L_corr_log_lambda] * log_lambda_raw); 
                                                                              // m_b~normal(0, 3); 
                                                                              // log_lambda_raw~normal(0, 1)

}


model {
  alpha ~ std_normal();
  beta ~ std_normal();

  m_b~normal(0, 3);
  m_log_lambda~normal(0, 3);


  b_raw~std_normal();
  log_lambda_raw~std_normal();

  L_corr_log_lambda ~ lkj_corr_cholesky(1);
  L_corr_b ~ lkj_corr_cholesky(1);
 


    vector[N3] u3;
    vector[N3] m3;
  for(i in 1:N3){
    real temp1 = 0;
    real temp2 = t3[i];
    for (j in 1:int3[i]){
      temp1 += exp(log_lambda[j] + alpha[study3[i]] + (beta[study3[i]] + b[j]) * x3[i]);
    }
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

  vector[J] log_surv_rate_x0;
    log_surv_rate_x0[1] = -exp(log_lambda[1] + mean(alpha));
  vector[J] log_surv_rate_x1;
    log_surv_rate_x1[1] = -exp(log_lambda[1] + mean(alpha) + b[1] + mean(beta));

  for (j in 2:J){
    log_surv_rate_x0[j] = (log_surv_rate_x0[j - 1])  -exp(log_lambda[j] + mean(alpha));
    log_surv_rate_x1[j] = (log_surv_rate_x1[j - 1])  -exp(log_lambda[j] + mean(alpha) + b[j] + mean(beta)) ; 
  }

  vector[J] surv_rate_x0 = exp(log_surv_rate_x0);
  vector[J] surv_rate_x1 = exp(log_surv_rate_x1);

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
