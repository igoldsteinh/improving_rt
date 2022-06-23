//current state of the art model as of 1/13/2022
//gamma model of incidence
//exponential seeded incidence
//observations as nb with mean as function of tests
functions{
 vector my_reverse(vector v) {
   int num = rows(v);
   vector[num] v2;
   for (i in 1:num) {
     v2[i] = v[num - i + 1];
   }
   
   return(v2);
 }
}
data {
  //data length (number of time steps)
  int<lower=0> n;
  //number of discretized gamma values (previous incidence days to use)
  int<lower=0> d; 
  //observed cases at time t
  int<lower=0> obs[n];
  //tests performed at time t
  int<lower=0> test[n];
  //discretized gamma values
  vector[d] w; 
  //incidences 3 through n
  //int<lower=0> num_weights[n-2];
  vector[d+1] delay_weights;
  int<lower=0> prev_vals;

  real log_incid_rate_mean;
  real log_incid_rate_sd;
    
  real log_sigma_mu;
  real log_sigma_sd;
  
  real log_rho_mu;
  real log_rho_sd;
  
  real log_r0_mu;
  real log_r0_sd;
  
  real kappa_mu;
  real kappa_sd;
  
  // int obs_pred;
  // int test_pred[obs_pred];
}

// The parameters accepted by the model.
parameters {
  real log_rt0_raw; //initial log r0
  vector[n-1] diff_rt; //vector of changes in log_rt for each time point
  real log_incid_rate_raw; // standard deviation for incidence
  real<lower=0> seed_incid_raw[prev_vals];
  real log_sigma; //variance on random walk parameter
  real <lower=0> log_rho; // raw underreporting parameter for cases as a function of tests
  real<lower=0> i_raw[n]; //incidence transformed so that the mean is always 1
  real<lower = 0> kappa;  // overdispersion parameter for cases
  real<lower=0> exp_rate;


  
}

transformed parameters{
  //kappa for beta binomial model
  //create a vector of log_rts
  vector[n] log_rt;
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd*log_sigma);
  real<lower=0> incid_rate = exp(log_incid_rate_raw*log_incid_rate_sd + 
                                 log_incid_rate_mean);
  vector<lower=0>[prev_vals] seed_incid;
  real<lower=0> incid[n]; //true unobserved incidence, transformed from raw
  real<lower=0> weighted_sum[n];
  real <lower=0> delay_sum[n];
  vector<lower=0>[prev_vals+n] aug_incid;
  real <lower=0> rho = exp(log_rho*log_rho_sd + log_rho_mu);


  //seed incid
  for (i in 1:prev_vals){
    seed_incid[i] = seed_incid_raw[i]*inv(exp_rate);
  }

  //actual standard deviation for incidence not equal to 1
  //real<lower=0>sigma;
  
  //first entry should just be initial log_rt
  log_rt[1] = log_r0_mu+log_r0_sd*log_rt0_raw;
  

  //next entry should be previous log_rt plus the scaled difference
  for (r in 2:n){
  
  log_rt[r] = log_rt[r-1] + diff_rt[r-1]*inv_sqrt(n) * sigma;
  } 
  

  
  //putting incidence fully into latent land
  //starting with aug_incid first 8 values being the seed_incids
  for (a in 1:prev_vals) {
      aug_incid[a] = seed_incid[a];
    }
  
  //true incidence modeled in terms of previous true incidence and seed incidence
  { int c;
    c = prev_vals;
    for (t in 1:n) {
      weighted_sum[t] = 0;
    
      if (c <= d ) {
        
        weighted_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(w, c)));

        }
    
      if (c > d) {

        weighted_sum[t] = dot_product(segment(aug_incid, c-d +1, d), 
                                      my_reverse(head(w, d)));
          
    
        }
    //create incid and aug incid
      incid[t] = exp(log_rt[t])*weighted_sum[t]*i_raw[t];
      aug_incid[t+prev_vals] = incid[t];
      c = c+1;
      }
    }
    
      { int c;
    //changing the indexing on this so that it allows for weight on the same day as 
    //the observation
    c = prev_vals +1;
    for (t in 1:n) {
      delay_sum[t] = 0;
    
      if (c <= d ) {
        
        delay_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(delay_weights, c)));

        }
    
      if (c > d) {

        delay_sum[t] = dot_product(segment(aug_incid, c-d +1, d), 
                                      my_reverse(head(delay_weights, d)));
          
    
        }
      c = c+1;
      }
    }
}
// The model to be estimated
model {


  //prior on underreporting
  target += normal_lpdf(log_rho| 0, 1);

  //truncated normal prior on kappa, truncated at zero
  target +=normal_lpdf(kappa|kappa_mu, kappa_sd) -
           log_diff_exp(0, normal_lcdf(0|kappa_mu, kappa_sd));

  //c should be 1 for time step 2, in loop starts at 0 then is updated
  // prior on random walk sigma
  target += normal_lpdf(log_sigma | 0, 1);

  //priors for log R0
  target += normal_lpdf(log_rt0_raw| 0, 1);
  //prior on sd of incidence
  target+= normal_lpdf(log_incid_rate_raw| 0, 1);
  
    //seed incidence
  target+=exponential_lpdf(exp_rate| 0.3);
  
  for (i in 1:prev_vals) {
    target+=exponential_lpdf(seed_incid_raw[i]| 1);

  }

  //unscaled difference in rts
  for (r in 1:n-1) {
    
    target+=normal_lpdf(diff_rt[r]| 0,1);

  }
  
  //seed incidence

  //raw incidence modeled in terms of previous true incidence and seed incidence
  for (t in 1:n) {
    
    target+=gamma_lpdf(i_raw[t]|exp(log_rt[t])*weighted_sum[t]*incid_rate, 
    exp(log_rt[t])*weighted_sum[t]*incid_rate);

  }


  //observed incidence model
      for (t in 1:n) {
        target+= neg_binomial_2_lpmf(obs[t] |rho*test[t]*delay_sum[t], kappa );

  }

  


}

generated quantities{
  // vector[n + obs_pred] log_gen_rt;
  int gen_obs[n];
  vector[n] log_lik;
  // 
  // //known values
  // for (a in 1:prev_vals){
  //   gen_incid[a] = seed_incid[a];
  // }
  // 
  // for (t in prev_vals+1:n+prev_vals) {
  //   gen_incid[t] = incid[t-prev_vals];
  // }
  // 
  for (t in 1:n) {
      gen_obs[t] = neg_binomial_2_rng(rho*test[t]*delay_sum[t], kappa);
  }
  for (t in 1:n) {
    log_lik[t] = neg_binomial_2_lpmf(obs[t] |rho*test[t]*delay_sum[t], kappa);
  }
  // 
  // for (t in 1:n) {
  //   log_gen_rt[t] = log_rt[t];
  // }
  // // 
  // // exp_log_gen_rt[1] = log_rt[1];
  // // 
  // // for (t in 1:n-1) {
  // //   gen_diff_rt[t] = normal_rng(0, 1);
  // // }
  // // for (r in 2:n){
  // // 
  // // exp_log_gen_rt[r] = exp_log_gen_rt[r-1] + gen_diff_rt[r-1]*inv_sqrt(n) * sigma;
  // // } 
  // 
  // //prediction time 
  // 
  // for (t in n+1:n+obs_pred) {
  //   log_gen_rt[t] = normal_rng(log_gen_rt[t-1], inv_sqrt(n) *sigma);
  // }
  // 
  // c = n+ prev_vals;
  // for (t in n+prev_vals+1:n + prev_vals + obs_pred) {
  //   
  //           gen_weightedsum= dot_product(segment(gen_incid, c-d +1, d), 
  //                                     reverse(head(w, d)));
  // 
  //  //print(gen_weightedsum);
  //   // print(log_gen_rt[t-prev_vals]);
  //   // print("t-prevvals");
  //   // print(t-prev_vals);
  //   // 
  //   // print("t"); 
  //   // print(t);
  //   gen_incid[t] = gamma_rng(gen_weightedsum*exp(log_gen_rt[t-prev_vals])*incid_rate,
  //                 incid_rate);
  //                 
  //   c = c + 1;
  // 
  // }
  // 
  // for (t in n+1:n + obs_pred) {
  //   gen_obs[t] = neg_binomial_2_rng(rho*test_pred[t-n]*gen_incid[t], kappa);
  // }
  // 
  // 

}


