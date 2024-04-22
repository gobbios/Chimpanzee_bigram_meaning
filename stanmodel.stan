functions {
  // function to calculate distance between two call types
  // works on two probability vectors, i.e. the 'context distributions'
  // (don't forget: these prob vectors are not simplexes (and don't have to be))
  real dist_foo(int calltype1, int calltype2, int n_contexts, matrix muvals) {
    real res;
    vector[cols(muvals)] out1;
    vector[cols(muvals)] out2;
    vector[cols(muvals)] temp3;

    for (i in 1:cols(muvals)) {
      vector[n_contexts] temp1 = to_vector(muvals[calltype1, ]);
      vector[n_contexts] temp2 = to_vector(muvals[calltype2, ]);
      out1[i] = temp1[i];
      out2[i] = temp2[i];
    }

    temp3 = out1 - out2;
    for (i in 1:cols(muvals)) {
      temp3[i] = temp3[i]^2;
    }
    res = sqrt(sum(temp3));

    return(res);
  }
  real multinomial_logit2_lpmf(array[] int y, vector mu) {
     return multinomial_lpmf(y | softmax(mu));
  }

  // simulate recordings given a simplex with probabilities for each context (for one call type)
  // also simulate the event length
  // length_mu for a given call: xintercept + callblups[i]
  // returns adjusted probability vector
  vector simulate_recordings_rng(vector probvec, int n_recordings, real length_mu) {
    // return the actual probability vector (already multiplied)
    int n = num_elements(probvec);
    array[n_recordings, n] int event_mat = rep_array(0, n_recordings, n);
    vector[n] events = rep_vector(0.0, n);

    for (i in 1:n_recordings) {
      int still_looking = 0;
      while (still_looking == 0) {
        int event_length = 0;
        while (event_length < 1) {
          event_length = poisson_rng(exp(length_mu));
        }
        array[n] int temp = multinomial_rng(probvec, event_length);
        if (max(temp) == 1) {
          event_mat[i, ] = temp;
          still_looking = 1;
        }
      }
    }

    for (i in 1:n) { // go over columns (events)
      events[i] = (sum(event_mat[, i]) * 1.0)/n_recordings;
    }
    return(events);
  }

}

data {
  int<lower=1> N; // total number of observations
  int<lower=2> ncat; // number of contexts (response categories as factor)  [minus the first one, which is the reference]
  array[N, ncat] int Ymat; // response variable: integer levels of context
  array[N] int Ytrials; // response variable: integer levels of context
  int<lower=1> n_subjects; // number of subjects (callers)
  array[N] int<lower=1> index_subject; // indicator for subject (caller)
  int<lower=1> n_calltypes; // number of grouping levels
  array[N] int<lower=1> index_calltype; // indicator for calltype

  int n_distcomps;
  array[n_distcomps, 2] int index_comps;

  int n_combined;
  array[n_combined, 3] int index_combined;

  // model multiplication factor for event lengths
  array[N] int<lower=1> length; // length of event

}

parameters {
  // for multinomial model
  vector<lower=0>[ncat - 1] sd_vals_subject;
  matrix[n_subjects, ncat - 1] blups_mat_subject_z;

  vector<lower=0>[ncat - 1] sd_vals_calltype;
  matrix[n_calltypes, ncat - 1] blups_mat_calltype_z;

  // for length model
  vector[n_subjects] callerblups_z;
  vector[n_calltypes] callblups_z;
  real<lower=0> caller_sd;
  real<lower=0> call_sd;
  real xintercept;

}
transformed parameters {
  matrix[n_subjects, ncat - 1] blups_mat_subject_actual;
  matrix[n_calltypes, ncat - 1] blups_mat_calltype_actual;
  vector[n_subjects] callerblups;
  vector[n_calltypes] callblups;

  for (i in 1:(ncat - 1)) {
    blups_mat_subject_actual[, i] = sd_vals_subject[i] * blups_mat_subject_z[, i];
    blups_mat_calltype_actual[, i] = sd_vals_calltype[i] * blups_mat_calltype_z[, i];
  }

  callerblups = callerblups_z * caller_sd;
  callblups = callblups_z * call_sd;

}
model {
  // for length model
  vector[N] mu_length = rep_vector(0.0, N);
  matrix[N, ncat - 1] mu_mat = rep_matrix(0.0, N, ncat - 1);
  // linear predictor matrix
  matrix[N, ncat] mu;
  for (n in 1:N) {
    for (i in 1:(ncat - 1)) {
      mu_mat[n, i] = blups_mat_subject_actual[, i][index_subject[n]] + blups_mat_calltype_actual[, i][index_calltype[n]];
    }
  }

  for (n in 1:N) {
    mu[n, 1] = 0.0;
    for (i in 1:(ncat - 1)) {
      mu[n, i + 1] = mu_mat[n, i];
    }
  }

  for (n in 1:N) {
    target += multinomial_logit2_lpmf(Ymat[n] | to_vector(mu[n, ]));
  }

  // priors:
  for (i in 1:(ncat - 1)) {
    blups_mat_subject_z[, i] ~ normal(0, 1);
    blups_mat_calltype_z[, i] ~ normal(0, 1);
    sd_vals_subject[i] ~ student_t(3, 0, 2.5);
    sd_vals_calltype[i] ~ student_t(3, 0, 2.5);
  }

  mu_length = mu_length + xintercept;
  mu_length = mu_length + callerblups[index_subject];
  mu_length = mu_length + callblups[index_calltype];
  // can't be vectorized??? "Outcomes in truncated distributions must be univariate."
  for (i in 1:N) {
    length[i] ~ poisson(exp(mu_length[i])) T[1, ];
  }

  callerblups_z ~ normal(0, 1);
  callblups_z ~ normal(0, 1);
  caller_sd ~ exponential(1);
  call_sd ~ exponential(1);
  xintercept ~ student_t(3, 0, 1);
}
generated quantities {
  array[N, ncat] int y_rep; // probs for context for each observation
  vector[n_distcomps] distvals; // Eucl. distance for desired pairs
  matrix<lower=0>[n_calltypes, ncat] pred_mat_rep_multiplied; // predicted probs for each combination of calltype and context
  matrix[n_combined, ncat] pred_mat_rep_multiplied_combined; // combined predicted probs for each combination of calltype and context for all bigrams (max(call1, call2))
  array[N] int y_length_rep = rep_array(0, N); // event length

  // simulated event lengths/context lengths per observation:
  {
    vector[N] mu_length_rep = rep_vector(0.0, N);
    mu_length_rep = mu_length_rep + xintercept;
    mu_length_rep = mu_length_rep + callerblups[index_subject];
    mu_length_rep = mu_length_rep + callblups[index_calltype];
    // can't be vectorized!: "Outcomes in truncated distributions must be univariate."
    for (i in 1:N) {
      while (y_length_rep[i] < 1) {
        y_length_rep[i] = poisson_rng(exp(mu_length_rep[i]));
      }
    }
  }

  // probs for context for each calltype
  {
    pred_mat_rep_multiplied[, 1] = rep_vector(0.0, n_calltypes);
    for (i in 2:(ncat)) {
      for (j in 1:(ncat - 1)) {
        pred_mat_rep_multiplied[, i] = (blups_mat_calltype_actual[, i - 1][1:n_calltypes]);
      }
    }
    for (i in 1:n_calltypes) {
      // get valid probability vector for a given call type (simplex)
      pred_mat_rep_multiplied[i, ] = to_row_vector(softmax(to_vector(pred_mat_rep_multiplied[i, ])));
      // generate 1000 simulated data sets (i.e. 1000 simulated recordings) for the current call type
      // requires simulated event lengths for that specific call (and an 'average' individual) [this is done inside the function]
      // obtain adjusted probabilities ('context distribution')
      // this is not a simplex anymore, but still the max value per element is 1
      pred_mat_rep_multiplied[i, ] = to_row_vector(simulate_recordings_rng(to_vector(pred_mat_rep_multiplied[i, ]), 1000, xintercept + callblups[i]));
    }
  }

  // combined probs for context for each calltype
  {
    for (i in 1:n_combined) {
      for (k in 1:ncat) {
        pred_mat_rep_multiplied_combined[i, k] = max(pred_mat_rep_multiplied[index_combined[i, 2:3], k]);
      }
    }
  }

  // probs for each context for each observation
  {
    matrix[N, ncat - 1] mu_mat_rep = rep_matrix(0.0, N, ncat - 1);
    // linear predictor matrix
    matrix[N, ncat] mu_rep;
    for (n in 1:N) {
      for (i in 1:(ncat - 1)) {
        mu_mat_rep[n, i] = blups_mat_subject_actual[, i][index_subject[n]] + blups_mat_calltype_actual[, i][index_calltype[n]];
      }
    }

    for (n in 1:N) {
      mu_rep[n, 1] = 0.0;
      for (i in 1:(ncat - 1)) {
        mu_rep[n, i + 1] = mu_mat_rep[n, i];
      }
    }

    // sample from multinomial, but such that each event occurs max once, even though there is more than one trial
    for (n in 1:N) {
      int xconstraint = 0;
      while(xconstraint == 0) {
        array[ncat] int aux;
        y_rep[n, ] = multinomial_rng(softmax(to_vector(mu_rep[n, ])), Ytrials[n]);
        if (max(y_rep[n, ]) == 1) {
          xconstraint = 1;
        }
      }
    }
  }

  // Eucl. distance values for each calltype pair
  {
    for (k in 1:n_distcomps) {
      distvals[k] = dist_foo(index_comps[k, 1], index_comps[k, 2], ncat, pred_mat_rep_multiplied);
    }
  }
}
