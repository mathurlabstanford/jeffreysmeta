functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei){

    // variables that will be recalculated for each observation
	  real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
    real sigma;
		// will just be set to 1
		int n;

    real Si;


		// this will be the TOTALS for all observations
		matrix[2,2] fishinfototal;
		fishinfototal[1,1] = 0;
  	fishinfototal[1,2] = 0;
  	fishinfototal[2,1] = 0;
  	fishinfototal[2,2] = 0;


		// build a Fisher info matrix for EACH observation
		for (i in 1:k) {

		  Si = sqrt( tau^2 + sei[i]^2 );
      
      kmm = -Si^(-2);
      kms = 0;
      kss = -2 * tau^2 * Si^(-4); 

  		fishinfo[1,1] = -kmm;
      fishinfo[1,2] = -kms;
      fishinfo[2,1] = -kms;
      fishinfo[2,2] = -kss;

  		// add the new fisher info to the total one
  		fishinfototal = fishinfototal + fishinfo;
		}

		return sqrt(determinant(fishinfototal));
	}
}

data{
	int<lower=0> k;
  real sei[k];
	real y[k];
}

parameters{
  real mu;
	real<lower=0> tau;
}


model{
  // see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
	target += log( jeffreys_prior(mu, tau, k, sei) );
	for(i in 1:k) {
      y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
	}
}

// this chunk doesn't actually affect the model that's being fit to the data;
//  it's just re-calculating the prior, lkl, and post to return to user
// Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values.
//  The block is executed only after a sample has been generated'
generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei) );
  real log_post;

  // calculate joint log-likelihood
  for ( i in 1:k ){
    log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
  }

  log_post = log_prior + log_lik;
}
