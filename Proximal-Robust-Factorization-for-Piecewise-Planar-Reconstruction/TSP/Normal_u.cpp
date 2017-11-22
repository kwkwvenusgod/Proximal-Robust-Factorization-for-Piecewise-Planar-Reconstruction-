#include "Normal_u.h"

#ifndef pi
#define pi 3.14159265
#endif

// --------------------------------------------------------------------------
// -- Normal_u
// --   constructor; initializes to empty
// --------------------------------------------------------------------------
Normal_u::Normal_u() : IW(), theta(NULL), post_theta(NULL), temp_theta(NULL)
{
}

// --------------------------------------------------------------------------
// -- Normal_u
// --   constructor; intializes to the right dimensions and allocates memory
// --------------------------------------------------------------------------
Normal_u::Normal_u(int newD) : IW(newD)
{
   theta = allocate_memory<double>(D);
   post_theta = allocate_memory<double>(D);
   temp_theta = allocate_memory<double>(D);
}

// --------------------------------------------------------------------------
// -- Normal_u
// --   copy constructor;
// --------------------------------------------------------------------------
Normal_u::Normal_u(const Normal_u& that) : IW(that), kappa(that.kappa),
   post_kappa(that.post_kappa), temp_kappa(that.temp_kappa)
{
   theta = allocate_memory<double>(D);
   post_theta = allocate_memory<double>(D);
   temp_theta = allocate_memory<double>(D);

   for (int d=0; d<D; d++)
   {
      theta[d] = that.theta[d];
      post_theta[d] = that.post_theta[d];
      temp_theta[d] = that.temp_theta[d];
   }
}

// --------------------------------------------------------------------------
// -- Normal_u
// --   constructor; intializes to all the values given
// --------------------------------------------------------------------------
Normal_u::Normal_u(int newD, arr(double) newtheta, double newkappa, arr(double) newDelta, double newnu) :
   IW(newD, newDelta, newnu), kappa(newkappa)
{
   theta = allocate_memory<double>(D);
   post_theta = allocate_memory<double>(D);
   temp_theta = allocate_memory<double>(D);

   for (int d=0; d<D; d++)
      theta[d] = newtheta[d];
}
// --------------------------------------------------------------------------
// -- Normal_u
// --   constructor; intializes to all the values given
// --------------------------------------------------------------------------
Normal_u::Normal_u(int newD, arr(double) newtheta, arr(double) newoffset, double newkappa, arr(double) newDelta, double newnu) :
   IW(newD, newoffset, newDelta, newnu), kappa(newkappa)
{
   theta = allocate_memory<double>(D);
   post_theta = allocate_memory<double>(D);
   temp_theta = allocate_memory<double>(D);

   for (int d=0; d<D; d++)
      theta[d] = newtheta[d];
}

// --------------------------------------------------------------------------
// -- operator=
// --   assignment operator
// --------------------------------------------------------------------------
Normal_u& Normal_u::operator=(const Normal_u& that)
{
   if (this != &that)
   {
      cleanup();
      *((IW*)this) = that;
      kappa = that.kappa;

      theta = allocate_memory<double>(D);
      post_theta = allocate_memory<double>(D);
      temp_theta = allocate_memory<double>(D);
      post_kappa = that.post_kappa;
      temp_kappa = that.temp_kappa;
      for (int d=0; d<D; d++)
      {
         theta[d] = that.theta[d];
         post_theta[d] = that.post_theta[d];
         temp_theta[d] = that.temp_theta[d];
      }
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- ~Normal_u
// --   destructor
// --------------------------------------------------------------------------
Normal_u::~Normal_u()
{
   cleanup();
}

// --------------------------------------------------------------------------
// -- copy
// --   returns a copy of this
// --------------------------------------------------------------------------
Normal_u* Normal_u::copy()
{
   return new Normal_u(*this);
}


// --------------------------------------------------------------------------
// -- ~cleanup
// --   deletes all the memory allocated by this
// --------------------------------------------------------------------------
void Normal_u::cleanup()
{
   if (D>0)
   {
      if (theta!=NULL)
         deallocate_memory(theta);
      if (post_theta!=NULL)
         deallocate_memory(post_theta);
      if (temp_theta!=NULL)
         deallocate_memory(temp_theta);
   }
}


// --------------------------------------------------------------------------
// -- get_
// --   functions to get parameters.  returns pointers to actual data
// --------------------------------------------------------------------------
arr(double) Normal_u::get_theta()
{
   return theta;
}
double Normal_u::get_kappa()
{
   return kappa;
}

// --------------------------------------------------------------------------
// -- get_XXXX_mode
// --   returns a pointer to the most likely parameter. returns a pointer to
// -- internal memory that should *not* be deallocated
// --
// --   parameters
// --     - update : whether or not to update the posterior hyperparameters
//--    before returning the mode
// --------------------------------------------------------------------------
arr(double) Normal_u::get_mean_mode(bool update, bool useTemp)
{
   if (!useTemp && update && !updated)
      update_posteriors();
   if (useTemp)
      for (int d=0; d<D; d++)
         mean[d] = temp_theta[d];
   else
      for (int d=0; d<D; d++)
         mean[d] = post_theta[d];
   return mean;
}



// --------------------------------------------------------------------------
// -- update_posteriors
// --   updates the internal posterior hyperparameters based on the stored
// -- sufficient statistics
// --------------------------------------------------------------------------
void Normal_u::update_posteriors()
{
   if (!updated)
   {
      if (N==0)
      {
         post_kappa = kappa;
         post_nu = nu;
         for (int d=0; d<D; d++)
            post_theta[d] = theta[d];
         for (int d=0; d<D*D; d++)
            post_Delta[d] = Delta[d];
      }
      else
      {
         post_kappa = kappa + N;
         for (int d=0; d<D; d++)
         {
            theta[d] = total[d] / N;
            post_theta[d] = theta[d];
         }

         post_nu = nu + N;
         for (int d1=0; d1<D; d1++)
         {
            double temp_prior_theta = theta[d1];
            for (int d2=d1; d2<D; d2++)
            {
               double value = (nu*Delta[d1+d2*D] + total2[d1+d2*D] - N*temp_prior_theta*theta[d2]) / post_nu;
               post_Delta[d1+d2*D] = value;
               post_Delta[d2+d1*D] = value;
            }
         }
      }
      updated = true;
   }
}
// --------------------------------------------------------------------------
// -- update_posteriors_new
// --   updates the internal posterior hyperparameters based on the stored
// -- sufficient statistics and the new data points
// --
// --   parameters:
// --     - data: the new point to consider
// --------------------------------------------------------------------------
void Normal_u::update_posteriors_new(arr(double) data)
{
   temp_kappa = kappa + N + 1;
   for (int d=0; d<D; d++)
      temp_theta[d] = (total[d] + data[d]) / (N+1);

   temp_nu = nu + N + 1;
   for (int d1=0; d1<D; d1++)
   {
      double temp_prior_theta = temp_theta[d1];
      double temp_data = data[d1];
      for (int d2=d1; d2<D; d2++)
      {
         double value = (nu*Delta[d1+d2*D] + total2[d1+d2*D]+temp_data*data[d2] - (N+1)*temp_prior_theta*temp_theta[d2] ) / temp_nu;
         temp_Delta[d1+d2*D] = value;
         temp_Delta[d2+d1*D] = value;
      }
   }
}
void Normal_u::update_posteriors_rem(arr(double) data)
{
   if (N==1)
   {
      temp_kappa = kappa;
      temp_nu = nu;
      for (int d=0; d<D; d++)
         temp_theta[d] = theta[d];
      for (int d=0; d<D*D; d++)
         temp_Delta[d] = Delta[d];
   }
   else
   {
      temp_kappa = kappa + N - 1;
      for (int d=0; d<D; d++)
         temp_theta[d] = (total[d] - data[d]) / (N-1);

      temp_nu = nu + N - 1;
      for (int d1=0; d1<D; d1++)
      {
         double temp_prior_theta = temp_theta[d1];
         double temp_data = data[d1];
         for (int d2=d1; d2<D; d2++)
         {
            double value = (nu*Delta[d1+d2*D] + total2[d1+d2*D]-temp_data*data[d2] - (N-1)*temp_prior_theta*temp_theta[d2]) / temp_nu;
            temp_Delta[d1+d2*D] = value;
            temp_Delta[d2+d1*D] = value;
         }
      }
   }
}

// --------------------------------------------------------------------------
// -- update_posteriors_new
// --   updates the internal posterior hyperparameters based on the stored
// -- sufficient statistics and the new data points
// --
// --   parameters:
// --     - other : another Normal_u to merge with and test the likelihood
// --------------------------------------------------------------------------
void Normal_u::update_posteriors_new(IW* &other)
{
   arr(double) o_total = other->get_total();
   arr(double) o_total2 = other->get_total2();
   double o_N = other->get_N();

   temp_kappa = post_kappa + o_N;
   for (int d=0; d<D; d++)
      temp_theta[d] = (total[d] + o_total[d]) / (N + o_N);

   temp_nu = post_nu + o_N;
   for (int d1=0; d1<D; d1++)
   {
      double temp_prior_theta = temp_theta[d1];
      for (int d2=d1; d2<D; d2++)
      {
         double value = (nu*Delta[d1+d2*D] + total2[d1+d2*D]+o_total2[d1+d2*D] - (N+o_N)*temp_prior_theta*temp_theta[d2]) / temp_nu;
         temp_Delta[d1+d2*D] = value;
         temp_Delta[d2+d1*D] = value;
      }
   }
}
void Normal_u::update_posteriors_new(IW* &other1, IW* &other2)
{
   update_posteriors();

   arr(double) o1_total = other1->get_total();
   arr(double) o1_total2 = other1->get_total2();
   arr(double) o2_total = other2->get_total();
   arr(double) o2_total2 = other2->get_total2();
   double o_N = other1->get_N() + other2->get_N();

   temp_kappa = post_kappa + o_N;
   for (int d=0; d<D; d++)
      temp_theta[d] = (total[d] + o1_total[d] + o2_total[d]) / (N+o_N);

   temp_nu = post_nu + o_N;
   for (int d1=0; d1<D; d1++)
   {
      double temp_prior_theta = temp_theta[d1];
      for (int d2=d1; d2<D; d2++)
      {
         double value = (nu*Delta[d1+d2*D] + total2[d1+d2*D]+o1_total2[d1+d2*D]+o2_total2[d1+d2*D] - (N+o_N)*temp_prior_theta*temp_theta[d2]) / temp_nu;
         temp_Delta[d1+d2*D] = value;
         temp_Delta[d2+d1*D] = value;
      }
   }
}






// --------------------------------------------------------------------------
// -- calc_log_data_test_point
// --   calculate log p(data | x,hyperparams) for the prior hyper parameters.
// --------------------------------------------------------------------------
double Normal_u::calc_log_data_test_point(arr(double) data)
{
   update_posteriors();

   for (int d=0; d<D; d++)
      mean[d] = (total[d] ) / (N) - data[d];
//      mean[d] = (total[d]+data[d] ) / (N+1) - data[d];  not sure if we should include the point or not
   for (int d1=0; d1<D; d1++)
   {
      for (int d2=d1; d2<D; d2++)
      {
         double value = post_nu*post_Delta[d1+d2*D]*(post_kappa+1) / (post_kappa*(post_nu-D+1));
         cov[d1+d2*D] = value;
         cov[d2+d1*D] = value;
      }
   }

   double this_nu = post_nu-D+1;

   return my_sf_lngamma(this_nu+D) - my_sf_lngamma(this_nu) - 0.5*D*log(this_nu*pi) - 0.5*log(det(cov,D)) - 0.5*(this_nu+D)*log(1.0 + 1.0/this_nu * xt_Ainv_x(mean, cov, D));
}
// --------------------------------------------------------------------------
// -- calc_log_data_test_point_rem
// --   calculate log p(data | x,hyperparams) for the prior hyper parameters.
// --------------------------------------------------------------------------
double Normal_u::calc_log_data_test_point_rem(arr(double) data)
{
   update_posteriors_rem(data);

   for (int d=0; d<D; d++)
      mean[d] = (total[d]-data[d]) / (N-1) - data[d];
//      mean[d] = total[d] / N - data[d];  not sure if we should include the point or not
   for (int d1=0; d1<D; d1++)
   {
      for (int d2=d1; d2<D; d2++)
      {
         double value = temp_nu*temp_Delta[d1+d2*D]*(temp_kappa+1) / (temp_kappa*(temp_nu-D+1));
         cov[d1+d2*D] = value;
         cov[d2+d1*D] = value;
      }
   }

   double this_nu = temp_nu-D+1;

   return my_sf_lngamma(this_nu+D) - my_sf_lngamma(this_nu) - 0.5*D*log(this_nu*pi) - 0.5*log(det(cov,D)) - 0.5*(this_nu+D)*log(1.0 + 1.0/this_nu * xt_Ainv_x(mean, cov, D));
}


double Normal_u::calc_logdata_internal(bool useTemp)
{
   double this_nu = (useTemp ? temp_nu : post_nu);
   arr(double) this_Delta = (useTemp ? temp_Delta : post_Delta);
   double this_kappa = (useTemp ? temp_kappa : post_kappa);
   double N = this_nu - nu;

   if (N==0) return 0;

   return -0.5*N*D*1.144729885849 + mylogmgamma(this_nu,D) - mylogmgamma(nu,D) + 0.5*nu*(D*log(nu) + log(det(Delta,D))) - 0.5*this_nu*(D*log(this_nu) + log(det(this_Delta,D))) + 0.5*D*log(kappa/this_kappa);
}

// --------------------------------------------------------------------------
// -- calc_logposterior_internal
// --   calculate log p(params | post_hyperparams) for the posterior hyper
// -- parameters that are stored.
// --------------------------------------------------------------------------
double Normal_u::calc_logposterior_internal(bool useTemp)
{
   // find the maximum likely parameters
   get_mean_mode(false, useTemp);
   get_cov_mode(false, useTemp);

   double this_nu = (useTemp ? temp_nu : post_nu);
   arr(double) this_Delta = (useTemp ? temp_Delta : post_Delta);
   double this_kappa = (useTemp ? temp_kappa : post_kappa);
   arr(double) this_theta = (useTemp ? temp_theta : post_theta);

   // log(|cov|)
   double logdeterminant = log(det(cov, D));

   // covariance
   //double logprob = 0.5*post_nu*log(post_nu*det(post_Delta,D)) - 0.5*post_nu*D*log(2) - logmgamma(post_nu/2, D) - 0.5*(post_nu+D+1)*logdeterminant - 0.5*D*(D+post_nu+1);
   double logprob = 0.5*this_nu*(D*log(this_nu) + log(det(this_Delta,D))) - 0.3465735903*this_nu*D - mylogmgamma(this_nu, D) - 0.5*(this_nu+D+1)*logdeterminant - 0.5*D*(D+this_nu+1);
   //logprob += -0.5*D*log(2*pi) -0.5*(logdeterminant - log(post_kappa));

   /*for (int d1=0; d1<D; d1++)
   {
      mexPrintf(" -> ");
      for (int d2=0; d2<D; d2++)
      {
         mexPrintf("%f\t", this_Delta[d1+d2*D]);
      }
      mexPrintf("\n");
   }
   mexPrintf("-----------\n");
   for (int d1=0; d1<D; d1++)
   {
      mexPrintf(" -> ");
      for (int d2=0; d2<D; d2++)
      {
         mexPrintf("%f\t", cov[d1+d2*D]);
      }
      mexPrintf("\n");
   }
   mexPrintf(" -------> nu=%f,%f,%f,%f,%f,%f \t prob=%f\n", this_nu, 0.5*this_nu*log(this_nu*det(this_Delta,D)), - 0.3465735903*this_nu*D, - mylogmgamma((int)(this_nu+0.5), D), - 0.5*(this_nu+D+1)*logdeterminant, - 0.5*D*(D+this_nu+1), logprob);
   //mexPrintf(" -> nu=%f\tlogmgamma=%f\n", this_nu, mylogmgamma((int)(this_nu+0.5), D));*/

   logprob += -0.9189385332*D -0.5*(logdeterminant - D*log(this_kappa));

   if (logprob!=logprob)
   {
      mexErrMsgTxt("AHHHHHHHHHHHH logprob is nan!\n");
   }

   //logprob -= (D-1)*this_nu/2*log(this_nu) + (D-1)*log(this_kappa)/2;

   return logprob;
   /*
   double logprob = IW::calc_logposterior_internal();
   logprob += -0.5*D*log(2*pi) -0.5*(logdeterminant - log(post_kappa));
   return logprob;*/
}

