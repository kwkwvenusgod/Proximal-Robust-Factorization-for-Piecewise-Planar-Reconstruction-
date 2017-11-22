#ifndef _Normal_u_H_INCLUDED_
#define _Normal_u_H_INCLUDED_

#include "matrix.h"
#include "mex.h"
#include <math.h>
#include "array.h"
#include "IW.h"

#include "helperMEX.h"
#include "debugMEX.h"

#include "linear_algebra.h"

class Normal_u : public IW
{
private:
   // prior hyperparameters for normal mean
   arr(double) theta;
   double kappa;

   // posterior hyperparameters
   arr(double) post_theta;
   double post_kappa;

   // temporary hyperparameters
   arr(double) temp_theta;
   double temp_kappa;

public:
   // --------------------------------------------------------------------------
   // -- Normal_u
   // --   constructor; initializes to empty
   // --------------------------------------------------------------------------
   Normal_u();
   // --------------------------------------------------------------------------
   // -- Normal_u
   // --   constructor; intializes to the right dimensions and allocates memory
   // --------------------------------------------------------------------------
   Normal_u(int newD);
   // --------------------------------------------------------------------------
   // -- Normal_u
   // --   copy constructor;
   // --------------------------------------------------------------------------
   Normal_u(const Normal_u& that);
   // --------------------------------------------------------------------------
   // -- Normal_u
   // --   constructor; intializes to all the values given
   // --------------------------------------------------------------------------
   Normal_u(int newD, arr(double) newtheta, double newkappa, arr(double) newDelta, double newnu);
   // --------------------------------------------------------------------------
   // -- Normal_u
   // --   constructor; intializes to all the values given
   // --------------------------------------------------------------------------
   Normal_u(int newD, arr(double) newtheta, arr(double) newoffset, double newkappa, arr(double) newDelta, double newnu);
   // --------------------------------------------------------------------------
   // -- operator=
   // --   assignment operator
   // --------------------------------------------------------------------------
   Normal_u& operator=(const Normal_u& that);
   // --------------------------------------------------------------------------
   // -- ~Normal_u
   // --   destructor
   // --------------------------------------------------------------------------
   ~Normal_u();
   // --------------------------------------------------------------------------
   // -- copy
   // --   returns a copy of this
   // --------------------------------------------------------------------------
   Normal_u* copy();

   // --------------------------------------------------------------------------
   // -- get_
   // --   functions to get parameters.  returns pointers to actual data
   // --------------------------------------------------------------------------
   arr(double) get_theta();
   double get_kappa();

   // --------------------------------------------------------------------------
   // -- get_XXXX_mode
   // --   returns a pointer to the most likely parameter. returns a pointer to
   // -- internal memory that should *not* be deallocated
   // --
   // --   parameters
   // --     - update : whether or not to update the posterior hyperparameters
   //--    before returning the mode
   // --------------------------------------------------------------------------
   arr(double) get_mean_mode(bool update, bool useTemp=false);

public:
   // --------------------------------------------------------------------------
   // -- ~cleanup
   // --   deletes all the memory allocated by this
   // --------------------------------------------------------------------------
   void cleanup();

   // --------------------------------------------------------------------------
   // -- update_posteriors
   // --   updates the internal posterior hyperparameters based on the stored
   // -- sufficiend statistics
   // --------------------------------------------------------------------------
   void update_posteriors();
   // --------------------------------------------------------------------------
   // -- update_posteriors_new
   // --   updates the internal posterior hyperparameters based on the stored
   // -- sufficient statistics and the new data points
   // --
   // --   parameters:
   // --     - data: the new point to consider
   // --------------------------------------------------------------------------
   void update_posteriors_new(arr(double) data);
   void update_posteriors_rem(arr(double) data);
   // --------------------------------------------------------------------------
   // -- update_posteriors_new
   // --   updates the internal posterior hyperparameters based on the stored
   // -- sufficient statistics and the new data points
   // --
   // --   parameters:
   // --     - other : another Normal_u to merge with and test the likelihood
   // --------------------------------------------------------------------------
   void update_posteriors_new(IW* &other);
   void update_posteriors_new(IW* &other1, IW* &other2);

   // --------------------------------------------------------------------------
   // -- calc_log_data_test_point
   // --   calculate log p(data | x,hyperparams) for the prior hyper parameters.
   // --------------------------------------------------------------------------
   double calc_log_data_test_point(arr(double) data);
   // --------------------------------------------------------------------------
   // -- calc_log_data_test_point_rem
   // --   calculate log p(data | x,hyperparams) for the prior hyper parameters.
   // --------------------------------------------------------------------------
   double calc_log_data_test_point_rem(arr(double) data);
   double calc_logdata_internal(bool useTemp);
   // --------------------------------------------------------------------------
   // -- calc_logposterior_internal
   // --   calculate log p(params | post_hyperparams) for the posterior hyper
   // -- parameters that are stored.
   // --------------------------------------------------------------------------
   double calc_logposterior_internal(bool useTemp=false);
};


#endif
