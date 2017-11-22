#ifndef _SP_H_INCLUDED_
#define _SP_H_INCLUDED_
#include "linkedList.cpp"
#include "array.h"
#include "matrix.h"
#include "mex.h"
#include <math.h>
#include "NormalD.h"


#ifndef pi
#define pi 3.14159265
#endif

#include "helperMEX.h"
#include "gsl/gsl_sf_gamma.h"

double calc_log_label(double N, double alpha, double epsilon);
std::pair<int, int> increment_neighbor_count(std::pair<int, int> n);
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n);
std::pair<int, int> increment_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2);
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2);

class SP
{
private:
   // Likelihood term
   NormalD* pos;
   NormalD* app;

   // probabilities
   double log_likelihood;
   double log_likelihood_empty;

   // SP things
   int N;
   linkedList<int> pixels;
   linkedList<int> borders;
   linkedList< std::pair<int, int> > neighbors;

   arr(double) prev_v;

   bool is_old;

   unsigned long UID;
   friend class IMG;

public:
   // --------------------------------------------------------------------------
   // -- SP
   // --   constructor; initializes to empty... probably shouldn't use
   // --------------------------------------------------------------------------
   SP();

   // --------------------------------------------------------------------------
   // -- initialize
   // --   initializes all the variables conditioned on pos and app being set
   // --------------------------------------------------------------------------
   void initialize();
   // --------------------------------------------------------------------------
   // -- SP
   // --   constructor; initializes to empty new super pixel
   // --------------------------------------------------------------------------
   SP(NormalD &new_pos, NormalD &new_app, unsigned long new_UID, bool isOld=false, arr(double) the_prev_v=NULL);
   SP(const SP& that);

   // --------------------------------------------------------------------------
   // -- operator=
   // --   assignment operator
   // --------------------------------------------------------------------------
   SP& operator=(const SP& that);
   // --------------------------------------------------------------------------
   // -- SP
   // --   destructor;
   // --------------------------------------------------------------------------
   ~SP();

   // --------------------------------------------------------------------------
   // -- empty
   // --   Empties out the SP
   // --
   // --   parameters:
   // --     - check_merged : checks to see if the linked list, pixels, is
   // --       empty. if it isn't, throws exception.
   // --------------------------------------------------------------------------
   void empty(bool check_merged=true);

   int get_N();
   double get_log_likelihood_empty();
   double get_log_likelihood();
   unsigned long get_UID();
   arr(double) get_prev_v();

   /*arr(double) get_flow();*/
   void set_flow(arr(double) flow);
   void set_flow(double flowx, double flowy);

   // --------------------------------------------------------------------------
   // -- isempty
   // --   returns whether or not the super pixel contains any pixels
   // --------------------------------------------------------------------------
   bool isempty();
   bool isold();
   void checkCount();

   // --------------------------------------------------------------------------
   // -- calculate_log_probs
   // --   Updates log_likelihood and log_label based on the current parameters
   // --------------------------------------------------------------------------
   void calculate_log_probs();

   // --------------------------------------------------------------------------
   // -- add_pixel_init
   // --   Adds a pixel to the super pixel, updates the appearance and position
   // -- parameters, and the linked lists. Does not update the likelihoods.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - is_border : indicator as to whether or not to add to border LL
   // --   return parameters:
   // --     - pixel_ptr : a pointer to the added linked list node pixel
   // --     - border_ptr : a pointer to the added linked list node border
   // --------------------------------------------------------------------------
   void add_pixel_init(arr(double) data, int index, bool is_border,
      linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp=true);
   // --------------------------------------------------------------------------
   // -- add_pixel
   // --   Adds a pixel to the super pixel, updates the appearance and position
   // -- parameters, the linked lists, and the probabilities.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - is_border : indicator as to whether or not to add to border LL
   // --   return parameters:
   // --     - pixel_ptr : a pointer to the added linked list node pixel
   // --     - border_ptr : a pointer to the added linked list node border
   // --------------------------------------------------------------------------
   void add_pixel(arr(double) data, int index, bool is_border,
      linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp=true);

   void merge_with(SP *other, arr(int) label, arr(linkedListNode<int>*) border_ptr, int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- rem_pixel
   // --   Removes a pixel from the super pixel, updates the appearance and
   // -- position parameters, and the linked lists.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - pixel_ptr : a poitner into the pixels linked list to be removed
   // --     - border_ptr : a pointer into the borders linked list to be removed
   // --------------------------------------------------------------------------
   void rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp=true);
   void rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr, bool doApp=true);


   // --------------------------------------------------------------------------
   // -- fix_borders
   // --   Fixes the border linked list and the border_ptr image for a single
   // -- super pixel.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - border_ptr : the border_ptr image
   // --     - (xdim, ydim) : the size of the image
   // --------------------------------------------------------------------------
   void fix_borders(arr(int) label, arr(linkedListNode<int>*) border_ptr, int xdim, int ydim);



   // --------------------------------------------------------------------------
   // -- update_neighbors_label_rem
   // --   Decrements the corresponding neighbor count for neighbor_label
   // --
   // --   parameters:
   // --     - neighbor_label : the label of the neighbor to decrement
   // --------------------------------------------------------------------------
   void update_neighbors_label_rem(int neighbor_label);

   // --------------------------------------------------------------------------
   // -- update_neighbors_label_rem_check
   // --   Checks to see if the neighbor count at index should be decremented
   // -- by removing one neighbor of label neighbor_label. If so, it decrements.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index bordering the removed pixel
   // --     - (xdim,ydim) : dimensions of image
   // --     - neighbor_label : the label of the neighbor to decrement
   // --------------------------------------------------------------------------
   void update_neighbors_label_rem_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label);

   // --------------------------------------------------------------------------
   // -- update_neighbors_add_self
   // --   Updates the neighbor lists and counts by adding one particular pixel
   // -- at index. Does not update the neighboring neighbor lists.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index of the added pixel
   // --     - (xdim,ydim) : dimensions of image
   // --------------------------------------------------------------------------
   void update_neighbors_add_self(arr (int) label, int index, int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- update_neighbors_self
   // --   Updates the neighbor lists and counts by looking at all borders.
   // -- Empties previous list. The borders list must be correct!
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - (xdim,ydim) : dimensions of image
   // --------------------------------------------------------------------------
   void update_neighbors_self(arr (int) label, int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- update_neighbors_label_add
   // --   Increments the corresponding neighbor count for neighbor_label
   // --
   // --   parameters:
   // --     - neighbor_label : the label of the neighbor to increment
   // --------------------------------------------------------------------------
   void update_neighbors_label_add(int neighbor_label);
   // --------------------------------------------------------------------------
   // -- update_neighbors_label_add_check
   // --   Checks to see if the neighbor count at index should be incremented
   // -- by adding one neighbor of label neighbor_label. If so, it increments.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index bordering the added pixel
   // --     - (xdim,ydim) : dimensions of image
   // --     - neighbor_label : the label of the neighbor to increment
   // --------------------------------------------------------------------------
   void update_neighbors_label_add_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label);

   bool has_no_neighbors();



   void switch_priors(SP* other);


   double log_likelihood_test_point_MM(arr(double) data);
   double log_likelihood_test_point_MM_pos(arr(double) data);
   double log_likelihood_test_point_app(arr(double) data);
   double log_likelihood_test_point_app_rem(arr(double) data);
   double log_likelihood_test_point_pos(arr(double) data);
   double log_likelihood_test_point_pos_rem(arr(double) data);

   double get_log_likelihood_app();
   double get_log_likelihood_pos();


   double log_likelihood_switch_app_prior(SP* new_prior);
   double log_likelihood_switch_pos_prior(SP* new_prior);
   double log_likelihood_switch_prior(SP* new_prior);

   // --------------------------------------------------------------------------
   // -- log_likelihood_test_point
   // --   finds the change in probability for adding this data point to the super
   // -- pixel
   // --
   // --   parameters
   // --     - data : a [1 5] vector of a test point to add
   // --------------------------------------------------------------------------
   double log_likelihood_test_point(arr(double) data, bool checkApp=true);
   double log_likelihood_test_point_rem(arr(double) data, bool checkApp=true);
   // --------------------------------------------------------------------------
   // -- log_likelihood_test_merge
   // --   finds the change in probability for adding the points in indices to
   // -- the super pixel
   // --
   // --   parameters
   // --     - other : another SP that we are testing for a merge
   // --------------------------------------------------------------------------
   double log_likelihood_test_merge(SP *other);
   double log_likelihood_test_merge_pos(SP *other);
   double log_likelihood_test_merge(SP *other1, SP *other2);

   // --------------------------------------------------------------------------
   // -- get_XXXX
   // --   returns a pointer to the most likely parameter. returns a pointer to
   // -- internal memory that should *not* be deallocated
   // --
   // --   parameters
   // --     - update : whether or not to update the posterior hyperparameters
   //--    before returning the mode
   // --------------------------------------------------------------------------
   arr(double) get_flow();
   arr(double) get_theta_pos();
   arr(double) get_theta_app();
   arr(double) get_Delta_pos();
   arr(double) get_Delta_app();
   arr(double) get_Sigma_pos();
   arr(double) get_Sigma_app();
   // iidmode == true indicates to find the mode of the mean based on current
   // pixels (ignoring possible neighbors)
   arr(double) get_mean_pos(bool iidmode=true);
   arr(double) get_mean_app();

   arr(double) get_total_app();
   arr(double) get_total_pos();
   arr(double) get_total2_pos();
   double get_sumlogDelta_div2_pos();

   void set_mean_app(arr(double) new_mean);
   void set_mean_pos(arr(double) new_mean);
   void set_meansum_pos(arr(double) new_mean1, arr(double) new_mean2);
   void set_flowsum(arr(double) flow, arr(double) other);
};

#endif
