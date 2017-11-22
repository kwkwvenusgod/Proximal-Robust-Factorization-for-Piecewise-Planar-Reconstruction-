#include "SP.h"


std::pair<int, int> increment_neighbor_count(std::pair<int, int> n)
{
   n.second++;
   return n;
}
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n)
{
   n.second--;
   return n;
}

std::pair<int, int> increment_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2)
{
   n.second+=n2.second;
   return n;
}
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2)
{
   n.second-=n2.second;
   return n;
}


// --------------------------------------------------------------------------
// -- SP
// --   constructor; initializes to empty
// --------------------------------------------------------------------------
SP::SP() : N(0), pos(NULL), app(NULL), is_old(false)
{
   prev_v = allocate_memory<double>(2,0);
}

// --------------------------------------------------------------------------
// -- initialize
// --   initializes all the variables conditioned on pos and app being set
// --------------------------------------------------------------------------
void SP::initialize()
{
   log_likelihood_empty = pos->calc_logposterior() + app->calc_logposterior();
   log_likelihood = log_likelihood_empty;
}

// --------------------------------------------------------------------------
// -- SP
// --   constructor; initializes to empty new super pixel
// --------------------------------------------------------------------------
SP::SP(NormalD &new_pos, NormalD &new_app, unsigned long new_UID, bool isOld, arr(double) the_prev_v) :
   N(0), UID(new_UID), is_old(isOld)
{
   pos = new NormalD(new_pos);
   app = new NormalD(new_app);

   if (the_prev_v == NULL)
      prev_v = allocate_memory<double>(2,0);
   else
   {
      prev_v = allocate_memory<double>(2);
      prev_v[0] = the_prev_v[0];
      prev_v[1] = the_prev_v[1];
   }

   initialize();
}


SP::SP(const SP& that) :
   pixels(that.pixels), borders(that.borders), neighbors(that.neighbors),
   log_likelihood(that.log_likelihood), log_likelihood_empty(that.log_likelihood_empty),
   UID(that.UID), is_old(that.is_old), N(that.N)
{
   pos = that.pos->copy();
   app = that.app->copy();
   prev_v[0] = that.prev_v[0];
   prev_v[1] = that.prev_v[1];
}


// --------------------------------------------------------------------------
// -- operator=
// --   assignment operator
// --------------------------------------------------------------------------
SP& SP::operator=(const SP& that)
{
   if (this != &that)
   {
      if (pos!=NULL)
         delete pos;
      if (app!=NULL)
         delete app;

      pos = that.pos->copy();
      app = that.app->copy();

      log_likelihood = that.log_likelihood;
      log_likelihood_empty = that.log_likelihood_empty;

      prev_v[0] = that.prev_v[0];
      prev_v[1] = that.prev_v[1];

      UID = that.UID;
      is_old = that.is_old;

      pixels.clear();
      borders.clear();
      neighbors.clear();
      mexPrintf("SP ASSIGNMENT OPERATOR NOT COPYING LINKED LISTS AND UID COPIED!\n");
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- SP
// --   destructor;
// --------------------------------------------------------------------------
SP::~SP()
{
   deallocate_memory(prev_v);
   if (pos!=NULL)
      delete pos;
   if (app!=NULL)
      delete app;
}


// --------------------------------------------------------------------------
// -- empty
// --   Empties out the SP
// --------------------------------------------------------------------------
void SP::empty(bool check_merged)
{
   pos->empty();
   app->empty();

   N = 0;
   //if (check_merged && !pixels.isempty())
   //   mexErrMsgTxt("Emptying a SP that hasn't been merged!");
   if (!borders.isempty())
      borders.clear();
   if (!pixels.isempty())
      pixels.clear();
   if (!neighbors.isempty())
      neighbors.clear();
   calculate_log_probs();
}


int SP::get_N()                        { return N;}
double SP::get_log_likelihood_empty()  { return log_likelihood_empty;}
double SP::get_log_likelihood()        { return log_likelihood;}
unsigned long SP::get_UID()            { return UID;}
arr(double) SP::get_prev_v()           { return prev_v;}

void SP::set_flow(arr(double) flow)    { pos->set_offset(flow);}

void SP::set_flow(double flowx, double flowy)
{
   arr(double) flow = allocate_memory<double>(2);
   flow[0] = flowx;
   flow[1] = flowy;
   pos->set_offset(flow);
   deallocate_memory(flow);
}


// --------------------------------------------------------------------------
// -- isempty
// --   returns whether or not the super pixel contains any pixels
// --------------------------------------------------------------------------
bool SP::isempty()
{
   return (N<=0);
}
bool SP::isold()
{
   return is_old;
}

void SP::checkCount()
{
   if (N!=pixels.getLength())
   {
      mexPrintf("N=%d, length=%d\n", N, pixels.getLength());
      mexErrMsgTxt("SHOOT\n");
   }
}


// --------------------------------------------------------------------------
// -- calculate_log_probs
// --   Updates log_likelihood and log_label based on the current parameters
// --------------------------------------------------------------------------
void SP::calculate_log_probs()
{
   if (N<=0)
      log_likelihood = log_likelihood_empty;
   else
      //log_likelihood = pos->calc_logposterior(false, is_old) + app->calc_logposterior();
      log_likelihood = pos->calc_logposterior() + app->calc_logposterior();
}

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
void SP::add_pixel_init(arr(double) data, int index, bool is_border,
   linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp)
{
   N++;
   pos->add_data(data + index*5);

   if (doApp)
      app->add_data(data + index*5+2);

   pixel_ptr = pixels.addNodeEnd(index);
   if (is_border)
      border_ptr = borders.addNodeEnd(index);
   else
      border_ptr = NULL;
}
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
void SP::add_pixel(arr(double) data, int index, bool is_border,
   linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp)
{
   add_pixel_init(data, index, is_border, pixel_ptr, border_ptr, doApp);
   calculate_log_probs();
}


void SP::merge_with(SP *other, arr(int) label, arr(linkedListNode<int>*) border_ptr, int xdim, int ydim)
{
   // update the MIW stuff
   N += other->N;
   pixels.append(other->pixels);
   app->merge_with(other->app);
   pos->merge_with(other->pos);
   other->empty();

   // update the border pixels
   fix_borders(label, border_ptr, xdim, ydim);

   calculate_log_probs();
}

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
void SP::rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp)
{
   if (N<=0)
      mexErrMsgTxt("Trying to remove a pixel from an empty super pixel!");

   N--;
   pos->rem_data(data + index*5);
   if (doApp)
      app->rem_data(data + index*5+2);

   if (border_ptr != NULL)
      borders.deleteNode(border_ptr);
   if (pixel_ptr != NULL)
      pixels.deleteNode(pixel_ptr);
   border_ptr = NULL;
   pixel_ptr = NULL;

   calculate_log_probs();
}
void SP::rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr, bool doApp)
{
   if (N<=0)
      mexErrMsgTxt("Trying to remove a pixel from an empty super pixel!");

   N--;
   pos->rem_data(data + index*5);
   if (doApp)
      app->rem_data(data + index*5+2);

   if (pixel_ptr != NULL)
      pixels.deleteNode(pixel_ptr);
   pixel_ptr = NULL;

   calculate_log_probs();
}



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
void SP::fix_borders(arr(int) label, arr(linkedListNode<int>*) border_ptr, int xdim, int ydim)
{
   // update the border pixels
   borders.clear();
   linkedListNode<int>* node = pixels.getFirst();
   while (node!=NULL)
   {
      int index = node->getData();
      int x = index%xdim;
      int y = index/xdim;
      int cur_label = label[index];

      bool border = false;
      border = border || (x>0 && label[index-1]!=cur_label);
      border = border || (y>0 && label[index-xdim]!=cur_label);
      border = border || (x<xdim-1 && label[index+1]!=cur_label);
      border = border || (y<ydim-1 && label[index+xdim]!=cur_label);

      if (border)
         border_ptr[index] = borders.addNodeEnd(index);
      else
         border_ptr[index] = NULL;

      node = node->getNext();
   }
}




// --------------------------------------------------------------------------
// -- update_neighbors_label_rem
// --   Decrements the corresponding neighbor count for neighbor_label
// --
// --   parameters:
// --     - neighbor_label : the label of the neighbor to decrement
// --------------------------------------------------------------------------
void SP::update_neighbors_label_rem(int neighbor_label)
{
   linkedListNode< std::pair<int, int> >* node = neighbors.getFirst();
   bool found = false;
   while (node!=NULL)
   {
      if (node->getData().first==neighbor_label)
      {
         found = true;
         break;
      }
      node = node->getNext();
   }
   if (found)
   {
      node->applyFunction( &decrement_neighbor_count );
      if (node->getData().second==0)
         neighbors.deleteNode(node);
   }
   else
      mexErrMsgTxt("Trying to remove a neighbor that was never added");
}
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
void SP::update_neighbors_label_rem_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label)
{
   // check to see how many of the neighbors have label neighbor_label
   int x = index%xdim;
   int y = index/xdim;

   int count = 0;
   if (x>0 && label[index-1]==neighbor_label) count++;
   if (y>0 && label[index-xdim]==neighbor_label) count++;
   if (x<xdim-1 && label[index+1]==neighbor_label) count++;
   if (y<ydim-1 && label[index+xdim]==neighbor_label) count++;

   if (count==0)
      update_neighbors_label_rem(neighbor_label);
}

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
void SP::update_neighbors_add_self(arr (int) label, int index, int xdim, int ydim)
{
   int x = index%xdim;
   int y = index/xdim;
   int cur_label = label[index];

   // don't do all neighbors, only the unique ones
   int llabel = (x>0 ? label[index-1] : -1);
   int ulabel = (y>0 ? label[index-xdim] : -1);
   int rlabel = (x<xdim-1 ? label[index+1] : -1);
   int dlabel = (y<ydim-1 ? label[index+xdim] : -1);
   if (llabel>=0 && llabel!=cur_label)
      update_neighbors_label_add(llabel);
   if (ulabel>=0 && ulabel!=cur_label && ulabel!=llabel)
      update_neighbors_label_add(ulabel);
   if (rlabel>=0 && rlabel!=cur_label && rlabel!=llabel && rlabel!=ulabel)
      update_neighbors_label_add(rlabel);
   if (dlabel>=0 && dlabel!=cur_label && dlabel!=llabel && dlabel!=ulabel && dlabel!=rlabel)
      update_neighbors_label_add(dlabel);
}

// --------------------------------------------------------------------------
// -- update_neighbors_self
// --   Updates the neighbor lists and counts by looking at all borders.
// -- Empties previous list. The borders list must be correct!
// --
// --   parameters:
// --     - label : the label image
// --     - (xdim,ydim) : dimensions of image
// --------------------------------------------------------------------------
void SP::update_neighbors_self(arr (int) label, int xdim, int ydim)
{
   linkedListNode< int >* node = borders.getFirst();
   while (node!=NULL)
   {
      int index = node->getData();
      update_neighbors_add_self(label, index, xdim, ydim);
      node = node->getNext();
   }
}


// --------------------------------------------------------------------------
// -- update_neighbors_label_add
// --   Increments the corresponding neighbor count for neighbor_label
// --
// --   parameters:
// --     - neighbor_label : the label of the neighbor to increment
// --------------------------------------------------------------------------
void SP::update_neighbors_label_add(int neighbor_label)
{
   linkedListNode< std::pair<int, int> >* node = neighbors.getFirst();
   bool found = false;
   while (node!=NULL)
   {
      if (node->getData().first==neighbor_label)
      {
         found = true;
         break;
      }
      node = node->getNext();
   }
   if (found)
      node->applyFunction( &increment_neighbor_count );
   else
      neighbors.addNodeEnd( std::pair<int, int>(neighbor_label, 1) );
}
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
void SP::update_neighbors_label_add_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label)
{
   // check to see how many of the neighbors have label neighbor_label
   int x = index%xdim;
   int y = index/xdim;

   int count = 0;
   if (x>0 && label[index-1]==neighbor_label) count++;
   if (y>0 && label[index-xdim]==neighbor_label) count++;
   if (x<xdim-1 && label[index+1]==neighbor_label) count++;
   if (y<ydim-1 && label[index+xdim]==neighbor_label) count++;

   if (count==1)
      update_neighbors_label_add(neighbor_label);
}

bool SP::has_no_neighbors()
{
   return neighbors.isempty();
}


void SP::switch_priors(SP* other)
{
   pos->switch_priors(other->pos);
   app->switch_priors(other->app);

   app->get_mean_mode();

   double tempDouble = log_likelihood_empty;
   log_likelihood_empty = other->log_likelihood_empty;
   other->log_likelihood_empty = tempDouble;

   unsigned long tempLong = UID;
   UID = other->UID;
   other->UID = tempLong;

   bool tempBool = is_old;
   is_old = other->is_old;
   other->is_old = tempBool;
   // MAKE SURE TO SWITC THE NEIGHBORS AFTERWARDS!!!
}




// --------------------------------------------------------------------------
// -- log_likelihood_test_point
// --   finds the change in probability for adding this data point to the super
// -- pixel
// --
// --   parameters
// --     - data : a [1 5] vector of a test point to add
// --------------------------------------------------------------------------
double SP::log_likelihood_test_point_MM(arr(double) data)
{
   // calculate the new log probability
   return pos->calc_logposterior_MM(data) + app->calc_logposterior_MM(data+2);
   //return pos->calc_logposterior_MM(data) + app->calc_logposterior_new(data+2);
}
double SP::log_likelihood_test_point_MM_pos(arr(double) data)
{
   // calculate the new log probability
   return pos->calc_logposterior_MM(data);
}

double SP::log_likelihood_test_point_app(arr(double) data)
{
   return app->calc_logposterior_new(data+2);
}
double SP::log_likelihood_test_point_app_rem(arr(double) data)
{
   return app->calc_logposterior_rem(data+2);
}
double SP::log_likelihood_test_point_pos(arr(double) data)
{
   return pos->calc_logposterior_new(data);
}
double SP::log_likelihood_test_point_pos_rem(arr(double) data)
{
   return pos->calc_logposterior_rem(data);
}

double SP::get_log_likelihood_app()
{
   return app->calc_logposterior(false);
}
double SP::get_log_likelihood_pos()
{
   return pos->calc_logposterior(false);
}


double SP::log_likelihood_switch_app_prior(SP* new_prior)
{
   arr(double) total = app->get_total();
   arr(double) total2 = app->get_total2();
   double logprob = -3*(N+0.5)*0.5*1.837877066409 - new_prior->app->get_sumlogDelta_div2();
   if (!(new_prior->is_old))
   {
      for (int d=0; d<3; d++)
         logprob += (total[d]*total[d]/N - total2[d]) / 2;
   }
   else
   {
      arr(double) Delta = new_prior->get_Delta_app();
      arr(double) theta = new_prior->get_theta_app();
      for (int d=0; d<3; d++)
      {
         double totald = total[d];
         double totald_sq = totald*totald;
         double thetad = theta[d];
         logprob += (Delta[d]*totald_sq + 2*totald*thetad - N*thetad*thetad) / (2*N*Delta[d]+2) - total2[d]/2;
      }
   }
   return logprob;

   /*if (is_uniform)
      return app->get_sumlogDelta_div2() - 0.5*1.837877066409;
   else
   {
      double logprob = app->get_sumlogDelta_div2() - 0.5*1.837877066409;
      arr(double) mean = get_mean_app();
      arr(double) theta = new_prior->get_theta_app();
      arr(double) Delta = new_prior->get_Delta_app();
      for (int d=0; d<3; d++)
      {
         double diff = (mean[d] - theta[d]);
         logprob -= 0.5*diff*diff / (Delta[d]);
      }
      return logprob;
   }*/
}

double SP::log_likelihood_switch_pos_prior(SP* new_prior)
{
   arr(double) total = pos->get_total();
   arr(double) total2 = pos->get_total2();
   double logprob = -(N+0.5)*1.837877066409 - new_prior->pos->get_sumlogDelta_div2();
   if (!(new_prior->is_old))
   {
      for (int d=0; d<2; d++)
         logprob += (total[d]*total[d]/N - total2[d]) / 2;
   }
   else
   {
      arr(double) Delta = new_prior->get_Delta_pos();
      arr(double) theta = new_prior->get_theta_pos();
      for (int d=0; d<2; d++)
      {
         double totald = total[d];
         double totald_sq = totald*totald;
         double thetad = theta[d];
         logprob += (Delta[d]*totald_sq + 2*totald*thetad - N*thetad*thetad) / (2*N*Delta[d]+2) - total2[d]/2;
      }
   }
   return logprob;
}

double SP::log_likelihood_switch_prior(SP* new_prior)
{
   return log_likelihood_switch_pos_prior(new_prior) + log_likelihood_switch_app_prior(new_prior);
}


// --------------------------------------------------------------------------
// -- log_likelihood_test_point
// --   finds the change in probability for adding this data point to the super
// -- pixel
// --
// --   parameters
// --     - data : a [1 5] vector of a test point to add
// --------------------------------------------------------------------------
double SP::log_likelihood_test_point(arr(double) data, bool checkApp)
{
   // calculate the new log probability
   if (checkApp)
      return pos->calc_logposterior_new(data) + app->calc_logposterior_new(data+2);
   else
      return pos->calc_logposterior_new(data) + app->calc_logposterior();
}
double SP::log_likelihood_test_point_rem(arr(double) data, bool checkApp)
{
   // calculate the new log probability
   if (checkApp)
      return pos->calc_logposterior_rem(data) + app->calc_logposterior_rem(data+2);
   else
      return pos->calc_logposterior_new(data) + app->calc_logposterior();
}

// --------------------------------------------------------------------------
// -- log_likelihood_test_merge
// --   finds the change in probability for adding the points in indices to
// -- the super pixel
// --
// --   parameters
// --     - other : another SP that we are testing for a merge
// --------------------------------------------------------------------------
double SP::log_likelihood_test_merge(SP *other)
{
   // calculate the new log probability
   //return pos->calc_logposterior_new(other->pos, is_old) + app->calc_logposterior_new(other->app);
   return pos->calc_logposterior_new(other->pos) + app->calc_logposterior_new(other->app);
}
double SP::log_likelihood_test_merge_pos(SP *other)
{
   // calculate the new log probability
   //return pos->calc_logposterior_new(other->pos, is_old);
   return pos->calc_logposterior_new(other->pos);
}
double SP::log_likelihood_test_merge(SP *other1, SP *other2)
{
   //return pos->calc_logposterior_new(other1->pos, other2->pos, is_old) + app->calc_logposterior_new(other1->app, other2->app);
   return pos->calc_logposterior_new(other1->pos, other2->pos) + app->calc_logposterior_new(other1->app, other2->app);
}

// --------------------------------------------------------------------------
// -- get_XXXX_mode
// --   returns a array to the most likely parameter. returns a static array
// -- that should be explicitly deleted
// --
// --   parameters
// --     - update : whether or not to update the posterior hyperparameters
//--    before returning the mode
// --------------------------------------------------------------------------
arr(double) SP::get_flow()
{
   return pos->get_offset();
   /*arr(double) x = allocate_memory<double>(2);
   x[0] = pos->x[0].get_offset();
   x[1] = pos->x[1].get_offset();
   return x;*/
}
arr(double) SP::get_theta_pos()
{
   return pos->get_theta();
   /*arr(double) x = allocate_memory<double>(2);
   x[0] = pos->x[0].get_theta();
   x[1] = pos->x[1].get_theta();
   return x;*/
}
arr(double) SP::get_theta_app()
{
   return app->get_theta();
   /*arr(double) x = allocate_memory<double>(3);
   x[0] = app->x[0].get_theta();
   x[1] = app->x[1].get_theta();
   x[2] = app->x[2].get_theta();
   return x;*/
}
arr(double) SP::get_Delta_pos()
{
   return pos->get_Delta();
   /*arr(double) x = allocate_memory<double>(2);
   x[0] = pos->x[0].get_Delta();
   x[1] = pos->x[1].get_Delta();
   return x;*/
}
arr(double) SP::get_Delta_app()
{
   return app->get_Delta();
   /*arr(double) x = allocate_memory<double>(3);
   x[0] = app->x[0].get_Delta();
   x[1] = app->x[1].get_Delta();
   x[2] = app->x[2].get_Delta();
   return x;*/
}
// iidmode == true indicates to find the mode of the mean based on current
// pixels (ignoring possible neighbors)
arr(double) SP::get_mean_pos(bool iidmode)
{
   if (iidmode)
      return pos->get_mean_mode();
   else
   {
      mexPrintf("SP shoudln't be here\n");
      return pos->get_mean();
   }
   /*arr(double) x = allocate_memory<double>(2);
   x[0] = pos->x[0].get_mean_mode();
   x[1] = pos->x[1].get_mean_mode();
   return x;*/
}
arr(double) SP::get_mean_app()
{
   return app->get_mean_mode();
   /*arr(double) x = allocate_memory<double>(3);
   x[0] = app->x[0].get_mean_mode();
   x[1] = app->x[1].get_mean_mode();
   x[2] = app->x[2].get_mean_mode();
   return x;*/
}

void SP::set_mean_app(arr(double) new_mean)
{
   app->set_mean(new_mean);
}
void SP::set_mean_pos(arr(double) new_mean)
{
   pos->set_mean(new_mean);
}
void SP::set_meansum_pos(arr(double) new_mean1, arr(double) new_mean2)
{
   pos->set_meansum(new_mean1, new_mean2);
}
void SP::set_flowsum(arr(double) flow, arr(double) other)
{
   pos->set_offsetsum(flow, other);
}

arr(double) SP::get_total_app()
{
   return app->get_total();
}
arr(double) SP::get_total_pos()
{
   return pos->get_total();
}
arr(double) SP::get_total2_pos()
{
   return pos->get_total2();
}
double SP::get_sumlogDelta_div2_pos()
{
   return pos->get_sumlogDelta_div2();
}