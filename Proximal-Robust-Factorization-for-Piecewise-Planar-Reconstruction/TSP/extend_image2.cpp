
#include "helperMEX.h"
#include "binTree.cpp"
#include "distancePixel.cpp"
#include "math.h"
#include "linkedList.cpp"

void repeatpad(int xdim, int ydim, arr(double) imageIn, arr(bool) maskIn,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double farthest_size,
   arr(double) farthest_x, arr(double) farthest_y, arr(double) imageOut);

void repeatpad2(int xdim, int ydim, arr(double) imageIn, arr(bool) maskIn,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double farthest_size,
   arr(double) farthest_x, arr(double) farthest_y, arr(double) imageOut);

void lineBresenham(int x0, int y0, int x1, int y1, int X, int Y, arr(bool) linear);

void calc_closest_curve(int xdim, int ydim, arr(bool) regMaskIn, double band_size,
   arr(double) distance, arr(double) closest_curve_x, arr(double) closest_curve_y,
   double &farthest_size, arr(double) farthest_x, arr(double) farthest_y);

void calc_dist(const int x, const int y, const int xNew, const int yNew,
   int xdim, int ydim, binTree<distancePixel> &unmarked, arr(bool) accepted,
   arr(binNode<distancePixel>*) tree_ptr,  arr(double) distance,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double band_size);


void repeatpad2(int xdim, int ydim, arr(double) imageIn, arr(bool) maskIn,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double farthest_size,
   arr(double) farthest_x, arr(double) farthest_y, arr(double) imageOut)
{
   int i,x,y;
   arr(bool) updatedAngle = allocate_memory<bool>(xdim*ydim);

   for (int i=0; i<xdim*ydim; i++)
   {
      if (maskIn[i])
      {
         imageOut[i] = imageIn[i];
         updatedAngle[i] = true;
      }
      else
      {
         imageOut[i] = 0;
         updatedAngle[i] = false;
      }
   }

   int patterni = 0;
   arr(double) pattern = allocate_memory<double>(20);

   for (i=0; i<farthest_size; i++)
   {
      x = (int)farthest_x[i];
      y = (int)farthest_y[i];
      int index = x + y*xdim;
      int cx = (int)closest_curve_x[index];
      int cy = (int)closest_curve_y[index];

      // check if we need to update it
      if (!updatedAngle[x*ydim+y])
      {
         int dx = x - cx;
         int dy = y - cy;

         double c1 = xdim+ydim;
         double c2 = xdim+ydim;
         if (dx>0) // cx - dx*c = 0;
            c1 = (double)(cx) / (double)(dx);
         else if (dx<0) // cx - dx*c = xdim-1
            c1 = (double)(cx-xdim+1) / (double)(dx);
         if (dy>0) // cy - dy*c = 0;
            c2 = (double)(cy) / (double)(dy);
         else if (dy<0)// cy - dy*c = ydim-1
            c2 = (double)(cy-ydim+1) / (double)(dy);
         double c = (c1 > c2 ? c2 : c1);

         //mexPrintf("%0.2f, %0.2f, %0.2f\n", c1, c2, c);
         int px = cx - dx*c;
         int py = cy - dy*c;

         // draw a line from the point to the closest point, and extend it by a lot... also make it thick
         patterni = 0;
         linkedList<pair> lineCoords;
         int x0 = x;
         int y0 = y;
         int x1 = px;
         int y1 = py;
         int line_dy = y1 - y0;
         int line_dx = x1 - x0;
         int stepx, stepy;

         if (line_dy < 0) { line_dy = -line_dy;  stepy = -1; } else { stepy = 1; }
         if (line_dx < 0) { line_dx = -line_dx;  stepx = -1; } else { stepx = 1; }
         line_dy <<= 1;                                                  // dy is now 2*dy
         line_dx <<= 1;                                                  // dx is now 2*dx

         lineCoords.addNodeEnd(pair(x0, y0));
         if (line_dx > line_dy)
         {
            int fraction = line_dy - (line_dx >> 1);                         // same as 2*dy - dx
            while (x0 != x1)
            {
               if (fraction >= 0)
               {
                  y0 += stepy;
                  fraction -= line_dx;                                // same as fraction -= 2*dx
               }
               x0 += stepx;
               fraction += line_dy;                                    // same as fraction -= 2*dy
               if (!maskIn[x0 + y0*xdim] && patterni==0)
                  lineCoords.addNodeEnd(pair(x0, y0));
               else if (maskIn[x0 + y0*xdim] && patterni<20)
                  pattern[patterni++] = imageIn[x0 + y0*xdim];
               else
                  break;
            }
         }
         else
         {
            int fraction = line_dx - (line_dy >> 1);
            while (y0 != y1)
            {
               if (fraction >= 0)
               {
                  x0 += stepx;
                  fraction -= line_dy;
               }
               y0 += stepy;
               fraction += line_dx;
               if (!maskIn[x0 + y0*xdim] && patterni==0)
                  lineCoords.addNodeEnd(pair(x0, y0));
               else if (maskIn[x0 + y0*xdim] && patterni<20)
                  pattern[patterni++] = imageIn[x0 + y0*xdim];
               else
                  break;
            }
         }

         if (patterni>1)
         {
            // find the target points by tracing through the lines beginning at the closest point
            linkedListNode<pair>* curLine = lineCoords.getFirst();

            int linei = 0;
            while (curLine != NULL)
            {
               pair temp = curLine->getData();
               int curx = temp.xpix;
               int cury = temp.ypix;
               if (!updatedAngle[curx + cury*xdim])
               {
                  updatedAngle[curx + cury*xdim] = true;
                  int ind = (patterni-1-linei)%patterni;
                  if (ind < 0)
                     ind += patterni;
                  imageOut[curx + cury*xdim] = pattern[ind];
                  //imageOut[index] = 255;
               }
               linei++;
               curLine = curLine->getNext();
            }
         }
         else
            imageOut[x + y*xdim] = imageIn[cx + cy*xdim];
      }
   }
   deallocate_memory(pattern);
   deallocate_memory(updatedAngle);
}





void lineBresenham(int x0, int y0, int x1, int y1, int X, int Y, arr(bool) linear)
{
   //if (x1<0 || x1>=X || y1<0 || y1>=Y)
      //mexPrintf("SHIT! (%d, %d) (%d, %d)\n", x0, y0, x1, y1);

   int dy = y1 - y0;
   int dx = x1 - x0;
   int stepx, stepy;

   if (dy < 0) { dy = -dy;  stepy = -1; } else { stepy = 1; }
   if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
   dy <<= 1;                                                  // dy is now 2*dy
   dx <<= 1;                                                  // dx is now 2*dx

   linear[x0 + y0*X] = true;
   if (dx > dy) {
      int fraction = dy - (dx >> 1);                         // same as 2*dy - dx
      while (x0 != x1) {
            if (fraction >= 0) {
               y0 += stepy;
               fraction -= dx;                                // same as fraction -= 2*dx
            }
            x0 += stepx;
            fraction += dy;                                    // same as fraction -= 2*dy
            linear[x0 + y0*X] = true;
      }
   } else {
      int fraction = dx - (dy >> 1);
      while (y0 != y1) {
            if (fraction >= 0) {
               x0 += stepx;
               fraction -= dy;
            }
            y0 += stepy;
            fraction += dx;
            linear[x0 + y0*X] = true;
      }
   }
}



void repeatpad(int xdim, int ydim, arr(double) imageIn, arr(bool) maskIn,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double farthest_size,
   arr(double) farthest_x, arr(double) farthest_y, arr(double) imageOut)
{
   int i,x,y;

   arr(bool) linear = allocate_memory<bool>(xdim*ydim);
   arr(bool) updatedAngle = allocate_memory<bool>(xdim*ydim);

   for (int i=0; i<xdim*ydim; i++)
   {
      linear[i] = false;
      if (maskIn[i])
      {
         imageOut[i] = imageIn[i];
         updatedAngle[i] = true;
      }
      else
      {
         imageOut[i] = 0;
         updatedAngle[i] = false;
      }
   }

   int patterni = 0;
   arr(double) pattern = allocate_memory<double>(20);

   for (i=0; i<farthest_size; i++)
   {
      x = (int)farthest_x[i];
      y = (int)farthest_y[i];
      int index = x + y*xdim;
      int cx = (int)closest_curve_x[index];
      int cy = (int)closest_curve_y[index];

      // check if we need to update it
      if (!updatedAngle[x*ydim+y])
      {
         memset(linear, false, sizeof(bool)*xdim*ydim);
         int dx = x - cx;
         int dy = y - cy;

         double c1 = xdim+ydim;
         double c2 = xdim+ydim;
         if (dx>0) // cx - dx*c = 0;
            c1 = (double)(cx) / (double)(dx);
         else if (dx<0) // cx - dx*c = xdim-1
            c1 = (double)(cx-xdim+1) / (double)(dx);
         if (dy>0) // cy - dy*c = 0;
            c2 = (double)(cy) / (double)(dy);
         else if (dy<0)// cy - dy*c = ydim-1
            c2 = (double)(cy-ydim+1) / (double)(dy);
         double c = (c1 > c2 ? c2 : c1);

         //mexPrintf("%0.2f, %0.2f, %0.2f\n", c1, c2, c);
         int px = cx - dx*c;
         int py = cy - dy*c;

         // draw a line from the point to the closest point, and extend it by a lot... also make it thick
         lineBresenham(x, y, px, py, xdim, ydim, linear);

         int checkx = (dx>0) ? -1 : 1;
         int checky = (dy>0) ? -1 : 1;

         linkedList<pair> lineCoords;
         int curx = x;
         int cury = y;
         bool closest_found = true;
         patterni = 0;

         // loop from the current point back to a point on the level set
         // keep track of these points in lineCoords
         while (!maskIn[curx + cury*xdim])
         {
            lineCoords.addNodeEnd(pair(curx,cury));

            bool cancheckx = (curx + checkx >= 0 && curx + checkx < xdim);
            bool canchecky = (cury + checky >= 0 && cury + checky < ydim);

            if (cancheckx && linear[curx+checkx + cury*xdim])
               curx += checkx;
            else if (canchecky && linear[curx + (cury+checky)*xdim])
               cury += checky;
            else if (cancheckx && canchecky && linear[curx+checkx + (cury+checky)*xdim])
            {
               curx += checkx;
               cury += checky;
            }
            else
            {
               // no closest point found... just copy over it
               imageOut[index] = imageIn[index];
               updatedAngle[index] = true;
               closest_found = false;
               break;
            }
         }
         // now go even further back and try to find a pattern to repeat
         while (closest_found && maskIn[curx + cury*xdim] && patterni<20)
         {
            pattern[patterni] = imageIn[curx + cury*xdim];
            patterni++;

            bool cancheckx = (curx + checkx >= 0 && curx + checkx < xdim);
            bool canchecky = (cury + checky >= 0 && cury + checky < ydim);
            if (cancheckx && linear[curx+checkx + cury*xdim])
               curx += checkx;
            else if (canchecky && linear[curx + (cury+checky)*xdim])
               cury += checky;
            else if (cancheckx && canchecky && linear[curx+checkx + (cury+checky)*xdim])
            {
               curx += checkx;
               cury += checky;
            }
            else
               break;
         }
         if (closest_found && patterni > 1)
         {
            // find the target points by tracing through the lines beginning at the closest point
            linkedListNode<pair>* curLine = lineCoords.getFirst();

            int linei = 0;
            while (curLine != NULL)
            {
               pair temp = curLine->getData();
               curx = temp.xpix;
               cury = temp.ypix;
               if (!updatedAngle[curx + cury*xdim])
               {
                  updatedAngle[curx + cury*xdim] = true;
                  int ind = (patterni-1-linei)%patterni;
                  if (ind < 0)
                     ind += patterni;
                  imageOut[curx + cury*xdim] = pattern[ind];
                  //imageOut[index] = 255;
               }
               linei++;
               curLine = curLine->getNext();
            }
         }
         else
            imageOut[curx + cury*xdim] = imageIn[cx + cy*xdim];
      }
   }
   deallocate_memory(pattern);
   deallocate_memory(updatedAngle);
   deallocate_memory(linear);
}








void calc_closest_curve(int xdim, int ydim, arr(bool) regMaskIn, double band_size,
   arr(double) distance, arr(double) closest_curve_x, arr(double) closest_curve_y,
   double &farthest_size, arr(double) farthest_x, arr(double) farthest_y)
{
   // =======================================================================
   // extend the image
   // =======================================================================
   binTree<distancePixel> unmarked;

   arr(bool) accepted = allocate_memory<bool>(xdim*ydim);
   arr(binNode<distancePixel>*) tree_ptr = allocate_memory< binNode<distancePixel>* >(xdim*ydim);

   for (int x=0; x<xdim; x++)
   {
      for (int y=0; y<ydim; y++)
      {
         int index = x + y*xdim;
         tree_ptr[index] = NULL;

         if (regMaskIn[index])
         {
            closest_curve_x[index] = x;
            closest_curve_y[index] = y;
            accepted[index] = true;
            distance[index] = 0;
         }
         else
         {
            closest_curve_x[index] = 0;
            closest_curve_y[index] = 0;
            accepted[index] = false;
            distance[index] = band_size;
         }
      }
   }

   for (int x=0; x<xdim; x++) for (int y=0; y<ydim; y++)
   {
      int index = x + y*xdim;
      if (regMaskIn[index] &&
          ( ((x-1 >= 0) && (!regMaskIn[index - 1])) ||
            ((x+1 < xdim) && (!regMaskIn[index+1])) ||
            ((y-1 >= 0) && (!regMaskIn[index-xdim])) ||
            ((y+1 < ydim) && (!regMaskIn[index+xdim])) ))
      {
         //tree_ptr[x][y] = unmarked.addLeaf( distancePixel(calc_dist_front(x,y,xdim,ydim,distance), x, y) );
         tree_ptr[index] = unmarked.addLeaf( distancePixel(0, x, y) );
         closest_curve_x[index] = x;
         closest_curve_y[index] = y;
      }
   }


   // recursively pick off the smallest distance until all pixels are done
   linkedList<pair> farthest;
   while (unmarked.treeNotEmpty())
   {
      distancePixel new_point = unmarked.pickOffRoot();
      int x = new_point.coords.xpix;
      int y = new_point.coords.ypix;
      int index = x + y*xdim;
      farthest.addNode(pair(x,y));
      accepted[index] = true;
      tree_ptr[index] = NULL;
      distance[index] = new_point.distance;

      if (x+1<xdim && !accepted[index+1] ) calc_dist(x+1,y,x,y, xdim, ydim, unmarked, accepted, tree_ptr, distance, closest_curve_x, closest_curve_y, band_size);
      if (x-1>=0   && !accepted[index-1] ) calc_dist(x-1,y,x,y, xdim, ydim, unmarked, accepted, tree_ptr, distance, closest_curve_x, closest_curve_y, band_size);
      if (y+1<ydim && !accepted[index+xdim] ) calc_dist(x,y+1,x,y, xdim, ydim, unmarked, accepted, tree_ptr, distance, closest_curve_x, closest_curve_y, band_size);
      if (y-1>=0   && !accepted[index-xdim] ) calc_dist(x,y-1,x,y, xdim, ydim, unmarked, accepted, tree_ptr, distance, closest_curve_x, closest_curve_y, band_size);
   }

   farthest_size = (double)farthest.getLength();
   linkedListNode<pair>* cur = farthest.getFirst();
   for (int i=0; i<farthest_size; i++)
   {
      pair curPair = cur->getData();
      farthest_x[i] = curPair.xpix;
      farthest_y[i] = curPair.ypix;
      cur = cur->getNext();
   }

   deallocate_memory(tree_ptr);
}

// --------------------------------------------------------------------------
// -- calc_dist
// --   calculates distance to zero level set at far point (x,y)
// --
// --   parameters:
// --     - (x,y) : the coordinate where the distance is needed
// --     - (xNew, yNew) : the newly accepted distance that touches (x,y)
// --     - unmarked : the tree of points that have a finite distance but
// --         have not been accepted
// --     - accepted : the 2D boolean matrix saying if the pixel is accepted
// --     - tree_ptr : a 2D matrix of pointers pointing to the node in the
// --         unmarked tree (if exists) or NULL (if doesn't exist)
// --     - band_size : the narrow band size
// --------------------------------------------------------------------------
void calc_dist(const int x, const int y, const int xNew, const int yNew,
   int xdim, int ydim, binTree<distancePixel> &unmarked, arr(bool) accepted,
   arr(binNode<distancePixel>*) tree_ptr,  arr(double) distance,
   arr(double) closest_curve_x, arr(double) closest_curve_y, double band_size)
{
   int index = x + y*xdim;

   // do x direction
   double a = 0;
   double b = 0;
   double c = -1;

   bool checkbck = (x-1>=0) && accepted[index-1];
   bool checkfwd = (x+1<xdim) && accepted[index+1];
   double dist = -1;
   if (checkbck && !checkfwd)
      dist = distance[index-1];
   else if (!checkbck && checkfwd)
      dist = distance[index+1];
   else if (checkbck && checkfwd)
      dist = (distance[index+1] > distance[index-1]) ? distance[index-1] : distance[index+1];
   if (dist != -1)
   {
      a += 1;
      b += -2*dist;
      c += dist*dist;
   }

   // do y direction
   checkbck = (y-1>=0) && accepted[index-xdim];
   checkfwd = (y+1<ydim) && accepted[index+xdim];
   dist = -1;
   if (checkbck && !checkfwd)
      dist = distance[index-xdim];
   else if (!checkbck && checkfwd)
      dist = distance[index+xdim];
   else if (checkbck && checkfwd)
      dist = (distance[index+xdim] > distance[index-xdim]) ? distance[index-xdim] : distance[index+xdim];
   if (dist != -1)
   {
      a += 1;
      b += -2*dist;
      c += dist*dist;
   }

   double sign = (distance[index] > 0) ? 1 : -1;
   dist = (-b + sign * sqrt(b*b - 4*a*c)) / (2*a);
   if (dist < 0) dist = -dist;

   if (dist < distance[index] && dist<=band_size)
   {
      // add it to the tree because it had never been added before
      if (tree_ptr[index] == NULL)
      {
         closest_curve_x[index] = closest_curve_x[xNew + yNew*xdim];
         closest_curve_y[index] = closest_curve_y[xNew + yNew*xdim];
         tree_ptr[index] = unmarked.addLeaf(distancePixel(dist, x, y));
      }
      // modify it from the tree
      else
         unmarked.modifyLeaf(tree_ptr[index], distancePixel(dist, x, y));
      distance[xdim] = dist;
   }
}
