/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   Copyright (C) 2010 William Hart
   
*****************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"

/*
   Determines a "Newton box" containing all the coefficients of a polynomial.
   
   Input: a list of coefficients for a polynomial, poly, and its length, len,
   along with a precision, prec, which all the coefficients are known to.

   Firstly an approximate log of the absolute value of each coefficient is 
   taken. Double precision is fine for this, as we need only an approximation. 
   Call these values a_0, a_1, ..., a_{n-1} where n = len.

   If we plot these values on a graph, we are looking for a set of two parallel
   lines which bound all the points between them, L1 bounding them above, L2
   bounding them below. If the two lines are vertically separated by <= prec
   then the function returns 1 and the slope and y-intercept of L1. If not, the
   function returns 0.

   Such lines are found as follows. Firstly, parallel horizontal lines L1 and L2
   are drawn which bound all the points. Both lines are then rotated (remaining)
   parallel, in the direction which allows them to come closer, until they can be
   rotated no more without one of the points crossing the line.

   ----------------------------------------------------- L1
                                                 X
                    X
   X                                   X
           X

		   X
   ----------------------------------------------------- L2

                                      ------------------ L1
                   ------------------/           X
   ---------------/ X
   X                                   X
           X                                        ---- L2
                                  -----------------/
		   X  -------------------/
  -----------/

*/
  
int _mpfr_poly_bound_newton(double * inter, double * slope, 
			                               mpfr * poly, long len, ulong prec)
{
   /*
	  Precision things can be stored to in a double (can't be less more than
	  prec!
   */
   ulong dprec = FLINT_MIN(53, prec);
	
   /*
	  First get the logarithms of all the points as doubles, rounding up, as
	  we want L1 to miss all the points.
   */
	
   long i, max_i, min_i;
   double slope1, slope2;
   double max_d, min_d;
   double k;
   mpfr_t t;
   mpfr_init2(t, dprec); // don't waste time working at a higher precision
   double * mag = (double *) malloc(sizeof(double)*len); 

   for (i = 0; i < len; i++) 
   {
	  mpfr_abs(t, poly + i, GMP_RNDU);
	  mpfr_log(t, t, GMP_RNDU);
	  mag[i] = mpfr_get_d(t, GMP_RNDU);
   }

   /*
	  Now figure out the max and min logs
   */
   max_i = min_i = 0;
   max_d = min_d = mag[0];
   for (i = 1; i < len; i++)
   {
	   if (mag[i] < min_d) { min_d = mag[i]; min_i = i; }
	   if (mag[i] > max_d) { max_d = mag[i]; max_i = i; }
   }

   /*
      We distinguish two cases:
	    (i)  The point touching the high line L1 is to the right of the point
		     touching the low line L2
		(ii) The point touching the high line L1 is to the left of the point
		     touching the low line L2
		 
		In case (i) we can drop the left of L1 and raise the right of L2.
		In case (ii) we can drop the right of L1 and raise the left of L2.
		If in neither case (i) nor case (ii) we do nothing.
   */
   slope1 = slope2 = 0.0;
   if (max_i > min_i) // case (i)
   {
      if (max_i) // slope to first point
		  slope1 = (max_d - mag[0])/((double) max_i); 
	  
	  if (max_i != len - 1) // slope to last point
		  slope2 = (mag[len-1] - min_d)/((double) (len - 1 - min_i)); 
	  
	  for (i = 1; i < max_i; i++)
	  {
	     k = (max_d - mag[i])/((double) (max_i - i)); // slope from current point
         if (k < slope1) slope1 = k;
	  }
	  
	  for (i = len - 1; i > min_i; i--)
	  {
	     k = (mag[i] - min_d)/((double) (i - min_i)); // slope to current point
         if (k < slope2) slope2 = k;
	  }

	  if (slope1 > slope2)
	     slope1 = slope2;
   } else if (max_i < min_i) // case (ii)
   {
      if (min_i) // slope to first point
		  slope2 = (min_d - mag[0])/((double) min_i); 
	  
	  if (max_i != len - 1) // slope to last point
		  slope1 = (mag[len-1] - max_d)/((double) (len - 1 - max_i)); 
	  
	  for (i = 1; i < min_i; i++)
	  {
	     k = (min_d - mag[i])/((double) (min_i - i)); // slope from current point
         if (k > slope2) slope2 = k;
	  }
	  
	  for (i = len - 1; i > max_i; i--)
	  {
	     k = (mag[i] - max_d)/((double) (i - max_i)); // slope to current point
         if (k > slope1) slope1 = k;
	  }

	  if (slope1 < slope2)
	     slope1 = slope2;
   }

   (*slope) = slope1;
   (*inter) = max_d - (double) max_i*slope1;

   // if inter(L2) < inter(L1) - prec 
   if (min_d - (double) min_i*slope2 <= *inter - prec) 
	  return 0;

   return 1;
   
   mpfr_clear(t);
}
